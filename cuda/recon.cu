/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "wfs.h"
#include "recon.h"
#include "pcg.h"
#include "curmat.h"
#include "cucmat.h"
#include "accphi.h"
#undef TIMING
#define TIMING 0
#if !TIMING
#define TIC_test
#define tic_test
#define toc_test(A)
#else
#define TIC_test TIC
#define tic_test tic
#define toc_test(A) toc2(A);tic
#endif
__global__ static void saloc2ptr_do(int (*restrict saptr)[2], float (*restrict saloc)[2], 
				    int nsa, float ox, float oy, float dx, float dy){
    const int step=blockDim.x * gridDim.x;
    const float dx1=1./dx;
    const float dy1=1./dy;
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	saptr[isa][0]=(int)roundf((saloc[isa][0]-ox)*dx1);
	saptr[isa][1]=(int)roundf((saloc[isa][1]-oy)*dy1);
    }
}
W01_T *gpu_get_W01(dsp *R_W0, dmat *R_W1){
    if(!R_W0 || !R_W1){
	error("R0, R1 must not be empty\n");
    }
    W01_T *W01=(W01_T*)calloc(1, sizeof(W01_T));
    cp2gpu(&W01->W1, R_W1);
    {
	/*W0 of partially illuminates subaps are stored as sparse matrix in
	  GPU. W0 of fully illuminated subaps are not.*/
	spint *pp=R_W0->p;
	spint *pi=R_W0->i;
	double *px=R_W0->x;
	dsp *W0new=spnew(R_W0->m, R_W0->n, R_W0->nzmax);
	spint *pp2=W0new->p;
	spint *pi2=W0new->i;
	double *px2=W0new->x;
	int *full;
	cudaMallocHost(&full, R_W0->n*sizeof(int));
	//#define W0_BW 1
	double W1max=dmax(R_W1);
	double thres=W1max*(1.f-1e-6);
	W01->W0v=(float)(W1max*4./9.);//max of W0 is 4/9 of max of W1. 
	info("W0v=%g\n", W01->W0v);
	int count=0;
	int count2=0;
	for(int ic=0; ic<R_W0->n; ic++){
	    pp2[ic]=count;
	    if(R_W1->p[ic]>thres){
		full[count2]=ic;
		count2++;
	    }else{
		int nv=pp[ic+1]-pp[ic];
		memcpy(pi2+count, pi+pp[ic], sizeof(spint)*nv);
		memcpy(px2+count, px+pp[ic], sizeof(double)*nv);
		count+=nv;
	    }
	}
	pp2[R_W0->n]=count;
	W0new->nzmax=count;
	cp2gpu(&W01->W0p, W0new);
	cp2gpu(&W01->W0f, full, count2);
	W01->nW0f=count2;
	spfree(W0new);
	cudaFreeHost(full);
    }
    return W01;
}
/*
  The caller must specify current GPU.
*/
static void gpu_setup_recon_do(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    if(!cudata->recon){
	cudata->recon=new curecon_t;
    }
    curecon_t *curecon=cudata->recon;
    if(parms->recon.mvm && (!parms->gpu.tomo || !parms->gpu.fit)){
	return; /*Use CPU to assemble MVM*/
    }
    if(parms->recon.alg!=0){
	error("Only MVR is implemented in GPU\n");
    }
    cuwloc_t *cupowfs=cudata->powfs;
    if((parms->gpu.tomo || parms->gpu.fit) && !parms->sim.idealfit){
	curecon->opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
	curecon->opdr_vec=curcellnew(recon->npsr, 1);
	for(int ips=0; ips<recon->npsr; ips++){
	    curecon->opdr_vec->p[ips]=curecon->opdr->p[ips]->ref(1);
	}
    }
    if(parms->gpu.tomo || parms->gpu.fit){
	curecon->amap=new cugrid_t[parms->ndm];
	curecon->xmap=new cugrid_t[recon->npsr];
	for(int idm=0; idm<parms->ndm; idm++){
	    curecon->amap[idm].init(recon->amap[idm]);
	}
	for(int ipsr=0; ipsr<recon->npsr; ipsr++){
	    curecon->xmap[ipsr].init(recon->xmap[ipsr]);
	}
	if(parms->fit.cachedm){
	    curecon->acmap=new cumap_t[parms->ndm];
	    for(int idm=0; idm<parms->ndm; idm++){
		curecon->acmap[idm].init(recon->acmap[idm]);
	    }
	}
	if(parms->fit.cachex){
	    curecon->xcmap=new cugrid_t[recon->npsr];
	    for(int ipsr=0; ipsr<recon->npsr; ipsr++){
		curecon->xcmap[ipsr].init(recon->xcmap[ipsr]);
	    }
	}
	curecon->pmap.init(recon->pmap);
	curecon->fmap.init(recon->fmap);
    }
    if(parms->gpu.tomo){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].skip) continue;
	    int nsa=powfs[ipowfs].pts->nsa;
	    cudaMalloc(&cupowfs[ipowfs].saptr, nsa*2*sizeof(int));
	    saloc2ptr_do<<<DIM(nsa,256)>>>
		(cupowfs[ipowfs].saptr, cupowfs[ipowfs].saloc, nsa, 
		 recon->pmap->ox, recon->pmap->oy, recon->pmap->dx, recon->pmap->dy);
	    if(recon->GP->p[ipowfs]){
		const int use_mat=parms->tomo.pos==2 ||parms->tomo.pos==1 ;
		if(use_mat){//normally true
		    dsp *GP=sptrans(recon->GP->p[ipowfs]);
		    spint *pp=GP->p;
		    spint *pi=GP->i;
		    double *px=GP->x;
		    //convert the max float to max 2 byte integer
		    double pxscale=floor(32767./maxabs(px, GP->nzmax));
		    int np1=parms->tomo.pos+1;
		    int np=np1*np1;
		    int zmax=parms->tomo.pos;
		    short2 *partxy=(short2*)calloc(sizeof(short2),np*nsa);//need to zero memory
		    int nsa=powfs[ipowfs].pts->nsa;
		    double dx1=1./recon->ploc->dx;
		    double dy1=1./recon->ploc->dy;
		    for(int ic=0; ic<GP->n; ic++){
			int isa=(ic<nsa)?ic:(ic-nsa);
			for(spint ir=pp[ic]; ir<pp[ic+1]; ir++){
			    int ix=pi[ir];
			    double lx=recon->ploc->locx[ix];
			    double ly=recon->ploc->locy[ix];
			    double sx=powfs[ipowfs].saloc->locx[isa];
			    double sy=powfs[ipowfs].saloc->locy[isa];
			    int zx=(int)round((lx-sx)*dx1);
			    int zy=(int)round((ly-sy)*dy1);
			    /**
			       When the points used to generate GP align well
			       with the subaperture edge, the coupled points are
			       confined within the subaperture.
			    */
			    if(zx<0 || zx>zmax || zy<0 || zy>zmax){
				warning("isa=%d, zxy=%d %d\n", isa, zx, zy);
			    }
			    if(zx<0) zx=0;
			    if(zx>zmax) zx=zmax;
			    if(zy<0) zy=0;
			    if(zy>zmax) zy=zmax;
			    if(ic<nsa){
				partxy[np*isa+zx+zy*np1].x+=(short)round(px[ir]*pxscale);
			    }else{
				partxy[np*isa+zx+zy*np1].y+=(short)round(px[ir]*pxscale);
			    }
			}
		    }
		    cupowfs[ipowfs].GPp=new cumat<int>(np, nsa);
		    cudaMemcpy(cupowfs[ipowfs].GPp->p, partxy, sizeof(int)*np*nsa, cudaMemcpyHostToDevice);
		    cupowfs[ipowfs].GPscale=1./pxscale;
		    free(partxy);
		    spfree(GP);
		}else{/*use sparse */
		    cp2gpu(&cupowfs[ipowfs].GP, recon->GP->p[ipowfs]);
		}
	    }else{
		error("GP is required\n");
	    }
	}
 
	curecon->l2c=(float*)calloc(recon->npsr, sizeof(float));
	for(int ips=0; ips<recon->npsr; ips++){
	    float tmp=laplacian_coef(recon->r0, recon->wt->p[ips], recon->xmap[ips]->dx)*0.25f;
	    curecon->l2c[ips]=tmp*tmp*TOMOSCALE;
	}
	if(parms->tomo.piston_cr){
	    curecon->zzi=(int*)calloc(recon->npsr, sizeof(int));
	    curecon->zzv=(float*)calloc(recon->npsr, sizeof(float));
	    for(int ips=0; ips<recon->npsr; ips++){
		double r0=recon->r0;
		double dx=recon->xloc[ips]->dx;
		double wt=recon->wt->p[ips];
		int icenter=loccenter(recon->xloc[ips]);
		curecon->zzi[ips]=icenter;
		curecon->zzv[ips]=pow(laplacian_coef(r0,wt,dx),2)*TOMOSCALE*1e-6;
	    }
	}
	curecon->neai=curcellnew(parms->nwfsr, 1);
	/*convert recon->saneai to our format. */
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    int nsa=powfs[ipowfs].pts->nsa;
	    int iwfs0=parms->recon.glao?iwfs:parms->powfs[ipowfs].wfs[0];/*first wfs in this group. */
	    if(iwfs!=iwfs0 && recon->saneai->p[iwfs+iwfs*parms->nwfsr]->p
	       ==recon->saneai->p[iwfs0+iwfs0*parms->nwfsr]->p){
		curecon->neai->p[iwfs]=curecon->neai->p[iwfs0]->ref();
	    }else{
		dsp *nea=recon->saneai->p[iwfs+iwfs*parms->nwfsr];
		spint *pp=nea->p;
		spint *pi=nea->i;
		double *px=nea->x;
		
		float (*neai)[3]=(float(*)[3])calloc(3*nsa, sizeof(float));
		if(nea->n!=2*nsa) error("nea doesn't have 2nsa x 2nsa dimension\n");
		for(int ic=0; ic<nea->n; ic++){
		    for(spint ir=pp[ic]; ir<pp[ic+1]; ir++){
			int ix=pi[ir];
			int isa=ic<nsa?ic:ic-nsa;
			float val=(float)px[ir]*TOMOSCALE;
			if(ix==ic){/*diagonal part. */
			    if(ic==isa){/*x */
				neai[isa][0]=val;
			    }else{/*y */
				neai[isa][1]=val;
			    }
			}else if(ix==ic-nsa || ix==ic+nsa){/*cross part. symmetric. */
			    neai[isa][2]=val;
			}else{
			    error("saneai has invalid format\n");
			}
		    }
		}
		curecon->neai->p[iwfs]=curnew(3, nsa);
		DO(cudaMemcpy(curecon->neai->p[iwfs]->p, neai, 3*nsa*sizeof(float), cudaMemcpyHostToDevice));
		free(neai);
	    }
	}/*for iwfs */
	CUDA_SYNC_DEVICE;
	if(recon->PTT && !curecon->PTT){//for t/t proj in 1)uplink t/t 2) recon
	    cp2gpu(&curecon->PTT, recon->PTT);
	}
	if(parms->tomo.precond==1){
	    gpu_setup_recon_fdpcg(parms, recon);
	}
	if(parms->tomo.alg==0){//CBS
	    chol_convert(recon->RL.C, 0);
	    cp2gpu(&curecon->RCl, recon->RL.C->Cl);
	    cp2gpu(&curecon->RCp, recon->RL.C->Cp, recon->RL.C->Cl->m);
	    if(recon->RL.Up){
		cp2gpu(&curecon->RUp, recon->RL.Up);
		cp2gpu(&curecon->RVp, recon->RL.Vp);
	    }
	}else if(parms->tomo.alg==2){//SVD
	    cp2gpu(&curecon->RMI, recon->RL.MI);
	}

	const int nwfs=parms->nwfsr;
	int nxp=recon->pmap->nx;
	int nyp=recon->pmap->ny;
  
	int nxpw[nwfs], nypw[nwfs], ngw[nwfs];
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip){
		nxpw[iwfs]=0;
		nypw[iwfs]=0;
		ngw[iwfs]=0;
	    }else{
		nxpw[iwfs]=nxp;
		nypw[iwfs]=nyp;
		ngw[iwfs]=powfs[ipowfs].pts->nsa*2;
	    }
	}
	curecon->opdwfs=curcellnew(nwfs, 1, nxpw, nypw);
	curecon->grad=curcellnew(nwfs, 1, ngw, (int*)NULL);
	curecon->ttf=curnew(3*nwfs, 1);

	GPU_PROP_GRID_T *hxdata=new GPU_PROP_GRID_T[nwfs*recon->npsr];
	DO(cudaMalloc(&curecon->hxdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr));
	for(int ips=0; ips<recon->npsr; ips++){ 
	    const float ht=recon->ht->p[ips]; 
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		const int ipowfs = parms->wfsr[iwfs].powfs;
		if(!parms->powfs[ipowfs].skip){
		    const float hs = parms->powfs[ipowfs].hs; 
		    const float scale = 1.f - ht/hs; 
		    float dispx=parms->wfsr[iwfs].thetax*ht; 
		    float dispy=parms->wfsr[iwfs].thetay*ht; 
		    cugrid_t pmapscale=curecon->pmap.scale(scale);
		    gpu_prop_grid_prep(hxdata+iwfs+ips*nwfs, pmapscale, curecon->xmap[ips],
				       dispx, dispy, NULL); 
		    {
			float tmp=laplacian_coef(recon->r0, recon->wt->p[ips], recon->xmap[ips]->dx)*0.25f;
			hxdata[iwfs+ips*nwfs].l2c=tmp*tmp*TOMOSCALE;
			if(parms->tomo.piston_cr){
			    hxdata[iwfs+ips*nwfs].zzi=loccenter(recon->xloc[ips]);
			    hxdata[iwfs+ips*nwfs].zzv=tmp*tmp*TOMOSCALE*1e-6;
			}else{
			    hxdata[iwfs+ips*nwfs].zzi=-1;
			}
		    }
		}
		hxdata[iwfs+ips*nwfs].togpu(&curecon->hxdata[iwfs+ips*nwfs]);
	    }
	}
	delete [] hxdata;
	GPU_GP_T *gpdata=new GPU_GP_T[nwfs];

	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    if(parms->powfs[ipowfs].wfs[0]!=0){
		error("Check this case. We had assumption that this powfs is the first group.\n");
	    }
	    gpdata[iwfs].ipowfs=ipowfs;
	    gpdata[iwfs].nwfs=parms->powfs[ipowfs].nwfsr;
	    gpdata[iwfs].jwfs=parms->powfs[ipowfs].wfsind[iwfs];//wfs index in this group
	    gpdata[iwfs].saptr=cupowfs[ipowfs].saptr;
	    gpdata[iwfs].dsa=powfs[ipowfs].pts->dsa;
	    gpdata[iwfs].GPp=(short2*)cupowfs[ipowfs].GPp->p;
	    gpdata[iwfs].GPscale=cupowfs[ipowfs].GPscale;
	    gpdata[iwfs].pos=parms->tomo.pos;
	    if(curecon->PTT){
		gpdata[iwfs].PTT=curecon->PTT->p[iwfs+iwfs*nwfs]->p;
	    }
	    if(parms->powfs[ipowfs].dfrs){
		/*We only use the first diagonal block for each powfs. The
		  off diagonal is simply -0.2 times the diagonal block*/
		int iwfs0=parms->powfs[ipowfs].wfs[0];//first wfs
		int iwfs1=parms->powfs[ipowfs].wfs[1];//second wfs
		if(!curecon->PDF){
		    curecon->PDF=curcellnew(nwfs, 1);
		}
		if(iwfs==iwfs0){//not the first one.
		    cp2gpu(&curecon->PDF->p[iwfs], recon->PDF->p[iwfs1*nwfs+iwfs1]);
		}
		gpdata[iwfs].PDF=curecon->PDF->p[iwfs0]->p;//every one in this group.
		if(curecon->PTT){
		    /*coupling between TT and DF modes. 
		      We desire (I-DF*PDF)(I-TT*PTT)g=(I-TT*PTT-DF*PDF+DF*PDF*TT*PTT)g
		      So we first compute tt=PTT*g; df=PDF*g; then
		      g2=(I-TT*tt-DF*(df-(PDF*TT)*tt))
		      Here we record the values of PDF*TT
		    */
		    dcell *pdftt=NULL;
		    dcellmm(&pdftt, recon->PDF, recon->TT, "nn", 1);
		    if(!curecon->PDFTT){
			curecon->PDFTT=curcellnew(nwfs, 1);
		    }
		    if(iwfs==iwfs0){
			cp2gpu(&curecon->PDFTT->p[iwfs], pdftt->p[iwfs1*nwfs+iwfs1]);
		    }
		    gpdata[iwfs].PDFTT=curecon->PDFTT->p[iwfs0]->p;
		    dcellfree(pdftt);
		}
	    }
	    gpdata[iwfs].neai=(const float(*)[3])curecon->neai->p[iwfs]->p;
	    gpdata[iwfs].nsa=powfs[ipowfs].pts->nsa;
	    gpdata[iwfs].nxp=recon->pmap->nx;
	    gpdata[iwfs].dxp=recon->pmap->dx;
	    gpdata[iwfs].dyp=recon->pmap->dy;
	    gpdata[iwfs].oxp=recon->pmap->ox;
	    gpdata[iwfs].oyp=recon->pmap->oy;
	}
	DO(cudaMalloc(&curecon->gpdata, sizeof(GPU_GP_T)*nwfs));
	DO(cudaMemcpy(curecon->gpdata, gpdata, sizeof(GPU_GP_T)*nwfs, cudaMemcpyHostToDevice));
	delete [] gpdata;
    }
    if(parms->gpu.fit){
	long npsr=recon->npsr;
	long ndir=parms->fit.nfit;
	long ndm=parms->ndm;
	if(parms->gpu.fit==1){ /*For fitting using sparse matrix*/
	    cp2gpu(&curecon->FR, &recon->FR);
	    cp2gpu(&curecon->FL, &recon->FL);
	    curecon->dmfit=curcellnew(ndm, 1, recon->anloc, (long*)NULL);
	    curecon->fitrhs=curcellnew(ndm, 1, recon->anloc, (long*)NULL);
	    curecon->dmfit_vec=curecon->dmfit;
	}else if(parms->gpu.fit==2){ /*For fitting using ray tracing*/
	    if(!recon->W0 || !recon->W1){
		error("W0, W1 is required\n");
	    }
	    if(parms->sim.idealfit){
		cp2gpu(&curecon->floc, recon->floc);
		curecon->nfloc=recon->floc->nloc;
	    }
	    curecon->W01=gpu_get_W01(recon->W0, recon->W1);
	    if(recon->fitNW){
		dmat *fitNW=dcell2m(recon->fitNW);
		cp2gpu(&curecon->fitNW, fitNW);
		dfree(fitNW);
		curecon->dotNW=curnew(curecon->fitNW->ny, 1);
	    }
	    if(recon->actslave){
		cp2gpu(&curecon->actslave, recon->actslave);
	    }
	    curecon->dmfit=curcellnew(ndm, 1, recon->anx, recon->any);
	    curecon->dmfit_vec=curcellnew(ndm, 1);
	    for(int idm=0; idm<ndm; idm++){
		curecon->dmfit_vec->p[idm]=curecon->dmfit->p[idm]->ref(1);
	    }
	    curecon->fitrhs=curcellnew(ndm, 1, recon->anx, recon->any);
	    if(parms->fit.cachedm){
		long acnx[ndm], acny[ndm];
		for(int idm=0; idm<ndm; idm++){
		    acnx[idm]=curecon->acmap[idm].nx;
		    acny[idm]=curecon->acmap[idm].ny;
		}
		curecon->dmcache=curcellnew(ndm, 1, acnx, acny);
	    }
	    if(parms->fit.cachex){
		long xcnx[npsr], xcny[npsr];
		for(int ips=0; ips<npsr; ips++){
		    xcnx[ips]=curecon->xcmap[ips].nx;
		    xcny[ips]=curecon->xcmap[ips].ny;
		}
		curecon->xcache=curcellnew(npsr, 1, xcnx, xcny);
	    } 
	    cp2gpu(&curecon->fitwt, recon->fitwt);
	    curecon->cubic_cc=curcellnew(ndm, 1);
	    for(int idm=0; idm<ndm; idm++){
		if(parms->dm[idm].cubic){
		    curecon->cubic_cc->p[idm]=gpu_dmcubic_cc(parms->dm[idm].iac);
		}
	    }
	    //xloc -> floc
	    GPU_PROP_GRID_T *hxpdata=new GPU_PROP_GRID_T[npsr*ndir];
	    //dm -> floc
	    GPU_PROP_GRID_T *hadata=new GPU_PROP_GRID_T[ndm*ndir];
	    DO(cudaMalloc(&curecon->hxpdata, sizeof(GPU_PROP_GRID_T)*npsr*ndir));
	    DO(cudaMalloc(&curecon->hadata, sizeof(GPU_PROP_GRID_T)*ndm*ndir));
	    //dm: amap->acmap
	    GPU_PROP_GRID_T *ha0data=NULL, *ha1data=NULL;
	    if(parms->fit.cachedm){
		ha0data=new GPU_PROP_GRID_T[ndm];
		ha1data=new GPU_PROP_GRID_T[ndm*ndir];
		DO(cudaMalloc(&curecon->ha0data, sizeof(GPU_PROP_GRID_T)*ndm));
		DO(cudaMalloc(&curecon->ha1data, sizeof(GPU_PROP_GRID_T)*ndm*ndir));
	    }
	    GPU_PROP_GRID_T *hxp0data=NULL, *hxp1data=NULL;
	    if(parms->fit.cachex){
		hxp0data=new GPU_PROP_GRID_T[npsr];
		hxp1data=new GPU_PROP_GRID_T[npsr*ndir];
		DO(cudaMalloc(&curecon->hxp0data, sizeof(GPU_PROP_GRID_T)*npsr));
		DO(cudaMalloc(&curecon->hxp1data, sizeof(GPU_PROP_GRID_T)*npsr*ndir));
	    }


	    for(int ipsr=0; ipsr<npsr; ipsr++){
		const float ht=recon->ht->p[ipsr];
		if(parms->fit.cachex){
		    gpu_prop_grid_prep(hxp0data+ipsr, curecon->xcmap[ipsr], curecon->xmap[ipsr],
				       0,0, NULL);
		    hxp0data[ipsr].togpu(&curecon->hxp0data[ipsr]);
		}
		for(int idir=0; idir<ndir; idir++){
		    const float hs=parms->fit.hs[idir];
		    const float thetax=(float)parms->fit.thetax[idir];
		    const float thetay=(float)parms->fit.thetay[idir];
		    const float scale=1.f-ht/hs;
		    cugrid_t fmapscale=curecon->fmap*scale;
		    gpu_prop_grid_prep(hxpdata+idir+ipsr*ndir, fmapscale, curecon->xmap[ipsr],
				       thetax*ht, thetay*ht, NULL);
		    hxpdata[idir+ipsr*ndir].togpu(&curecon->hxpdata[idir+ipsr*ndir]);
		    if(parms->fit.cachex){
			gpu_prop_grid_prep(hxp1data+idir+ipsr*ndir, fmapscale, curecon->xcmap[ipsr],
					   thetax*ht, thetay*ht, NULL);
			hxp1data[idir+ipsr*ndir].togpu(&curecon->hxp1data[idir+ipsr*ndir]);
		    }
		}
	    }
	    for(int idm=0; idm<ndm; idm++){
		const float ht=parms->dm[idm].ht;
		if(parms->fit.cachedm){
		    gpu_prop_grid_prep(ha0data+idm, curecon->acmap[idm], curecon->amap[idm],
				       0, 0, curecon->cubic_cc->p[idm]);
		    ha0data[idm].togpu(&curecon->ha0data[idm]);
		}
		for(int idir=0; idir<ndir; idir++){
		    const float hs=parms->fit.hs[idir];
		    const float thetax=(float)parms->fit.thetax[idir];
		    const float thetay=(float)parms->fit.thetay[idir];
		    const float scale=1.f-ht/hs;
		    cugrid_t fmapscale=curecon->fmap*scale;
		    if(parms->fit.cachedm){
			gpu_prop_grid_prep(ha1data+idir+idm*ndir, fmapscale, curecon->acmap[idm],
					   thetax*ht, thetay*ht, NULL);
			ha1data[idir+idm*ndir].togpu(&curecon->ha1data[idir+idm*ndir]);
		    }
		    gpu_prop_grid_prep(hadata+idir+idm*ndir, fmapscale, curecon->amap[idm],
				       thetax*ht, thetay*ht, curecon->cubic_cc->p[idm]);	
		    hadata[idir+idm*ndir].togpu(&curecon->hadata[idir+idm*ndir]);
		}
	    }
	    delete [] hxpdata;
	    delete [] hxp0data;
	    delete [] hxp1data;
	    delete [] hadata;
	    delete [] ha0data;
	    delete [] ha1data;
	}

	if(parms->fit.alg==0){
	    chol_convert(recon->FL.C, 0);
	    cp2gpu(&curecon->FCl, recon->FL.C->Cl);
	    cp2gpu(&curecon->FCp, recon->FL.C->Cp, recon->FL.C->Cl->m);
	    if(recon->FL.Up){
		cp2gpu(&curecon->FUp, recon->FL.Up);
		cp2gpu(&curecon->FVp, recon->FL.Vp);
	    }
	}else if(parms->fit.alg==2){
	    cp2gpu(&curecon->FMI, recon->FL.MI);
	}
	const int nfit=parms->fit.nfit;
	long fnx[nfit],fny[nfit];
	for(int ifit=0; ifit<nfit; ifit++){
	    fnx[ifit]=recon->fmap->nx;
	    fny[ifit]=recon->fmap->ny;
	}
	curecon->opdfit=curcellnew(nfit, 1, fnx, fny);
	curecon->opdfit2=curcellnew(nfit, 1, fnx, fny);
	curecon->opdfitv=curnew(recon->fmap->nx*recon->fmap->ny, nfit, curecon->opdfit->m->p, 0);
	curecon->opdfit2v=curnew(recon->fmap->nx*recon->fmap->ny, nfit, curecon->opdfit2->m->p, 0);
	curecon->pis=curnew(1, parms->fit.nfit);
    }
 
    if(recon->RFlgsx){
	cp2gpu(&curecon->RFlgsx, recon->RFlgsx);
    }
    if(recon->RFngsx){
	cp2gpu(&curecon->RFngsx, recon->RFngsx);
    }
    if(recon->RFdfx){
	cp2gpu(&curecon->RFdfx, recon->RFdfx);
    }
    if(recon->GXL){
	cp2gpu(&curecon->GXL, recon->GXL);
    }
    gpu_print_mem("recon init");
}
void gpu_setup_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    if(parms->recon.mvm && parms->gpu.tomo && parms->gpu.fit && !parms->load.mvm){
	for(int igpu=0; igpu<NGPU; igpu++){
	    gpu_set(igpu);
	    gpu_setup_recon_do(parms, powfs, recon);
	}
    }else{
	gpu_set(gpu_recon);
	gpu_setup_recon_do(parms, powfs, recon);
    }
}
/*Copy FDPCG information to GPU. This may be done every time step if prediction is used.*/
static void gpu_setup_recon_fdpcg_do(const PARMS_T *parms, RECON_T *recon){
    if(!parms->tomo.precond==1) return;
    curecon_t *curecon=cudata->recon;
    FDPCG_T *fdpcg=recon->fdpcg;
    if(!fdpcg){
	return;
    }
    int bs=fdpcg->bs;
    int nb=(fdpcg->nbx/2+1)*fdpcg->nby;//half frequency range
    if(curecon->fdpcg){//already initialized. Just update Mb
	int nxsave=fdpcg->Mbinv->nx;
	fdpcg->Mbinv->nx=nb;
	cp2gpu(&curecon->fdpcg->Mb, fdpcg->Mbinv);
	fdpcg->Mbinv->nx=nxsave;
	return;
    }
    cufdpcg_t *cufd=curecon->fdpcg=new cufdpcg_t;
    cufd->scale=fdpcg->scale;
    cp2gpu(&cufd->perm, fdpcg->permhf, nb*bs);
    //copy only needed blocks to gpu
    int nxsave=fdpcg->Mbinv->nx;
    fdpcg->Mbinv->nx=nb;
    cp2gpu(&cufd->Mb, fdpcg->Mbinv);
    fdpcg->Mbinv->nx=nxsave;
    int nps=recon->npsr;
    int count=0;
    int osi=-1;
    int start[nps];
    for(int ips=0; ips<nps; ips++){
	/*group layers with the same os together in a batch fft.*/
	if(osi != parms->atmr.os[ips]){
	    start[count]=ips;
	    osi = parms->atmr.os[ips];
	    count++;
	}
    }
    cufd->fft=(cufftHandle*)calloc(count, sizeof(cufftHandle));
    cufd->ffti=(cufftHandle*)calloc(count, sizeof(cufftHandle));
    cufd->fftnc=count;
    cufd->fftips=(int*)calloc(count+1, sizeof(int));
    for(int ic=0; ic<count; ic++){
	cufd->fftips[ic]=start[ic];
    }
    cufd->fftips[count]=nps;
    for(int ic=0; ic<count; ic++){
	int ncomp[2];
	/*Notice the reverse in specifying dimensions. THe first element is outmost rank.*/
	ncomp[0]=recon->xny[start[ic]];
	ncomp[1]=recon->xnx[start[ic]];

	int nembed[2];
	nembed[0]=recon->xnx[start[ic]]*recon->xny[start[ic]];
	nembed[1]=recon->xnx[start[ic]];
	DO(cufftPlanMany(&cufd->fft[ic], 2, ncomp, 
			 nembed, 1, ncomp[0]*ncomp[1], 
			 nembed, 1, ncomp[0]*ncomp[1], 
			 CUFFT_R2C, cufd->fftips[ic+1]-cufd->fftips[ic]));
	DO(cufftPlanMany(&cufd->ffti[ic], 2, ncomp, 
			 nembed, 1, ncomp[0]*ncomp[1], 
			 nembed, 1, ncomp[0]*ncomp[1],
			 CUFFT_C2R, cufd->fftips[ic+1]-cufd->fftips[ic]));
	DO(cufftSetStream(cufd->ffti[ic], curecon->cgstream));

	DO(cufftSetStream(cufd->fft[ic], curecon->cgstream));
    }
    cufd->xhat1=cuccellnew(recon->npsr, 1, recon->xnx, recon->xny);
    {
	int nby=256/bs;//number of blocks in each grid
	int nbz=nb/nby;//number of grids to launch.
	while(nb!=nbz*nby){
	    nby--;
	    nbz=nb/nby;
	}
	cufd->nby=nby;
	cufd->nbz=nbz;
    }
    /* notice: performance may be improved by using
       R2C FFTs instead of C2C. Need to update perm
       and Mbinv to use R2C.*/
    GPU_FDPCG_T *fddata=new GPU_FDPCG_T[nps];
    for(int ips=0; ips<nps; ips++){
	fddata[ips].nx=recon->xnx[ips];
	fddata[ips].ny=recon->xny[ips];
	if(cufd->scale){
	    fddata[ips].scale=1.f/sqrtf((float)(recon->xnx[ips]*recon->xny[ips]));
	}else{
	    fddata[ips].scale=1.f;
	}
    }
    cudaMalloc(&curecon->fddata, sizeof(GPU_FDPCG_T)*nps);
    cudaMemcpy(curecon->fddata, fddata, sizeof(GPU_FDPCG_T)*nps, cudaMemcpyHostToDevice);
    delete [] fddata;
}
void gpu_setup_recon_fdpcg(const PARMS_T *parms, RECON_T *recon){
    if(parms->recon.mvm && parms->gpu.tomo && parms->gpu.fit && !parms->load.mvm){
	for(int igpu=0; igpu<NGPU; igpu++){
	    gpu_set(igpu);
	    gpu_setup_recon_fdpcg_do(parms, recon);
	}
    }else{
	gpu_set(gpu_recon);
	gpu_setup_recon_fdpcg_do(parms, recon);
    }
}
void gpu_recon_free_do(){
    curecon_t *curecon=cudata->recon;
    if(!curecon) return;
    curcellfree(curecon->neai);
    curcellfree(curecon->opdwfs);
    curcellfree(curecon->grad); 
    curcellfree(curecon->opdr); 
    curcellfree(curecon->opdr_vec); 
    delete curecon->fdpcg;
    if(curecon->dmfit_vec!=curecon->dmfit){
	curcellfree(curecon->dmfit_vec);
    }else{
	curecon->dmfit_vec=NULL;
    }
    curcellfree(curecon->dmfit);
    free(curecon->l2c);
    free(curecon->zzi);
    free(curecon->zzv);
    if(curecon->W01){
	W01_T *W01=curecon->W01;
	curfree(W01->W1);
	delete W01->W0p;
	cudaFree(W01->W0f);
    }
    curcellfree(curecon->opdfit);
    curcellfree(curecon->opdfit2);
    curfree(curecon->pis);
    cudaFree(curecon->floc);
    curcellfree(curecon->fitNW);
    delete curecon->actslave;
    delete curecon->RCl;
    cudaFree(curecon->RCp);
    curfree(curecon->RUp);
    curfree(curecon->RVp);
    curfree(curecon->RMI);
    delete curecon->FCl;
    cudaFree(curecon->FCp);
    curfree(curecon->FUp);
    curfree(curecon->FVp);
    curfree(curecon->FMI);
    curcellfree(curecon->GXL);
    //    delete curecon;
}
void gpu_recon_free(){
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	if(cudata->recon){
	    gpu_recon_free_do();
	}
    }
}
void gpu_setup_recon_mvm(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    /*The following routine assemble MVM and put in recon->MVM*/
    if(parms->recon.mvm==1){
	gpu_setup_recon_mvm_trans(parms, recon, powfs);
    }else{
	gpu_setup_recon_mvm_direct(parms, recon, powfs);
    }
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	gpu_recon_free_do();
	CUDA_SYNC_DEVICE;
    }///for GPU
    if(!parms->sim.mvmport){
	gpu_set(gpu_recon);
	curecon_t *curecon=cudata->recon;
	cp2gpu(&curecon->MVM, recon->MVM);
    }
    gpu_print_mem("MVM");
}
void gpu_setup_recon_predict_do(const PARMS_T *parms, RECON_T *recon){
    if(!parms->gpu.tomo || !parms->tomo.predict){
	return;
    }
    curecon_t *curecon=cudata->recon;
    const int nwfs=parms->nwfsr;
    GPU_PROP_GRID_T *hxdata=new GPU_PROP_GRID_T[nwfs*recon->npsr];
    for(int ips=0; ips<recon->npsr; ips++){ 
	const float ht=recon->ht->p[ips]; 
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    const float hs = parms->powfs[ipowfs].hs; 
	    const float scale = 1.f - ht/hs; 
	    float dispx=parms->wfsr[iwfs].thetax*ht; 
	    float dispy=parms->wfsr[iwfs].thetay*ht; 
	    if(parms->tomo.predict){ 
		int ips0=parms->atmr.indps[ips]; 
		dispx+=cudata->atm[ips0].vx*parms->sim.dt*2; 
		dispy+=cudata->atm[ips0].vy*parms->sim.dt*2; 
	    } 
	    cugrid_t pmapscale=curecon->pmap*scale;
	    gpu_prop_grid_prep(hxdata+iwfs+ips*nwfs,  pmapscale, curecon->xmap[ips],
			       dispx, dispy, NULL); 
	    {
		float tmp=laplacian_coef(recon->r0, recon->wt->p[ips], recon->xmap[ips]->dx)*0.25f;
		hxdata[iwfs+ips*nwfs].l2c=tmp*tmp*TOMOSCALE;
		if(parms->tomo.piston_cr){
		    hxdata[iwfs+ips*nwfs].zzi=loccenter(recon->xloc[ips]);
		    hxdata[iwfs+ips*nwfs].zzv=tmp*tmp*TOMOSCALE*1e-6;
		}else{
		    hxdata[iwfs+ips*nwfs].zzi=-1;
		}
	    }
	}
    }
    if(!curecon->hxdata){
	DO(cudaMalloc(&curecon->hxdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr));
    }
    DO(cudaMemcpy(curecon->hxdata, hxdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr, cudaMemcpyHostToDevice));
    delete [] hxdata;
}
void gpu_setup_recon_predict(const PARMS_T *parms, RECON_T *recon){
    if(parms->recon.mvm && parms->gpu.tomo && parms->gpu.fit && !parms->load.mvm){
	for(int igpu=0; igpu<NGPU; igpu++){
	    gpu_set(igpu);
	    gpu_setup_recon_predict_do(parms, recon);
	}
    }else{
	gpu_set(gpu_recon);
	gpu_setup_recon_predict_do(parms, recon);
    }
}
/*update reconstruction parameters after slodar.*/
void gpu_update_recon(const PARMS_T *parms, RECON_T *recon){
    gpu_set(gpu_recon);
    TO_IMPLEMENT;//copy to GPU struct (hxtdata)
    curecon_t *curecon=cudata->recon;
    for(int ips=0; ips<recon->npsr; ips++){
	float tmp=laplacian_coef(recon->r0, recon->wt->p[ips], recon->xmap[ips]->dx)*0.25f;
	curecon->l2c[ips]=tmp*tmp*TOMOSCALE;
    }
    if(parms->tomo.piston_cr){
	for(int ips=0; ips<recon->npsr; ips++){
	    double r0=recon->r0;
	    double dx=recon->xloc[ips]->dx;
	    double wt=recon->wt->p[ips];
	    int icenter=loccenter(recon->xloc[ips]);
	    curecon->zzi[ips]=icenter;
	    curecon->zzv[ips]=pow(laplacian_coef(r0,wt,dx),2)*TOMOSCALE*1e-6;
	}
    }
}
void gpu_recon_reset(const PARMS_T *parms){/*reset warm restart.*/
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    curcellzero(curecon->opdr, 0);
    curcellzero(curecon->dmfit, 0);
    if(curecon->dm_wfs){
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    curcellzero(curecon->dm_wfs[iwfs], 0);
	}
    }
    if(curecon->dm_evl){
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    curcellzero(curecon->dm_evl[ievl], 0);
	}
    }
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	if(cudata->dm_wfs){
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		cudata->dm_wfs[iwfs].p.zero();
	    }
	}
	if(cudata->dm_evl){
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		cudata->dm_evl[ievl].p.zero();
	    }
	}
	CUDA_SYNC_DEVICE;
    }
}

void gpu_tomo(SIM_T *simu){
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    TIC_test;tic_test;
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    curecon->reconisim=simu->reconisim;
#if 0
    gpu_tomo_test(simu);
#endif
    toc_test("Before gradin");
    cp2gpu(&curecon->gradin, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
    toc_test("Gradin");
    cudaProfilerStart();
    curcell *opdrsave=NULL;
    warning_once("remove opdrsave after debugging\n");
    curcellcp(&opdrsave, curecon->opdr, curecon->cgstream);
    simu->cgres->p[0]->p[simu->reconisim]=
	gpu_tomo_do(parms, recon, curecon->opdr, NULL, curecon->gradin, curecon->cgstream);
    //Sanity check the result
    float opdrmax=curcellmax(curecon->opdr, curecon->cgstream);
    if(curecon->opdr->nx>1 && opdrmax>4e-6){
	simu->status->warning=2;
	info("opdrmax=%g\n", opdrmax);
	curcellwrite(curecon->gradin, "dbg_gradin_%d", simu->reconisim);
	curcellwrite(opdrsave, "dbg_opdrlast_%d", simu->reconisim);
	curcellwrite(curecon->opdr, "dbg_opdr_%d", simu->reconisim);
	curcellcp(&curecon->opdr, opdrsave, curecon->cgstream);
	extern int pcg_save;
	pcg_save=1;
	double newres=gpu_tomo_do(parms, recon, curecon->opdr, NULL, curecon->gradin, curecon->cgstream);
	pcg_save=0;
	curcellwrite(curecon->opdr, "dbg_opdrredo_%d", simu->reconisim);
	info("oldres=%g. newres=%g\n", simu->cgres->p[0]->p[simu->reconisim], newres);
    }
    curfree(opdrsave);
    if(!parms->gpu.fit || parms->save.opdr || (recon->moao && !parms->gpu.moao)){
	cp2cpu(&simu->opdr, 0, curecon->opdr_vec, 1, curecon->cgstream);
    }
    if(parms->recon.split==2){
	curcell *gngsmvst=NULL;
	curcellmm(&gngsmvst, 1, curecon->GXL, curecon->opdr_vec, "nn", 1./parms->sim.dtrat_lo, curecon->cgstream);
	add2cpu(&simu->gngsmvst, gngsmvst, curecon->cgstream);
	curfree(gngsmvst);
    }
    if(parms->dbg.deltafocus){
	curcell *tmp1=NULL;
	curcellmm(&tmp1, 1, curecon->RFdfx, curecon->opdr_vec, "nn", 1, curecon->cgstream);
	scell *tmp=NULL;
	cp2cpu(&tmp, tmp1, curecon->cgstream);
	curcellfree(tmp1);
	if(tmp->nx!=1 || tmp->ny!=1 || tmp->p[0]->nx!=1 || tmp->p[0]->ny!=1){
	    error("Wrong format");
	}
	curecon->cgstream.sync();
	simu->deltafocus=tmp->p[0]->p[0];
	scellfree(tmp); 
    }else{
	curecon->cgstream.sync();
    }
    cudaProfilerStop();
    toc_test("Tomo");
}

void gpu_fit(SIM_T *simu){
    TIC_test;tic_test;
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(!parms->gpu.tomo){
	cp2gpu(&curecon->opdr_vec, simu->opdr);
    }
#if 0
    gpu_fit_test(simu);
#endif
    toc_test("Before FitR");
    curcell *fitsave=NULL;//for debugging purpose.
    warning_once("remove fitsave after debugging\n");
    curcellcp(&fitsave, curecon->dmfit, curecon->cgstream);
    simu->cgres->p[1]->p[simu->reconisim]=
	gpu_fit_do(parms, recon, curecon->dmfit, curecon->opdr, curecon->cgstream);
    curecon->cgstream.sync();
    if(simu->reconisim>0 && simu->cgres->p[1]->p[simu->reconisim]>simu->cgres->p[1]->p[simu->reconisim-1]*2){
	simu->status->warning=2;
	curcellwrite(curecon->gradin, "dbg_gradin_%d", simu->reconisim);
	curcellwrite(curecon->dmfit, "dbg_dmfit_%d", simu->reconisim);
	curcellwrite(curecon->opdr, "dbg_opdr_%d", simu->reconisim);
	curcellwrite(fitsave, "dbg_dmfitlast_%d", simu->reconisim);
	curcellcp(&curecon->dmfit, fitsave, curecon->cgstream);
	double newres=gpu_fit_do(parms, recon, curecon->dmfit, curecon->opdr, curecon->cgstream);
	curcellwrite(curecon->dmfit, "dbg_dmfitredo_%d", simu->reconisim);
	info("oldres=%g newres=%g\n", simu->cgres->p[1]->p[simu->reconisim], newres);
    }
    cp2cpu(&simu->dmfit, 0, curecon->dmfit_vec, 1, curecon->cgstream);
    curfree(fitsave);
    /*Don't free opdr. Needed for warm restart in tomo.*/
    toc_test("Fit");
}
void gpu_recon_mvm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    cp2gpu(&curecon->gradin, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
    curcellmm(&curecon->dmfit_vec, 0., curecon->MVM, curecon->gradin,"nn", 1., curecon->cgstream);
    cp2cpu(&simu->dmerr, 0., curecon->dmfit_vec, 1., curecon->cgstream);
    curecon->cgstream.sync();
}
