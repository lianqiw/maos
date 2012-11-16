/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
				    int nsa, float ox, float oy, float dx){
    const int step=blockDim.x * gridDim.x;
    const float dx1=1./dx;
    for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
	saptr[isa][0]=(int)roundf((saloc[isa][0]-ox)*dx1);
	saptr[isa][1]=(int)roundf((saloc[isa][1]-oy)*dx1);
    }
}
W01_T *gpu_get_W01(dsp *R_W0, dmat *R_W1){
    if(!R_W0 || !R_W1){
	error("R0, R1 must not be empty\n");
    }
    W01_T *W01=(W01_T*)calloc(1, sizeof(W01_T));
    cp2gpu(&W01->W1, R_W1);
    if(0){warning("change to 0\n");
	W01->nW0f=0;
	cp2gpu(&W01->W0p, R_W0);
    }else{
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
	dsp *W0new2=sptrans(W0new);
	cp2gpu(&W01->W0p, W0new2);
	cp2gpu(&W01->W0f, full, count2);
	W01->nW0f=count2;
	spfree(W0new);
	spfree(W0new2);
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
    
    curecon->cgstream =new stream_t;
    if((parms->gpu.tomo || parms->gpu.fit) && !parms->sim.idealfit){
	curecon->opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
	curecon->opdr_vec=curcellnew(recon->npsr, 1);
	for(int ips=0; ips<recon->npsr; ips++){
	    curecon->opdr_vec->p[ips]=curref(curecon->opdr->p[ips]);
	    curecon->opdr_vec->p[ips]->nx=curecon->opdr->p[ips]->nx*curecon->opdr->p[ips]->ny;
	    curecon->opdr_vec->p[ips]->ny=1;
	}
    }
    if(parms->gpu.tomo){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].skip) continue;
	    int nsa=powfs[ipowfs].pts->nsa;
	    cudaMalloc(&cupowfs[ipowfs].saptr, nsa*2*sizeof(int));
	    saloc2ptr_do<<<DIM(nsa,256)>>>
		(cupowfs[ipowfs].saptr, cupowfs[ipowfs].saloc, nsa, 
		 recon->pmap->ox, recon->pmap->oy, recon->pmap->dx);
	    if(recon->GP->p[ipowfs]){
		const int use_mat=parms->tomo.pos==2 ||parms->tomo.pos==1 ;
		if(use_mat){//normally true
		    dsp *GP=sptrans(recon->GP->p[ipowfs]);
		    spint *pp=GP->p;
		    spint *pi=GP->i;
		    double *px=GP->x;
		    dmat *partxy=NULL;
		    int np1=parms->tomo.pos+1;
		    int np=np1*np1;
		    int zmax=parms->tomo.pos;
		    partxy=dnew(np*2, nsa);
		    int nsa=powfs[ipowfs].pts->nsa;
		    double dx1=1./recon->ploc->dx;
		
		    for(int ic=0; ic<GP->n; ic++){
			int isa=(ic<nsa)?ic:(ic-nsa);
			for(spint ir=pp[ic]; ir<pp[ic+1]; ir++){
			    int ix=pi[ir];
			    double lx=recon->ploc->locx[ix];
			    double ly=recon->ploc->locy[ix];
			    double sx=powfs[ipowfs].saloc->locx[isa];
			    double sy=powfs[ipowfs].saloc->locy[isa];
			    int zx=(int)round((lx-sx)*dx1);
			    int zy=(int)round((ly-sy)*dx1);
			    /**
			       Squeeze the weights to closed points in 3x3 in the subaperture.
			       Does not work well. simply drop these points doesn't work either.
			    */
			    if(zx<0 || zx>zmax || zy<0 || zy>zmax){
				warning("isa=%d, zxy=%d %d\n", isa, zx, zy);
			    }
			    if(zx<0) zx=0;
			    if(zx>zmax) zx=zmax;
			    if(zy<0) zy=0;
			    if(zy>zmax) zy=zmax;
			    partxy->p[np*2*isa+zx+zy*np1+(ic<nsa?0:np)]+=px[ir];
			}
		    }
		    cp2gpu(&cupowfs[ipowfs].GPp, partxy);
		    dfree(partxy);
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
	    int iwfs0=parms->powfs[ipowfs].wfs[0];/*first wfs in this group. */
	    if(iwfs!=iwfs0 && recon->saneai->p[iwfs+iwfs*parms->nwfsr]->p
	       ==recon->saneai->p[iwfs0+iwfs0*parms->nwfsr]->p){
		curecon->neai->p[iwfs]=curref(curecon->neai->p[iwfs0]);
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
	if(recon->PTT && !curecon->PTT){
	    cp2gpu(&curecon->PTT, recon->PTT);
	}
	if(parms->tomo.precond==1){/*fdpcg*/
	    FDPCG_T *fdpcg=recon->fdpcg;
	    cufdpcg_t *cufd=curecon->fdpcg=new cufdpcg_t;
	    cufd->scale=fdpcg->scale;
	    cufd->half=fdpcg->half;
	    int nb=fdpcg->Mbinv->nx;
	    int bs=fdpcg->Mbinv->p[0]->nx;
	    cp2gpu(&cufd->perm, fdpcg->perm, nb*bs);//not bs*bs
	    cp2gpu(&cufd->Mb, fdpcg->Mbinv);
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
		DO(cufftSetStream(cufd->ffti[ic], curecon->cgstream[0]));

		DO(cufftSetStream(cufd->fft[ic], curecon->cgstream[0]));
	    }
	    cufd->xhat1=cuccellnew(recon->npsr, 1, recon->xnx, recon->xny);
	    {
		int nby=256/bs;
		int nbz=nb/nby;
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

	const float oxp=recon->pmap->ox;
	const float oyp=recon->pmap->oy;
	const float dxp=recon->pmap->dx;
	GPU_PROP_GRID_T *hxdata=new GPU_PROP_GRID_T[nwfs*recon->npsr];
	GPU_PROP_GRID_T *hxtdata=new GPU_PROP_GRID_T[nwfs*recon->npsr];
	for(int ips=0; ips<recon->npsr; ips++){ 
	    const float ht=recon->ht->p[ips]; 
	    const float oxx=recon->xmap[ips]->ox; 
	    const float oyx=recon->xmap[ips]->oy; 
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		const int ipowfs = parms->wfsr[iwfs].powfs;
		if(parms->powfs[ipowfs].skip) continue;
		const float hs = parms->powfs[ipowfs].hs; 
		const float scale = 1.f - ht/hs; 
		float dispx=parms->wfsr[iwfs].thetax*ht; 
		float dispy=parms->wfsr[iwfs].thetay*ht; 
		gpu_prop_grid_prep(hxdata+iwfs+ips*nwfs, curecon->opdwfs->p[iwfs], oxp*scale, oyp*scale, dxp*scale, 
				   curecon->opdr->p[ips], oxx, oyx, recon->xmap[ips]->dx, 
				   dispx, dispy, 'n'); 
		gpu_prop_grid_prep(hxtdata+iwfs+ips*nwfs, curecon->opdwfs->p[iwfs], oxp*scale, oyp*scale, dxp*scale, 
				   curecon->opdr->p[ips], oxx, oyx, recon->xmap[ips]->dx, 
				   dispx, dispy, 't'); 
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
	DO(cudaMalloc(&curecon->hxdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr));
	DO(cudaMemcpy(curecon->hxdata, hxdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr, cudaMemcpyHostToDevice));
	DO(cudaMalloc(&curecon->hxtdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr));
	DO(cudaMemcpy(curecon->hxtdata, hxtdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr, cudaMemcpyHostToDevice));
	delete [] hxdata;
	delete [] hxtdata;
	GPU_GP_T *gpdata=new GPU_GP_T[nwfs];
	if(recon->PDF){
	    curecon->PDF=(curcell**)calloc(sizeof(curcell*), nwfs);
	}
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    gpdata[iwfs].saptr=cupowfs[ipowfs].saptr;
	    gpdata[iwfs].dsa=powfs[ipowfs].pts->dsa;
	    gpdata[iwfs].GPp=cupowfs[ipowfs].GPp->p;
	    gpdata[iwfs].pos=parms->tomo.pos;
	    if(curecon->PTT){
		gpdata[iwfs].PTT=curecon->PTT->p[iwfs+iwfs*nwfs]->p;
	    }
	    if(recon->PDF){
		dcell *tmp=dcellnew(recon->PDF->nx, 1);
		for(int i=0; i<tmp->nx; i++){
		    tmp->p[i]=dref(recon->PDF->p[i+nwfs*iwfs]);
		}
		cp2gpu(&curecon->PDF[iwfs], tmp);
		dcellfree(tmp);
		gpdata[iwfs].PDF=curecon->PDF[iwfs]->pm;
	    }
	    gpdata[iwfs].neai=(const float(*)[3])curecon->neai->p[iwfs]->p;
	    gpdata[iwfs].nsa=powfs[ipowfs].pts->nsa;
	    gpdata[iwfs].nxp=recon->pmap->nx;
	    gpdata[iwfs].dxp=recon->pmap->dx;
	    gpdata[iwfs].oxp=recon->pmap->ox;
	    gpdata[iwfs].oyp=recon->pmap->oy;
	}
	DO(cudaMalloc(&curecon->gpdata, sizeof(GPU_GP_T)*nwfs));
	DO(cudaMemcpy(curecon->gpdata, gpdata, sizeof(GPU_GP_T)*nwfs, cudaMemcpyHostToDevice));
	delete [] gpdata;
    }
    if(parms->gpu.fit){
	curecon->fitstream=new stream_t[parms->fit.nfit];
	curecon->psstream =new stream_t[recon->npsr];
	curecon->dmstream =new stream_t[recon->ndm];
	if(parms->gpu.fit==1){ /*For fitting using sparse matrix*/
	    cp2gpu(&curecon->FR, &recon->FR);
	    cp2gpu(&curecon->FL, &recon->FL);
	    curecon->FR.fitstream=curecon->fitstream;
	    curecon->FR.dmstream=curecon->dmstream;
	    curecon->FL.fitstream=curecon->fitstream;
	    curecon->FL.dmstream=curecon->dmstream;
	    curecon->dmfit=curcellnew(parms->ndm, 1, recon->anloc, (long*)NULL);
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
	    cp2gpu(&curecon->fitNW, recon->fitNW);
	    cp2gpu(&curecon->actslave, recon->actslave);
	    curecon->dmfit=curcellnew(parms->ndm, 1);
	    curecon->dmfit_vec=curcellnew(parms->ndm, 1);
	    int ntotact=0;
	    for(int idm=0; idm<parms->ndm; idm++){
		ntotact+=recon->amap[idm]->nx*recon->amap[idm]->ny;
	    }
	    curecon->dmfit->m=curnew(ntotact, 1);
	    int ct=0;
	    for(int idm=0; idm<parms->ndm; idm++){
		curecon->dmfit->p[idm]=new cumap_t(recon->amap[idm]->nx,recon->amap[idm]->ny,
						   curecon->dmfit->m->p+ct, 0);
		ct+=recon->amap[idm]->nx*recon->amap[idm]->ny;
		curecon->dmfit_vec->p[idm]=curref(curecon->dmfit->p[idm]);
		curecon->dmfit_vec->p[idm]->nx=curecon->dmfit_vec->p[idm]->nx
		    *curecon->dmfit_vec->p[idm]->ny;
		curecon->dmfit_vec->p[idm]->ny=1;
	    }
	}
	curecon->cubic_cc=new float *[parms->ndm];
	for(int idm=0; idm<parms->ndm; idm++){
	    curecon->cubic_cc[idm]=gpu_dmcubic_cc(parms->dm[idm].iac);
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
	curecon->opdfit=curcellnew(nfit, 1);
	curecon->opdfit2=curcellnew(nfit, 1);
	for(int ifit=0; ifit<nfit; ifit++){
	    curecon->opdfit->p[ifit] =curnew(recon->fmap->nx, recon->fmap->ny);
	    curecon->opdfit2->p[ifit]=curnew(recon->fmap->nx, recon->fmap->ny);
	}
	curecon->pis=curnew(parms->fit.nfit, 1);
    }
 
    if(recon->RFlgsx){
	cp2gpu(&curecon->RFlgsx, recon->RFlgsx);
    }
    if(recon->RFngsx){
	cp2gpu(&curecon->RFngsx, recon->RFngsx);
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
static void gpu_recon_free_do(){
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
    //    delete curecon;
}
void gpu_recon_free(){
    gpu_set(gpu_recon);
    gpu_recon_free_do();
}

typedef struct MVM_IGPU_T{
    const PARMS_T *parms;
    RECON_T *recon;
    POWFS_T *powfs;
    curcell *mvmig; /*intermediate TomoL result*/
    curcell *mvmfg; /*intermediate FitR result*/
    smat *mvmt;     /*result: tranpose of MVM calculated by this GPU.*/
    float *FLI;
    smat *residual;
    long (*curp)[2];
    int ntotact;
    int ntotgrad;
    int load_mvmf; /*intermediate FitR result is for 1) loading, 0) saving.*/
}MVM_IGPU_T;
void gpu_setup_recon_mvm_igpu(thread_t *info){
    TIC;tic;
    double tk_prep=0, tk_fitL=0, tk_fitR=0, tk_TomoL=0, tk_TomoR=0, tk_cp=0;
    MVM_IGPU_T *data=(MVM_IGPU_T*)info->data;
    const PARMS_T *parms=data->parms;
    RECON_T *recon=data->recon;
    smat *residual=data->residual;
    long (*curp)[2]=data->curp;
    const int ntotact=data->ntotact;
    const int ntotgrad=data->ntotgrad;
    const int load_mvmf=data->load_mvmf;
    int igpu=info->ithread;
    gpu_set(igpu);
    curecon_t *curecon=cudata->recon;
    curmat *mvmi=data->mvmig?data->mvmig->p[igpu]:NULL;/*Tomography output, for warm restart*/
    curmat *mvmf=data->mvmfg?data->mvmfg->p[igpu]:NULL;/*loaded FitR output.*/
    /*Tomography*/
    G_PREFUN prefun=NULL;
    void *predata=NULL;
    if(parms->tomo.precond==1){
	prefun=gpu_Tomo_fdprecond;
	predata=(void*)recon;
    }
    G_CGFUN cg_fun;
    void *cg_data;
    curcell *eyec=NULL;/* Only use eyec for CG.*/
    float eye2[2]={0,1.};
    float eye1[1]={1.};
    //const int nwfs=parms->nwfs;
    const int ndm=parms->ndm;
    if(parms->gpu.fit==1){//sparse matrix
	cg_fun=(G_CGFUN) cumuv;
	cg_data=&curecon->FL;
    }else{
	cg_fun=(G_CGFUN) gpu_FitL;
	cg_data=recon;
    }
    const float *FLI=data->FLI;
    if(!FLI && !load_mvmf){
	if(parms->fit.square){
	    eyec=curcellnew(ndm, 1, recon->anx, recon->any);
	}else{
	    eyec=curcellnew(ndm, 1, recon->anloc, (long*)0);
	}
    }
 
    curcell *dmfit=load_mvmf?NULL:curcellnew(curecon->dmfit);
    curcell *opdx=curcellnew(recon->npsr, 1, recon->xnx, recon->xny, (float*)(mvmf?1:0));
    curcell *opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny, (float*)(mvmi?1:0));
    curcell *grad=curcellnew(parms->nwfs, 1, recon->ngrad, (long*)0, (float*)1);
    if(ntotact==0){
	error("ntotact=0;\n");
    }
    curmat *mvmt=curnew(ntotgrad, info->end-info->start);/*contains result*/
    tk_prep+=toc3;tic;
    stream_t &stream=curecon->cgstream[0];
    for(int iact=info->start; iact<info->end; iact++){
	int curdm=curp[iact][0];
	int curact=curp[iact][1];
	if(info->ithread==0){
	    if(!detached){
		info2("%6d of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", iact*NGPU, ntotact);
	    }else if(iact % 100==0){
		info2("%6d of %6d\n", iact*NGPU, ntotact);
	    }
	}
	if(eyec){
	    if(iact){
		cudaMemcpyAsync(eyec->m->p+iact-1, eye2, 2*sizeof(float),
				cudaMemcpyHostToDevice, stream);
	    }else{
		cudaMemcpyAsync(eyec->m->p+iact, eye1, sizeof(float), 
				cudaMemcpyHostToDevice, stream);
	    }
	}
	if(!recon->actcpl || recon->actcpl->p[curdm]->p[curact]>EPS){
	    if(mvmf) opdx->replace(mvmf->p+(iact-info->start)*mvmf->nx, 0, stream);
	    if(!load_mvmf){
		if(eyec){
		    /*Fitting operator*/
		    curcellzero(dmfit, stream);//temp
		    if(gpu_pcg(&dmfit, (G_CGFUN)cg_fun, cg_data, NULL, NULL, eyec, &curecon->cgtmp_fit,
			       parms->recon.warm_restart, parms->fit.maxit, stream)>1.){
			warning("Fit CG not converge.\n");
		    }
		}else{
		    cudaMemcpyAsync(dmfit->m->p, FLI+iact*ntotact, sizeof(float)*ntotact, 
				    cudaMemcpyHostToDevice, stream);
		}
		cudaStreamSynchronize(stream);
    		tk_fitL+=toc3; tic;
		/*Transpose of fitting operator*/
		if(parms->gpu.fit==1){//sparse matrix
		    cumuv_trans(&opdx, 0, &curecon->FR, dmfit, 1);
		}else{
		    gpu_FitRt(&opdx, 0, recon, dmfit, 1);
		}
	    }
	    tk_fitR+=toc3; tic;
	    switch(parms->tomo.alg){
	    case 0:
		if(!opdr->m || !opdx->m){
		    error("opdr and opdx must be continuous\n");
		}
		cuchol_solve(opdr->m->p, curecon->RCl, curecon->RCp, opdx->m->p, stream);
		if(curecon->RUp){
		    curmat *tmp=curnew(curecon->RVp->ny, 1);
		    curmv(tmp->p, 0, curecon->RVp, opdx->m->p, 't', -1, stream);
		    curmv(opdr->m->p, 1, curecon->RUp, tmp->p, 'n', 1, stream);
		    curfree(tmp);
		}
		break;
	    case 1:{
		if(mvmi){
		    opdr->replace(mvmi->p+(iact-info->start)*mvmi->nx, 0, stream);
		}
		if(parms->recon.mvm==2){//disable warm restart in CG using neighboring act.
		    curcellzero(opdr, stream); 
		}
		int disablelrt=curecon->disablelrt;
		curecon->disablelrt=1;
		/*disable the t/t removal lrt in split tomo that creats problem in fdpcg mode*/
		if((residual->p[iact]=gpu_pcg(&opdr, gpu_TomoL, recon, prefun, predata, opdx, &curecon->cgtmp_tomo,
					      parms->recon.warm_restart, parms->tomo.maxit,
					      stream, parms->tomo.cgthres))>1){
		    warning2("Tomo CG residual is %.2f for %d\n", residual->p[iact], iact);
		}
		curecon->disablelrt=disablelrt;
	    }
		break;
	    case 2:
		curmv(opdr->m->p, 0, curecon->RMI, opdx->m->p, 'n', 1, stream);
		break;
	    default:
		error("Invalid");
	    }
	    tk_TomoL+=toc3; tic;
	    /*Right hand side. output directly to mvmt*/
	    grad->replace(mvmt->p+(iact-info->start)*ntotgrad, 0, stream);
	    gpu_TomoRt(&grad, 0, recon, opdr, 1, stream);
	    tk_TomoR+=toc3; tic;
	}
    }//for iact
    int nn=ntotgrad*(info->end-info->start)*sizeof(float);
    float *mvmtc=data->mvmt->p+info->start*ntotgrad;
    cudaMemcpyAsync(mvmtc, mvmt->p, nn, cudaMemcpyDeviceToHost, curecon->cgstream[0]);
    cudaStreamSynchronize(curecon->cgstream[0]);
    curcellfree(dmfit);
    curcellfree(opdx);
    curcellfree(opdr);
    curcellfree(grad);
    curcellfree(eyec);
    curfree(mvmt);
    tk_cp+=toc3;tic;
    info2("GPU %d: Prep %.2f FitL %.2f FitR %.2f TomoL %.1f TomoR %.1f cp %.2f\n", 
	  igpu, tk_prep, tk_fitL, tk_fitR, tk_TomoL, tk_TomoR, tk_cp);
}
void gpu_setup_recon_mvm(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    TIC;tic;
    if(parms->recon.alg!=0){
	error("Please adept to LSR\n");
    } 
    if(!parms->load.mvm){
	info2("Assembling MVR MVM in GPU\n");
	int ntotact=0;
	int ntotgrad=0;
	int ntotxloc=0;
	const int ndm=parms->ndm;
	for(int idm=0; idm<ndm; idm++){
	    ntotact+=recon->anloc[idm];
	} 
	for(int ips=0; ips<recon->npsr; ips++){
	    ntotxloc+=recon->xloc[ips]->nloc;
	}
	long (*curp)[2]=(long(*)[2])malloc(ntotact*2*sizeof(long));
	int nact=0;
	for(int idm=0; idm<ndm; idm++){
	    for(int iact=0; iact<recon->anloc[idm]; iact++){
		curp[nact+iact][0]=idm;
		curp[nact+iact][1]=iact;
	    }
	    nact+=recon->anloc[idm];
	}   
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    ntotgrad+=recon->ngrad[iwfs];
	}
	
	smat *residual=NULL;
	if(parms->tomo.alg==1){
	    residual=snew(ntotact, 1);
	}
	dmat *FLId=NULL; /* MI is inv(FL) for direct methods*/
	float *FLI=NULL;

	/* Loading or saving intermediate TomoL result. */
	smat *mvmi=NULL; 
	if(parms->load.mvmi){
	    mvmi=sread("%s", parms->load.mvmi);
	    if(mvmi->nx!=ntotxloc || mvmi->ny!=ntotact){
		error("loaded mvmi has dimension (%ld, %ld) but we expect (%d, %d)",
		      mvmi->nx, mvmi->ny, ntotxloc, ntotact);
	    }
       	}else if(parms->save.mvmi){
	    mvmi=snew(ntotxloc, ntotact);
	}
	curcell *mvmig=NULL;
	if(mvmi){
	    mvmig=curcellnew(NGPU, 1);
	}

	/* Loading or saving intermediate FitR Result.*/
	smat *mvmf=NULL;
	if(parms->load.mvmf){
	    /*Load FitR FitL results from file. Resembling warm restart case
	      where mvmf is kept in memory*/
	    mvmf=sread("%s", parms->load.mvmf);
	    if(mvmf->nx!=ntotxloc || mvmf->ny!=ntotact){
		error("loaded mvmf has dimension (%ld, %ld) but we expect (%d, %d)",
		      mvmf->nx, mvmf->ny, ntotxloc, ntotact);
	    }
	}else if(parms->save.mvmf){
	    /*save FitR FitL resutls to file, for later loading.*/
	    mvmf=snew(ntotxloc, ntotact);
	}
	curcell *mvmfg=NULL;
	if(mvmf){
	    mvmfg=curcellnew(NGPU, 1);
	}
	if(!parms->load.mvmf){
	    /*Prepare FitR, FitL is don't load fitting results using mvmf*/
	    switch(parms->fit.alg){
	    case 0:{
		dmat *eye=dnew(ntotact, ntotact);
		daddI(eye, 1);
		FLId=dnew(ntotact, ntotact);
		muv_direct_solve(&FLId, &recon->FL, eye);
		dfree(eye);
		toc("Fit CBS");tic;
	    }
		break;
	    case 1:
		break;
	    case 2:
		FLId=dref(recon->FL.MI);
		break;
	    default:
		error("Invalid fit.alg=%d\n", parms->fit.alg);
	    }
	    if(FLId){
		FLI=(float*)malloc4async(sizeof(float)*ntotact*ntotact);
		for(long i=0; i<ntotact*ntotact; i++){
		    FLI[i]=(float)FLId->p[i];
		}
		dwrite(FLId, "FLId");
		dfree(FLId);
	    }
	}
    	smat *mvmt=snew(ntotgrad, ntotact);
	MVM_IGPU_T data={parms, recon, powfs, mvmig, mvmfg, mvmt, FLI, residual, curp, ntotact, ntotgrad, parms->load.mvmf?1:0};
	int nthread=NGPU;
	thread_t info[NGPU];
	thread_prep(info, 0, ntotact, nthread, gpu_setup_recon_mvm_igpu, &data);

	/*Initialyze intermediate TomoL result array in GPU. Send intermediate
	  TomoL results to GPU if load.mvmi is set.*/
	if(mvmi){
	    TIC;tic;
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		mvmig->p[i]=curnew(ntotxloc, info[i].end-info[i].start);
		if(parms->load.mvmi){
		    cudaMemcpy(mvmig->p[i]->p, mvmi->p+info[i].start*ntotxloc, 
			       sizeof(float)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyHostToDevice);
		}
	    }
	    if(parms->load.mvmi){
		toc2("copy mvmi to gpu");
	    }
	}
	/*Initialyze intermediate FitL/FitR result array in GPU. Send
	  intermediate FitL/FitR results to GPU if load.mvmf is set.*/
	if(mvmf){
	    TIC;tic;
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		mvmfg->p[i]=curnew(ntotxloc, info[i].end-info[i].start);
		if(parms->load.mvmf){
		    cudaMemcpy(mvmfg->p[i]->p, mvmf->p+info[i].start*ntotxloc, 
			       sizeof(float)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyHostToDevice);
		}
	    }
	    if(parms->load.mvmf){
		toc2("copy mvmf to gpu");
	    }
	}
	/*Do real MVM control matrix assemble in multiply CPU/GPU*/
	CALL_THREAD(info, nthread, 1);
	/*Copy MVM control matrix results back*/
	{
	    TIC;tic;
	    int ndm=parms->ndm;
	    int nwfs=parms->nwfs;
	    recon->MVM=dcellnew(ndm, nwfs);
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(!parms->powfs[ipowfs].skip){
		    for(int idm=0; idm<ndm; idm++){
			recon->MVM->p[idm+ndm*iwfs]=dnew(recon->anloc[idm], powfs[ipowfs].saloc->nloc*2);
		    }
		}
	    }
	    dmat *mvmtt=dnew(mvmt->ny, mvmt->nx);
	    for(int iy=0; iy<mvmtt->ny; iy++){
		for(int ix=0; ix<mvmtt->nx; ix++){
		    mvmtt->p[ix+iy*mvmtt->nx]=(double)mvmt->p[iy+ix*mvmt->nx];
		}
	    }
	    toc2("MVM Reshape in CPU 1");
	    sfree(mvmt);
	    d2cell(&recon->MVM, mvmtt, NULL);
	    dfree(mvmtt);
	    toc2("MVM Reshape in CPU 2");
	}
	if(parms->save.setup || parms->save.mvm){
	    dcellwrite(recon->MVM, "MVM.bin");
	}
	swrite(residual, "MVM_RL_residual");
	
	if(parms->save.mvmi){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		cudaMemcpy(mvmi->p+info[i].start*ntotxloc, mvmig->p[i]->p,  
			   sizeof(float)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyDeviceToHost);
	    }
	    swrite(mvmi, "MVM_Tomo.bin");
	}
	if(parms->save.mvmf){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		cudaMemcpy(mvmf->p+info[i].start*ntotxloc, mvmfg->p[i]->p,  
			   sizeof(float)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyDeviceToHost);
	    }
	    swrite(mvmf, "MVM_FitL.bin");
	}
	if(mvmig){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		curfree(mvmig->p[i]);
	    }
	    curcellfree(mvmig);
	}
	if(mvmfg){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		curfree(mvmfg->p[i]);
	    }
	    curcellfree(mvmfg);
	}
	sfree(mvmi);
	sfree(mvmf);
	sfree(residual);

	free(curp);
	if(FLI) free4async(FLI);
    }//if assemble in gpu
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	gpu_recon_free_do();
	CUDA_SYNC_DEVICE;
    }///for GPU
    if(!parms->sim.mvmport){
	gpu_set(gpu_recon);
	curecon_t *curecon=cudata->recon;
	cp2gpu(&curecon->MVM, recon->MVM);
	dcellfree(recon->MVM);
    }
    toc("MVM Final");
    gpu_print_mem("MVM");
}
void gpu_setup_recon_predict(const PARMS_T *parms, RECON_T *recon){
    if(!parms->gpu.tomo || !parms->tomo.predict){
	return;
    }
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    const int nwfs=parms->nwfs;
    const float oxp=recon->pmap->ox;
    const float oyp=recon->pmap->oy;
    const float dxp=recon->pmap->dx;
    GPU_PROP_GRID_T *hxdata=new GPU_PROP_GRID_T[nwfs*recon->npsr];
    GPU_PROP_GRID_T *hxtdata=new GPU_PROP_GRID_T[nwfs*recon->npsr];
    for(int ips=0; ips<recon->npsr; ips++){ 
	const float ht=recon->ht->p[ips]; 
	const float oxx=recon->xmap[ips]->ox; 
	const float oyx=recon->xmap[ips]->oy; 
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    const float hs = parms->powfs[ipowfs].hs; 
	    const float scale = 1.f - ht/hs; 
	    float dispx=parms->wfsr[iwfs].thetax*ht; 
	    float dispy=parms->wfsr[iwfs].thetay*ht; 
	    if(parms->tomo.predict){ 
		int ips0=parms->atmr.indps[ips]; 
		dispx+=cudata->atm[ips0]->vx*parms->sim.dt*2; 
		dispy+=cudata->atm[ips0]->vy*parms->sim.dt*2; 
	    } 
	    gpu_prop_grid_prep(hxdata+iwfs+ips*nwfs, curecon->opdwfs->p[iwfs], oxp*scale, oyp*scale, dxp*scale, 
			       curecon->opdr->p[ips], oxx, oyx, recon->xmap[ips]->dx, 
			       dispx, dispy, 'n'); 
	    gpu_prop_grid_prep(hxtdata+iwfs+ips*nwfs, curecon->opdwfs->p[iwfs], oxp*scale, oyp*scale, dxp*scale, 
			       curecon->opdr->p[ips], oxx, oyx, recon->xmap[ips]->dx, 
			       dispx, dispy, 't'); 
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
    if(!curecon->hxtdata){
	DO(cudaMalloc(&curecon->hxtdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr));
    }
    DO(cudaMemcpy(curecon->hxdata, hxdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr, cudaMemcpyHostToDevice));
    DO(cudaMemcpy(curecon->hxtdata, hxtdata, sizeof(GPU_PROP_GRID_T)*nwfs*recon->npsr, cudaMemcpyHostToDevice));
    delete [] hxdata;
    delete [] hxtdata;
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
	curcellzero(cudata->dm_wfs,0);
	curcellzero(cudata->dm_evl,0);
	CUDA_SYNC_DEVICE;
    }
}

void gpu_tomo(SIM_T *simu){
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    TIC_test;tic_test;
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    /*if(parms->tomo.pos!=2){
      TO_IMPLEMENT;
      }*/
    /*if(curecon->PDF){
	TO_IMPLEMENT;
	}*/
    /*first send gradients to GPU. can be skipped if keep grad in gpu. fast though. */
    /*Create temporary memory */
    curecon->reconisim=simu->reconisim;
#if 0
    gpu_tomo_test(simu);
#endif
    toc_test("Before gradin");
    cp2gpu(&curecon->gradin, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
    toc_test("Gradin");
    curcell *rhs=NULL;
    gpu_TomoR(&rhs, 0, recon, curecon->gradin, 1, curecon->cgstream[0]);
    toc_test("TomoR");
    switch(parms->tomo.alg){
    case 0:
	if(!curecon->opdr->m){
	    error("opdr must be continuous\n");
	}
	if(!rhs->m){
	    error("rhs must be continuous\n");
	}
	cuchol_solve(curecon->opdr->m->p, curecon->RCl, curecon->RCp, rhs->m->p, curecon->cgstream[0]);
	if(curecon->RUp){
	    curmat *tmp=curnew(curecon->RVp->ny, 1);
	    curmv(tmp->p, 0, curecon->RVp, rhs->m->p, 't', -1, curecon->cgstream[0]);
	    curmv(curecon->opdr->m->p, 1, curecon->RUp, tmp->p, 'n', 1, curecon->cgstream[0]);
	    cudaStreamSynchronize(curecon->cgstream[0]);
	    curfree(tmp);
	}
	/*{
	  curcellwrite(rhs, "GPU_RHS");
	  curcellwrite(curecon->opdr, "GPU_OPDR");
	  muv_solve(&simu->opdr, &recon->RL, &recon->RR, simu->gradlastol);
	  dcellwrite(simu->opdr, "CPU_OPDR");
	  exit(1);
	  }*/
	break;
    case 1:{
	G_PREFUN prefun=NULL;
	void *predata=NULL;
	if(parms->tomo.precond==1){
	    prefun=gpu_Tomo_fdprecond;
	    predata=(void*)recon;
	}
	if(gpu_pcg(&curecon->opdr, gpu_TomoL, recon, prefun, predata, rhs, &curecon->cgtmp_tomo, 
		   simu->parms->recon.warm_restart, parms->tomo.maxit, curecon->cgstream[0])>1){
	    warning("Tomo CG not converge.\n");
	}
	toc_test("TomoL CG");
    }break;
    case 2:
	curmv(curecon->opdr->m->p, 0, curecon->RMI, rhs->m->p, 'n', 1, curecon->cgstream[0]);
	break;
    default:
	error("Invalid");
    }
    curcellfree(rhs); rhs=NULL;
    if(!parms->gpu.fit || parms->save.opdr || parms->recon.split==2 || (recon->moao && !parms->gpu.moao)){
	cp2cpu(&simu->opdr, 0, curecon->opdr_vec, 1, curecon->cgstream[0]);
    }
    if(curecon->RFlgsx){
	curcell *focus=NULL;
	curcellmm(&focus, 0, curecon->RFlgsx, curecon->opdr_vec, "nn", 1, curecon->cgstream[0]);
	cp2cpu(&simu->focuslgsx, 0, focus, 1, curecon->cgstream[0]);
    }
    if(curecon->RFngsx){
	curcell *focus=NULL;
	curcellmm(&focus, 0, curecon->RFngsx, curecon->opdr_vec, "nn", 1, curecon->cgstream[0]);
	cp2cpu(&simu->focusngsx, 0, focus, 1, curecon->cgstream[0]);
    }
    cudaStreamSynchronize(curecon->cgstream[0]);
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
    curcell *rhs=NULL;
    G_CGFUN cg_fun;
    void *cg_data;
    if(parms->gpu.fit==1){//sparse matrix
	cumuv(&rhs, 0, &curecon->FR, curecon->opdr, 1);
	cg_fun=(G_CGFUN) cumuv;
	cg_data=&curecon->FL;
    }else{
	gpu_FitR(&rhs, 0, recon, curecon->opdr, 1);
	cg_fun=(G_CGFUN) gpu_FitL;
	cg_data=(void*)recon;
    }
    toc_test("FitR");
    switch(parms->fit.alg){
    case 0:
	cuchol_solve(curecon->dmfit->m->p, curecon->FCl, curecon->FCp, rhs->m->p, curecon->cgstream[0]);
	if(curecon->FUp){
	    curmat *tmp=curnew(curecon->FVp->ny, 1);
	    curmv(tmp->p, 0, curecon->FVp, rhs->m->p, 't', -1, curecon->cgstream[0]);
	    curmv(curecon->dmfit->m->p, 1, curecon->FUp, tmp->p, 'n', 1, curecon->cgstream[0]);
	    cudaStreamSynchronize(curecon->cgstream[0]);
	    curfree(tmp);
	}
	break;
    case 1:{
	double res;
	if((res=gpu_pcg(&curecon->dmfit, (G_CGFUN)cg_fun, cg_data, NULL, NULL, rhs, &curecon->cgtmp_fit,
			simu->parms->recon.warm_restart, parms->fit.maxit, curecon->cgstream[0]))>1){
	    warning("DM Fitting PCG not converge. res=%g\n", res);
	}
    }
	break;
    case 2:
	curmv(curecon->dmfit->m->p, 0, curecon->FMI, rhs->m->p, 'n', 1, curecon->cgstream[0]);
	break;
    default:
	error("Invalid");
    }
    cp2cpu(&simu->dmfit, 0, curecon->dmfit_vec, 1, curecon->cgstream[0]);
    toc_test("FitL CG");
    cudaStreamSynchronize(curecon->cgstream[0]);
    /*Don't free opdr. */
    curcellfree(rhs); rhs=NULL;
    toc_test("Fit");
}
void gpu_recon_mvm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    cp2gpu(&curecon->gradin, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
    curcellmm(&curecon->dmfit_vec, 0., curecon->MVM, curecon->gradin,"nn", 1., curecon->cgstream[0]);
    cp2cpu(&simu->dmerr, 0., curecon->dmfit_vec, 1., curecon->cgstream[0]);
}
