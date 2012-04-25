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
curecon_t *curecon=NULL;
#undef TIMING
#define TIMING 0
#if !TIMING
#undef TIC
#undef tic
#undef toc
#define TIC
#define tic
#define toc(A)
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
void gpu_setup_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    gpu_set(gpu_recon);
    if(!curecon){
	curecon=(curecon_t*)calloc(1, sizeof(curecon_t));
    }
    STREAM_NEW(curecon->cgstream);
    HANDLE_NEW(curecon->cghandle, curecon->cgstream);
    if(parms->recon.mvm && (!parms->gpu.tomo || !parms->gpu.fit)){
	return; /*Use CPU to assemble MVM*/
    }
    if(parms->recon.alg!=0){
	error("Only MVR is implemented in GPU\n");
    }
    cuwloc_t *cupowfs=cudata->powfs;
    cudaStream_t stream;
    STREAM_NEW(stream);
    cublasHandle_t handle;
    DO(cublasCreate(&handle));
    DO(cublasSetStream(handle, stream));
    curecon->wfsstream=(cudaStream_t*)calloc(parms->nwfsr, sizeof(cudaStream_t));
    curecon->wfshandle=(cublasHandle_t*)calloc(parms->nwfsr, sizeof(cublasHandle_t));
    curecon->wfssphandle=(cusparseHandle_t*)calloc(parms->nwfsr, sizeof(cusparseHandle_t));
    curecon->wfsevent=(cudaEvent_t*)calloc(parms->nwfsr, sizeof(cudaEvent_t));

    curecon->fitstream=(cudaStream_t*)calloc(parms->fit.nfit, sizeof(cudaStream_t));
    curecon->fitsphandle=(cusparseHandle_t*)calloc(parms->fit.nfit, sizeof(cusparseHandle_t));
    curecon->fithandle=(cublasHandle_t*)calloc(parms->fit.nfit, sizeof(cublasHandle_t));
    
    curecon->psstream=(cudaStream_t*)calloc(recon->npsr, sizeof(cudaStream_t));
    curecon->dmstream=(cudaStream_t*)calloc(recon->ndm, sizeof(cudaStream_t));
    curecon->dmhandle=(cublasHandle_t*)calloc(recon->ndm, sizeof(cublasHandle_t));
    curecon->dmsphandle=(cusparseHandle_t*)calloc(recon->ndm, sizeof(cusparseHandle_t));
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	STREAM_NEW(curecon->wfsstream[iwfs]);
	HANDLE_NEW(curecon->wfshandle[iwfs], curecon->wfsstream[iwfs]);
	SPHANDLE_NEW(curecon->wfssphandle[iwfs], curecon->wfsstream[iwfs]);
	DO(cudaEventCreateWithFlags(&curecon->wfsevent[iwfs], cudaEventDisableTiming));
    }
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	STREAM_NEW(curecon->fitstream[ifit]);
	HANDLE_NEW(curecon->fithandle[ifit], curecon->fitstream[ifit]);
	SPHANDLE_NEW(curecon->fitsphandle[ifit], curecon->fitstream[ifit]);
    }
    for(int ips=0; ips<recon->npsr; ips++){
	STREAM_NEW(curecon->psstream[ips]);
    }
    for(int idm=0; idm<recon->ndm; idm++){
	STREAM_NEW(curecon->dmstream[idm]);
	HANDLE_NEW(curecon->dmhandle[idm], curecon->dmstream[idm]);
	SPHANDLE_NEW(curecon->dmsphandle[idm],curecon->dmstream[idm]);
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
		const int use_mat=parms->tomo.pos==2;
		if(use_mat){//normally true
		    dsp *GP=sptrans(recon->GP->p[ipowfs]);
		    spint *pp=GP->p;
		    spint *pi=GP->i;
		    double *px=GP->x;
		    dmat *partx=NULL;
		    dmat *party=NULL;
		    partx=dnew(9, nsa);
		    party=dnew(9, nsa);

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
			    if(zx<0 || zx>2 || zy<0 || zy>2){
				warning("isa=%d, zxy=%d %d\n", isa, zx, zy);
			    }
			    if(zx<0) zx=0;
			    if(zx>2) zx=2;
			    if(zy<0) zy=0;
			    if(zy>2) zy=2;
			    if(ic<nsa){/*x */
				partx->p[9*isa+zx+zy*3]+=px[ir];
			    }else{/*y */
				party->p[9*isa+zx+zy*3]+=px[ir];
			    }
			}
		    }
		    cp2gpu(&cupowfs[ipowfs].GPpx, partx);
		    cp2gpu(&cupowfs[ipowfs].GPpy, party);
		    dfree(partx);
		    dfree(party);
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
		info2("wfs %d is reusing nea for %d\n", iwfs, iwfs0);
		curecon->neai->p[iwfs]=curref(curecon->neai->p[iwfs0]);
	    }else{
		info2("wfs %d has different nea than %d\n", iwfs, iwfs0);
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
	if(recon->PDF && !curecon->PDF){
	    cp2gpu(&curecon->PDF, recon->PDF);
	}
	if(parms->tomo.precond==1){/*fdpcg*/
	    FDPCG_T *fdpcg=recon->fdpcg;
	    int nb=fdpcg->Mbinv->nx;
	    int bs=fdpcg->Mbinv->p[0]->nx;
	    cp2gpu(&curecon->fd_perm, fdpcg->perm, nb*bs);//not bs*bs
	    //curecon->fd_Mb=cuccellnew(nb, 1, bs, bs);
	    cp2gpu(&curecon->fd_Mb, fdpcg->Mbinv);
	    curecon->fd_nxtot=nb*bs;
	    int nps=recon->npsr;
	    int count=0;
	    int osi=-1;
	    int start[nps];
	    for(int ips=0; ips<nps; ips++){
		if(osi != parms->atmr.os[ips]){
		    start[count]=ips;
		    osi = parms->atmr.os[ips];
		    count++;
		}
	    }
	    curecon->fd_fft=(cufftHandle*)calloc(count, sizeof(cufftHandle));
	    curecon->fd_fftnc=count;
	    curecon->fd_fftips=(int*)calloc(count+1, sizeof(int));
	    for(int ic=0; ic<count; ic++){
		curecon->fd_fftips[ic]=start[ic];
	    }
	    curecon->fd_fftips[count]=nps;
	    for(int ic=0; ic<count; ic++){
		int ncomp[2];
		ncomp[0]=recon->xnx[start[ic]];
		ncomp[1]=recon->xny[start[ic]];
		cufftPlanMany(&curecon->fd_fft[ic], 2, ncomp, NULL, 1, 0, NULL, 1, 0, 
			      CUFFT_C2C, curecon->fd_fftips[ic+1]-curecon->fd_fftips[ic]);
		cufftSetStream(curecon->fd_fft[ic], curecon->cgstream);
	    }
	    /* notice: performance may be improved by using
	       R2C FFTs instead of C2C. Need to update perm
	       and Mbinv to use R2C.*/
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
	curecon->opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
	for(int ips=0; ips<recon->npsr; ips++){
	    info("opdr[%d] is %dx%d\n", ips, curecon->opdr->p[ips]->nx, curecon->opdr->p[ips]->ny);
	}
	const int nwfs=parms->nwfsr;
	int nxp=recon->pmap->nx;
	int nyp=recon->pmap->ny;
  
	curecon->opdwfs=curcellnew(nwfs, 1);
	curecon->grad=curcellnew(nwfs, 1);/*intermediate. */
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    const int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip) continue;
	    curecon->opdwfs->p[iwfs]=curnew(nxp, nyp);
	    const int nsa=powfs[ipowfs].pts->nsa;
	    curecon->grad->p[iwfs]=curnew(nsa*2,1);
	}
	if(parms->tomo.precond==1){
	    curecon->fd_xhat1=cuccellnew(recon->npsr, 1, recon->xnx, recon->xny);
	    curecon->fd_xhat2=cuccellnew(recon->npsr, 1, recon->xnx, recon->xny);
	}
    }
    if(parms->gpu.fit){
	if(parms->gpu.fit==1){ /*For fitting using sparse matrix*/
	    cp2gpu(&curecon->FR, &recon->FR);
	    cp2gpu(&curecon->FL, &recon->FL);
	    curecon->FR.fitstream=curecon->fitstream;
	    curecon->FR.fithandle=curecon->fithandle;
	    curecon->FR.fitsphandle=curecon->fitsphandle;
	    curecon->FR.dmstream=curecon->dmstream;
	    curecon->FR.dmhandle=curecon->dmhandle;
	    curecon->FR.dmsphandle=curecon->dmsphandle;
	    curecon->FL.fitstream=curecon->fitstream;
	    curecon->FL.fithandle=curecon->fithandle;
	    curecon->FL.fitsphandle=curecon->fitsphandle;
	    curecon->FL.dmstream=curecon->dmstream;
	    curecon->FL.dmhandle=curecon->dmhandle;
	    curecon->FL.dmsphandle=curecon->dmsphandle;
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
    STREAM_DONE(stream);
    cublasDestroy(handle);
    gpu_print_mem("recon init");
}
static __global__ void 
copy_mvm_do(float *dest, long stride, float *in, long n){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	dest[i*stride]=in[i];
    }
}
void gpu_setup_recon_mvm(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    gpu_set(gpu_recon);
    if(parms->recon.alg!=0){
	error("Please adept to LSR\n");
    } 
    if(!parms->gpu.tomo || !parms->gpu.fit){
	setup_recon_mvr_mvm(recon, parms, powfs);
	cp2gpu(&curecon->MVM, recon->MVM);
	dcellfree(recon->MVM);
    }else{
	info2("Assembling MVR MVM in GPU\n");
	int ntotact=0;
	int ntotgrad=0;
	const int nwfs=parms->nwfs;
	const int ndm=parms->ndm;
	for(int idm=0; idm<ndm; idm++){
	    ntotact+=recon->anloc[idm];
	}    
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    ntotgrad+=recon->ngrad[iwfs];
	}
	curcell *eyec;
	if(parms->fit.square){
	    eyec=curcellnew(ndm, 1, recon->anx, recon->any);
	}else{
	    eyec=curcellnew(ndm, 1, recon->anloc, (long*)0);
	}
	float eye2[2]={0,1.};
	float eye1[1]={1.};
	curecon->MVM=curcellnew(ndm, nwfs);
	for(int idm=0; idm<ndm; idm++){
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(!parms->powfs[ipowfs].skip){
		    curecon->MVM->p[idm+iwfs*ndm]
			=curnew(recon->anloc[idm], powfs[ipowfs].saloc->nloc*2);
		}
	    }
	}
	curcell *dmfit=curcellnew(curecon->dmfit);
	curcell *opdx=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
	curcell *opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
	curcell *grad=curcellnew(parms->nwfs, 1, recon->ngrad, (long*)0);
	G_CGFUN cg_fun;
	void *cg_data;
	if(parms->gpu.fit==1){//sparse matrix
	    cg_fun=(G_CGFUN) cumuv;
	    cg_data=&curecon->FL;
	}else{
	    cg_fun=(G_CGFUN) gpu_FitL;
	    cg_data=recon;
	}
	int curdm=0, curact=0;
	for(int iact=0; iact<ntotact; iact++){
	    TIC;tic;
	    if(!detached){
		info2("%6d of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", iact, ntotact);
	    }
	    if(iact){
		cudaMemcpyAsync(eyec->m->p+iact-1, eye2, 2*sizeof(float),
				cudaMemcpyHostToDevice, curecon->cgstream);
	    }else{
		cudaMemcpyAsync(eyec->m->p+iact, eye1, sizeof(float), 
				cudaMemcpyHostToDevice, curecon->cgstream);
	    }
	    toc("copy");
	    /*Fitting operator*/
	    switch(parms->fit.alg){
	    case 0:
		cuchol_solve(dmfit->m->p, curecon->FCl, curecon->FCp, eyec->m->p, curecon->cgstream);
		if(curecon->FUp){
		    curmat *tmp=curnew(curecon->FVp->ny, 1);
		    curmv(tmp->p, 0, curecon->FVp, eyec->m->p, 't', -1, curecon->cghandle);
		    curmv(curecon->dmfit->m->p,1,curecon->FUp,tmp->p,'n',1,curecon->cghandle);
		    cudaStreamSynchronize(curecon->cgstream);
		    curfree(tmp);
		}
		break;
	    case 1:
		curcellzero(dmfit, curecon->cgstream);//temp
		if(gpu_pcg(&dmfit, (G_CGFUN)cg_fun, cg_data, NULL, NULL, eyec,
			   parms->recon.warm_restart, parms->fit.maxit, curecon->cgstream)){
		    error("Fit CG failed\n");
		}
		break;
	    case 2:
		curmv(dmfit->m->p, 0, curecon->FMI, eyec->m->p, 'n', 1, curecon->cghandle);
		break;
	    default:
		error("Invalid");
	    }
	    toc("FitL");
	    /*Transpose of fitting operator*/
	    if(parms->gpu.fit==1){//sparse matrix
		cumuv_trans(&opdx, 0, &curecon->FR, dmfit, 1);
	    }else{
		gpu_FitRt(&opdx, 0, recon, dmfit, 1);
	    }
	    
	    toc("FitRt");
	    /*Tomography*/
	    G_PREFUN prefun=NULL;
	    void *predata=NULL;
	    if(parms->tomo.precond==1){
		prefun=gpu_Tomo_fdprecond;
		predata=(void*)recon;
	    }
	    curcellzero(opdr, curecon->cgstream);
	    switch(parms->tomo.alg){
	    case 0:
		if(!opdr->m || !opdx->m){
		    error("opdr and opdx must be continuous\n");
		}
		cuchol_solve(opdr->m->p, curecon->RCl, curecon->RCp, opdx->m->p, curecon->cgstream);
		if(curecon->RUp){
		    curmat *tmp=curnew(curecon->RVp->ny, 1);
		    curmv(tmp->p, 0, curecon->RVp, opdx->m->p, 't', -1, curecon->cghandle);
		    curmv(opdr->m->p, 1, curecon->RUp, tmp->p, 'n', 1, curecon->cghandle);
		    cudaStreamSynchronize(curecon->cgstream);
		    curfree(tmp);
		}
		break;
	    case 1:{
		int disablelrt=curecon->disablelrt;
		curecon->disablelrt=1;
		/*disable the t/t removal lrt in split tomo that creats problem in fdpcg mode*/
		if(gpu_pcg(&opdr, gpu_TomoL, recon, prefun, predata, opdx, 
			   parms->recon.warm_restart, parms->tomo.maxit, curecon->cgstream)){
		    error("Tomo CG failed\n");
		}
		curecon->disablelrt=disablelrt;
	    }
		break;
	    case 2:
		curmv(opdr->m->p, 0, curecon->RMI, opdx->m->p, 'n', 1, curecon->cghandle);
		break;
	    default:
		error("Invalid");
	    }
	    toc("TomoL");
	    /*Right hand size*/
	    gpu_TomoRt(&grad, 0, recon, opdr, 1);
	    toc("TomoR");
	    //copy grad out and transpose
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		curmat *to=curecon->MVM->p[curdm+ndm*iwfs];
		if(to){
		    copy_mvm_do<<<DIM(ntotgrad, 256), 0, curecon->cgstream>>>
			(to->p+curact, recon->anloc[curdm], grad->p[iwfs]->p, recon->ngrad[iwfs]);
		}
	    }
	    /*{
	      curcellwrite(eyec, "eyec_%d", iact);
	      curcellwrite(dmfit, "dmfit_%d", iact);
	      curcellwrite(opdx, "opdx_%d", iact);
	      curcellwrite(opdr, "opdr_%d", iact);
	      curcellwrite(grad, "grad_%d", iact);
	      }*/
	    curact++;
	    if(curact>=recon->anloc[curdm]){
		curdm++;
		curact=0;
	    }
	    toc("Copy");
	}//for iact
	warning("\n\nFree the various data for GPU tomo/fit\n\n");
	curcellfree(dmfit);
	curcellfree(opdx);
	curcellfree(opdr);
	curcellfree(grad);
    }
    if(parms->save.setup){
	curcellwrite(curecon->MVM, "%s/MVM.bin", dirsetup);
    }
}
/*update reconstruction parameters after slodar.*/
void gpu_update_recon(const PARMS_T *parms, RECON_T *recon){
    gpu_set(gpu_recon);
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
    curcellzero(curecon->opdr, 0);
    curcellzero(curecon->dmfit, 0);
    if(curecon->moao_wfs){
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    curcellzero(curecon->moao_wfs[iwfs], 0);
	}
    }
    if(curecon->moao_evl){
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    curcellzero(curecon->moao_evl[ievl], 0);
	}
    }
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	curcellzero(cudata->moao_wfs,0);
	curcellzero(cudata->moao_evl,0);
    }
    cudaStreamSynchronize(0);
}

void gpu_tomo(SIM_T *simu){
    gpu_set(gpu_recon);
    TIC;tic;
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    if(parms->tomo.pos!=2){
	TO_IMPLEMENT;
    }
    if(curecon->PDF){
	TO_IMPLEMENT;
    }
    /*first send gradients to GPU. can be skipped if keep grad in gpu. fast though. */
    /*Create temporary memory */
    curecon->reconisim=simu->reconisim;
#if 0
    {
	/*Debugging. */
	dcell *rhsc=NULL;
	dcell *lc=NULL;
	muv(&rhsc, &recon->RR, simu->gradlastol, 1);
	dcellwrite(rhsc, "CPU_TomoR");
	muv(&lc, &recon->RL, rhsc, 1);
	dcellwrite(lc, "CPU_TomoL");
	muv(&lc, &recon->RL, rhsc, -1);
	dcellwrite(lc, "CPU_TomoL2");
	dcellzero(lc);
	muv(&lc, &recon->RL, rhsc, 1);
	dcellwrite(lc, "CPU_TomoL3");
	if(parms->tomo.alg==1 && parms->tomo.precond==1){
	    dcell *lp=NULL;
	    fdpcg_precond(&lp, recon, lc);
	    dcellwrite(lp, "CPU_TomoP");
	    fdpcg_precond(&lp, recon, lc);
	    dcellwrite(lp, "CPU_TomoP2");
	}
	cp2gpu(&curecon->gradin, simu->gradlastol);
	curcell *rhsg=NULL;
	curcell *lg=NULL;
	gpu_TomoR(&rhsg, 0, recon, curecon->gradin, 1);
	curcellwrite(rhsg, "GPU_TomoR");
	gpu_TomoL(&lg, 0, recon, rhsg, 1);
	curcellwrite(lg, "GPU_TomoL");
	gpu_TomoL(&lg, 1, recon, rhsg, -1);
	curcellwrite(lg, "GPU_TomoL2");
	gpu_TomoL(&lg, 0, recon, rhsg, 1);
	curcellwrite(lg, "GPU_TomoL3");
	if(parms->tomo.alg==1 && parms->tomo.precond==1){
	    curcell *lp=NULL;
	    gpu_Tomo_fdprecond(&lp, recon, lg, curecon->cgstream);
	    curcellwrite(lp, "GPU_TomoP");
	    gpu_Tomo_fdprecond(&lp, recon, lg, curecon->cgstream);
	    curcellwrite(lp, "GPU_TomoP2");
	}
	CUDA_SYNC_DEVICE;
	exit(0);
    }
#endif
    toc("Before gradin");
    cp2gpu(&curecon->gradin, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
    toc("Gradin");
    curcell *rhs=NULL;
    gpu_TomoR(&rhs, 0, recon, curecon->gradin, 1);
    toc("TomoR");
    switch(parms->tomo.alg){
    case 0:
	if(!curecon->opdr->m){
	    error("opdr must be continuous\n");
	}
	if(!rhs->m){
	    error("rhs must be continuous\n");
	}
	cuchol_solve(curecon->opdr->m->p, curecon->RCl, curecon->RCp, rhs->m->p, curecon->cgstream);
	if(curecon->RUp){
	    curmat *tmp=curnew(curecon->RVp->ny, 1);
	    curmv(tmp->p, 0, curecon->RVp, rhs->m->p, 't', -1, curecon->cghandle);
	    curmv(curecon->opdr->m->p, 1, curecon->RUp, tmp->p, 'n', 1, curecon->cghandle);
	    cudaStreamSynchronize(curecon->cgstream);
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
	if(gpu_pcg(&curecon->opdr, gpu_TomoL, recon, prefun, predata, rhs, 
		   simu->parms->recon.warm_restart, parms->tomo.maxit, curecon->cgstream)){
	    error("Tomo CG failed\n");
	}
	toc("TomoL CG");
    }break;
    case 2:
	curmv(curecon->opdr->m->p, 0, curecon->RMI, rhs->m->p, 'n', 1, curecon->cghandle);
	break;
    default:
	error("Invalid");
    }
    curcellfree(rhs); rhs=NULL;
    if(!parms->gpu.fit || parms->save.opdr || parms->recon.split==2 || (recon->moao && !parms->gpu.moao)){
	cp2cpu(&simu->opdr, 0, curecon->opdr, 1, curecon->cgstream);
	cudaStreamSynchronize(curecon->cgstream);
	for(int i=0; i<simu->opdr->nx; i++){
	    simu->opdr->p[i]->nx=simu->opdr->p[i]->nx*simu->opdr->p[i]->ny;
	    simu->opdr->p[i]->ny=1;
	}
    }
    toc("Tomo");
}
void gpu_fit(SIM_T *simu){
    TIC;tic;
    gpu_set(gpu_recon);
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(!parms->gpu.tomo){
	cp2gpu(&curecon->opdr, simu->opdr);
    }
#if 0
    {
	/*Debugging. */
	dcell *rhsc=NULL;
	dcell *lc=NULL;
	if(!simu->opdr){
	    cp2cpu(&simu->opdr, 0, curecon->opdr, 1, 0);
	    for(int i=0; i<simu->opdr->nx; i++){
		simu->opdr->p[i]->nx=simu->opdr->p[i]->nx*simu->opdr->p[i]->ny;
		simu->opdr->p[i]->ny=1;
	    }
	}
	dcellwrite(simu->opdr, "opdr");
	muv(&rhsc, &recon->FR, simu->opdr, 1);
	dcellwrite(rhsc, "CPU_FitR");
	muv(&lc, &recon->FL, rhsc, 1);
	dcellwrite(lc, "CPU_FitL");
	muv(&lc, &recon->FL, rhsc, -1);
	dcellwrite(lc, "CPU_FitL2");
	dcellzero(lc);
	muv(&lc, &recon->FL, rhsc, 1);
	dcellwrite(lc, "CPU_FitL3");
	dcellzero(lc);
	muv_solve(&lc, &recon->FL, NULL, rhsc);
	dcellwrite(lc, "CPU_FitCG");
	/*dcell *lhs=NULL;
	muv_trans(&lhs, &recon->FR, rhsc, 1);
	dcellwrite(lhs, "CPU_FitRt");*/
	curcell *rhsg=NULL;
	curcell *lg=NULL;
	gpu_FitR(&rhsg, 0, recon, curecon->opdr, 1);
	curcellwrite(rhsg, "GPU_FitR");
	gpu_FitL(&lg, 0, recon, rhsg, 1);
	curcellwrite(lg, "GPU_FitL");
	gpu_FitL(&lg, 1, recon, rhsg, -1);
	curcellwrite(lg, "GPU_FitL2");
	gpu_FitL(&lg, 0, recon, rhsg, 1);
	curcellwrite(lg, "GPU_FitL3");
	/*curcell *lhsg=NULL;
	gpu_FitRt(&lhsg, 0, recon, rhsg, 1);
	curcellwrite(lhsg, "GPU_FitRt");*/
	curcellzero(lg, curecon->cgstream);
	gpu_pcg(&lg, (G_CGFUN)gpu_FitL, (void*)recon, NULL, NULL, rhsg,
		simu->parms->recon.warm_restart, parms->fit.maxit, curecon->cgstream);
	curcellwrite(lg, "GPU_FitCG");
	CUDA_SYNC_DEVICE;
	exit(0);
    }
#endif
    toc("Before FitR");
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
    toc("FitR");
    switch(parms->fit.alg){
    case 0:
	cuchol_solve(curecon->dmfit->m->p, curecon->FCl, curecon->FCp, rhs->m->p, curecon->cgstream);
	if(curecon->FUp){
	    curmat *tmp=curnew(curecon->FVp->ny, 1);
	    curmv(tmp->p, 0, curecon->FVp, rhs->m->p, 't', -1, curecon->cghandle);
	    curmv(curecon->dmfit->m->p, 1, curecon->FUp, tmp->p, 'n', 1, curecon->cghandle);
	    cudaStreamSynchronize(curecon->cgstream);
	    curfree(tmp);
	}
	break;
    case 1:
	if(gpu_pcg(&curecon->dmfit, (G_CGFUN)cg_fun, cg_data, NULL, NULL, rhs,
		   simu->parms->recon.warm_restart, parms->fit.maxit, curecon->cgstream)){
	    error("DM Fitting PCG failed\n");
	}
	break;
    case 2:
	curmv(curecon->dmfit->m->p, 0, curecon->FMI, rhs->m->p, 'n', 1, curecon->cghandle);
	break;
    default:
	error("Invalid");
    }
    cp2cpu(&simu->dmfit, 0, curecon->dmfit_vec, 1, curecon->cgstream);
    toc("FitL CG");
    cudaStreamSynchronize(curecon->cgstream);
    /*Don't free opdr. */
    curcellfree(rhs); rhs=NULL;
    toc("Fit");
}
void gpu_recon_mvm(SIM_T *simu){
    gpu_set(gpu_recon);
    const PARMS_T *parms=simu->parms;
    cp2gpu(&curecon->gradin, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
    curcellmm(&curecon->dmfit_vec, 0., curecon->MVM, curecon->gradin,"nn", 1., curecon->cghandle);
    cp2cpu(&simu->dmerr, 0., curecon->dmfit_vec, 1., curecon->cgstream);
    if(parms->tomo.psol){
	dcelladd(&simu->dmerr, 1, simu->dmint->mint[parms->dbg.psol?0:1], -1);
    }
}
