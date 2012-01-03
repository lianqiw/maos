/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
curecon_t *curecon;
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

void gpu_setup_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    gpu_set(0);
    cuwloc_t *cupowfs=cudata->powfs;
    cudaStream_t stream;
    STREAM_NEW(stream);
    cublasHandle_t handle;
    DO(cublasCreate(&handle));
    DO(cublasSetStream(handle, stream));

    curecon=(curecon_t*)calloc(1, sizeof(curecon_t));
    curecon->wfsstream=(cudaStream_t*)calloc(parms->nwfsr, sizeof(cudaStream_t));
    curecon->wfshandle=(cublasHandle_t*)calloc(parms->nwfsr, sizeof(cublasHandle_t));
    curecon->wfssphandle=(cusparseHandle_t*)calloc(parms->nwfsr, sizeof(cusparseHandle_t));
    curecon->wfsevent=(cudaEvent_t*)calloc(parms->nwfsr, sizeof(cudaEvent_t));

    curecon->fitstream=(cudaStream_t*)calloc(parms->fit.nfit, sizeof(cudaStream_t));
    curecon->fitsphandle=(cusparseHandle_t*)calloc(parms->fit.nfit, sizeof(cusparseHandle_t));
    curecon->fithandle=(cublasHandle_t*)calloc(parms->fit.nfit, sizeof(cublasHandle_t));

    curecon->psstream=(cudaStream_t*)calloc(recon->npsr, sizeof(cudaStream_t));
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	STREAM_NEW(curecon->wfsstream[iwfs]);
	DO(cublasCreate(&curecon->wfshandle[iwfs]));
	DO(cublasSetStream(curecon->wfshandle[iwfs], curecon->wfsstream[iwfs]));
	DO(cusparseCreate(&curecon->wfssphandle[iwfs]));
	DO(cusparseSetKernelStream(curecon->wfssphandle[iwfs], curecon->wfsstream[iwfs]));
	DO(cudaEventCreateWithFlags(&curecon->wfsevent[iwfs], cudaEventDisableTiming));
    }
    for(int ifit=0; ifit<parms->fit.nfit; ifit++){
	STREAM_NEW(curecon->fitstream[ifit]);
	DO(cusparseCreate(&curecon->fitsphandle[ifit]));
	DO(cusparseSetKernelStream(curecon->fitsphandle[ifit], curecon->fitstream[ifit]));
	DO(cublasCreate(&curecon->fithandle[ifit]));
	DO(cublasSetStream(curecon->fithandle[ifit], curecon->fitstream[ifit]));
    }
    for(int ips=0; ips<recon->npsr; ips++){
	STREAM_NEW(curecon->psstream[ips]);
    }
    STREAM_NEW(curecon->cgstream);
    DO(cublasCreate(&curecon->cghandle));
    DO(cublasSetStream(curecon->cghandle, curecon->cgstream));
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
		if(use_mat){
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
			for(int ir=pp[ic]; ir<pp[ic+1]; ir++){
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
		    gpu_dmat2cu(&cupowfs[ipowfs].GPpx, partx);
		    gpu_dmat2cu(&cupowfs[ipowfs].GPpy, party);
		    dfree(partx);
		    dfree(party);
		    spfree(GP);
		}else{/*use sparse */
		    gpu_sp2dev(&cupowfs[ipowfs].GP, recon->GP->p[ipowfs]);
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
		curecon->zzv[ips]=pow(laplacian_coef(r0,wt,dx),2)*TOMOSCALE;
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
		    for(int ir=pp[ic]; ir<pp[ic+1]; ir++){
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
		DO(cudaMemcpy(curecon->neai->p[iwfs]->p, neai, 3*nsa*sizeof(float), cudaMemcpyDefault));
		free(neai);
	    }
	}/*for iwfs */
	CUDA_SYNC_DEVICE;
	if(recon->PTT && !curecon->PTT){
	    gpu_dcell2cu(&curecon->PTT, recon->PTT);
	}
	if(recon->PDF && !curecon->PDF){
	    gpu_dcell2cu(&curecon->PDF, recon->PDF);
	}
	if(parms->tomo.precond==1){/*fdpcg*/
	    FDPCG_T *fdpcg=recon->fdpcg;
	    int nb=fdpcg->Mbinv->nx;
	    int bs=fdpcg->Mbinv->p[0]->nx;
	    gpu_long2dev(&curecon->fd_perm, fdpcg->perm, nb*bs);//not bs*bs
	    curecon->fd_Mb=cuccellnew(nb, 1, bs, bs);
	    fcomplex *tmp;
	    cudaMallocHost(&tmp, nb*bs*bs*sizeof(fcomplex));
	    for(int ib=0; ib<nb; ib++){
		dcomplex *in=(dcomplex*)fdpcg->Mbinv->p[ib]->p;
		fcomplex *tmp2=tmp+ib*(bs*bs);
		for(int i=0; i<bs*bs; i++){
		    tmp2[i]=make_cuFloatComplex((float)cuCreal(in[i]), (float)cuCimag(in[i]));
		}
	    }
	    DO(cudaMemcpy(curecon->fd_Mb->p[0]->p, tmp, nb*bs*bs*sizeof(fcomplex), cudaMemcpyHostToDevice));
	    curecon->fd_nxtot=nb*bs;
	    cudaFreeHost(tmp);
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
	    }
	    /* notice: performance may be improved by using
	       R2C FFTs instead of C2C. Need to update perm
	       and Mbinv to use R2C.*/
	}
    }
    if(parms->gpu.fit){
	/*For fitting */
	/*
	  gpu_dmat2cu(&curecon->W1, recon->W1);
	  if(recon->W0){
	  spint *pp=recon->W0->p;
	  spint *pi=recon->W0->i;
	  double *px=recon->W0->x;
	  dsp *W0new=spnew(recon->W0->m, recon->W0->n, recon->W0->nzmax);
	  spint *pp2=W0new->p;
	  spint *pi2=W0new->i;
	  double *px2=W0new->x;
	  int *full;
	  cudaMallocHost(&full, recon->W0->n*sizeof(int));
	  double W1max=dmax(recon->W1);
	  double thres=W1max*(1.f-1e-6);
	  curecon->W0v=(float)(W1max*4./9.);//max of W0 is 4/9 of max of W1. 
	  info("W0v=%g\n", curecon->W0v);
	  int count=0;
	  int count2=0;
	  for(int ic=0; ic<recon->W0->n; ic++){
	  pp2[ic]=count;
	  if(recon->W1->p[ic]>thres){
	  full[count2]=ic;
	  count2++;
	  }else{
	  int nv=pp[ic+1]-pp[ic];
	  memcpy(pi2+count, pi+pp[ic], sizeof(spint)*nv);
	  memcpy(px2+count, px+pp[ic], sizeof(double)*nv);
	  count+=nv;
	  }
	  }
	  pp2[recon->W0->n]=count;
	  W0new->nzmax=count;
	  dsp *W0new2=sptrans(W0new);
	  gpu_sp2dev(&curecon->W0p, W0new2);
	  gpu_int2dev(&curecon->W0f, full, count2);
	  curecon->nW0f=count2;
	  spfree(W0new);
	  spfree(W0new2);
	  cudaFreeHost(full);
	  }else{
	  error("W0, W1 is required\n");
	  }*/
	gpu_muv2dev(&curecon->FR, &recon->FR);
	gpu_muv2dev(&curecon->FL, &recon->FL);
    }
    STREAM_DONE(stream);
    cublasDestroy(handle);
    gpu_print_mem("recon init");
}
void gpu_recon_reset(){
    curcellfree(curecon->opdr); curecon->opdr=NULL;
    curcellfree(curecon->dmfit); curecon->dmfit=NULL;
}
void cumuv(curcell **out, float beta, cumuv_t *A, const curcell *in, float alpha){
    if(!A->Mt) error("A->M Can not be empty\n");
    if(A->U->ny>1 || A->V->ny>1) error("Not handled yet\n");
    if(!*out){
	*out=curcellnew(A->Mt->ny, 1);
	for(int ix=0; ix<A->Mt->ny; ix++){
	    (*out)->p[ix]=curnew(A->Mt->p[ix*A->Mt->nx]->ny, 1);
	}
    }else if(fabsf(beta-1.)>1e-5){
	curcellscale(*out, beta, curecon->fitstream[0]);
    }
    cusparseHandle_t sphandle=curecon->fitsphandle[0];
    cublasHandle_t handle=curecon->fithandle[0];
    for(int iy=0; iy<A->Mt->nx; iy++){
	for(int ix=0; ix<A->Mt->ny; ix++){
	    cusptmul((*out)->p[ix]->p, A->Mt->p[iy+ix*A->Mt->nx], in->p[iy]->p, alpha, sphandle);
	}
    }
    curmat *tmp=curnew(A->V->p[0]->ny, 1);
    for(int iy=0; iy<A->V->nx; iy++){
	curmv(&tmp, 1.f, A->V->p[iy], in->p[iy], 't', 1, handle);
    }
    for(int ix=0; ix<A->U->nx; ix++){
	curmv(&(*out)->p[ix], 1.f, A->U->p[ix], tmp, 'n', -alpha, handle);
    }
    cudaStreamSynchronize(curecon->fitstream[0]);
    curfree(tmp);
}
__global__ void nothing(void){

}
void gpu_tomo(SIM_T *simu){
    gpu_set(0);
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
    int nxp=recon->pmap->nx;
    int nyp=recon->pmap->ny;
   
    const int nwfs=parms->nwfsr;
    /*Create temporary memory */
    curecon->reconisim=simu->reconisim;
    curecon->opdwfs=curcellnew(nwfs, 1);
    curecon->grad=curcellnew(nwfs, 1);/*intermediate. */
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	curecon->opdwfs->p[iwfs]=curnew(nxp, nyp);
	const int nsa=simu->powfs[ipowfs].pts->nsa;
	curecon->grad->p[iwfs]=curnew(nsa*2,1);
    }
    if(0){
	/*Debugging. */
	dcell *rhsc=NULL;
	dcell *lc=NULL;
	muv(&rhsc, &recon->RR, simu->gradlastol, 1);
	muv(&lc, &recon->RL, rhsc, -1);
	dcellwrite(rhsc, "CPU_TomoR");
	dcellwrite(lc, "CPU_TomoL");

	gpu_dcell2cu(&curecon->gradin, simu->gradlastol);
	curcell *rhsg=NULL;
	curcell *lg=NULL;
	gpu_TomoR(&rhsg, simu, curecon->gradin, 1);
	curcellwrite(rhsg, "GPU_TomoR");
	gpu_TomoL(&lg, 1, simu, rhsg, -1);
	curcellwrite(lg, "GPU_TomoL");
	CUDA_SYNC_DEVICE;
	exit(0);
    }
    toc("Before gradin");tic;
    gpu_dcell2cu(&curecon->gradin, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
    toc("Gradin");tic;
    curcell *rhs=NULL;
    gpu_TomoR(&rhs, simu, curecon->gradin, 1);
    toc("TomoR");
    G_PREFUN prefun=NULL;
    void *predata=NULL;
    if(parms->tomo.precond==1){
	curecon->fd_xhat1=cuccellnew(recon->npsr, 1, recon->xnx, recon->xny);
	curecon->fd_xhat2=cuccellnew(recon->npsr, 1, recon->xnx, recon->xny);
	prefun=gpu_Tomo_fdprecond;
	predata=(void*)simu;
    }
    tic;
    if(gpu_pcg(&curecon->opdr, gpu_TomoL, simu, prefun, predata, rhs, 
	       simu->parms->recon.warm_restart, parms->tomo.maxit, curecon->cgstream)){
	dcellwrite(simu->gradlastol, "tomo_gradlastol_%d", simu->reconisim);
	curcellwrite(curecon->gradin, "tomo_gradin_%d", simu->reconisim);
	curcellwrite(rhs, "tomo_rhs_%d", simu->reconisim);
	exit(1);
    }
    toc("TomoL CG");tic;
    if(parms->tomo.precond==1){
	cuccellfree(curecon->fd_xhat1);
	cuccellfree(curecon->fd_xhat2);
    }
    curcellfree(rhs); rhs=NULL;
    curcellfree(curecon->opdwfs); curecon->opdwfs=NULL;
    curcellfree(curecon->grad);   curecon->grad=NULL;
    curcellfree(curecon->gradin); curecon->gradin=NULL;
    if(!parms->gpu.fit || parms->save.opdr || parms->recon.split==2){
	gpu_curcell2d(&simu->opdr, curecon->opdr, curecon->cgstream);
	cudaStreamSynchronize(curecon->cgstream);
	for(int i=0; i<simu->opdr->nx; i++){
	    simu->opdr->p[i]->nx=simu->opdr->p[i]->nx*simu->opdr->p[i]->ny;
	    simu->opdr->p[i]->ny=1;
	}
    }
}
void gpu_fit(SIM_T *simu){
    TIC;tic;
    gpu_set(0);
    const PARMS_T *parms=simu->parms;
    if(!parms->gpu.tomo){
	gpu_dcell2cu(&curecon->opdr, simu->opdr);
    }
    toc("Before FitR");
    curcell *rhs=NULL;
    cumuv(&rhs, 0, &curecon->FR, curecon->opdr, 1);
    toc("FitR");
    if(gpu_pcg(&curecon->dmfit, (G_CGFUN)cumuv, &curecon->FL, NULL, NULL, rhs,
	       simu->parms->recon.warm_restart, parms->fit.maxit, curecon->cgstream)){
	curcellwrite(curecon->opdr, "fit_opdr_%d", simu->reconisim);
	curcellwrite(rhs, "fit_rhs_%d", simu->reconisim);
	curcellwrite(curecon->dmfit, "fit_dmfit_%d", simu->reconisim);
	dcellwrite(simu->gradlastol, "fit_gradlastol_%d", simu->reconisim);
	exit(1);
    }
    toc("FitL CG");
    gpu_curcell2d(&simu->dmfit_hi, curecon->dmfit, curecon->cgstream);
    cudaStreamSynchronize(curecon->cgstream);
    /*Don't free opdr. */
    curcellfree(rhs); rhs=NULL;
}
