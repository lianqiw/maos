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
#include "accphi.h"
#include "cucmat.h"
#include "kernel.h"
#include "cudata.h"
static int *cunembed=NULL;
static int *cupsfsize=NULL;
static float *cuwvls=NULL;    
pthread_mutex_t *evlmutex=NULL;
int *evlgpu=NULL;
static cudaStream_t *evlstream=NULL;
static cublasHandle_t *evlhandle=NULL;
static cufftHandle *evlplan=NULL;
/** 
    save aper_locs, aper_amp to GPU.
*/

__global__ static void calc_ptt_do( float *cc,
				    const float (*restrict loc)[2], 
				    const int nloc,
				    const float *restrict phi,
				    const float *restrict amp){
    __shared__ float ccb[4];/*for each block. */
    if(threadIdx.x<4){
	ccb[threadIdx.x]=0.f;
    }
    __syncthreads();
    float cci[4]={0.f,0.f,0.f,0.f};/*for each thread */
    int step=blockDim.x * gridDim.x; 
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float tmp=phi[i]*amp[i];
	cci[0]+=tmp;
	cci[1]+=tmp*loc[i][0];
	cci[2]+=tmp*loc[i][1];
	cci[3]+=tmp*phi[i];
    }
    /*Add results to shared value in each block. */
    atomicAdd(&ccb[0], cci[0]);
    atomicAdd(&ccb[1], cci[1]);
    atomicAdd(&ccb[2], cci[2]);
    atomicAdd(&ccb[3], cci[3]);
    __syncthreads();/*Wait until all threads in this block is done. */
    if(threadIdx.x<4){/*This is the first thread of a block. add block result to global. */
	atomicAdd(&cc[threadIdx.x], ccb[threadIdx.x]);
    }
}
__global__ static void calc_ngsmod_do( float *cc,
				       const float (*restrict loc)[2], 
				       const int nloc,
				       const float *restrict phi,
				       const float *restrict amp){
    int step=blockDim.x * gridDim.x; 
    float cci[7]={0,0,0,0,0,0,0};/*for each thread */
    __shared__ float ccb[7];/*for each block. */
    if(threadIdx.x<7){
	ccb[threadIdx.x]=0.f;
    }
    __syncthreads();
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	float tmp=phi[i]*amp[i];
	cci[0]+=tmp;
	cci[1]+=tmp*loc[i][0];
	cci[2]+=tmp*loc[i][1];
	cci[3]+=tmp*loc[i][0]*loc[i][0];
	cci[4]+=tmp*loc[i][1]*loc[i][1];
	cci[5]+=tmp*loc[i][0]*loc[i][1];
	cci[6]+=tmp*phi[i];
    }
    /*Add results to shared value in each block. */
    atomicAdd(&ccb[0], cci[0]);
    atomicAdd(&ccb[1], cci[1]);
    atomicAdd(&ccb[2], cci[2]);
    atomicAdd(&ccb[3], cci[3]);
    atomicAdd(&ccb[4], cci[4]);
    atomicAdd(&ccb[5], cci[5]);
    atomicAdd(&ccb[6], cci[6]);
    __syncthreads();/*Wait until all threads in this block is done. */
    if(threadIdx.x<7){/*This is the first thread of a block. add block result to global. */
	atomicAdd(&cc[threadIdx.x], ccb[threadIdx.x]);
    }
}
/*
  Let M be the modal matrix of pistion/tip/tilt. Calculate M'*diag(amp)*phi
  where amp is the amptliude weighting.  */
void calc_ptt(double *rmsout, double *coeffout, 
		  const double ipcc, const dmat *imcc,
		  const float (*restrict loc)[2], 
		  const int nloc,
		  const float *restrict phi,
		  const float *restrict amp,
		  cudaStream_t stream
		  ){
        /*sum with 16 blocks, each with 256 threads. */
    float *cc;
#if CUDAVER >=20
    float ccb[4];
#else
    float *ccb;
    cudaMallocHost(&ccb, 4*sizeof(float));
#endif
    cudaCalloc(cc, 4*sizeof(float), stream);
    calc_ptt_do<<<DIM(nloc, 256), 0, stream>>>(cc, loc, nloc, phi, amp);
    cudaMemcpyAsync(ccb, cc, 4*sizeof(float), cudaMemcpyDeviceToHost, stream);
    cudaFree(cc); cc=NULL;//this synchronizes
    double coeff[3], tot;
    coeff[0]=ccb[0]; coeff[1]=ccb[1]; coeff[2]=ccb[2]; tot=ccb[3];
    if(coeffout){
	dmulvec3(coeffout, imcc, coeff);
    }
    if(rmsout){
	double pis=ipcc*coeff[0]*coeff[0];/*piston mode variance */
	double ptt=dwdot3(coeff, imcc, coeff);/*p/t/t mode variance. */
	rmsout[0]=tot-pis;/*PR */
	rmsout[1]=ptt-pis;/*TT */
	rmsout[2]=tot-ptt;/*PTTR	 */
    }
#if CUDAVER <20
    cudaFreeHost(ccb);
#endif
}
void calc_ngsmod(double *pttr_out, double *pttrcoeff_out,
		 double *ngsmod_out, int nmod,
		 double MCC_fcp, double ht, double scale,
		 double thetax, double thetay,
		 const double ipcc, const dmat *imcc,
		 const float (*restrict loc)[2], 
		 const int nloc,
		 const float *restrict phi,
		 const float *restrict amp,
		 const PARMS_T *parms,
		 cudaStream_t stream){
    float *cc;
    double tot=0;
    cudaCalloc(cc, 7*sizeof(float), stream);
    if(nmod==2){/*single DM. */
	calc_ptt_do<<<DIM(nloc,256),0,stream>>>(cc, loc, nloc, phi, amp);
    }else if(nmod>=5){/*AHST mode */
	calc_ngsmod_do<<<DIM(nloc,256),0,stream>>>(cc, loc, nloc, phi, amp);
    }else{
	TO_IMPLEMENT;
    }
    float ccb[7];
    cudaMemcpyAsync(ccb, cc, 7*sizeof(float), cudaMemcpyDeviceToHost, stream);
    cudaFree(cc); 
    tot=ccb[nmod==2?3:6];
    
    double coeff[6];/*convert to double*/
    coeff[0]=ccb[0]; coeff[1]=ccb[1]; 
    coeff[2]=ccb[2]; coeff[3]=ccb[3];
    coeff[4]=ccb[4]; coeff[5]=ccb[5];
    
    if(pttrcoeff_out){/*p/t/t*/
	memset(pttrcoeff_out, 0, sizeof(double)*3);
	dmulvec(pttrcoeff_out, imcc, coeff, 1);
    }
    if(pttr_out){
	/*compute TT removed wavefront variance as a side product */
	double pis=ipcc*coeff[0]*coeff[0];
	double ptt=dwdot3(coeff, imcc, coeff);
	pttr_out[0]=tot-pis;/*PR */
	pttr_out[1]=ptt-pis;/*TT */
	pttr_out[2]=tot-ptt;/*PTTR */
    }
    /*don't use +=. need locking */
    ngsmod_out[0]=coeff[1];
    ngsmod_out[1]=coeff[2];
    const double scale1=1.-scale;
    if(nmod>=5){
	if(parms->sim.ahstfocus){
	    ngsmod_out[2]=(-2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
	}else{
	    ngsmod_out[2]=(scale1*(coeff[3]+coeff[4]-coeff[0]*MCC_fcp)
			   -2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
	}
	ngsmod_out[3]=(scale1*(coeff[3]-coeff[4])
		       -2*scale*ht*(thetax*coeff[1]-thetay*coeff[2]));
	ngsmod_out[4]=(scale1*(coeff[5])
		       -scale*ht*(thetay*coeff[1]+thetax*coeff[2]));
	if(nmod>5){
	    ngsmod_out[5]=(coeff[3]+coeff[4]-coeff[0]*MCC_fcp);
	}
    }
}
/**
   Convert NGS mode vector to aperture grid for science directions.  */
void gpu_ngsmod2science(curmat *opd, float (*restrict loc)[2],
			const NGSMOD_T *ngsmod, const double *mod, 
			double thetax, double thetay, 
			double alpha, cudaStream_t stream){
    if(ngsmod->nmod==2){
	curaddptt(opd, loc, 0, mod[0]*alpha, mod[1]*alpha, stream);
    }else{
	const float ht=ngsmod->ht;
	const float scale=ngsmod->scale;
	float focus;
	if(ngsmod->nmod>5){
	    focus=mod[5];
	    if(!ngsmod->ahstfocus){
		focus+=mod[2]*(1.-scale);
	    }
	}else{
	    focus=mod[2]*(1.-scale);
	}
	add_ngsmod_do<<<DIM(opd->nx*opd->ny, 256), 0, stream>>>
	    (opd->p, loc, opd->nx*opd->ny, 
	     mod[0], mod[1], mod[2], mod[3], mod[4], focus,
	     thetax, thetay, scale, ht, alpha);
    }
}
/**
   Initialize perfevl
*/
void gpu_perfevl_init(const PARMS_T *parms, APER_T *aper){
    const int nevl=parms->evl.nevl;
    const int nwvl=parms->evl.nwvl;
    evlgpu=(int*)calloc(nevl, sizeof(int));
    for(int ievl=0; ievl<nevl; ievl++){
	evlgpu[ievl]=gpu_next();
	if(NGPU>2 && evlgpu[ievl]==gpu_recon){
	    evlgpu[ievl]=gpu_next();
	}
    }
    /*The following lives in CPU. */
    if(parms->evl.psfmean || parms->evl.psfhist){
	cunembed =(int*)  calloc(nwvl, sizeof(int));
	cupsfsize=(int*)  calloc(nwvl, sizeof(int));
	cuwvls   =(float*)calloc(nwvl, sizeof(float));
    
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cunembed[iwvl]=(int)aper->nembed[iwvl];
	    cupsfsize[iwvl]=parms->evl.psfsize[iwvl];
	    cuwvls[iwvl]=parms->evl.wvl[iwvl];
	}
    }
    evlmutex=(pthread_mutex_t*)calloc(NGPU, sizeof(pthread_mutex_t));
    /*The following lives in GPU. */
    for(int im=0; im<NGPU; im++){
	pthread_mutex_init(&evlmutex[im], 0);
	gpu_set(im);
	cudata->plocs=new culoc_t(aper->locs);
	cp2gpu(&cudata->pamp, aper->amp);
	if(parms->evl.psfmean || parms->evl.psfhist){
	    cudata->embed    = (int**) calloc(nwvl, sizeof(int*));
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		cp2gpu(&cudata->embed[iwvl], aper->embed[iwvl], aper->locs->nloc, 1);
	    }
	}
    }/*for igpu */
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(evlgpu[ievl]);
	if(!cudata->plocs_dm){
	    cudata->plocs_dm=new culoc_t**[nevl];
	}
	cudata->plocs_dm[ievl]=new culoc_t*[parms->ndm];
	for(int idm=0; idm<parms->ndm; idm++){
	    loc_t *loc_dm;
	    if(aper->locs_dm && aper->locs_dm[ievl+idm*nevl]){
		loc_dm=aper->locs_dm[ievl+idm*nevl];
	    }else{
		loc_dm=aper->locs;
	    }
	    cudata->plocs_dm[ievl][idm]=new culoc_t(loc_dm);
	}
    }
    evlstream=(cudaStream_t*)calloc(nevl, sizeof(cudaStream_t));
    evlhandle=(cublasHandle_t*)calloc(nevl, sizeof(cublasHandle_t));
    if(parms->evl.psfmean || parms->evl.psfhist){
	evlplan  = (cufftHandle*)calloc(nwvl*nevl, sizeof(cufftHandle));
    }
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(evlgpu[ievl]);
	STREAM_NEW(evlstream[ievl]);
	cublasCreate(&evlhandle[ievl]);
	cublasSetStream(evlhandle[ievl], evlstream[ievl]);
	if(parms->evl.psfmean || parms->evl.psfhist){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(iwvl>0 && cunembed[iwvl]==cunembed[0]){
		    evlplan[iwvl+nwvl*ievl]=evlplan[0+nwvl*ievl];
		}else{
		    DO(cufftPlan2d(&evlplan[iwvl+nwvl*ievl],cunembed[iwvl],cunembed[iwvl],CUFFT_C2C));
		    DO(cufftSetStream(evlplan[iwvl+nwvl*ievl], evlstream[ievl]));
		}
	    }/*for iwvl */
	}
    }
    gpu_print_mem("perf init");
}
/*
  Initialize simulation data. Seed dependent.
*/
void gpu_perfevl_init_sim(const PARMS_T *parms, APER_T *aper){
    const int nevl=parms->evl.nevl;
    const int nwvl=parms->evl.nwvl;
    int nloc=aper->locs->nloc;
    if(!parms->gpu.evl){
	return;
    }
    /*first open loop ones are on every GPU.*/
    if(parms->evl.psfol){
	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    if(parms->evl.opdcov && parms->gpu.psf){ /*do OL opd cov*/
		initzero(&cudata->evlopdcovol, nloc, nloc);
		initzero(&cudata->evlopdmeanol, nloc, 1);
	    }
	    if(parms->evl.psfmean || parms->evl.psfhist){
		if(cudata->evlpsfol){
		    curcellzero(cudata->evlpsfol);
		}else{
		    cudata->evlpsfol=curcellnew(nwvl,1);
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			cudata->evlpsfol->p[iwvl]=curnew(cupsfsize[iwvl], cupsfsize[iwvl]);
		    }
		}
	    }
	}
    }

    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cudata->evlopd=curcellnew(nevl,1);
    }
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(evlgpu[ievl]);
	if(!cudata->evlopd->p[ievl]){
	    cudata->evlopd->p[ievl]=curnew(nloc, 1);
	}
    }
    
    if(parms->evl.opdcov && parms->gpu.psf && !parms->sim.evlol){
	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    if(!cudata->evlopdcov_ngsr){
		cudata->evlopdcov_ngsr=curcellnew(nevl, 1);
		cudata->evlopdmean_ngsr=curcellnew(nevl, 1);
		cudata->evlopdcov=curcellnew(nevl, 1);
		cudata->evlopdmean=curcellnew(nevl, 1);
	    }
	}
	for(int ievl=0; ievl<nevl; ievl++){
	    if(parms->evl.psf[ievl]==0){
		continue;
	    }
	    gpu_set(evlgpu[ievl]);
	    if(parms->evl.psfngsr[ievl]){
		initzero(&cudata->evlopdcov_ngsr->p[ievl], nloc,nloc);
		initzero(&cudata->evlopdmean_ngsr->p[ievl], nloc,1);
	    }
	    if(parms->evl.psfngsr[ievl]!=2){
		initzero(&cudata->evlopdcov->p[ievl],nloc,nloc);
		initzero(&cudata->evlopdmean->p[ievl],nloc,1);
	    }
	}
    }
	
    if((parms->evl.psfmean || parms->evl.psfhist) && !parms->sim.evlol){
	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    if(!cudata->evlwvf){/*temporary*/
		cudata->evlwvf=cuccellnew(nwvl, 1);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    if(!parms->evl.psfhist && iwvl>0 && cunembed[iwvl] == cunembed[iwvl-1]){
			cudata->evlwvf->p[iwvl]=cucref(cudata->evlwvf->p[iwvl-1]);
		    }else{
			cudata->evlwvf->p[iwvl]=cucnew(cunembed[iwvl], cunembed[iwvl]);
		    }
		}
	    }
	    if(!cudata->evlpsfcl){
		cudata->evlpsfcl = curcellnew(nwvl, parms->evl.nevl);
		cudata->evlpsfcl_ngsr = curcellnew(nwvl, parms->evl.nevl);
	    }
	}
	for(int ievl=0; ievl<nevl; ievl++){
	    if(parms->evl.psf[ievl]==0){
		continue;
	    }
	    gpu_set(evlgpu[ievl]);
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(parms->evl.psfngsr[ievl]){
		    initzero(&cudata->evlpsfcl_ngsr->p[iwvl+nwvl*ievl], 
			     cupsfsize[iwvl], cupsfsize[iwvl]);
		}
		if(parms->evl.psfngsr[ievl]!=2){
		    initzero(&cudata->evlpsfcl->p[iwvl+nwvl*ievl],
			     cupsfsize[iwvl], cupsfsize[iwvl]);
		}
	    }	
	}
    }
    CUDA_SYNC_DEVICE;
}
/**
   Add surface to surfevl;
*/
void gpu_evlsurf2gpu(APER_T *aper){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	if(aper->opdadd){
	    cp2gpu(&cudata->surfevl, aper->opdadd);
	}
    }
}


/**
   Compute complex PSF and return.
*/
static cuccell *psfcomp(curmat *iopdevl, int nwvl, int ievl, int nloc, cudaStream_t stream){
    LOCK(evlmutex[evlgpu[ievl]]);/*wvf is allocated per GPU.*/
    cuccell *psfs=cuccellnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	cucmat *wvf=cudata->evlwvf->p[iwvl];
	cuczero(wvf, stream);
	embed_wvf_do<<<DIM(iopdevl->nx,256),0,stream>>>
	    (wvf->p, iopdevl->p, cudata->pamp, cudata->embed[iwvl], nloc, cuwvls[iwvl]);
	CUFFT(evlplan[iwvl+nwvl*ievl], wvf->p, CUFFT_FORWARD);
	cucmat *psf=NULL;
	if(cupsfsize[iwvl]<cunembed[iwvl]){
	    psf=cucnew(cupsfsize[iwvl], cupsfsize[iwvl]);
	    corner2center_do<<<DIM2(psf->nx,psf->ny,16),0,stream>>>
		(psf->p, psf->nx, psf->ny, wvf->p, wvf->nx, wvf->ny);
	}else{
	    psf=cucref(wvf);
	    fftshift_do<<<DIM2(psf->nx,psf->ny,16),0,stream>>>
		(psf->p, psf->nx, psf->ny);
	}
	psfs->p[iwvl]=psf;
    }
    UNLOCK(evlmutex[evlgpu[ievl]]);
    return psfs;
}
/**
   Compute only PSF and add to result.
*/
static void psfcomp_r(curmat **psf, curmat *iopdevl, int nwvl, int ievl, int nloc, int atomic, cudaStream_t stream){
    LOCK(evlmutex[evlgpu[ievl]]);/*wvf is allocated per GPU.*/
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	cucmat *wvf=cudata->evlwvf->p[iwvl];
	cuczero(wvf, stream);
	embed_wvf_do<<<DIM(iopdevl->nx,256),0,stream>>>
	    (wvf->p, iopdevl->p, cudata->pamp, cudata->embed[iwvl], nloc, cuwvls[iwvl]);
	CUFFT(evlplan[iwvl+nwvl*ievl], wvf->p, CUFFT_FORWARD);
	if(!psf[iwvl]) psf[iwvl]=curnew(cupsfsize[iwvl], cupsfsize[iwvl]);
	if(atomic){
	    corner2center_abs2_atomic_do<<<DIM2((psf[iwvl])->nx,(psf[iwvl])->ny,16),0,stream>>>
		((psf[iwvl])->p, (psf[iwvl])->nx, (psf[iwvl])->ny, wvf->p, wvf->nx, wvf->ny);
	}else{
	    corner2center_abs2_do<<<DIM2((psf[iwvl])->nx,(psf[iwvl])->ny,16),0,stream>>>
		((psf[iwvl])->p, (psf[iwvl])->nx, (psf[iwvl])->ny, wvf->p, wvf->nx, wvf->ny);
	}
    }
    UNLOCK(evlmutex[evlgpu[ievl]]);
}
#define PERFEVL_WFE(pclep, pclmp, cleNGSmp)				\
    if(parms->recon.split){						\
	if(parms->ndm<=2){						\
	    PDMAT(cleNGSmp->p[ievl], pcleNGSmp);			\
	    if(nmod==3){						\
		calc_ngsmod(pclep[isim], pclmp[isim], pcleNGSmp[isim],recon->ngsmod->nmod, \
			    recon->ngsmod->aper_fcp, recon->ngsmod->ht, \
			    recon->ngsmod->scale, thetax, thetay,	\
			    aper->ipcc, aper->imcc,			\
			    cudata->plocs->p, nloc, iopdevl->p, cudata->pamp, parms, stream); \
	    }else{							\
		calc_ngsmod(NULL, NULL, pcleNGSmp[isim],recon->ngsmod->nmod, \
			    recon->ngsmod->aper_fcp, recon->ngsmod->ht, \
			    recon->ngsmod->scale, thetax, thetay,	\
			    aper->ipcc, aper->imcc,			\
			    cudata->plocs->p, nloc, iopdevl->p, cudata->pamp, parms, stream); \
		TO_IMPLEMENT;/*mode decomposition. */			\
	    }								\
	}								\
    }else{								\
	if(nmod==3){							\
	    calc_ptt(pclep[isim], pclmp[isim], aper->ipcc, aper->imcc,	\
		     cudata->plocs->p, nloc, iopdevl->p, cudata->pamp, stream); \
	}else{								\
	    TO_IMPLEMENT;						\
	}								\
    }

/**
   Performance evaluation. Designed to replace perfevl_ievl in maos/perfevl.c
*/
void gpu_perfevl(thread_t *info){
    const int ievl=info->start;
    /*lock the mutex because iopdevl, evlwvf is allocated per GPU.*/
    gpu_set(evlgpu[ievl]);
    SIM_T *simu=(SIM_T*)info->data;
    assert(info->end==info->start+1);/*only one evl. */
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int isim=simu->isim;
    const int nmod=parms->evl.nmod;
    const int imoao=parms->evl.moao;
    const double dt=simu->dt;
    const int do_psf_cov=(parms->evl.psfmean || parms->evl.psfhist || parms->evl.opdcov) 
	&& isim>=parms->evl.psfisim && parms->evl.psf[ievl]!=0;
    const int save_evlopd=parms->save.evlopd>0 && ((isim+1)%parms->save.evlopd)==0;
    const int nloc=aper->locs->nloc;
    const int nwvl=parms->evl.nwvl;
    const double thetax=parms->evl.thetax[ievl];
    const double thetay=parms->evl.thetay[ievl];
    /*Setup pointers for easy usage */
    PDMAT(simu->olmp->p[ievl],polmp);/*OL mode for each dir */
    PDMAT(simu->olep->p[ievl],polep);/*OL error for each dir */
    PDMAT(simu->clmp->p[ievl],pclmp);
    PDMAT(simu->clep->p[ievl],pclep);
    cudaStream_t stream=evlstream[ievl];
    cublasHandle_t handle=evlhandle[ievl];
    curmat *iopdevl=cudata->evlopd->p[ievl];
    curzero(iopdevl, stream);
    /* iopdevl must be in device memory. 6 times slower if in host memory.*/
    if(cudata->surfevl && cudata->surfevl->p[ievl]){
	curcp(&iopdevl, cudata->surfevl->p[ievl], stream);
    }else{
	curset(iopdevl, 0, stream);
    }

    if(parms->sim.idealevl){
	gpu_dm2loc(iopdevl->p, cudata->plocs_dm[ievl], cudata->dmproj, cudata->ndm,
		   parms->evl.hs[ievl], thetax, thetay,
		   0,0,1, stream);
    }else if(simu->atm && !parms->sim.wfsalias){
	gpu_atm2loc(iopdevl->p, cudata->plocs, parms->evl.hs[ievl], thetax, thetay, 
		    0,0,isim*dt, 1, stream);
    }
    
    if(simu->telws){/*Wind shake */
	float tt=simu->telws->p[isim];
	float angle=simu->winddir?simu->winddir->p[0]:0;
	curaddptt(iopdevl, cudata->plocs->p, 0, tt*cosf(angle), tt*sinf(angle), stream);
    }
    if(save_evlopd){
	cellarr_cur(simu->save->evlopdol[ievl], isim, iopdevl, stream);
    }
    if(parms->plot.run>1){
	dmat *tmp=NULL;
	cp2cpu(&tmp, iopdevl, stream);
	drawopdamp("OL", aper->locs, tmp->p, aper->amp1->p, NULL,
		   "Science Open Loop OPD", "x (m)", "y (m)", "OL %d", ievl);
	dfree(tmp);
    }
    PERFEVL_WFE(polep, polmp, simu->oleNGSmp);
    if((parms->evl.psfmean  || parms->evl.opdcov)
       && isim>=parms->evl.psfisim 
       &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
	  ||(parms->evl.psfol==2 && parms->evl.psf[ievl]))){
	/*calculate Openloop PSF. we also test psfisim to synchronize with psfcl.*/
	curmat *opdcopy=NULL;
	if(parms->evl.psfpttr[ievl]){
	    curcp(&opdcopy, iopdevl, stream);
	    curaddptt(opdcopy, cudata->plocs->p, -polmp[isim][0], -polmp[isim][1], -polmp[isim][2], stream);
	}else if(parms->evl.opdcov){/* remove piston*/
	    curcp(&opdcopy, iopdevl, stream);
	    curadd(opdcopy, -polmp[isim][0], stream);
	}else{
	    opdcopy=iopdevl;
	}
	if(parms->evl.opdcov){
	    if(parms->gpu.psf){
		curmm(&cudata->evlopdcovol, 1, opdcopy, opdcopy, "nt", 1, handle);
		curadd(&cudata->evlopdmeanol, 1, opdcopy, 1, stream);
	    }else{
		dmat *tmp=NULL;
		cp2cpu(&tmp, opdcopy, stream);
		dmm(&simu->evlopdcovol, tmp, tmp, "nt", 1);
		dadd(&simu->evlopdmeanol, 1, tmp, 1);
		dfree(tmp);
	    }
	}
	if(parms->evl.psfmean){
	    psfcomp_r(cudata->evlpsfol->p, opdcopy, nwvl, ievl, nloc, parms->evl.psfol==2?1:0, stream);
	    if(opdcopy!=iopdevl){
		curfree(opdcopy);
	    }
	    if(!parms->gpu.psf){ /*need to move psf from GPU to CPU for accumulation.*/
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    add2cpu(&simu->evlpsfolmean->p[iwvl], cudata->evlpsfol->p[iwvl], stream);
		    curfree(cudata->evlpsfol->p[iwvl]); cudata->evlpsfol->p[iwvl]=NULL;
		}
	    }
	}
    }
    if(parms->sim.evlol) goto end;
    
    if(parms->evl.tomo){
	TO_IMPLEMENT;
    }else{
	gpu_dm2loc(iopdevl->p, cudata->plocs_dm[ievl], cudata->dmreal, cudata->ndm, 
		   parms->evl.hs[ievl], thetax, thetay,
		   0,0,-1, stream);
	if(imoao!=-1){
	    gpu_dm2loc(iopdevl->p, &cudata->plocs, cudata->dm_evl[ievl], 1,
		       INFINITY, 0, 0, 0, 0, -1, stream);
	}
    }
    if(save_evlopd){
	cellarr_cur(simu->save->evlopdcl[ievl], isim, iopdevl, stream);
    }
    if(parms->plot.run>1){
	dmat *tmp=NULL;
	cp2cpu(&tmp, iopdevl, stream);
	drawopdamp("CL", aper->locs,tmp->p , aper->amp1->p, NULL,
		   "Science Closed loop OPD", "x (m)", "y (m)", "CL %d", ievl);
	dfree(tmp);
    }
    PERFEVL_WFE(pclep, pclmp, simu->cleNGSmp);
    if(do_psf_cov){
	if(parms->evl.psfngsr[ievl]!=2){/*also do normal one.*/
	    if(parms->evl.psfpttr[ievl]){
		curaddptt(iopdevl, cudata->plocs->p, -pclmp[isim][0], -pclmp[isim][1], -pclmp[isim][2], stream);
	    }else{
		curadd(iopdevl, -pclmp[isim][0], stream);
	    }
	    if(parms->evl.opdcov){
		if(parms->gpu.psf){
		    curmm(&cudata->evlopdcov->p[ievl], 1, iopdevl, iopdevl, "nt", 1, handle);
		    curadd(&cudata->evlopdmean->p[ievl], 1, iopdevl, 1, stream);
		}else{
		    dmat *tmp=NULL;
		    cp2cpu(&tmp, iopdevl, stream);
		    dmm(&simu->evlopdcov->p[ievl], tmp, tmp, "nt", 1);
		    dadd(&simu->evlopdmean->p[ievl], 1, tmp, 1);
		    dfree(tmp);
		}
	    }/*opdcov */
	    if(parms->evl.psfhist || parms->evl.psfmean){
		if(parms->evl.psfhist){
		    /*Compute complex. */
		    cuccell *psfs=psfcomp(iopdevl, nwvl, ievl, nloc, stream);
		    cellarr_cuccell(simu->save->evlpsfhist[ievl], isim, psfs, stream);
		    if(parms->evl.psfmean){
			for(int iwvl=0; iwvl<nwvl; iwvl++){
			    curaddcabs2(cudata->evlpsfcl->p+iwvl+nwvl*ievl, 1, psfs->p[iwvl], 1, stream);
			}
		    }
		    curcellfree(psfs);
		}else if(parms->evl.psfmean){
		    psfcomp_r(cudata->evlpsfcl->p+nwvl*ievl, iopdevl, nwvl, ievl, nloc, 0, stream);
		}
		if(!parms->gpu.psf){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			add2cpu(&simu->evlpsfmean->p[iwvl+ievl*nwvl], cudata->evlpsfcl->p[iwvl+ievl*nwvl], stream);
			curfree(cudata->evlpsfcl->p[iwvl+ievl*nwvl]); cudata->evlpsfcl->p[iwvl+ievl*nwvl]=NULL;
		    }
		}
	    }
	}
    }
 end:
    CUDA_SYNC_STREAM;
}
/**
   Compute the PSF or OPDCOV for NGS mode removed opd.
*/
void gpu_perfevl_ngsr(SIM_T *simu, double *cleNGSm){
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const int nloc=aper->locs->nloc;
    const int nwvl=parms->evl.nwvl;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(parms->evl.psfngsr[ievl]==0){
	    continue;
	}
	warning("Compare with CPU code to verify accuracy. Need to verify focus mode\n");
	gpu_set(evlgpu[ievl]);
	curmat *iopdevl=cudata->evlopd->p[ievl];
	cudaStream_t stream=evlstream[ievl];
	cublasHandle_t handle=evlhandle[ievl];
	gpu_ngsmod2science(iopdevl, cudata->plocs->p, simu->recon->ngsmod, cleNGSm, 
			   parms->evl.thetax[ievl], parms->evl.thetay[ievl],
			   -1, stream);
	if(parms->evl.psfpttr[ievl]){
	    double ptt[3];
	    calc_ptt(NULL, ptt,  aper->ipcc, aper->imcc,
		     cudata->plocs->p, nloc, iopdevl->p, cudata->pamp, stream);
	    curaddptt(iopdevl, cudata->plocs->p, -ptt[0], -ptt[1], -ptt[2], stream);
	}
	if(parms->evl.opdcov){
	    if(parms->gpu.psf){
		curmm(&cudata->evlopdcov_ngsr->p[ievl], 1, iopdevl, iopdevl, "nt", 1, handle);
		curadd(&cudata->evlopdmean_ngsr->p[ievl], 1, iopdevl, 1, stream);
	    }else{
		dmat *tmp=NULL;
		cp2cpu(&tmp, iopdevl, stream);
		dmm(&simu->evlopdcov_ngsr->p[ievl], tmp, tmp, "nt", 1);
		dadd(&simu->evlopdmean_ngsr->p[ievl], 1, tmp, 1);
		dfree(tmp);
	    }
	}/*opdcov */
	if(parms->evl.psfhist||parms->evl.psfmean){
	    if(parms->evl.psfhist){
		/*Compute complex. */
		cuccell *psfs=psfcomp(iopdevl, nwvl, ievl, nloc, stream);
		cellarr_cuccell(simu->save->evlpsfhist_ngsr[ievl], simu->isim, psfs, stream);
		if(parms->evl.psfmean){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			curaddcabs2(cudata->evlpsfcl_ngsr->p+iwvl+nwvl*ievl, 1, psfs->p[iwvl], 1, stream);
		    }
		}
	    }else if(parms->evl.psfmean){
		psfcomp_r(cudata->evlpsfcl_ngsr->p+nwvl*ievl, iopdevl, nwvl, ievl, nloc, 0, stream);
	    }
	    if(!parms->gpu.psf){
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    add2cpu(&simu->evlpsfmean_ngsr->p[iwvl+ievl*nwvl], cudata->evlpsfcl_ngsr->p[iwvl+ievl*nwvl], stream);
		    curfree(cudata->evlpsfcl_ngsr->p[iwvl+ievl*nwvl]); cudata->evlpsfcl_ngsr->p[iwvl+ievl*nwvl]=NULL;
		}
	    }
	}
	CUDA_SYNC_STREAM;
   }
}
void gpu_perfevl_save(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(!parms->evl.nevl) return;
    const int isim=simu->isim;
    if(parms->evl.psfmean && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.psfmean)){
	info2("Step %d: Output PSF\n", isim);
	const int nwvl=parms->evl.nwvl;
	if(cudata->evlpsfol){
	    /*copy the PSF accumulated in all the GPUs to CPU.*/
	    scell *temp=scellnew(nwvl, 1);
	    scell *temp2=scellnew(nwvl, 1);
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    if(parms->evl.psfol==2){
		scale=scale/parms->evl.npsf;
	    }
	    for(int im=0; im<NGPU; im++){
		gpu_set(im);
		cp2cpu(&temp2, cudata->evlpsfol, 0);
		cudaStreamSynchronize(0);
		scelladd(&temp, 1, temp2, scale);
	    }
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(!temp || !temp->p[iwvl]) continue;
		temp->p[iwvl]->header=evl_header(simu->parms, simu->aper, -1, iwvl);
		cellarr_smat(simu->save->evlpsfolmean, isim, temp->p[iwvl]);
		free(temp->p[iwvl]->header); temp->p[iwvl]->header=NULL;
	    }
	    scellfree(temp);
	    scellfree(temp2);
	}
	if(cudata->evlpsfcl){
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->evl.psf[ievl] || parms->evl.psfngsr[ievl]==2) continue;
		gpu_set(evlgpu[ievl]);
		cudaStream_t stream=evlstream[ievl];
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    curmat *pp=cudata->evlpsfcl->p[iwvl+nwvl*ievl];
		    curscale(pp, scale, stream);
		    if(!pp->header){
			pp->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		    }
		    cellarr_cur(simu->save->evlpsfmean[ievl], isim, pp, stream);
		    curscale(pp, 1.f/scale, stream);
		}
	    }
	}
	if(cudata->evlpsfcl_ngsr){
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->evl.psf[ievl] || !parms->evl.psfngsr[ievl]) continue;
		gpu_set(evlgpu[ievl]);
		cudaStream_t stream=evlstream[ievl];
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    curmat *pp=cudata->evlpsfcl_ngsr->p[iwvl+nwvl*ievl];
		    curscale(pp, scale, stream);
		    if(!pp->header){
			pp->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		    }
		    cellarr_cur(simu->save->evlpsfmean_ngsr[ievl], isim, pp, stream);
		    curscale(pp, 1.f/scale, stream);
		}
	    }
	}
    }
    if(parms->evl.opdcov && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.opdcov)){
	info2("Step %d: Output opdcov\n", isim);
	double scale=1./(isim+1-parms->evl.psfisim);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!parms->evl.psf[ievl]|| parms->evl.psfngsr[ievl]==2) continue;
	    gpu_set(evlgpu[ievl]);
	    cudaStream_t stream=evlstream[ievl];
	    curmat *pp;
	    {
		pp=cudata->evlopdcov->p[ievl];
		curscale(pp, scale, stream);
		cellarr_cur(simu->save->evlopdcov[ievl], isim, pp, stream);
		curscale(pp, 1./scale, stream);
	    }
	    {
		pp=cudata->evlopdmean->p[ievl];
		curscale(pp, scale, stream);
		cellarr_cur(simu->save->evlopdmean[ievl], isim, pp, stream);
		curscale(pp, 1./scale, stream);
	    }
	}
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!parms->evl.psf[ievl]|| !parms->evl.psfngsr[ievl]) continue;
	    gpu_set(evlgpu[ievl]);
	    cudaStream_t stream=evlstream[ievl];
	    curmat *pp;
	    {
		pp=cudata->evlopdcov_ngsr->p[ievl];
		curscale(pp, scale, stream);
		cellarr_cur(simu->save->evlopdcov_ngsr[ievl], isim, pp, stream);
		curscale(pp, 1./scale, stream);
	    }
	    {
		pp=cudata->evlopdmean_ngsr->p[ievl];
		curscale(pp, scale, stream);
		cellarr_cur(simu->save->evlopdmean_ngsr[ievl], isim, pp, stream);
		curscale(pp, 1./scale, stream);
	    }
	}
	if(parms->evl.psfol){
	    if(parms->evl.psfol==2){
		scale=scale/parms->evl.npsf;
	    }
	    {
		smat *temp=NULL;
		smat *temp2=NULL;
		for(int im=0; im<NGPU; im++){
		    gpu_set(im);
		    cp2cpu(&temp2, cudata->evlopdcovol, 0);
		    cudaStreamSynchronize(0);
		    sadd(&temp, 1, temp2, scale);
		}
		cellarr_smat(simu->save->evlopdcovol, isim, temp);
		sfree(temp);
		sfree(temp2);
	    }
	    {
		smat *temp=NULL;
		smat *temp2=NULL;
		for(int im=0; im<NGPU; im++){
		    gpu_set(im);
		    cp2cpu(&temp2, cudata->evlopdmeanol, 0);
		    cudaStreamSynchronize(0);
		    sadd(&temp, 1, temp2, scale);
		}
		cellarr_smat(simu->save->evlopdmeanol, isim, temp);
		sfree(temp);
		sfree(temp2);
	    }
	}
    }
}
