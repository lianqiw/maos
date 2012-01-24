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
#include "accphi.h"
#include "cucmat.h"
#include "kernel.h"
static int *cunembed=NULL;
static int *cupsfsize=NULL;
static float *cuwvls=NULL;    
static int evlfft_lock=0;
static int *evlgpu=NULL;
static cudaStream_t *evlstream=NULL;
static cublasHandle_t *evlhandle=NULL;
static cufftHandle *evlplan=NULL;
/** 
    save aper_locs, aper_amp to GPU.
*/
static void gpu_plocs2gpu(loc_t *loc, dmat *amp){
    if(cudata->plocs) error("plocs is already Copied.\n");
    int nloc=loc->nloc;
    gpu_loc2dev(&cudata->plocs, loc);
    gpu_dbl2dev(&cudata->pamp, amp->p, nloc);
}

/**
   Convert NGS mode vector to aperture grid for science directions.  */
void gpu_ngsmod2science(curmat *opd, const PARMS_T *parms,
			const RECON_T *recon, const APER_T *aper,
			const double *mod, int ievl, double alpha, cudaStream_t stream){
    if(recon->ngsmod->nmod==2){
	curaddptt(opd, cudata->plocs, 0, mod[0]*alpha, mod[1]*alpha, stream);
    }else{
	const float MCC_fcp=recon->ngsmod->aper_fcp;
	const float ht=recon->ngsmod->ht;
	const float scale=recon->ngsmod->scale;
	const float thetax=parms->evl.thetax[ievl]; 
	const float thetay=parms->evl.thetay[ievl];
	add_ngsmod_do<<<DIM(opd->nx*opd->ny, 256), 0, stream>>>
	    (opd->p, cudata->plocs, opd->nx*opd->ny, 
	     mod[0], mod[1], mod[2], mod[3], mod[4],
	     thetax, thetay, scale, ht, MCC_fcp, alpha);
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
    }
    /*The following lives in CPU. */
    if(parms->evl.psfmean || parms->evl.psfhist){
	cunembed =(int*)  calloc(nwvl, sizeof(int));
	cupsfsize=(int*)  calloc(nwvl, sizeof(int));
	cuwvls   =(float*)calloc(nwvl, sizeof(float));
    
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cunembed[iwvl]=(int)aper->nembed[iwvl];
	    if(cunembed[iwvl]>1024 && nevl>2){
		evlfft_lock=1;/**FFT will take a lot of memory, so do one direction at a time */
	    }
	    cupsfsize[iwvl]=parms->evl.psfsize[iwvl];
	    cuwvls[iwvl]=parms->evl.wvl[iwvl];
	}
    }
    /*The following lives in GPU. */
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	gpu_plocs2gpu(aper->locs, aper->amp);
	if(parms->evl.psfmean || parms->evl.psfhist){
	    cudata->embed    = (int**) calloc(nwvl, sizeof(int*));
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		gpu_long2dev(&cudata->embed[iwvl], aper->embed[iwvl], aper->locs->nloc);
	    }
	}
    }/*for igpu */
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
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	if(parms->evl.opdcov){
	    curcellfree(cudata->evlopdcov);
	    curcellfree(cudata->evlopdcov_ngsr);
	    cudata->evlopdcov     =curcellnew(nevl, 1);
	    cudata->evlopdcov_ngsr=curcellnew(nevl, 1);
	}
	curcellfree(cudata->evlopd);
	cudata->evlopd=curcellnew(nevl,1);
	if(parms->evl.psfmean || parms->evl.psfhist){
	    curcellfree(cudata->evlpsfcl);
	    curcellfree(cudata->evlpsfcl_ngsr);
	    curcellfree(cudata->evlpsfol);
	    cudata->evlpsfcl = curcellnew(nwvl, parms->evl.nevl);
	    cudata->evlpsfcl_ngsr = curcellnew(nwvl, parms->evl.nevl);
	    cudata->evlpsfol = curcellnew(nwvl, 1);
	}
	CUDA_SYNC_DEVICE;
    }
}
/**
   Add surface to surfevl;
*/
void gpu_evlsurf2gpu(APER_T *aper){
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	if(aper->opdadd){
	    gpu_dcell2cu(&cudata->surfevl, aper->opdadd);
	    dcellfree(aper->opdadd);
	}
    }
}

__global__ static void evl_embed_wvf_do(fcomplex *restrict wvf, 
					const float *restrict opd, const float *restrict amp, 
					const int *embed, const int nloc, const float wvl){
    const float pi2l=2.f*M_PI/wvl;
    for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nloc; ix+=blockDim.x*gridDim.x){
	float s,c;
	sincosf(pi2l*opd[ix], &s, &c);
	wvf[embed[ix]]=make_cuComplex(amp[ix]*c, amp[ix]*s);
    }
}

/**
   Embed or crop an array to another array. Preserve corner.
*/
__global__ static void evl_corner2center_do(fcomplex *restrict out, int noutx,  int nouty,
					    const fcomplex *restrict in, int ninx, int niny){
    int nx,ny;
    ny=MIN(niny, nouty)>>1;
    nx=MIN(ninx, noutx)>>1;
    int noutx2=noutx>>1;
    int nouty2=nouty>>1;
    for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
	for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
	    out[(iy+nouty2 )*noutx+(ix+noutx2)  ] = in[iy*ninx+ix];
	    out[(iy+nouty2 )*noutx+(noutx2-1-ix)] = in[iy*ninx+(ninx-1-ix)];
	    out[(nouty2-1-iy)*noutx+(noutx2-1-ix)] = in[(niny-1-iy)*ninx+(ninx-1-ix)];
	    out[(nouty2-1-iy)*noutx+(ix+noutx2)  ] = in[(niny-1-iy)*ninx+(ix)];
	}
    }
}

/**
   Embed or crop an array to another array. Preserve corner.
*/
__global__ static void evl_corner2center_abs2_do(float *restrict out, int noutx,  int nouty,
					    const fcomplex *restrict in, int ninx, int niny){
    int nx,ny;
    ny=MIN(niny, nouty)>>1;
    nx=MIN(ninx, noutx)>>1;
    int noutx2=noutx>>1;
    int nouty2=nouty>>1;
    for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
	for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
	    out[(iy+nouty2)*noutx+(ix+noutx2)]+=     CABS2(in[iy*ninx+ix]);
	    out[(iy+nouty2)*noutx+(noutx2-1-ix)]+=   CABS2(in[iy*ninx+(ninx-1-ix)]);
	    out[(nouty2-1-iy)*noutx+(noutx2-1-ix)]+= CABS2(in[(niny-1-iy)*ninx+(ninx-1-ix)]);
	    out[(nouty2-1-iy)*noutx+(ix+noutx2)]+=   CABS2(in[(niny-1-iy)*ninx+(ix)]);
	}
    }
}
/**
   Embed or crop an array to another array. Preserve corner.
*/
__global__ static void evl_corner2center_abs2_atomic_do(float *restrict out, int noutx,  int nouty,
					    const fcomplex *restrict in, int ninx, int niny){
    int nx,ny;
    ny=MIN(niny, nouty)>>1;
    nx=MIN(ninx, noutx)>>1;
    int noutx2=noutx>>1;
    int nouty2=nouty>>1;
    for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
	for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
	    atomicAdd(&out[(iy+nouty2)*noutx+(ix+noutx2)],     CABS2(in[iy*ninx+ix]));
	    atomicAdd(&out[(iy+nouty2)*noutx+(noutx2-1-ix)],   CABS2(in[iy*ninx+(ninx-1-ix)]));
	    atomicAdd(&out[(nouty2-1-iy)*noutx+(noutx2-1-ix)], CABS2(in[(niny-1-iy)*ninx+(ninx-1-ix)]));
	    atomicAdd(&out[(nouty2-1-iy)*noutx+(ix+noutx2)],   CABS2(in[(niny-1-iy)*ninx+(ix)]));
	}
    }
}
/**
   FFT Shift.
*/
__global__ static void evl_fftshift_do(fcomplex *wvf, const int nx, const int ny){
    int nx2=nx>>1;
    int ny2=ny>>1;
    for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny2; iy+=blockDim.y*gridDim.y){
	for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx2; ix+=blockDim.x*gridDim.x){
	    fcomplex tmp;
	    tmp=wvf[ix+iy*nx];
	    wvf[ix+iy*nx]=wvf[(ix+nx2)+(iy+ny2)*nx];
	    wvf[(ix+nx2)+(iy+ny2)*nx]=tmp;
	    tmp=wvf[ix+(iy+ny2)*nx];
	    wvf[ix+(iy+ny2)*nx]=wvf[(ix+nx2)+iy*nx];
	    wvf[(ix+nx2)+iy*nx]=tmp;
	}
    }
}
/**
   Compute complex PSF and return.
*/
static cuccell *psfcomp(curmat *iopdevl, int nwvl, int ievl, int nloc, cudaStream_t stream){
    static pthread_mutex_t evl_mutex=PTHREAD_MUTEX_INITIALIZER;
    if(evlfft_lock) LOCK(evl_mutex);
    cuccell *psfs=cuccellnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	cucmat *wvf=cucnew(cunembed[iwvl], cunembed[iwvl],stream);
	evl_embed_wvf_do<<<DIM(iopdevl->nx,256),0,stream>>>
	    (wvf->p, iopdevl->p, cudata->pamp, cudata->embed[iwvl], nloc, cuwvls[iwvl]);
	CUFFT(evlplan[iwvl+nwvl*ievl], wvf->p, CUFFT_FORWARD);
	cucmat *psf=NULL;
	if(cupsfsize[iwvl]<cunembed[iwvl]){
	    psf=cucnew(cupsfsize[iwvl], cupsfsize[iwvl]);
	    evl_corner2center_do<<<DIM2(psf->nx,psf->ny,16),0,stream>>>
		(psf->p, psf->nx, psf->ny, wvf->p, wvf->nx, wvf->ny);
	}else{
	    psf=wvf;
	    evl_fftshift_do<<<DIM2(psf->nx,psf->ny,16),0,stream>>>
	    (psf->p, psf->nx, psf->ny);
	}
	if(psf!=wvf) cucfree(wvf);
	psfs->p[iwvl]=psf;
    }
    if(evlfft_lock) UNLOCK(evl_mutex);
    return psfs;
}
/**
   Compute only PSF and add to result.
*/
static void psfcomp_r(curmat **psf, curmat *iopdevl, int nwvl, int ievl, int nloc, int atomic, cudaStream_t stream){
    static pthread_mutex_t evl_mutex=PTHREAD_MUTEX_INITIALIZER;
    if(evlfft_lock) LOCK(evl_mutex);
    cucmat *wvf=NULL;
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	if(cunembed[iwvl]!=cunembed[0]){
	    CUDA_SYNC_STREAM;
	    cucfree(wvf); wvf=NULL;
	}
	if(!wvf){
	    wvf=cucnew(cunembed[iwvl], cunembed[iwvl], stream);
	}else{
	    cuczero(wvf, stream);
	}
	evl_embed_wvf_do<<<DIM(iopdevl->nx,256),0,stream>>>
	    (wvf->p, iopdevl->p, cudata->pamp, cudata->embed[iwvl], nloc, cuwvls[iwvl]);
	CUFFT(evlplan[iwvl+nwvl*ievl], wvf->p, CUFFT_FORWARD);
	if(!psf[iwvl]) psf[iwvl]=curnew(cupsfsize[iwvl], cupsfsize[iwvl], stream);
	if(atomic){
	    evl_corner2center_abs2_atomic_do<<<DIM2((psf[iwvl])->nx,(psf[iwvl])->ny,16),0,stream>>>
		((psf[iwvl])->p, (psf[iwvl])->nx, (psf[iwvl])->ny, wvf->p, wvf->nx, wvf->ny);
	}else{
	    evl_corner2center_abs2_do<<<DIM2((psf[iwvl])->nx,(psf[iwvl])->ny,16),0,stream>>>
		((psf[iwvl])->p, (psf[iwvl])->nx, (psf[iwvl])->ny, wvf->p, wvf->nx, wvf->ny);
	}
    }
    CUDA_SYNC_STREAM;
    cucfree(wvf);
    if(evlfft_lock) UNLOCK(evl_mutex);
}


/**
   Performance evaluation. Designed to replace perfevl_ievl in maos/perfevl.c
*/
void gpu_perfevl(thread_t *info){
    const int ievl=info->start;
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

    curmat *iopdevl=curnew(aper->locs->nloc, 1);
    cudaStream_t stream=evlstream[ievl];
    cublasHandle_t handle=evlhandle[ievl];
    /* iopdevl must be in device memory. 6 times slower if in host memory.*/
    if(cudata->surfevl && cudata->surfevl->p[ievl]){
	curcp(&iopdevl, cudata->surfevl->p[ievl], stream);
    }else{
	curset(iopdevl, 0, stream);
    }

    if(parms->sim.idealevl){
	gpu_dm2loc(iopdevl->p, cudata->plocs, nloc, cudata->dmproj, parms->evl.hs[ievl], thetax, thetay,
		   parms->evl.misreg[0], parms->evl.misreg[1], 1, stream);
    }else if(simu->atm && !parms->sim.wfsalias){
	gpu_atm2loc(iopdevl->p, cudata->plocs, nloc, parms->evl.hs[ievl], thetax, thetay, 
		    parms->evl.misreg[0], parms->evl.misreg[1], isim*dt, 1, stream);
    }
    
    if(simu->telws){/*Wind shake */
	float tt=simu->telws->p[isim];
	float angle=simu->winddir?simu->winddir->p[0]:0;
	curaddptt(iopdevl, cudata->plocs, 0, tt*cosf(angle), tt*sinf(angle), stream);
    }
    if(save_evlopd){
	cellarr_cur(simu->save->evlopdol[ievl], iopdevl, stream);
    }
    if(parms->plot.run){
	/*TO_IMPLEMENT; */
    }
    if(nmod==3){
	gpu_calc_ptt(polep[isim], polmp[isim], aper->ipcc, aper->imcc,
		     cudata->plocs, nloc, iopdevl->p, cudata->pamp, stream);
    }else{
	TO_IMPLEMENT;
    }
    if(parms->evl.psfmean &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
			     ||(parms->evl.psfol==2 && parms->evl.psf[ievl]))){
	/*calculate Openloop PSF. */
	curmat *opdcopy=NULL;
	if(parms->evl.psfpttr[ievl]){
	    curcp(&opdcopy, iopdevl, stream);
	    curaddptt(opdcopy, cudata->plocs, -polmp[isim][0], -polmp[isim][1], -polmp[isim][2], stream);
	}else{
	    opdcopy=iopdevl;
	}
	psfcomp_r(cudata->evlpsfol->p, opdcopy, nwvl, ievl, nloc, parms->evl.psfol==2?1:0, stream);
	if(opdcopy!=iopdevl){
	    CUDA_SYNC_STREAM;
	    curfree(opdcopy);
	}
    }
    
    if(parms->sim.evlol) goto end;
    
    if(parms->evl.tomo){
	TO_IMPLEMENT;
    }else{
	gpu_dm2loc(iopdevl->p, cudata->plocs, nloc, cudata->dmreal, parms->evl.hs[ievl], thetax, thetay,
	       parms->evl.misreg[0], parms->evl.misreg[1], -1, stream);
	if(imoao>-1){
	    TO_IMPLEMENT;
	}
    }
    if(save_evlopd){
	cellarr_cur(simu->save->evlopdcl[ievl], iopdevl, stream);
    }
    if(parms->plot.run){
	/*TO_IMPLEMENT; */
    }
    if(parms->recon.split){
	if(parms->ndm<=2){
	    PDMAT(simu->cleNGSmp->p[ievl], pcleNGSmp);
	    if(nmod==3){
		gpu_calc_ngsmod(pclep[isim], pclmp[isim], pcleNGSmp[isim],recon->ngsmod->nmod,
				recon->ngsmod->aper_fcp, recon->ngsmod->ht,
				recon->ngsmod->scale, thetax, thetay,
				aper->ipcc, aper->imcc,
				cudata->plocs, nloc, iopdevl->p, cudata->pamp, stream);
	    }else{
		gpu_calc_ngsmod(NULL, NULL, pcleNGSmp[isim],recon->ngsmod->nmod,
				recon->ngsmod->aper_fcp, recon->ngsmod->ht,
				recon->ngsmod->scale, thetax, thetay,
				aper->ipcc, aper->imcc,
				cudata->plocs, nloc, iopdevl->p, cudata->pamp, stream);	
		TO_IMPLEMENT;/*mode decomposition. */
	    }
	}
    }else{
	if(nmod==3){
	    gpu_calc_ptt(pclep[isim], pclmp[isim], aper->ipcc, aper->imcc,
			 cudata->plocs, nloc, iopdevl->p, cudata->pamp, stream);
	}else{
	    TO_IMPLEMENT;
	}
    }
    if(do_psf_cov){
	if(parms->evl.psfngsr[ievl]!=0){/*do after ngs mode removal */
	    if(cudata->evlopd->p[ievl]) error("cudata->evlopd->p[%d] should be NULL\n", ievl);
	    cudata->evlopd->p[ievl]=iopdevl;/*record. */
	}
	if(parms->evl.psfngsr[ievl]!=2){
	    if(parms->evl.psfpttr[ievl]){
		curaddptt(iopdevl, cudata->plocs, -pclmp[isim][0], -pclmp[isim][1], -pclmp[isim][2], stream);
	    }else{
		curadds(iopdevl, -pclmp[isim][0], stream);
	    }
	    if(parms->evl.opdcov){
		curmm(&cudata->evlopdcov->p[ievl], 1, iopdevl, iopdevl, "nt", 1, handle);
	    }/*opdcov */
	    if(parms->evl.psfhist){
		/*Compute complex. */
		cuccell *psfs=psfcomp(iopdevl, nwvl, ievl, nloc, stream);
		cellarr_cuccell(simu->save->evlpsfhist[ievl], psfs, stream);
		if(parms->evl.psfmean){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			curaddcabs2(cudata->evlpsfcl->p+iwvl+nwvl*ievl, 1, psfs->p[iwvl], 1, stream);
		    }
		}
	    }else if(parms->evl.psfmean){
		psfcomp_r(cudata->evlpsfcl->p+nwvl*ievl, iopdevl, nwvl, ievl, nloc, 0, stream);
	    }
	}
    }
 end:
    CUDA_SYNC_STREAM;
    if(cudata->evlopd->p[ievl]!=iopdevl) curfree(iopdevl);
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
	gpu_set(evlgpu[ievl]);
	curmat *iopdevl=cudata->evlopd->p[ievl];
	if(!iopdevl) continue;
	cudaStream_t stream=evlstream[ievl];
	cublasHandle_t handle=evlhandle[ievl];
	gpu_ngsmod2science(iopdevl, parms, simu->recon, aper, cleNGSm, ievl, -1, stream);
	if(parms->evl.psfpttr[ievl]){
	    double ptt[3];
	    gpu_calc_ptt(NULL, ptt,  aper->ipcc, aper->imcc,
			 cudata->plocs, nloc, iopdevl->p, cudata->pamp, stream);
	    curaddptt(iopdevl, cudata->plocs, -ptt[0], -ptt[1], -ptt[2], stream);
	}
	if(parms->evl.opdcov){
	    curmm(&cudata->evlopdcov_ngsr->p[ievl], 1, iopdevl, iopdevl, "nt", 1, handle);
	}/*opdcov */
	if(parms->evl.psfhist||parms->evl.psfmean){
	    if(parms->evl.psfhist){
		/*Compute complex. */
		cuccell *psfs=psfcomp(iopdevl, nwvl, ievl, nloc, stream);
		cellarr_cuccell(simu->save->evlpsfhist_ngsr[ievl], psfs, stream);
		if(parms->evl.psfmean){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			curaddcabs2(cudata->evlpsfcl_ngsr->p+iwvl+nwvl*ievl, 1, psfs->p[iwvl], 1, stream);
		    }
		}
	    }else if(parms->evl.psfmean){
		psfcomp_r(cudata->evlpsfcl_ngsr->p+nwvl*ievl, iopdevl, nwvl, ievl, nloc, 0, stream);
	    }
	}
	CUDA_SYNC_STREAM;
	curfree(cudata->evlopd->p[ievl]);
	cudata->evlopd->p[ievl]=NULL;
    }
}
void gpu_perfevl_save(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(!parms->evl.nevl) return;
    const int isim=simu->isim;
    if(parms->evl.psfmean && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.psfmean)){
	info2("Output PSF\n");
	const int nwvl=parms->evl.nwvl;
	if(cudata->evlpsfol){
	    scell *temp=scellnew(nwvl, 1);
	    scell *temp2=scellnew(nwvl, 1);
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    if(parms->evl.psfol==2){
		scale=scale/parms->evl.npsf;
	    }
	    for(int im=0; im<NGPU; im++){
		gpu_set(im);
		gpu_curcell2s(&temp2, cudata->evlpsfol, 0);
		cudaStreamSynchronize(0);
		scelladd(&temp, 1, temp2, scale);
	    }
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		temp->p[iwvl]->header=evl_header(simu->parms, simu->aper, -1, iwvl);
		cellarr_smat(simu->save->evlpsfolmean, temp->p[iwvl]);
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
		curmat *temp=NULL;
		cudaStream_t stream=evlstream[ievl];
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    curadd(&temp, 0, cudata->evlpsfcl->p[iwvl+nwvl*ievl], scale, stream);
		    temp->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		    cellarr_cur(simu->save->evlpsfmean[ievl], temp, stream);
		    free(temp->header); temp->header=NULL;
		}
		curfree(temp);
	    }
	}
	if(cudata->evlpsfcl_ngsr){
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->evl.psf[ievl] || !parms->evl.psfngsr[ievl]) continue;
		gpu_set(evlgpu[ievl]);
		curmat *temp=NULL;
		cudaStream_t stream=evlstream[ievl];
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    curadd(&temp, 0, cudata->evlpsfcl_ngsr->p[iwvl+nwvl*ievl], scale, stream);
		    temp->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		    cellarr_cur(simu->save->evlpsfmean_ngsr[ievl], temp, stream);
		    free(temp->header); temp->header=NULL;
		}
		curfree(temp);
	    }
	}
    }
    if(parms->evl.opdcov && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.opdcov)){
	info("Output opdcov\n");
	double scale=1./(isim+1-parms->evl.psfisim);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!parms->evl.psf[ievl]|| parms->evl.psfngsr[ievl]==2) continue;
	    gpu_set(evlgpu[ievl]);
	    curmat *temp=NULL;
	    cudaStream_t stream=evlstream[ievl];
	    curadd(&temp, 0, cudata->evlopdcov->p[ievl], scale, stream);
	    cellarr_cur(simu->save->evlopdcov[ievl], temp, stream);
	    CUDA_SYNC_STREAM;
	    curfree(temp);
	}
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!parms->evl.psf[ievl]|| !parms->evl.psfngsr[ievl]) continue;
	    gpu_set(evlgpu[ievl]);
	    curmat *temp=NULL;
	    cudaStream_t stream=evlstream[ievl];
	    curadd(&temp, 0, cudata->evlopdcov_ngsr->p[ievl], scale, stream);
	    cellarr_cur(simu->save->evlopdcov_ngsr[ievl], temp, stream);
	    CUDA_SYNC_STREAM;
	    curfree(temp);
	}
    }
}
