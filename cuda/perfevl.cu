extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#include "cufft.h"
#include "cucmat.h"

static float (*cuplocs)[2]=NULL;
static float *cupamp=NULL;
static curcell *cusurfevl;
static int *cunembed=NULL;
static int *cupsfsize=NULL;
static int **cuembed=NULL;
static float *cuwvls=NULL;
static curcell *cuevlpsfol=NULL;
static curcell *cuevlpsfcl=NULL;
static cufftHandle *cuevlplan=NULL;
/**
  save aper_locs, aper_amp to GPU.
*/
static void gpu_plocs2gpu(loc_t *loc, dmat *amp){
    if(cuplocs) error("plocs is already Copied.\n");
    int nloc=loc->nloc;
    gpu_loc2dev(&cuplocs, loc);
    gpu_dbl2dev(&cupamp, amp->p, nloc);
}
/**
   Initialize perfevl
*/
void gpu_perfevl_init(const PARMS_T *parms, APER_T *aper){
    (void)parms;
    gpu_plocs2gpu(aper->locs, aper->amp);
    if(parms->evl.psfmean || parms->evl.psfhist){
	const int nwvl=parms->evl.nwvl;
	cunembed =(int*)  calloc(nwvl, sizeof(int));
	cupsfsize=(int*)  calloc(nwvl, sizeof(int));
	cuembed  =(int**) calloc(nwvl, sizeof(int*));
	cuwvls   =(float*)calloc(nwvl, sizeof(float));
	cuevlplan  = (cufftHandle*)calloc(nwvl, sizeof(cufftHandle));
	cuevlpsfcl = curcellnew(nwvl, parms->evl.nevl);
	cuevlpsfol = curcellnew(nwvl, 1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cunembed[iwvl]=(int)aper->nembed[iwvl];
	    cupsfsize[iwvl]=parms->evl.psfsize[iwvl];
	    cuwvls[iwvl]=parms->evl.wvl[iwvl];
	    gpu_long2dev(&cuembed[iwvl], aper->embed[iwvl], aper->locs->nloc);
	    DO(cufftPlan2d(&cuevlplan[iwvl], cunembed[iwvl], cunembed[iwvl], CUFFT_C2C));
	}//for iwvl
    }
}
/**
   Add surface to surfevl;
*/
void gpu_evlsurf2gpu(APER_T *aper){
    if(aper->opdadd){
	gpu_dcell2cu(&cusurfevl, aper->opdadd);
	dcellfree(aper->opdadd);
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
static cuccell *psfcomp(curmat *iopdevl, int nwvl, int nloc, cudaStream_t stream){
    cuccell *psfs=cuccellnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	cucmat *wvf=cucnew(cunembed[iwvl], cunembed[iwvl]);
	evl_embed_wvf_do<<<DIM(iopdevl->nx,256),0,stream>>>
	    (wvf->p, iopdevl->p, cupamp, cuembed[iwvl], nloc, cuwvls[iwvl]);
	CUFFT(cuevlplan[iwvl], wvf->p, CUFFT_FORWARD);
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
    return psfs;
}
/**
   Compute only PSF and add to result.
*/
static void psfcomp_r(curmat **psf, curmat *iopdevl, int iwvl, int nloc, int atomic, cudaStream_t stream){
    cucmat *wvf=cucnew(cunembed[iwvl], cunembed[iwvl]);
    evl_embed_wvf_do<<<DIM(iopdevl->nx,256),0,stream>>>
	(wvf->p, iopdevl->p, cupamp, cuembed[iwvl], nloc, cuwvls[iwvl]);
    CUFFT(cuevlplan[iwvl], wvf->p, CUFFT_FORWARD);
    if(!*psf) *psf=curnew(cupsfsize[iwvl], cupsfsize[iwvl]);
    if(atomic){
	evl_corner2center_abs2_atomic_do<<<DIM2((*psf)->nx,(*psf)->ny,16),0,stream>>>
	    ((*psf)->p, (*psf)->nx, (*psf)->ny, wvf->p, wvf->nx, wvf->ny);
    }else{
	evl_corner2center_abs2_do<<<DIM2((*psf)->nx,(*psf)->ny,16),0,stream>>>
	    ((*psf)->p, (*psf)->nx, (*psf)->ny, wvf->p, wvf->nx, wvf->ny);
    }
    cucfree(wvf);
}


/**
   Performance evaluation. Designed to replace perfevl_ievl in maos/perfevl.c
*/
void gpu_perfevl(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    const int ievl=info->start;
    assert(info->end==info->start+1);//only one evl.
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int isim=simu->isim;
    const int nmod=parms->evl.nmod;
    const int imoao=parms->evl.moao;
    const double dt=simu->dt;
    const int do_psf=(parms->evl.psfmean || parms->evl.psfhist) && isim>=parms->evl.psfisim;
    const int save_evlopd=parms->save.evlopd>0 && ((isim+1)%parms->save.evlopd)==0;
    const int nloc=aper->locs->nloc;
    const int nwvl=parms->evl.nwvl;
    const double thetax=parms->evl.thetax[ievl];
    const double thetay=parms->evl.thetay[ievl];
    //Setup pointers for easy usage
    PDMAT(simu->olmp->p[ievl],polmp);//OL mode for each dir
    PDMAT(simu->olep->p[ievl],polep);//OL error for each dir
    PDMAT(simu->clmp->p[ievl],pclmp);
    PDMAT(simu->clep->p[ievl],pclep);

    curmat *iopdevl=curnew(aper->locs->nloc, 1);
    cudaStream_t stream;
    STREAM_NEW(stream);
    /* iopdevl must be in device memory. 6 times slower if in host memory.*/
    if(cusurfevl && cusurfevl->p[ievl]){
	curcp(&iopdevl, cusurfevl->p[ievl], stream);
    }else{
	curset(iopdevl, 0, stream);
    }
    if(parms->sim.idealevl){
	gpu_dm2loc(iopdevl->p, cuplocs, nloc, cudmproj, parms->evl.hs[ievl], thetax, thetay,
		   parms->evl.misreg[0], parms->evl.misreg[1], 1, stream);
    }else if(simu->atm && !parms->sim.wfsalias){
	gpu_atm2loc(iopdevl->p, cuplocs, nloc, parms->evl.hs[ievl], thetax, thetay, 
		    parms->evl.misreg[0], parms->evl.misreg[1], isim*dt, 1, stream);
    }
    
    if(simu->telws){//Wind shake
	float tt=simu->telws->p[isim];
	float angle=simu->winddir?simu->winddir->p[0]:0;
	curaddptt(iopdevl, cuplocs, tt*cosf(angle), tt*sinf(angle), stream);
    }
    if(save_evlopd){
	cellarr_cur(simu->save->evlopdol[ievl], iopdevl, stream);
    }
    if(parms->plot.run){
	TO_IMPLEMENT;
    }
    if(nmod==3){
	gpu_calc_ptt(polep[isim], polmp[isim], aper->ipcc, aper->imcc,
		     cuplocs, nloc, iopdevl->p, cupamp, stream);
    }else{
	TO_IMPLEMENT;
    }

    if(parms->evl.psfmean &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
			     ||(parms->evl.psfol==2 && parms->evl.psf[ievl]))){
	//calculate Openloop PSF.
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    psfcomp_r(cuevlpsfol->p+iwvl, iopdevl, iwvl, nloc, parms->evl.psfol==2?1:0, stream);
	}
    }
    
    if(parms->sim.evlol) goto end;
    
    if(parms->evl.tomo){
	TO_IMPLEMENT;
    }else{
	gpu_dm2loc(iopdevl->p, cuplocs, nloc, cudmreal, parms->evl.hs[ievl], thetax, thetay,
	       parms->evl.misreg[0], parms->evl.misreg[1], -1, stream);
	if(imoao>-1){
	    TO_IMPLEMENT;
	}
    }
    if(save_evlopd){
	cellarr_cur(simu->save->evlopdcl[ievl], iopdevl, stream);
    }
    if(parms->plot.run){
	TO_IMPLEMENT;
    }
    if(parms->recon.split){
	if(parms->ndm<=2){
	    PDMAT(simu->cleNGSmp->p[ievl], pcleNGSmp);
	    if(nmod==3){
		gpu_calc_ngsmod(pclep[isim], pclmp[isim], pcleNGSmp[isim],recon->ngsmod->nmod,
				recon->ngsmod->aper_fcp, recon->ngsmod->ht,
				recon->ngsmod->scale, thetax, thetay,
				aper->ipcc, aper->imcc,
				cuplocs, nloc, iopdevl->p, cupamp, stream);
	    }else{
		gpu_calc_ngsmod(NULL, NULL, pcleNGSmp[isim],recon->ngsmod->nmod,
				recon->ngsmod->aper_fcp, recon->ngsmod->ht,
				recon->ngsmod->scale, thetax, thetay,
				aper->ipcc, aper->imcc,
				cuplocs, nloc, iopdevl->p, cupamp, stream);	
		TO_IMPLEMENT;//mode decomposition.
	    }
	}
    }else{
	if(nmod==3){
	    gpu_calc_ptt(pclep[isim], pclmp[isim], aper->ipcc, aper->imcc,
			 cuplocs, nloc, iopdevl->p, cupamp, stream);
	}else{
	    TO_IMPLEMENT;
	}
    }
    if(parms->evl.psf[ievl] && isim>=parms->evl.psfisim && do_psf){
	if(parms->evl.psfhist){
	    //Compute complex.
	    cuccell *psfs=psfcomp(iopdevl, nwvl, nloc, stream);
	    zcell *temp=NULL;
	    gpu_cuccell2z(&temp, psfs, stream);
	    CUDA_SYNC_STREAM;
	    cellarr_zcell(simu->save->evlpsfhist[ievl], temp);
	    zcellfree(temp);
	    if(parms->evl.psfmean){
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    curaddcabs2(cuevlpsfcl->p+iwvl+nwvl*ievl, 1, psfs->p[iwvl], 1, stream);
		}
	    }
	}else if(parms->evl.psfmean){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		psfcomp_r(cuevlpsfcl->p+iwvl+nwvl*ievl, iopdevl, iwvl, nloc, 0, stream);
	    }
	}
    }
 end:
    CUDA_SYNC_STREAM;
    STREAM_DONE(stream);
    curfree(iopdevl);
}

void gpu_perfevl_save(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    if(parms->evl.psfmean && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.psfmean)){
	info2("Output PSF\n");
	cudaStream_t stream;
	STREAM_NEW(stream);

	const int nwvl=parms->evl.nwvl;
	curcell *temp=curcellnew(nwvl, 1);
	if(cuevlpsfol){
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    if(parms->evl.psfol==2){
		scale=scale/parms->evl.npsf;
	    }
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		curadd(&temp->p[iwvl], 0, cuevlpsfol->p[iwvl], scale, stream);
		cellarr_cur(simu->save->evlpsfolmean, temp->p[iwvl], stream);
	    }
	}
	if(cuevlpsfcl){
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		    curadd(&temp->p[iwvl], 0, cuevlpsfcl->p[iwvl+nwvl*ievl], 
			   scale, stream);
		    cellarr_cur(simu->save->evlpsfmean[ievl], 
				temp->p[iwvl], stream);
		}
	    }
	}
	CUDA_SYNC_STREAM;
	curcellfree(temp);
	STREAM_DONE(stream);
    }
    if(parms->evl.opdcov && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.opdcov)){
	TO_IMPLEMENT;
    }
}
