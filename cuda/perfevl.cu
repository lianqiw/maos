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
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#include "cucmat.h"
#include "kernel.h"
#include "cudata.h"
#include "perf.h"

/** 
    save aper_locs, aper_amp to GPU.
*/

__global__ static void calc_ptt_do( Real *cc,
				    const Real (*restrict loc)[2], 
				    const int nloc,
				    const Real *restrict phi,
				    const Real *restrict amp){
    extern __shared__ Real ccb0[];
    Real *ccb[4];
    for(int i=0; i<4; i++){
	ccb[i]=ccb0+blockDim.x*i;
	ccb[i][threadIdx.x]=0.f;
    }
    int step=blockDim.x * gridDim.x; 
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	Real tmp=phi[i]*amp[i];
	ccb[0][threadIdx.x]+=tmp;
	ccb[1][threadIdx.x]+=tmp*loc[i][0];
	ccb[2][threadIdx.x]+=tmp*loc[i][1];
	ccb[3][threadIdx.x]+=tmp*phi[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    for(int i=0; i<4; i++){
		ccb[i][threadIdx.x]+=ccb[i][threadIdx.x+step];
	    }
	}
    }
    if(threadIdx.x<4){/*This is the first thread of a block. add block result to global. */
	atomicAdd(&cc[threadIdx.x], ccb[threadIdx.x][0]);
    }
}
__global__ static void calc_ngsmod_do( Real *cc,
				       const Real (*restrict loc)[2], 
				       const int nloc,
				       const Real *restrict phi,
				       const Real *restrict amp){
    extern __shared__ Real ccb0[];
    Real *ccb[7];
    for(int i=0; i<7; i++){
	ccb[i]=ccb0+blockDim.x*i;
	ccb[i][threadIdx.x]=0.f;
    }
    int step=blockDim.x * gridDim.x; 
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=step){
	const Real tmp=phi[i]*amp[i];
	const Real x=loc[i][0];
	const Real y=loc[i][1];
	ccb[0][threadIdx.x]+=tmp;
	ccb[1][threadIdx.x]+=tmp*x;
	ccb[2][threadIdx.x]+=tmp*y;
	ccb[3][threadIdx.x]+=tmp*x*x;
	ccb[4][threadIdx.x]+=tmp*y*y;
	ccb[5][threadIdx.x]+=tmp*x*y;
	ccb[6][threadIdx.x]+=tmp*phi[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    for(int i=0; i<7; i++){
		ccb[i][threadIdx.x]+=ccb[i][threadIdx.x+step];
	    }
	}
    }
    if(threadIdx.x<7){/*This is the first thread of a block. add block result to global. */
	atomicAdd(&cc[threadIdx.x], ccb[threadIdx.x][0]);
    }
}
/*
  Let M be the modal matrix of pistion/tip/tilt. Calculate M'*diag(amp)*phi
  where amp is the amptliude weighting.  */
static void calc_ptt(double *rmsout, double *coeffout, 
		     const double ipcc, const dmat *imcc,
		     const Real (*restrict loc)[2], 
		     const int nloc,
		     const Real *restrict phi,
		     const Real *restrict amp,
		     Real *cc, Real *ccb,
		     cudaStream_t stream
    ){
    /*sum with 16 blocks, each with 256 threads. */
    cudaMemsetAsync(cc, 0, sizeof(Real)*4, stream);
    calc_ptt_do<<<DIM(nloc, 128), 128*4*sizeof(Real), stream>>>(cc, loc, nloc, phi, amp);
    CUDA_SYNC_STREAM;
    cudaMemcpy(ccb, cc, 4*sizeof(Real), cudaMemcpyDeviceToHost);
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
	rmsout[2]=tot-ptt;/*PTTR*/
    }
}
static void calc_ngsmod(double *pttr_out, double *pttrcoeff_out,
			double *ngsmod_out, int nmod,
			double MCC_fcp, double ht, double scale,
			double thetax, double thetay,
			const double ipcc, const dmat *imcc,
			const Real (*restrict loc)[2], 
			const int nloc,
			const Real *restrict phi,
			const Real *restrict amp,
			const PARMS_T *parms,
			Real *cc, Real *ccb,
			cudaStream_t stream){
    double tot=0;
    CUDA_SYNC_STREAM;
    cudaMemsetAsync(cc, 0, 7*sizeof(Real), stream);
    CUDA_SYNC_STREAM;
    if(nmod==2){/*single DM. */
	calc_ptt_do<<<DIM(nloc,128),128*4*sizeof(Real),stream>>>(cc, loc, nloc, phi, amp);
    }else if(nmod>=5){/*AHST mode */
	calc_ngsmod_do<<<DIM(nloc,128),128*7*sizeof(Real),stream>>>(cc, loc, nloc, phi, amp);
    }else{
	info("nmod=%d\n", nmod);
	TO_IMPLEMENT;
    }
    CUDA_SYNC_STREAM;
    cudaMemcpy(ccb, cc, 7*sizeof(Real), cudaMemcpyDeviceToHost);
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



__global__ static void 
strehlcomp_do(Comp *strehlc, 
	      const Real *opd, const Real *amp, const int nloc, const Real kk){
    extern __shared__ Real sbx[];
    Real *sby=sbx+blockDim.x;
    sbx[threadIdx.x]=0;
    sby[threadIdx.x]=0;
    int skip=blockDim.x * gridDim.x ;
    Real s,c;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=skip){
	Z(sincos)(kk*opd[i], &s, &c);
	sbx[threadIdx.x]+=amp[i]*c;
	sby[threadIdx.x]+=amp[i]*s;
    }
    for(int step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    sbx[threadIdx.x]+=sbx[threadIdx.x+step];
	    sby[threadIdx.x]+=sby[threadIdx.x+step];
	}
    }
    if(threadIdx.x==0){
	if(strehlc){
	    atomicAdd((Real*)strehlc, sbx[0]);
	    atomicAdd((Real*)strehlc+1, sby[0]);
	}
	//donot try to accumuate x*x+y*y. that is not correct because of many blocks.
    }
}
/**
   Compute complex PSF and return.
*/
static void psfcomp(curmat *iopdevl, int nwvl, int ievl, int nloc, cudaStream_t stream){
    LOCK(cudata->perf->mutex);/*wvf is allocated per GPU.*/
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	cucmat *psf=cudata->perf->psfs->p[iwvl];
	if(cuperf_t::psfsize[iwvl]==1){
	    strehlcomp_do<<<REDUCE(nloc), DIM_REDUCE*sizeof(Comp),stream>>>
		(psf->p, iopdevl->p, cudata->perf->amp, nloc, 2.*M_PI/cuperf_t::wvls[iwvl]);
	}else{
	    cucmat *wvf=cudata->perf->wvf->p[iwvl];
	    cuczero(wvf, stream);
	    embed_wvf_do<<<DIM(iopdevl->nx,256),0,stream>>>
		(wvf->p, iopdevl->p, cudata->perf->amp, cudata->perf->embed[iwvl], nloc, cuperf_t::wvls[iwvl]);
	    CUFFT(cuperf_t::plan[iwvl+nwvl*ievl], wvf->p, CUFFT_FORWARD);
	    if(cuperf_t::psfsize[iwvl]<cuperf_t::nembed[iwvl]){
		corner2center_do<<<DIM2(psf->nx,psf->ny,16),0,stream>>>
		    (psf->p, psf->nx, psf->ny, wvf->p, wvf->nx, wvf->ny);
	    }else{
		fftshift_do<<<DIM2(psf->nx,psf->ny,16),0,stream>>>
		    (psf->p, psf->nx, psf->ny);
	    }
	}
    }
    UNLOCK(cudata->perf->mutex);
}
/**
   Compute only PSF and add to result.
*/
static void psfcomp_r(curmat **psf, curmat *iopdevl, int nwvl, int ievl, int nloc, int atomic, cudaStream_t stream){
    LOCK(cudata->perf->mutex);/*wvf is allocated per GPU.*/
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	cucmat *wvf=cudata->perf->wvf->p[iwvl];
	cuczero(wvf, stream);
	if(!psf[iwvl]) psf[iwvl]=curnew(cuperf_t::psfsize[iwvl], cuperf_t::psfsize[iwvl]);
	if(cuperf_t::psfsize[iwvl]==1){
	    strehlcomp_do<<<REDUCE(nloc), DIM_REDUCE*sizeof(Real)*2,stream>>>
		(wvf->p, iopdevl->p, cudata->perf->amp, nloc, 2.*M_PI/cuperf_t::wvls[iwvl]);
	    //do abs2.
	    addcabs2_do<<<1,1,0,stream>>>(psf[iwvl]->p, 1.f, wvf->p, 1.f, 1);
	}else{
	    embed_wvf_do<<<DIM(iopdevl->nx,256),0,stream>>>
		(wvf->p, iopdevl->p, cudata->perf->amp, cudata->perf->embed[iwvl], nloc, cuperf_t::wvls[iwvl]);
	    CUFFT(cuperf_t::plan[iwvl+nwvl*ievl], wvf->p, CUFFT_FORWARD);
	    if(atomic){
		corner2center_abs2_atomic_do<<<DIM2((psf[iwvl])->nx,(psf[iwvl])->ny,16),0,stream>>>
		    ((psf[iwvl])->p, (psf[iwvl])->nx, (psf[iwvl])->ny, wvf->p, wvf->nx, wvf->ny);
	    }else{
		corner2center_abs2_do<<<DIM2((psf[iwvl])->nx,(psf[iwvl])->ny,16),0,stream>>>
		    ((psf[iwvl])->p, (psf[iwvl])->nx, (psf[iwvl])->ny, wvf->p, wvf->nx, wvf->ny);
	    }
	}
    }
    UNLOCK(cudata->perf->mutex);
}
#define PERFEVL_WFE(pclep, pclmp, cleNGSmp)				\
    if(parms->recon.split){						\
	if(parms->ndm<=2){						\
	    PDMAT(cleNGSmp->p[ievl], pcleNGSmp);			\
	    calc_ngsmod(nmod==3?pclep[isim]:0, nmod==3?pclmp[isim]:0,	\
			pcleNGSmp[isim],recon->ngsmod->nmod,		\
			recon->ngsmod->aper_fcp, recon->ngsmod->ht,	\
			recon->ngsmod->scale, thetax, thetay,		\
			aper->ipcc, aper->imcc,				\
			cudata->perf->locs->p, nloc, iopdevl->p, cudata->perf->amp, parms, \
			cuperf_t::cc->p[ievl]->p,cuperf_t::ccb->p[ievl]->p, stream); \
	    if(nmod!=3){						\
		TO_IMPLEMENT;/*mode decomposition. */			\
	    }								\
	}								\
    }else{								\
	if(nmod==3){							\
	    calc_ptt(pclep[isim], pclmp[isim], aper->ipcc, aper->imcc,	\
		     cudata->perf->locs->p, nloc, iopdevl->p, cudata->perf->amp,	\
		     cuperf_t::cc->p[ievl]->p,cuperf_t::ccb->p[ievl]->p,stream); \
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
    gpu_set(cudata_t::evlgpu[ievl]);
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
    cudaStream_t stream=cuperf_t::stream[ievl];
    cublasHandle_t handle=cuperf_t::handle[ievl];
    curmat *iopdevl=cuperf_t::opd->p[ievl];
    curzero(iopdevl, stream);
    /* iopdevl must be in device memory. 6 times slower if in host memory.*/
    if(cuperf_t::surf && cuperf_t::surf->p[ievl]){
	curcp(&iopdevl, cuperf_t::surf->p[ievl], stream);
    }else{
	curset(iopdevl, 0, stream);
    }
    if(parms->sim.idealevl){
	gpu_dm2loc(iopdevl->p, cudata->perf->locs_dm[ievl], cudata->dmproj, cudata->ndm,
		   parms->evl.hs[ievl], thetax, thetay,
		   0,0,1, stream);
    }else if(simu->atm && !parms->sim.wfsalias){
	gpu_atm2loc(iopdevl->p, cudata->perf->locs, parms->evl.hs[ievl], thetax, thetay, 
		    0,0,isim*dt, 1, stream);
    }
    if(simu->telws){/*Wind shake */
	Real tt=simu->telws->p[isim];
	Real angle=simu->winddir?simu->winddir->p[0]:0;
	curaddptt(iopdevl, cudata->perf->locs->p, 0, tt*cosf(angle), tt*sinf(angle), stream);
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
	    curaddptt(opdcopy, cudata->perf->locs->p, -polmp[isim][0], -polmp[isim][1], -polmp[isim][2], stream);
	}else if(parms->evl.opdcov){/* remove piston*/
	    curcp(&opdcopy, iopdevl, stream);
	    curadd(opdcopy, -polmp[isim][0], stream);
	}else{
	    opdcopy=iopdevl;
	}
	if(parms->evl.opdcov){
	    if(parms->gpu.psf){
		curmm(&cudata->perf->opdcovol, 1, opdcopy, opdcopy, "nt", 1, handle);
		curadd(&cudata->perf->opdmeanol, 1, opdcopy, 1, stream);
	    }else{
		dmat *tmp=NULL;
		cp2cpu(&tmp, opdcopy, stream);
		dmm(&simu->evlopdcovol, 1,tmp, tmp, "nt", 1);
		dadd(&simu->evlopdmeanol, 1, tmp, 1);
		dfree(tmp);
	    }
	}
	if(parms->evl.psfmean){
	    psfcomp_r(cudata->perf->psfol->p, opdcopy, nwvl, ievl, nloc, parms->evl.psfol==2?1:0, stream);
	    if(opdcopy!=iopdevl){
		curfree(opdcopy);
	    }
	    if(!parms->gpu.psf){ /*need to move psf from GPU to CPU for accumulation.*/
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    add2cpu(&simu->evlpsfolmean->p[iwvl], 1, cudata->perf->psfol->p[iwvl], 1, stream);
		    curzero(cudata->perf->psfol->p[iwvl]); //do not accumulate in gpu.
		}
	    }
	}
    }
    if(parms->sim.evlol) goto end;
    if(parms->evl.tomo){
	TO_IMPLEMENT;
    }else{
	gpu_dm2loc(iopdevl->p, cudata->perf->locs_dm[ievl], cudata->dmreal, cudata->ndm, 
		   parms->evl.hs[ievl], thetax, thetay,
		   0,0,-1, stream);
	if(simu->ttmreal){
	    curaddptt(iopdevl, cudata->perf->locs->p, 0, -simu->ttmreal->p[0], -simu->ttmreal->p[1], stream);
	}
	if(imoao!=-1){
	    gpu_dm2loc(iopdevl->p, &cudata->perf->locs, cudata->dm_evl[ievl], 1,
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
		curaddptt(iopdevl, cudata->perf->locs->p, -pclmp[isim][0], -pclmp[isim][1], -pclmp[isim][2], stream);
	    }else{
		curadd(iopdevl, -pclmp[isim][0], stream);
	    }
	    if(parms->evl.opdcov){
		if(parms->gpu.psf){
		    curmm(&cuperf_t::opdcov->p[ievl], 1, iopdevl, iopdevl, "nt", 1, handle);
		    curadd(&cuperf_t::opdmean->p[ievl], 1, iopdevl, 1, stream);
		}else{
		    dmat *tmp=NULL;
		    cp2cpu(&tmp, iopdevl, stream);
		    dmm(&simu->evlopdcov->p[ievl], 1,tmp, tmp, "nt", 1);
		    dadd(&simu->evlopdmean->p[ievl], 1, tmp, 1);
		    dfree(tmp);
		}
	    }/*opdcov */
	    if(parms->evl.psfhist || parms->evl.psfmean){
		if(parms->evl.psfhist){
		    /*Compute complex. */
		    psfcomp(iopdevl, nwvl, ievl, nloc, stream);
		    cellarr_cuccell(simu->save->evlpsfhist[ievl], isim, cudata->perf->psfs, stream);
		    if(parms->evl.psfmean){
			for(int iwvl=0; iwvl<nwvl; iwvl++){
			    curaddcabs2(cuperf_t::psfcl->p+iwvl+nwvl*ievl, 1, 
					cudata->perf->psfs->p[iwvl], 1, stream);
			}
		    }
		}else if(parms->evl.psfmean){
		    psfcomp_r(cuperf_t::psfcl->p+nwvl*ievl, iopdevl, nwvl, ievl, nloc, 0, stream);
		}
		if(!parms->gpu.psf){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			add2cpu(&simu->evlpsfmean->p[iwvl+ievl*nwvl], 1, cuperf_t::psfcl->p[iwvl+ievl*nwvl], 1, stream);
			curzero(cuperf_t::psfcl->p[iwvl+ievl*nwvl]); 
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
	gpu_set(cudata_t::evlgpu[ievl]);
	curmat *iopdevl=cuperf_t::opd->p[ievl];
	cudaStream_t stream=cuperf_t::stream[ievl];
	cublasHandle_t handle=cuperf_t::handle[ievl];
	gpu_ngsmod2science(iopdevl, cudata->perf->locs->p, simu->recon->ngsmod, cleNGSm, 
			   parms->evl.thetax[ievl], parms->evl.thetay[ievl],
			   -1, stream);
	if(parms->evl.psfpttr[ievl]){
	    double ptt[3];
	    calc_ptt(NULL, ptt,  aper->ipcc, aper->imcc,
		     cudata->perf->locs->p, nloc, iopdevl->p, cudata->perf->amp, 
		     cuperf_t::cc->p[ievl]->p,cuperf_t::ccb->p[ievl]->p, stream);
	    curaddptt(iopdevl, cudata->perf->locs->p, -ptt[0], -ptt[1], -ptt[2], stream);
	}
	if(parms->evl.opdcov){
	    if(parms->gpu.psf){
		curmm(&cuperf_t::opdcov_ngsr->p[ievl], 1, iopdevl, iopdevl, "nt", 1, handle);
		curadd(&cuperf_t::opdmean_ngsr->p[ievl], 1, iopdevl, 1, stream);
	    }else{
		dmat *tmp=NULL;
		cp2cpu(&tmp, iopdevl, stream);
		dmm(&simu->evlopdcov_ngsr->p[ievl], 1,tmp, tmp, "nt", 1);
		dadd(&simu->evlopdmean_ngsr->p[ievl], 1, tmp, 1);
		dfree(tmp);
	    }
	}/*opdcov */
	if(parms->evl.psfhist||parms->evl.psfmean){
	    if(parms->evl.psfhist){
		/*Compute complex. */
		psfcomp(iopdevl, nwvl, ievl, nloc, stream);
		cellarr_cuccell(simu->save->evlpsfhist_ngsr[ievl], simu->isim, cudata->perf->psfs, stream);
		if(parms->evl.psfmean){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			curaddcabs2(cuperf_t::psfcl_ngsr->p+iwvl+nwvl*ievl, 1, 
				    cudata->perf->psfs->p[iwvl], 1, stream);
		    }
		}
	    }else if(parms->evl.psfmean){
		psfcomp_r(cuperf_t::psfcl_ngsr->p+nwvl*ievl, iopdevl, nwvl, ievl, nloc, 0, stream);
	    }
	    if(!parms->gpu.psf){
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    add2cpu(&simu->evlpsfmean_ngsr->p[iwvl+ievl*nwvl], 1, cuperf_t::psfcl_ngsr->p[iwvl+ievl*nwvl], 1, stream);
		    curzero(cuperf_t::psfcl_ngsr->p[iwvl+ievl*nwvl]); 
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
	if(cudata->perf->psfol){
	    /*copy the PSF accumulated in all the GPUs to CPU.*/
	    X(cell) *temp=X(cellnew)(nwvl, 1);
	    X(cell) *temp2=X(cellnew)(nwvl, 1);
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    if(parms->evl.psfol==2){
		scale=scale/parms->evl.npsf;
	    }
	    for(int im=0; im<NGPU; im++){
		gpu_set(im);
		cp2cpu(&temp2, cudata->perf->psfol, 0);
		cudaStreamSynchronize(0);
		X(celladd)(&temp, 1, temp2, scale);
	    }
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(!temp || !temp->p[iwvl]) continue;
		temp->p[iwvl]->header=evl_header(simu->parms, simu->aper, -1, iwvl);
		cellarr_mat(simu->save->evlpsfolmean, isim*nwvl+iwvl, temp->p[iwvl]);
		free(temp->p[iwvl]->header); temp->p[iwvl]->header=NULL;
	    }
	    X(cellfree)(temp);
	    X(cellfree)(temp2);
	}
	if(cuperf_t::psfcl){
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->evl.psf[ievl] || parms->evl.psfngsr[ievl]==2) continue;
		gpu_set(cudata_t::evlgpu[ievl]);
		cudaStream_t stream=cuperf_t::stream[ievl];
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    curmat *pp=cuperf_t::psfcl->p[iwvl+nwvl*ievl];
		    curscale(pp, scale, stream);
		    if(!pp->header){
			pp->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		    }
		    cellarr_cur(simu->save->evlpsfmean[ievl], isim*nwvl+iwvl, pp, stream);
		    curscale(pp, 1.f/scale, stream);
		}
	    }
	}
	if(cuperf_t::psfcl_ngsr){
	    double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->evl.psf[ievl] || !parms->evl.psfngsr[ievl]) continue;
		gpu_set(cudata_t::evlgpu[ievl]);
		cudaStream_t stream=cuperf_t::stream[ievl];
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    curmat *pp=cuperf_t::psfcl_ngsr->p[iwvl+nwvl*ievl];
		    curscale(pp, scale, stream);
		    if(!pp->header){
			pp->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		    }
		    cellarr_cur(simu->save->evlpsfmean_ngsr[ievl], isim*nwvl+iwvl, pp, stream);
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
	    gpu_set(cudata_t::evlgpu[ievl]);
	    cudaStream_t stream=cuperf_t::stream[ievl];
	    curmat *pp;
	    {
		pp=cuperf_t::opdcov->p[ievl];
		curscale(pp, scale, stream);
		cellarr_cur(simu->save->evlopdcov[ievl], isim, pp, stream);
		curscale(pp, 1./scale, stream);
	    }
	    {
		pp=cuperf_t::opdmean->p[ievl];
		curscale(pp, scale, stream);
		cellarr_cur(simu->save->evlopdmean[ievl], isim, pp, stream);
		curscale(pp, 1./scale, stream);
	    }
	}
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!parms->evl.psf[ievl]|| !parms->evl.psfngsr[ievl]) continue;
	    gpu_set(cudata_t::evlgpu[ievl]);
	    cudaStream_t stream=cuperf_t::stream[ievl];
	    curmat *pp;
	    {
		pp=cuperf_t::opdcov_ngsr->p[ievl];
		curscale(pp, scale, stream);
		cellarr_cur(simu->save->evlopdcov_ngsr[ievl], isim, pp, stream);
		curscale(pp, 1./scale, stream);
	    }
	    {
		pp=cuperf_t::opdmean_ngsr->p[ievl];
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
		X(mat) *temp=NULL;
		X(mat) *temp2=NULL;
		for(int im=0; im<NGPU; im++){
		    gpu_set(im);
		    cp2cpu(&temp2, cudata->perf->opdcovol, 0);
		    cudaStreamSynchronize(0);
		    X(add)(&temp, 1, temp2, scale);
		}
		cellarr_mat(simu->save->evlopdcovol, isim, temp);
		X(free)(temp);
		X(free)(temp2);
	    }
	    {
		X(mat) *temp=NULL;
		X(mat) *temp2=NULL;
		for(int im=0; im<NGPU; im++){
		    gpu_set(im);
		    cp2cpu(&temp2, cudata->perf->opdmeanol, 0);
		    cudaStreamSynchronize(0);
		    X(add)(&temp, 1, temp2, scale);
		}
		cellarr_mat(simu->save->evlopdmeanol, isim, temp);
		X(free)(temp);
		X(free)(temp2);
	    }
	}
    }
}
