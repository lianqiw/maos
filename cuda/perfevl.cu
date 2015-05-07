/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "gpu.h"
#include "utils.h"
#include "accphi.h"
#include "cucmat.h"
#include "kernel.h"
#include "cudata.h"
#include "perf.h"
#undef TIMING
#define TIMING 0
#if !TIMING
#undef TIC
#undef tic
#undef toc
#define TIC
#define tic
#define ctoc(A)
#else
#define ctoc(A) toc2(A)
#endif
/** 
    save aper_locs, aper_amp to GPU.
*/
const int TT_NBX=128;//Number of thread in a block. (for reduction).
__global__ static void calc_ptt_do( Real *cc,
				    const Real (*restrict loc)[2], 
				    const int nloc,
				    const Real *restrict phi,
				    const Real *restrict amp){
    __shared__ Real ccb[4][TT_NBX];
    for(int i=0; i<4; i++){
	ccb[i][threadIdx.x]=0.f;
    }
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=blockDim.x * gridDim.x){
	const Real tmp=phi[i]*amp[i];
	ccb[0][threadIdx.x]+=tmp*phi[i];
	ccb[1][threadIdx.x]+=tmp;
	ccb[2][threadIdx.x]+=tmp*loc[i][0];
	ccb[3][threadIdx.x]+=tmp*loc[i][1];
    }
    for(int step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    for(int i=0; i<4; i++){
		ccb[i][threadIdx.x]+=ccb[i][threadIdx.x+step];
	    }
	}
    }
    if(threadIdx.x<4){
	atomicAdd(&cc[threadIdx.x], ccb[threadIdx.x][0]);
    }
}

__global__ static void calc_ngsmod_do( Real *cc,
				       const Real (*restrict loc)[2], 
				       const int nloc,
				       const Real *restrict phi,
				       const Real *restrict amp){
    __shared__ Real ccb[7][TT_NBX];
    for(int i=0; i<7; i++){
	ccb[i][threadIdx.x]=0.f;
    }
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=blockDim.x * gridDim.x){
	const Real tmp=phi[i]*amp[i];
	const Real x=loc[i][0];
	const Real y=loc[i][1];
	ccb[0][threadIdx.x]+=tmp*phi[i];
	ccb[1][threadIdx.x]+=tmp;
	ccb[2][threadIdx.x]+=tmp*x;
	ccb[3][threadIdx.x]+=tmp*y;
	ccb[4][threadIdx.x]+=tmp*x*x;
	ccb[5][threadIdx.x]+=tmp*y*y;
	ccb[6][threadIdx.x]+=tmp*x*y;
    }
    for(int step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
#pragma unroll
	    for(int i=0; i<7; i++){
		ccb[i][threadIdx.x]+=ccb[i][threadIdx.x+step];
	    }
	}
    }
    __syncthreads();
    if(threadIdx.x<7){
	atomicAdd(&cc[threadIdx.x], ccb[threadIdx.x][0]);
    }
}
/*
  Let M be the modal matrix of pistion/tip/tilt. Calculate M'*diag(amp)*phi
  where amp is the amptliude weighting.  */
static int calc_ptt_post(double *rmsout, double *coeffout, 
		     const double ipcc, const dmat *imcc,
		     Real *ccb){
    double coeff[3];
    double tot=ccb[0];
    coeff[0]=ccb[1]; coeff[1]=ccb[2]; coeff[2]=ccb[3]; 
    if(coeffout){
	dmulvec3(coeffout, imcc, coeff);
    }
    int ans=0;
    if(rmsout){
	double pis=ipcc*coeff[0]*coeff[0];/*piston mode variance */
	double ptt=dwdot3(coeff, imcc, coeff);/*p/t/t mode variance. */
	rmsout[0]=tot-pis;/*PR */
	rmsout[1]=ptt-pis;/*TT */
	rmsout[2]=tot-ptt;/*PTTR*/
	if(tot*1.01<pis || tot*1.01<ptt){//sanity check. allow round off error
	    warning("tot=%g, pis=%g, ptt=%g\n", tot, pis, ptt);
	    ans=1;
	}
    }
    return ans;
}
static int calc_ngsmod(double *pttr_out, double *pttrcoeff_out,
		       double *ngsmod_out, int nmod,
		       double MCC_fcp, double ht, double scale,
		       double thetax, double thetay,
		       const double ipcc, const dmat *imcc,
		       const PARMS_T *parms,
		       Real *ccb){
    double tot=(double)ccb[0];
    double coeff[6];//convert to double
    coeff[0]=ccb[1]; coeff[1]=ccb[2]; 
    coeff[2]=ccb[3]; coeff[3]=ccb[4];
    coeff[4]=ccb[5]; coeff[5]=ccb[6];
    
    if(pttrcoeff_out){//p/t/t
	memset(pttrcoeff_out, 0, sizeof(double)*3);
	dmulvec(pttrcoeff_out, imcc, coeff, 1);
    }
    int ans=0;
    if(pttr_out){
	//compute TT removed wavefront variance as a side product 
	double pis=ipcc*coeff[0]*coeff[0];
	double ptt=dwdot3(coeff, imcc, coeff);
	pttr_out[0]=tot-pis;//PR
	pttr_out[1]=ptt-pis;//TT
	pttr_out[2]=tot-ptt;//PTTR
	if(tot<pis || tot<ptt || ptt<pis || pis<0){
	    warning("tot=%g, pis=%g, ptt=%g\n", tot, pis, ptt);
	    ans=1;
	}
    }
    //don't use +=. need locking
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
    return ans;
}



__global__ static void 
strehlcomp_do(Comp *strehlc, 
	      const Real *opd, const Real *amp, const int nloc, const Real kk){
    extern __shared__ Real sbx[];
    Real *sby=sbx+blockDim.x;
    sbx[threadIdx.x]=0;
    sby[threadIdx.x]=0;
    Real s,c;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nloc; i+=blockDim.x * gridDim.x){
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
#define PERFEVL_WFE_GPU(cc,ccb)						\
    if((parms->recon.split && recon->ngsmod->nmod==2)			\
       || (!parms->recon.split && parms->evl.nmod==3)){			\
	cudaMemsetAsync(cc, 0, 4*sizeof(Real), stream);			\
	calc_ptt_do<<<DIM(nloc, TT_NBX), 0, stream>>>			\
	    (cc, cudata->perf->locs->p, nloc, iopdevl->p, cudata->perf->amp); \
	cudaMemcpyAsync(ccb, cc, 4*sizeof(Real), cudaMemcpyDeviceToHost, stream); \
    }else if(parms->recon.split){					\
	cudaMemsetAsync(cc, 0, 7*sizeof(Real), stream);			\
	calc_ngsmod_do<<<DIM(nloc,TT_NBX),0,stream>>>			\
	    (cc, cudata->perf->locs->p, nloc, iopdevl->p, cudata->perf->amp);\
	cudaMemcpyAsync(ccb, cc, 7*sizeof(Real), cudaMemcpyDeviceToHost, stream); \
    }

#define PERFEVL_WFE_CPU(ans, pclep, pclmp, cleNGSmp, ccb)		\
    if(nmod!=3){							\
	TO_IMPLEMENT;/*mode decomposition. */				\
    }									\
    int ans=0;								\
    if(parms->recon.split){						\
	double *pcleNGSmp=COLUMN(cleNGSmp->p[ievl], isim);		\
	ans=calc_ngsmod(nmod==3?pclep:0, nmod==3?pclmp:0,		\
			pcleNGSmp,recon->ngsmod->nmod,			\
			recon->ngsmod->aper_fcp, recon->ngsmod->ht,	\
			recon->ngsmod->scale, thetax, thetay,		\
			aper->ipcc, aper->imcc,	parms, ccb);		\
    }else{								\
	ans=calc_ptt_post(pclep, pclmp, aper->ipcc, aper->imcc, ccb);	\
    }

/**
   Performance evaluation. Designed to replace perfevl_ievl in maos/perfevl.c
*/
void gpu_perfevl_queue(thread_t *info){
    TIC;tic;
    SIM_T *simu=(SIM_T*)info->data;
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int isim=simu->isim;
    const int imoao=parms->evl.moao;
    const int nloc=aper->locs->nloc;
    const int nwvl=parms->evl.nwvl;
    for(int ievl=info->start; ievl<info->end; ievl++){
	gpu_set(cudata_t::evlgpu[ievl]);
	//info2("thread %ld gpu %d ievl %d start\n", thread_id(), cudata->igpu, ievl);
	const int do_psf_cov=(parms->evl.psfmean || parms->evl.psfhist || parms->evl.cov) 
	    && isim>=parms->evl.psfisim && parms->evl.psf->p[ievl]!=0;
	const int save_evlopd=parms->save.evlopd>0 && ((isim+1)%parms->save.evlopd)==0;
	const double thetax=parms->evl.thetax->p[ievl];
	const double thetay=parms->evl.thetay->p[ievl];
 
	cudaStream_t stream=cuperf_t::stream[ievl];
	cublasHandle_t handle=cuperf_t::handle[ievl];
	curmat *iopdevl=cuperf_t::opd->p[ievl];
	// iopdevl must be in device memory. 6 times slower if in host memory.
	if(cuperf_t::surf && cuperf_t::surf->p[ievl]){
	    curcp(&iopdevl, cuperf_t::surf->p[ievl], stream);
	}else{
	    curset(iopdevl, 0, stream);
	}
	if(parms->sim.idealevl){
	    gpu_dm2loc(iopdevl->p, cudata->perf->locs_dm[ievl], cudata->dmproj, cudata->ndm,
		       parms->evl.hs->p[ievl], thetax, thetay, 0,0,1, stream);
	}else if(simu->atm && !parms->sim.wfsalias){
	    gpu_atm2loc(iopdevl->p, cudata->perf->locs, parms->evl.hs->p[ievl], thetax, thetay, 
			0,0,parms->sim.dt,isim, 1, stream);
	}
	if(simu->telws){//Wind shake 
	    Real tt=simu->telws->p[isim];
	    Real angle=simu->winddir?simu->winddir->p[0]:0;
	    curaddptt(iopdevl, cudata->perf->locs->p, 0, tt*cosf(angle), tt*sinf(angle), stream);
	}
	if(save_evlopd){
	    cellarr_cur(simu->save->evlopdol[ievl], isim, iopdevl, stream);
	}
	if(parms->plot.run){
	    drawopdamp_gpu("OL", aper->locs, iopdevl, stream, aper->amp1->p, NULL,
			   "Science Open Loop OPD", "x (m)", "y (m)", "OL %d", ievl);
	}
	PERFEVL_WFE_GPU(cuperf_t::cc_ol->p[ievl]->p, cuperf_t::ccb_ol[ievl]);
	if((parms->evl.psfmean  || parms->evl.cov)
	   && isim>=parms->evl.psfisim 
	   &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
	      ||(parms->evl.psfol==2 && parms->evl.psf->p[ievl]))){
	    //calculate Openloop PSF. we also test psfisim to synchronize with psfcl.
	    curmat *opdcopy=NULL;
	    curmv(cuperf_t::coeff->p[ievl]->p, 0, cudata->perf->imcc, 
		  cuperf_t::cc_ol->p[ievl]->p, 'n', 1, handle);
	    curcp(&opdcopy, iopdevl, stream);
	    if(parms->evl.pttr->p[ievl]){//remove piston/tip/tilt
		curaddptt(opdcopy, cudata->perf->locs->p, cuperf_t::coeff->p[ievl]->p, -1,-1,-1,stream);
	    }else{//remove piston only
		curaddptt(opdcopy, cudata->perf->locs->p, cuperf_t::coeff->p[ievl]->p, -1, 0, 0, stream);
	    }
	    if(parms->evl.cov){
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
		if(!parms->gpu.psf){ //need to move psf from GPU to CPU for accumulation.
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			add2cpu(&simu->evlpsfolmean->p[iwvl], 1, cudata->perf->psfol->p[iwvl], 1, stream);
			curzero(cudata->perf->psfol->p[iwvl]); //do not accumulate in gpu.
		    }
		}
	    }
	}
	if(parms->sim.evlol) continue;
	if(parms->evl.tomo){
	    TO_IMPLEMENT;
	}else{
	    gpu_dm2loc(iopdevl->p, cudata->perf->locs_dm[ievl], cudata->dmreal, cudata->ndm, 
		       parms->evl.hs->p[ievl], thetax, thetay,
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

	if(parms->plot.run){
	    drawopdamp_gpu("CL", aper->locs,iopdevl, stream , aper->amp1->p, NULL,
			   "Science Closed loop OPD", "x (m)", "y (m)", "CL %d", ievl);
	}
	PERFEVL_WFE_GPU(cuperf_t::cc_cl->p[ievl]->p, cuperf_t::ccb_cl[ievl]);
	if(do_psf_cov && parms->evl.psfngsr->p[ievl]!=2){//also do normal one.
	    curmv(cuperf_t::coeff->p[ievl]->p, 0, cudata->perf->imcc, 
		  cuperf_t::cc_cl->p[ievl]->p, 'n', 1, handle);
	    if(parms->evl.pttr->p[ievl]){
		curaddptt(iopdevl, cudata->perf->locs->p, cuperf_t::coeff->p[ievl]->p, -1, -1, -1, stream);
	    }else{
		curaddptt(iopdevl, cudata->perf->locs->p, cuperf_t::coeff->p[ievl]->p, -1, 0, 0, stream);
	    }
	    if(parms->evl.cov){
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
	    }//opdcov 
	    if(parms->evl.psfhist || parms->evl.psfmean){
		if(parms->evl.psfhist){
		    //Compute complex. 
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
	//info2("thread %ld gpu %d ievl %d queued\n", thread_id(), cudata->igpu, ievl);
	ctoc("queued");
    }//for ievl
}
void gpu_perfevl_sync(thread_t *info){
    TIC;tic;
    SIM_T *simu=(SIM_T*)info->data;
    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int nmod=parms->evl.nmod;
    for(int ievl=info->start; ievl<info->end; ievl++){
	/*lock the mutex because iopdevl, evlwvf is allocated per GPU.*/
	gpu_set(cudata_t::evlgpu[ievl]);
	cudaStream_t stream=cuperf_t::stream[ievl];
	const double thetax=parms->evl.thetax->p[ievl];
	const double thetay=parms->evl.thetay->p[ievl];
	/*Setup pointers for easy usage */
	double *polmp=COLUMN(simu->olmp->p[ievl], isim);
	double *pclmp=COLUMN(simu->clmp->p[ievl], isim);
	double *polep=COLUMN(simu->olep->p[ievl], isim);
	double *pclep=COLUMN(simu->clep->p[ievl], isim);
	CUDA_SYNC_STREAM;
	PERFEVL_WFE_CPU(ans1, polep, polmp, simu->oleNGSmp, cuperf_t::ccb_ol[ievl]);
	PERFEVL_WFE_CPU(ans2, pclep, pclmp, simu->cleNGSmp, cuperf_t::ccb_cl[ievl]);
	if(ans1 || ans2){
	    warning("Perfevl fails, redo\n");
	    gpu_perfevl_queue(info);
	    gpu_perfevl_sync(info);
	}
	//info2("thread %ld gpu %d ievl %d end\n", thread_id(), cudata->igpu, ievl);
    }//for ievl
    ctoc("done");
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
	if(parms->evl.psfngsr->p[ievl]==0){
	    continue;
	}
	warning("Compare with CPU code to verify accuracy. Need to verify focus mode\n");
	gpu_set(cudata_t::evlgpu[ievl]);
	curmat *iopdevl=cuperf_t::opd->p[ievl];
	cudaStream_t stream=cuperf_t::stream[ievl];
	cublasHandle_t handle=cuperf_t::handle[ievl];
	gpu_ngsmod2science(iopdevl, cudata->perf->locs->p, simu->recon->ngsmod, cleNGSm, 
			   parms->evl.thetax->p[ievl], parms->evl.thetay->p[ievl],
			   -1, stream);
	if(parms->evl.pttr->p[ievl]){
	    double ptt[3];
	    calc_ptt_do<<<DIM(nloc,TT_NBX), 0, stream>>>	
		(cuperf_t::cc_cl->p[ievl]->p, cudata->perf->locs->p, nloc, iopdevl->p, cudata->perf->amp); 
	    cudaMemcpyAsync(cuperf_t::ccb_cl[ievl], cuperf_t::cc_cl->p[ievl]->p, 
			    4*sizeof(Real), cudaMemcpyDeviceToHost, stream); 
	    calc_ptt_post(NULL, ptt,  aper->ipcc, aper->imcc, cuperf_t::ccb_cl[ievl]);
	    curaddptt(iopdevl, cudata->perf->locs->p, -ptt[0], -ptt[1], -ptt[2], stream);
	}
	if(parms->evl.cov){
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
		if(!parms->evl.psf->p[ievl] || parms->evl.psfngsr->p[ievl]==2) continue;
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
		if(!parms->evl.psf->p[ievl] || !parms->evl.psfngsr->p[ievl]) continue;
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
    if(parms->evl.cov && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.cov)){
	info2("Step %d: Output opdcov\n", isim);
	double scale=1./(isim+1-parms->evl.psfisim);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!parms->evl.psf->p[ievl]|| parms->evl.psfngsr->p[ievl]==2) continue;
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
	    if(!parms->evl.psf->p[ievl]|| !parms->evl.psfngsr->p[ievl]) continue;
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
