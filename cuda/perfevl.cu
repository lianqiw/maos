/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "utils.h"
#include "accphi.h"
#include "cucmat.h"
#include "kernel.h"
#include "cudata.h"
#include "perf.h"
#if !USE_CPP
extern "C"{
#endif
#include "../maos/utils.h"
#include "../maos/ahst.h"
#if !USE_CPP
}
#endif
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
//wraps calc_ptt_do
static void calc_ptt(Real *cc,
		     const Real (*restrict loc)[2], 
		     const int nloc,
		     const Real *restrict phi,
		     const Real *restrict amp, stream_t &stream){
    cudaMemsetAsync(cc, 0, 4*sizeof(Real), stream);			
    calc_ptt_do<<<DIM(nloc, TT_NBX), 0, stream>>>			
	(cc, loc, nloc, phi, amp);
}
/*
  Let M be the modal matrix of pistion/tip/tilt. Calculate M'*diag(amp)*phi
  where amp is the amptliude weighting.  */
static int calc_ptt_post(double *rmsout, double *coeffout, 
			 const double ipcc, const dmat *imcc,
			 const Real *ccb){
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
	if(tot+1e-18<pis || tot+1e-18<ptt){//sanity check. allow round off error
	    warning("tot=%g, pis=%g, ptt=%g\n", tot, pis, ptt);
	    ans=1;
	}
    }
    return ans;
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
static void calc_ngsmod(Real *cc,
			const Real (*restrict loc)[2], 
			const int nloc,
			const Real *restrict phi,
			const Real *restrict amp,
			stream_t&stream){
    cudaMemsetAsync(cc, 0, 7*sizeof(Real), stream);			
    calc_ngsmod_do<<<DIM(nloc,TT_NBX),0,stream>>>
	(cc, loc, nloc, phi, amp);
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
static void psfcomp(cuccell psfs, const curmat &iopdevl, int nwvl, int ievl, int nloc, cudaStream_t stream){
    cucmat wvf;
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	cucmat &psf=psfs[iwvl];
	if(cuglobal->perf.psfsize[iwvl]==1){
	    strehlcomp_do<<<REDUCE(nloc), DIM_REDUCE*sizeof(Comp),stream>>>
		(psf(), iopdevl(), cudata->perf.amp, nloc, 2.*M_PI/cuglobal->perf.wvls[iwvl]);
	}else{
	    if(wvf.Nx()!=cuglobal->perf.nembed[iwvl]){
		wvf=cucmat(cuglobal->perf.nembed[iwvl], cuglobal->perf.nembed[iwvl]);
	    }else{
		cuzero(wvf, stream);
	    }
	    embed_wvf_do<<<DIM(iopdevl.Nx(),256),0,stream>>>
		(wvf(), iopdevl(), cudata->perf.amp, cudata->perf.embed[iwvl], nloc, cuglobal->perf.wvls[iwvl]);
	    CUFFT(cuglobal->perf.plan[iwvl+nwvl*ievl], wvf(), CUFFT_FORWARD);
	    if(cuglobal->perf.psfsize[iwvl]<cuglobal->perf.nembed[iwvl]){
		corner2center_do<<<DIM2(psf.Nx(),psf.Ny(),16),0,stream>>>
		    (psf(), psf.Nx(), psf.Ny(), wvf(), wvf.Nx(), wvf.Ny());
	    }else{
		fftshift_do<<<DIM2(psf.Nx(),psf.Ny(),16),0,stream>>>
		    (psf(), psf.Nx(), psf.Ny());
	    }
	}
    }
}
/**
   Compute only PSF and add to result.
*/
static void psfcomp_r(curmat *psf, const curmat &iopdevl, int nwvl, int ievl, int nloc, int atomic, cudaStream_t stream){
    cucmat wvf;
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	//cucmat &wvf=cudata->perf.wvf[iwvl];
	if(wvf.Nx()!=cuglobal->perf.nembed[iwvl]){
	    wvf=cucmat(cuglobal->perf.nembed[iwvl], cuglobal->perf.nembed[iwvl]);
	}else{
	    cuzero(wvf, stream);
	}
	if(!psf[iwvl]) psf[iwvl]=curmat(cuglobal->perf.psfsize[iwvl], cuglobal->perf.psfsize[iwvl]);
	if(cuglobal->perf.psfsize[iwvl]==1){
	    strehlcomp_do<<<REDUCE(nloc), DIM_REDUCE*sizeof(Real)*2,stream>>>
		(wvf(), iopdevl(), cudata->perf.amp, nloc, 2.*M_PI/cuglobal->perf.wvls[iwvl]);
	    //do abs2.
	    addcabs2_do<<<1,1,0,stream>>>(psf[iwvl](), 1.f, wvf(), 1.f, 1);
	}else{
	    embed_wvf_do<<<DIM(iopdevl.Nx(),256),0,stream>>>
		(wvf(), iopdevl(), cudata->perf.amp, cudata->perf.embed[iwvl], nloc, cuglobal->perf.wvls[iwvl]);
	    CUFFT(cuglobal->perf.plan[iwvl+nwvl*ievl], wvf(), CUFFT_FORWARD);
	    if(atomic){
		corner2center_abs2_atomic_do<<<DIM2((psf[iwvl]).Nx(),(psf[iwvl]).Ny(),16),0,stream>>>
		    ((psf[iwvl])(), (psf[iwvl]).Nx(), (psf[iwvl]).Ny(), wvf(), wvf.Nx(), wvf.Ny());
	    }else{
		corner2center_abs2_do<<<DIM2((psf[iwvl]).Nx(),(psf[iwvl]).Ny(),16),0,stream>>>
		    ((psf[iwvl])(), (psf[iwvl]).Nx(), (psf[iwvl]).Ny(), wvf(), wvf.Nx(), wvf.Ny());
	    }
	}
    }
}
#define PERFEVL_WFE_GPU(cc,ccb)						\
    if((parms->recon.split && recon->ngsmod->nmod==2)			\
       || (!parms->recon.split && parms->evl.nmod==3)){			\
	calc_ptt(cc, cudata->perf.locs(), nloc, iopdevl(), cudata->perf.amp, stream); \
	cudaMemcpyAsync(ccb, cc, 4*sizeof(Real), cudaMemcpyDeviceToHost, stream); \
    }else if(parms->recon.split){					\
	calc_ngsmod(cc, cudata->perf.locs(), nloc, iopdevl(), cudata->perf.amp, stream); \
	cudaMemcpyAsync(ccb, cc, 7*sizeof(Real), cudaMemcpyDeviceToHost, stream); \
    }

#define PERFEVL_WFE_CPU(ans, pclep, pclmp, cleNGSmp, ccb)		\
    if(nmod!=3){							\
	TO_IMPLEMENT;/*mode decomposition. */				\
    }									\
    int ans=0;								\
    if(parms->recon.split){						\
	double *pcleNGSmp=PCOL(cleNGSmp->p[ievl], isim);		\
	double coeff[6];						\
	coeff[0]=ccb[1]; coeff[1]=ccb[2];				\
	coeff[2]=ccb[3]; coeff[3]=ccb[4];				\
	coeff[4]=ccb[5]; coeff[5]=ccb[6];				\
	calc_ngsmod_post(nmod==3?pclep:0, nmod==3?pclmp:0,		\
			 pcleNGSmp,ccb[0],coeff,recon->ngsmod, aper,thetax,thetay); \
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
    const int isim=simu->perfisim;
    const int imoao=parms->evl.moao;
    const int nloc=aper->locs->nloc;
    const int nwvl=parms->evl.nwvl;
    for(int ievl=info->start; ievl<info->end; ievl++){
	gpu_set(cuglobal->evlgpu[ievl]);
	//info("thread %ld gpu %d ievl %d start\n", thread_id(), cudata->igpu, ievl);
	const int do_psf_cov=(parms->evl.psfmean || parms->evl.psfhist || parms->evl.cov) 
	    && isim>=parms->evl.psfisim && parms->evl.psf->p[ievl]!=0;
	const int save_evlopd=parms->save.evlopd>0 && ((isim+1)%parms->save.evlopd)==0;
	const double thetax=parms->evl.thetax->p[ievl];
	const double thetay=parms->evl.thetay->p[ievl];
 
	stream_t &stream=cuglobal->perf.stream[ievl];
	curmat &iopdevl=cuglobal->perf.opd[ievl];
	// iopdevl must be in device memory. 6 times slower if in host memory.
	if(cuglobal->perf.surf && cuglobal->perf.surf[ievl]){
	    curcp(iopdevl, cuglobal->perf.surf[ievl], stream);
	}else{
	    curset(iopdevl, 0, stream);
	}
	if(parms->sim.idealevl){
	    gpu_dm2loc(iopdevl(), cudata->perf.locs_dm[ievl], cudata->dmproj, parms->ndm,
		       parms->evl.hs->p[ievl], 0, thetax, thetay, 0,0,1, stream);
	}else if(simu->atm && !parms->sim.wfsalias){
	    gpu_atm2loc(iopdevl(), cudata->perf.locs, parms->evl.hs->p[ievl], 0, thetax, thetay, 
			0,0,parms->sim.dt,isim, 1, stream);
	}
	if(simu->telws){//Wind shake 
	    Real tt=simu->telws->p[isim];
	    Real angle=simu->winddir?simu->winddir->p[0]:0;
	    curaddptt(iopdevl, cudata->perf.locs(), 0, tt*cosf(angle), tt*sinf(angle), stream);
	}
	if(simu->telfocusreal){
	    Real focus=-simu->telfocusreal->p[0]->p[0];
	    add_focus_do<<<DIM(nloc, 256),0,stream>>>(iopdevl, cudata->perf.locs(), nloc, focus);
	}
	if(save_evlopd){
	    zfarr_push(simu->save->evlopdol[ievl], isim, iopdevl, stream);
	}
	if(parms->plot.run){
	    drawopdamp_gpu("Evlol", aper->locs, iopdevl, stream, aper->amp1->p, NULL,
			   "Science Open Loop OPD", "x (m)", "y (m)", "OL %d", ievl);
	}
	PERFEVL_WFE_GPU(cuglobal->perf.cc_ol[ievl](), cuglobal->perf.ccb_ol[ievl]);
	if((parms->evl.psfmean  || parms->evl.cov)
	   && isim>=parms->evl.psfisim 
	   &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
	      ||(parms->evl.psfol==2 && parms->evl.psf->p[ievl]))){
	    //calculate Openloop PSF. we also test psfisim to synchronize with psfcl.
	    curmat opdcopy;
	    curmv(cuglobal->perf.coeff[ievl](), 0, cudata->perf.imcc, 
		  cuglobal->perf.cc_ol[ievl](), 'n', 1, stream);
	    curcp(opdcopy, iopdevl, stream);
	    if(parms->evl.pttr->p[ievl]){//remove piston/tip/tilt
		curaddptt(opdcopy, cudata->perf.locs(), cuglobal->perf.coeff[ievl](), -1,-1,-1,stream);
		warning_once("Removing piston/tip/tilt from OPD.\n");
	    }else{//remove piston only
		curaddptt(opdcopy, cudata->perf.locs(), cuglobal->perf.coeff[ievl](), -1, 0, 0, stream);
	    }
	    if(parms->evl.cov){
		if(parms->gpu.psf){
		    curmm(cudata->perf.opdcovol, 1, opdcopy, opdcopy, "nt", 1, stream);
		    curadd(cudata->perf.opdmeanol, 1, opdcopy, 1, stream);
		}else{
		    dmat *tmp=NULL;
		    cp2cpu(&tmp, opdcopy, stream);
		    dmm(&simu->evlopdcovol, 1,tmp, tmp, "nt", 1);
		    dadd(&simu->evlopdmeanol, 1, tmp, 1);
		    dfree(tmp);
		}
	    }
	    if(parms->evl.psfmean){
		psfcomp_r(cudata->perf.psfol(), opdcopy, nwvl, ievl, nloc, parms->evl.psfol==2?1:0, stream);
		if(parms->plot.run){
		    int count=parms->gpu.psf?(simu->perfisim+1-parms->evl.psfisim):1;
		    if(parms->evl.psfol==2){
			count*=lsum(parms->evl.psf);
		    }
		    
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			drawpsf_gpu("PSFol", cudata->perf.psfol[iwvl], count, stream,
				    parms->plot.psf,  "Science Open Loop PSF", 
				    "x", "y", "OL%2d %.2f", ievl, parms->evl.wvl->p[iwvl]*1e6);
		    }
		}
		if(!parms->gpu.psf){ //need to move psf from GPU to CPU for accumulation.
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			add2cpu(&simu->evlpsfolmean->p[iwvl], 1, cudata->perf.psfol[iwvl], 1, stream);
			cuzero(cudata->perf.psfol[iwvl]); //do not accumulate in gpu.
		    }
		}
	    }
	}
	if(parms->sim.evlol) continue;
	if(parms->evl.tomo){
	    TO_IMPLEMENT;
	}else{
	    wait_dmreal(simu, simu->perfisim);
	    gpu_dm2loc(iopdevl(), cudata->perf.locs_dm[ievl], cudata->dmreal, parms->ndm, 
		       parms->evl.hs->p[ievl], 0, thetax, thetay,
		       0,0,-1, stream);
	    if(simu->ttmreal){
		curaddptt(iopdevl, cudata->perf.locs(), 0, -simu->ttmreal->p[0], -simu->ttmreal->p[1], stream);
	    }
	    if(imoao!=-1){
		gpu_dm2loc(iopdevl(), cudata->perf.locs, cudata->dm_evl[ievl], 1,
			   INFINITY, 0, 0, 0, 0, 0, -1, stream);
	    }
	}
	if(save_evlopd){
	    zfarr_push(simu->save->evlopdcl[ievl], isim, iopdevl, stream);
	}

	if(parms->plot.run){
	    drawopdamp_gpu("Evlcl", aper->locs,iopdevl, stream , aper->amp1->p, NULL,
			   "Science Closed loop OPD", "x (m)", "y (m)", "CL %d", ievl);
	}
	PERFEVL_WFE_GPU(cuglobal->perf.cc_cl[ievl](), cuglobal->perf.ccb_cl[ievl]);
	if(do_psf_cov && parms->evl.psfngsr->p[ievl]!=2){//also do normal one.
	    curmv(cuglobal->perf.coeff[ievl](), 0, cudata->perf.imcc, 
		  cuglobal->perf.cc_cl[ievl](), 'n', 1, stream);
	    if(parms->evl.pttr->p[ievl]){
		curaddptt(iopdevl, cudata->perf.locs(), cuglobal->perf.coeff[ievl], -1, -1, -1, stream);
	    }else{
		curaddptt(iopdevl, cudata->perf.locs(), cuglobal->perf.coeff[ievl], -1, 0, 0, stream);
	    }
	    if(parms->evl.cov){
		if(parms->gpu.psf){
		    curmm(cuglobal->perf.opdcov[ievl], 1, iopdevl, iopdevl, "nt", 1, stream);
		    curadd(cuglobal->perf.opdmean[ievl], 1, iopdevl, 1, stream);
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
		    cuccell psfs(nwvl,1);
		    psfcomp(psfs, iopdevl, nwvl, ievl, nloc, stream);
		    zfarr_push(simu->save->evlpsfhist[ievl], isim, psfs, stream);
		    if(parms->evl.psfmean){
			for(int iwvl=0; iwvl<nwvl; iwvl++){
			    curaddcabs2(cuglobal->perf.psfcl[iwvl+nwvl*ievl], 1, 
					psfs[iwvl], 1, stream);
			}
		    }
		}else if(parms->evl.psfmean){
		    psfcomp_r(cuglobal->perf.psfcl+nwvl*ievl, iopdevl, nwvl, ievl, nloc, 0, stream);
		}
		if(parms->plot.run){
		    int count=parms->gpu.psf?(simu->perfisim+1-parms->evl.psfisim):1;
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			drawpsf_gpu("PSFcl",  cuglobal->perf.psfcl[iwvl+nwvl*ievl], count, stream,
				    parms->plot.psf,  "Science Closed Loop PSF", 
				    "x", "y", "CL%2d %.2f", ievl, parms->evl.wvl->p[iwvl]*1e6);
		    }
		}
		if(!parms->gpu.psf){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			add2cpu(&simu->evlpsfmean->p[iwvl+ievl*nwvl], 1, cuglobal->perf.psfcl[iwvl+ievl*nwvl], 1, stream);
			cuzero(cuglobal->perf.psfcl[iwvl+ievl*nwvl]); 
		    }
		}
	    }
	}
	//info("thread %ld gpu %d ievl %d queued\n", thread_id(), cudata->igpu, ievl);
	ctoc("queued");
    }//for ievl
}
void gpu_perfevl_sync(thread_t *info){
    TIC;tic;
    SIM_T *simu=(SIM_T*)info->data;
    const PARMS_T *parms=simu->parms;
    const int isim=simu->perfisim;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int nmod=parms->evl.nmod;
    for(int ievl=info->start; ievl<info->end; ievl++){
	gpu_set(cuglobal->evlgpu[ievl]);
	cudaStream_t stream=cuglobal->perf.stream[ievl];
	const double thetax=parms->evl.thetax->p[ievl];
	const double thetay=parms->evl.thetay->p[ievl];
	/*Setup pointers for easy usage */
	double *polmp=PCOL(simu->olmp->p[ievl], isim);
	double *pclmp=PCOL(simu->clmp->p[ievl], isim);
	double *polep=PCOL(simu->olep->p[ievl], isim);
	double *pclep=PCOL(simu->clep->p[ievl], isim);
	CUDA_SYNC_STREAM;
	PERFEVL_WFE_CPU(ans1, polep, polmp, simu->oleNGSmp, cuglobal->perf.ccb_ol[ievl]);
	PERFEVL_WFE_CPU(ans2, pclep, pclmp, simu->cleNGSmp, cuglobal->perf.ccb_cl[ievl]);
	if(ans1 || ans2){
	    warning("Perfevl fails, redo\n");
	    gpu_perfevl_queue(info);
	    gpu_perfevl_sync(info);
	}
	//info("thread %ld gpu %d ievl %d end\n", thread_id(), cudata->igpu, ievl);
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
    for(int ievl=0; ievl<parms->evl.nevl; ievl++)
	if(parms->evl.psfngsr->p[ievl] && parms->evl.psf->p[ievl])
#pragma omp task
	{

	    //warning_once("Compare with CPU code to verify accuracy. Need to verify focus mode\n");
	    gpu_set(cuglobal->evlgpu[ievl]);
	    curmat &iopdevl=cuglobal->perf.opd[ievl];
	    stream_t &stream=cuglobal->perf.stream[ievl];
	    gpu_ngsmod2science(iopdevl, cudata->perf.locs(), simu->recon->ngsmod, cleNGSm, 
			       parms->evl.thetax->p[ievl], parms->evl.thetay->p[ievl],
			       -1, stream);
	    if(parms->plot.run){
		drawopdamp_gpu("Evlcl", aper->locs,iopdevl, stream , aper->amp1->p, NULL,
			       "Science Closed loop OPD", "x (m)", "y (m)", "ngsr %d", ievl);
	    }
	    if(parms->evl.pttr->p[ievl]){
		calc_ptt(cuglobal->perf.cc_cl[ievl](), cudata->perf.locs(), nloc, iopdevl(), cudata->perf.amp, stream); 
		cudaMemcpyAsync(cuglobal->perf.ccb_cl[ievl], cuglobal->perf.cc_cl[ievl](), 
				4*sizeof(Real), cudaMemcpyDeviceToHost, stream); 
		CUDA_SYNC_STREAM;
		double ptt[3]={0,0,0};
		calc_ptt_post(NULL, ptt,  aper->ipcc, aper->imcc, cuglobal->perf.ccb_cl[ievl]);
		curaddptt(iopdevl, cudata->perf.locs(), -ptt[0], -ptt[1], -ptt[2], stream);
	    }
	    if(parms->evl.cov){
		if(parms->gpu.psf){
		    curmm(cuglobal->perf.opdcov_ngsr[ievl], 1, iopdevl, iopdevl, "nt", 1, stream);
		    curadd(cuglobal->perf.opdmean_ngsr[ievl], 1, iopdevl, 1, stream);
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
		    cuccell psfs(nwvl, 1);
		    psfcomp(psfs, iopdevl, nwvl, ievl, nloc, stream);
		    zfarr_push(simu->save->evlpsfhist_ngsr[ievl], simu->perfisim, psfs, stream);
		    if(parms->evl.psfmean){
			for(int iwvl=0; iwvl<nwvl; iwvl++){
			    curaddcabs2(cuglobal->perf.psfcl_ngsr[iwvl+nwvl*ievl], 1, 
					psfs[iwvl], 1, stream);
			}
		    }
		}else if(parms->evl.psfmean){
		    psfcomp_r(cuglobal->perf.psfcl_ngsr+nwvl*ievl, iopdevl, nwvl, ievl, nloc, 0, stream);
		}
		if(parms->plot.run){
		    int count=parms->gpu.psf?(simu->perfisim+1-parms->evl.psfisim):1;
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			drawpsf_gpu("PSFngsr",  cuglobal->perf.psfcl_ngsr[iwvl+nwvl*ievl], count, stream,
				    parms->plot.psf,  "Science Closed Loop PSF", 
				    "x", "y", "CL%2d %.2f", ievl, parms->evl.wvl->p[iwvl]*1e6);
		    }
		}
		if(!parms->gpu.psf){
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			add2cpu(&simu->evlpsfmean_ngsr->p[iwvl+ievl*nwvl], 1, cuglobal->perf.psfcl_ngsr[iwvl+ievl*nwvl], 1, stream);
			cuzero(cuglobal->perf.psfcl_ngsr[iwvl+ievl*nwvl]); 
		    }
		}
	    }
	    CUDA_SYNC_STREAM;
	}
#pragma omp taskwait
}
void gpu_perfevl_save(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(!parms->evl.nevl) return;
    const int isim=simu->perfisim;
    if(parms->evl.psfmean && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.psfmean)){
	info("Step %d: Output PSF\n", isim);
	const int nwvl=parms->evl.nwvl;
	int nacc=(simu->perfisim+1-parms->evl.psfisim);//total accumulated.
	const double scale=1./(double)nacc;
	if(cudata->perf.psfol){
	    const double scaleol=(parms->evl.psfol==2)?(scale/parms->evl.npsf):(scale);
	    /*copy the PSF accumulated in all the GPUs to CPU.*/
	    X(cell) *temp=X(cellnew)(nwvl, 1);
	    X(cell) *temp2=X(cellnew)(nwvl, 1);
	    for(int im=0; im<NGPU; im++){
		gpu_set(im);
		cp2cpu(&temp2, cudata->perf.psfol, 0);
		cudaStreamSynchronize(0);
		X(celladd)(&temp, 1, temp2, scaleol);
	    }
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(!temp || !temp->p[iwvl]) continue;
		temp->p[iwvl]->header=evl_header(simu->parms, simu->aper, -1, iwvl, isim);
		zfarr_push(simu->save->evlpsfolmean, isim*nwvl+iwvl, temp->p[iwvl]);
		free(temp->p[iwvl]->header); temp->p[iwvl]->header=NULL;
	    }
	    X(cellfree)(temp);
	    X(cellfree)(temp2);
	}
	if(cuglobal->perf.psfcl){
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->evl.psf->p[ievl] || parms->evl.psfngsr->p[ievl]==2) continue;
		gpu_set(cuglobal->evlgpu[ievl]);
		cudaStream_t stream=cuglobal->perf.stream[ievl];
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    curmat &pp=cuglobal->perf.psfcl[iwvl+nwvl*ievl];
		    curscale(pp, scale, stream);
		    pp.header=evl_header(simu->parms, simu->aper, ievl, iwvl, isim);
		    zfarr_push(simu->save->evlpsfmean[ievl], isim*nwvl+iwvl, pp, stream);
		    curscale(pp, 1./scale, stream);
		}
	    }
	}
	if(cuglobal->perf.psfcl_ngsr){
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->evl.psf->p[ievl] || !parms->evl.psfngsr->p[ievl]) continue;
		gpu_set(cuglobal->evlgpu[ievl]);
		cudaStream_t stream=cuglobal->perf.stream[ievl];
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    curmat &pp=cuglobal->perf.psfcl_ngsr[iwvl+nwvl*ievl];
		    curscale(pp, scale, stream);
		    pp.header=evl_header(simu->parms, simu->aper, ievl, iwvl, isim);
		    zfarr_push(simu->save->evlpsfmean_ngsr[ievl], isim*nwvl+iwvl, pp, stream);
		    curscale(pp, 1./scale, stream);
		}
	    }
	}
    }
    if(parms->evl.cov && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.cov)){
	info("Step %d: Output opdcov\n", isim);
	int nacc=(simu->perfisim+1-parms->evl.psfisim);//total accumulated.
	const double scale=1./(double)nacc;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!parms->evl.psf->p[ievl]) continue;
	    gpu_set(cuglobal->evlgpu[ievl]);
	    cudaStream_t stream=cuglobal->perf.stream[ievl];
	    if(parms->evl.psfngsr->p[ievl]!=2){
		{
		    curmat &pp=cuglobal->perf.opdcov[ievl];
		    curscale(pp, scale, stream);
		    zfarr_push(simu->save->evlopdcov[ievl], isim, pp, stream);
		    curscale(pp, 1./scale, stream);
		}
		{
		    curmat &pp=cuglobal->perf.opdmean[ievl];
		    curscale(pp, scale, stream);
		    zfarr_push(simu->save->evlopdmean[ievl], isim, pp, stream);
		    curscale(pp, 1./scale, stream);
		}
	    }
	    if(parms->evl.psfngsr->p[ievl]){
		{
		    curmat &pp=cuglobal->perf.opdcov_ngsr[ievl];
		    curscale(pp, scale, stream);
		    zfarr_push(simu->save->evlopdcov_ngsr[ievl], isim, pp, stream);
		    curscale(pp, 1./scale, stream);
		}
		{
		    curmat &pp=cuglobal->perf.opdmean_ngsr[ievl];
		    curscale(pp, scale, stream);
		    zfarr_push(simu->save->evlopdmean_ngsr[ievl], isim, pp, stream);
		    curscale(pp, 1./scale, stream);
		}
	    }
	}
	if(parms->evl.psfol){
	    const double scaleol=(parms->evl.psfol==2)?(scale/parms->evl.npsf):(scale);
	    {
		X(mat) *temp=NULL;
		X(mat) *temp2=NULL;
		for(int im=0; im<NGPU; im++){
		    gpu_set(im);
		    cp2cpu(&temp2, cudata->perf.opdcovol, 0);
		    cudaStreamSynchronize(0);
		    X(add)(&temp, 1, temp2, scaleol);
		}
		zfarr_push(simu->save->evlopdcovol, isim, temp);
		X(free)(temp);
		X(free)(temp2);
	    }
	    {
		X(mat) *temp=NULL;
		X(mat) *temp2=NULL;
		for(int im=0; im<NGPU; im++){
		    gpu_set(im);
		    cp2cpu(&temp2, cudata->perf.opdmeanol, 0);
		    cudaStreamSynchronize(0);
		    X(add)(&temp, 1, temp2, scaleol);
		}
		zfarr_push(simu->save->evlopdmeanol, isim, temp);
		X(free)(temp);
		X(free)(temp2);
	    }
	}
    }
}
