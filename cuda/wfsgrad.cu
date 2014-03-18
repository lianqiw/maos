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
#include "../maos/sim.h"
}
#include "utils.h"
#include "accphi.h"
#include <curand_kernel.h>
#include <cusparse.h>
#include <cufft.h>
#include "wfs.h"
#include "cudata.h"
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
#define ctoc(A) CUDA_SYNC_STREAM; toc2(A)
#endif
extern const char *dirskysim;
/*
  Notice that both blocks and threads are partitioning isa
 */
__global__ static void add_geom_noise_do(float *restrict g, const float *restrict nea, 
				      int nsa, curandState *restrict rstat){
    const int id=threadIdx.x + blockIdx.x * blockDim.x;
    curandState lstat=rstat[id];
    const int nstep=blockDim.x * gridDim.x;
    for(int i=id; i<nsa; i+=nstep){
	float n1=curand_normal(&lstat);
	float n2=curand_normal(&lstat);
	g[i]+=n1*nea[i];
	g[i+nsa]+=n2*nea[i+nsa]+n1*nea[i+nsa*2];/*cross term. */
    }
    rstat[id]=lstat;
}

/**
   Compute ztilt over a square subaperture.
*/
static __global__ void 
cuztilt_do(float *restrict g, float *restrict opd, 
	   const int nsa, const float dx, const int nx, float (**imcc)[3],
	   const float (*orig)[2], const float*restrict amp, float alpha){
    extern __shared__ float a0[];
    int idx=threadIdx.x+threadIdx.y*blockDim.x;
    float *a[3];
    for(int i=0; i<3; i++){
	a[i]=a0+blockDim.x*blockDim.y*i;
	a[i][idx]=0;
    }
    const int isa=blockIdx.x;
    const int skip=isa*nx*nx;
    const float ox=orig[isa][0];
    const float oy=orig[isa][1];
    for(int iy=threadIdx.y; iy<nx; iy+=blockDim.y){
	const int skip2=skip+iy*nx;
	const float y=iy*dx+oy;
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    const int ind=skip2+ix;
	    const float tmp=amp[ind]*opd[ind];
	    a[0][idx]+=tmp;
	    a[1][idx]+=tmp*(dx*ix+ox);
	    a[2][idx]+=tmp*y;
	}
    }
    for(int step=(blockDim.x*blockDim.y)>>1;step>0;step>>=1){
	__syncthreads();
	if(idx<step){
	    for(int i=0; i<3; i++){
		a[i][idx]+=a[i][idx+step];
	    }
	}
    }
    __syncthreads();
    if(threadIdx.x<2 && threadIdx.y==0){
	float (*restrict A)[3]=imcc[isa];
	atomicAdd(&g[isa+threadIdx.x*nsa], 
		  alpha*(+a[0][0]*A[0][threadIdx.x+1]
			 +a[1][0]*A[1][threadIdx.x+1]
			 +a[2][0]*A[2][threadIdx.x+1]));
	
    }
}
void cuztilt(float *restrict g, float *restrict opd, 
	     const int nsa, const float dx, const int nx, float (**imcc)[3],
	     const float (*orig)[2], const float*restrict amp, float alpha, cudaStream_t stream){
    const int tx=4;
    cuztilt_do<<<nsa, dim3(tx,tx), tx*tx*3*sizeof(float), stream>>>
	(g, opd, nsa, dx, nx, imcc, orig, amp, alpha);
    CUDA_SYNC_STREAM;
}
/**
   Apply matched filter. \todo this implementation relies on shared variable. It
is probably causing competition.  */
__global__ static void mtche_do(float *restrict grad, float (*restrict *restrict mtches)[2], 
				const float *restrict ints, const float *restrict i0sum, 
				int pixpsa, int nsa){
    extern __shared__ float g0[];
    float *g[3];
    for(int i=0; i<3; i++){
	g[i]=g0+blockDim.x*i;
	g[i][threadIdx.x]=0;
    }
    int isa=blockIdx.x;
    ints+=isa*pixpsa;
    const float (*const restrict mtche)[2]=mtches[isa];
 
    for (int ipix=threadIdx.x; ipix<pixpsa; ipix+=blockDim.x){
	g[0][threadIdx.x]+=mtche[ipix][0]*ints[ipix];
	g[1][threadIdx.x]+=mtche[ipix][1]*ints[ipix];
	g[2][threadIdx.x]+=ints[ipix];
    }
    for(int step=(blockDim.x)>>1;step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    for(int i=0; i<3; i++){
		g[i][threadIdx.x]+=g[i][threadIdx.x+step];
	    }
	}
    }
    __syncthreads();
    if(threadIdx.x<2){
	if(i0sum){/*normalize gradients according to siglev.*/
	    g[threadIdx.x][0]*=i0sum[isa]/g[2][0];
	}
	grad[isa+nsa*threadIdx.x]=g[threadIdx.x][0];
    }
}

static void mtche(float *restrict grad, float (*restrict *restrict mtches)[2], 
		  const float *restrict ints, const float *restrict i0sum, 
		  int pixpsa, int nsa, int msa, cudaStream_t stream){
    for(int isa=0; isa<nsa; isa+=msa){
	int ksa=MIN(msa, nsa-isa);
	mtche_do<<<ksa, 16, ksa*sizeof(float)*3, stream>>>
	    (grad+isa, mtches+isa, ints+pixpsa*isa, i0sum?i0sum+isa:0, pixpsa, nsa);
    }
}
/**
   Apply tCoG.
*/
__global__ static void tcog_do(float *grad, const float *restrict ints, 
			       int nx, int ny, float pixthetax, float pixthetay, int nsa, float (*cogcoeff)[2], float *srot){
    __shared__ float sum[3];
    if(threadIdx.x<3 && threadIdx.y==0) sum[threadIdx.x]=0.f;
    __syncthreads();//is this necessary?
    int isa=blockIdx.x;
    ints+=isa*nx*ny;
    float cogthres=cogcoeff[isa][0];
    float cogoff=cogcoeff[isa][1];
    for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    float im=ints[ix+iy*nx]-cogoff;
	    if(im>cogthres){
		atomicAdd(&sum[0], im);
		atomicAdd(&sum[1], im*ix);
		atomicAdd(&sum[2], im*iy);
	    }
	}
    }
    __syncthreads();
    if(threadIdx.x==0 && threadIdx.y==0){
	if(fabsf(sum[0])>0){
	    float gx=(sum[1]/sum[0]-(nx-1)*0.5)*pixthetax;
	    float gy=(sum[2]/sum[0]-(ny-1)*0.5)*pixthetay;
	    if(srot){
		float s,c;
		sincos(srot[isa], &s, &c);
		float tmp=gx*c-gy*s;
		gy=gx*s+gy*c;
		gx=tmp;
	    }
	    grad[isa]=gx;
	    grad[isa+nsa]=gy;
	}else{
	    grad[isa]=0;
	    grad[isa+nsa]=0;
	}
    }
}
/**
   Poisson random generator.
*/
__device__ static float curandp(curandState *rstat, float xm){
    float g, t, xmu;
    int x=0, xu;
    if(xm>200){
	x=(int)round(xm+curand_normal(rstat)*sqrt(xm));
    }else{
	while(xm>0){
	    xmu = xm > 12.f ? 12.f : xm;
	    xm -= xmu;
	    g   = __expf(-xmu);
	    xu  = -1;
	    t   = 1.f;
	    while(t>g){
		xu++;
		t *= curand_uniform(rstat);
	    }
	    x += xu;
	}
    }
    return x;
}
/**
   Add noise to pix images.
*/
__global__ static void addnoise_do(float *restrict ints0, int nsa, int pixpsa, float bkgrnd, float bkgrndc, 
				   float *const restrict *restrict bkgrnd2s, float *const restrict *restrict bkgrnd2cs,
				   float rne, curandState *rstat){
    const int id=threadIdx.x + blockIdx.x * blockDim.x;
    const int nstep=blockDim.x * gridDim.x;
    curandState lstat=rstat[id];
    for(int isa=id; isa<nsa; isa+=nstep){
	float *restrict ints=ints0+isa*pixpsa;
	const float *restrict bkgrnd2=bkgrnd2s?bkgrnd2s[isa]:NULL;
	const float *restrict bkgrnd2c=bkgrnd2cs?bkgrnd2cs[isa]:NULL;
	for(int ipix=0; ipix<pixpsa; ipix++){
	    float corr=bkgrnd2c?(bkgrnd2c[ipix]+bkgrndc):bkgrndc;
	    if(bkgrnd2){
		ints[ipix]=curandp(&lstat, ints[ipix]+bkgrnd+bkgrnd2[ipix])+rne*curand_normal(&lstat)-corr;
	    }else{
		ints[ipix]=curandp(&lstat, ints[ipix]+bkgrnd)+rne*curand_normal(&lstat)-corr;
	    }
	}
    }
    rstat[id]=lstat;
}
void gpu_fieldstop(curmat *opd, float *amp, int *embed, int nembed, 
		   curmat* fieldstop, float wvl, cufftHandle fftplan, cudaStream_t stream){
    cucmat wvf(nembed, nembed);
    embed_wvf_do<<<DIM(opd->nx, 256), 0, stream>>>
	(wvf.p, opd->p, amp, embed, opd->nx, wvl);
    CUFFT(fftplan, wvf.p, CUFFT_FORWARD);
    cwm_do<<<DIM(wvf.nx*wvf.ny, 256),0,stream>>>
      (wvf.p, fieldstop->p, wvf.nx*wvf.ny);
    CUFFT(fftplan, wvf.p, CUFFT_INVERSE);
    unwrap_phase_do<<<DIM2(wvf.nx, wvf.ny, 16),0,stream>>>
	(wvf.p, opd->p, embed, opd->nx, wvl);
}
__global__ static void
dither_acc_do(float *restrict *im0, float *restrict *imx, float *restrict *imy, 
	      float *restrict const *pints, float cd, float sd, int pixpsa, int nsa){
    for(int isa=blockIdx.x; isa<nsa; isa+=gridDim.x){
	const float *ints=pints[isa];
	float *restrict acc_ints=im0[isa];
	float *restrict acc_intsx=imx[isa];
	float *restrict acc_intsy=imy[isa];
	for(int ipix=threadIdx.x; ipix<pixpsa; ipix+=blockDim.x){
	    acc_ints[ipix]+=ints[ipix];
	    acc_intsx[ipix]+=ints[ipix]*cd;
	    acc_intsy[ipix]+=ints[ipix]*sd;
	}
    }
}
dither_t::dither_t(int nsa, int pixpsax, int pixpsay):imc(0){
    im0=curcellnew(nsa,1,pixpsax,pixpsay);
    imx=curcellnew(nsa,1,pixpsax,pixpsay);
    imy=curcellnew(nsa,1,pixpsax,pixpsay);
}
void dither_t::reset(){
    imc=0;
    curcellzero(im0);
    curcellzero(imx);
    curcellzero(imy);
}
/**Accumulate for matched filter updating*/
void dither_t::acc(curcell *ints, float angle, cudaStream_t stream){
    const int nsa=ints->nx*ints->ny;
    const int pixpsa=ints->p[0]->nx*ints->p[0]->ny;
    dither_acc_do<<<ints->nx, pixpsa, 0, stream>>>
	(im0->pm, imx->pm, imy->pm, ints->pm, cosf(angle), sinf(angle), pixpsa, nsa);
    imc++;
}
/**Output for matched filter updating*/
void dither_t::output(float a2m, int iwfs, int isim, cudaStream_t stream){
    curcellscale(im0, 1./(imc), stream);
    curcellscale(imx, 2./(a2m*imc), stream);
    curcellscale(imy, 2./(a2m*imc), stream);
    CUDA_SYNC_STREAM;
    curcellwrite(im0, "wfs%d_i0_%d", iwfs, isim);
    curcellwrite(imx, "wfs%d_gx_%d", iwfs, isim);
    curcellwrite(imy, "wfs%d_gy_%d", iwfs, isim);
    reset();
}

/**
   Ray tracing and gradient computation for WFS. \todo Expand to do gradients in GPU without transfering
   data back to CPU.
*/
void gpu_wfsgrad_iwfs(SIM_T *simu, int iwfs){
    gpu_set(cudata_t::wfsgpu[iwfs]);
    cuwloc_t *cupowfs=cudata->powfs;
    cuwfs_t *cuwfs=cudata_t::wfs;
    TIC;tic;
    const PARMS_T *parms=simu->parms;
    assert(iwfs<parms->nwfs);
    const POWFS_T *powfs=simu->powfs;
    const RECON_T *recon=simu->recon;
    /*output */
    const int CL=parms->sim.closeloop;
    const int isim=simu->isim;
    /*The following are truly constants for this powfs */
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int imoao=parms->powfs[ipowfs].moao;
    const int nsa=powfs[ipowfs].pts->nsa;
    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
    const float hs=parms->powfs[ipowfs].hs;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const int save_gradgeom=parms->save.gradgeom[iwfs];
    const int save_opd =parms->save.wfsopd[iwfs];
    const int save_ints=parms->save.ints[iwfs];
    const int noisy=parms->powfs[ipowfs].noisy;
    /*The following depends on isim */
    const int dtrat_output=((isim+1)%dtrat==0);
    const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
    const int do_pistatout=parms->powfs[ipowfs].pistatout&&isim>=parms->powfs[ipowfs].pistatstart;
    const int do_geom=!do_phy || save_gradgeom || do_pistatout;
    const float thetax=parms->wfs[iwfs].thetax;
    const float thetay=parms->wfs[iwfs].thetay;
    const float dtisim=parms->sim.dt*isim;
    float (*loc)[2]=cupowfs[ipowfs].loc->p;
    const int nloc=cupowfs[ipowfs].loc->nloc;
    /*Out to host for now. \todo : keep grad in device when do reconstruction on device. */
    stream_t &stream=*cuwfs[iwfs].stream;
    dmat *gradcl=simu->gradcl->p[iwfs];
    curmat *phiout=cuwfs[iwfs].phiout;
    curmat *gradacc=cuwfs[iwfs].gradacc;
    curmat *gradcalc=cuwfs[iwfs].gradcalc;
    if(cuwfs[iwfs].opdadd){ /*copy to phiout. */
	curcp(&phiout, cuwfs[iwfs].opdadd, stream);
    }else{
	curzero(phiout, stream);
    }
    curzero(cuwfs[iwfs].gradcalc, stream);
    if(parms->sim.idealwfs){
	gpu_dm2loc(phiout->p, cuwfs[iwfs].loc_dm, cudata->dmproj, cudata->ndm,
		   hs, thetax, thetay, 0, 0, 1, stream);
    }else{
	if(simu->atm){
	    gpu_atm2loc(phiout->p, cuwfs[iwfs].loc_tel, hs, thetax, thetay, 0, 0, dtisim, 1, stream);
	}
	if(parms->sim.wfsalias){
	    gpu_dm2loc(phiout->p, cuwfs[iwfs].loc_dm, cudata->dmproj, cudata->ndm,
		       hs, thetax, thetay, 0, 0, -1, stream);
	}
    }
    if(simu->telws){
	float tt=simu->telws->p[isim];
	float angle=simu->winddir?simu->winddir->p[0]:0;
	curaddptt(phiout, loc, 0, tt*cosf(angle), tt*sinf(angle), stream);
    }
    if(save_opd){
	cellarr_cur(simu->save->wfsopdol[iwfs], simu->isim, phiout, stream);
    }
    if(CL){
	gpu_dm2loc(phiout->p, cuwfs[iwfs].loc_dm, cudata->dmreal, cudata->ndm,
		   hs, thetax, thetay, 0, 0, -1, stream);
	if(simu->ttmreal){
	    curaddptt(phiout, loc, 0, -simu->ttmreal->p[0], -simu->ttmreal->p[1], stream);
	}
    }
    if(parms->tomo.ahst_idealngs && parms->powfs[ipowfs].skip){
	const double *cleNGSm=simu->cleNGSm->p+isim*recon->ngsmod->nmod;
	gpu_ngsmod2science(phiout, cupowfs[ipowfs].loc->p, recon->ngsmod, cleNGSm, 
			   parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay, 
			   -1, stream);
    }
    /*CUDA_SYNC_STREAM; */
    if(imoao>-1){
	gpu_dm2loc(phiout->p, cuwfs[iwfs].loc_dm, cudata->dm_wfs[iwfs], 1,
		   INFINITY, 0, 0, 0, 0, -1, stream);
    }
  
    if(parms->powfs[ipowfs].llt){
	float focus=(float)wfsfocusadj(simu, iwfs);
	if(fabsf(focus)>1e-20){
	    add_focus_do<<<DIM(nloc, 256), 0, stream>>>(phiout->p, loc, nloc, focus);
	}
    }
    if(parms->powfs[ipowfs].fieldstop){
	gpu_fieldstop(phiout, cuwfs[iwfs].amp, cupowfs[ipowfs].embed, cupowfs[ipowfs].nembed, 
		      cupowfs[ipowfs].fieldstop, parms->powfs[ipowfs].wvl[0], cuwfs[iwfs].plan_fs, stream);
    }
    if(save_opd){
	cellarr_cur(simu->save->wfsopd[iwfs], simu->isim, phiout, stream);
    }
    if(parms->plot.run>1){
	const double *realamp=powfs[ipowfs].realamp->p[wfsind]->p;
	dmat *tmp=NULL;
	cp2cpu(&tmp, phiout, stream);
	drawopdamp("wfsopd",powfs[ipowfs].loc,tmp->p,realamp,NULL,
		   "WFS OPD","x (m)", "y (m)", "WFS %d", iwfs);
	dfree(tmp);
    }
    curmat *gradref=0;
    if(do_geom){
	double ratio;
	if(do_pistatout && dtrat>1){
	    gradref=gradcalc; 
	    ratio=1;
	}else{
	    gradref=gradacc;
	    ratio=1.f/(float)dtrat;
	}

	if(parms->powfs[ipowfs].gtype_sim==1){
	    cuztilt(gradref->p, phiout->p, 
		    cupowfs[ipowfs].pts->nloc, 
		    cupowfs[ipowfs].pts->dxsa, 
		    cupowfs[ipowfs].pts->nxsa, cuwfs[iwfs].imcc,
		    cupowfs[ipowfs].pts->p, cuwfs[iwfs].amp, ratio, stream);
	}else{
	    cusp *GS0=cuwfs[iwfs].GS0;
	    cuspmul(gradref->p, GS0, phiout->p, 1, 'n', ratio, stream);
	}
	if(gradacc!=gradref){
	    curadd(&gradacc, 1, gradref, 1.f/(float)dtrat, stream);
	}
    }   
    if(parms->powfs[ipowfs].psfout){
	cellarr_cur(simu->save->ztiltout[iwfs], simu->isim, gradcalc, stream);
    }
    if(do_phy || parms->powfs[ipowfs].psfout || do_pistatout){/*physical optics */
	gpu_wfsints(simu, phiout->p, gradref, iwfs, isim, stream);
    }/*do phy */
    ctoc("grad");
    if(dtrat_output){
	if(do_phy){
	    /*signal level was already multiplied in ints. */
	    curcell *ints=cuwfs[iwfs].ints;
	    const int pixpsa=powfs[ipowfs].pixpsax*powfs[ipowfs].pixpsay;
	    if(save_ints){
		cellarr_curcell(simu->save->intsnf[iwfs], simu->isim, ints, stream);
	    }
	    /*CUDA_SYNC_STREAM; */
	    ctoc("mtche");
	    float rne=0, bkgrnd=0;
	    if(noisy){
		rne=parms->powfs[ipowfs].rne;
		bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
		addnoise_do<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt, 0, stream>>>
		    (ints->p[0]->p, nsa, pixpsa, bkgrnd, bkgrnd*parms->powfs[ipowfs].bkgrndc,
		     cuwfs[iwfs].bkgrnd2, cuwfs[iwfs].bkgrnd2c, 
		     rne, cuwfs[iwfs].custat);
		ctoc("noise");
	    }
	    if(parms->powfs[ipowfs].dither && isim>=parms->powfs[ipowfs].dither_nskip){
		double angle=M_PI*0.5*isim/parms->powfs[ipowfs].dtrat;
		angle+=simu->dither[iwfs]->deltam;
		cuwfs[iwfs].dither->acc(ints, angle, stream);
		int nstat=parms->powfs[ipowfs].dither_nstat;
		int dtrat=parms->powfs[ipowfs].dtrat;
		if((isim-parms->powfs[ipowfs].dither_nskip+1)%(nstat*dtrat)==0){
		    warning2("Dither step%d, wfs%d: output statistics\n", isim, iwfs);
		    cuwfs[iwfs].dither->output((float)simu->dither[iwfs]->a2m, iwfs, isim, stream);
		}
	    }
	    curzero(gradcalc, stream);
	    switch(parms->powfs[ipowfs].phytypesim){
	    case 1:
		mtche(gradcalc->p, cuwfs[iwfs].mtche, ints->p[0]->p, 
		      parms->powfs[ipowfs].mtchscl?cuwfs[iwfs].i0sum:NULL,
		      pixpsa, nsa, cuwfs[iwfs].msa, stream);
		break;
	    case 2:{
		float pixthetax=(float)parms->powfs[ipowfs].radpixtheta;
		float pixthetay=(float)parms->powfs[ipowfs].pixtheta;
		int pixpsax=powfs[ipowfs].pixpsax;
		int pixpsay=powfs[ipowfs].pixpsay;
		float *srot=parms->powfs[ipowfs].radpix?cuwfs[iwfs].srot:NULL;
		tcog_do<<<nsa, dim3(pixpsax, pixpsay),0,stream>>>
		    (gradcalc->p, ints->p[0]->p, 
		     pixpsax, pixpsay, pixthetax, pixthetay, nsa, (float(*)[2])cuwfs[iwfs].cogcoeff, srot);
	    }
		break;
	    case 3:{
		dcell *cints=NULL;
		cp2cpu(&cints, ints, stream);
		CUDA_SYNC_STREAM;
		double geach[3];
		for(int isa=0; isa<nsa; isa++){
		    geach[0]=gradcl->p[isa];
		    geach[1]=gradcl->p[isa+nsa];
		    geach[2]=1;
		    maxapriori(geach, cints->p[isa], parms, powfs, iwfs, isa, 1, bkgrnd, rne);
		    gradcl->p[isa]=geach[0];
		    gradcl->p[isa+nsa]=geach[1];
		}
		dcellfree(cints);
	    }
		break;
	    default:
		TO_IMPLEMENT;
	    }

	    if(save_ints && noisy){
		cellarr_curcell(simu->save->intsny[iwfs], simu->isim, ints, stream);
	    }
	    ctoc("mtche");
	    /*send grad to CPU. */
	    if(parms->powfs[ipowfs].phytypesim!=3){//3 is handled in cpu.
		cp2cpu(&gradcl, gradcalc, stream);
	    }
	    ctoc("dev2dbl");
	    curcellzero(ints, stream);
	    if(save_gradgeom){//also do geom grad during phy grad sims
		CUDA_SYNC_STREAM;/*necessary. */
		cellarr_cur(simu->save->gradgeom[iwfs], simu->isim, gradacc, stream);
	    }
	}else{
	    if(noisy && !parms->powfs[ipowfs].usephy){
		add_geom_noise_do<<<cuwfs[iwfs].custatb, cuwfs[iwfs].custatt, 0, stream>>>
		    (gradacc->p, cuwfs[iwfs].neasim, nsa,cuwfs[iwfs].custat);
		ctoc("noise");
	    }
	    cp2cpu(&gradcl, gradacc, stream);
	    if(save_gradgeom){
		cellarr_cur(simu->save->gradgeom[iwfs], simu->isim, NULL, stream);
	    }
	    ctoc("dev2dbl");
	}
	if(do_geom){
	    curzero(gradacc, stream);
	}
    }/*dtrat_output */
    ctoc("done");
    CUDA_SYNC_STREAM;
}
void gpu_wfsgrad_save(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    if((isim % 50 ==0) || isim+1==parms->sim.end){
	for(int iwfs=0; iwfs<simu->parms->nwfs; iwfs++){
	    gpu_set(cudata_t::wfsgpu[iwfs]);
	    cuwfs_t *cuwfs=cudata_t::wfs;
	    const PARMS_T *parms=simu->parms;
	    const int ipowfs=simu->parms->wfs[iwfs].powfs;
	    stream_t &stream=*cuwfs[iwfs].stream;
	    if(cuwfs[iwfs].pistatout){
		int nstep=isim+1-parms->powfs[ipowfs].pistatstart;
		if(nstep>0){
		    curcell* tmp=cuwfs[iwfs].pistatout;
		    curcellscale(tmp, 1.f/(float)nstep, stream);
		    if(parms->sim.skysim){
			curcellwrite(tmp, "%s/pistat/pistat_seed%d_sa%d_x%g_y%g.bin",
				     dirskysim,simu->seed,
				     parms->powfs[ipowfs].order,
				     parms->wfs[iwfs].thetax*206265,
				     parms->wfs[iwfs].thetay*206265);
		    }else{
			curcellwrite(tmp,"pistat_seed%d_wfs%d.bin", simu->seed,iwfs);
		    }
		    curcellscale(tmp, nstep, stream);
		}
	    }
	}
    }
}
