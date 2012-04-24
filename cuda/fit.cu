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
#include "recon.h"
#include "accphi.h"
#define SYNC_PS  for(int ips=0; ips<recon->npsr; ips++){ cudaStreamSynchronize(curecon->psstream[ips]); }
#define SYNC_FIT  for(int ifit=0; ifit<parms->fit.nfit; ifit++){ cudaStreamSynchronize(curecon->fitstream[ifit]); }
#define SYNC_DM  for(int idm=0; idm<parms->ndm; idm++){ cudaStreamSynchronize(curecon->dmstream[idm]); }

#define TIMING 0
#if !TIMING
#undef TIC
#undef tic
#undef toc
#define TIC
#define tic
#define toc(A)
#define ctoc(A)
#else
#define ctoc(A) CUDA_SYNC_STREAM; toc2(A);tic
#endif
/*
  Todo: share the ground layer which is both matched and same.
*/
/**
   Apply W for fully illuminated points. fully illuminated points are surrounded
   by partially illuminated points. So unexpected wrapover won't happen.  */
__global__ void apply_W_do(float *restrict out, const float *restrict in, const int *W0f, 
				  float alpha, int nx, int n){
    const int step=blockDim.x * gridDim.x;
    for(int k=blockIdx.x * blockDim.x + threadIdx.x; k<n; k+=step){
	int i=W0f[k];
	out[i]+=alpha*(in[i]
		       +0.25f*(in[i-1]+in[i+1]+in[i-nx]+in[i+nx])
		       +0.0625*(in[i-nx-1]+in[i-nx+1]+in[i+nx-1]+in[i+nx+1]));
    }
}
/*
__global__ static void prop_grid_cubic_os2_do(float *restrict out, int nxo,
				  const float *in, int nxi,
				  float dispx, float dispy, float dispx2, float dispy2,
				  int nx, int ny,
				  const float *cc, float alpha){
    (void) dispx2;
    (void) dispy2;
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<ny; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<nx; ix+=stepx){

	}
    }
    }*/
/*
  Raytracing from in to out (dir=1) or in from out (dir=-1)
  in  has origin at (nxi, nyi), sampling dxi, normally DM
  out has origin at (nxo, nyo), sampling dxo, normally other OPD
*/
/*
void gpu_prop_grid_cubic_os2(curmat *in, float oxi, float oyi, float dxi,
			 curmat *out, float oxo, float oyo, float dxo,
			 float dispx, float dispy,  float *cubic_cc,
			 float alpha, int dir, cudaStream_t stream){
    //offset of origin. 
    dispx=(dispx-oxi+oxo);
    dispy=(dispy-oyi+oyo);
    int offx1=0, offy1=0;
    const double dxo1=1./dxo;
    //if output is outside of input from left/bottom, starting at an offset. 
    if(dispx<0){
	offx1=(int)ceilf(-dispx*dxo1);
	dispx+=offx1*dxo;
    }
    if(dispy<0){
	offy1=(int)ceilf(-dispy*dxo1);
	dispy+=offy1*dxo;
    }
    int nxi=in->nx;
    int nyi=in->ny;
    int nxo=out->nx;
    int nyo=out->ny;
    int nx=(int)floorf(((nxi-1)*dxi-dispx)*dxo1)+1;
    int ny=(int)floorf(((nyi-1)*dxi-dispy)*dxo1)+1;
    //if output is outside of input from right/top, limit number of points to do.
    if(nx>nxo-offx1) nx=nxo-offx1;
    if(ny>nyo-offy1) ny=nyo-offy1;
  
    if(fabsf(2*dxo1-dxi)<EPS){
	dispx=dispx/dxi;
	dispy=dispy/dxi;
	int offx2=floor(dispx); dispx-=offx2;
	int offy2=floor(dispy); dispy-=offy2;
	float dispx2=dispx+(dispx>0.5)?-0.5:0.5;
	float dispy2=dispy+(dispy>0.5)?-0.5:0.5;
	TO_IMPLEMENT;
	if(cubic_cc){
	    if(dir==1){
		prop_grid_cubic_do_os2<<<DIM2(nx,ny,8),0,stream>>>
		    (out->p+offy1*nxo+offx1, nxo,
		     in->p+offy2*nxi+offx2, nxi,
		     dispx, dispy, dispx2, dispy2,
		     nx, ny, cubic_cc, alpha);
	    }else if(dir==-1){
		TO_IMPLEMENT;
	    }
	}else{
	    TO_IMPLEMENT;
	}
    }else{
	TO_IMPLEMENT;
    }
    }*/

#define DO_HA(opdfit, xin)	/*opdfit = HA *xin*/			\
    curzero(opdfit->p[ifit], curecon->fitstream[ifit]);			\
    for(int idm=0; idm<recon->ndm; idm++){				\
	const float ht = (float)parms->dm[idm].ht;			\
	const float scale=1.f-ht/hs;					\
	const float dispx=thetax*ht;					\
	const float dispy=thetay*ht;					\
	if(curecon->cubic_cc[idm]){					\
	    gpu_prop_grid_cubic(opdfit->p[ifit], oxp*scale, oyp*scale, dxp*scale, \
				xin->p[idm], recon->amap[idm]->ox, recon->amap[idm]->oy, recon->amap[idm]->dx, \
				dispx, dispy, curecon->cubic_cc[idm], 1.f, 'n', curecon->fitstream[ifit]); \
	}else{								\
	    gpu_prop_grid(opdfit->p[ifit], oxp*scale, oyp*scale, dxp*scale, \
			  xin->p[idm], recon->amap[idm]->ox, recon->amap[idm]->oy, recon->amap[idm]->dx, \
			  dispx, dispy, 1.f, 'n', curecon->fitstream[ifit]); \
	}								\
    }

#define DO_HAT(xout, opdfit2)	/*xout = beta*xout + alpha * HA' *opdfit2 */ \
    for(int idm=0; idm<recon->ndm; idm++){				\
	if(fabsf(beta)<EPS) curzero((*xout)->p[idm], curecon->dmstream[idm]); \
	else if(fabsf(beta-1)>EPS)					\
	    curscale((*xout)->p[idm], beta, curecon->dmstream[idm]);	\
	const float ht = (float)parms->dm[idm].ht;			\
	for(int ifit=0; ifit<nfit; ifit++){				\
	    const float hs = (float)parms->fit.ht[ifit];		\
	    const float scale=1.f-ht/hs;				\
	    const float dispx=(float)parms->fit.thetax[ifit]*ht;	\
	    const float dispy=(float)parms->fit.thetay[ifit]*ht;	\
	    if(curecon->cubic_cc[idm]){					\
		gpu_prop_grid_cubic(opdfit2->p[ifit], oxp*scale, oyp*scale, dxp*scale, \
				    (*xout)->p[idm], recon->amap[idm]->ox, \
				    recon->amap[idm]->oy, recon->amap[idm]->dx, \
				    dispx, dispy, curecon->cubic_cc[idm], \
				    alpha, 't', curecon->dmstream[idm]); \
	    }else{							\
		gpu_prop_grid(opdfit2->p[ifit], oxp*scale, oyp*scale, dxp*scale, \
			      (*xout)->p[idm], recon->amap[idm]->ox,	\
			      recon->amap[idm]->oy, recon->amap[idm]->dx, \
			      dispx, dispy,				\
			      alpha, 't', curecon->dmstream[idm]);	\
	    }								\
	}								\
    }

#define DO_HX(opdfit, xin) /*opdfit = HX * xin */			\
    curzero(opdfit->p[ifit], curecon->fitstream[ifit]);			\
    if(xin){ /*do HX operation, from xin to opdfit.*/			\
	for(int ips=0; ips<npsr; ips++){				\
	    const float ht = (float)recon->ht->p[ips];			\
	    const float scale=1.f-ht/hs;				\
	    assert(xin->p[ips]->nx==recon->xmap[ips]->nx		\
		   &&xin->p[ips]->ny==recon->xmap[ips]->ny);		\
	    gpu_prop_grid(opdfit->p[ifit], oxp*scale, oyp*scale, dxp*scale, \
			  xin->p[ips], recon->xmap[ips]->ox, recon->xmap[ips]->oy, \
			  recon->xmap[ips]->dx, thetax*ht, thetay*ht,	\
			  1.f,'n', curecon->fitstream[ifit]);		\
	}								\
    }else{ /*propagate from atmosphere*/				\
	SIM_T *simu=recon->simu;					\
	gpu_atm2loc(opdfit->p[ifit]->p, curecon->floc, curecon->nfloc, hs, \
		    thetax, thetay, 0, 0, parms->sim.dt*simu->isim,	\
		    1, curecon->fitstream[ifit]);			\
    }

#define DO_HXT(xout, opdfit) /* *xout=*xout*beta + alpha*HX' * opdfit */ \
    for(int ips=0; ips<recon->npsr; ips++){				\
	if(fabsf(beta)<EPS) curzero((*xout)->p[ips], curecon->psstream[ips]); \
	else if(fabsf(beta-1)>EPS)					\
	    curscale((*xout)->p[ips], beta, curecon->psstream[ips]);	\
	const float ht = (float)recon->ht->p[ips];			\
	for(int ifit=0; ifit<nfit; ifit++){				\
	    const float hs = (float)parms->fit.ht[ifit];		\
	    const float scale=1.f-ht/hs;				\
	    float thetax=(float)parms->fit.thetax[ifit];		\
	    float thetay=(float)parms->fit.thetay[ifit];		\
	    gpu_prop_grid(opdfit->p[ifit], oxp*scale, oyp*scale, dxp*scale, \
			  (*xout)->p[ips], recon->xmap[ips]->ox, recon->xmap[ips]->oy, \
			  recon->xmap[ips]->dx, thetax*ht, thetay*ht,	\
			  alpha,'t', curecon->psstream[ips]);		\
	}								\
    }

#define DO_W(opdfit2, opdfit)  /*opdfit2 = (W0-W1*W1')*opdfit */	\
    curzero(opdfit2->p[ifit], curecon->fitstream[ifit]);		\
    cudaMemsetAsync(&curecon->pis->p[ifit], 0, sizeof(float), curecon->fitstream[ifit]); \
    inn_wrap(&curecon->pis->p[ifit], curecon->W01->W1->p,opdfit->p[ifit]->p, np,  curecon->fitstream[ifit]); \
    add_do<<<DIM(np, 256), 0, curecon->fitstream[ifit]>>>		\
	(opdfit2->p[ifit]->p, curecon->W01->W1->p, &curecon->pis->p[ifit], -recon->fitwt->p[ifit], np); \
    cuspmul(opdfit2->p[ifit]->p, curecon->W01->W0p, opdfit->p[ifit]->p,	\
	    recon->fitwt->p[ifit], curecon->fitsphandle[ifit]);		\
    if(curecon->W01->nW0f){						\
	apply_W_do<<<DIM(np, 256), 0, curecon->fitstream[ifit]>>>	\
	    (opdfit2->p[ifit]->p, opdfit->p[ifit]->p, curecon->W01->W0f, \
	     curecon->W01->W0v*recon->fitwt->p[ifit], nxp, curecon->W01->nW0f);	\
    }

/*
  Right hand size operator. 
*/
void gpu_FitR(curcell **xout, float beta, const void *A, const curcell *xin, float alpha){
    TIC;tic;
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=recon->parms;
    const int nfit=parms->fit.nfit;
    const int npsr=recon->npsr;
    float oxp=recon->fmap->ox;
    float oyp=recon->fmap->oy;
    float dxp=recon->fmap->dx;
    const int nxp=recon->fmap->nx;
    const int nyp=recon->fmap->ny;
    const int np=nxp*nyp;
    curcell *opdfit=curecon->opdfit;
    curcell *opdfit2=curecon->opdfit2;
    if(!*xout){
	*xout=curcellnew(recon->ndm, 1, recon->anx, recon->any);
    }
    for(int ifit=0; ifit<nfit; ifit++){
	double hs=parms->fit.ht[ifit];
	float thetax=(float)parms->fit.thetax[ifit];
	float thetay=(float)parms->fit.thetay[ifit];
   	DO_HX(opdfit, xin);    /*opdfit = HX * xin */
	DO_W(opdfit2, opdfit); /*opdfit2 = (W0-W1*W1')*opdfit */
    }
    SYNC_FIT; 
    toc("HX");tic;//1ms
    DO_HAT(xout, opdfit2);/*xout = beta*xout + alpha * HA' *opdfit2 */
    SYNC_DM;   
    toc("HAT");//2.4ms
    /*{
	curcellwrite(xin, "gpur_xin");
	curcellwrite(opdfit, "gpur_opdfit");
	curwrite(curecon->pis, "gpur_pis");
	curcellwrite(opdfit2, "gpur_opdfit2");
	curcellwrite(*xout, "gpur_xout");
	}*/
}
void gpu_FitRt(curcell **xout, float beta, const void *A, const curcell *xin, float alpha){
    TIC;tic;
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=recon->parms;
    const int nfit=parms->fit.nfit;
    float oxp=recon->fmap->ox;
    float oyp=recon->fmap->oy;
    float dxp=recon->fmap->dx;
    const int nxp=recon->fmap->nx;
    const int nyp=recon->fmap->ny;
    const int np=nxp*nyp;
    curcell *opdfit=curecon->opdfit;
    curcell *opdfit2=curecon->opdfit2;
    if(!*xout){
	*xout=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
    }
    for(int ifit=0; ifit<nfit; ifit++){
	float hs    =(float)parms->fit.ht[ifit];
	float thetax=(float)parms->fit.thetax[ifit];
	float thetay=(float)parms->fit.thetay[ifit];
	DO_HA(opdfit,xin); /*do HA operation, opdfit = HA *xin*/
	DO_W(opdfit2,opdfit);  /*do W operation, opdfit2 = (W0-W1*W1')*opdfit*/
    }
    SYNC_FIT; toc("HA");tic;//0.8ms for cubic or linear
    DO_HXT(xout, opdfit2); /*do HXT operation, *xout=*xout*beta + alpha*HX' * opdfit2*/
    SYNC_PS; toc("HXT");
}
void gpu_FitL(curcell **xout, float beta, const void *A, const curcell *xin, float alpha){
    TIC;tic;
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=recon->parms;
    const int nfit=parms->fit.nfit;
    float oxp=recon->fmap->ox;
    float oyp=recon->fmap->oy;
    float dxp=recon->fmap->dx;
    const int nxp=recon->fmap->nx;
    const int nyp=recon->fmap->ny;
    const int np=nxp*nyp;
    curcell *opdfit=curecon->opdfit;
    curcell *opdfit2=curecon->opdfit2;
    if(!*xout){
	*xout=curcellnew(recon->ndm, 1, recon->anx, recon->any);
    }
    for(int ifit=0; ifit<nfit; ifit++){
	float hs    =(float)parms->fit.ht[ifit];
	float thetax=(float)parms->fit.thetax[ifit];
	float thetay=(float)parms->fit.thetay[ifit];
	DO_HA(opdfit,xin); /*do HA operation, opdfit = HA *xin*/
	DO_W(opdfit2,opdfit);  /*do W operation, opdfit2 = (W0-W1*W1')*opdfit*/
    }
    SYNC_FIT;
    toc("HA");tic;//0.8ms for cubic or linear
    /*do HAT operation, from opdfit2 to xout*/
    DO_HAT(xout, opdfit2);/*xout = beta*xout + alpha * HA' *opdfit2 */
#if TIMING
    SYNC_DM; toc("HAT");tic;//2.5ms for cubic, 0.2 ms for linear
#endif
    /*{
	curcellwrite(xin, "gpu_xin");
	curcellwrite(opdfit, "gpu_opdfit");
	curwrite(curecon->pis, "gpu_pis");
	curcellwrite(opdfit2, "gpu_opdfit2");
	curcellwrite(*xout, "gpu_xout");
	}*/
    if(curecon->fitNW){
	curcell *tmp=curcellnew(recon->ndm, 1);
	for(int idm=0; idm<recon->ndm; idm++){
	    tmp->p[idm]=curnew(curecon->fitNW->p[idm]->nx, 1);
	    curmv(tmp->p[idm]->p, 0, curecon->fitNW->p[idm], xin->p[idm]->p, 
		  't', 1, curecon->dmhandle[idm]);
	    curmv((*xout)->p[idm]->p, 1, curecon->fitNW->p[idm], tmp->p[idm]->p, 
		  'n', alpha, curecon->dmhandle[idm]);
	}
	curcellfree(tmp);
    }
    if(curecon->actslave){
	for(int idm=0; idm<recon->ndm; idm++){
	    cuspmul((*xout)->p[idm]->p, curecon->actslave->p[idm*(1+recon->ndm)],
		    xin->p[idm]->p, alpha, curecon->dmsphandle[idm]);
	}
    }
    SYNC_DM;
    toc("fitNW");//0ms.
}
