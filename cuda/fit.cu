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
#include "pcg.h"
#define SYNC_PS  for(int ips=0; ips<recon->npsr; ips++){ curecon->psstream[ips].sync(); }
#define SYNC_FIT  for(int ifit=0; ifit<parms->fit.nfit; ifit++){ curecon->fitstream[ifit].sync(); }
#define SYNC_DM  for(int idm=0; idm<parms->ndm; idm++){ curecon->dmstream[idm].sync(); }

#define TIMING 1

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
	    const float hs = (float)parms->fit.hs[ifit];		\
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
	    const float hs = (float)parms->fit.hs[ifit];		\
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
	    recon->fitwt->p[ifit], curecon->fitstream[ifit]);		\
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
    curecon_t *curecon=cudata->recon;
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
	double hs=parms->fit.hs[ifit];
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
}
void gpu_FitRt(curcell **xout, float beta, const void *A, const curcell *xin, float alpha){
    TIC;tic;
    curecon_t *curecon=cudata->recon;
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
	float hs    =(float)parms->fit.hs[ifit];
	float thetax=(float)parms->fit.thetax[ifit];
	float thetay=(float)parms->fit.thetay[ifit];
	DO_HA(opdfit,xin); /*do HA operation, opdfit = HA *xin*/
	DO_W(opdfit2,opdfit);  /*do W operation, opdfit2 = (W0-W1*W1')*opdfit*/
    }
    SYNC_FIT; toc("HA");tic;//0.8ms for cubic or linear
    DO_HXT(xout, opdfit2); /*do HXT operation, *xout=*xout*beta + alpha*HX' * opdfit2*/
    SYNC_PS; toc("HXT");
}
void gpu_FitL(curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream){
    TIC;tic;
    curecon_t *curecon=cudata->recon;
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
    CUDA_SYNC_STREAM;
    for(int ifit=0; ifit<nfit; ifit++){
	float hs    =(float)parms->fit.hs[ifit];
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
  
    if(curecon->fitNW){
	curcell *tmp=curcellnew(recon->ndm, 1);
	for(int idm=0; idm<recon->ndm; idm++){
	    tmp->p[idm]=curnew(curecon->fitNW->p[idm]->nx, 1);
	    curmv(tmp->p[idm]->p, 0, curecon->fitNW->p[idm], xin->p[idm]->p, 
		  't', 1, curecon->dmstream[idm]);
	    curmv((*xout)->p[idm]->p, 1, curecon->fitNW->p[idm], tmp->p[idm]->p, 
		  'n', alpha, curecon->dmstream[idm]);
	}
	curcellfree(tmp);
    }
    if(curecon->actslave){
	for(int idm=0; idm<recon->ndm; idm++){
	    cuspmul((*xout)->p[idm]->p, curecon->actslave->p[idm*(1+recon->ndm)],
		    xin->p[idm]->p, alpha, curecon->dmstream[idm]);
	}
    }
    SYNC_DM;
    toc("fitNW");//0ms.
}
/**
   Wrap of the DM fitting operation

   opdr is the OPD input.
   fitx is the right hand side vector computed from opdr. Allow NULL.
   fitr is the DM fitting result.
*/
double gpu_fit_do(const PARMS_T *parms,const RECON_T *recon, curcell *fitr, curcell *fitx, curcell *opdr, stream_t &stream){
    G_CGFUN cg_fun;
    void *cg_data;
    curcell *rhs=NULL;
    double res=0;
    curecon_t *curecon=cudata->recon;
    if(fitx){
	rhs=fitx;
    }
#if TIMING 
    EVENT_INIT(3);
    EVENT_TIC(0);
#endif
    curmat *tmp=NULL;
    if(parms->gpu.fit==1){//sparse matrix
	cumuv(&rhs, 0, &curecon->FR, opdr, 1);
	cg_fun=(G_CGFUN) cumuv;
	cg_data=&curecon->FL;
    }else{
	gpu_FitR(&rhs, 0, recon, opdr, 1);
	cg_fun=(G_CGFUN) gpu_FitL;
	cg_data=(void*)recon;
    }
#if TIMING 
    EVENT_TIC(1);
#endif
    switch(parms->fit.alg){
    case 0:
	cuchol_solve(fitr->m->p, curecon->FCl, curecon->FCp, rhs->m->p, stream);
	if(curecon->FUp){
	    tmp=curnew(curecon->FVp->ny, 1);
	    curmv(tmp->p, 0, curecon->FVp, rhs->m->p, 't', -1, stream);
	    curmv(fitr->m->p, 1, curecon->FUp, tmp->p, 'n', 1, stream);
	}
	break;
    case 1:{
	if((res=gpu_pcg(&fitr, (G_CGFUN)cg_fun, cg_data, NULL, NULL, rhs, &curecon->cgtmp_fit,
			parms->recon.warm_restart, parms->fit.maxit, stream))>1){
	    warning("DM Fitting PCG not converge. res=%g\n", res);
	}
    }
	break;
    case 2:
	curmv(fitr->m->p, 0, curecon->FMI, rhs->m->p, 'n', 1, stream);
	break;
    default:
	error("Invalid");
    }
#if TIMING 
    EVENT_TIC(2);
    EVENT_TOC;
    info2("Fit RHS: %6.0f, LHS: %6.0f, Tot: %5.0f\n", times[1], times[2], times[0]);
#endif
    if(!fitx){
	curcellfree(rhs);
    }
    if(tmp) curfree(tmp);
    return res;
}

void gpu_fit_test(SIM_T *simu){	/*Debugging. */
    gpu_set(gpu_recon);
    stream_t stream;
    curecon_t *curecon=cudata->recon;
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
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
    gpu_FitL(&lg, 0, recon, rhsg, 1, stream);
    curcellwrite(lg, "GPU_FitL");
    gpu_FitL(&lg, 1, recon, rhsg, -1, stream);
    curcellwrite(lg, "GPU_FitL2");
    gpu_FitL(&lg, 0, recon, rhsg, 1, stream);
    curcellwrite(lg, "GPU_FitL3");
    /*curcell *lhsg=NULL;
      gpu_FitRt(&lhsg, 0, recon, rhsg, 1);
      curcellwrite(lhsg, "GPU_FitRt");*/
    curcellzero(lg, curecon->cgstream[0]);
    gpu_pcg(&lg, (G_CGFUN)gpu_FitL, (void*)recon, NULL, NULL, rhsg, &curecon->cgtmp_fit,
	    simu->parms->recon.warm_restart, parms->fit.maxit, curecon->cgstream[0]);
    curcellwrite(lg, "GPU_FitCG");
    CUDA_SYNC_DEVICE;
    exit(0);
}
