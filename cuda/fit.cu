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
#include "recon.h"
#include "accphi.h"
#include "pcg.h"

#define TIMING 0
#if TIMING <1
#undef EVENT_INIT
#undef EVENT_TIC
#undef EVENT_TOC
#undef EVENT_PRINT
#define EVENT_INIT(A)
#define EVENT_TIC(A)
#define EVENT_TOC
#define EVENT_PRINT(A...)
#else
#define EVENT_PRINT(A...) fprintf(stderr, A);EVENT_DEINIT
#endif

#if TIMING <2
#define EVENT2_INIT(A)
#define EVENT2_TIC(A)
#define EVENT2_TOC
#define EVENT2_PRINT(A...)
#else
#define EVENT2_INIT EVENT_INIT
#define EVENT2_TIC EVENT_TIC
#define EVENT2_TOC EVENT_TOC
#define EVENT2_PRINT(A...) fprintf(stderr, A);EVENT_DEINIT
#endif


/*
  Todo: share the ground layer which is both matched and same.
*/
/**
   Apply W for fully illuminated points. fully illuminated points are surrounded
   by partially illuminated points. So unexpected wrapover won't happen.  */
__global__ void 
apply_W_do(float *outs, const float *ins, const int *W0f, float alpha, int nx, int n){
    float *out=outs+blockIdx.y*nx*nx;
    const float *in=ins+blockIdx.y*nx*nx;
    const int step=blockDim.x * gridDim.x;
    for(int k=blockIdx.x * blockDim.x + threadIdx.x; k<n; k+=step){
	int i=W0f[k];
	out[i]+=alpha*(in[i]
		       +0.25f*(in[i-1]+in[i+1]+in[i-nx]+in[i+nx])
		       +0.0625*(in[i-nx-1]+in[i-nx+1]+in[i+nx-1]+in[i+nx+1]));
    }
}
/**
do HXp operation, opdfit+=Hxp*xin*/
static inline void
fit_do_hxp(curecon_t *curecon, const curcell *xin, stream_t &stream){
    const int nfit=curecon->opdfit->nx;
    const int npsr=xin->nx;
    EVENT2_INIT(4);
    EVENT2_TIC(0);
    curecon->opdfit->m->zero(stream);
    EVENT2_TIC(1);
    if(curecon->xcache){
	curecon->xcache->m->zero(stream);
	gpu_prop_grid_do<<<dim3(3,3,npsr),dim3(16,16),0,stream>>>
	    (curecon->hxp0data, curecon->xcache->pm, xin->pm, 1, npsr, 1, NULL, 'n');
	EVENT2_TIC(2);
	gpu_prop_grid_do<<<dim3(3,3,nfit),dim3(16,16),0,stream>>>
	    (curecon->hxp1data, curecon->opdfit->pm, curecon->xcache->pm, nfit, npsr, 1, NULL, 'n');
	EVENT2_TIC(3);
    }else{
	EVENT2_TIC(2);
	gpu_prop_grid_do<<<dim3(3,3,nfit),dim3(16,16),0,stream>>>
	    (curecon->hxpdata, curecon->opdfit->pm, xin->pm, nfit, npsr, 1, NULL, 'n');
	EVENT2_TIC(3);
    }
    EVENT2_TOC;
    EVENT2_PRINT("do_hxp: zero: %.3f cache: %.3f prop: %.3f tot: %.3f\n", times[1], times[2], times[3], times[0]);
}
/**
do HXp' operation, xout+=alpha*Hxp'*opdfit2*/
static inline void
fit_do_hxpt(curecon_t *curecon, const curcell *xout, float alpha, stream_t &stream){
    const int nfit=curecon->opdfit->nx;
    const int npsr=xout->nx;
    EVENT2_INIT(4);
    EVENT2_TIC(0);
    if(curecon->xcache){
	curecon->xcache->m->zero(stream);
	gpu_prop_grid_do<<<dim3(3,3,npsr),dim3(16,16),0,stream>>>
	    (curecon->hxp1data, curecon->opdfit2->pm, curecon->xcache->pm, nfit, npsr, alpha, curecon->fitwt->p, 't');
	EVENT2_TIC(2);
	gpu_prop_grid_do<<<dim3(3,3,npsr),dim3(16,16),0,stream>>>
	    (curecon->hxp0data, curecon->xcache->pm, xout->pm, 1, npsr, alpha, NULL, 't');
	EVENT2_TIC(3);
    }else{
	gpu_prop_grid_do<<<dim3(3,3,npsr), dim3(16,16),0,stream>>>
	    (curecon->hxpdata, curecon->opdfit2->pm, xout->pm, nfit, npsr, alpha, curecon->fitwt->p, 't');
	EVENT2_TIC(3);
    }
    EVENT2_TOC;
    EVENT2_PRINT("do_hxpt: zero: %.3f prop: %.3f cache: %.3f tot: %.3f\n", times[1], times[2], times[3], times[0]);
}
/**
   res_vec[i]+=sum_j(as[i][j]*b[j]);
 */
__global__ void 
inn_multi_do(float *res_vec, const float *as, const float *b, const int n){
    extern __shared__ float sb[];
    sb[threadIdx.x]=0;
    const float *a=as+n*blockIdx.y;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	sb[threadIdx.x]+=a[i]*b[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    sb[threadIdx.x]+=sb[threadIdx.x+step];
	}
    }
    if(threadIdx.x==0){
	atomicAdd(&res_vec[blockIdx.y], sb[0]);
    }
}
/**
   a[i][j]=b[j]*beta1[i]*beta2;
 */
__global__ void 
assign_multi_do(float *as, const float *b, float *beta1, float beta2, int n){
    float *a=as+blockIdx.y*n;
    float beta=beta1[blockIdx.y]*beta2;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=b[i]*beta;
    }
}
/**
   opdfit2 = (W0-W1*W1')*opdfit */
static inline void
fit_do_w(curecon_t *curecon, stream_t &stream){
    EVENT2_INIT(6);
    EVENT2_TIC(0);
    EVENT2_TIC(1);
    curecon->pis->zero(stream);
    const int nfit=curecon->opdfit->nx;
    const int nxf=curecon->opdfit->p[0]->nx;
    const int nf=curecon->opdfit->p[0]->nx*curecon->opdfit->p[0]->ny;
    //inn_multi_do takes only 28 us while curmv takes 121 us.
    inn_multi_do<<<dim3(32, nfit), dim3(DIM_REDUCE), DIM_REDUCE*sizeof(float), stream>>>
	(curecon->pis->p, curecon->opdfit->m->p, curecon->W01->W1->p, nf);
    EVENT2_TIC(2);
    //assign_multi_do takes 8 us while curmm takes 28.6 us
    assign_multi_do<<<dim3(32, nfit), dim3(256), 0, stream>>>
	(curecon->opdfit2->m->p, curecon->W01->W1->p, curecon->pis->p, -1, nf);
    EVENT2_TIC(3);
    if(curecon->W01->nW0f){ //19us
	apply_W_do<<<dim3(16, nfit), dim3(256,1), 0, stream>>> 		
	    (curecon->opdfit2->m->p, curecon->opdfit->m->p, curecon->W01->W0f, 
	     curecon->W01->W0v, nxf, curecon->W01->nW0f); 
    }
    EVENT2_TIC(4);
    if(curecon->W01->W0p){ //53 us *********Can we use B/W map to avoid this matrix?***********
	cuspmul(curecon->opdfit2->m->p, curecon->W01->W0p, 
		curecon->opdfit->m->p, nfit, 'n', 1.f, stream); //80us
    }
    EVENT2_TIC(5);
    EVENT2_TOC;
    EVENT2_PRINT("do_w: zero: %.3f W1_dot: %.3f W1_add: %.3f W0f: %.3f W0p: %.3f tot: %.3f\n",
		times[1], times[2], times[3], times[4], times[5], times[0]);
}

/**
   opdfit+=Ha*xin;
*/
static inline void
fit_do_ha(curecon_t *curecon, const curcell *xin, stream_t &stream){
    curecon->opdfit->m->zero(stream); 
    const int nfit=curecon->opdfit->nx;
    const int ndm=xin->nx;
    EVENT2_INIT(3);
    EVENT2_TIC(0);
    if(curecon->dmcache){
	/*xout->dmcache*/ 
	curecon->dmcache->m->zero(stream); 
	gpu_prop_grid_do<<<dim3(8,8,ndm),dim3(16,16),0,stream>>> 
	    (curecon->ha0data, curecon->dmcache->pm, xin->pm, 
	     1, ndm, 1.f,NULL,'n');//25 us
	EVENT2_TIC(1);
	/*dmcache->opdfit*/ 
	gpu_prop_grid_do<<<dim3(8,8,nfit),dim3(16,16),0,stream>>> 
	    (curecon->ha1data, curecon->opdfit->pm, curecon->dmcache->pm, 
	     nfit, ndm, 1.f,NULL, 'n'); // 87 us
	EVENT2_TIC(2);
    }else{ 
	EVENT2_TIC(1);
	/*xout->opfit*/ 
	gpu_prop_grid_do<<<dim3(8,8,nfit),dim3(16,16),0,stream>>> 
	    (curecon->hadata, curecon->opdfit->pm, xin->pm, 
	     nfit, ndm, 1.f,NULL, 'n'); 
	EVENT2_TIC(2);
    }
    EVENT2_TOC;
    EVENT2_PRINT("do_ha: cache: %.3f prop: %.3f tot: %.3f\n", times[1], times[2], times[0]);
}

/**
   xout+=alpha*HA'*opdfit2*/
static inline void 
fit_do_hat(curecon_t *curecon, curcell *xout,  float alpha, stream_t &stream){
    const int nfit=curecon->opdfit->nx;
    const int ndm=xout->nx;
    EVENT2_INIT(3);
    EVENT2_TIC(0);
    if(curecon->dmcache){ 
	/*opdfit2->dmcache*/ 
	curecon->dmcache->m->zero(stream); 
	gpu_prop_grid_do<<<dim3(8,8,ndm),dim3(16,16),0,stream>>> 
	    (curecon->ha1data, curecon->opdfit2->pm, curecon->dmcache->pm, 
	     nfit, ndm, 1.,curecon->fitwt->p, 't'); //105 us
	EVENT2_TIC(1);
	/*dmcache->xout*/ 
	gpu_prop_grid_do<<<dim3(8,8,ndm),dim3(16,16),0,stream>>> 
	    (curecon->ha0data, curecon->dmcache->pm, xout->pm, 
	     1, ndm, alpha, NULL,'t'); //278 us. ******NEED Improvement****
	EVENT2_TIC(2);
    }else{ 
	/*opfit2->xout	*/ 
	gpu_prop_grid_do<<<dim3(8,8,ndm),dim3(16,16),0,stream>>> 
	    (curecon->hadata, curecon->opdfit2->pm, xout->pm, 
	     nfit, ndm, alpha, curecon->fitwt->p, 't'); 
	EVENT2_TIC(1);
	EVENT2_TIC(2);
    } 
    EVENT2_TOC;
    EVENT2_PRINT("do_hat: prop: %.3f cache: %.3f tot: %.3f\n", times[1], times[2], times[0]);
}

/*
  Right hand size operator. 
*/
void gpu_FitR(curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream){
    curecon_t *curecon=cudata->recon;
    const RECON_T *recon=(const RECON_T *)A;
    if(!*xout){
	*xout=curcellnew(recon->ndm, 1, recon->anx, recon->any);
    }else{
	curscale((*xout)->m, beta, stream);
    }
    fit_do_hxp(curecon, xin, stream);//153 us
    fit_do_w(curecon, stream);//123 us
    fit_do_hat(curecon, *xout, alpha, stream);//390 us
}
void gpu_FitRt(curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream){
    curecon_t *curecon=cudata->recon;
    const RECON_T *recon=(const RECON_T *)A;
    if(!*xout){
	*xout=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
    }else{
	curscale((*xout)->m, beta, stream);
    }
    fit_do_ha(curecon, xin, stream);
    fit_do_w(curecon, stream);
    fit_do_hxpt(curecon, *xout, alpha, stream);
}
void gpu_FitL(curcell **xout, float beta, const void *A, const curcell *xin, float alpha, stream_t &stream){
    curecon_t *curecon=cudata->recon;
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=recon->parms;
    const int ndm=parms->ndm;
    EVENT_INIT(6);
    EVENT_TIC(0);
    if(!*xout){
	*xout=curcellnew(ndm, 1, recon->anx, recon->any);
    }else{
	curscale((*xout)->m, beta, stream);
    }   
    fit_do_ha(curecon, xin, stream);//112 us
    EVENT_TIC(1);
    fit_do_w(curecon, stream);//108 us
    EVENT_TIC(2);
    fit_do_hat(curecon, *xout, alpha, stream);//390 us
    EVENT_TIC(3);
    if(curecon->fitNW){
	curmv(curecon->dotNW->p, 0, curecon->fitNW, xin->m->p, 't', 1, stream);
	curmv((*xout)->m->p, 1, curecon->fitNW, curecon->dotNW->p, 'n', alpha, stream);
    }
    EVENT_TIC(4);
    if(curecon->actslave){
	cuspmul((*xout)->m->p, curecon->actslave, xin->m->p, 1,'n', alpha, stream);
    }
    EVENT_TIC(5);
    EVENT_TOC;
    EVENT_PRINT("FitL HA: %.3f W: %.3f HAT: %.3f NW: %.3f SL: %.3f tot: %.3f\n",
	  times[1],times[2],times[3],times[4],times[5],times[0]);
}
/**
   Wrap of the DM fitting operation

   opdr is the OPD input.
   fitx is the right hand side vector computed from opdr for output. Allow NULL.
   fitr is the DM fitting result.
*/
double gpu_fit_do(const PARMS_T *parms,const RECON_T *recon, curcell *fitr, curcell *opdr, stream_t &stream){
    G_CGFUN cg_fun;
    void *cg_data;
    double res=0;
    curecon_t *curecon=cudata->recon;
    curcell *rhs=curecon->fitrhs;
    EVENT_INIT(3);
    EVENT_TIC(0);
    curmat *tmp=NULL;
    if(parms->gpu.fit==1){//sparse matrix
	cumuv(&curecon->fitrhs, 0, &curecon->FR, opdr, 1, stream);
	cg_fun=(G_CGFUN) cumuv;
	cg_data=&curecon->FL;
    }else{
	gpu_FitR(&curecon->fitrhs, 0, recon, opdr, 1, stream);
	cg_fun=(G_CGFUN) gpu_FitL;
	cg_data=(void*)recon;
    }
    EVENT_TIC(1);
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
	if((res=gpu_pcg(fitr, (G_CGFUN)cg_fun, cg_data, NULL, NULL, rhs, &curecon->cgtmp_fit,
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
    EVENT_TIC(2);
    EVENT_TOC;
    EVENT_PRINT("Fit RHS: %6.0f, LHS: %6.0f, Tot: %5.0f\n", times[1], times[2], times[0]);
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
    dcell *rb=NULL;
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
    muv_trans(&rb, &recon->FR, rhsc, 1);
    dcellwrite(rb, "CPU_FitRt");
    muv_trans(&rb, &recon->FR, rhsc, -1);
    dcellwrite(rb, "CPU_FitRt2");
    muv(&lc, &recon->FL, rhsc, 1);
    dcellwrite(lc, "CPU_FitL");
    muv(&lc, &recon->FL, rhsc, -1);
    dcellwrite(lc, "CPU_FitL2");
    dcellzero(lc);
    muv(&lc, &recon->FL, rhsc, 1);
    dcellwrite(lc, "CPU_FitL3");
    dcellzero(lc);
    for(int i=0; i<5; i++){
	muv_solve(&lc, &recon->FL, NULL, rhsc);
	dcellwrite(lc, "CPU_FitCG%d", i);
    }
   
    curcell *rhsg=NULL;
    curcell *lg=NULL;
    curcell *rg=NULL;
    if(parms->gpu.fit==2){
	gpu_FitR(&rhsg, 0, recon, curecon->opdr, 1, stream);
	curcellwrite(rhsg, "GPU_FitR");
	gpu_FitRt(&rg, 0, recon, rhsg, 1, stream);
	curcellwrite(rg, "GPU_FitRt");
	gpu_FitRt(&rg, 1, recon, rhsg, -1, stream);
	curcellwrite(rg, "GPU_FitRt2");
	gpu_FitL(&lg, 0, recon, rhsg, 1, stream);
	curcellwrite(lg, "GPU_FitL");
	gpu_FitL(&lg, 1, recon, rhsg, -1, stream);
	curcellwrite(lg, "GPU_FitL2");
	gpu_FitL(&lg, 0, recon, rhsg, 1, stream);
	curcellwrite(lg, "GPU_FitL3");
  
	curcellzero(lg, curecon->cgstream);
	for(int i=0; i<5; i++){
	    gpu_pcg(lg, (G_CGFUN)gpu_FitL, (void*)recon, NULL, NULL, rhsg, &curecon->cgtmp_fit,
		    simu->parms->recon.warm_restart, parms->fit.maxit, curecon->cgstream);
	    curcellwrite(lg, "GPU_FitCG%d",i);
	}
    }else{
	cumuv(&rhsg, 0, &curecon->FR, curecon->opdr, 1, stream);
	curcellwrite(rhsg, "GPU_FitR");
	cumuv_trans(&rg, 0, &curecon->FR, rhsg, 1, stream);
	curcellwrite(rg, "GPU_FitRt");
	cumuv_trans(&rg, 1, &curecon->FR, rhsg, -1, stream);
	curcellwrite(rg, "GPU_FitRt2");
	cumuv(&lg, 0, &curecon->FL, rhsg, 1, stream);
	curcellwrite(lg, "GPU_FitL");
	cumuv(&lg, 1, &curecon->FL, rhsg, -1, stream);
	curcellwrite(lg, "GPU_FitL2");
	cumuv(&lg, 0, &curecon->FL, rhsg, 1, stream);
	curcellwrite(lg, "GPU_FitL3");
  
	curcellzero(lg, curecon->cgstream);
	for(int i=0; i<5; i++){
	    gpu_pcg(lg, (G_CGFUN)cumuv, &curecon->FL, NULL, NULL, rhsg, &curecon->cgtmp_fit,
		    simu->parms->recon.warm_restart, parms->fit.maxit, curecon->cgstream);
	    curcellwrite(lg, "GPU_FitCG%d",i);
	}
    }
    CUDA_SYNC_DEVICE;
    exit(0);
}
