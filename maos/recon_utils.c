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

#include "recon_utils.h"
#include "setup_recon.h"
#include "ahst.h"
#include "sim.h"
#undef TIMING
#define TIMING 0
#if !TIMING
#define TIC_tm
#define tic_tm
#define toc_tm(A)
#else
#define TIC_tm TIC
#define tic_tm tic
#define toc_tm(A) toc2(A);tic
#endif
/**
   \file recon_utils.c
   Reusable utilities for wavefront reconstruction and DM fitting.
*/

/**
   Apply Laplacian2 to xin and accumulates to xout.
*/
void apply_L2(dcell **xout, const dspcell *L2, const dcell *xin, 
	      double alpha){
    dcell *xx=NULL;
    dcellmm(&xx, L2, xin, "nn", 1.);
    dcellmm(xout, L2, xx, "tn", alpha);
    dcellfree(xx);
}
/**
   Apply turbulence invpsd to xin in Fourier space, scaled by alpha and add to xout.
   do nothing if xb != yb, since we apply to diagonal only.
   if xb==-1, do all blocks. 
   \todo: this is not thread safe.
*/
void apply_invpsd(dcell **xout, const void *A, const dcell *xin, double alpha, int xb, int yb){
    if(xb!=yb) return;
    const INVPSD_T *extra=(const INVPSD_T*)A;
    dcell *invpsd=extra->invpsd;
    ccell *fftxopd=extra->fftxopd;
    int ips1, ips2;
    if(xb<0){/*do all cells */
	ips1=0; 
	ips2=xin->nx*xin->ny;
    }else{/*do the specified cell */
	ips1=xb;
	ips2=xb+1;
    }
    for(int ips=ips1; ips<ips2; ips++){
	/*for(int ips=0; ips<xin->nx*xin->ny; ips++){ */
	long nx=fftxopd->p[ips]->nx;
	long ny=fftxopd->p[ips]->ny;
	if(extra->square){
	    dmat *xini=dref_reshape(xin->p[ips], nx, ny);
	    ccpd(&fftxopd->p[ips], xini);
	    dfree(xini);
	}else{
	    czero(fftxopd->p[ips]);
	    cembed_locstat(&fftxopd->p[ips], 0, extra->xloc->p[ips], xin->p[ips]->p, alpha, 0);
	}
	cfft2(fftxopd->p[ips],-1);
	ccwmd(fftxopd->p[ips], invpsd->p[ips], 1);
	cfft2(fftxopd->p[ips],1);
	if(extra->square){
	    dmat *xouti=NULL;
	    xouti=dref_reshape((*xout)->p[ips], nx, ny);
	    creal2d(&xouti,1,fftxopd->p[ips],alpha);
	    dfree(xouti);
	}else{
	    cembed_locstat(&fftxopd->p[ips], 1, extra->xloc->p[ips], (*xout)->p[ips]->p, 1, 1);
	}
    }
}

/**
   Apply fractal regularization to x, scaled by alpha.
   do nothing if xb != yb, since we apply to diagonal only.
   if xb==-1, do all blocks. 
   \todo: this is not thread safe.
*/
void apply_fractal(dcell **xout, const void *A, const dcell *xin, double alpha, int xb, int yb){
    if(xb!=yb) return;
    const FRACTAL_T *extra=(const FRACTAL_T*)A;
    int ips1, ips2;
    if(xb<0){/*do all cells */
	ips1=0; 
	ips2=xin->nx*xin->ny;
    }else{/*do the specified cell */
	ips1=xb;
	ips2=xb+1;
    }
    for(int ips=ips1; ips<ips2; ips++){
	/*for(int ips=0; ips<xin->nx*xin->ny; ips++){ */
	dzero(extra->xopd->p[ips]);
	double r0i=extra->r0*pow(extra->wt[ips], -3./5.);
	dembed_locstat(&extra->xopd->p[ips], 0, extra->xloc->p[ips], xin->p[ips]->p, 
		       alpha*extra->scale, 0);
	fractal_inv(extra->xopd->p[ips],
		    extra->xloc->p[ips]->dx, r0i, extra->L0, extra->ninit);
	fractal_inv_trans(extra->xopd->p[ips],
			  extra->xloc->p[ips]->dx, r0i, extra->L0, extra->ninit);
	dembed_locstat(&extra->xopd->p[ips], 1, extra->xloc->p[ips], (*xout)->p[ips]->p, 1, 1);
    }
}

/**
   Removing Tip/Tilt/Focus from LGS grads. TTF is the Tip/tilt/focus modes, and
   PTTF is the pseudo inverse of it, weighted by subaperture noise.  It is not
   diagonal if differential focus is removed. So cannot be in Tomo_nea*/
void TTFR(dcell* x, const dcell *TTF, const dcell *PTTF){
    if(!TTF || !PTTF){
	return;
    }
    dcell *junk=NULL;
    dcellmm(&junk, PTTF, x, "nn", 1);
    dcellmm(&x, TTF, junk, "nn", -1);
    dcellfree(junk);
}
/**
   Removing Tip/Tilt/Focus from LGS grads. TTF is the Tip/tilt/focus modes, and
   PTTF is the pseudo inverse of it, weighted by subaperture noise.  It is not
   diagonal if differential focus is removed. So cannot be in Tomo_nea*/
void TTFRt(dcell* x, const dcell *TTF, const dcell *PTTF){
    if(!TTF || !PTTF){
	return;
    }
    dcell *junk=NULL;
    dcellmm(&junk, TTF, x, "tn", 1);
    dcellmm(&x, PTTF, junk, "tn", -1);
    dcellfree(junk);
}

/**
   Apply weighting W0/W1 to a vector. W0*x-W1*(W1'*x)
*/
static void applyWeach(dmat *xin, const dsp *W0, const dmat *W1, const double wt){
    if(!W0 || !W1) {
	warning("W0 or W1 is NULL\n");
	return;
    }
    dmat *xout=NULL;
    dmat *tmp=NULL;
    dspmm(&xout, W0, xin, "nn", wt);
    dmm(&tmp, 0, W1, xin, "tn", -1);
    dmm(&xout, 1, W1, tmp, "nn", wt);
    dcp(&xin, xout);
    dfree(xout); dfree(tmp);
}
/**
   apply weighting W0/W1 with weighting wt for each block.
   W0*xin-W1*(W1'*xin);
*/
void applyW(dcell *xin, const dsp *W0, const dmat *W1, const double *wt){
    const int nevl=xin->nx;
    for(int iy=0; iy<xin->ny; iy++){
	for(int ievl=0; ievl<nevl; ievl++){
	    int ind=iy*nevl+ievl;
	    applyWeach(xin->p[ind], W0, W1, wt[ievl]);
	}
    }
}

/**
   Compute W0/W1 weighted dot product: \f$A^T(W0 B-W1 (W1^T B))\f$
*/
dcell* calcWmcc(const dcell *A, const dcell *B, const dsp *W0, 
		const dmat *W1, const dmat *wt){

    assert(wt->nx==B->nx && wt->ny==1 && A->nx == B->nx);
    const int nevl=B->nx;
    dcell *res=dcellnew(A->ny, B->ny);
    for(int iy=0; iy<B->ny; iy++){
	for(int ievl=0; ievl<nevl; ievl++){
	    int ind=iy*nevl+ievl;
	    dmat *xout=NULL;
	    dmat *tmp=NULL;
	    dspmm(&xout, W0, B->p[ind], "nn", wt->p[ievl]);
	    dmm(&tmp, 0, W1, B->p[ind], "tn", -1);
	    dmm(&xout, 1, W1, tmp, "nn", wt->p[ievl]);
	    for(int ix=0; ix<A->ny; ix++){
		dmm(&res->p[ix+iy*res->nx], 1, A->p[ix*nevl+ievl], xout, "tn", 1);
	    }
	    dfree(xout);
	    dfree(tmp);
	}
    }
    return res;
}

typedef struct Tomo_T{
    const RECON_T *recon;
    const double alpha;
    const dcell *xin;/*input */
    dcell *gg;/*intermediate gradient */
    dcell *xout;/*output */
    int isR;/*Mark that we are right hand side.*/
}Tomo_T;

/**
   Speed up TomoL by gathering the first part of operations (GP*HXW) belonging
   to each WFS to facilitate threading. 

   gg  = GP * HXW * xin
*/
static void Tomo_prop_do(thread_t *info){
    Tomo_T *data=(Tomo_T*)info->data;
    const RECON_T *recon=data->recon;
    const PARMS_T *parms=global->parms;
    SIM_T *simu=global->simu;
    const int nps=recon->npsr;
    map_t xmap;/*make a temporary xmap for thread safety.*/
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
	int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	const double delay=parms->sim.dt*(parms->powfs[ipowfs].dtrat+1+parms->sim.alhi);
	dmat *xx=dnew(recon->ploc->nloc, 1);
	const double hs=parms->wfs[iwfs].hs;
	const double hc=parms->powfs[ipowfs].hc;
	for(int ips=0; ips<nps; ips++){
	    if(parms->tomo.square && !parms->dbg.tomo_hxw){
		/*Do the ray tracing instead of using HXW. */
		double ht=recon->ht->p[ips]-hc;
		double displace[2];
		displace[0]=parms->wfsr[iwfs].thetax*ht+parms->wfsr[iwfs].misreg_x;
		displace[1]=parms->wfsr[iwfs].thetay*ht+parms->wfsr[iwfs].misreg_y;
		if(parms->tomo.predict){
		    int ips0=parms->atmr.indps->p[ips];
		    displace[0]+=simu->atm->p[ips0]->vx*delay;
		    displace[1]+=simu->atm->p[ips0]->vy*delay;
		}
		double scale=1. - ht/hs;
		if(scale<0) continue;
		memcpy(&xmap, recon->xmap->p[ips], sizeof(map_t));
		xmap.p=data->xin->p[ips]->p;
		prop_grid_stat(&xmap, recon->ploc->stat, xx->p, 1, 
			       displace[0],displace[1], scale, 0, 0, 0);
	    }else{
		dspcell* HXW=recon->HXWtomo/*PDSPCELL*/;
		dspmm(&xx, P(HXW,iwfs,ips), data->xin->p[ips], "nn", 1);
	    }
	}
	/*Apply the gradient operation */
	dspmm(&data->gg->p[iwfs], recon->GP->p[iwfs], xx, "nn", 1);
	dfree(xx);
	/* For each wfs, Ray tracing takes 1.5 ms.  GP takes 0.7 ms. */
    }
}
/**
   Wrapper of Tomo_prop_do
*/
void Tomo_prop(Tomo_T *data, int nthread){
    thread_t info[nthread];
    thread_prep(info, 0, data->gg->nx, nthread, Tomo_prop_do, data);
    CALL_THREAD(info, 1);
}
/**
   Speed up TomoL by gathering the second part of operations (GP') belonging to
   each WFS to facilitate threading. 
   
   gg = GP' * NEAI * gg;
*/
static void Tomo_nea_gpt_do(thread_t *info){
    Tomo_T *data=(Tomo_T*)info->data;
    const RECON_T *recon=data->recon;
    dspcell*  NEAI=recon->saneai/*PDSPCELL*/;
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
	dmat *gg2=NULL;
	/*Apply the gradient operation */
	dspmm(&gg2, P(NEAI,iwfs,iwfs), data->gg->p[iwfs], "nn", 1);
	dfree(data->gg->p[iwfs]); /*We reuse gg. */
	dspmm(&data->gg->p[iwfs], recon->GP->p[iwfs], gg2, "tn", data->alpha);
	dfree(gg2);
    }
}

static void Tomo_nea_do(thread_t *info){
    Tomo_T *data=(Tomo_T*)info->data;
    const RECON_T *recon=data->recon;
    dspcell*  NEAI=recon->saneai/*PDSPCELL*/;
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
	dmat *gg2=NULL;
	/*Apply the gradient operation */
	dspmm(&gg2, P(NEAI,iwfs,iwfs), data->gg->p[iwfs], "nn", 1);
	dcp(&data->gg->p[iwfs], gg2);
	dfree(gg2);
    }
}

/*Wrapp of Tomo_nea_do for multi-threads*/
void Tomo_nea(Tomo_T *data, int nthread, int gpt){
    thread_t info[nthread];
    thread_prep(info, 0, data->gg->nx, nthread, gpt?Tomo_nea_gpt_do:Tomo_nea_do, data);
    CALL_THREAD(info, 1);
}
/**
   Speed up TomoL by gathering the third part of operations (GP') belonging to
   each WFS to facilitate threading. gg->xout.

   xout = Cxx^-1 * xin + HXW * gg;
*/
static void Tomo_iprop_do(thread_t *info){
    Tomo_T *data=(Tomo_T*)info->data;
    const RECON_T *recon=data->recon;
    const PARMS_T *parms=global->parms;
    SIM_T *simu=global->simu;
    const int nps=recon->npsr;
    map_t xmap;
    for(int ips=info->start; ips<info->end; ips++){
	if(parms->tomo.square && !parms->dbg.tomo_hxw){
	    /*Do the ray tracing instead of using HXW. */
	    if(!data->xout->p[ips]){
		data->xout->p[ips]=dnew(recon->xloc->p[ips]->nloc, 1);
	    }
	    memcpy(&xmap, recon->xmap->p[ips], sizeof(map_t));
	    xmap.p=data->xout->p[ips]->p;
	    double ht=recon->ht->p[ips];
	    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		if(!data->gg->p[iwfs]) continue;
		const double hs=parms->wfs[iwfs].hs;
		const int ipowfs=parms->wfs[iwfs].powfs;
		const double hc=parms->powfs[ipowfs].hc;
		double displace[2];
		displace[0]=parms->wfsr[iwfs].thetax*(ht-hc)+parms->wfsr[iwfs].misreg_x;
		displace[1]=parms->wfsr[iwfs].thetay*(ht-hc)+parms->wfsr[iwfs].misreg_y;
		if(parms->tomo.predict){
		    const double delay=parms->sim.dt*(parms->powfs[ipowfs].dtrat+1+parms->sim.alhi);
		    int ips0=parms->atmr.indps->p[ips];
		    displace[0]+=simu->atm->p[ips0]->vx*delay;
		    displace[1]+=simu->atm->p[ips0]->vy*delay;
		}
		double scale=1. - ht/hs;
		if(scale<0) continue;
		prop_grid_stat_transpose(&xmap, recon->ploc->stat, data->gg->p[iwfs]->p, 1, 
					 displace[0],displace[1], scale, 0, 0, 0);
	    }
	}else{
	    dspcell* HXW=recon->HXWtomo/*PDSPCELL*/;
	    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		dspmm(&data->xout->p[ips], P(HXW,iwfs,ips), data->gg->p[iwfs], "tn", 1);
	    }
	}
	if(data->xin){/*data->xin is empty when called from TomoR */
	    switch(recon->cxxalg){
	    case 0:{/*L2 */
		dmat *xx=NULL;
		dspmm(&xx, recon->L2->p[ips+ips*nps], data->xin->p[ips],"nn", 1);
		dspmm(&data->xout->p[ips], recon->L2->p[ips+ips*nps], xx, "tn", data->alpha);
		dfree(xx);
	    }
		break;
	    case 1:
		apply_invpsd(&data->xout, recon->invpsd, data->xin, data->alpha, ips, ips);
		break;
	    case 2:
		apply_fractal(&data->xout, recon->fractal, data->xin, data->alpha, ips, ips);
		break;
	    }
	    if(recon->ZZT){
		dspmm(&data->xout->p[ips], recon->ZZT->p[ips+ips*nps], data->xin->p[ips], "tn", data->alpha);
	    }
	}
    }
}
/**
   Wrapper of Tomo_iprop_do
 */
void Tomo_iprop(Tomo_T *data, int nthread){
    thread_t info[nthread];
    thread_prep(info, 0, data->xout->nx, nthread, Tomo_iprop_do, data);
    CALL_THREAD(info, 1);
}
/**
   Apply tomography right hand operator without using assembled matrix. Fast and
   saves memory.  The operation is the same as the Tomo_nea and Tomo_iprop in
   TomoL, so merge the implemenations.

   xout=HXW'*GP'*NEAI*(1-TTF*PTTF)*gin.
*/

void TomoR(dcell **xout, const void *A, 
	   const dcell *gin, const double alpha){
    TIC_tm;tic_tm;
    const RECON_T *recon=(const RECON_T *)A;
    dcell *gg=NULL;
    dcellcp(&gg, gin);/*copy to gg so we don't touch the input. */
    TTFR(gg, recon->TTF, recon->PTTF);
    if(!*xout){
	*xout=dcellnew(recon->npsr, 1);
    }
    Tomo_T data={recon, alpha, NULL, gg, *xout, 1};
    Tomo_nea(&data, recon->nthread, 1);
    Tomo_iprop(&data, recon->nthread);
    dcellfree(gg);
    toc_tm("TomoR");
}

/**
   Transpose operation of TomoRt. From xout -> gin
 */
void TomoRt(dcell **gout, const void *A, 
	    const dcell *xin, const double alpha){
    const RECON_T *recon=(const RECON_T *)A;
    if(!*gout){
	*gout=dcellnew(recon->saneai->nx, 1);
    }
    Tomo_T data={recon, alpha, xin, *gout, NULL, 1};
    Tomo_prop(&data, recon->nthread);
    /*Using Tomo_nea followed by TTFRt is equilvaent as using TTFR followed by Tomo_nea.*/
    TTFR(*gout, recon->TTF, recon->PTTF);
    Tomo_nea(&data, recon->nthread, 0);
}
/**
   Apply tomography left hand side operator without using assembled matrix. Fast
   and saves memory. Only useful in CG. Accumulates to xout;

   Some timing information: 
   Generated in T410s with Intel Core i5 520 M @ 2.4 GHz (1 T is 1 thread, 2 T is 2 thread)
   and Poweredge with dual Xeon X5355 2.66Ghz (1 P is 1 thread, 8 P is 8 thread)
   Time  1 T  2 T  1 P  8 P
   HXW:  7.8  5.4  7.7  5.0
   GP:   4.3  2.6  3.6  1.4
   TTFR: 0.3  0.3  0.4  0.6
   neai: 0.5  0.3  0.5  0.4
   GP':  3.2  2.1  3.4  1.0
   HXW': 5.5  4.7  7.4  4.9
   cxx:  3.2  2.5  3.5  2.7
*/

#define test_TomoL 0
void TomoL(dcell **xout, const void *A, 
	   const dcell *xin, const double alpha){
    TIC_tm;tic_tm;
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=global->parms;
    assert(xin->ny==1);/*modify the code for ny>1 case. */
    dcell *gg=dcellnew(parms->nwfsr, 1);
    if(!*xout){
	*xout=dcellnew(recon->npsr, 1);
    }
    Tomo_T data={recon, alpha, xin, gg, *xout, 0};
  
    Tomo_prop(&data, recon->nthread);  

    if(!parms->recon.split || parms->tomo.splitlrt){
	/*Remove global Tip/Tilt, differential focus only in integrated
	  tomography to limit noise propagation (?).*/
	TTFR(gg, recon->TTF, recon->PTTF);
    }
    Tomo_nea(&data, recon->nthread, 1);
    Tomo_iprop(&data, recon->nthread);
    /* 
       square=1  square=0 (1 thread on T410s)
       prop:  takes 6 ms 13 ms
       nea:   takes 4 ms 4 ms
       iprop: takes 6 ms 9 ms
    */
    dcellfree(gg);
    /*Tikhonov regularization is not added because it is not necessary in CG mode.*/
    toc_tm("TomoL");
}

/**
   Apply fit right hand side matrix in CG mode without using assembled matrix.
   Slow. don't use. Assembled matrix is faster because of multiple directions.
*/
void FitR(dcell **xout, const void *A, 
	  const dcell *xin, const double alpha){
    const FIT_T *fit=(const FIT_T *)A;
    const int nfit=fit->thetax->nx;
    dcell *xp=dcellnew(nfit,1);

    if(!xin){/*xin is empty. We will trace rays from atmosphere directly */
	const PARMS_T *parms=global->parms;
	SIM_T *simu=global->simu;
	int isim=fit->notrecon?simu->wfsisim:simu->reconisim;
	const double atmscale=simu->atmscale?simu->atmscale->p[isim]:1;
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=fit->hs->p[ifit];
	    xp->p[ifit]=dnew(fit->floc->nloc,1);
	    for(int ips=0; ips<parms->atm.nps; ips++){
		const double ht = parms->atm.ht->p[ips]-fit->floc->ht;
		double scale=1-ht/hs;
		double displace[2];
		if(scale<0) continue;
		displace[0]=fit->thetax->p[ifit]*ht-simu->atm->p[ips]->vx*isim*parms->sim.dt;
		displace[1]=fit->thetay->p[ifit]*ht-simu->atm->p[ips]->vy*isim*parms->sim.dt;
		prop_grid(simu->atm->p[ips], fit->floc, xp->p[ifit]->p, 
			  atmscale, displace[0], displace[1], scale, 1, 0, 0);
	    }
	}
    }else if(fit->HXF){
	dcellmm(&xp, fit->HXF, xin, "nn", 1.);
    }else{/*Do the ray tracing from xloc to ploc */
	const int npsr=fit->xloc->nx;
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=fit->hs->p[ifit];
	    xp->p[ifit]=dnew(fit->floc->nloc,1);
	    for(int ips=0; ips<npsr; ips++){
		const double ht = fit->xloc->p[ips]->ht-fit->floc->ht;
		double scale=1-ht/hs;
		double displace[2];
		if(scale<0) continue;
		displace[0]=fit->thetax->p[ifit]*ht;
		displace[1]=fit->thetay->p[ifit]*ht;
		prop_nongrid(fit->xloc->p[ips], xin->p[ips]->p, fit->floc, 
			     xp->p[ifit]->p, 1, displace[0], displace[1], scale, 0, 0);
	    }
	}
    }
    //writebin(xp, "CPU_FitR_x1");
    applyW(xp, fit->W0, fit->W1, fit->wt->p);
    //writebin(xp, "CPU_FitR_x2");
    dcellmm(xout, fit->HA, xp, "tn", alpha);
    //writebin(*xout, "CPU_FitR_x3");
    dcellfree(xp);
}
/**
   Apply fit left hand side matrix in CG mode without using assembled
   matrix. Slow. don't use. Assembled matrix is faster because of multiple
   directions.  */
void FitL(dcell **xout, const void *A, 
	  const dcell *xin, const double alpha){
    const FIT_T *fit=(const FIT_T *)A;
    dcell *xp=NULL;
    dcellmm(&xp, fit->HA, xin, "nn", 1.);
    applyW(xp, fit->W0, fit->W1, fit->wt->p);
    dcellmm(xout, fit->HA, xp, "tn", alpha);
    dcellfree(xp);xp=NULL;
    dcellmm(&xp,fit->NW, xin, "tn", 1);
    dcellmm(xout,fit->NW, xp, "nn", alpha);
    dcellfree(xp);
    if(fit->actslave){
	dcellmm(xout, fit->actslave, xin, "nn", alpha);
    }
}


/**
   Convert block of 2x2 neas to sparse matrix.
 */
dsp *nea2sp(dmat **nea, long nsa){
    dsp *sanea=dspnew(nsa*2, nsa*2, 4*nsa);
    spint *pp=sanea->p;
    spint *pi=sanea->i;
    double *px=sanea->x;
    long count=0;
    for(long isa=0; isa<nsa; isa++){
	/*Cxx */
	pp[isa]=count;
	pi[count]=isa;
	px[count]=nea[isa]->p[0];
	count++;
	/*Cyx */
	pi[count]=isa+nsa;
	px[count]=nea[isa]->p[1];
	count++;
    }
    for(long isa=0; isa<nsa; isa++){
	/*Cxy */
	pp[isa+nsa]=count;
	pi[count]=isa;
	px[count]=nea[isa]->p[2];
	count++;
	/*Cyy */
	pi[count]=isa+nsa;
	px[count]=nea[isa]->p[3];
	count++;
    }
    pp[nsa*2]=count;
    return sanea;
}

/**
   Compute and save PSF reconstruction telemetry.  For pseudo open loop
   estimations, like high order loop, we add opdr to, and subtract dmpsol from
   the OPD.  For closed loop estimations, like ahst low order or lsr, we add
   dmerr_lo, and dmerr to the OPD.*/
void psfr_calc(SIM_T *simu, dcell *opdr, dcell *dmpsol, dcell *dmerr, dcell *dmerr_lo){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    /* The tomography estimates, opdr is pseudo open loop estimates. We need to
       subtract the constribution of the added DM command to form closed loop
       estimates.
    */
    dcell *dmadd=NULL;
    if(opdr && parms->dbg.useopdr){
	/* The original formulation using Hx*x-Ha*a.  Changed from ploc to plocs
	   on July 18, 2011. Ploc is too low sampled. Make sure the sampling of
	   plocs is not too big.  Deprecated. Use dm space instead.
	*/
	if(dmpsol){/*Pseudo OL estimates */
	    if(parms->recon.split==1){
		/* We will remove NGS modes from dmlast which is in NULL modes of
		   tomography reconstructor (is this 100% true)?  SHould we remove NGS
		   modes from final OPD, xx, instead?*/
		dcell *tmp = dcelldup(dmpsol);/*The DM command used for high order. */
		remove_dm_ngsmod(simu, tmp);/*remove NGS modes as we do in ahst. */
		dcelladd(&dmadd, 1, tmp, -1);
		dcellfree(tmp);
	    }else{
		dcelladd(&dmadd, 1, dmpsol, -1);
	    }
	}

	loc_t *locs=simu->aper->locs;
	dmat *xx = dnew(locs->nloc, 1);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    double hs = parms->evl.hs->p[ievl];
	    if(parms->evl.psfr->p[ievl]){
		dzero(xx);
		if(opdr){
		    map_t xmap;
		    const int npsr=recon->npsr;
		    /*First compute residual opd: Hx*x-Ha*a*/
		    for(int ips=0; ips<npsr; ips++){
			const double ht = recon->ht->p[ips];
			double scale=1-ht/hs;
			double dispx=parms->evl.thetax->p[ievl]*ht;
			double dispy=parms->evl.thetay->p[ievl]*ht;
			if(scale<0) continue;
			if(parms->tomo.square){/*square xloc */
			    memcpy(&xmap, recon->xmap->p[ips], sizeof(map_t));
			    xmap.p=opdr->p[ips]->p;
			    prop_grid_stat(&xmap, locs->stat, xx->p, 1, 
					   dispx, dispy, scale, 0, 0, 0);
			}else{
			    prop_nongrid(recon->xloc->p[ips], opdr->p[ips]->p, locs,
					 xx->p, 1, dispx, dispy, scale, 0, 0);
			}
		    }
		}
		if(dmadd){
		    for(int idm=0; idm<parms->ndm; idm++){
			const double ht = parms->dm[idm].ht;
			double scale=1.-ht/hs;
			double dispx=parms->evl.thetax->p[ievl]*ht;
			double dispy=parms->evl.thetay->p[ievl]*ht;
			if(scale<0) continue;
			prop_nongrid(recon->aloc->p[idm], dmadd->p[idm]->p, locs,
				     xx->p, 1, dispx, dispy, scale, 0, 0);
		    }
		}
		dmm(&simu->ecov->p[ievl], 1, xx, xx, "nt", 1);
		if(parms->dbg.ecovxx){
		    zfarr_push(simu->save->ecovxx[ievl], simu->reconisim, xx);
		}
	    }/*if psfr[ievl] */
	}/*ievl */
	dfree(xx);
    }else{/*Do Ha in postproc, so just do a. */
	if(dmerr){/*high order closed loop estimates. (lsr)*/
	    dcelladd(&dmadd, 1, dmerr, 1);
	}
	if(dmerr_lo){/*In AHST, dmerr_lo is CL Estimation.*/
	    addlow2dm(&dmadd, simu, dmerr_lo, 1);
	}
	dcellmm(&simu->ecov, dmadd, dmadd, "nt", 1);
    }
    dcellfree(dmadd);
}

/**
   Shift gradient when new gradients are ready (in the end of parallel section
   in sim in CL or wfsgrad in OL). Do not execute in parallel with other
   routines. In GLAO mode, also averaged gradients from the same type of powfs.
*/
void shift_grad(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(parms->sim.evlol || parms->sim.idealfit || parms->sim.idealtomo) return;
    if(PARALLEL==2){
	if(simu->wfsisim>0){
	    while(simu->wfsgrad_count<1){//not being consumed yet
		//dbg("waiting: wfsgrad_count is %d, need %d\n", simu->wfsgrad_count, 1);
		pthread_cond_wait(&simu->wfsgrad_condw, &simu->wfsgrad_mutex);
		pthread_mutex_unlock(&simu->wfsgrad_mutex);
	    }
	    //dbg("ready: wfsgrad_count is ready: %d\n", simu->wfsgrad_count);
	}
    }
    if(parms->recon.glao){
	/* Average the gradients in GLAO mode. */
	if(simu->gradlastcl){
	    dcellzero(simu->gradlastcl);
	}else{
	    long nnx[parms->nwfsr];
	    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		int ipowfs=parms->wfsr[iwfs].powfs;
		nnx[iwfs]=simu->powfs[ipowfs].saloc->nloc*2;
	    }
	    simu->gradlastcl=dcellnew3(parms->nwfsr, 1, nnx, NULL);
	}
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    const double scale=1./parms->powfs[ipowfs].nwfs;
	    for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[indwfs];
		dadd(&simu->gradlastcl->p[ipowfs], 1., simu->gradcl->p[iwfs], scale);
		if(simu->gradoff->p[iwfs]){
		    //Gradient offset due to mainly NCPA calibration. Must be after gain adjustment.
		    dadd(&simu->gradlastcl->p[ipowfs], 1, simu->gradoff->p[iwfs], -parms->dbg.gradoff_scale*scale);
		}
		if(parms->dbg.gradoff){
		    info_once("Add injected gradient offset vector\n");
		    int icol=(simu->wfsisim+1)%parms->dbg.gradoff->ny;
		    dadd(&simu->gradlastcl->p[ipowfs], 1, P(parms->dbg.gradoff, iwfs, icol), -1*scale);
		}
	    }
	}
    }else{
	dcellcp(&simu->gradlastcl, simu->gradcl); 
	//Gradient offset due to mainly NCPA calibration. Must be after gain adjustment.
	//Add gradient offset to gradlastcl, so that lpfocus on gradcl does not affect it.
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    if(simu->gradoff->p[iwfs]){
		dadd(&simu->gradlastcl->p[iwfs], 1, simu->gradoff->p[iwfs], -parms->dbg.gradoff_scale);
	    }
	    if(parms->dbg.gradoff){
		info_once("Add injected gradient offset vector\n");
		int icol=(simu->wfsisim+1)%parms->dbg.gradoff->ny;
		dadd(&simu->gradlastcl->p[iwfs], 1, P(parms->dbg.gradoff, iwfs, icol), -1);
	    }
	}
    }
    if(PARALLEL==2){
	//Signal recon wfsgrad is ready/
	simu->wfsgrad_isim=simu->wfsisim;
	simu->wfsgrad_count=0;//reset the counter
	pthread_cond_broadcast(&simu->wfsgrad_condr);
	//dbg("wfsgrad_isim is set to %d\n", simu->wfsgrad_isim);
    }
}

/**
   Parse the input dead actuator location to actuator indices based on aloc.
   2015-03-30: build a mask for dead actuators based on coordinate.
*/
lmat* loc_coord2ind(loc_t *aloc,       /**<[in] Aloc*/
		    const char *fndead /**<[in] File containing dead actuators*/
    ){
    dmat *dead=dread("%s", fndead);
    if(dead->ny!=2 && dead->ny!=3){
	error("%s must contain 2 or 3 columns of data\n", fndead);
    }
    loc_create_map(aloc);
    map_t *map=aloc->map;
    double ox=aloc->map->ox;
    double oy=aloc->map->oy;
    double dx1=1./aloc->dx;
    dmat*  ps=dead/*PDMAT*/;
    lmat *out=lnew(aloc->nloc, 1);
    for(long jact=0; jact<dead->nx; jact++){
	long mapx=(long)round((P(ps,jact,0)-ox)*dx1);
	long mapy=(long)round((P(ps,jact,1)-oy)*dx1);
	long iact=loc_map_get(map, mapx, mapy)-1;
	if(iact>=0){
	    out->p[iact]=(dead->ny==3?P(ps,jact,2)*1e9:1);//integer in nm.
	}
    }
    dfree(dead);
    return out;
}

/**
   Prepare arrays for cn2 estimation. Multiple pairs can be used to do Cn2
Estimation. The result will be an average of them.  */


cn2est_t *cn2est_prepare(const PARMS_T *parms, const POWFS_T *powfs){
    dmat *pair=parms->cn2.pair;
    int npair=pair->nx*pair->ny;
    int ipowfs=-1;
    for(int ind=0; ind<npair; ind++){
	int iwfs=(int)pair->p[ind];
	if(ipowfs==-1){
	    ipowfs=parms->wfs[iwfs].powfs;
	}else if(ipowfs!=parms->wfs[iwfs].powfs){
	    error("All wfs in parms->cn2.pair do not belong to the same powfs\n");
	}
    }
    dmat *wfstheta=dnew(parms->nwfs, 2);
    dmat*  ptheta=wfstheta/*PDMAT*/;
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	P(ptheta,iwfs,0)=parms->wfs[iwfs].thetax;
	P(ptheta,iwfs,1)=parms->wfs[iwfs].thetay;
    }
    dmat *ht=0;
    if(parms->cn2.keepht){
	ht=dnew(parms->atmr.nps, 1);
	for(int iht=0; iht<ht->nx; iht++){
	    ht->p[iht]=parms->atmr.ht->p[iht];
	}
    }else{
	int nht=parms->cn2.nhtomo;
	ht=dnew(nht,1);
	double dht=parms->cn2.hmax/(double)(nht-1);
	if(dht<=0) error("Invalid hmax or nhtomo\n");
	for(int iht=0; iht<nht; iht++){
	    ht->p[iht]=dht*iht;
	}
    }
    dmat *hs=dnew(parms->nwfs, 1);
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	hs->p[iwfs]=parms->wfs[iwfs].hs;
    }
    cn2est_t *cn2est=cn2est_new(pair, wfstheta, powfs[ipowfs].saloc, powfs[ipowfs].saa, parms->cn2.saat, 
			     hs, ht, parms->cn2.keepht, parms->atmr.L0);
    cn2est->os=dnew(ht->nx, 1);
    if(!parms->cn2.keepht){
	/*preserve the number of over sampled layers. */
	int osc=0;
	for(int ips=0; ips<parms->atmr.nps; ips++){
	    if(parms->atmr.os->p[ips]>1){
		osc++;
		if(parms->atmr.os->p[ips]!=2){
		    error("os is not 2. adept this code to it.\n");
		}	
	    }
	}
	for(int iht=0; iht<ht->nx; iht++){
	    if(iht<osc){
		cn2est->os->p[iht]=2;
	    }else{
		cn2est->os->p[iht]=1;
	    }
	}
    }else{
	for(int iht=0; iht<ht->nx; iht++){
	    cn2est->os->p[iht]=parms->atmr.os->p[iht];
	}
    }
    if(parms->cn2.verbose){
	cn2est->dmht=dnew(parms->ndm, 1);
	for(int idm=0; idm<parms->ndm; idm++){
	    cn2est->dmht->p[idm]=parms->dm[idm].ht;
	}
    }
    if(parms->save.setup){
	writebin(cn2est->overlapi, "cn2_overlapi");
	writebin(cn2est->iPnk,"cn2_iPnk");
	writebin(cn2est->Pnk,"cn2_Pnk");
	writebin(cn2est->ht,"cn2_ht");
	writebin(cn2est->wtconvert,"cn2_wtconvert");
    }
    dfree(wfstheta);
    dfree(ht);
    dfree(hs);
    return cn2est;
}
/**
   Implemented mechanism to move height of layers.
 */
static void cn2est_moveht(RECON_T *recon){
    (void)recon;
    /*cn2est_t *cn2est=recon->cn2est; */
    /*
      Implemented mechanism to move height of layers. Need to redo HXF, GX, etc.
    */
    error("moveht not yet implemented");
}

/**
   Wrapper of Cn2 Estimation operations in recon.c
*/
void cn2est_isim(RECON_T *recon, const PARMS_T *parms, dcell *grad, int *tomo_update){
    cn2est_t *cn2est=recon->cn2est;
    cn2est_push(cn2est, grad);
    static int icn2=-1;
    if(cn2est->count%parms->cn2.step == 0){
	icn2++;
	int nset=cn2est->count/parms->cn2.step;
	int reset=parms->cn2.reset && (nset%parms->cn2.reset)==0;
	cn2est_est(cn2est, parms->cn2.verbose, reset);/*do the CN2 estimation */
	if(global->simu->cn2est){
	    global->simu->cn2est->p[0]->p[icn2]=cn2est->r0m;
	    memcpy(PCOL(global->simu->cn2est->p[1], icn2),cn2est->wtrecon->p[0]->p, 
		   cn2est->htrecon->nx*sizeof(double));
	}
	if(parms->cn2.tomo){
	    if(parms->cn2.moveht){
		cn2est_moveht(recon);
	    }

	    if(parms->cn2.verbose){
		info("Updating tomography weights\n");
	    }
	    /*Changes recon parameters. cannot be parallel with tomofit(). */
	    /*wtrecon is referenced so should be updated automaticaly. */
	    if(recon->wt->p!=cn2est->wtrecon->p[0]->p){
		dfree(recon->wt);
		recon->wt=dref(cn2est->wtrecon->p[0]);
	    }
	    recon->r0=cn2est->r0m;
	    recon->L0=cn2est->L0;
	
	    *tomo_update=1;
	}
    }
}
