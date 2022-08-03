/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
 * \file recon_utils.h
   Reusable utilities for wavefront reconstruction and DM fitting, supporting recon.c
*/

/**
   Apply Laplacian2 to xin and accumulates to xout.
*/
void apply_L2(dcell** xout, const dspcell* L2, const dcell* xin,
	real alpha){
	dcell* xx=NULL;
	dspcellmm(&xx, L2, xin, "nn", 1.);
	dspcellmm(xout, L2, xx, "tn", alpha);
	dcellfree(xx);
}
/**
   Apply turbulence invpsd to xin in Fourier space, scaled by alpha and add to xout.
   do nothing if xb != yb, since we apply to diagonal only.
   if xb==-1, do all blocks.
   \todo: this is not thread safe.
*/
void apply_invpsd(dcell** xout, const void* A, const dcell* xin, real alpha, int xb, int yb){
	if(xb!=yb) return;
	const invpsd_t* extra=(const invpsd_t*)A;
	dcell* invpsd=extra->invpsd;
	ccell* fftxopd=extra->fftxopd;
	int ips1, ips2;
	if(xb<0){/*do all cells */
		ips1=0;
		ips2=NX(xin)*NY(xin);
	} else{/*do the specified cell */
		ips1=xb;
		ips2=xb+1;
	}
	for(int ips=ips1; ips<ips2; ips++){
	/*for(int ips=0; ips<NX(xin)*NY(xin); ips++){ */
		long nx=P(fftxopd, ips)->nx;
		long ny=P(fftxopd, ips)->ny;
		if(extra->square){
			dmat* xini=dref_reshape(P(xin, ips), nx, ny);
			ccpd(&P(fftxopd, ips), xini);
			dfree(xini);
		} else{
			czero(P(fftxopd, ips));
			cembed_locstat(&P(fftxopd, ips), 0, P(extra->xloc, ips), P(P(xin, ips)), alpha, 0);
		}
		cfft2(P(fftxopd, ips), -1);
		ccwmd(P(fftxopd, ips), P(invpsd, ips), 1);
		cfft2(P(fftxopd, ips), 1);
		if(extra->square){
			dmat* xouti=NULL;
			xouti=dref_reshape(P(*xout, ips), nx, ny);
			creal2d(&xouti, 1, P(fftxopd, ips), alpha);
			dfree(xouti);
		} else{
			cembed_locstat(&P(fftxopd, ips), 1, P(extra->xloc, ips), P(P(*xout, ips)), 1, 1);
		}
	}
}

/**
   Apply fractal regularization to x, scaled by alpha.
   do nothing if xb != yb, since we apply to diagonal only.
   if xb==-1, do all blocks.
   \todo: this is not thread safe.
*/
void apply_fractal(dcell** xout, const void* A, const dcell* xin, real alpha, int xb, int yb){
	if(xb!=yb) return;
	const fractal_t* extra=(const fractal_t*)A;
	int ips1, ips2;
	if(xb<0){/*do all cells */
		ips1=0;
		ips2=NX(xin)*NY(xin);
	} else{/*do the specified cell */
		ips1=xb;
		ips2=xb+1;
	}
	for(int ips=ips1; ips<ips2; ips++){
	/*for(int ips=0; ips<NX(xin)*NY(xin); ips++){ */
		dzero(P(extra->xopd, ips));
		real r0i=extra->r0*pow(extra->wt[ips], -3./5.);
		dembed_locstat(&P(extra->xopd, ips), 0, P(extra->xloc, ips), P(P(xin, ips)),
			alpha*extra->scale, 0);
		fractal_inv(P(extra->xopd, ips),
			P(extra->xloc, ips)->dx, r0i, extra->L0, extra->ninit);
		fractal_inv_trans(P(extra->xopd, ips),
			P(extra->xloc, ips)->dx, r0i, extra->L0, extra->ninit);
		dembed_locstat(&P(extra->xopd, ips), 1, P(extra->xloc, ips), P(P(*xout, ips)), 1, 1);
	}
}

/**
   Removing Tip/Tilt/Focus from LGS grads. TTF is the Tip/tilt/focus modes, and
   PTTF is the pseudo inverse of it, weighted by subaperture noise.  It is not
   diagonal if differential focus is removed. So cannot be in Tomo_nea*/
void TTFR(dcell* x, const dcell* TTF, const dcell* PTTF){
	if(!TTF||!PTTF){
		return;
	}
	dcell* junk=NULL;
	dcellmm(&junk, PTTF, x, "nn", 1);
	dcellmm(&x, TTF, junk, "nn", -1);
	dcellfree(junk);
}
/**
   Removing Tip/Tilt/Focus from LGS grads. TTF is the Tip/tilt/focus modes, and
   PTTF is the pseudo inverse of it, weighted by subaperture noise.  It is not
   diagonal if differential focus is removed. So cannot be in Tomo_nea*/
void TTFRt(dcell* x, const dcell* TTF, const dcell* PTTF){
	if(!TTF||!PTTF){
		return;
	}
	dcell* junk=NULL;
	dcellmm(&junk, TTF, x, "tn", 1);
	dcellmm(&x, PTTF, junk, "tn", -1);
	dcellfree(junk);
}

/**
   Apply weighting W0/W1 to a vector. W0*x-W1*(W1'*x)
*/
static void applyWeach(dmat* xin, const dsp* W0, const dmat* W1, const real wt){
	if(!W0||!W1){
		warning("W0 or W1 is NULL\n");
		return;
	}
	dmat* xout=NULL;
	dmat* tmp=NULL;
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
void applyW(dcell* xin, const dsp* W0, const dmat* W1, const real* wt){
	const int nevl=NX(xin);
	for(int iy=0; iy<NY(xin); iy++){
		for(int ievl=0; ievl<nevl; ievl++){
			int ind=iy*nevl+ievl;
			applyWeach(P(xin, ind), W0, W1, wt[ievl]);
		}
	}
}

/**
   Compute W0/W1 weighted dot product: \f$A^T(W0 B-W1 (W1^T B))\f$
*/
dcell* calcWmcc(const dcell* A, const dcell* B, const dsp* W0,
	const dmat* W1, const dmat* wt){

	assert(NX(wt)==NX(B)&&NY(wt)==1&&NX(A)==NX(B));
	const int nevl=NX(B);
	dcell* res=dcellnew(NY(A), NY(B));
	for(int iy=0; iy<NY(B); iy++){
		for(int ievl=0; ievl<nevl; ievl++){
			int ind=iy*nevl+ievl;
			dmat* xout=NULL;
			dmat* tmp=NULL;
			dspmm(&xout, W0, P(B, ind), "nn", P(wt, ievl));
			dmm(&tmp, 0, W1, P(B, ind), "tn", -1);
			dmm(&xout, 1, W1, tmp, "nn", P(wt, ievl));
			for(int ix=0; ix<NY(A); ix++){
				dmm(&P(res, ix, iy), 1, P(A, ievl, ix), xout, "tn", 1);
			}
			dfree(xout);
			dfree(tmp);
		}
	}
	return res;
}

typedef struct Tomo_T{
	const recon_t* recon;
	const real alpha;
	const dcell* xin;/*input */
	dcell* gg;/*intermediate gradient */
	dcell* xout;/*output */
	int isR;/*Mark that we are right hand side.*/
	unsigned int ips;//counter for CALL
	unsigned int iwfs;//counter for CALL
}Tomo_T;

/**
   Speed up TomoL by gathering the first part of operations (GP*HXW) belonging
   to each WFS to facilitate threading.

   gg  = GP * HXW * xin
*/
static void Tomo_prop_do(Tomo_T *data){
	const recon_t* recon=data->recon;
	const parms_t* parms=global->parms;
	sim_t* simu=global->simu;
	const int nps=recon->npsr;
	int iwfs=0;
	while((iwfs=atomic_fetch_add(&data->iwfs, 1))<parms->nwfsr){
		int ipowfs=parms->wfsr[iwfs].powfs;
		if(parms->powfs[ipowfs].skip) continue;
		const real delay=parms->sim.dt*(parms->powfs[ipowfs].dtrat+1+parms->sim.alhi);
		dmat* xx=dnew(recon->ploc->nloc, 1);
		const real hs=parms->wfs[iwfs].hs;
		const real hc=parms->wfs[iwfs].hc;
		for(int ips=0; ips<nps; ips++){
			if(parms->tomo.square&&!parms->dbg.tomo_hxw){
			/*Do the ray tracing instead of using HXW. */
				real ht=P(recon->ht, ips);
				real displace[2];
				displace[0]=parms->wfsr[iwfs].thetax*ht+parms->wfsr[iwfs].misreg_x;
				displace[1]=parms->wfsr[iwfs].thetay*ht+parms->wfsr[iwfs].misreg_y;
				if(parms->tomo.predict){
					int ips0=P(parms->atmr.indps, ips);
					displace[0]+=P(simu->atm, ips0)->vx*delay;
					displace[1]+=P(simu->atm, ips0)->vy*delay;
				}
				real scale=1.-(ht-hc)/hs;
				if(scale<0) continue;
				map_t xmap;/*make a temporary xmap for thread safety.*/
				memcpy(&xmap, P(recon->xmap, ips), sizeof(map_t));
				xmap.p=P(P(data->xin, ips));
				prop_grid_stat(&xmap, recon->ploc->stat, P(xx), 1,
					displace[0], displace[1], scale, 0, 0, 0);
			} else{
				dspcell* HXW=recon->HXWtomo/*PDSPCELL*/;
				dspmm(&xx, P(HXW, iwfs, ips), P(data->xin, ips), "nn", 1);
			}
		}
		/*Apply the gradient operation */
		dspmm(&P(data->gg, iwfs), P(recon->GP, iwfs), xx, "nn", 1);
		dfree(xx);
		/* For each wfs, Ray tracing takes 1.5 ms.  GP takes 0.7 ms. */
	}
}
/**
   Wrapper of Tomo_prop_do
*/
void Tomo_prop(Tomo_T* data, int nthread){
	//thread_t info[nthread];
	//thread_prep(info, 0, NX(data->gg), nthread, Tomo_prop_do, data);
	//CALL_THREAD(info, 1);
	data->iwfs=0;
	nthread=1;
	CALL(Tomo_prop_do, data, nthread, 1);
}
/**
   Speed up TomoL by gathering the second part of operations (GP') belonging to
   each WFS to facilitate threading.

   gg = GP' * NEAI * gg;
*/
static void Tomo_nea_gpt_do(Tomo_T *data){
	const recon_t* recon=data->recon;
	const parms_t *parms=global->parms;
	dspcell* NEAI=recon->saneai/*PDSPCELL*/;
	int iwfs=0;
	while((iwfs=atomic_fetch_add(&data->iwfs, 1))<parms->nwfsr){
		dmat* gg2=NULL;
		/*Apply the gradient operation */
		dspmm(&gg2, P(NEAI, iwfs, iwfs), P(data->gg, iwfs), "nn", 1);
		dfree(P(data->gg, iwfs)); /*We reuse gg. */
		dspmm(&P(data->gg, iwfs), P(recon->GP, iwfs), gg2, "tn", data->alpha);
		dfree(gg2);
	}
}

static void Tomo_nea_do(Tomo_T *data){
	const recon_t* recon=data->recon;
	const parms_t *parms=global->parms;
	dspcell* NEAI=recon->saneai/*PDSPCELL*/;
	int iwfs=0;
	while((iwfs=atomic_fetch_add(&data->iwfs, 1))<parms->nwfsr){
		dmat* gg2=NULL;
		/*Apply the gradient operation */
		dspmm(&gg2, P(NEAI, iwfs, iwfs), P(data->gg, iwfs), "nn", 1);
		dcp(&P(data->gg, iwfs), gg2);
		dfree(gg2);
	}
}

/*Wrapp of Tomo_nea_do for multi-threads*/
void Tomo_nea(Tomo_T* data, int nthread, int gpt){
	//thread_t info[nthread];
	//thread_prep(info, 0, NX(data->gg), nthread, gpt?Tomo_nea_gpt_do:Tomo_nea_do, data);
	//CALL_THREAD(info, 1);
	data->iwfs=0;
	nthread=1;
	if(gpt){
		CALL(Tomo_nea_gpt_do, data, nthread, 1);
	}else{
		CALL(Tomo_nea_do, data, nthread, 1);
	}
}
/**
   Speed up TomoL by gathering the third part of operations (GP') belonging to
   each WFS to facilitate threading. gg->xout.

   xout = Cxx^-1 * xin + HXW * gg;
*/
static void Tomo_iprop_do(Tomo_T *data){
	const recon_t* recon=data->recon;
	const parms_t* parms=global->parms;
	sim_t* simu=global->simu;
	map_t xmap;
	int ips;
	while((ips=atomic_fetch_add(&data->ips, 1))<NX(data->xout)){
		if(parms->tomo.square&&!parms->dbg.tomo_hxw){
			/*Do the ray tracing instead of using HXW. */
			if(!P(data->xout, ips)){
				P(data->xout, ips)=dnew(P(recon->xloc, ips)->nloc, 1);
			}
			memcpy(&xmap, P(recon->xmap, ips), sizeof(map_t));
			xmap.p=P(P(data->xout, ips));
			real ht=P(recon->ht, ips);
			for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
				if(!P(data->gg, iwfs)) continue;
				const real hs=parms->wfs[iwfs].hs;
				const real hc=parms->wfs[iwfs].hc;
				const int ipowfs=parms->wfs[iwfs].powfs;
				real displace[2];
				displace[0]=parms->wfsr[iwfs].thetax*ht+parms->wfsr[iwfs].misreg_x;
				displace[1]=parms->wfsr[iwfs].thetay*ht+parms->wfsr[iwfs].misreg_y;
				if(parms->tomo.predict){
					const real delay=parms->sim.dt*(parms->powfs[ipowfs].dtrat+1+parms->sim.alhi);
					int ips0=P(parms->atmr.indps, ips);
					displace[0]+=P(simu->atm, ips0)->vx*delay;
					displace[1]+=P(simu->atm, ips0)->vy*delay;
				}
				real scale=1.-(ht-hc)/hs;
				if(scale<0) continue;
				prop_grid_stat_transpose(&xmap, recon->ploc->stat, P(P(data->gg, iwfs)), 1,
					displace[0], displace[1], scale, 0, 0, 0);
			}
		} else{
			dspcell* HXW=recon->HXWtomo/*PDSPCELL*/;
			for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
				dspmm(&P(data->xout, ips), P(HXW, iwfs, ips), P(data->gg, iwfs), "tn", 1);
			}
		}
		if(data->xin){/*data->xin is empty when called from TomoR */
			switch(recon->cxxalg){
			case 0:{/*L2 */
				dmat* xx=NULL;
				dspmm(&xx, P(recon->L2, ips, ips), P(data->xin, ips), "nn", 1);
				dspmm(&P(data->xout, ips), P(recon->L2, ips, ips), xx, "tn", data->alpha);
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
				dspmm(&P(data->xout, ips), P(recon->ZZT, ips, ips), P(data->xin, ips), "tn", data->alpha);
			}
		}
	}
}
/**
   Wrapper of Tomo_iprop_do
 */
void Tomo_iprop(Tomo_T* data, int nthread){
	//thread_t info[nthread];
	//thread_prep(info, 0, NX(data->xout), nthread, Tomo_iprop_do, data);
	//CALL_THREAD(info, 1);
	data->ips=0;
	nthread=1;
	CALL(Tomo_iprop_do, data, nthread, 1);
}
/**
   Apply tomography right hand operator without using assembled matrix. Fast and
   saves memory.  The operation is the same as the Tomo_nea and Tomo_iprop in
   TomoL, so merge the implemenations.

   xout=HXW'*GP'*NEAI*(1-TTF*PTTF)*gin.
*/

void TomoR(dcell** xout, const void* A,
	const dcell* gin, const real alpha){
	TIC_tm;tic_tm;
	const recon_t* recon=(const recon_t*)A;
	dcell* gg=NULL;
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
void TomoRt(dcell** gout, const void* A,
	const dcell* xin, const real alpha){
	const recon_t* recon=(const recon_t*)A;
	if(!*gout){
		*gout=dcellnew(NX(recon->saneai), 1);
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

void TomoL(dcell** xout, const void* A,
	const dcell* xin, const real alpha){
	TIC_tm;tic_tm;
	const recon_t* recon=(const recon_t*)A;
	const parms_t* parms=global->parms;
	assert(NY(xin)==1);/*modify the code for ny>1 case. */
	dcell* gg=dcellnew(parms->nwfsr, 1);
	if(!*xout){
		*xout=dcellnew(recon->npsr, 1);
	}
	Tomo_T data={recon, alpha, xin, gg, *xout, 0};

	Tomo_prop(&data, recon->nthread);

	if(!parms->recon.split||parms->tomo.splitlrt){
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
void FitR(dcell** xout, const void* A,
	const dcell* xin, const real alpha){
	const fit_t* fit=(const fit_t*)A;
	const int nfit=NX(fit->thetax);
	dcell* xp=dcellnew(nfit, 1);

	if(!xin){/*xin is empty. We will trace rays from atmosphere directly */
		const parms_t* parms=global->parms;
		sim_t* simu=global->simu;
		int isim=fit->notrecon?simu->wfsisim:simu->reconisim;
		const real atmscale=simu->atmscale?P(simu->atmscale, isim):1;
		for(int ifit=0; ifit<nfit; ifit++){
			real hs=P(fit->hs, ifit);
			P(xp, ifit)=dnew(fit->floc->nloc, 1);
			for(int ips=0; ips<parms->atm.nps; ips++){
				const real ht=P(parms->atm.ht, ips)-fit->floc->ht;
				real scale=1-ht/hs;
				real displace[2];
				if(scale<0) continue;
				displace[0]=P(fit->thetax, ifit)*ht-P(simu->atm, ips)->vx*isim*parms->sim.dt;
				displace[1]=P(fit->thetay, ifit)*ht-P(simu->atm, ips)->vy*isim*parms->sim.dt;
				prop_grid(P(simu->atm, ips), fit->floc, P(P(xp, ifit)),
					atmscale, displace[0], displace[1], scale, 1, 0, 0);
			}
		}
	} else if(fit->HXF){
		dspcellmm(&xp, fit->HXF, xin, "nn", 1.);
	} else{/*Do the ray tracing from xloc to ploc */
		const int npsr=NX(fit->xloc);
		for(int ifit=0; ifit<nfit; ifit++){
			real hs=P(fit->hs, ifit);
			P(xp, ifit)=dnew(fit->floc->nloc, 1);
			for(int ips=0; ips<npsr; ips++){
				const real ht=P(fit->xloc, ips)->ht-fit->floc->ht;
				real scale=1-ht/hs;
				real displace[2];
				if(scale<0) continue;
				displace[0]=P(fit->thetax, ifit)*ht;
				displace[1]=P(fit->thetay, ifit)*ht;
				prop_nongrid(P(fit->xloc, ips), P(P(xin, ips)), fit->floc,
					P(P(xp, ifit)), 1, displace[0], displace[1], scale, 0, 0);
			}
		}
	}
	//writebin(xp, "CPU_FitR_x1");
	applyW(xp, fit->W0, fit->W1, P(fit->wt));
	//writebin(xp, "CPU_FitR_x2");
	//dcellzero(xp); P(P(xp,0), PN(xp,0)-1)=1e-7;
	dspcellmm(xout, fit->HA, xp, "tn", alpha);
	//writebin(*xout, "CPU_FitR_x3");
	dcellfree(xp);
}
/**
   Apply fit left hand side matrix in CG mode without using assembled
   matrix. Slow. don't use. Assembled matrix is faster because of multiple
   directions.  */
void FitL(dcell** xout, const void* A,
	const dcell* xin, const real alpha){
	const fit_t* fit=(const fit_t*)A;
	dcell* xp=NULL;
	dspcellmm(&xp, fit->HA, xin, "nn", 1.);
	applyW(xp, fit->W0, fit->W1, P(fit->wt));
	dspcellmm(xout, fit->HA, xp, "tn", alpha);
	dcellfree(xp);xp=NULL;
	dcellmm(&xp, fit->NW, xin, "tn", 1);
	dcellmm(xout, fit->NW, xp, "nn", alpha);
	dcellfree(xp);
	if(fit->actslave){
		dspcellmm(xout, fit->actslave, xin, "nn", alpha);
	}
}


/**
 * Convert column vector nsax2 or nsax3 to sparse matrix. The third column of nea is coupling.
 * ll: controls lower left
 * ur: controls upper right 
 */
dsp* nea2sp(dmat* nea, int ll, int ur, int ng){
	if(NY(nea) != ng && NY(nea)!=3){
		error("nea has wrong format:%ldx%ld. Expect 2 or 3 columns.\n", NX(nea), NY(nea));
	}
	if(NY(nea)==ng){//off-diagonal not available
		ll=0;
		ur=0;
	}
	const long nsa=NX(nea);
	dsp* sanea=dspnew(nsa*ng, nsa*ng, (ng+(ll?1:0)+(ur?1:0))*nsa);
	spint* pp=sanea->pp;
	spint* pi=sanea->pi;
	real* px=sanea->px;
	long count=0;
	for(long ig=0; ig<ng; ig++){
		for(long isa=0; isa<nsa; isa++){
			pp[isa+nsa*ig]=count;
			if(ig==1 && ur){/*Cxy */
				pi[count]=isa;
				px[count]=P(nea, isa, 2);
				count++;
			}
			{/*Cxx */
				pi[count]=isa+nsa*ig;
				px[count]=P(nea, isa, ig);
				count++;
			}
			if(ig==0 && ll){/*Cyx */
				pi[count]=isa+nsa;
				px[count]=P(nea, isa, 2);
				count++;
			}
		}
	}
	pp[nsa*ng]=count;
	if(count>sanea->nzmax){
		error("memory overflow\n");
	}
	return sanea;
}

/**
   Compute and save PSF reconstruction telemetry.  For pseudo open loop
   estimations, like high order loop, we add opdr to, and subtract dmpsol from
   the OPD.  For closed loop estimations, like ahst low order or lsr, we add
   dmerr_lo, and dmerr to the OPD.*/
void psfr_calc(sim_t* simu, dcell* opdr, dcell* dmpsol, dcell* dmerr, dcell* dmerr_lo){
	const parms_t* parms=simu->parms;
	recon_t* recon=simu->recon;
	/* The tomography estimates, opdr is pseudo open loop estimates. We need to
	   subtract the constribution of the added DM command to form closed loop
	   estimates.
	*/
	dcell* dmadd=NULL;
	if(opdr&&parms->dbg.useopdr){
	/* The original formulation using Hx*x-Ha*a.  Changed from ploc to plocs
	   on July 18, 2011. Ploc is too low sampled. Make sure the sampling of
	   plocs is not too big.  Deprecated. Use dm space instead.
	*/
		if(dmpsol){/*Pseudo OL estimates */
			if(parms->recon.split==1){
			/* We will remove NGS modes from dmlast which is in NULL modes of
			   tomography reconstructor (is this 100% true)?  SHould we remove NGS
			   modes from final OPD, xx, instead?*/
				dcell* tmp=dcelldup(dmpsol);/*The DM command used for high order. */
				remove_dm_ngsmod(simu, tmp);/*remove NGS modes as we do in ahst. */
				dcelladd(&dmadd, 1, tmp, -1);
				dcellfree(tmp);
			} else{
				dcelladd(&dmadd, 1, dmpsol, -1);
			}
		}

		loc_t* locs=simu->aper->locs;
		dmat* xx=dnew(locs->nloc, 1);
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			real hs=P(parms->evl.hs, ievl);
			if(P(parms->evl.psfr, ievl)){
				dzero(xx);
				if(opdr){
					map_t xmap;
					const int npsr=recon->npsr;
					/*First compute residual opd: Hx*x-Ha*a*/
					for(int ips=0; ips<npsr; ips++){
						const real ht=P(recon->ht, ips);
						real scale=1-ht/hs;
						real dispx=P(parms->evl.thetax, ievl)*ht;
						real dispy=P(parms->evl.thetay, ievl)*ht;
						if(scale<0) continue;
						if(parms->tomo.square){/*square xloc */
							memcpy(&xmap, P(recon->xmap, ips), sizeof(map_t));
							xmap.p=P(P(opdr, ips));
							prop_grid_stat(&xmap, locs->stat, P(xx), 1,
								dispx, dispy, scale, 0, 0, 0);
						} else{
							prop_nongrid(P(recon->xloc, ips), P(P(opdr, ips)), locs,
								P(xx), 1, dispx, dispy, scale, 0, 0);
						}
					}
				}
				if(dmadd){
					for(int idm=0; idm<parms->ndm; idm++){
						const real ht=parms->dm[idm].ht;
						real scale=1.-ht/hs;
						real dispx=P(parms->evl.thetax, ievl)*ht;
						real dispy=P(parms->evl.thetay, ievl)*ht;
						if(scale<0) continue;
						prop_nongrid(P(recon->aloc, idm), P(P(dmadd, idm)), locs,
							P(xx), 1, dispx, dispy, scale, 0, 0);
					}
				}
				dmm(&P(simu->ecov, ievl), 1, xx, xx, "nt", 1);
				if(parms->dbg.ecovxx){
					zfarr_push(simu->save->ecovxx[ievl], simu->reconisim, xx);
				}
			}/*if psfr[ievl] */
		}/*ievl */
		dfree(xx);
	} else{/*Do Ha in postproc, so just do a. */
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
void shift_grad(sim_t* simu){
	const parms_t* parms=simu->parms;
	if(parms->sim.evlol||parms->sim.idealfit||parms->sim.idealtomo) return;
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
		} else{
			long nnx[parms->nwfsr];
			for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
				int ipowfs=parms->wfsr[iwfs].powfs;
				nnx[iwfs]=simu->powfs[ipowfs].saloc->nloc*2;
			}
			simu->gradlastcl=dcellnew3(parms->nwfsr, 1, nnx, NULL);
		}
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			const real scale=1./parms->powfs[ipowfs].nwfs;
			for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs, indwfs);
				dadd(&P(simu->gradlastcl, ipowfs), 1., P(simu->gradcl, iwfs), scale);
			}
		}
	} else{
		dcellcp(&simu->gradlastcl, simu->gradcl);
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
lmat* loc_coord2ind(loc_t* aloc,       /**<[in] Aloc*/
	const char* fndead /**<[in] File containing dead actuators*/
){
	dmat* dead=dread("%s", fndead);
	if(NY(dead)!=2&&NY(dead)!=3){
		error("%s must contain 2 or 3 columns of data\n", fndead);
	}
	loc_create_map(aloc);
	map_t* map=aloc->map;
	real ox=aloc->map->ox;
	real oy=aloc->map->oy;
	real dx1=1./aloc->dx;
	dmat* ps=dead/*PDMAT*/;
	lmat* out=lnew(aloc->nloc, 1);
	for(long jact=0; jact<NX(dead); jact++){
		long mapx=(long)round((P(ps, jact, 0)-ox)*dx1);
		long mapy=(long)round((P(ps, jact, 1)-oy)*dx1);
		long iact=loc_map_get(map, mapx, mapy)-1;
		if(iact>=0){
			P(out, iact)=(NY(dead)==3?P(ps, jact, 2)*1e9:1);//integer in nm.
		}
	}
	dfree(dead);
	return out;
}

/**
   Prepare arrays for cn2 estimation. Multiple pairs can be used to do Cn2
Estimation. The result will be an average of them.  */


cn2est_t* cn2est_prepare(const parms_t* parms, const powfs_t* powfs){
	dmat* pair=parms->cn2.pair;
	int npair=NX(pair)*NY(pair);
	int ipowfs=-1;
	for(int ind=0; ind<npair; ind++){
		int iwfs=(int)P(pair, ind);
		if(ipowfs==-1){
			ipowfs=parms->wfs[iwfs].powfs;
		} else if(ipowfs!=parms->wfs[iwfs].powfs){
			error("All wfs in parms->cn2.pair do not belong to the same powfs\n");
		}
	}
	dmat* wfstheta=dnew(parms->nwfs, 2);
	dmat* ptheta=wfstheta/*PDMAT*/;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		P(ptheta, iwfs, 0)=parms->wfs[iwfs].thetax;
		P(ptheta, iwfs, 1)=parms->wfs[iwfs].thetay;
	}
	dmat* ht=0;
	if(parms->cn2.keepht){
		ht=dnew(parms->atmr.nps, 1);
		for(int iht=0; iht<NX(ht); iht++){
			P(ht, iht)=P(parms->atmr.ht, iht);
		}
	} else{
		int nht=parms->cn2.nhtomo;
		ht=dnew(nht, 1);
		real dht=parms->cn2.hmax/(real)(nht-1);
		if(dht<=0) error("Invalid hmax or nhtomo\n");
		for(int iht=0; iht<nht; iht++){
			P(ht, iht)=dht*iht;
		}
	}
	dmat* hs=dnew(parms->nwfs, 1);
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		P(hs, iwfs)=parms->wfs[iwfs].hs;
	}
	cn2est_t* cn2est=cn2est_new(pair, wfstheta, powfs[ipowfs].saloc, powfs[ipowfs].saa, parms->cn2.saat,
		hs, ht, parms->cn2.keepht, parms->atmr.L0);
	dfree(hs);
	dfree(wfstheta);
	if(cn2est){
		cn2est->os=dnew(NX(ht), 1);
		if(!parms->cn2.keepht){
			/*preserve the number of over sampled layers. */
			int osc=0;
			for(int ips=0; ips<parms->atmr.nps; ips++){
				if(P(parms->atmr.os, ips)>1){
					osc++;
					if(P(parms->atmr.os, ips)!=2){
						error("os is not 2. adept this code to it.\n");
					}
				}
			}
			for(int iht=0; iht<NX(ht); iht++){
				if(iht<osc){
					P(cn2est->os, iht)=2;
				} else{
					P(cn2est->os, iht)=1;
				}
			}
		} else{
			for(int iht=0; iht<NX(ht); iht++){
				P(cn2est->os, iht)=P(parms->atmr.os, iht);
			}
		}
		if(parms->cn2.verbose){
			cn2est->dmht=dnew(parms->ndm, 1);
			for(int idm=0; idm<parms->ndm; idm++){
				P(cn2est->dmht, idm)=parms->dm[idm].ht;
			}
		}
		if(parms->save.setup){
			writebin(cn2est->overlapi, "cn2_overlapi");
			writebin(cn2est->iPnk, "cn2_iPnk");
			writebin(cn2est->Pnk, "cn2_Pnk");
			writebin(cn2est->ht, "cn2_ht");
			writebin(cn2est->wtconvert, "cn2_wtconvert");
		}
	}
	dfree(ht);
	return cn2est;
}
/**
   Implemented mechanism to move height of layers.
 */
static void cn2est_moveht(recon_t* recon){
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
void cn2est_isim(dcell* cn2res, recon_t* recon, const parms_t* parms, const dcell* grad, int* tomo_update){
	cn2est_t* cn2est=recon->cn2est;
	cn2est_push(cn2est, grad);
	static int icn2=-1;
	if(cn2est->count%parms->cn2.step==0){
		icn2++;
		int nset=cn2est->count/parms->cn2.step;
		cn2est_est(cn2est, parms->cn2.verbose);/*do the CN2 estimation */
		if(parms->cn2.reset&&(nset%parms->cn2.reset)==0){
			cn2est_reset(cn2est);
		}

		if(cn2res){
			P(P(cn2res, 0), icn2)=cn2est->r0m;
			memcpy(PCOL(P(cn2res, 1), icn2), P(P(cn2est->wtrecon, 0)),
				NX(cn2est->htrecon)*sizeof(real));
		}
		if(parms->cn2.tomo){
			if(parms->cn2.moveht){
				cn2est_moveht(recon);
			}

			if(parms->cn2.verbose){
				info2("Updating tomography weights\n");
			}
			/*Changes recon parameters. cannot be parallel with tomofit(). */
			/*wtrecon is referenced so should be updated automaticaly. */
			if(P(recon->wt)!=P(P(cn2est->wtrecon, 0))){
				dfree(recon->wt);
				recon->wt=dref(P(cn2est->wtrecon, 0));
			}
			recon->r0=cn2est->r0m;
			recon->L0=cn2est->L0;

			*tomo_update=1;
		}
	}
}
