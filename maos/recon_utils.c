/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "ahst.h"
#include "sim.h"
/**
   \file recon_utils.c
   Reusable utilities for wavefront reconstruction and DM fitting.
*/

/**
   Apply Laplacian2 to xin and accumulates to xout.
*/
void apply_L2(dcell **xout, const spcell *L2, const dcell *xin, 
	      double alpha, int nthread){
    dcell *xx=NULL;
    spcellmulmat_thread(&xx, L2, xin, 1.);
    sptcellmulmat_thread(xout, L2, xx, alpha);
    dcellfree(xx);
}
/**
   Apply turbulence invpsd to xin in Fourier space, scaled by alpha and add to xout.
   do nothing if xb != yb, since we apply to diagonal only.
   if xb==-1, do all blocks. 
*/
void apply_invpsd(dcell **xout, const void *A, const dcell *xin, double alpha, int xb, int yb){
    if(xb!=yb) return;
    const INVPSD_T *extra=A;
    dcell *invpsd=extra->invpsd;
    ccell *fftxopd=extra->fftxopd;
    int ips1, ips2;
    if(xb<0){//do all cells
	ips1=0; 
	ips2=xin->nx*xin->ny;
    }else{//do the specified cell
	ips1=xb;
	ips2=xb+1;
    }
    for(int ips=ips1; ips<ips2; ips++){
	//for(int ips=0; ips<xin->nx*xin->ny; ips++){
	long nx=fftxopd->p[ips]->nx;
	long ny=fftxopd->p[ips]->ny;
	if(extra->square){
	    dmat *xini=dref_reshape(xin->p[ips], nx, ny);
	    ccpd(&fftxopd->p[ips], xini);
	    dfree(xini);
	}else{
	    czero(fftxopd->p[ips]);
	    cembed_locstat(&fftxopd->p[ips], 0, extra->xloc[ips], xin->p[ips]->p, alpha, 0);
	}
	cfft2(fftxopd->p[ips],-1);
	//cwrite(fftxopd->p[ips],"fftxopd_%d",ips);
	ccwmd(fftxopd->p[ips], invpsd->p[ips], 1);
	cfft2(fftxopd->p[ips],1);
	if(extra->square){
	    dmat *xouti=NULL;
	    xouti=dref_reshape((*xout)->p[ips], nx, ny);
	    creal2d(&xouti,1,fftxopd->p[ips],alpha);
	    dfree(xouti);
	}else{
	    cembed_locstat(&fftxopd->p[ips], 1, extra->xloc[ips], (*xout)->p[ips]->p, 1, 1);
	}
    }
}

/**
   Apply fractal regularization to x, scaled by alpha.
   do nothing if xb != yb, since we apply to diagonal only.
   if xb==-1, do all blocks. 
*/
void apply_fractal(dcell **xout, const void *A, const dcell *xin, double alpha, int xb, int yb){
    if(xb!=yb) return;
    const FRACTAL_T *extra=A;
    int ips1, ips2;
    if(xb<0){//do all cells
	ips1=0; 
	ips2=xin->nx*xin->ny;
    }else{//do the specified cell
	ips1=xb;
	ips2=xb+1;
    }
    for(int ips=ips1; ips<ips2; ips++){
	//for(int ips=0; ips<xin->nx*xin->ny; ips++){
	dzero(extra->xopd->p[ips]);
	double r0i=extra->r0*pow(extra->wt[ips], -3./5.);
	dembed_locstat(&extra->xopd->p[ips], 0, extra->xloc[ips], xin->p[ips]->p, 
		       alpha*extra->scale, 0);
	fractal_inv(extra->xopd->p[ips]->p, extra->xopd->p[ips]->nx, extra->xopd->p[ips]->ny, 
		    extra->xloc[ips]->dx, r0i, extra->l0, extra->ninit);
	fractal_inv_trans(extra->xopd->p[ips]->p, extra->xopd->p[ips]->nx, extra->xopd->p[ips]->ny, 
			  extra->xloc[ips]->dx, r0i, extra->l0, extra->ninit);
	dembed_locstat(&extra->xopd->p[ips], 1, extra->xloc[ips], (*xout)->p[ips]->p, 1, 1);
    }
}

/**
   Removing Tip/Tilt/Focus from LGS grads. TTF is the Tip/tilt/focus modes, and
   PTTF is the pseudo inverse of it, weighted by subaperture noise.  */
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
   Apply weighting W0/W1 to a vector. W0*x-W1*(W1'*x)
*/
static void applyWeach(dmat *xin, const dsp *W0, const dmat *W1, const double wt){
    if(!W0 || !W1) {
	warning("W0 or W1 is NULL\n");
	return;
    }
    dmat *xout=NULL;
    dmat *tmp=NULL;
    spmulmat(&xout, W0, xin, wt);
    dmm(&tmp, W1, xin, "tn", -1);
    dmm(&xout,W1, tmp, "nn", wt);
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
   Compute W0/W1 weighting dot product: \f$A^T(W0 B-W1 (W1^T B))\f$
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
	    spmulmat(&xout, W0, B->p[ind], wt->p[ievl]);
	    dmm(&tmp, W1, B->p[ind], "tn", -1);
	    dmm(&xout, W1, tmp, "nn", wt->p[ievl]);
	    for(int ix=0; ix<A->ny; ix++){
		dmm(&res->p[ix+iy*res->nx], A->p[ix*nevl+ievl], xout, "tn", 1);
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
    const dcell *xin;//input
    dcell *gg;//intermediate gradient
    dcell *xout;//output
}Tomo_T;

/**
   Speed up TomoL by gathering the first part of operations (GP*HXW) belonging
   to each WFS to facilitate threading. 

   gg  = GP * HXW * xin
*/
#define USE_PROP 1
static void Tomo_prop(thread_t *info){
    Tomo_T *data=info->data;
    const RECON_T *recon=data->recon;
    const PARMS_T *parms=recon->parms;
    PSPCELL(recon->HXWtomo,HXW);
    const int nps=recon->HXWtomo->ny;
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
	int ipowfs = parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	dmat *xx=dnew(recon->ploc->nloc, 1);
	const double hs=parms->powfs[ipowfs].hs;
	for(int ips=0; ips<nps; ips++){
	    if(parms->tomo.square && !parms->dbg.tomo_hxw){
		//Do the ray tracing instead of using HXW.
		double ht=recon->ht->p[ips];
		double displace[2];
		displace[0]=parms->wfsr[iwfs].thetax*ht;
		displace[1]=parms->wfsr[iwfs].thetay*ht;
		double scale=1. - ht/hs;
		prop_grid_stat(recon->xmap[ips], recon->ploc->stat, xx->p, 1, 
			       displace[0],displace[1], scale, 0, 0, 0);
	    }else{
		spmulmat(&xx, HXW[ips][iwfs], data->xin->p[ips], 1);
	    }
	}
	//Apply the gradient operation
	spmulmat(&data->gg->p[iwfs], recon->GP->p[iwfs], xx, 1);
	dfree(xx);
	/*
	  For each wfs, Ray tracing takes 1.5 ms.  GP takes 0.7 ms.
	*/
    }
}

/**
   Speed up TomoL by gathering the second part of operations (GP') belonging to
   each WFS to facilitate threading. 
   
   gg = GP' * NEAI * gg;
*/
static void Tomo_nea(thread_t *info){
    Tomo_T *data=info->data;
    const RECON_T *recon=data->recon;
    PSPCELL(recon->saneai, NEAI);
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
	dmat *gg2=NULL;
	//Apply the gradient operation
	spmulmat(&gg2, NEAI[iwfs][iwfs], data->gg->p[iwfs], 1);
	dfree(data->gg->p[iwfs]); //We reuse gg.
	sptmulmat(&data->gg->p[iwfs], recon->GP->p[iwfs], gg2, data->alpha);
	dfree(gg2);
    }
}
/**
   Speed up TomoL by gathering the third part of operations (GP') belonging to
   each WFS to facilitate threading. gg->xout.

   xout = Cxx^-1 * xin + HXW * gg;
*/
static void Tomo_iprop(thread_t *info){
    Tomo_T *data=info->data;
    const RECON_T *recon=data->recon;
    const PARMS_T *parms=recon->parms;
    const int nps=recon->HXWtomo->ny;
    PSPCELL(recon->HXWtomo,HXW);
    //for(int ips=0; ips<recon->HXWtomo->ny; ips++){
    for(int ips=info->start; ips<info->end; ips++){
	if(parms->tomo.square && !parms->dbg.tomo_hxw){
	    //Do the ray tracing instead of using HXW.
	    if(!data->xout->p[ips]){
		data->xout->p[ips]=dnew(recon->xloc[ips]->nloc, 1);
	    }
	    recon->xmap[ips]->p=data->xout->p[ips]->p;
	    double ht=recon->ht->p[ips];
	    for(int iwfs=0; iwfs<recon->HXWtomo->nx; iwfs++){
		if(!data->gg->p[iwfs]) continue;
		int ipowfs = parms->wfsr[iwfs].powfs;
		const double hs=parms->powfs[ipowfs].hs;
		double displace[2];
		displace[0]=parms->wfsr[iwfs].thetax*ht;
		displace[1]=parms->wfsr[iwfs].thetay*ht;
		double scale=1. - ht/hs;
		prop_grid_stat_transpose(recon->xmap[ips], recon->ploc->stat, data->gg->p[iwfs]->p, 1, 
					 displace[0],displace[1], scale, 0, 0, 0);
	    }
	}else{
	    for(int iwfs=0; iwfs<recon->HXWtomo->nx; iwfs++){
		sptmulmat(&data->xout->p[ips], HXW[ips][iwfs], data->gg->p[iwfs], 1);
	    }
	}
	if(data->xin){//data->xin is empty when called from TomoR
	    switch(recon->cxx){
	    case 0:{//L2
		dmat *xx=NULL;
		spmulmat(&xx, recon->L2->p[ips+ips*nps], data->xin->p[ips], 1);
		sptmulmat(&data->xout->p[ips], recon->L2->p[ips+ips*nps], xx, data->alpha);
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
	}
    }
}

/**
   Apply tomography right hand operator without using assembled matrix. Fast and
   saves memory.  The operation is the same as the Tomo_nea and Tomo_iprop in
   TomoL, so merge the implemenations.

   xout=HXW'*GP'*NEAI*(1-TTF*PTTF)*gin.
*/

void TomoR(dcell **xout, const void *A, 
	   const dcell *gin, const double alpha){
    
    const RECON_T *recon=(const RECON_T *)A;
    dcell *gg=NULL;
    dcellcp(&gg, gin);//copy to gg so we don't touch the input.
    TTFR(gg, recon->TTF, recon->PTTF);

    if(!*xout){
	*xout=dcellnew(recon->npsr, 1);
    }
    Tomo_T data={recon, alpha, NULL, gg, *xout};
    thread_t info_iwfs2[recon->nthread];
    thread_prep(info_iwfs2, 0, gg->nx, recon->nthread, Tomo_nea, &data);
    thread_t info_ips[recon->nthread];
    thread_prep(info_ips, 0, recon->npsr, recon->nthread, Tomo_iprop, &data);
    
    CALL_THREAD(info_iwfs2, recon->nthread, 1);
    CALL_THREAD(info_ips, recon->nthread, 1);
    dcellfree(gg);
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

void TomoL(dcell **xout, const void *A, 
	   const dcell *xin, const double alpha){
    const RECON_T *recon=(const RECON_T *)A;
    const PARMS_T *parms=recon->parms;
    assert(xin->ny==1);//modify the code for ny>1 case.
    dcell *gg=dcellnew(parms->nwfsr, 1);
    const int nps=recon->npsr;
    if(!*xout){
	*xout=dcellnew(recon->npsr, 1);
    }
    Tomo_T data={recon, alpha, xin, gg, *xout};

    if(parms->tomo.square){//do ray tracing instead of HXW.
	for(int ips=0; ips<nps; ips++){
	    recon->xmap[ips]->p=xin->p[ips]->p; //replace the vector.
	}
    }
    thread_t info_iwfs1[recon->nthread];
    thread_prep(info_iwfs1, 0, gg->nx, recon->nthread, Tomo_prop, &data);
    thread_t info_iwfs2[recon->nthread];
    thread_prep(info_iwfs2, 0, gg->nx, recon->nthread, Tomo_nea, &data);
    thread_t info_ips[recon->nthread];
    thread_prep(info_ips, 0, nps, recon->nthread, Tomo_iprop, &data);
    CALL_THREAD(info_iwfs1, recon->nthread, 1);
    if(!parms->tomo.split || parms->dbg.splitlrt){
	/*Remove global Tip/Tilt, differential focus only in integrated
	  tomography to limit noise propagation.*/
	TTFR(gg, recon->TTF, recon->PTTF);
    }
    CALL_THREAD(info_iwfs2, recon->nthread, 1);
    CALL_THREAD(info_ips, recon->nthread, 1);
   
    /* 
       square=1  square=0 (1 thread on T410s)
       iwfs1: takes 6 ms 13 ms
       iwfs2: takes 4 ms 4 ms
       ips:   takes 6 ms 9 ms
    */
    dcellfree(gg);
    
    if(recon->ZZT){//single point piston constraint. fast if any
	sptcellmulmat_thread(xout, recon->ZZT, xin, alpha);
    }
    /*Tikhonov regularization is not added because it is not necessary in CG
      mode.*/
}

/**
   Apply fit right hand side matrix in CG mode without using assembled matrix.
   Slow. don't use. Assembled matrix is faster because of multiple directions.
*/
void FitR(dcell **xout, const void *A, 
	  const dcell *xin, const double alpha){
    const RECON_T *recon=(const RECON_T *)A;
    dcell *xp=NULL;
    if(!xin){//xin is empty. We will trace rays from atmosphere directly
	const PARMS_T *parms=recon->parms;
	SIM_T *simu=recon->simu;
	int isim=parms->sim.closeloop?simu->isim-1:simu->isim;
	const int nfit=parms->fit.nfit;
	xp=dcellnew(nfit,1);
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=parms->fit.ht[ifit];
	    xp->p[ifit]=dnew(recon->ploc->nloc,1);
	    for(int ips=0; ips<parms->atm.nps; ips++){
		const double ht = parms->atm.ht[ips];
		double scale=1-ht/hs;
		double displace[2];
		displace[0]=parms->fit.thetax[ifit]*ht-simu->atm[ips]->vx*isim*simu->dt;
		displace[1]=parms->fit.thetay[ifit]*ht-simu->atm[ips]->vy*isim*simu->dt;
		prop_grid(simu->atm[ips], recon->ploc, xp->p[ifit]->p, 
			  alpha, displace[0], displace[1], scale, 1, 0, 0);
	    }
	}
    }else if(recon->HXF){
	spcellmulmat_thread(&xp, recon->HXF, xin, 1.);
    }else{//Do the ray tracing from xloc to ploc
	const PARMS_T *parms=recon->parms;
	const int nfit=parms->fit.nfit;
	const int npsr=recon->npsr;
	xp=dcellnew(nfit,1);
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=parms->fit.ht[ifit];
	    xp->p[ifit]=dnew(recon->ploc->nloc,1);
	    for(int ips=0; ips<npsr; ips++){
		const double ht = recon->ht->p[ips];
		double scale=1-ht/hs;
		double displace[2];
		displace[0]=parms->fit.thetax[ifit]*ht;
		displace[1]=parms->fit.thetay[ifit]*ht;
		prop_nongrid(recon->xloc[ips], xin->p[ips]->p, recon->ploc, NULL, 
			     xp->p[ifit]->p, alpha, displace[0], displace[1], scale, 0, 0);
	    }
	}
    }
    applyW(xp, recon->W0, recon->W1, recon->fitwt->p);
    sptcellmulmat_thread(xout, recon->HA, xp, alpha);
    dcellfree(xp);
}
/**
   Apply fit left hand side matrix in CG mode without using assembled
   matrix. Slow. don't use. Assembled matridx is faster because of multiple
   directions.  */
void FitL(dcell **xout, const void *A, 
	  const dcell *xin, const double alpha){
    const RECON_T *recon=(const RECON_T *)A;
    dcell *xp=NULL;
    spcellmulmat_thread(&xp, recon->HA, xin, 1.);
    applyW(xp, recon->W0, recon->W1, recon->fitwt->p);
    sptcellmulmat_thread(xout, recon->HA, xp, alpha);
    dcellfree(xp);xp=NULL;
    dcellmm(&xp,recon->fitNW, xin, "tn", 1);
    dcellmm(xout,recon->fitNW, xp, "nn", alpha);
    dcellfree(xp);
    if(recon->actslave){
	spcellmulmat(xout, recon->actslave, xin, alpha);
    }
}


/**
   Update LGS focus or reference vector.
   \fixme: what if dtrat!=1
   \verbatim
   LGS: CL grads works, but not as good.
   LGS: CL grads + DM grads does not work
   LGS: CL grads + DM grads - Xhat grad is the choice.

   mffocus is 0: no focus tracking.
   mffocus is 1: use CL grads + DM grads - Xhat grad for LGS and NGS
   mffocus is 2: use CL grads + DM grads - Xhat grad for LGS; 
   CL grads for NGS.
   \endverbatim
*/
void focus_tracking(SIM_T*simu){

    if(!simu->recon->RFlgs){
	warning("There is no LGS. No need to do focus tracking\n");
	return;
    }
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    dcell *graduse=dcellnew(parms->nwfs,1);
    PDSPCELL(recon->GXfocus,GX);
    int ngs_psol=0;
    int ngs_x=0;
    int lgs_psol=0;
    int lgs_x=0;
    lgs_psol=1;
    lgs_x=1;
    if(parms->sim.mffocus==1){
	ngs_psol=1;
	ngs_x=1;
    }else if(parms->sim.mffocus!=2){
	error("Invalid mffocus: %d\n", parms->sim.mffocus);
    }
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	int gs_psol, gs_x;
	if(parms->powfs[ipowfs].hasllt){//LGS
	    gs_psol=lgs_psol;
	    gs_x=lgs_x;
	}else{//NGS
	    if(parms->powfs[ipowfs].order==1) continue;//no bother with TT NGS.
	    gs_psol=ngs_psol;
	    gs_x=ngs_x;
	}
	if(gs_psol){//psol
	    if(simu->gradlastol->p[iwfs]){
		graduse->p[iwfs]=ddup(simu->gradlastol->p[iwfs]);
	    }else{
		error("Require PSOL grads for wfs %d\n",iwfs);
	    }
	}else{//cl
	    graduse->p[iwfs]=ddup(simu->gradlastcl->p[iwfs]);
	}
	if(gs_x){
	    info("Subtracing tomo grad from wfs %d\n",iwfs);
	    for(int ips=0; ips<simu->recon->npsr; ips++){
		if(!GX[ips][iwfs]){
		    error("GX[%d][%d] is empty\n",ips,iwfs);
		}
		spmulmat(&graduse->p[iwfs],GX[ips][iwfs],simu->opdr->p[ips],-1);
	    }
	}
    }
    dcell *NGSfocus=NULL;
    dcellmm(&NGSfocus, recon->RFngs, graduse, "nn", 1);
    info("NGSfocus is %g\n", NGSfocus->p[0]->p[0]);
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	if(!parms->powfs[ipowfs].hasllt) 
	    continue;
	dmat *LGSfocus=NULL;
	dmm (&LGSfocus,recon->RFlgs->p[iwfs],graduse->p[iwfs],"nn",1);
	//Use NGS focus - LGS focus to drive the zoom optics/reference vector
	dadd(&LGSfocus,-1, NGSfocus->p[0], 1);
	dadd(&simu->focuslpf->p[iwfs], 1.-parms->sim.lpfocus, 
	     LGSfocus, parms->sim.lpfocus);
	dadd(&simu->focusint->p[iwfs], 1, 
	     simu->focuslpf->p[iwfs], parms->sim.epfocus);
	dfree(LGSfocus);
    }
    dcellfree(NGSfocus);
    dcellfree(graduse);
}


/**
   Experimental routine to do wind estimation using correlation tracking of
   tomography outputs. Not finished. Not verified.  */
void windest(SIM_T *simu){
    long nps=simu->opdr->nx;
    if(!simu->opdrhat){
	simu->opdrhat=ccellnew(nps,1);
	simu->opdrhat=ccellnew(nps,1);
	for(long ips=0; ips<nps; ips++){
	    simu->opdrhat->p[ips]=cnew(simu->recon->xloc_nx[ips],
				       simu->recon->xloc_ny[ips]);
	    simu->opdrhatlast->p[ips]=cnew(simu->recon->xloc_nx[ips],
					   simu->recon->xloc_ny[ips]);
	    cfft2plan(simu->opdrhat->p[ips], -1);
	    cfft2plan(simu->opdrhatlast->p[ips], -1);
	}
	simu->windest=dnew(2, nps);
    }
    ccell *temp;
    temp=simu->opdrhatlast;
    simu->opdrhatlast=simu->opdrhat;
    simu->opdrhat=temp;
    for(long ips=0; ips<nps; ips++){
	dcomplex *restrict dest=simu->opdrhat->p[ips]->p;
	double *restrict source=simu->opdr->p[ips]->p;
	for(long ipix=0; ipix<simu->opdr->p[ips]->nx*simu->opdr->p[ips]->nx; ipix++){
	    dest[ipix]=source[ipix];
	}
	cfftshift(simu->opdrhat->p[ips]);//may be able to remove
	cfft2(simu->opdrhat->p[ips], -1);
	cfftshift(simu->opdrhat->p[ips]);
    }
    if(simu->isim > simu->parms->sim.start){
	//Compute Wind. not done yet.
	ccellwrite(simu->opdrhat,"opdrhat");
	ccellwrite(simu->opdrhatlast,"opdrhatlast");
	exit(0);
    }
}

/**
   Convert block of 2x2 neas to sparse matrix.
 */
dsp *nea2sp(dmat **nea, long nsa){
    dsp *sanea=spnew(nsa*2, nsa*2, 4*nsa);
    spint *pp=sanea->p;
    spint *pi=sanea->i;
    double *px=sanea->x;
    long count=0;
    for(long isa=0; isa<nsa; isa++){
	//Cxx
	pp[isa]=count;
	pi[count]=isa;
	px[count]=nea[isa]->p[0];
	count++;
	//Cyx
	pi[count]=isa+nsa;
	px[count]=nea[isa]->p[1];
	count++;
    }
    for(long isa=0; isa<nsa; isa++){
	//Cxy
	pp[isa+nsa]=count;
	pi[count]=isa;
	px[count]=nea[isa]->p[2];
	count++;
	//Cyy
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
   dmerr_lo, and dmerr_hi to the OPD.*/
void psfr_calc(SIM_T *simu, dcell *opdr, dcell *dmpsol, dcell *dmerr_hi, dcell *dmerr_lo){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    /*
      
      Do we do PSFR on evl directions, or separate directions?

    */
    /*
      The tomography estimates, opdr is pseudo open loop
      estimates. We need to subtract the constribution of the
      added DM command to form closed loop estimates. 
    */
    dcell *dmadd=NULL;

    if(parms->tomo.split==1){
	/*
	  We will remove NGS modes from dmlast which is in NULL
	  modes of tomography reconstructor (is this 100% true)?
	  SHould we remove NGS modes from final OPD, xx,
	  instead?*/
	dcell *tmp = dcelldup(dmpsol);//The DM command used for high order.
	remove_dm_ngsmod(simu, tmp);//remove NGS modes as we do in ahst.
	dcelladd(&dmadd, 1, tmp, -1);
	dcellfree(tmp);
    }else{
	dcelladd(&dmadd, 1, dmpsol, -1);
    }
    if(dmerr_hi){/*high order closed loop estimates. (lsr)*/
	dcelladd(&dmadd, 1, dmerr_hi, 1);
    }
    if(dmerr_lo){/*In AHST, dmerr_lo is CL Estimation.*/
	addlow2dm(&dmadd, simu, dmerr_lo, 1);
    }
    dmat *xx = dnew(recon->ploc->nloc, 1);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	double hs = parms->evl.ht[ievl];
	if(parms->evl.psfr[ievl]){
	    dzero(xx);
	    if(opdr){
		const int npsr=recon->npsr;
		/*First compute residual opd: Hx*x-Ha*a*/
		for(int ips=0; ips<npsr; ips++){
		    const double ht = recon->ht->p[ips];
		    double scale=1-ht/hs;
		    double dispx=parms->evl.thetax[ievl]*ht;
		    double dispy=parms->evl.thetay[ievl]*ht;
		    if(parms->tomo.square){//square xloc
			recon->xmap[ips]->p=opdr->p[ips]->p;
			prop_grid_stat(recon->xmap[ips], recon->ploc->stat, xx->p, 1, 
				       dispx, dispy, scale, 0, 0, 0);
		    }else{
			prop_nongrid(recon->xloc[ips], opdr->p[ips]->p, recon->ploc, NULL,
				     xx->p, 1, dispx, dispy, scale, 0, 0);
		    }
		}
	    }
	    if(dmadd){
		for(int idm=0; idm<parms->ndm; idm++){
		    const double ht = parms->dm[idm].ht;
		    double scale=1-ht/hs;
		    double dispx=parms->evl.thetax[ievl]*ht;
		    double dispy=parms->evl.thetay[ievl]*ht;
		    prop_nongrid(recon->aloc[idm], dmadd->p[idm]->p, recon->ploc, NULL,
				 xx->p, 1, dispx, dispy, scale, 0, 0);
		}
	    }
	    dmm(&simu->ecov->p[ievl], xx, xx, "nt", 1);
	}//if psfr[ievl]
    }//ievl
    dfree(xx);
    dcellfree(dmadd);
}

/**
   Shift gradient when new gradients are ready (in the end of parallel section
   in sim in CL or wfsgrad in OL). Do not execute in parallel with other
   routines. In GLAO mode, also averaged gradients from the same type of powfs.
*/
void shift_grad(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(parms->sim.glao){
	if(simu->gradlastcl){
	    dcellzero(simu->gradlastcl);
	}else{
	    simu->gradlastcl=dcellnew(parms->nwfsr, 1);
	}
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    const double scale=1./parms->powfs[ipowfs].nwfs;
	    for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
		int iwfs=parms->powfs[ipowfs].wfs[indwfs];
		dadd(&simu->gradlastcl->p[ipowfs], 1., simu->gradcl->p[iwfs], scale);
	    }
	}
    }else{
	dcellcp(&simu->gradlastcl, simu->gradcl); 
    }
    dcellcp(&simu->dmcmdlast, simu->dmcmd); 
    simu->reconisim = simu->isim;
}

/**
   Parse the input dead actuator location to actuator indices based on aloc.
*/
imat* act_coord2ind(loc_t *aloc,       /**<[in] Aloc*/
		    const char *fndead /**<[in] File containing dead actuators*/
		    ){
    dmat *dead=dread("%s", fndead);
    if(dead->ny!=2){
	error("%s must contain 2 columns of data\n", fndead);
    }
    loc_create_map_npad(aloc,1);
    long (*map)[aloc->map->nx]=(void*)aloc->map->p;
    double ox=aloc->map->ox;
    double oy=aloc->map->oy;
    double dx1=1./aloc->dx;
    PDMAT(dead, ps);
    imat *out=inew(dead->nx, 1);
    for(long jact=0; jact<dead->nx; jact++){
	long mapx=(long)round((ps[0][jact]-ox)*dx1);
	long mapy=(long)round((ps[1][jact]-oy)*dx1);
	long iact=map[mapy][mapx]-1;
	out->p[jact]=iact;
    }
    dfree(dead);
    return out;
}
