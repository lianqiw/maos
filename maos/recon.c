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

#include "maos.h"
#include "recon.h"
#include "sim_utils.h"
#include "fdpcg.h"
#include "sim.h"
#include "recon_utils.h"
#include "ahst.h"
#include "moao.h"

/**
   \file recon.c
   Wavefront reconstruction and DM fitting routines
*/

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
	int ipowfs = parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].skip) continue;
	dmat *xx=dnew(recon->ploc->nloc, 1);
	const double hs=parms->powfs[ipowfs].hs;
	for(int ips=0; ips<nps; ips++){
	    if(parms->tomo.square && !parms->dbg.tomo_hxw){
		//Do the ray tracing instead of using HXW.
		double ht=recon->ht->p[ips];
		double displace[2];
		displace[0]=parms->wfs[iwfs].thetax*ht;
		displace[1]=parms->wfs[iwfs].thetay*ht;
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
		int ipowfs = parms->wfs[iwfs].powfs;
		const double hs=parms->powfs[ipowfs].hs;
		double displace[2];
		displace[0]=parms->wfs[iwfs].thetax*ht;
		displace[1]=parms->wfs[iwfs].thetay*ht;
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
		apply_invpsd(&data->xout, recon->invpsd, data->xin, data->alpha, ips);
		break;
	    case 2:
		apply_fractal(&data->xout, recon->fractal, data->xin, data->alpha, ips);
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
#define USE_NEW 0
#if USE_NEW
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
#else
    dcell *g3=NULL;
    spcellmulmat_thread(&g3,recon->saneai, gg, 1,recon->nthread);
    dcellfree(gg);
    spcellmulmat_each(&gg, recon->GP, g3, alpha, 1, recon->nthread);
    sptcellmulmat_thread(xout, recon->HXWtomo, gg, 1, recon->nthread);
    dcellfree(g3);
#endif
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
    dcell *gg=dcellnew(parms->nwfs, 1);
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
    //Remove global Tip/Tilt, differential focus. Need all wfs gradients.
    TTFR(gg, recon->TTF, recon->PTTF);
    CALL_THREAD(info_iwfs2, recon->nthread, 1);
    CALL_THREAD(info_ips, recon->nthread, 1);
   
    /*      square=1  square=0 (1 thread on T410s)
	    iwfs1: takes 6 ms 13 ms
	    iwfs2: takes 4 ms 4 ms
	    ips:   takes 6 ms 9 ms
    */
    dcellfree(gg);
    
    if(recon->ZZT){//single point piston constraint. fast if any
	sptcellmulmat_thread(xout, recon->ZZT, xin, alpha, recon->nthread);
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
	spcellmulmat(&xp, recon->HXF, xin, 1.);
    }else{//Do the ray tracing directly.
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
    sptcellmulmat(xout, recon->HA, xp, alpha);
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
    spcellmulmat(&xp, recon->HA, xin, 1.);
    applyW(xp, recon->W0, recon->W1, recon->fitwt->p);
    sptcellmulmat(xout, recon->HA, xp, alpha);
    dcellfree(xp);xp=NULL;
    dcellmm(&xp,recon->NW, xin, "tn", 1);
    dcellmm(xout,recon->NW, xp, "nn", alpha);
    dcellfree(xp);
}

/**
   Carry out tomography. (Precondition) CG with pcg() and Cholesky
   Backsubtitution with muv_direct_solve() is implemented.
*/
void tomo(dcell **opdr, const PARMS_T *parms, const RECON_T *recon, const dcell *grad, int maxit){
    dcell *rhs=NULL;
    if(parms->tomo.assemble){
	muv(&rhs, &(recon->RR), grad, 1);
    }else{
	TomoR(&rhs, recon, grad, 1);
    }
    switch(parms->tomo.alg){
    case 0:
	muv_direct_solve_cell(opdr,&(recon->RL), rhs);
	break;
    case 1:{
	PREFUN pfun;
	const void *pdata;
	CGFUN cgfun;
	const void *cgdata;
	switch(parms->tomo.precond){
	case 0: /*This is no preconditioner*/
	    pfun=NULL;
	    pdata=NULL;
	    break;
	case 1: /*Fourier Domain preconditioner*/
	    pfun=fdpcg_precond;
	    pdata=(void*)recon;
	    break;
	default:
	    pfun=NULL;
	    pdata=NULL;
	    error("Invalid tomo.precond\n");
	}

	if(parms->tomo.assemble){
	    cgfun=(CGFUN)muv;
	    cgdata=&(recon->RL);
	}else{
	    cgfun=TomoL;
	    cgdata=recon;
	}
	pcg(opdr, cgfun, cgdata, pfun, pdata, rhs, 
	    recon->warm_restart, maxit);
    }
	break;
    case 2:
	muv_direct_solve_cell(opdr,&(recon->RL), rhs);
	break;
    case 3: {//Block Gaussian Seidel
	
    }
	break;
    default:
	error("Not implemented: alg=%d\n",parms->tomo.alg);
    }
    dcellfree(rhs);
}
/**
   Carry out DM fit. Un-Precondition CG and Cholesky Backsubtitution is
   implemented.  */
void fit(dcell **adm, const PARMS_T *parms, 
	 const RECON_T *recon, const dcell *opdr){
    if(parms->ndm==0) return;
    dcell *rhs=NULL;
    if(recon->FR.M){
	muv(&rhs, &(recon->FR), opdr, 1);
    }else{
	FitR(&rhs, recon, opdr, 1);
    }
    switch(parms->fit.alg){
    case 0:
	muv_direct_solve_cell(adm,&(recon->FL),rhs);
	break;
    case 1:{
	PREFUN pfun;
	const void *pdata;
	CGFUN cgfun;
	const void *cgdata;
	switch(parms->fit.precond){
	case 0:
	    pfun=NULL;
	    pdata=NULL;
	    break;
	default:
	    pfun=NULL;
	    pdata=NULL;
	    error("Invalid fit.precond\n");
	}
	cgfun=(CGFUN)muv;
	cgdata=&(recon->FL);
	pcg(adm, cgfun, cgdata, pfun, pdata, rhs, 
	    recon->warm_restart, parms->fit.maxit);
	break;
    }
    case 2:
	muv_direct_solve_cell(adm,&(recon->FL),rhs);
	break;
    default:
	error("Not implemented: alg=%d\n",parms->fit.alg);
    }
    dcellfree(rhs);
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
static void windest(SIM_T *simu){
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
   Calls tomo() and fit() to do the tomography and DM fit. Do error signal and
   split tomography.  */
void tomofit(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    //2010-12-16: replaced isim+1 by isim since recon is delayed from wfsgrad by 1 frame.
    if(!parms->sim.closeloop || parms->sim.fitonly || simu->dtrat_hi==1 || (simu->isim)%simu->dtrat_hi==0){
	if(parms->sim.fitonly){
	    dcellfree(simu->opdr);
	    //simu->opdr=atm2xloc(simu);
	}else{
	    int maxit=parms->tomo.maxit;
	    if(parms->dbg.ntomo_maxit){
		if(simu->isim<parms->dbg.ntomo_maxit){
		    maxit=parms->dbg.tomo_maxit[simu->isim];
		    info2("Running tomo.maxit=%d\n",maxit);
		}else{
		    error("Out of range\n");
		}
	    }
	    //computes simu->opdr
	    tomo(&simu->opdr,parms,recon,simu->gradlastol,maxit);
	}
	if(parms->tomo.windest){
	    info2("Estimating wind direction and speed using FFT method\n");
	    windest(simu);
	    //Update wind, and interpolation matrix.
	}
	if(parms->tomo.windshift){
	    int factor=parms->tomo.windshift;
	    if(!simu->windshift){
		simu->windshift=spcellnew(recon->npsr, 1);
		for(int ips=0; ips<recon->npsr; ips++){
		    double dispx=simu->dt*simu->atm[ips]->vx*factor;//2 is two cycle delay.
		    double dispy=simu->dt*simu->atm[ips]->vy*factor;
		    info("ips=%d: dispx=%g, dispy=%g\n", ips, dispx, dispy);
		    simu->windshift->p[ips]=mkhb(recon->xloc[ips], recon->xloc[ips], NULL,
						 dispx,dispy,1,0,0);
		}
		spcellwrite(simu->windshift,"windshift");
	    }
	    info2("Using wind information to shift opdr by %d v*dt.\n", factor);
	    for(int ips=0; ips<recon->npsr; ips++){
		dmat *tmp=simu->opdr->p[ips];
		simu->opdr->p[ips]=NULL;
		spmulmat(&simu->opdr->p[ips], simu->windshift->p[ips], tmp, 1);
		dfree(tmp);
	    }
	}

	fit(&simu->dmfit_hi,parms,recon,simu->opdr);

	dcellcp(&simu->dmerr_hi, simu->dmfit_hi);//keep dmfit_hi for warm restart
    
	/*
	  Forming LGS error signal.
	  2010-01-07: changed dmreal_hi to dmreal to comply with
	  the block diagram. This is before NGS mode removal
	  keep dmfit for warm restart. 
	*/
 
	if(parms->sim.fuseint){
	    if(parms->dbg.psol){
		warning("Using dm for next step to form err signal\n");
		dcelladd(&simu->dmerr_hi, 1., simu->dmint[0], -1);
	    }else{
		dcelladd(&simu->dmerr_hi, 1., simu->dmint[1], -1);
	    }
	}else{
	    /**
	       2010-07-23: Moved remove_dm_ngsmod to after forming error signal. 
	    */
	    if(parms->dbg.psol){
		warning("Using dm for next step to form err signal\n");
		dcelladd(&simu->dmerr_hi, 1., simu->dmint_hi[0], -1);
	    }else{
		dcelladd(&simu->dmerr_hi, 1., simu->dmint_hi[1], -1);
	    }
	}

	if(!parms->sim.fitonly && parms->tomo.split==1){//ahst
	    remove_dm_ngsmod(simu, simu->dmerr_hi);
	}
	if(parms->tomo.ahst_rtt && parms->tomo.split){
	    remove_dm_tt(simu, simu->dmerr_hi);
	}

    }//if high order has output

    if(!parms->sim.fitonly && parms->tomo.split){
	if(parms->tomo.split==2){
	    dcelladd(&simu->opdrmvst, 1, simu->opdr, 1./simu->dtrat_lo);
	}
	//Low order has output
	//2010-12-16: replaces isim+1 by isim.
	if(!parms->sim.closeloop || simu->dtrat_lo==1 || (simu->isim)%simu->dtrat_lo==0){
	    dcellzero(simu->Merr_lo);
	    switch(parms->tomo.split){
	    case 1:{
		NGSMOD_T *ngsmod=recon->ngsmod;
		if(!parms->tomo.ahst_idealngs){//Low order NGS recon.
		    dcellmm(&simu->Merr_lo,ngsmod->Rngs,simu->gradlastcl,"nn",1);
		}//else: there is ideal NGS correction done in perfevl.
	    }
		break;
	    case 2:{
		dcellmm(&simu->gradlastol, recon->GXL, simu->opdrmvst, "nn",-1);
		dcellmm(&simu->Merr_lo, recon->MVRngs, simu->gradlastol, "nn",1);
		if(parms->sim.fuseint){
		    dcelladd(&simu->Merr_lo, 1., simu->dmint[1], -1);
		    error("This mode is not finished\n");
		}else{//form error signal
		    dcelladd(&simu->Merr_lo, 1., simu->Mint_lo[1], -1);
		}
		dcellzero(simu->opdrmvst);
	    }
		break;
	    default:
		error("Invalid parms->tomo.split: %d",parms->tomo.split);
	    }
	}else{
	    dcellfree(simu->Merr_lo);//don't have output.
	}
    }
   
}
/**
   least square reconstructor
*/
void lsr(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int do_hi=(!parms->sim.closeloop || parms->sim.fitonly || 
		     simu->dtrat_hi==1 || (simu->isim)%simu->dtrat_hi==0);
    const int do_low=parms->tomo.split && (!parms->sim.closeloop || simu->dtrat_lo==1
					   || (simu->isim)%simu->dtrat_lo==0);
    dcell *graduse=NULL;
    if(do_hi || do_low){
	if(parms->sim.recon==2){
	    graduse=dcellnew(parms->npowfs, 1);
	    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		double scale=1./parms->powfs[ipowfs].nwfs;
		for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs[indwfs];
		    dadd(&graduse->p[ipowfs], 1, simu->gradlastcl->p[iwfs], scale);
		}
	    }
	}else{
	    graduse=simu->gradlastcl;
	}
    }
    //2010-12-16: replaced isim+1 by isim since recon is delayed from wfsgrad by 1 frame.
    if(do_hi){
	dcell *rhs=NULL;
	muv(&rhs, &(recon->LR), graduse, 1);
	switch(parms->tomo.alg){
	case 0://CBS
	    muv_direct_solve_cell(&simu->dmerr_hi, &recon->LL, rhs);
	    break;
	case 1://CG
	    pcg(&simu->dmerr_hi, (CGFUN)muv, &recon->LL, NULL, NULL, rhs,
		recon->warm_restart, parms->tomo.maxit);
	    break;
	case 2://SVD inversion
	    muv_direct_solve_cell(&simu->dmerr_hi, &recon->LL, rhs);
	    break;
	default:
	    error("Not implemented\n");
	}
	dcellfree(rhs);
	if(!parms->sim.fitonly && parms->tomo.split==1){//ahst
	    remove_dm_ngsmod(simu, simu->dmerr_hi);
	}
    }//if high order has output
    if(parms->tomo.split){
	//Low order has output
	//2010-12-16: replaced isim+1 by isim since recon is delayed from wfsgrad by 1 frame.
	if(do_low){
	    dcellzero(simu->Merr_lo);
	    NGSMOD_T *ngsmod=recon->ngsmod;
	    dcellmm(&simu->Merr_lo,ngsmod->Rngs,graduse,"nn",1);
	}else{
	    dcellfree(simu->Merr_lo);//don't have output.
	}
    } 
    if(graduse != simu->gradlastcl){
	dcellfree(graduse);
    }
}
/**
   Deformable mirror control. call tomofit() to do tomo()/fit() or lsr() to do
   least square reconstruction. */
void reconstruct(SIM_T *simu){
    double tk_start=myclockd();
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(simu->gradlastol || simu->gradlastcl){
	switch(parms->sim.recon){//mv
	case 0:
	    tomofit(simu);//tomography and fitting.
	    break;
	case 1:
	case 2:
	    lsr(simu);
	    break;
	}
	if(recon->moao){
	    moao_recon(simu);
	}
	if(parms->sim.mffocus){
	    focus_tracking(simu);
	}
    }
    if(parms->plot.run){
	if(parms->sim.recon==0){
	    for(int i=0; simu->opdr && i<simu->opdr->nx; i++){
		drawopd("Recon", recon->xloc[i], simu->opdr->p[i]->p, NULL,
			"Reconstructed Atmosphere","x (m)","y (m)","opdr %d",i);
	    }
	    for(int i=0; simu->dmfit_hi && i<simu->dmfit_hi->nx; i++){
		drawopd("DM", recon->aloc[i], simu->dmfit_hi->p[i]->p,NULL,
			"DM Fitting Output","x (m)", "y (m)","Fit %d",i);
	    }
	}
	for(int idm=0; simu->dmerr_hi && idm<parms->ndm; idm++){
	    drawopd("DM",recon->aloc[idm], simu->dmerr_hi->p[idm]->p,NULL,
		    "DM Error Signal (Hi)","x (m)","y (m)",
		    "Err Hi %d",idm);
	}
    }
    if(parms->plot.run && simu->Merr_lo){
	dcell *dmlo=NULL;
	switch(simu->parms->tomo.split){
	case 1:
	    ngsmod2dm(&dmlo, recon, simu->Merr_lo, 1);
	    break;
	case 2:
	    dcellmm(&dmlo, simu->recon->MVModes, simu->Merr_lo, "nn", 1);
	    break;
	}
	for(int idm=0; dmlo && idm<parms->ndm; idm++){
	    drawopd("DM",recon->aloc[idm], dmlo->p[idm]->p,NULL,
		    "DM Error Signal (Lo)","x (m)","y (m)",
		    "Err Lo %d",idm);
	}
	dcellfree(dmlo);
    }
    if(parms->sim.recon==0){
	if(parms->save.opdr){
	    cellarr_dcell(simu->save->opdr, simu->opdr);
	}
	if(parms->save.dm){
	    cellarr_dcell(simu->save->dmfit_hi, simu->dmfit_hi);
	}
	if(parms->save.opdx || parms->plot.opdx){
	    dcell *opdx;
	    if(parms->sim.fitonly){
		opdx=simu->opdr;
	    }else{
		opdx=atm2xloc(simu);
	    }
	    if(parms->save.opdx){
		cellarr_dcell(simu->save->opdx, opdx);
	    }
	    if(parms->plot.opdx){ //draw opdx
		for(int i=0; i<opdx->nx; i++){
		    drawopd("Recon", recon->xloc[i], opdx->p[i]->p, NULL,
			    "Atmosphere Projected to XLOC","x (m)","y (m)","opdx %d",i);
		}
	    }
	    if(!parms->sim.fitonly){
		dcellfree(opdx);
	    }
	}
    }
    if(parms->save.dm){
	cellarr_dcell(simu->save->dmerr_hi, simu->dmerr_hi);
	if(parms->sim.fuseint){
	    cellarr_dcell(simu->save->dmint, simu->dmint[0]);
	}else{
	    cellarr_dcell(simu->save->dmint_hi, simu->dmint_hi[0]);
	}
	if(simu->save->Merr_lo){
	    cellarr_dcell(simu->save->Merr_lo, simu->Merr_lo);
	}
	if(simu->Mint_lo){
	    cellarr_dcell(simu->save->Mint_lo, simu->Mint_lo[0]);
	}
    }
    simu->tk_recon=myclockd()-tk_start;
}
