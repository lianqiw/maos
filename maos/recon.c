/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
/**
   \file recon.c
   Wavefront reconstruction and DM fitting routines
*/

typedef struct{
    int ic;
#if USE_PTHREAD > 0
    pthread_mutex_t ilock;
#endif
    int nc;
    dcell *xout;
    spcell *A;
    dcell *xin; 
    double alpha;
    int trans;
}EACH_T;
/**
   added temporarily to test speed
 */
static void spcellmulmat_each_do(EACH_T *info){
    int ic;
    while(LOCK(info->ilock),ic=info->ic++,UNLOCK(info->ilock),ic<info->nc){
	if(info->trans){
	    sptmulmat(&info->xout->p[ic], info->A->p[ic], info->xin->p[ic], info->alpha);
	}else{
	    spmulmat(&info->xout->p[ic], info->A->p[ic], info->xin->p[ic], info->alpha);
	}
    }
}
static void spcellmulmat_each(dcell **xout, spcell *A, dcell *xin, double alpha, int trans, int nthread){
    if(!*xout){
	*xout=dcellnew(xin->nx, xin->ny);
    }
    assert(xin->ny==1);
    EACH_T info;
    info.xout=*xout;
    info.A=A;
    info.xin=xin;
    info.alpha=alpha;
    info.trans=trans;
    info.ic=0;
    info.nc=info.xout->nx;
    PINIT(info.ilock);
    CALL(spcellmulmat_each_do, &info, nthread);
}

/**
   Apply tomography right hand operator without using assembled matrix. Fast and
   saves memory.  */
void TomoR(dcell **xout, const void *A, 
	   const dcell *gin, const double alpha){
    
    const RECON_T *recon=(const RECON_T *)A;
    dcell *g2=NULL;
    dcellcp(&g2, gin);
    TTFR(g2, recon->TTF, recon->PTTF);
    dcell *g3=NULL;
    spcellmulmat_thread(&g3,recon->saneai, g2, 1,recon->nthread);
    dcell *xx2=NULL;
    spcellmulmat_each(&xx2, recon->GG, g3, alpha, 1, recon->nthread);
    sptcellmulmat_thread(xout, recon->H0tomo, xx2, 1, recon->nthread);
    dcellfree(xx2);
    dcellfree(g2);
    dcellfree(g3);
}

/**
   Apply tomography left hand side operator without using assembled matrix. Fast
   and saves memory. Only useful in CG. Accumulates to xout; */
void TomoL(dcell **xout, const void *A, 
	   const dcell *xin, const double alpha){

    const RECON_T *recon=(const RECON_T *)A;
    dcell *gg=NULL;
    dcell *xx=NULL;
    spcellmulmat_thread(&xx, recon->H0tomo, xin, 1., recon->nthread);
    spcellmulmat_each(&gg, recon->GG, xx, 1., 0, recon->nthread);
    dcellfree(xx);
    TTFR(gg, recon->TTF, recon->PTTF);
    dcell *gg2=NULL;
    spcellmulmat_thread(&gg2, recon->saneai, gg,1,recon->nthread);
    dcellfree(gg);
    dcell *xx2=NULL;
    spcellmulmat_each(&xx2, recon->GG, gg2, alpha, 1, recon->nthread);
    sptcellmulmat_thread(xout, recon->H0tomo, xx2, 1, recon->nthread);
    dcellfree(xx2);

    dcellfree(gg2);
    if(recon->L2){
	apply_L2(xout, recon->L2, xin, alpha, recon->nthread);
    }else{
	apply_invpsd(xout, recon->RL.extra, xin, alpha);
    }
    if(recon->ZZT){//single point piston constraint
	sptcellmulmat_thread(xout, recon->ZZT, xin, alpha,
			     recon->nthread);
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
    spcellmulmat(&xp, recon->HX, xin, 1.);
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
   Backsubtitution with muv_chol_solve() is implemented.
*/
void tomo(dcell **opdr, const PARMS_T *parms, const RECON_T *recon, const dcell *grad){
    dcell *rhs=NULL;
    if(parms->tomo.assemble){
	muv(&rhs, &(recon->RR), grad, 1);
    }else{
	TomoR(&rhs, recon, grad, 1);
    }
    switch(parms->tomo.alg){
    case 0:
	muv_chol_solve_cell(opdr,&(recon->RL), rhs);
	break;
    case 1:{
	PREFUN pfun=NULL;
	const void *pdata=NULL;
	switch(parms->tomo.precond){
	case 0:
	    pfun=NULL;
	    pdata=NULL;
	    break;
	case 1:
	    pfun=fdpcg_precond;
	    pdata=(void*)recon;
	    break;
	default:
	    error("Invalid tomo.precond\n");
	}
	if(parms->tomo.assemble){
	    pcg(opdr, (CGFUN)muv, &(recon->RL), pfun, pdata,rhs, 1, parms->tomo.maxit);
	}else{
	    pcg(opdr, TomoL, recon, pfun,pdata,rhs, 1, parms->tomo.maxit);
	}
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
    muv(&rhs, &(recon->FR), opdr, 1);
    switch(parms->fit.alg){
    case 0:
	muv_chol_solve_cell(adm,&(recon->FL),rhs);
	break;
    case 1:
	pcg(adm, (CGFUN)muv, &(recon->FL),NULL,NULL, rhs, 1, parms->fit.maxit);
	break;
    default:
	error("Not implemented: alg=%d\n",parms->fit.alg);
    }
    dcellfree(rhs);
}
/**
   Update LGS focus or reference vector.

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
    //Updating LGS focus or reference vector.
    //fixme: what if strat!=1
    if(!simu->recon->RFlgs){
	warning("There is no LGS. No need to do focus tracking\n");
	return;
    }
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    dcell *graduse=dcellnew(parms->nwfs,1);
    PSPCELL(recon->GA,GA);
    PSPCELL(recon->G0focus,G0);
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
	    if(simu->gradpsol->p[iwfs]){
		graduse->p[iwfs]=ddup(simu->gradpsol->p[iwfs]);
	    }else{
		info("Forming PSOL grads for wfs %d\n",iwfs);
		graduse->p[iwfs]=ddup(simu->gradcl->p[iwfs]);
		if(simu->dmreal){
		    for(int idm=0; idm<parms->ndm; idm++){
			spmulmat(&graduse->p[iwfs],GA[idm][iwfs],
				 simu->dmreal->p[idm],1);
		    }
		}
	    }
	}else{//cl
	    graduse->p[iwfs]=ddup(simu->gradcl->p[iwfs]);
	}
	if(gs_x){
	    info("Subtracing tomo grad from wfs %d\n",iwfs);
	    for(int ips=0; ips<simu->recon->npsr; ips++){
		if(!G0[ips][iwfs]){
		    error("G0[%d][%d] is empty\n",ips,iwfs);
		}
		spmulmat(&graduse->p[iwfs],G0[ips][iwfs],simu->opdr->p[ips],-1);
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
   remove NGS modes from LGS DM commands
   if split_wt==1
   Rngs*GA*dmerr is zero
   if split_wt==2
   Doesn't perturb NGS modes in science direction.
*/
static inline void remove_dm_ngsmod(SIM_T *simu, dcell *dmerr){
    const RECON_T *recon=simu->recon;
    dcell *Mngs=NULL;
    dcellmm(&Mngs, recon->ngsmod->Pngs, dmerr, "nn",1);
    ngsmod2dm(&dmerr,recon, Mngs,-1);
    /*{
      info("NGSmod before removal\n");
      dshow(Mngs->p[0]);
      dcellzero(Mngs);
      dcellmm(&Mngs, recon->ngsmod->Pngs, dmerr, "nn",1);
      info("NGSmod after removal\n");
      dshow(Mngs->p[0]);
      }*/
    dcellfree(Mngs);
}
/**
   Removal tip/tilt on invidual DMs. Be careful about the roll off near the
edge.  */
static inline void remove_dm_tt(SIM_T *simu, dcell *dmerr){
    const RECON_T *recon=simu->recon;
    for(int idm=0; idm<simu->parms->ndm; idm++){
	dmat *utt=NULL;
	dmm(&utt, recon->ngsmod->Ptt->p[idm], dmerr->p[idm], "nn", -1);
	double *ptt;
	if(utt->nx==2){
	    ptt=alloca(3*sizeof(double));
	    ptt[0]=0; ptt[1]=utt->p[0]; ptt[2]=utt->p[1];
	}else{
	    ptt=utt->p;
	}
	loc_add_ptt(dmerr->p[idm]->p, ptt, recon->aloc[idm]);
	info("Adding P/T/T %g m %f %f mas to dm %d\n",
	     ptt[0],ptt[1]*206265000,ptt[2]*206265000,idm);
	dfree(utt);
    }
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
   cyclic shift the dmats.  */
static void shift_ring(int nap, dmat **ring, dmat *new){
    dmat *keep=ring[nap-1];
    for(int iap=nap-1; iap>=0; iap--){
	if(!iap){
	    ring[iap]=dref(new);
	}else{
	    ring[iap]=ring[iap-1];
	}
    }
    dfree(keep);
}

/**
   Apply fit right hand side matrix in CG mode without using assembled matrix
   for MOAO. subtract contributions from DMs that are in common path. Be careful
   which time step the dmcommon is. The DM common should use the commands on the
   step that you are going to apply the MOAO command for. That is the integrator
   output after this computation.  */

static void 
moao_FitR(dcell **xout, const RECON_T *recon, const PARMS_T *parms, int imoao, 
	  double thetax, double thetay, double hs, 
	  const dcell *opdr, const dcell *dmcommon, dcell **rhsout, const double alpha){
  
    //LOC_T *maloc=recon->moao[imoao].aloc;
    dcell *xp=dcellnew(1,1);
    xp->p[0]=dnew(recon->ploc->nloc,1);
    
    for(int ipsr=0; ipsr<recon->npsr; ipsr++){
	const double ht = parms->atmr.ht[ipsr];
	double scale=1.-ht/hs;
	prop_nongrid(recon->xloc[ipsr], opdr->p[ipsr]->p,
		     recon->ploc, NULL, xp->p[0]->p, 1, 
		     thetax*ht, thetay*ht, scale, 
		     0, 0);
    }
    for(int idm=0; idm<recon->ndm; idm++){
	const double ht = parms->dm[idm].ht;
	double scale=1.-ht/hs;
	prop_nongrid_cubic(recon->aloc[idm], dmcommon->p[idm]->p,
			   recon->ploc, NULL, xp->p[0]->p, -1, 
			   thetax*ht, thetay*ht, scale, 
			   parms->dm[idm].iac, 0, 0);
    }
    if(rhsout){
	*rhsout=dcelldup(xp);
    }
    double wt=1;
    applyW(xp, recon->moao[imoao].W0, recon->moao[imoao].W1, &wt);
    sptcellmulmat(xout, recon->moao[imoao].HA, xp, alpha);
    dcellfree(xp);
}
/**
   Apply fit left hand side matrix in CG mode
   without using assembled matrix. Slow. don't
   use. Assembled matridx is faster because of multiple
   directions.
*/

static void 
moao_FitL(dcell **xout, const void *A, 
	  const dcell *xin, const double alpha){
    const MOAO_T *moao=(const MOAO_T *)A;
    dcell *xp=NULL;
    double wt=1;
    spcellmulmat(&xp, moao->HA, xin, 1.);
    applyW(xp, moao->W0, moao->W1, &wt);
    sptcellmulmat(xout, moao->HA, xp, alpha);
    dcellfree(xp);xp=NULL;
    dcellmm(&xp, moao->NW, xin, "tn", 1);
    dcellmm(xout,moao->NW, xp, "nn", alpha);
    dcellfree(xp);
    spcellmulmat(xout, moao->actslave, xin, alpha);
}
/**
   mao_recon happens after the common DM fitting and its integrator output
   to take into account the delay in DM commands. there is no close loop
   filtering in MOAO DM commands, but there is still a time delay of 2
   cycles.
*/

void moao_recon(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    dcell *dmcommon=NULL;
    if(1){//Take OL result
	dcellcp(&dmcommon, simu->dmfit_hi);
    }else{//CL result
	if(parms->sim.closeloop){
	    if(parms->sim.fuseint){
		dcellcp(&dmcommon, simu->dmint[0]);
		if(parms->tomo.split){
		    remove_dm_ngsmod(simu, dmcommon);
		}
	    }else{
		dcellcp(&dmcommon, simu->dmint_hi[0]);
		//addlow2dm(&dmcommon, simu,simu->Mint_lo[0], 1);
	    }
	}else{
	    dcellcp(&dmcommon, simu->dmerr_hi);
	}
    }
    dcell *rhs=NULL;
    if(simu->moao_wfs){
	PDCELL(simu->moao_wfs, dmwfs);
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    dcell *dmmoao=NULL;
	    dcell *rhsout=NULL;
	    if(imoao>-1){
		double hs=parms->powfs[ipowfs].hs;
		dcellzero(rhs);
		moao_FitR(&rhs, recon, parms,  imoao, 
			  parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay, 
			  hs, simu->opdr, dmcommon, &rhsout, 1);
		pcg(&dmmoao, moao_FitL, &recon->moao[imoao], NULL, NULL, rhs, 
		    1, parms->fit.maxit);
		if(!isinf(parms->moao[imoao].stroke)){
		    int nclip=dclip(dmmoao->p[0],
				    -parms->moao[imoao].stroke,
				    parms->moao[imoao].stroke);
		    if(nclip>0){
			info("wfs %d: %d actuators clipped\n", iwfs, nclip);
		    }
		}
		shift_ring(simu->moao_wfs->nx, dmwfs[iwfs], dmmoao->p[0]);
		if(parms->plot.run){
		    drawopd("MOAO WFS RHS", recon->ploc, rhsout->p[0]->p, 
			    "MOAO for WFS","x (m)", "y(m)", "Wfs rhs %d", iwfs);
		    drawopd("MOAO WFS", recon->moao[imoao].aloc, dmmoao->p[0]->p,
			    "MOAO for WFS","x (m)", "y(m)", "Wfs %d", iwfs);
		}
		dcellfree(rhsout);
		dcellfree(dmmoao);
	    }//if imoao
	}//if wfs

	if(parms->save.run){
	    dcellwrite(simu->moao_wfs,"moao_%d/moao_wfs_%d",simu->seed, simu->isim);
	}
    }
    if(simu->moao_evl){
	PDCELL(simu->moao_evl, dmevl);
	int imoao=parms->evl.moao;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    dcell *dmmoao=dcellnew(1,1);
	    dmmoao->p[0]=ddup(dmevl[ievl][0]);//warm restart.
	    dcell *rhsout=NULL;
	    dcellzero(rhs);
	    moao_FitR(&rhs, recon, parms, imoao, 
		      parms->evl.thetax[ievl], parms->evl.thetay[ievl],
		      INFINITY, simu->opdr, dmcommon, &rhsout, 1);
	    
	    if(parms->save.run){
		dwrite(rhsout->p[0],"moao_%d/moao_rhs_evl%d_%d",simu->seed,ievl,simu->isim);
	    }
	    pcg(&dmmoao, moao_FitL, &recon->moao[imoao], NULL, NULL, rhs,
		1, parms->fit.maxit);
	    if(!isinf(parms->moao[imoao].stroke)){
		int nclip=dclip(dmmoao->p[0],
				-parms->moao[imoao].stroke,
				parms->moao[imoao].stroke);
		if(nclip>0){
		    info("evl %d: %d actuators clipped\n", ievl, nclip);
		}
	    }
	    shift_ring(simu->moao_evl->nx, dmevl[ievl], dmmoao->p[0]);
	    if(parms->plot.run){
		drawopd("MOAO EVL RHS", recon->ploc, rhsout->p[0]->p, 
			"MOAO for WFS","x (m)", "y(m)", "Evl %d", ievl);
		drawopd("MOAO EVL", recon->moao[imoao].aloc, dmevl[ievl][0]->p,
			"MOAO for EVL","x (m)", "y(m)", "Evl %d", ievl);
	    }
	    if(parms->save.dm){
		dwrite(rhsout->p[0], "moao_evlrhs_%d_%d", imoao, simu->isim);
		dwrite(dmmoao->p[0], "moao_evl_%d_%d", imoao, simu->isim);
	    }
	    dcellfree(dmmoao);
	    dcellfree(rhsout);
	}//ievl
	if(parms->save.run){
	    dcellwrite(simu->moao_evl,"moao_%d/moao_evl_%d",simu->seed,simu->isim);
	}
    }
    dcellfree(dmcommon);
    dcellfree(rhs);
}
/**
   Calls tomo() and fit() to do the tomography and DM fit. Do error signal and
   split tomography.  */
void tomofit(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(!parms->sim.closeloop || parms->dbg.fitonly || 
       simu->dtrat_hi==1 || (simu->isim+1)%simu->dtrat_hi==0){
	if(parms->dbg.fitonly){
	    dcellfree(simu->opdr);
	    simu->opdr=atm2xloc(simu);
	}else{
	    tomo(&simu->opdr,parms,recon,simu->gradpsol);
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
	if(parms->save.opdr){
	    dcellwrite(simu->opdr,"opdr_%d/opdr_%d.bin", simu->seed,simu->isim);
	}
	if(parms->save.opdx){
	    dcell *opdx=atm2xloc(simu);
	    dcellwrite(opdx,"opdx_%d/opdx_%d.bin.gz",simu->seed,simu->isim);
	    dcellfree(opdx);
	}
	if(parms->save.dm){
	    dcellwrite(simu->dmfit_hi,"dmfit_hi_%d/dmfit_hi_%d.bin.gz",
		       simu->seed,simu->isim);
	}
	//Ploting.
	if(parms->plot.run){
	    for(int i=0; simu->opdr && i<simu->opdr->nx; i++){
		drawopd("Recon", recon->xloc[i], simu->opdr->p[i]->p, 
			"Reconstructed Atmosphere","x (m)","y (m)","opdr %d",i);
	    }
	    for(int i=0; simu->dmfit_hi && i<simu->dmfit_hi->nx; i++){
		drawopd("DM", recon->aloc[i], simu->dmfit_hi->p[i]->p,
			"DM Fitting Output","x (m)", "y (m)","Fit %d",i);
	    }
	    for(int idm=0; simu->dmreal && idm<simu->parms->ndm; idm++){
		drawopd("DM", simu->recon->aloc[idm], simu->dmreal->p[idm]->p,
			"Actual DM Actuator Commands","x (m)", "y (m)",
			"Real %d",idm);
	    }
	}
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

	if(parms->save.dm){
	    dcellwrite(simu->dmerr_hi,"dmerr_hi_%d/dmerr_hi_%d_raw.bin.gz",
		       simu->seed,simu->isim);
	    if(parms->sim.fuseint){
		dcellwrite(simu->dmint[0],"dmint_%d/dmint_%d.bin",
			   simu->seed,simu->isim);
	    }else{
		dcellwrite(simu->dmint_hi[0],"dmint_hi_%d/dmint_hi_%d.bin",
			   simu->seed,simu->isim);   
	    }
	}

	if(!parms->dbg.fitonly && parms->tomo.split==1){//ahst
	    remove_dm_ngsmod(simu, simu->dmerr_hi);
	}
	if(parms->tomo.split_rtt && parms->tomo.split){
	    remove_dm_tt(simu, simu->dmerr_hi);
	}
	
	if(parms->save.dm){
	    dcellwrite(simu->dmerr_hi,"dmerr_hi_%d/dmerr_hi_%d_net.bin.gz",
		       simu->seed,simu->isim);
	}
	if(parms->plot.run){
	    for(int idm=0; idm<parms->ndm; idm++){
		drawopd("DM",recon->aloc[idm], simu->dmerr_hi->p[idm]->p,
			"DM Error Signal (Hi)","x (m)","y (m)",
			"Err Hi %d",idm);
	    }
	}
    }//if high order has output

    if(!parms->dbg.fitonly && parms->tomo.split){
	if(parms->tomo.split==2){
	    info("accumulating opdrmvst\n");
	    dcelladd(&simu->opdrmvst, 1, simu->opdr, 1./simu->dtrat_lo);
	}
	//Low order
	if(!parms->sim.closeloop || simu->dtrat_lo==1 || (simu->isim+1)%simu->dtrat_lo==0){
	    dcellzero(simu->Merr_lo);
	    switch(parms->tomo.split){
	    case 1:{
		NGSMOD_T *ngsmod=recon->ngsmod;
		if(!parms->tomo.split_idealngs){//Low order NGS recon.
		    dcellmm(&simu->Merr_lo,ngsmod->Rngs,simu->gradcl,"nn",1);
		}//else: there is ideal NGS correction done in perfevl.
	    }
		break;
	    case 2:{
		dcellmm(&simu->gradpsol, recon->G0L, simu->opdrmvst, "nn",-1);
		dcellmm(&simu->Merr_lo, recon->MVRngs, simu->gradpsol, "nn",1);
		if(parms->sim.fuseint){
		    dcelladd(&simu->Merr_lo, 1., simu->dmint[1], -1);
		    error("This mode is not finished\n");
		}else{//form error signal
		    dcelladd(&simu->Merr_lo, 1., simu->Mint_lo[1], -1);
		}
		dcellzero(simu->opdrmvst);info("zero opdrmvst\n");
	    }
		break;
	    default:
		error("Invalid parms->tomo.split: %d",parms->tomo.split);
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
		for(int idm=0; idm<parms->ndm; idm++){
		    drawopd("DM",recon->aloc[idm], dmlo->p[idm]->p,
			    "DM Error Signal (Lo)","x (m)","y (m)",
			    "Err Lo %d",idm);
		}
		dcellfree(dmlo);
	    }
	}else{
	    dcellfree(simu->Merr_lo);//don't have output.
	}

	if(parms->save.dm){
	    dcellwrite(simu->Merr_lo, "Merr_lo_%d/Merr_lo_%d", simu->seed,simu->isim);
	    if(simu->Mint_lo){
		dcellwrite(simu->Mint_lo[0], "Mint_lo_%d/Mint_lo_%d", simu->seed,simu->isim);
	    }
	}
    }
    if(parms->sim.mffocus){
	focus_tracking(simu);
    }
}
