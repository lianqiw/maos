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
#include "common.h"
#include "sim.h"
#include "sim_utils.h"
#include "ahst.h"
#include "mtch.h"
#include "save.h"
#include "setup_recon.h"
#include "recon_utils.h"
#include "setup_powfs.h"
#include "pywfs.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file wfsgrad.c
   contains functions that computes WFS gradients in geometric or physical optics mode.
*/
#define TIMING 0
#if TIMING == 1
#define TIM(A) double tk##A=myclockd()
#else
#define TIM(A)
#endif

/**
   Propagate atm onto WFS subaperture grid, and then to fine lenslet grid.
 */
void wfs_ideal_atm(SIM_T *simu, dmat *opd, int iwfs, double alpha){
    const PARMS_T *parms=simu->parms;
    POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int jwfs=parms->powfs[ipowfs].wfsind->p[iwfs];
    const double hs=parms->wfs[iwfs].hs;
    const double hc=parms->powfs[ipowfs].hc;
    if(parms->sim.wfsalias==2 || parms->sim.idealwfs==2){
	loc_t *aloc=powfs[ipowfs].fit[jwfs].aloc->p[0];
	dcell *wfsopd=dcellnew(1,1); wfsopd->p[0]=dnew(aloc->nloc, 1);
	FIT_T *fit=&powfs[ipowfs].fit[jwfs];
	muv_solve(&wfsopd, &fit->FL, &fit->FR, 0);
	prop_nongrid(aloc, wfsopd->p[0]->p, powfs[ipowfs].loc, opd->p, alpha, 0, 0, 1,  0, 0);
	dcellfree(wfsopd);
    }else{
	const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
	for(int idm=0; idm<parms->ndm; idm++){
	    loc_t *loc=powfs[ipowfs].loc_dm?powfs[ipowfs].loc_dm->p[wfsind+idm*parms->nwfs]:powfs[ipowfs].loc;
	    const double ht = parms->dm[idm].ht+parms->dm[idm].vmisreg-hc;
	    double dispx=ht*parms->wfs[iwfs].thetax;
	    double dispy=ht*parms->wfs[iwfs].thetay;
	    double scale=1.-ht/hs;
	    if(scale<0) continue;
	    prop_grid(simu->dmprojsq->p[idm], loc, opd->p, 
		      alpha, dispx, dispy, scale, 0, 0, 0);
	}
    }
}

/**
   computes close loop and pseudo open loop gradidents for both gometric and
   physical optics WFS. Calls wfsints() to accumulate WFS subapertures images in
   physical optics mode.  */

void wfsgrad_iwfs(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    const int isim=simu->isim;
    const int iwfs=info->start;
    const PARMS_T *parms=simu->parms;
    const int ipowfs=parms->wfs[iwfs].powfs;
    //if(isim<parms->powfs[ipowfs].step) return;
    assert(iwfs<parms->nwfs);
    /*
      simu->gradcl is CL grad output (also for warm-restart of maxapriori
      simu->gradacc is internal, to accumulate geometric grads.
      do not accumulate opd. accumate ints for phy, g for GS
    */
    /*input */
    
    mapcell *atm=simu->atm;
    const RECON_T *recon=simu->recon;
    const POWFS_T *powfs=simu->powfs;
    /*output */
    const int CL=parms->sim.closeloop;
    const int nps=parms->atm.nps;
    const double atmscale=simu->atmscale?simu->atmscale->p[isim]:1;
    const double dt=parms->sim.dt;
    TIM(0);
    /*The following are truly constants for this powfs */
    const int imoao=parms->powfs[ipowfs].moao;
    const int nsa=powfs[ipowfs].saloc->nloc;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    const double hs=parms->wfs[iwfs].hs;
    const int dtrat=parms->powfs[ipowfs].dtrat;
    const int save_gradgeom=parms->save.gradgeom->p[iwfs];
    const int save_opd =parms->save.wfsopd->p[iwfs];
    const int save_ints=parms->save.ints->p[iwfs];
    const int noisy=parms->powfs[ipowfs].noisy;
    /*The following depends on isim */
    /*const int dtrat_reset=(isim%dtrat==0); */
    const int dtrat_output=(isim+1)%dtrat==0;
    const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
    const int do_pistatout=parms->powfs[ipowfs].pistatout&&isim>=parms->powfs[ipowfs].pistatstart;
    const int do_geom=(!do_phy || save_gradgeom || do_pistatout) && parms->powfs[ipowfs].type==0;
    const double *realamp=powfs[ipowfs].realamp?powfs[ipowfs].realamp->p[wfsind]->p:0;
    dmat *gradcalc=NULL;
    dmat **gradacc=&simu->gradacc->p[iwfs];
    dmat **gradout=&simu->gradcl->p[iwfs];
    dcell *ints=simu->ints->p[iwfs];
    dmat  *opd=simu->wfsopd->p[iwfs];
    dzero(opd);
    if(isim%dtrat==0){
	dcellzero(ints);
	dzero(*gradacc);
    }
    /* Now begin ray tracing. */
    if(atm && ((!parms->sim.idealwfs && !parms->powfs[ipowfs].lo)
	       || (!parms->sim.wfsalias && parms->powfs[ipowfs].lo))){
	for(int ips=0; ips<nps; ips++){
	    thread_t *wfs_prop=simu->wfs_prop_atm[iwfs+parms->nwfs*ips];
	    PROPDATA_T *wfs_propdata=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
	    wfs_propdata->phiout=opd->p;
	    wfs_propdata->displacex1=-atm->p[ips]->vx*dt*isim;
	    wfs_propdata->displacey1=-atm->p[ips]->vy*dt*isim;
	    wfs_propdata->alpha=atmscale;
	    /* have to wait to finish before another phase screen. */
	    CALL_THREAD(wfs_prop, 0);
	}
    }
    /* 
       Propagate controllable component of atm (within range of DM) to wfs.
       wfsalias: atm - controllable.
       idealwfs: just controllable.
    */
    /* timing: most expensive 0.10 per LGS for*/
    if(!parms->powfs[ipowfs].lo && (parms->sim.wfsalias || parms->sim.idealwfs)){
	double alpha=parms->sim.idealwfs?1:-1;
	wfs_ideal_atm(simu, opd, iwfs, alpha);
    }


    if(simu->telws){/*Wind shake */
	double tmp=simu->telws->p[isim];
	double angle=simu->winddir?simu->winddir->p[0]:0;
	double ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
	loc_add_ptt(opd->p, ptt, powfs[ipowfs].loc);
    }
    
    double focus=wfsfocusadj(simu, iwfs);
    if(fabs(focus)>1e-20){
	loc_add_focus(opd->p, powfs[ipowfs].loc, focus);
    }
    
    /* Add surface error*/
    if(powfs[ipowfs].opdadd && powfs[ipowfs].opdadd->p[wfsind]){
	dadd(&opd,1, powfs[ipowfs].opdadd->p[wfsind],1);
    }

    if(save_opd){
	zfarr_push(simu->save->wfsopdol[iwfs], isim, opd);
    }
    TIM(1);
    if(CL){
	for(int idm=0; idm<parms->ndm; idm++){
	    thread_t *wfs_prop=simu->wfs_prop_dm[iwfs+parms->nwfs*idm];
	    PROPDATA_T *wfs_propdata=&simu->wfs_propdata_dm[iwfs+parms->nwfs*idm];
	    wfs_propdata->phiout=opd->p;
	    CALL_THREAD(wfs_prop, 0);
	}/*idm */
	double ptt[3]={0,0,0};
	if(simu->ttmreal){
	    ptt[1]-=simu->ttmreal->p[0];
	    ptt[2]-=simu->ttmreal->p[1];
	}
	//For dithering with downlink instead of uplink FSM
	if(simu->fsmreal && simu->fsmreal->p[iwfs] && !powfs[ipowfs].llt){
	    ptt[1]-=simu->fsmreal->p[iwfs]->p[0];
	    ptt[2]-=simu->fsmreal->p[iwfs]->p[1];
	}
	if(ptt[1] || ptt[2]){
	    loc_add_ptt(opd->p, ptt, powfs[ipowfs].loc);
	}
    }
    if(parms->powfs[ipowfs].skip && parms->tomo.ahst_idealngs==1){
	//apply ideal NGS modes to NGS WFS
	ngsmod2science(opd, powfs[ipowfs].loc, recon->ngsmod, 
		       parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay,
		       simu->cleNGSm->p+isim*recon->ngsmod->nmod, -1);
    }
    if(imoao>-1){
	dmat **dmwfs=simu->dm_wfs->p;
	if(dmwfs[iwfs]){
	    /* No need to do mis registration here since the MOAO DM is attached
	       to close to the WFS.*/
	    prop_nongrid_pts(recon->moao[imoao].aloc->p[0], dmwfs[iwfs]->p,
			     powfs[ipowfs].pts, opd->p, -1, 0, 0, 1, 0, 0);
	}
    }
 
    if(parms->powfs[ipowfs].fieldstop>0 && parms->powfs[ipowfs].type==0){
	locfft_fieldstop(powfs[ipowfs].fieldstop, opd, parms->powfs[ipowfs].wvlwts);
    }

    if(save_opd){
	zfarr_push(simu->save->wfsopd[iwfs], isim, opd);
    }
    if(parms->plot.run){
	drawopdamp("wfsopd",powfs[ipowfs].loc,opd->p,realamp,parms->dbg.draw_opdmax->p,
		   "WFS OPD","x (m)", "y (m)", "WFS %d", iwfs);
    }
    if(do_geom){
	/* Now Geometric Optics gradient calculations. if dtrat==1, we compute
	   gradients directly to gradacc, which is the same as gradcalc. If
	   dtrat>1, we compute gradients to gradcalc, and accumulate to
	   gradacc. gradcalc is used to shift pistat. We DONOT include gradoff
	   adjustment to gradref, but only do it on gradcl. This will make the
	   pistat always peak in center no matter what NCPA is present.
	*/
	if(!do_pistatout || parms->powfs[ipowfs].pistatstc || dtrat==1){
	    //we do not need separate gradcalc.
	    gradcalc=dref(*gradacc);
	}//else: calculate first to gradcalc then add to gradacc
	if(parms->powfs[ipowfs].gtype_sim==1){ /*compute ztilt. */
	    pts_ztilt(&gradcalc,powfs[ipowfs].pts,
		      powfs[ipowfs].saimcc->p[powfs[ipowfs].nsaimcc>1?wfsind:0], 
		      realamp, opd->p);
	}else{/*G tilt */
	    dspmm(&gradcalc,PR(powfs[ipowfs].GS0,wfsind,0),opd,"nn",1);
	}
	if(gradcalc->p!=(*gradacc)->p){
	    dadd(gradacc, 1, gradcalc, 1);
	}
    }

    ccell *psfout=NULL;
    zfarr *psfoutzfarr=NULL;
    zfarr *ztiltoutzfarr=NULL;
    if(parms->powfs[ipowfs].psfout){
	psfout=simu->wfspsfout->p[iwfs];
	psfoutzfarr=simu->save->wfspsfout[iwfs];
	ztiltoutzfarr=simu->save->ztiltout[iwfs];
    }
    TIM(2);
    /* Now begin Physical Optics Intensity calculations */
    if(do_phy || psfout || do_pistatout || parms->powfs[ipowfs].dither==1){
	dmat *lltopd=NULL;
	if(powfs[ipowfs].llt){//If there is LLT, apply FSM onto LLT
	    if(powfs[ipowfs].llt->ncpa){
		lltopd=ddup(PR(powfs[ipowfs].llt->ncpa, wfsind, 0));
	    }else{
		lltopd=dnew(powfs[ipowfs].llt->pts->nx, powfs[ipowfs].llt->pts->nx);
	    }
	    const long illt=parms->powfs[ipowfs].llt->i->p[wfsind];
	    if(atm){/*LLT OPD */
		for(int ips=0; ips<nps; ips++){
		    const double hl=atm->p[ips]->h;
		    const double scale=1.-hl/hs;
		    if(scale<0) continue;
		    const double thetax=parms->wfs[iwfs].thetax-parms->powfs[ipowfs].llt->ox->p[illt]/hs;
		    const double thetay=parms->wfs[iwfs].thetay-parms->powfs[ipowfs].llt->oy->p[illt]/hs;
		    const double displacex=-atm->p[ips]->vx*isim*dt+thetax*hl+parms->powfs[ipowfs].llt->misreg->p[0];
		    const double displacey=-atm->p[ips]->vy*isim*dt+thetay*hl+parms->powfs[ipowfs].llt->misreg->p[1];
		    prop_grid_pts(atm->p[ips],powfs[ipowfs].llt->pts,
				  lltopd->p,atmscale,displacex,displacey,
				  scale, 1., 0, 0);
		}
	    }
	    double ttx=0, tty=0;//FSM + wind shake induced jitter
	    if((simu->fsmreal && simu->fsmreal->p[iwfs]) ||do_pistatout||parms->sim.idealfsm){
		if(do_pistatout||parms->sim.idealfsm){
		    /* remove tip/tilt completely */
		    dmat *lltg=dnew(2,1);
		    pts_ztilt(&lltg,powfs[ipowfs].llt->pts,
			      powfs[ipowfs].llt->imcc,
			      powfs[ipowfs].llt->amp->p,
			      lltopd->p);
		    simu->fsmreal->p[iwfs]->p[0]=-lltg->p[0];
		    simu->fsmreal->p[iwfs]->p[1]=-lltg->p[1];
		    dfree(lltg);
		}
		ttx=simu->fsmreal->p[iwfs]->p[0];
		tty=simu->fsmreal->p[iwfs]->p[1];
	    }
	    if(simu->telws){
		double tmp=simu->telws->p[isim]*parms->powfs[ipowfs].llt->ttrat;
		double angle=simu->winddir?simu->winddir->p[0]:0;
		ttx+=tmp*cos(angle);
		tty+=tmp*sin(angle);
	    }
	    if(simu->llt_tt && simu->llt_tt->p[iwfs]){
		ttx+=simu->llt_tt->p[iwfs]->p[isim];//put all to x direction.
	    }
	    if(ttx !=0 || tty != 0){ /* add tip/tilt to llt opd */
		double ptt[3]={0, ttx, tty};
		loc_add_ptt(lltopd->p, ptt, powfs[ipowfs].llt->loc);
	    }
	    if(save_opd){
		zfarr_push(simu->save->wfslltopd[iwfs], isim, lltopd);
	    }
	}
	if(parms->powfs[ipowfs].type==0){
	    WFSINTS_T *intsdata=simu->wfs_intsdata+iwfs;
	    intsdata->ints=ints;
	    intsdata->psfout=psfout;
	    intsdata->pistatout=simu->pistatout->p[iwfs];
	    if(parms->powfs[ipowfs].pistatout){
		intsdata->gradref=gradcalc;
	    }
	    intsdata->opd=opd;
	    intsdata->lltopd=lltopd;
	    intsdata->isim=isim;
	    CALL_THREAD(simu->wfs_ints[iwfs], 0);
	    dfree(lltopd);
	    intsdata->opd=0;
	    intsdata->lltopd=0;
	    if(psfout){
		zfarr_push(psfoutzfarr, isim, psfout);
		zfarr_push(ztiltoutzfarr, isim, *gradacc);
	    }
	}else{//Pywfs
	    pywfs_fft(&ints->p[0], powfs[ipowfs].pywfs, opd);
	    dscale(ints->p[0], parms->wfs[iwfs].siglevsim);
	}
    }
    TIM(3);
    if(dtrat_output){
	const double rne=parms->powfs[ipowfs].rne;
	const double bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
	if(do_phy){
	    /* In Physical optics mode, do integration and compute
	       gradients. The matched filter are in x/y coordinate even if
	       radpix=1. */
	    if(save_ints){
		zfarr_push(simu->save->intsnf[iwfs], isim, ints);
	    }
	    if(noisy){/*add noise */
		const double bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
		dmat **bkgrnd2=NULL;
		dmat **bkgrnd2c=NULL;
		if(powfs[ipowfs].bkgrnd){
		    if(powfs[ipowfs].bkgrnd->ny==1){
			bkgrnd2=powfs[ipowfs].bkgrnd->p;
		    }else{
			bkgrnd2=powfs[ipowfs].bkgrnd->p+nsa*wfsind;
		    }
		}
		if(powfs[ipowfs].bkgrndc){
		    if(powfs[ipowfs].bkgrndc->ny==1){
			bkgrnd2c=powfs[ipowfs].bkgrndc->p;
		    }else{
			bkgrnd2c=powfs[ipowfs].bkgrndc->p+nsa*wfsind;
		    }
		}
		for(int isa=0; isa<ints->nx; isa++){
		    dmat *bkgrnd2i=(bkgrnd2)?bkgrnd2[isa]:NULL;
		    dmat *bkgrnd2ic=(bkgrnd2c)?bkgrnd2c[isa]:NULL;
		    addnoise(ints->p[isa], &simu->wfs_rand[iwfs],
			     bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, parms->powfs[ipowfs].qe, rne, 1.);
		}
		if(save_ints){
		    zfarr_push(simu->save->intsny[iwfs], isim, ints);
		}
	    }
	    if(parms->powfs[ipowfs].i0save==2){
		dcelladd(&simu->ints->p[iwfs], 1, ints, 1);
	    }
	    if(parms->powfs[ipowfs].dither==1 && isim>=parms->powfs[ipowfs].dither_ogskip
	       && parms->powfs[ipowfs].type==0 && parms->powfs[ipowfs].phytype_sim2==1){
		/*Collect statistics with dithering*/
		DITHER_T *pd=simu->dither[iwfs];
		double cs, ss;
		dither_position(&cs, &ss, parms, ipowfs, isim, pd->deltam);
		//accumulate for matched filter
		dcelladd(&pd->imb, 1, ints, 1.);
		dcelladd(&pd->imx, 1, ints, cs);
		dcelladd(&pd->imy, 1, ints, ss);
	    }
	}
	if(do_phy){

	    if(parms->powfs[ipowfs].type==0){
		calc_phygrads(gradout, ints->p, parms, powfs, iwfs, parms->powfs[ipowfs].phytype_sim);
	    }else{
		pywfs_grad(gradout, powfs[ipowfs].pywfs, ints->p[0]);
	    }
	}else{
	    /* geomtric optics accumulation mode. scale and copy results to output. */
	    dcp(gradout,*gradacc);
	    if(dtrat!=1){
		dscale(*gradout,1./dtrat);/*average */
	    }
	    if(noisy && !parms->powfs[ipowfs].usephy){
		const dmat *nea=PR(powfs[ipowfs].neasim,wfsind,0);//neasim is the LL' decomposition
		const double *neax=nea->p;
		const double *neay=nea->p+nsa;
		const double *neaxy=nea->p+nsa*2;
		double *restrict ggx=(*gradout)->p;
		double *restrict ggy=(*gradout)->p+nsa;
		for(int isa=0; isa<nsa; isa++){
		    /*Preserve the random sequence. */
		    double n1=randn(&simu->wfs_rand[iwfs]);
		    double n2=randn(&simu->wfs_rand[iwfs]);
		    double errx=neax[isa]*n1;
		    double erry=neay[isa]*n2+neaxy[isa]*n1;/*cross term. */
		    ggx[isa]+=errx;
		    ggy[isa]+=erry;
		}
	    }
	}
	if(save_gradgeom){
	    dmat *gradtmp=NULL;
	    dadd(&gradtmp, 1, *gradacc, 1./dtrat);
	    zfarr_push(simu->save->gradgeom[iwfs], isim, gradtmp);/*noise free. */
	    dfree(gradtmp);
	}
    }//dtrat_out
    dfree(gradcalc);
    TIM(4);
#if TIMING==1
    info("WFS %d grad timing: atm %.2f dm %.2f ints %.2f grad %.2f\n",iwfs,tk1-tk0,tk2-tk1,tk3-tk2,tk4-tk3);
#endif
}
/**
   Calculate the amplitude of dithering signal or response based on frequency of the dithering signal.
*/
static double calc_dither_amp(dmat *dithersig, /**<array of data. nmod*nsim */
			      long dtrat,   /**<skip columns due to wfs/sim dt ratio*/
			      long npoint,  /**<number of points during dithering*/
			      int detrend   /**<flag for detrending (remove linear dithersig)*/
    ){
    const long nmod=dithersig->nx;
    long nframe=(dithersig->ny-1)/dtrat+1;//number of wavefront frames
    double slope=0;
    long offset=(nframe/npoint-1)*npoint;//number of WFS frame separations between first and last cycle
    if(detrend && offset){//detrending
	for(long ip=0; ip<npoint; ip++){
	    for(long im=0; im<nmod; im++){
		long i0=ip*dtrat*nmod+im;
		long i1=(ip+offset)*dtrat*nmod+im;
		slope+=dithersig->p[i1]-dithersig->p[i0];
	    }
	}
	slope/=(npoint*nmod*offset);
	//dbg("slope=%g. npoint=%ld, nmod=%ld, nframe=%ld, offset=%ld\n", slope, npoint, nmod, nframe, offset);
    }
    double anglei=M_PI*2/npoint;
    double ipv=0, qdv=0;
    if(nmod==2){//tip and tilt
	for(int iframe=0; iframe<nframe; iframe++){
	    double angle=anglei*iframe;//position of dithering
	    double cs=cos(angle);
	    double ss=sin(angle);
	    double ttx=dithersig->p[iframe*dtrat*2]-slope*iframe;
	    double tty=dithersig->p[iframe*dtrat*2+1]-slope*iframe;
	    ipv+=(ttx*cs+tty*ss);
	    qdv+=(ttx*ss-tty*cs);
	}
    }else if(nmod==1){//single mode dithering
	for(int iframe=0; iframe<nframe; iframe++){
	    double angle=anglei*iframe;//position of dithering
	    double cs=cos(angle);
	    double ss=sin(angle);
	    double mod=dithersig->p[iframe*dtrat]-slope*iframe;
	    ipv+=(mod*cs);
	    qdv+=(mod*ss);
	}
    }
    double a2m=sqrt(ipv*ipv+qdv*qdv)/nframe;
    return a2m;
}

/*Fast steering mirror for each WFS*/
void wfsgrad_fsm(SIM_T *simu, int iwfs){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int isim=simu->isim;
    /*Uplink FSM*/
    const int ind=parms->recon.glao?(ipowfs+ipowfs*parms->npowfs):(iwfs+iwfs*parms->nwfs);
    dmat *PTT=recon->PTT?(recon->PTT->p[ind]):0;
    if(!PTT){
	error("powfs %d has FSM, but PTT is empty\n", ipowfs);
    }
    /* Compute FSM error. */
    simu->fsmerr=simu->fsmerr_store;
    dmm(&simu->fsmerr->p[iwfs], 0, PTT, simu->gradcl->p[iwfs], "nn", 1);
    //Save data
    P(simu->fsmerrs->p[iwfs], 0, isim)=simu->fsmerr->p[iwfs]->p[0];
    P(simu->fsmerrs->p[iwfs], 1, isim)=simu->fsmerr->p[iwfs]->p[1];
}

/**
   Postprocessing for dithering dithersig extraction:
   1. Every step: accumulate signal for phase detection.
   2. At PLL output: determine input/output amplitude of dithering signal.
   3. At Gain output:determine matched filter i0, gx, gy, or CoG gain.
   4. Subtract t/t from gradients for non-comon-path (TT) dithering.

*/
static void wfsgrad_dither(SIM_T *simu, int iwfs){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    POWFS_T *powfs=simu->powfs;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int isim=simu->isim;
    const int npll=parms->powfs[ipowfs].dither_pllrat;
    if(!parms->powfs[ipowfs].dither || isim<parms->powfs[ipowfs].dither_pllskip){
	return;
    }
    DITHER_T *pd=simu->dither[iwfs];
    if(parms->powfs[ipowfs].dither==1){ //T/T dithering.
	double cs, ss;
	//Compute gradient derivative for every subaperture
	dither_position(&cs, &ss, parms, ipowfs, isim, pd->deltam);
	if(parms->powfs[ipowfs].type==0 && parms->powfs[ipowfs].phytype_sim2!=1){
	    const int nsa=powfs[ipowfs].saloc->nloc;
	    if(!pd->ggm){
		pd->ggm=dnew(nsa*2,1);
	    }
	    for(int isa=0; isa<nsa; isa++){
		pd->ggm->p[isa]+=cs*simu->gradcl->p[iwfs]->p[isa];
		pd->ggm->p[isa+nsa]+=ss*simu->gradcl->p[iwfs]->p[isa+nsa];
	    }
	}
	    
	/*Use delay locked loop to determine the phase of actual dithering
	  dithersig from WFS error dithersig, i.e., average position of expected
	  spot during integration. Uplink propagation is accounted for in
	  LGS.*/
	double err, cd, sd;
	dither_position(&cd, &sd, parms, ipowfs, isim, pd->deltam);
	err=(-sd*(simu->fsmerr->p[iwfs]->p[0])
	     +cd*(simu->fsmerr->p[iwfs]->p[1]))/(parms->powfs[ipowfs].dither_amp);
	pd->delta+=parms->powfs[ipowfs].dither_gpll*(err/npll);
    }else if(parms->powfs[ipowfs].dither>1){ //DM dithering.
	//Compute dither signal strength in input (DM) and output (gradients) by correlation.
	dmat *tmp=0;
	const int idm=parms->idmground;
	dmm(&tmp, 0, P(recon->dither_ra, idm, idm), simu->dmreal->p[idm], "nn", 1);
	simu->dither[iwfs]->mr->p[0]->p[isim]=tmp->p[0];
	dmm(&tmp, 0, P(recon->dither_rg, iwfs, iwfs), simu->gradcl->p[iwfs], "nn", 1);
	simu->dither[iwfs]->mr->p[1]->p[isim]=tmp->p[0];
	dfree(tmp);
    }
    
    int npllacc=(simu->isim-parms->powfs[ipowfs].dither_pllskip+1)/parms->powfs[ipowfs].dtrat;
    if(npllacc >0 && npllacc%npll==0){
	//Synchronous detection of dither dithersig in input (DM) and output
	//(gradients) dithersig. The ratio between the two is used for optical gain adjustment.
	const int npoint=parms->powfs[ipowfs].dither_npoint;
	const int ncol=(npll-1)*parms->powfs[ipowfs].dtrat+1;
	if(parms->powfs[ipowfs].dither==1){//TT
	    pd->deltam=pd->delta+pd->deltao;
	    dmat *tmp=0;
	    const double norm=1./parms->powfs[ipowfs].dither_amp;
	    const int detrend=parms->powfs[ipowfs].llt?0:1;
	    tmp=drefcols(simu->fsmcmds->p[iwfs], simu->isim-ncol+1, ncol);
	    pd->a2m=calc_dither_amp(tmp, parms->powfs[ipowfs].dtrat, npoint, detrend)*norm;
	    dfree(tmp);
	    tmp=drefcols(simu->fsmerrs->p[iwfs], simu->isim-ncol+1, ncol);
	    pd->a2me=calc_dither_amp(tmp, parms->powfs[ipowfs].dtrat, npoint, detrend)*norm;
	    dfree(tmp);
	}else if(parms->powfs[ipowfs].dither>1){//DM
	    dmat *tmp=0;
	    tmp=drefcols(pd->mr->p[0], simu->isim-ncol+1, ncol);//DM
	    pd->a2m=calc_dither_amp(tmp, parms->powfs[ipowfs].dtrat, npoint, 1);
	    dfree(tmp);
	    tmp=drefcols(pd->mr->p[1], simu->isim-ncol+1, ncol);//Grad
	    pd->a2me=calc_dither_amp(tmp, parms->powfs[ipowfs].dtrat, npoint, 1);
	    dfree(tmp);
	}
	    
	if(iwfs==parms->powfs[ipowfs].wfs->p[0]){
	    const double anglei=(2*M_PI/parms->powfs[ipowfs].dither_npoint);
	    info("Step %d wfs%d PLL: delay=%.2f frame, dither amplitude=%.2f, estimate=%.2f\n",
		  isim, iwfs, pd->deltam/anglei, pd->a2m, pd->a2me);
	}
	if(simu->resdither){
	    int ic=(npllacc-1)/(npll);
	    P(simu->resdither->p[iwfs], 0, ic)=pd->deltam;
	    P(simu->resdither->p[iwfs], 1, ic)=pd->a2m;
	    P(simu->resdither->p[iwfs], 2, ic)=pd->a2me;
	}
    }
    int nogacc=(simu->isim-parms->powfs[ipowfs].dither_ogskip+1)/parms->powfs[ipowfs].dtrat;
    if(nogacc>0 && nogacc%npll==0){//Gain update output
	if(parms->powfs[ipowfs].dither==1){//TT Dither
	    double scale1=1./npll;
	    double amp=pd->a2m*parms->powfs[ipowfs].dither_amp;
	    double scale2=scale1*2./(amp);
	    if(pd->imb){//matched filter
		dcellscale(pd->imb, scale1);
		dmat *ibgrad=0;
		calc_phygrads(&ibgrad, pd->imb->p, parms, powfs, iwfs, 2);
		if(parms->powfs[ipowfs].trs){//tip/tilt drift dithersig
		    /* Update drift mode computation. Only useful when wfs t/t is removed*/
		    dmat *tt=dnew(2,1);
		    const int ind=parms->recon.glao?(ipowfs+ipowfs*parms->npowfs):(iwfs+iwfs*parms->nwfs);
		    dmat *PTT=recon->PTT?(recon->PTT->p[ind]):0;
		    dmm(&tt, 0, PTT, ibgrad, "nn", 1);
		    simu->fsmerr->p[iwfs]->p[0]+=tt->p[0];
		    simu->fsmerr->p[iwfs]->p[1]+=tt->p[1];
		    dfree(tt);
		}
		//Smooth trombone movement by provide continuous err.
		if(parms->powfs[ipowfs].llt){
		    dmat *focus=dnew(1,1);
		    dmat *RFlgsg=P(recon->RFlgsg, iwfs, iwfs);
		    dmm(&focus, 0, RFlgsg, ibgrad, "nn", 1);
		    //zoomerr is adjust by the gain, and scaling due to zoh for npll frames.
		    simu->zoomerr->p[iwfs]=focus->p[0]*(parms->powfs[ipowfs].zoomgain/npll);
		    //if(iwfs==0) info("zoomerr_ib[%d]=%g nm\n", isim, focus->p[0]*1e9);
		    dfree(focus);
		}
		dfree(ibgrad);
		//Accumulate data for matched filter
		dcelladd(&pd->i0, 1, pd->imb, 1);//imb was already scaled
		dcelladd(&pd->gx, 1, pd->imx, scale2);
		dcelladd(&pd->gy, 1, pd->imy, scale2);
		dcellzero(pd->imb);
		dcellzero(pd->imx);
		dcellzero(pd->imy);
	    }else if(pd->ggm){//cog
		dadd(&pd->gg0, 1, pd->ggm, scale2);
		dzero(pd->ggm);
	    }
	}//Drift mode computation
    }

    if(!parms->powfs[ipowfs].trs){
	/*when WFS t/t is used for reconstruction, do not close FSM
	 * loop. Subtract actual dithering dithersig.*/
	if(parms->powfs[ipowfs].dither==1 && simu->gradscale->p[iwfs]){
	    double amp,cs,ss; 
	    dither_position(&cs, &ss, parms, ipowfs, isim, pd->deltam);
	    amp=pd->a2me*parms->powfs[ipowfs].dither_amp;
	    double ptt[2]={-cs*amp, -ss*amp};
	    dmulvec(simu->gradcl->p[iwfs]->p, recon->TT->p[iwfs+iwfs*parms->nwfsr], ptt, 1);
	}
    }
}

/**
   Accomplish Two tasks:
   1) Use LPF'ed LGS focus measurement to drive the trombone.
   2) HPF lgs focus on gradients to remove sodium range variation effact.
   
   We trust the focus measurement of the LGS WFS at high temporal frequency
   which NGS cannot provide due to low frame rate. After the HPF on lgs
   gradients, our system is NO LONGER affected by sodium layer variation.

   if sim.mffocus==1: The HPF is applied to each LGS WFS indepently. This largely
   removes the effect of differential focus. powfs.dfrs is no longer needed. (preferred).

   if sim.mffocus==2: We apply a LPF on the average focus from six LGS WFS, and
   then remove this value from all LGS WFS. The differential focus is still
   present and powfs.dfrs need to be set to 1 to handle it in tomography. This
   is the original focus tracking method, and is no longer recommended.
*/
static void wfsgrad_lgsfocus(SIM_T* simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    
    /*New plate mode focus offset for LGS WFS. Not really needed*/
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].llt && parms->tomo.ahst_focus==2 
	   && simu->Mint_lo && simu->Mint_lo->mint->p[1]
	   && simu->isim+1>parms->powfs[ipowfs].step
	   && (simu->isim+1)%parms->powfs[ipowfs].dtrat==0){
	    /*When tomo.ahst_focus>0, the first plate scale mode contains focus for
	      lgs. But it turns out to be not necessary to remove it because the
	      HPF in the LGS path removed the influence of this focus mode. set
	      tomo.ahst_focus=2 to enable adjust gradients.*/

	    double scale=simu->recon->ngsmod->scale;
	    int indps=simu->recon->ngsmod->indps;
	    dmat *mint=simu->Mint_lo->mint->p[0]->p[0];//2018-12-11: changed first p[1] to p[0]
	    double focus=mint->p[indps]*(scale-1); 
	    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[iwfs], focus);
	    }
	}
    }


    dcell *LGSfocus=simu->LGSfocus;//computed in wfsgrad_post from gradcl.
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(simu->isim<parms->powfs[ipowfs].step) continue;
	if(parms->sim.closeloop && (simu->isim+1)%parms->sim.dtrat_hi!=0) continue;
	if(!parms->powfs[ipowfs].llt) continue;
	const int do_phy=(parms->powfs[ipowfs].usephy && simu->isim>=parms->powfs[ipowfs].phystep);
	if(!do_phy) continue;
	double lgsfocusm=0;
	if(parms->powfs[ipowfs].zoomshare || parms->sim.mffocus==2){
	    lgsfocusm=0;
	    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		lgsfocusm+=LGSfocus->p[iwfs]->p[0];
	    }
	    lgsfocusm/=parms->powfs[ipowfs].nwfs;
	}
	if(parms->powfs[ipowfs].zoomset && simu->isim==parms->sim.start && parms->powfs[ipowfs].phytype_sim!=1){
	    /*Here we set trombone position according to focus in the first
	      measurement. And adjust the focus content of this
	      measurement. This simulates the initial focus acquisition
	      step. No need if start with pre-built matched filter.*/
	    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		if(parms->powfs[ipowfs].zoomshare){
		    simu->zoomint->p[iwfs]=lgsfocusm;
		}else{
		    simu->zoomint->p[iwfs]=LGSfocus->p[iwfs]->p[0]; 
		}
		LGSfocus->p[iwfs]->p[0]-=simu->zoomint->p[iwfs];
		dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[iwfs], -simu->zoomint->p[iwfs]);
		info("wfs %d: Set trombone position to %g.\n", iwfs, simu->zoomint->p[iwfs]);
	    }
	}
	//In RTC. LPF can be put after using the value to put it off critical path.
	double lpfocus=parms->sim.lpfocushi;
	double infocus=0;
	//For those WFS that dither and run mtch, use focus error from ib instead
	const int zoom_from_grad=(parms->powfs[ipowfs].dither!=1 || parms->powfs[ipowfs].phytype_sim2!=1);
	const int zoomavg_output=((simu->reconisim+1)%parms->powfs[ipowfs].zoomdtrat==0);
	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    if(parms->sim.mffocus){//Focus HPF
		if(parms->sim.mffocus==1){//remove LPF focus from each lgs
		    infocus=LGSfocus->p[iwfs]->p[0];
		}else{//remove powfs averaged focus form each lgs.
		    infocus=lgsfocusm;
		}
		simu->lgsfocuslpf->p[iwfs]=simu->lgsfocuslpf->p[iwfs]*(1-lpfocus)+infocus*lpfocus;
		dadd(&simu->gradcl->p[iwfs], 1, recon->GFall->p[iwfs], -simu->lgsfocuslpf->p[iwfs]);
	    }
	    if(zoom_from_grad){//Trombone averager
		if(parms->powfs[ipowfs].zoomshare){//all lgs share the same trombone.
		    simu->zoomavg->p[iwfs]+=lgsfocusm;
		}else{//each lgs has individual zoom mechanism.
		    simu->zoomavg->p[iwfs]+=LGSfocus->p[iwfs]->p[0];
		}
		if(zoomavg_output){
		    /*zoom error is zero order hold even if no output from averager*/
		    simu->zoomerr->p[iwfs]=simu->zoomavg->p[iwfs]*parms->powfs[ipowfs].zoomgain*pow(parms->powfs[ipowfs].zoomdtrat,-2);
		    //info("zoomerr_grad[%d]=%g nm\n", simu->isim, simu->zoomerr->p[0]/parms->powfs[ipowfs].zoomgain*parms->powfs[ipowfs].zoomdtrat*1e9);
		    simu->zoomavg->p[iwfs]=0;
		}
	    }
	}
    }
    
    /*The zoomerr is prescaled, and ZoH. moved from filter.c as it only relates to WFS.*/
    dcp(&simu->zoomreal, simu->zoomint);
    dadd(&simu->zoomint, 1, simu->zoomerr, 1);
}

/**
   Every operation here should be in the Simulator not the Controller 
*/
void wfsgrad_post(thread_t *info){
    SIM_T *simu=(SIM_T*)info->data;
    const PARMS_T *parms=simu->parms;
    //Postprocessing gradients
    const int isim=simu->isim;
    for(int iwfs=info->start; iwfs<info->end; iwfs++){
#if USE_CUDA
	if(parms->gpu.wfs){
	    gpu_wfsgrad_sync(simu, iwfs);
	}
#endif
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int dtrat=parms->powfs[ipowfs].dtrat;
	if(isim<parms->powfs[ipowfs].step) continue;
	const int dtrat_output=((isim+1)%dtrat==0);
	const int do_phy=(parms->powfs[ipowfs].usephy && isim>=parms->powfs[ipowfs].phystep);
	dmat **gradout=&simu->gradcl->p[iwfs];
	/* copy fsmreal to output  */
	if(simu->fsmreal && simu->fsmreal->p[iwfs]){
	    P(simu->fsmcmds->p[iwfs], 0, isim)=simu->fsmreal->p[iwfs]->p[0];
	    P(simu->fsmcmds->p[iwfs], 1, isim)=simu->fsmreal->p[iwfs]->p[1];
	}
	if(dtrat_output){
	    if(parms->plot.run>1 && (simu->gradoff->p[iwfs]|| parms->dbg.gradoff)){
		drawopd("Gclx", simu->powfs[ipowfs].saloc, simu->gradcl->p[iwfs]->p,
			parms->dbg.draw_gmax->p,
			"WFS Closeloop Gradients (x)","x (m)", "y (m)",
			"raw x %d",  iwfs);
		drawopd("Gcly", simu->powfs[ipowfs].saloc, simu->gradcl->p[iwfs]->p+
			simu->powfs[ipowfs].saloc->nloc,
			parms->dbg.draw_gmax->p,
			"WFS Closeloop Gradients (y)","x (m)", "y (m)",
			"raw y %d",  iwfs);
	    }
	    //scaling on gradout
	    if(simu->gradscale->p[iwfs]){
		dcwm(*gradout, simu->gradscale->p[iwfs]);
	    }else{
		dscale(*gradout, parms->powfs[ipowfs].gradscale);
	    }
	    //Gradient offset due to mainly NCPA calibration. Must be after gain adjustment.
	    if(simu->gradoff->p[iwfs]){
		dadd(gradout, 1, simu->gradoff->p[iwfs], -parms->dbg.gradoff_scale);
	    }
	    //Injected gradient offset for testing.
	    if(parms->dbg.gradoff){
		info_once("Add injected gradient offset vector\n");
		int icol=(isim+1)%parms->dbg.gradoff->ny;
		dadd(gradout, 1, P(parms->dbg.gradoff, iwfs, icol), -1);
	    }
	    if(do_phy){
		if(simu->fsmerr_store->p[iwfs]){
		    wfsgrad_fsm(simu, iwfs);
		}
		if(parms->powfs[ipowfs].dither){
		    wfsgrad_dither(simu, iwfs);
		}
		if(!parms->powfs[ipowfs].trs && parms->powfs[ipowfs].skip!=2 && simu->fsmerr){
		    dzero(simu->fsmerr->p[iwfs]);//do not close fsm loop
		}
		if(parms->powfs[ipowfs].llt){
		    dmm(PP(simu->LGSfocus, iwfs), 0, P(simu->recon->RFlgsg, iwfs, iwfs), P(simu->gradcl, iwfs), "nn", 1);
		}
	    }
	    if(parms->save.grad->p[iwfs]){
		zfarr_push(simu->save->gradcl[iwfs], isim, simu->gradcl->p[iwfs]);
	    }
	    if(parms->plot.run){
		drawopd("Gclx", simu->powfs[ipowfs].saloc, simu->gradcl->p[iwfs]->p,
			parms->dbg.draw_gmax->p,
			"WFS Closeloop Gradients (x)","x (m)", "y (m)",
			"x %d",  iwfs);
		drawopd("Gcly", simu->powfs[ipowfs].saloc, simu->gradcl->p[iwfs]->p+
			simu->powfs[ipowfs].saloc->nloc,
			parms->dbg.draw_gmax->p,
			"WFS Closeloop Gradients (y)","x (m)", "y (m)",
			"y %d",  iwfs);
		if(do_phy){
		    if(simu->ints->p[iwfs]->nx==1){
			ddraw("Ints", simu->ints->p[iwfs]->p[0], NULL,NULL, "WFS Subaperture Images",
			      "x", "y", "wfs %d", iwfs);
		    }else if(simu->ints->p[iwfs]->nx==4){
			cellreshape(simu->ints->p[iwfs], 2, 2);
			dmat *ints2=dcell2m(simu->ints->p[iwfs]);
			cellreshape(simu->ints->p[iwfs], 4, 1);
			ddraw("Ints", ints2, NULL, NULL, "WFS Subaperture Images",
			      "x", "y", "wfs %d", iwfs);
			dfree(ints2);
		    }
		}
	    }
	}
    }//for iwfs
}

/**
   Dither update: zoom corrector, matched filter, gain ajustment, TWFS.
*/
static void wfsgrad_dither_post(SIM_T *simu){
    POWFS_T *powfs=simu->powfs;
    const PARMS_T *parms=simu->parms;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(!parms->powfs[ipowfs].dither) continue;
	if(simu->isim<parms->powfs[ipowfs].step) continue;
	if((simu->isim+1)%parms->powfs[ipowfs].dtrat!=0) continue;
	const int npllacc=(simu->isim-parms->powfs[ipowfs].dither_pllskip+1)/parms->powfs[ipowfs].dtrat;
	const int nwfs=parms->powfs[ipowfs].nwfs;
	const int npll=parms->powfs[ipowfs].dither_pllrat;
	
	if(parms->powfs[ipowfs].zoomshare && parms->powfs[ipowfs].llt //this is LLT
	   && npllacc>0 && npllacc % npll==0){//There is drift mode computation
	    //Average zoomerr computed from wfsgrad_dither.
	    double sum=0;
	    for(int jwfs=0; jwfs<nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		sum+=simu->zoomerr->p[iwfs];
	    }
	    sum/=(double)nwfs;
	    for(int jwfs=0; jwfs<nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		simu->zoomerr->p[iwfs]=sum;
	    }
	}
	int ngrad=parms->powfs[ipowfs].dither_ograt;
	const int nogacc=(simu->isim-parms->powfs[ipowfs].dither_ogskip+1)/parms->powfs[ipowfs].dtrat;
	if(nogacc>0 && nogacc % ngrad==0){//This is matched filter or cog update
	    const int nsa=powfs[ipowfs].saloc->nloc;
	    double scale1=(double)parms->powfs[ipowfs].dither_pllrat/(double)parms->powfs[ipowfs].dither_ograt;
	    if(parms->powfs[ipowfs].phytype_sim2==1 && parms->powfs[ipowfs].type==0){
		info("Step %d: Update matched filter for powfs %d\n", simu->isim, ipowfs);
		//For matched filter
		if(!powfs[ipowfs].intstat){
		    powfs[ipowfs].intstat=mycalloc(1,INTSTAT_T);
		}
		parms->powfs[ipowfs].radgx=0;//ensure derivate is interpreted as along x/y.
		if(!powfs[ipowfs].intstat->i0 || powfs[ipowfs].intstat->i0->ny!=nwfs){
		    dcellfree(powfs[ipowfs].intstat->i0);
		    dcellfree(powfs[ipowfs].intstat->gx);
		    dcellfree(powfs[ipowfs].intstat->gy);

		    powfs[ipowfs].intstat->i0=dcellnew(nsa, nwfs);
		    powfs[ipowfs].intstat->gx=dcellnew(nsa, nwfs);
		    powfs[ipowfs].intstat->gy=dcellnew(nsa, nwfs);
		}

		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    DITHER_T *pd=simu->dither[iwfs];
		    //Scale the output due to accumulation
		    for(int isa=0; isa<nsa; isa++){
			dadd(powfs[ipowfs].intstat->i0->p+isa+jwfs*nsa, 0, pd->i0->p[isa], scale1);
			dadd(powfs[ipowfs].intstat->gx->p+isa+jwfs*nsa, 0, pd->gx->p[isa], scale1);
			dadd(powfs[ipowfs].intstat->gy->p+isa+jwfs*nsa, 0, pd->gy->p[isa], scale1);
		    }
		    dcellzero(pd->i0);
		    dcellzero(pd->gx);
		    dcellzero(pd->gy);
		    dmat *goff=0;
		    if(1){//Always use this
			/*Compute the gradient of i0 using old gradient algorithm
			  and subtract from the gradient offset to prevent sudden
			  jump of gradient measurement.*/
			calc_phygrads(&goff, powfs[ipowfs].intstat->i0->p+jwfs*nsa, parms, powfs, iwfs, parms->powfs[ipowfs].phytype_sim);
			dadd(&simu->gradoff->p[iwfs], 1, goff, -1);
		    }else{
			dzero(simu->gradoff->p[iwfs]);
		    }
		    if(parms->save.dither){
			writebin(simu->gradoff->p[iwfs], "wfs%d_gradoff_%d", iwfs, simu->isim);
		    }
		    if(parms->powfs[ipowfs].dither_gdrift>0){
			//outer loop to prevent i0 from drifting in COG
			dzero(goff);
			//Compute CoG of i0 + goff and drive it toward gradncpa with low gain (0.1)
			calc_phygrads(&goff, powfs[ipowfs].intstat->i0->p+jwfs*nsa, parms, powfs, iwfs, 2);
			dadd(&goff, 1, simu->gradoff->p[iwfs], 1);
			if(powfs[ipowfs].gradncpa){
			    dadd(&goff, 1, powfs[ipowfs].gradncpa->p[jwfs], -1);
			}
			if(parms->powfs[ipowfs].llt){
			    //Remove focus drift control in LGS WFS as it is fixed using HFP and trombone.
			    dmat *focus=dnew(1,1);
			    dmm(&focus, 0, P(simu->recon->RFlgsg, iwfs, iwfs), goff, "nn", 1);
			    info("Step %d, wfs %d: removing focus=%g\n", simu->isim, iwfs, focus->p[0]);
			    dadd(&goff, 1, simu->recon->GFall->p[iwfs], -focus->p[0]);
			    dfree(focus);
			}	
			dadd(&simu->gradoff->p[iwfs], 1, goff, -parms->powfs[ipowfs].dither_gdrift);
			//outer loop to prevent i0 derivative from drifting.
			dmat *i0sx=0, *i0sy=0;
			double theta=0;
			const double gshift=parms->powfs[ipowfs].pixtheta*0.1;
			for(int isa=0; isa<nsa; isa++){
			    double g0[3], gx[3], gy[3];
			    dcp(&i0sx, PR(powfs[ipowfs].intstat->i0, isa, jwfs));
			    dcog(g0, i0sx,0.,0., powfs[ipowfs].cogcoeff->p[jwfs]->p[isa*2], powfs[ipowfs].cogcoeff->p[jwfs]->p[isa*2+1], 0);
			    dcp(&i0sy, PR(powfs[ipowfs].intstat->i0, isa, jwfs));
			    dadd(&i0sx, 1, PR(powfs[ipowfs].intstat->gx, isa, jwfs), gshift);
			    dadd(&i0sy, 1, PR(powfs[ipowfs].intstat->gy, isa, jwfs), gshift);
			    dcog(gx, i0sx,0.,0., powfs[ipowfs].cogcoeff->p[jwfs]->p[isa*2], powfs[ipowfs].cogcoeff->p[jwfs]->p[isa*2+1], 0);
			    dcog(gy, i0sy,0.,0., powfs[ipowfs].cogcoeff->p[jwfs]->p[isa*2], powfs[ipowfs].cogcoeff->p[jwfs]->p[isa*2+1], 0);
			    //Works in both x/y and r/a coordinate.
			    theta+=(atan2(gx[1]-g0[1], gx[0]-g0[0])+atan2(gy[1]-g0[1], gy[0]-g0[0]));
			}
			theta*=0.5/nsa;
			pd->deltao=-theta;
			info("wfs[%d] deltao is %g.\n", iwfs, pd->deltao);
		    }
		    dfree(goff);
		    if(parms->save.dither){
			writebin(simu->gradoff->p[iwfs], "wfs%d_gradoff_drift_%d", iwfs, simu->isim);
		    }
		}
		if(parms->save.dither){
		    writebin(powfs[ipowfs].intstat->i0, "powfs%d_i0_%d", ipowfs, simu->isim);
		    writebin(powfs[ipowfs].intstat->gx, "powfs%d_gx_%d", ipowfs, simu->isim);
		    writebin(powfs[ipowfs].intstat->gy, "powfs%d_gy_%d", ipowfs, simu->isim);
		}
	    }else{
		//For CoG gain
		for(int jwfs=0; jwfs<nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    DITHER_T *pd=simu->dither[iwfs];
		    const int ng=powfs[ipowfs].saloc->nloc*2;
		    if(!simu->gradscale->p[iwfs]){
			simu->gradscale->p[iwfs]=dnew(ng,1);
			dset(simu->gradscale->p[iwfs], parms->powfs[ipowfs].gradscale);
		    }
		    double mgold=dsum(simu->gradscale->p[iwfs])/ng;
		    double mgnew;
		    const char *ogtype=0;
		    //gg0 is output/input of dither dithersig.
		    if(!pd->gg0 || parms->powfs[ipowfs].dither_ogsingle){//single gain for all subapertures. For Pyramid WFS
			ogtype="globally";
			double gerr=pd->a2me/pd->a2m;
#define HIA_G_UPDATE 0
#if HIA_G_UPDATE //HIA method.
			double adj=parms->powfs[ipowfs].dither_gog*mgold*(1-gerr);
			dadds(simu->gradscale->p[iwfs], adj);
			mgnew=mgold+adj;
			while(mgnew<0){//prevent negative gain
			    adj*=0.5;
			    dadds(simu->gradscale->p[iwfs], -adj);
			    mgnew+=-adj;
			}
#else
			double adj=pow(gerr, -parms->powfs[ipowfs].dither_gog);
			dscale(simu->gradscale->p[iwfs], adj);
			mgnew=mgold*adj;
#endif
		    }else{//separate gain for each gradient. For shwfs.
			ogtype="on average";
			dscale(pd->gg0, scale1); //Scale value at end of accumulation
			for(long ig=0; ig<ng; ig++){
#if HIA_G_UPDATE
			    double adj=parms->powfs[ipowfs].dither_gog*mgold*(1.-pd->gg0->p[ig]);
			    simu->gradscale->p[iwfs]->p[ig]+=adj;
			    while(simu->gradscale->p[iwfs]->p[ig]<0){
				adj*=0.5;
				simu->gradscale->p[iwfs]->p[ig]+=-adj;
			    }
#else
			    if(pd->gg0->p[ig]>0.01){//skip weakly determined subapertures.
				simu->gradscale->p[iwfs]->p[ig]*=pow(pd->gg0->p[ig], -parms->powfs[ipowfs].dither_gog);
			    }
#endif
			}
			mgnew=dsum(simu->gradscale->p[iwfs])/ng;
			dzero(pd->gg0);
		    }
		    info("Step %5d wfs %d CoG gain adjusted from %g to %g %s.\n", 
			 simu->isim, iwfs, mgold, mgnew, ogtype);
		    if(simu->resdither){
			int ic=(npllacc-1)/(npll);
			P(simu->resdither->p[iwfs], 3, ic)=mgnew;
		    }			   
		    //adjust WFS measurement dither dithersig by gain adjustment. used for dither t/t removal from gradients.
		    //dbg("a2me=%g, mgold=%g, mgnew=%g\n", pd->a2me, mgold, mgnew);
		    pd->a2me*=(mgnew/mgold);
		    dcellscale(powfs[ipowfs].sanea, pow(mgnew/mgold, 2));
		    if(parms->save.dither){
			writebin(simu->gradscale->p[iwfs], "wfs%d_gradscale_%d", iwfs, simu->isim);
		    }
		    /*if(parms->powfs[ipowfs].skip==2 && parms->recon.fnsphpsd){//TWFS. Update TWFS gain.
			simu->eptwfs=twfs_gain_optim(parms, simu->recon, powfs);
			}*/
		}
	    }
	    if(parms->powfs[ipowfs].phytype_sim != parms->powfs[ipowfs].phytype_sim2){
		parms->powfs[ipowfs].phytype_sim=parms->powfs[ipowfs].phytype_sim2;
		parms->powfs[ipowfs].phytype_recon=parms->powfs[ipowfs].phytype_sim;
		info("Step %5d: powfs %d changed to %s\n", simu->isim, ipowfs, 
		      parms->powfs[ipowfs].phytype_sim==1?"matched filter":"CoG");
	    }
	    if(parms->powfs[ipowfs].phytype_sim==1){//Matched filter
		if(parms->powfs[ipowfs].neareconfile || parms->powfs[ipowfs].phyusenea){
		    warning("Disable neareconfile and phyusenea\n");
		    parms->powfs[ipowfs].neareconfile=NULL;
		    parms->powfs[ipowfs].phyusenea=0;
		}
		parms->powfs[ipowfs].phytype_recon=1;//Make sure MF is used for reconstruction.
		genmtch(parms, powfs, ipowfs);
		if(parms->save.dither==1){
		    writebin(powfs[ipowfs].intstat->mtche, "powfs%d_mtche_%d", ipowfs, simu->isim);
		    writebin(powfs[ipowfs].intstat->i0sum, "powfs%d_i0sum_%d", ipowfs, simu->isim);
		    writebin(powfs[ipowfs].sanea, "powfs%d_sanea_%d", ipowfs, simu->isim);
		}
#if USE_CUDA
		if(parms->gpu.wfs){
		    info("Update matched filter in GPU\n");
		    gpu_wfsgrad_update_mtche(parms, powfs);
		}
#endif
	    

		if(!parms->powfs[ipowfs].lo && parms->recon.alg==0){//no need to update LSR.
		    simu->tomo_update=2;
		}
	    }
	}
    }
}
/**
   TWFS has output. Accumulate result to simu->gradoff. It is put in wfsgrad.c
   instead of recon.c to avoid race condition because it updates simu->gradoff.
 */
void wfsgrad_twfs_recon(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int itpowfs=parms->itpowfs;
    const int ntstep=(simu->isim-parms->powfs[itpowfs].step+1);
    if(ntstep>0 && ntstep%parms->powfs[itpowfs].dtrat==0){
	info("Step %5d: TWFS[%d] has output with gain %g\n", simu->isim, itpowfs, simu->eptwfs);
	dcell *Rmod=0;
	//Build radial mode error using closed loop TWFS measurements from this time step.
	dcellmm(&Rmod, simu->recon->RRtwfs, simu->gradcl, "nn", 1);
	memcpy(PCOL(simu->restwfs, ntstep/parms->powfs[itpowfs].dtrat-1),
	       Rmod->p[0]->p, sizeof(double)*simu->restwfs->nx);
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    if(parms->powfs[ipowfs].llt){
		dmm(&simu->gradoff->p[iwfs], 1, simu->recon->GRall->p[iwfs], Rmod->p[0], "nn", -simu->eptwfs); 
		
		if(parms->plot.run){
		    const int nsa=simu->powfs[ipowfs].saloc->nloc;
		    drawopd("Goffx",simu->powfs[ipowfs].saloc, simu->gradoff->p[iwfs]->p,
			    parms->dbg.draw_gmax->p,
			    "WFS Offset (x)","x (m)", "y (m)", "x %d",  iwfs);
		    drawopd("Goffy",simu->powfs[ipowfs].saloc, simu->gradoff->p[iwfs]->p+nsa,
			    parms->dbg.draw_gmax->p,
			    "WFS Offset (y)","x (m)", "y (m)", "y %d",  iwfs);
		}
	    }
	    if(parms->save.dither){
		writebin(simu->gradoff->p[iwfs], "wfs%d_gradoff_twfs_%d", iwfs, simu->isim);
	    }
	}
	if(parms->save.dither){
	    writebin(Rmod, "Rmod_%d", simu->isim);
	}
	dcellfree(Rmod);
    
	if(parms->recon.psd){
	    const int ntacc=ntstep/parms->powfs[itpowfs].dtrat;
	    const int dtrat=parms->recon.psddtrat_twfs;
	    if(ntacc % dtrat==0){//output
		dbg("Step %5d: TWFS output psd\n", simu->isim);
		dmat *ts=dsub(simu->restwfs, 0, 0, ntacc-dtrat, dtrat);
		dmat *tts=dtrans(ts);dfree(ts);
		const double dt=parms->sim.dt*parms->powfs[itpowfs].dtrat;	
		dmat *psd=psd1dt(tts, parms->recon.psdnseg, dt);
		dfree(tts);
		//Sum all the PSDs
		psd_sum(psd, 1);
		dmat *psdol=servo_rej2ol(psd, parms->sim.dt, parms->powfs[itpowfs].dtrat, simu->eptwfs, 0);
		writebin(psd, "psdcl_twfs_%d", ntacc);
		writebin(psdol, "psdol_twfs_%d", ntacc);
		dcell *coeff=servo_optim(psdol, parms->sim.dt, parms->powfs[itpowfs].dtrat, M_PI*0.25, 0, 1);
		const double g=0.5;
		simu->eptwfs=simu->eptwfs*(1-g)+coeff->p[0]->p[0]*g;
		info("Step %5d New gain (twfs): %.3f\n", simu->isim, simu->eptwfs);
		dfree(psdol);
		cellfree(coeff);
		dfree(psd);
	    }
	}
    }
}
/**
   Calls wfsgrad_iwfs() to computes WFS gradient in parallel.
   It also includes operations on Gradients before tomography.
*/
void wfsgrad(SIM_T *simu){
    double tk_start=myclockd();
    const PARMS_T *parms=simu->parms;
    if(parms->sim.idealfit || parms->sim.evlol || parms->sim.idealtomo) return;
    // call the task in parallel and wait for them to finish. It may be done in CPU or GPU.
    extern int PARALLEL;
    if(!PARALLEL || parms->tomo.ahst_idealngs==1 || !parms->gpu.wfs){
	CALL_THREAD(simu->wfsgrad_pre, 0);
    }//else: already called by sim.c
    CALL_THREAD(simu->wfsgrad_post, 0);
    wfsgrad_dither_post(simu);//must be before wfsgrad_lgsfocus because wfsgrad_lgsfocus runs zoom integrator.
    if(parms->itpowfs!=-1){
	wfsgrad_twfs_recon(simu);
    }
    if(parms->nlgspowfs){ 
	//high pass filter lgs focus to remove sodium range variation effect
	wfsgrad_lgsfocus(simu);
    }
    if(1+simu->isim==parms->sim.end){
#if USE_CUDA
	if(parms->gpu.wfs){
	    gpu_save_pistat(simu);
	}else
#endif
	    save_pistat(simu);
    }
    simu->tk_wfs=myclockd()-tk_start;
}
