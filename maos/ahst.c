/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

/*
  The 5 NGS mode for split tomography with 2DM
  I work in not normalized zernike space. So the result is in radians,
  multiply to 2R to get zernike mode.
  2010-07-23: Added tikholnov regularization to Wa. 
*/
#include "common.h"
#include "ahst.h"

/**
   \file maos/ahst.c Contains functions to setup NGS modes and reconstructor
   using AHST for one or two DMs.  Use parms->wfsr instead of parms->wfs for wfs
   information, which hands GLAO mode correctly.
   */

static TIC;

/**
   Compute number of ahst modes from number of DMs.
 */
static int ngsmod_nmod(int ndm, double hs){
    int nmod=0;
    if(ndm==1 || !isfinite(hs)){
	nmod=2;
    }else if(ndm>=2){
	nmod=5;
    }else{
	error("Invalid ndm: %d\n",ndm);
    }
    return nmod;
}

/**
   computes the cross-coupling of NGS modes in science field.
   MCC=(M'*Ha'*W*Ha*M) where M is ngsmod on DM, Ha is propagator from DM to
   science. W is weighting in science.
*/
static dcell* ngsmod_mcc(const PARMS_T *parms, RECON_T *recon, const APER_T *aper, const double *wt){
    NGSMOD_T *ngsmod=recon->ngsmod;
    const loc_t *plocs=aper->locs;
    double *x=plocs->locx;
    double *y=plocs->locy;
    int nloc=plocs->nloc;
    double *amp=aper->amp->p;
    const int nmod=ngsmod->nmod;
    dcell *mcc=cellnew(parms->evl.nevl,1);
    PDMAT(aper->mcc,aMCC);
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	mcc->p[ievl]=dnew(nmod,nmod);
	PDMAT(mcc->p[ievl],MCC);
	MCC[0][0]=aMCC[1][1];
	MCC[1][1]=aMCC[2][2];
	MCC[1][0]=MCC[0][1]=aMCC[1][2];
    }
    if(ngsmod->nmod>=5){
	double *mod[nmod];
	mod[0]=x;
	mod[1]=y;
	for(int imod=2; imod<nmod; imod++){
	    mod[imod]=malloc(nloc*sizeof(double));
	}
	/*dc component of the focus mod. subtract during evaluation. */
	/*this is not precisely R^2/2 due to obscuration */
	const double MCC_fcp=aper->fcp;

	const double ht=ngsmod->ht;
	const double scale=ngsmod->scale;
	const double scale1=1.-scale;

	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    PDMAT(mcc->p[ievl],MCC);
	    if(fabs(wt[ievl])<1.e-12) continue;
	    double thetax=parms->evl.thetax->p[ievl];
	    double thetay=parms->evl.thetay->p[ievl];
	    for(int iloc=0; iloc<nloc; iloc++){
		double xx=x[iloc]*x[iloc];
		double xy=x[iloc]*y[iloc];
		double yy=y[iloc]*y[iloc];
		/*remove piston in focus */
		if(parms->sim.ahstfocus){
		    mod[2][iloc]=/*mod[2] has no focus effect on science*/
			-2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]);
		}else{
		    mod[2][iloc]=scale1*(xx+yy-MCC_fcp)
			-2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]);
		}
		mod[3][iloc]=scale1*(xx-yy)
		    -2.*ht*scale*(thetax*x[iloc]-thetay*y[iloc]);
		mod[4][iloc]=scale1*(xy)
		    -ht*scale*(thetay*x[iloc]+thetax*y[iloc]);
		if(nmod>5){
		    mod[5][iloc]=(xx+yy-MCC_fcp);//for focus tracking
		}
	    }
	    
	    for(int jmod=0; jmod<nmod; jmod++){
		for(int imod=jmod; imod<nmod; imod++){
		    if(imod<2&&jmod<2) continue;
		    double tmp=dotdbl(mod[imod],mod[jmod],amp,nloc);
		    MCC[imod][jmod]=tmp;
		    if(imod!=jmod){
			MCC[jmod][imod]=MCC[imod][jmod];
		    }
		}
	    }
	}
	for(int imod=2; imod<nmod; imod++){
	    free(mod[imod]);
	}
    }
 
    return mcc;
}
/**
   Compute NGS mode aperture weighting using science field.  Wa=Ha'*W*Ha where
   Ha is from ALOC to PLOCS and W is the amplitude weighting in PLOCS when
   use_ploc==0. Otherwise, Ha is from ALOC to PLOC and W is the amplitude
   weighting in PLOC.  */
static dspcell *ngsmod_Wa(const PARMS_T *parms, RECON_T *recon, 
			 const APER_T *aper, int use_ploc){
    const double *wt=parms->evl.wt->p;
    const int ndm=parms->ndm;
    loc_t *loc;
    double *amp=NULL;
    if(use_ploc){
	loc=recon->floc;
	amp=calloc(loc->nloc,sizeof(double));
	prop_nongrid_bin(aper->locs,aper->amp->p, loc,NULL,amp,1,0,0,1);
	normalize_sum(amp,loc->nloc,1);
    }else{
	amp=aper->amp->p;
	loc=aper->locs;
    }
    dspcell *Wa=NULL;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(fabs(wt[ievl])<1.e-12) continue;
	double thetax=parms->evl.thetax->p[ievl];
	double thetay=parms->evl.thetay->p[ievl];
	dspcell *Hat=cellnew(ndm,1);
	for(int idm=0; idm<ndm; idm++){
	    double hc = parms->dm[idm].ht;
	    double displacex=thetax*hc;
	    double displacey=thetay*hc;
	    /*from DM to ploc (plocs) science beam */
	    Hat->p[idm]=mkhb(recon->aloc->p[idm], loc, NULL, displacex,displacey,1.,0,0);
	    dspmuldiag(Hat->p[idm], amp, wt[ievl]);
	}
	dcellmm(&Wa, Hat, Hat, "nt", 1);
	dspcellfree(Hat);
    }
    if(use_ploc){
	free(amp);
    }
    return Wa;
}
/**
   compute NGS mode removal Pngs from LGS commands using aperture
   weighting. Pngs=(MCC)^-1 (Hm'*W*Ha).

   2012-05-25: The NGS mode removal should be based on old five modes even if now focus on PS1 is merged with defocus mode
*/
static dcell* ngsmod_Pngs_Wa(const PARMS_T *parms, RECON_T *recon, 
		     const APER_T *aper, int use_ploc){

    NGSMOD_T *ngsmod=recon->ngsmod;
    const double ht=ngsmod->ht;
    const double scale=ngsmod->scale;
    const double scale1=1.-scale;
    const double *wt=parms->evl.wt->p;
    const int ndm=parms->ndm;
    const int nmod=ngsmod->nmod;
    loc_t *loc;
    double *x, *y;
    int nloc;
    double *amp=NULL;
    if(use_ploc){
	loc=recon->floc;
	amp=calloc(loc->nloc,sizeof(double));
	prop_nongrid_bin(aper->locs,aper->amp->p,loc,NULL,amp,1,0,0,1);
	normalize_sum(amp,loc->nloc,1);
    }else{
	amp=aper->amp->p;
	loc=aper->locs;
    }
    x=loc->locx;
    y=loc->locy;
    nloc=loc->nloc;

    dcell *modc=cellnew(1,1);/*W*Hm*M */
    modc->p[0]=dnew(nloc,nmod);
    PDMAT(modc->p[0],mod);
    for(int iloc=0; iloc<nloc; iloc++){
	mod[0][iloc]=x[iloc]*amp[iloc];
	mod[1][iloc]=y[iloc]*amp[iloc];
    }
    const double MCC_fcp=ngsmod->aper_fcp;
    /*dc component of the focus mod. subtract during evaluation. */
    /*this is not precisely R^2/2 due to obscuration */
    dcell *HatWHmt=cellnew(ndm,1);
    dsp *HatGround=NULL;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(fabs(wt[ievl])<1.e-12) continue;
	double thetax=parms->evl.thetax->p[ievl];
	double thetay=parms->evl.thetay->p[ievl];
	if(nmod>=5){
	    for(int iloc=0; iloc<nloc; iloc++){
		double xx=x[iloc]*x[iloc];
		double xy=x[iloc]*y[iloc];
		double yy=y[iloc]*y[iloc];
		/*remove piston in focus */
		if(parms->sim.ahstfocus){
		    mod[2][iloc]=amp[iloc]
			*(-2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]));
		}else{
		    mod[2][iloc]=amp[iloc]
			*(scale1*(xx+yy-MCC_fcp)
			  -2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]));
		}
		mod[3][iloc]=amp[iloc]
		    *(scale1*(xx-yy)
		      -2.*ht*scale*(thetax*x[iloc]-thetay*y[iloc]));
		mod[4][iloc]=amp[iloc]
		    *(scale1*(xy)
		      -ht*scale*(thetay*x[iloc]+thetax*y[iloc]));
		if(nmod>5){
		    mod[5][iloc]=amp[iloc]*(xx+yy-MCC_fcp);
		}
	    }
	}
	dspcell *Hat=cellnew(ndm,1);
	for(int idm=0; idm<ndm; idm++){
	    double hc = parms->dm[idm].ht;
	    double displacex=thetax*hc;
	    double displacey=thetay*hc;
	    if(parms->dm[idm].isground && HatGround){
		info("Reusing HatGround\n");
		Hat->p[idm]=dspref(HatGround);
	    }else{
		/*from DM to ploc (plocs) science beam */
		Hat->p[idm]=mkhb(recon->aloc->p[idm], loc, NULL, displacex,displacey,1.,0,0);
		if(parms->dm[idm].isground){
		    HatGround=dspref(Hat->p[idm]);
		}
	    }
	}
	dcellmm(&HatWHmt,Hat,modc,"nn",wt[ievl]);
	dspcellfree(Hat);
    }
    dspfree(HatGround);
    dcell *IMCC=cellnew(1,1);
    IMCC->p[0]=dref(ngsmod->IMCC);
    dcell *Pngs=NULL;
    dcellmm(&Pngs,IMCC,HatWHmt,"nt",1);
    dcellfree(IMCC);
    dcellfree(modc);
    if(use_ploc){
	free(amp);
    }
    dcellfree(HatWHmt);
    return Pngs;
}
/**
   compute tip/tilt mode removal from each DM commands using aperture
   weighting. Ptt=(MCC_TT)^-1 *(Hmtt * W * Ha)
*/
static dcell* ngsmod_Ptt_Wa(const PARMS_T *parms, RECON_T *recon, 
			    const APER_T *aper, int use_ploc){
    NGSMOD_T *ngsmod=recon->ngsmod;
    const double *wt=parms->evl.wt->p;
    const int ndm=parms->ndm;
    const int nmod=2;
    loc_t *loc;
    double *x, *y;
    int nloc;
    double *amp=NULL;
    if(use_ploc){
	loc=recon->floc;
	amp=calloc(loc->nloc,sizeof(double));
	prop_nongrid_bin(aper->locs, aper->amp->p,loc,NULL,amp,1,0,0,1);
	normalize_sum(amp,loc->nloc,1);
    }else{
	amp=aper->amp->p;
	loc=aper->locs;
    }
    x=loc->locx;
    y=loc->locy;
    nloc=loc->nloc;

    dcell *modc=cellnew(1,1);/*W*Hm*M */
    modc->p[0]=dnew(nloc,nmod);
    PDMAT(modc->p[0],mod);
    for(int iloc=0; iloc<nloc; iloc++){
	mod[0][iloc]=x[iloc]*amp[iloc];
	mod[1][iloc]=y[iloc]*amp[iloc];
    }
    dcell *HatWHmt=cellnew(ndm,1);
    dsp *HatGround=NULL;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(fabs(wt[ievl])<1.e-12) continue;
	double thetax=parms->evl.thetax->p[ievl];
	double thetay=parms->evl.thetay->p[ievl];
	dspcell *Hat=cellnew(ndm,1);
	for(int idm=0; idm<ndm; idm++){
	    double hc = parms->dm[idm].ht;
	    double displacex=thetax*hc;
	    double displacey=thetay*hc;
	    if(!parms->dm[idm].isground || !HatGround){
		/*from DM to ploc (plocs) science beam */
		Hat->p[idm]=mkhb(recon->aloc->p[idm], loc, NULL, displacex,displacey,1.,0,0);
		if(parms->dm[idm].isground){
		    HatGround=dspref(Hat->p[idm]);
		}
	    }else{
		Hat->p[idm]=dspref(HatGround);
	    }
	    dspmm(&HatWHmt->p[idm],Hat->p[idm],modc->p[0],"nn",wt[ievl]);
	}
	dspcellfree(Hat);
    }
    dcell *IMCC=cellnew(1,1);
    IMCC->p[0]=dref(ngsmod->IMCC_TT);
    dcell *Ptt=cellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	dmm(&(Ptt->p[idm]),0,IMCC->p[0],HatWHmt->p[idm],"nt",1);
    }
    dcellfree(IMCC);
    dcellfree(modc);
    if(use_ploc){
	free(amp);
    }
    dcellfree(HatWHmt);
    return Ptt;
}
/**
   DM modes for all the low order modes, defined on DM grid. It uses ngsmod2dm
   to define the modes*/
static dcell *ngsmod_m(const PARMS_T *parms, RECON_T *recon){
    NGSMOD_T *ngsmod=recon->ngsmod; 
    int ndm=parms->ndm;
    int nmod=ngsmod->nmod;
    dcell *M=cellnew(1,1);
    M->p[0]=dnew(nmod,1);
    dcell *mod=cellnew(ndm,1);
    dcell *dmt=cellnew(ndm,1);
    loc_t **aloc=recon->aloc->p;
    for(int idm=0; idm<ndm; idm++){
	dmt->p[idm]=dnew(aloc[idm]->nloc,1);
	mod->p[idm]=dnew(aloc[idm]->nloc,nmod);
    }
    for(int imod=0; imod<nmod; imod++){
	dcellzero(dmt);
	M->p[0]->p[imod]=1;
	ngsmod2dm(&dmt,recon,M,1.);
	M->p[0]->p[imod]=0;
	for(int idm=0; idm<ndm; idm++){
	    memcpy(mod->p[idm]->p+imod*aloc[idm]->nloc, dmt->p[idm]->p, 
		   sizeof(double)*aloc[idm]->nloc);
	}
    }
    dcellfree(M);
    dcellfree(dmt);
    return mod;
}

/**
   Compute NGS modes Ha*M in the science directon using ray tracing. Not used
*/
dcell *ngsmod_hm_accphi(const PARMS_T *parms, RECON_T *recon, const APER_T *aper){
    /*Build NGS mod in science direction using accphi */
    NGSMOD_T *ngsmod=recon->ngsmod;
    loc_t **aloc=recon->aloc->p;
    const int ndm=parms->ndm;
    dcell *dmt=cellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	dmt->p[idm]=dnew(aloc[idm]->nloc,1);
    }
    dcell *M=cellnew(1,1);
    const int nmod=ngsmod->nmod;
    M->p[0]=dnew(nmod,1);
    dcell *HMC=cellnew(parms->evl.nevl,nmod);/*fine with SCAO or GLAO */
    PDCELL(HMC,HM);
    int nloc=aper->locs->nloc;
    for(int imod=0; imod<nmod; imod++){
	dzero(M->p[0]);
	M->p[0]->p[imod]=1;
	dcellzero(dmt);
	ngsmod2dm(&dmt,recon,M,1.);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    HM[imod][ievl]=dnew(nloc,1);
	    for(int idm=0; idm<parms->ndm; idm++){
		double ht=parms->dm[idm].ht;
		double displacex=parms->evl.thetax->p[ievl]*ht;
		double displacey=parms->evl.thetay->p[ievl]*ht;
		prop_nongrid(aloc[idm], 
			     dmt->p[idm]->p,
			     aper->locs, NULL, HM[imod][ievl]->p,1,
			     displacex, displacey, 1.,0,0);
	    }
	}
    }
    return HMC;
}
/**
   Compute NGS modes Ha*M in the science directon using analytic formula. Not used
 */
dcell *ngsmod_hm_ana(const PARMS_T *parms, RECON_T *recon, const APER_T *aper){
    /*confirmed to agree with ngsmod_hm_accphi except DM artifacts */
    NGSMOD_T *ngsmod=recon->ngsmod;
    const double hs=ngsmod->hs;
    const double ht=ngsmod->ht;
    const double scale=pow(1.-ht/hs,-2);
    const double scale1=1.-scale;
    const int nmod=ngsmod->nmod;
    const double MCC_fcp=recon->ngsmod->aper_fcp;
    dcell *HMC=cellnew(parms->evl.nevl,nmod);
    PDCELL(HMC,HM);
    double *x=aper->locs->locx;
    double *y=aper->locs->locy;
    int nloc=aper->locs->nloc;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	double sx=parms->evl.thetax->p[ievl];
	double sy=parms->evl.thetay->p[ievl];
	for(int imod=0; imod<nmod; imod++){
	    HM[imod][ievl]=dnew(nloc,1);
	}
	if(nmod==2){
	    for(int iloc=0; iloc<nloc; iloc++){
		HM[0][ievl]->p[iloc]=x[iloc];
		HM[1][ievl]->p[iloc]=y[iloc];
	    }
	}else if(nmod>=5){
	    for(int iloc=0; iloc<nloc; iloc++){
		double xx=pow(x[iloc],2);
		double yy=pow(y[iloc],2);
		double xy=x[iloc]*y[iloc];
		HM[0][ievl]->p[iloc]=x[iloc];
		HM[1][ievl]->p[iloc]=y[iloc];
		if(parms->sim.ahstfocus){
		    HM[2][ievl]->p[iloc]=
			-2*ht*(sx*x[iloc]+sy*y[iloc])*scale;
		}else{
		    HM[2][ievl]->p[iloc]=(xx+yy-MCC_fcp)*scale1
			-2*ht*(sx*x[iloc]+sy*y[iloc])*scale;
		}
		HM[3][ievl]->p[iloc]=(xx-yy)*scale1-2*ht*(sx*x[iloc]-sy*y[iloc])*scale;
		HM[4][ievl]->p[iloc]=(xy)*scale1-(sx*y[iloc]+sy*x[iloc])*ht*scale;
		if(nmod>5){
		    HM[5][ievl]->p[iloc]=(xx+yy-MCC_fcp);
		}
	    }
	}else{
	    error("Invalid nmod\n");
	}
    }
    return HMC;
}
/**
   setup NGS modes and reconstructor using AHST for one or two DMs.
 */
void setup_ngsmod(const PARMS_T *parms, RECON_T *recon, 
		  const APER_T *aper, POWFS_T *powfs){
    ngsmod_free(recon->ngsmod);
    NGSMOD_T *ngsmod=recon->ngsmod=calloc(1, sizeof(NGSMOD_T));
    ngsmod->ahstfocus=parms->sim.ahstfocus;
    const int ndm=parms->ndm;	
    ngsmod->aper_fcp=aper->fcp;
    if(ndm==2 && fabs(parms->dm[0].ht)>1.e-10){
	error("Error configuration. First DM is not on ground\n");
    }
    double hs=NAN;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].lo || parms->powfs[ipowfs].skip) continue;
	if(isnan(hs)){
	    hs=parms->powfs[ipowfs].hs;
	}else{
	    if(isfinite(hs) || isfinite(parms->powfs[ipowfs].hs)){
		if(fabs(hs-parms->powfs[ipowfs].hs)>1000){
		    //High order GS at different altitude.
		    hs=INFINITY; break;
		}
	    }
	}
    }
    ngsmod->nmod=ngsmod_nmod(ndm, hs);
    if(isfinite(hs) && (parms->sim.mffocus || parms->tomo.ahst_idealngs || parms->sim.ahstfocus)){
	if(parms->ndm==1){
	    error("Not implemented\n");
	}
	ngsmod->nmod++;/*add a global focus mode.*/
    }
    info2("ngsmod nmod=%d, ahstfocus=%d\n", ngsmod->nmod, parms->sim.ahstfocus);
    ngsmod->hs=hs;
    if(ndm>1){
	ngsmod->ht=parms->dm[ndm-1].ht;//last DM.
    }else{
	ngsmod->ht=0;
    }
    ngsmod->scale=pow(1.-ngsmod->ht/ngsmod->hs,-2);
    /*modal cross coupling matrix along science*/
    ngsmod->MCCP=ngsmod_mcc(parms,recon,aper, parms->evl.wt->p);
    if(ngsmod->MCCP->nx==1){
	ngsmod->MCC=dref(ngsmod->MCCP->p[0]);
    }else{
	ngsmod->MCC=NULL;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    dadd(&ngsmod->MCC, 1, ngsmod->MCCP->p[ievl], parms->evl.wt->p[ievl]);
	}
    }
    if(parms->save.setup){
	writebin(recon->ngsmod->MCC, "ahst_MCC");
    }
    ngsmod->IMCC=dinvspd(ngsmod->MCC);
    PDMAT(ngsmod->MCC,MCC);
    ngsmod->IMCC_TT=dnew(2,2);
    ngsmod->IMCC_TT->p[0]=MCC[0][0];
    ngsmod->IMCC_TT->p[1]=MCC[0][1];
    ngsmod->IMCC_TT->p[2]=MCC[1][0];
    ngsmod->IMCC_TT->p[3]=MCC[1][1];
    dinvspd_inplace(ngsmod->IMCC_TT);
    /*the ngsmodes defined on the DM.*/
    ngsmod->Modes=ngsmod_m(parms,recon);
    if(recon->actstuck){
	warning2("Apply stuck actuators to ngs modes\n");
	act_zero(recon->aloc, recon->ngsmod->Modes, recon->actstuck);
    }
    /*if(recon->actfloat){
      We do extrapolation to float actuators, so no need to modify Pngs/Ptt.
      warning2("Apply float actuators to Pngs, Ptt\n");
      act_zero(recon->aloc, recon->ngsmod->Modes, recon->actfloat);
      }*/
    /*
       W is recon->saneai;
       Rngs=(M'*G'*W*G*M)^-1*M'*G'*W
       Pngs=Rngs*GA
     */
    dspcell *saneai=recon->saneai;
    if(parms->recon.split==1 && !parms->sim.skysim && parms->ntipowfs){
	/*we disabled GA for low order wfs in skysim mode. */
	dcellmm(&ngsmod->GM, recon->GAlo, ngsmod->Modes, "nn", 1);
	if(parms->nlowfs==1 && ngsmod->nmod>5){
	    /*There is only one wfs, remove first plate scale mode*/
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].lo){
		    int nx=ngsmod->GM->p[iwfs]->nx;
		    memset(ngsmod->GM->p[iwfs]->p+nx*2, 0, nx*sizeof(double));
		}
	    }
	}
	ngsmod->Rngs=dcellpinv(ngsmod->GM,saneai);
    }
    if(parms->tomo.ahst_wt==1){
	/*Use gradient weighting. */
	//dcellmulsp(&ngsmod->Pngs, ngsmod->Rngs, recon->GAlo, 1);
	dcellmm(&ngsmod->Pngs, ngsmod->Rngs, recon->GAlo, "nn", 1);
    }else if(parms->tomo.ahst_wt==2){
	/*Use science based weighting. */
	if(parms->dbg.wamethod==0){
	    info("Wa using DM mode\n");
	    tic;
	    ngsmod->Wa=ngsmod_Wa(parms,recon,aper,0);
	    /*
	      Add tikhonov regularization. H is from aloc to some other loc. 
	      the eigen value of H'*amp*H is about 4/aloc->nloc.
	    */
	    int nact=0;
	    for(int idm=0; idm<parms->ndm; idm++){
		nact+=recon->aloc->p[idm]->nloc;
	    }
	    double maxeig=4./nact;
	    dcelladdI(ngsmod->Wa, 1e-9*maxeig);
	    
	    toc("Wa");
	    ngsmod->Pngs=dcellpinv(ngsmod->Modes,ngsmod->Wa);
	    toc("Pngs");
	}else{
	    info("Wa using science mode\n");
	    tic;
	    ngsmod->Pngs=ngsmod_Pngs_Wa(parms,recon,aper,0);
	    toc("Pngs_Wa");
	}
    }else if(parms->tomo.ahst_wt==3){/*Identity weighting. */
	ngsmod->Pngs=dcellpinv(ngsmod->Modes, NULL);
    }else{
	error("Invalid parms->tomo.ahst_wt=%d\n", parms->tomo.ahst_wt);
    }
   
 
    if(parms->save.setup){
	/*ahst stands for ad hoc split tomography */
    	writebin(recon->ngsmod->GM,  "ahst_GM");
	writebin(recon->ngsmod->Rngs,"ahst_Rngs");
	writebin(recon->ngsmod->Pngs,"ahst_Pngs");
	writebin(recon->ngsmod->Modes, "ahst_Modes");
	writebin(recon->ngsmod->Wa, "ahst_Wa");
    }
    
}
/**
   used in performance evaluation on science opds. accumulate to out*/
void calc_ngsmod_dot(double *pttr_out, double *pttrcoeff_out,
		     double *ngsmod_out,
		     const PARMS_T *parms,
		     const RECON_T *recon, const APER_T *aper, 
		     const double *opd, int ievl){

    const double *restrict amp=aper->amp->p;
    const double *restrict locx=aper->locs->locx;
    const double *restrict locy=aper->locs->locy;
    const double MCC_fcp=recon->ngsmod->aper_fcp;
    const double ht=recon->ngsmod->ht;
    const double scale=recon->ngsmod->scale;
    const double scale1=1.-scale;
    double coeff[6]={0,0,0,0,0,0};
    double tot=0;
    const double thetax=parms->evl.thetax->p[ievl]; 
    const double thetay=parms->evl.thetay->p[ievl]; 
    if(recon->ngsmod->nmod==2){
	for(int iloc=0; iloc<aper->locs->nloc; iloc++){
	    const double junk=amp[iloc]*opd[iloc];
	    tot+=junk*opd[iloc];
	    const double junkx=locx[iloc]*junk;
	    const double junky=locy[iloc]*junk;
	    coeff[0]+=junk;
	    coeff[1]+=junkx;
	    coeff[2]+=junky;
	}
    }else if(recon->ngsmod->nmod>=5){
	for(int iloc=0; iloc<aper->locs->nloc; iloc++){
	    const double junk=amp[iloc]*opd[iloc];
	    tot+=junk*opd[iloc];
	    coeff[0]+=junk;
	    const double junkx=locx[iloc]*junk;
	    const double junky=locy[iloc]*junk;
	    coeff[1]+=junkx;
	    coeff[2]+=junky;
	    coeff[3]+=locx[iloc]*junkx;//xx
	    coeff[4]+=locy[iloc]*junky;//yy
	    coeff[5]+=locx[iloc]*junky;//xy
	}
    }else{
	error("Invalid nmod\n");
    }
    if(pttrcoeff_out){
	memset(pttrcoeff_out, 0, sizeof(double)*3);
	dmulvec(pttrcoeff_out, aper->imcc, coeff, 1);
    }
    if(pttr_out){
	/*compute TT removed wavefront variance as a side product */
	double pis=aper->ipcc*coeff[0]*coeff[0];
	double ptt=dwdot3(coeff, aper->imcc, coeff);
	pttr_out[0]=tot-pis;/*PR */
	pttr_out[1]=ptt-pis;/*TT */
	pttr_out[2]=tot-ptt;/*PTTR */
    }
    /*don't use +=. need locking */
    ngsmod_out[0]=coeff[1];
    ngsmod_out[1]=coeff[2];
    if(recon->ngsmod->nmod>=5){
	if(parms->sim.ahstfocus){
	    ngsmod_out[2]=(-2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
	}else{
	    ngsmod_out[2]=(scale1*(coeff[3]+coeff[4]-coeff[0]*MCC_fcp)
			   -2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
	}
	ngsmod_out[3]=(scale1*(coeff[3]-coeff[4])
		       -2*scale*ht*(thetax*coeff[1]-thetay*coeff[2]));
	ngsmod_out[4]=(scale1*(coeff[5])
		       -scale*ht*(thetay*coeff[1]+thetax*coeff[2]));
	if(recon->ngsmod->nmod>5){
	    ngsmod_out[5]=(coeff[3]+coeff[4]-coeff[0]*MCC_fcp);
	}
    }
}
/**
   Convert NGS modes to DM actuator commands using analytical expression. For >2 DMs, we only put NGS modes on ground and top-most DM.
*/
void ngsmod2dm(dcell **dmc, const RECON_T *recon, const dcell *M, double gain){
    if(!M || !M->p[0]) return;
    assert(M->nx==1 && M->ny==1);
    double scale=recon->ngsmod->scale;
    /*The MCC_fcp depends weakly on the aperture sampling. */
    double MCC_fcp=recon->ngsmod->aper_fcp;
    loc_t **aloc=recon->aloc->p;
    /*convert mode vector and add to dm commands */
    const int ndm=recon->ndm;
    if(!*dmc){
	*dmc=cellnew(ndm,1);
    }
    for(int idm=0; idm<ndm; idm++){
	if(!(*dmc)->p[idm]){
	    (*dmc)->p[idm]=dnew(recon->aloc->p[idm]->nloc, 1);
	}
    }

    /*first dm */
    double *pm=M->p[0]->p;
    if(recon->ngsmod->nmod==2){
	if(M->p[0]->nx!=2) error("Invalid mode\n");
	int idm=0;
	double *p=(*dmc)->p[idm]->p;
	double *xloc=aloc[idm]->locx;
	double *yloc=aloc[idm]->locy;
	unsigned long nloc=aloc[idm]->nloc;
	for(unsigned long iloc=0; iloc<nloc; iloc++){
	    p[iloc]+=gain*(pm[0]*xloc[iloc]+pm[1]*yloc[iloc]);
	}
    }else if(recon->ngsmod->nmod>=5){
	double scale2=-scale*gain;
	for(int idm=0; idm<ndm; idm++){
	    double *p=(*dmc)->p[idm]->p;
	    unsigned long nloc=aloc[idm]->nloc;
	    double *xloc=aloc[idm]->locx;
	    double *yloc=aloc[idm]->locy;
	    if(idm==0){
		double pm_focus=pm[2];
		if(recon->ngsmod->nmod>5){
		    if(recon->ngsmod->ahstfocus){
			/*The net focus put in should have no focus effect on
			  science. global focus is controlled independently*/
			pm_focus=pm[2]*scale+pm[5];
		    }else{
			pm_focus=pm[2]+pm[5];
		    }
		}
		for(unsigned long iloc=0; iloc<nloc; iloc++){
		    double xx=xloc[iloc]*xloc[iloc];
		    double xy=xloc[iloc]*yloc[iloc];
		    double yy=yloc[iloc]*yloc[iloc];
		    p[iloc]+=gain*(pm[0]*xloc[iloc]
				   +pm[1]*yloc[iloc]
				   +pm_focus*(xx+yy-MCC_fcp)
				   +pm[3]*(xx-yy)
				   +pm[4]*(xy));
		}
	    }else if(idm+1==ndm){
		for(unsigned long iloc=0; iloc<nloc; iloc++){
		    double xx=xloc[iloc]*xloc[iloc];
		    double xy=xloc[iloc]*yloc[iloc];
		    double yy=yloc[iloc]*yloc[iloc];
		    p[iloc]+=scale2*(pm[2]*(xx+yy-MCC_fcp)
				     +pm[3]*(xx-yy)
				     +pm[4]*(xy));
		}
	    }	
	}
    }
}
/**
   Convert NGS mode vector to aperture grid for science directions.  */
void ngsmod2science(dmat *iopd, loc_t *loc, const NGSMOD_T *ngsmod, 
		    double thetax, double thetay,
		    const double *mod, double alpha){
    const double *locx=loc->locx;
    const double *locy=loc->locy;
    const int nmod=ngsmod->nmod;
    if(nmod==2){
	for(int iloc=0; iloc<loc->nloc; iloc++){
	    double tmp=locx[iloc]*mod[0]+locy[iloc]*mod[1];
	    iopd->p[iloc]+=tmp*alpha;
	}
    }else{
	const double ht=ngsmod->ht;
	const double scale=ngsmod->scale;
	const double scale1=1.-scale;
	double focus;
	if(nmod>5){
	    warning_once("Check accuracy with ahstfocus=0,1\n");
	    focus=mod[5];
	    if(!ngsmod->ahstfocus){
		focus+=mod[2]*scale1;
	    }
	}else{
	    focus=mod[2]*scale1;
	}
	for(int iloc=0; iloc<loc->nloc; iloc++){
	    double x=locx[iloc];
	    double y=locy[iloc];
	    double xy=x*y;
	    double x2=x*x;
	    double y2=y*y;
	    double tmp= locx[iloc]*mod[0]
		+locy[iloc]*mod[1]
		+focus*(x2+y2)
		+mod[2]*(-2*scale*ht*(thetax*x+thetay*y))
		+mod[3]*((x2-y2)*scale1 - 2*scale*ht*(thetax*x-thetay*y))
		+mod[4]*(xy*scale1-scale*ht*(thetay*x+thetax*y));
	    iopd->p[iloc]+=tmp*alpha;
	}
    }
}
void ngsmod_free(NGSMOD_T *ngsmod){
    if(!ngsmod) return;
    dcellfree(ngsmod->GM);
    dcellfree(ngsmod->Rngs);
    dcellfree(ngsmod->Pngs);
    dcellfree(ngsmod->Modes);
    dfree(ngsmod->MCC);
    dcellfree(ngsmod->MCCP);
    dfree(ngsmod->IMCC);
    dfree(ngsmod->IMCC_TT);
    free(ngsmod);
}

/**
   remove NGS modes from LGS DM commands
   if ahst_wt==1
   Rngs*GA*dmerr is zero
   if ahst_wt==2
   Doesn't perturb NGS modes in science direction.
   if ahst_wt==3
   Identity weighting.

   if nmod==6: make sure the global focus mode is not removed from LGS result.
*/
void remove_dm_ngsmod(SIM_T *simu, dcell *dmerr){
    if(!dmerr) return;
    const RECON_T *recon=simu->recon;
    dcellzero(simu->Mngs);
    dcellmm(&simu->Mngs, recon->ngsmod->Pngs, dmerr, "nn",1);
    //zero out global focus mode if any.
    if(recon->ngsmod->nmod>5){
	if(simu->parms->sim.ahstfocus){
	    const double scale=recon->ngsmod->scale;
	    simu->Mngs->p[0]->p[5]=simu->Mngs->p[0]->p[2]*(1-scale);
	}else{
	    simu->Mngs->p[0]->p[5]=0;
	}
    }
    dcellmm(&dmerr, recon->ngsmod->Modes, simu->Mngs, "nn", -1);
    //ngsmod2dm(&dmerr,recon, simu->Mngs,-1);
}
