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

/*
  The 5 NGS mode for split tomography with 2DM
  I work in not normalized zernike space. So the result is in radians,
  multiply to 2R to get zernike mode.
  2010-07-23: Added tikholnov regularization to Wa. 
*/
#include "maos.h"
#include "ahst.h"

/**
   \file maos/ahst.c Contains functions to setup NGS modes and reconstructor
   using AHST for one or two DMs.  Use parms->wfsr instead of parms->wfs for wfs
   information, which hands GLAO mode correctly.
   */

TIC;
/**
   Compute number of ahst modes from number of DMs.
 */
static int ngsmod_nmod(int ndm){
    int nmod=0;
    if(ndm==1)
	nmod=2;
    else if(ndm==2)
	nmod=5;
    else
	error("Invalid ndm: %d\n",ndm);
    return nmod;
}

/**
   computes the cross-coupling of NGS modes in science field.
   MCC=(M'*Ha'*W*Ha*M) where M is ngsmod on DM, Ha is propagator from DM to
   science. W is weighting in science.
*/
static dcell* ngsmod_mcc(const PARMS_T *parms, RECON_T *recon, APER_T *aper, const double *wt){

    NGSMOD_T *ngsmod=recon->ngsmod;
    double *x, *y;
    int nloc;
    double *amp=NULL;
  
    const loc_t *plocs=aper->locs;
    x=plocs->locx;
    y=plocs->locy;
    nloc=plocs->nloc;
    amp=aper->amp->p;
    
    dcell *mcc=NULL;
    if(parms->ndm==1){
	//Single conjugate. low order WFS only controls Tip/tilt
	mcc=dcellnew(1,1);
	mcc->p[0]=dnew(2,2);
	PDMAT(mcc->p[0],MCC);
	PDMAT(aper->mcc,aMCC);
	//not yet available in ngsmod
	MCC[0][0]=aMCC[1][1];
	MCC[1][1]=aMCC[2][2];
	MCC[1][0]=MCC[0][1]=aMCC[1][2];
    }else if(parms->ndm==2){
	double *mod[5];
	mcc=dcellnew(parms->evl.nevl, 1);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    mcc->p[ievl]=dnew(5,5);
	    PDMAT(mcc->p[ievl],MCC);
	    PDMAT(aper->mcc,aMCC);
	    MCC[0][0]=aMCC[1][1];
	    MCC[1][1]=aMCC[2][2];
	    MCC[1][0]=MCC[0][1]=aMCC[1][2];
	}
	mod[0]=x;
	mod[1]=y;
	mod[2]=malloc(nloc*sizeof(double));
	mod[3]=malloc(nloc*sizeof(double));
	mod[4]=malloc(nloc*sizeof(double));
	const double MCC_fcp=aper->fcp;
	//dc component of the focus mod. subtract during evaluation.
	//this is not precisely R^2/2 due to obscuration
	const double ht=ngsmod->ht;
	const double scale=ngsmod->scale;
	const double scale1=1.-scale;

	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    PDMAT(mcc->p[ievl],MCC);
	    if(fabs(wt[ievl])<1.e-12) continue;
	    double thetax=parms->evl.thetax[ievl];
	    double thetay=parms->evl.thetay[ievl];
	
	    for(int iloc=0; iloc<nloc; iloc++){
		double xx=x[iloc]*x[iloc];
		double xy=x[iloc]*y[iloc];
		double yy=y[iloc]*y[iloc];
		//remove piston in focus
		mod[2][iloc]=scale1*(xx+yy-MCC_fcp)
		    -2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]);
		mod[3][iloc]=scale1*(xx-yy)
		    -2.*ht*scale*(thetax*x[iloc]-thetay*y[iloc]);
		mod[4][iloc]=scale1*(xy)
		    -ht*scale*(thetay*x[iloc]+thetax*y[iloc]);
	    }
	    for(int jmod=0; jmod<5; jmod++){
		for(int imod=jmod; imod<5; imod++){
		    if(imod<2&&jmod<2) continue;
		    double tmp=dotdbl(mod[imod],mod[jmod],amp,nloc);
		    MCC[imod][jmod]=tmp;
		    if(imod!=jmod){
			MCC[jmod][imod]=MCC[imod][jmod];
		    }
		}
	    }
	}
	for(int imod=2; imod<5; imod++){
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
static spcell *ngsmod_Wa(const PARMS_T *parms, RECON_T *recon, 
			 APER_T *aper, int use_ploc){
    const double *wt=parms->evl.wt;
    const int ndm=parms->ndm;
    loc_t *loc;
    double *amp=NULL;
    if(use_ploc){
	loc=recon->floc;
	amp=calloc(loc->nloc,sizeof(double));
	prop_nongrid_bin(aper->locs,aper->amp->p, loc,NULL,amp,1,0,0,1);
	normalize(amp,loc->nloc,1);
    }else{
	amp=aper->amp->p;
	loc=aper->locs;
    }
    spcell *Wa=NULL;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(fabs(wt[ievl])<1.e-12) continue;
	double thetax=parms->evl.thetax[ievl];
	double thetay=parms->evl.thetay[ievl];

	spcell *Hat=spcellnew(ndm,1);
	spcell *Ha=spcellnew(1,ndm);
	for(int idm=0; idm<ndm; idm++){
	    double hc = parms->dm[idm].ht;
	    double displacex=thetax*hc;
	    double displacey=thetay*hc;
	    //from DM to ploc (plocs) science beam
	    Hat->p[idm]=mkhb(recon->aloc[idm], loc, NULL, displacex,displacey,1.,0,0);
	    Ha->p[idm]=sptrans(Hat->p[idm]);
	    spmuldiag(Hat->p[idm], amp, wt[ievl]);
	}
	spcell *HatHai=spcellmulspcell(Hat,Ha,1);
	spcelladd(&Wa, HatHai);
	/*
	for(int idm=0; idm<ndm*ndm; idm++){
	    spfull(&(Wa->p[idm]), HatHai->p[idm], 1);
	    }*/
	spcellfree(Hat);
	spcellfree(Ha);
	spcellfree(HatHai);
    }
    if(use_ploc){
	free(amp);
    }
    return Wa;
}
/**
   compute NGS mode removal Pngs from LGS commands using aperture weighting. Pngs=(MCC)^-1 (Hm'*W*Ha)

*/
static dcell* ngsmod_Pngs_Wa(const PARMS_T *parms, RECON_T *recon, 
		     APER_T *aper, int use_ploc){

    NGSMOD_T *ngsmod=recon->ngsmod;
    const double ht=ngsmod->ht;
    const double scale=ngsmod->scale;
    const double scale1=1.-scale;
    const double *wt=parms->evl.wt;
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
	normalize(amp,loc->nloc,1);
    }else{
	amp=aper->amp->p;
	loc=aper->locs;
    }
    x=loc->locx;
    y=loc->locy;
    nloc=loc->nloc;

    dcell *modc=dcellnew(1,1);//W*Hm*M
    modc->p[0]=dnew(nloc,nmod);
    PDMAT(modc->p[0],mod);
    for(int iloc=0; iloc<nloc; iloc++){
	mod[0][iloc]=x[iloc]*amp[iloc];
	mod[1][iloc]=y[iloc]*amp[iloc];
    }
    const double MCC_fcp=ngsmod->aper_fcp;
    //dc component of the focus mod. subtract during evaluation.
    //this is not precisely R^2/2 due to obscuration
    dcell *HatWHmt=dcellnew(ndm,1);
    dsp *HatGround=NULL;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(fabs(wt[ievl])<1.e-12) continue;
	double thetax=parms->evl.thetax[ievl];
	double thetay=parms->evl.thetay[ievl];
	if(nmod==5){
	    for(int iloc=0; iloc<nloc; iloc++){
		double xx=x[iloc]*x[iloc];
		double xy=x[iloc]*y[iloc];
		double yy=y[iloc]*y[iloc];
		//remove piston in focus
		mod[2][iloc]=amp[iloc]
		    *(scale1*(xx+yy-MCC_fcp)
		      -2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]));
		mod[3][iloc]=amp[iloc]
		    *(scale1*(xx-yy)
		      -2.*ht*scale*(thetax*x[iloc]-thetay*y[iloc]));
		mod[4][iloc]=amp[iloc]
		    *(scale1*(xy)
		      -ht*scale*(thetay*x[iloc]+thetax*y[iloc]));
	    }
	}
	spcell *Hat=spcellnew(ndm,1);
	for(int idm=0; idm<ndm; idm++){
	    double hc = parms->dm[idm].ht;
	    double displacex=thetax*hc;
	    double displacey=thetay*hc;
	    if(parms->dm[idm].isground && HatGround){
		info("Reusing HatGround\n");
		Hat->p[idm]=spref(HatGround);
	    }else{
		//from DM to ploc (plocs) science beam
		Hat->p[idm]=mkhb(recon->aloc[idm], loc, NULL, displacex,displacey,1.,0,0);
		if(parms->dm[idm].isground){
		    HatGround=spref(Hat->p[idm]);
		}
	    }
	}
	spcellmulmat(&HatWHmt,Hat,modc,wt[ievl]);
	spcellfree(Hat);
    }
    spfree(HatGround);
    dcell *IMCC=dcellnew(1,1);
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
			    APER_T *aper, int use_ploc){
    NGSMOD_T *ngsmod=recon->ngsmod;
    const double *wt=parms->evl.wt;
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
	normalize(amp,loc->nloc,1);
    }else{
	amp=aper->amp->p;
	loc=aper->locs;
    }
    x=loc->locx;
    y=loc->locy;
    nloc=loc->nloc;

    dcell *modc=dcellnew(1,1);//W*Hm*M
    modc->p[0]=dnew(nloc,nmod);
    PDMAT(modc->p[0],mod);
    for(int iloc=0; iloc<nloc; iloc++){
	mod[0][iloc]=x[iloc]*amp[iloc];
	mod[1][iloc]=y[iloc]*amp[iloc];
    }
    dcell *HatWHmt=dcellnew(ndm,1);
    dsp *HatGround=NULL;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(fabs(wt[ievl])<1.e-12) continue;
	double thetax=parms->evl.thetax[ievl];
	double thetay=parms->evl.thetay[ievl];
	spcell *Hat=spcellnew(ndm,1);
	for(int idm=0; idm<ndm; idm++){
	    double hc = parms->dm[idm].ht;
	    double displacex=thetax*hc;
	    double displacey=thetay*hc;
	    if(!parms->dm[idm].isground || !HatGround){
		//from DM to ploc (plocs) science beam
		Hat->p[idm]=mkhb(recon->aloc[idm], loc, NULL, displacex,displacey,1.,0,0);
		if(parms->dm[idm].isground){
		    HatGround=spref(Hat->p[idm]);
		}
	    }else{
		Hat->p[idm]=spref(HatGround);
	    }
	    spmulmat(&HatWHmt->p[idm],Hat->p[idm],modc->p[0],wt[ievl]);
	}
	spcellfree(Hat);
    }
    dcell *IMCC=dcellnew(1,1);
    IMCC->p[0]=dref(ngsmod->IMCC_TT);
    dcell *Ptt=dcellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	dmm(&(Ptt->p[idm]),IMCC->p[0],HatWHmt->p[idm],"nt",1);
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
   DM modes for all the low order modes.*/
static dcell *ngsmod_m(const PARMS_T *parms, RECON_T *recon){
    NGSMOD_T *ngsmod=recon->ngsmod; 
    int ndm=parms->ndm;
    int nmod=ngsmod->nmod;
    dcell *M=dcellnew(1,1);
    M->p[0]=dnew(nmod,1);
    dcell *mod=dcellnew(ndm,1);
    dcell *dmt=dcellnew(ndm,1);
    loc_t **aloc=recon->aloc;
    for(int idm=0; idm<ndm; idm++){
	dmt->p[idm]=dnew(aloc[idm]->nloc,1);
	mod->p[idm]=dnew(aloc[idm]->nloc,nmod);
    }
    for(int imod=0; imod<nmod; imod++){
	dcellzero(M);
	M->p[0]->p[imod]=1;
	dcellzero(dmt);
	ngsmod2dm(&dmt,recon,M,1.);
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
   from ngsmod to ztilt/gtilt gradients using GA. this GM depends on aperture
   sampling before somehow GA is sensitive to piston due to partial
   illumination. GP has no gradient on piston.

*/

static dcell *ngsmod_g(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    NGSMOD_T *ngsmod=recon->ngsmod;
    int ndm=parms->ndm;
    loc_t **aloc=recon->aloc;
    int nmod=ngsmod->nmod;
    dcell *ZSN=dcellnew(parms->nwfsr,1);
    //NGS mode vector
    dcell *M=dcellnew(1,1);
    M->p[0]=dnew(nmod,1);
    //DM command sfor NGS mod
    dcell *dmt=dcellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	dmt->p[idm]=dnew(aloc[idm]->nloc,1);
    }
    PSPCELL(recon->GA,GA);
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].skip)
	    continue;

	int nsa=powfs[ipowfs].pts->nsa;
	if(ZSN->p[iwfs]) continue;
	ZSN->p[iwfs]=dnew(nsa*2,nmod);
	PDMAT(ZSN->p[iwfs],grad);
	dmat *grad2=calloc(1, sizeof(dmat));
	grad2->nx=nsa*2; 
	grad2->ny=1;
	
	for(int imod=0; imod<nmod; imod++){
	    grad2->p=grad[imod];
	    dcellzero(M);
	    M->p[0]->p[imod]=1;
	    dcellzero(dmt);
	    //convert NGS mode vector to DM commands
	    ngsmod2dm(&dmt,recon,M,1.);
	    dzero(grad2);
	    for(int idm=0; idm<parms->ndm; idm++){
		spmulmat(&grad2,GA[idm][iwfs],dmt->p[idm],1);
	    }
	}
	free(grad2);
    }
    dcellfree(M);
    dcellfree(dmt);
    return ZSN;
}
/**
   Compute NGS modes Ha*M in the science directon using ray tracing. Not used
*/
dcell *ngsmod_hm_accphi(const PARMS_T *parms, RECON_T *recon, APER_T *aper){
    //Build NGS mod in science direction using accphi
    NGSMOD_T *ngsmod=recon->ngsmod;
    loc_t **aloc=recon->aloc;
    const int ndm=parms->ndm;
    dcell *dmt=dcellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	dmt->p[idm]=dnew(aloc[idm]->nloc,1);
    }
    dcell *M=dcellnew(1,1);
    const int nmod=ngsmod->nmod;
    M->p[0]=dnew(nmod,1);
    dcell *HMC=dcellnew(parms->evl.nevl,nmod);//fine with SCAO or GLAO
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
		double displacex=parms->evl.thetax[ievl]*ht;
		double displacey=parms->evl.thetay[ievl]*ht;
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
dcell *ngsmod_hm_ana(const PARMS_T *parms, RECON_T *recon, 
			    APER_T *aper){
    //if(parms->tomo.split!=1) error("Only in split mode 1.");
    //confirmed to agree with ngsmod_hm_accphi except DM artifacts
    NGSMOD_T *ngsmod=recon->ngsmod;
    const double hs=ngsmod->hs;
    const double ht=ngsmod->ht;
    const double scale=pow(1.-ht/hs,-2);
    const double scale1=1.-scale;
    const int nmod=ngsmod->nmod;
    const double MCC_fcp=recon->ngsmod->aper_fcp;
    dcell *HMC=dcellnew(parms->evl.nevl,nmod);
    PDCELL(HMC,HM);
    double *x=aper->locs->locx;
    double *y=aper->locs->locy;
    int nloc=aper->locs->nloc;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	double sx=parms->evl.thetax[ievl];
	double sy=parms->evl.thetay[ievl];
	for(int imod=0; imod<nmod; imod++){
	    HM[imod][ievl]=dnew(nloc,1);
	}
	if(nmod==2){
	    for(int iloc=0; iloc<nloc; iloc++){
		HM[0][ievl]->p[iloc]=x[iloc];
		HM[1][ievl]->p[iloc]=y[iloc];
	    }
	}else if(nmod==5){
	    for(int iloc=0; iloc<nloc; iloc++){
		double xx=pow(x[iloc],2);
		double yy=pow(y[iloc],2);
		double xy=x[iloc]*y[iloc];
		HM[0][ievl]->p[iloc]=x[iloc];
		HM[1][ievl]->p[iloc]=y[iloc];
		HM[2][ievl]->p[iloc]=(xx+yy-MCC_fcp)*scale1
		    -2*ht*(sx*x[iloc]+sy*y[iloc])*scale;
		HM[3][ievl]->p[iloc]=(xx-yy)*scale1-2*ht*(sx*x[iloc]-sy*y[iloc])*scale;
		HM[4][ievl]->p[iloc]=(xy)*scale1-(sx*y[iloc]+sy*x[iloc])*ht*scale;
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
		  APER_T *aper, POWFS_T *powfs){
    //if(parms->tomo.split!=1) error("Only work in split mode 1.");

    NGSMOD_T *ngsmod=recon->ngsmod=calloc(1, sizeof(NGSMOD_T));
    const int ndm=parms->ndm;	
    ngsmod->aper_fcp=aper->fcp;
    ngsmod->nmod=ngsmod_nmod(ndm);
    if(ndm==2 && fabs(parms->dm[0].ht)>1.e-10){
	error("Error configuration. First DM is not on ground\n");
    }
    double hs=INFINITY;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(!isinf(parms->powfs[ipowfs].hs)){
	    if(!isinf(hs) && fabs(hs-parms->powfs[ipowfs].hs)>1.e-5) 
		error("There are multiple LGS type with different hs.\n");
	    hs=parms->powfs[ipowfs].hs;
	}
    }
    if(isinf(hs)){
	warning("No LGS found. No Need split tomography\n");
    }
    ngsmod->hs=hs;
    if(ndm>1){
	ngsmod->ht=parms->dm[1].ht;
    }else{
	ngsmod->ht=0;
    }
    ngsmod->scale=pow(1.-ngsmod->ht/ngsmod->hs,-2);
    ngsmod->MCCP=ngsmod_mcc(parms,recon,aper, parms->evl.wt);
    if(ngsmod->MCCP->nx==1){
	ngsmod->MCC=dref(ngsmod->MCCP->p[0]);
    }else{
	ngsmod->MCC=NULL;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    dadd(&ngsmod->MCC, 1, ngsmod->MCCP->p[ievl], parms->evl.wt[ievl]);
	}
    }
    ngsmod->IMCC=dinvspd(ngsmod->MCC);
    PDMAT(ngsmod->MCC,MCC);
    ngsmod->IMCC_TT=dnew(2,2);
    ngsmod->IMCC_TT->p[0]=MCC[0][0];
    ngsmod->IMCC_TT->p[1]=MCC[0][1];
    ngsmod->IMCC_TT->p[2]=MCC[1][0];
    ngsmod->IMCC_TT->p[3]=MCC[1][1];
    dinvspd_inplace(ngsmod->IMCC_TT);
    if(parms->save.setup){
	dwrite(recon->ngsmod->MCC, "%s/ahst_MCC", dirsetup);
    }
  
    ngsmod->Modes=ngsmod_m(parms,recon);
    /*
       W is recon->saneai;
       Rngs=(M'*G'*W*G*M)^-1*M'*G'*W
       Pngs=Rngs*GA
     */
    spcell *saneai=recon->saneai;
    if(parms->tomo.split==1 && !parms->sim.skysim){
	//we disabled GA for low order wfs in skysim mode.
	ngsmod->GM=ngsmod_g(parms,recon,powfs);
	ngsmod->Rngs=dcellpinv(ngsmod->GM,NULL,saneai);
    }
    if(parms->tomo.ahst_wt==1){
	//Use gradient weighting.
	dcellmulsp(&ngsmod->Pngs, ngsmod->Rngs, recon->GAlo, 1);
	if(parms->tomo.ahst_rtt){
	    ngsmod->Ptt=dcellnew(parms->ndm,1);
	    for(int idm=0; idm<ndm; idm++){
		dcell *Matt=dcellnew(ndm,1);
		dcell *GaM=NULL;
		Matt->p[idm]=loc2mat(recon->aloc[idm],0);
		spcellmulmat(&GaM, recon->GAlo, Matt, 1);
		dcell *tmp=dcellpinv(GaM, NULL,saneai);
		dcell *tmp2=NULL;
		dcellmulsp(&tmp2,tmp,recon->GAlo, 1);
		ngsmod->Ptt->p[idm]=dref(tmp2->p[idm]);
		dcellfree(tmp);
		dcellfree(tmp2);
		dcellfree(GaM);
		dcellfree(Matt);
	    }
	}
    }else if(parms->tomo.ahst_wt==2){
	//Use aperture weighting.
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
		nact+=recon->aloc[idm]->nloc;
	    }
	    double maxeig=4./nact;
	    spcelladdI(ngsmod->Wa, 1e-9*maxeig);
	    
	    toc("Wa");
	    ngsmod->Pngs=dcellpinv(ngsmod->Modes, NULL,ngsmod->Wa);
	    toc("Pngs");
	    if(parms->tomo.ahst_rtt){
		ngsmod->Ptt=dcellnew(parms->ndm,1);
		for(int idm=0; idm<ndm; idm++){
		    dmat *Matt=loc2mat(recon->aloc[idm],0);
		    ngsmod->Ptt->p[idm]=dpinv(Matt, NULL, ngsmod->Wa->p[idm+idm*ndm]);
		    dfree(Matt);
		}
	    }
	}else{
	    info("Wa using science mode\n");
	    tic;
	    ngsmod->Pngs=ngsmod_Pngs_Wa(parms,recon,aper,0);
	    toc("Pngs_Wa");
	    if(parms->tomo.ahst_rtt){
		tic;
		ngsmod->Ptt=ngsmod_Ptt_Wa(parms,recon,aper,0);
		toc("Pngs_Ptt");
	    }
	}
    }else if(parms->tomo.ahst_wt==3){//Identity weighting.
	ngsmod->Pngs=dcellpinv(ngsmod->Modes, NULL,NULL);
	if(parms->tomo.ahst_rtt){
	    ngsmod->Ptt=dcellnew(parms->ndm,1);
	    for(int idm=0; idm<ndm; idm++){
		dmat *Matt=loc2mat(recon->aloc[idm],0);
		ngsmod->Ptt->p[idm]=dpinv(Matt, NULL,NULL);
		dfree(Matt);
	    }
	}
    }else{
	error("Invalid parms->tomo.ahst_wt=%d\n", parms->tomo.ahst_wt);
    }
    if(recon->actstuck){
	warning2("Apply stuck actuators to Pngs, Ptt\n");
	act_stuck(recon->aloc, NULL, recon->ngsmod->Pngs, recon->actstuck);
	act_stuck(recon->aloc, NULL, recon->ngsmod->Ptt, recon->actstuck);
	act_zero(recon->aloc, recon->ngsmod->Modes, recon->actstuck);
    }
    if(recon->actfloat){
	warning2("Apply float actuators to Pngs, Ptt\n");
	act_float(recon->aloc, NULL, recon->ngsmod->Pngs, recon->actfloat);
	act_float(recon->aloc, NULL, recon->ngsmod->Ptt, recon->actfloat);
	act_zero(recon->aloc, recon->ngsmod->Modes, recon->actfloat);
    }
    if(parms->save.setup){
	//ahst stands for ad hoc split tomography
    	dcellwrite(recon->ngsmod->GM,  "%s/ahst_GM",  dirsetup);
	dcellwrite(recon->ngsmod->Rngs,"%s/ahst_Rngs",dirsetup);
	dcellwrite(recon->ngsmod->Pngs,"%s/ahst_Pngs",dirsetup);
	dcellwrite(recon->ngsmod->Ptt, "%s/ahst_Ptt", dirsetup);
	dcellwrite(recon->ngsmod->Modes, "%s/ahst_Modes", dirsetup);
	spcellwrite(recon->ngsmod->Wa, "%s/ahst_Wa", dirsetup);
    }
  
}
/**
   used in performance evaluation on science opds. accumulate to out*/
void calc_ngsmod_dot(double *pttr_out, double *pttrcoeff_out,
		     double *ngsmod_out,
		     const PARMS_T *parms,
		     const RECON_T *recon, const APER_T *aper, 
		     const double *opd, int ievl){

    const double *amp=aper->amp->p;
    const double *locx=aper->locs->locx;
    const double *locy=aper->locs->locy;
    const double MCC_fcp=recon->ngsmod->aper_fcp;
    const double ht=recon->ngsmod->ht;
    const double scale=recon->ngsmod->scale;
    const double scale1=1.-scale;
    double coeff[6]={0,0,0,0,0,0};
    double tot=0;
    const double thetax=parms->evl.thetax[ievl]; 
    const double thetay=parms->evl.thetay[ievl]; 
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
    }else if(recon->ngsmod->nmod==5){
	for(int iloc=0; iloc<aper->locs->nloc; iloc++){
	    const double junk=amp[iloc]*opd[iloc];
	    tot+=junk*opd[iloc];
	    const double junkx=locx[iloc]*junk;
	    const double junky=locy[iloc]*junk;
	    coeff[0]+=junk;
	    coeff[1]+=junkx;
	    coeff[2]+=junky;
	    coeff[3]+=locx[iloc]*junkx;
	    coeff[4]+=locy[iloc]*junky;
	    coeff[5]+=locx[iloc]*junky;
	}
    }else{
	error("Invalid nmod\n");
    }
    if(pttrcoeff_out){
	memset(pttrcoeff_out, 0, sizeof(double)*3);
	dmulvec(pttrcoeff_out, aper->imcc, coeff, 1);
    }
    if(pttr_out){
	//compute TT removed wavefront variance as a side product
	double pis=aper->ipcc*coeff[0]*coeff[0];
	double ptt=dwdot3(coeff, aper->imcc, coeff);
	pttr_out[0]=tot-pis;//PR
	pttr_out[1]=ptt-pis;//TT
	pttr_out[2]=tot-ptt;//PTTR
    }
    //don't use +=. need locking
    ngsmod_out[0]=coeff[1];
    ngsmod_out[1]=coeff[2];
    if(recon->ngsmod->nmod==5){
	ngsmod_out[2]=(scale1*(coeff[3]+coeff[4]-coeff[0]*MCC_fcp)
		       -2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
	ngsmod_out[3]=(scale1*(coeff[3]-coeff[4])
		       -2*scale*ht*(thetax*coeff[1]-thetay*coeff[2]));
	ngsmod_out[4]=(scale1*(coeff[5])
		       -scale*ht*(thetay*coeff[1]+thetax*coeff[2]));
    }
}
/**
   Convert NGS modes to DM actuator commands using analytical expression.
   \todo switch to use mode vectors
*/
void ngsmod2dm(dcell **dmc, const RECON_T *recon, const dcell *M, double gain){
    if(!M || !M->p[0]) return;
    assert(M->nx==1 && M->ny==1);
    double scale=recon->ngsmod->scale;
    //The MCC_fcp depends weakly on the aperture sampling.
    double MCC_fcp=recon->ngsmod->aper_fcp;
    loc_t **aloc=recon->aloc;
    //convert mode vector and add to dm commands
    const int ndm=recon->ndm;
    if(!*dmc){
	*dmc=dcellnew(ndm,1);
    }
    for(int idm=0; idm<ndm; idm++){
	if(!(*dmc)->p[idm]){
	    (*dmc)->p[idm]=dnew(recon->aloc[idm]->nloc, 1);
	}
    }

    if(ndm>2) error("Error Usage\n");
    //first dm
    double *pm=M->p[0]->p;
    if(ndm==1){
	if(M->p[0]->nx!=2) error("Invalid mode\n");
	int idm=0;
	double *p=(*dmc)->p[idm]->p;
	double *xloc=aloc[idm]->locx;
	double *yloc=aloc[idm]->locy;
	unsigned long nloc=aloc[idm]->nloc;
	for(unsigned long iloc=0; iloc<nloc; iloc++){
	    p[iloc]+=gain*(pm[0]*xloc[iloc]+pm[1]*yloc[iloc]);
	}
    }else if(ndm==2){
	if(M->p[0]->nx!=5) error("Invalid mode\n");
	double scale2=-scale*gain;
	for(int idm=0; idm<ndm; idm++){
	    double *p=(*dmc)->p[idm]->p;
	    unsigned long nloc=aloc[idm]->nloc;
	    double *xloc=aloc[idm]->locx;
	    double *yloc=aloc[idm]->locy;
	    if(idm==0){
		for(unsigned long iloc=0; iloc<nloc; iloc++){
		    double xx=xloc[iloc]*xloc[iloc];
		    double xy=xloc[iloc]*yloc[iloc];
		    double yy=yloc[iloc]*yloc[iloc];
		    p[iloc]+=gain*(pm[0]*xloc[iloc]
				   +pm[1]*yloc[iloc]
				   +pm[2]*(xx+yy-MCC_fcp)
				   +pm[3]*(xx-yy)
				   +pm[4]*(xy));
		}
	    }else if(idm==1){
		for(unsigned long iloc=0; iloc<nloc; iloc++){
		    double xx=xloc[iloc]*xloc[iloc];
		    double xy=xloc[iloc]*yloc[iloc];
		    double yy=yloc[iloc]*yloc[iloc];
		    p[iloc]+=scale2*(pm[2]*(xx+yy-MCC_fcp)
				     +pm[3]*(xx-yy)
				     +pm[4]*(xy));
		}
	    }else{
		error("Invalid\n");
	    }	
	}
    }
}
/**
   Convert NGS mode vector to aperture grid for science directions.  */
void ngsmod2science(dmat *iopdevl, const PARMS_T *parms,
		    const RECON_T *recon, const APER_T *aper,
		    const dcell *M, int ievl, double gain){
    const double *locx=aper->locs->locx;
    const double *locy=aper->locs->locy;
    const double *mod=M->p[0]->p;
    if(recon->ngsmod->nmod==2){
	for(int iloc=0; iloc<aper->locs->nloc; iloc++){
	    double tmp=locx[iloc]*mod[0]+locy[iloc]*mod[1];
	    iopdevl->p[iloc]+=tmp*gain;
	}
    }else{
	const double MCC_fcp=recon->ngsmod->aper_fcp;
	const double ht=recon->ngsmod->ht;
	const double scale=recon->ngsmod->scale;
	const double scale1=1.-scale;
	const double thetax=parms->evl.thetax[ievl]; 
	const double thetay=parms->evl.thetay[ievl];
	for(int iloc=0; iloc<aper->locs->nloc; iloc++){
	    double x=locx[iloc];
	    double y=locy[iloc];
	    double xy=x*y;
	    double x2=x*x;
	    double y2=y*y;
	    double tmp= locx[iloc]*mod[0]
		+locy[iloc]*mod[1]
		+mod[2]*((x2+y2-MCC_fcp)*scale1-2*scale*ht*(thetax*x+thetay*y))
		+mod[3]*((x2-y2)*scale1 - 2*scale*ht*(thetax*x-thetay*y))
		+mod[4]*(xy*scale1-scale*ht*(thetay*x+thetax*y));
	    iopdevl->p[iloc]+=tmp*gain;
	}
    }
}

/**
   remove NGS modes from LGS DM commands
   if ahst_wt==1
   Rngs*GA*dmerr is zero
   if ahst_wt==2
   Doesn't perturb NGS modes in science direction.
*/
void remove_dm_ngsmod(SIM_T *simu, dcell *dmerr){
    const RECON_T *recon=simu->recon;
    dcell *Mngs=NULL;
    dcellmm(&Mngs, recon->ngsmod->Pngs, dmerr, "nn",1);
    ngsmod2dm(&dmerr,recon, Mngs,-1);
    dcellfree(Mngs);
}
/**
   Removal tip/tilt on invidual DMs. Be careful about the roll off near the
   edge.  */
void remove_dm_tt(SIM_T *simu, dcell *dmerr){
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
