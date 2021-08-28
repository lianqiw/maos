/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/**
 * \file ahst.h
   Contains functions to setup NGS modes and reconstructor
   using AHST for one or more DMs.  Use parms->wfsr instead of parms->wfs for wfs
   information, which hands GLAO mode correctly.

   Notice that update of this file may require GPU code update accordingly
*/

/*
  The 5 NGS mode for split tomography with 2DM
  I work in not normalized zernike space. So the result is in radians,
  multiply to 2R to get zernike mode.
  2010-07-23: Added tikholnov regularization to Wa.
*/
#include "common.h"
#include "ahst.h"
#include "pywfs.h"
/*
   2017-09-11: When there is misregistration/distortion between the DM and
   science pupil, the assumed NGS mode on DM remain intact, but the influence on
   Science OPD needs to use ray tracing.
*/
static void ngsmod2dm(dcell** dmc, const recon_t* recon, const dcell* M, real gain);
static TIC;

/**
   computes the cross-coupling of NGS modes in science field.
   MCC=(M'*Ha'*W*Ha*M) where M is ngsmod on DM, Ha is propagator from DM to
   science. W is weighting in science.
*/
static dcell* ngsmod_mcc(const parms_t* parms, recon_t* recon, const aper_t* aper, const real* wt){
	ngsmod_t* ngsmod=recon->ngsmod;
	const loc_t* plocs=aper->locs;
	real* x=plocs->locx;
	real* y=plocs->locy;
	int nloc=plocs->nloc;
	real* amp=aper->amp->p;
	const int nmod=ngsmod->nmod;
	dcell* mcc=dcellnew(parms->evl.nevl, 1);
	dmat* aMCC=aper->mcc;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		dmat* MCC=P(mcc,ievl)=dnew(nmod, nmod);
		P(MCC, 0, 0)=P(aMCC, 1, 1);
		P(MCC, 1, 1)=P(aMCC, 2, 2);
		P(MCC, 0, 1)=P(MCC, 1, 0)=P(aMCC, 2, 1);
	}

	if(ngsmod->nmod>2){
		tic;
		real* mod[nmod];
		mod[0]=x;
		mod[1]=y;
		for(int imod=2; imod<nmod; imod++){
			mod[imod]=mymalloc(nloc, real);
		}
		/*dc component of the focus mod. subtract during evaluation. */
		/*this is not precisely R^2/2 due to obscuration */

		const real MCC_fcp=aper->fcp;
		const real ht=ngsmod->ht;
		const real scale=ngsmod->scale;
		const real scale1=1.-scale;
		//dmat *modvec=dnew(nmod, 1);
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			dmat* MCC=P(mcc,ievl);
			if(fabs(wt[ievl])<1.e-12) continue;
			real thetax=P(parms->evl.thetax,ievl);
			real thetay=P(parms->evl.thetay,ievl);

			for(int iloc=0; iloc<nloc; iloc++){
				real xx=x[iloc]*x[iloc];
				real xy=x[iloc]*y[iloc];
				real yy=y[iloc]*y[iloc];
				//remove piston in focus 
				if(ngsmod->indps){
					if(ngsmod->ahstfocus){
						mod[ngsmod->indps][iloc]=//mod[2] has no focus effect on science
							-2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]);
					} else{
						mod[ngsmod->indps][iloc]=scale1*(xx+yy-MCC_fcp)
							-2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]);
					}
					mod[ngsmod->indps+1][iloc]=scale1*(xx-yy)
						-2.*ht*scale*(thetax*x[iloc]-thetay*y[iloc]);
					mod[ngsmod->indps+2][iloc]=scale1*(xy)
						-ht*scale*(thetay*x[iloc]+thetax*y[iloc]);
				}
				if(ngsmod->indastig){
					mod[ngsmod->indastig][iloc]=(xx-yy);
					mod[ngsmod->indastig+1][iloc]=(xy);
				}
				if(ngsmod->indfocus){
					mod[ngsmod->indfocus][iloc]=(xx+yy-MCC_fcp);//for focus tracking
				}
			}

			for(int jmod=0; jmod<nmod; jmod++){
				for(int imod=jmod; imod<nmod; imod++){
					if(imod<2&&jmod<2) continue;
					real tmp=dvecdot(mod[imod], mod[jmod], amp, nloc);
					P(MCC, imod, jmod)=P(MCC, jmod, imod)=tmp;
				}
			}
		}//for(ievl)
		//dfree(modvec);
		for(int imod=2; imod<nmod; imod++){
			free(mod[imod]);
		}
		toc2("mcc");
	}

	return mcc;
}
/**
   Compute NGS mode aperture weighting using science field.  Wa=Ha'*W*Ha where
   Ha is from ALOC to PLOCS and W is the amplitude weighting in PLOCS when
   use_ploc==0. Otherwise, Ha is from ALOC to PLOC and W is the amplitude
   weighting in PLOC.  */
static dspcell* ngsmod_Wa(const parms_t* parms, recon_t* recon,
	const aper_t* aper, int use_ploc){
	const real* wt=parms->evl.wt->p;
	const int ndm=parms->ndm;
	loc_t* loc;
	real* amp=NULL;
	if(use_ploc){
		loc=recon->floc;
		amp=mycalloc(loc->nloc, real);
		prop_nongrid_bin(aper->locs, aper->amp->p, loc, amp, 1, 0, 0, 1);
		dnormalize_sumabs(amp, loc->nloc, 1);
	} else{
		amp=aper->amp->p;
		loc=aper->locs;
	}
	dspcell* Wa=NULL;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(fabs(wt[ievl])<1.e-12) continue;
		real thetax=P(parms->evl.thetax,ievl);
		real thetay=P(parms->evl.thetay,ievl);
		dspcell* Hat=dspcellnew(ndm, 1);
		for(int idm=0; idm<ndm; idm++){
			real hc=parms->dm[idm].ht;
			real displacex=thetax*hc;
			real displacey=thetay*hc;
			/*from DM to ploc (plocs) science beam */
			P(Hat,idm)=mkhb(P(recon->aloc,idm), loc, displacex, displacey, 1.);
			dspmuldiag(P(Hat,idm), amp, wt[ievl]);
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
static dcell* ngsmod_Pngs_Wa(const parms_t* parms, recon_t* recon,
	const aper_t* aper, int use_ploc){

	ngsmod_t* ngsmod=recon->ngsmod;
	const real ht=ngsmod->ht;
	const real scale=ngsmod->scale;
	const real scale1=1.-scale;
	const real* wt=parms->evl.wt->p;
	const int ndm=parms->ndm;
	const int nmod=ngsmod->nmod;
	loc_t* loc;
	real* x, * y;
	int nloc;
	real* amp=NULL;
	if(use_ploc){
		loc=recon->floc;
		amp=mycalloc(loc->nloc, real);
		prop_nongrid_bin(aper->locs, aper->amp->p, loc, amp, 1, 0, 0, 1);
		dnormalize_sumabs(amp, loc->nloc, 1);
	} else{
		amp=aper->amp->p;
		loc=aper->locs;
	}
	x=loc->locx;
	y=loc->locy;
	nloc=loc->nloc;

	dcell* WHm=dcellnew(1, 1);/* W_amp*Hm */
	P(WHm,0)=dnew(nloc, nmod);
	dmat* mod=P(WHm,0);
	for(int iloc=0; iloc<nloc; iloc++){
		P(mod, iloc, 0)=x[iloc]*amp[iloc];
		P(mod, iloc, 1)=y[iloc]*amp[iloc];
	}
	const real MCC_fcp=ngsmod->aper_fcp;
	/*dc component of the focus mod. subtract during evaluation. */
	/*this is not precisely R^2/2 due to obscuration */
	dcell* HatWHm=dcellnew(ndm, 1);
	dsp* HatGround=NULL;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(fabs(wt[ievl])<1.e-12) continue;
		real thetax=P(parms->evl.thetax,ievl);
		real thetay=P(parms->evl.thetay,ievl);
		if(nmod>2){
			for(int iloc=0; iloc<nloc; iloc++){
				real xx=x[iloc]*x[iloc];
				real xy=x[iloc]*y[iloc];
				real yy=y[iloc]*y[iloc];
				if(ngsmod->indps){
					if(ngsmod->ahstfocus){//no focus mode at ps
						P(mod, iloc, ngsmod->indps)=amp[iloc]
							*(-2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]));
					} else{
						P(mod, iloc, ngsmod->indps)=amp[iloc]
							*(scale1*(xx+yy-MCC_fcp)
								-2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]));
					}
					P(mod, iloc, ngsmod->indps+1)=amp[iloc]
						*(scale1*(xx-yy)
							-2.*ht*scale*(thetax*x[iloc]-thetay*y[iloc]));
					P(mod, iloc, ngsmod->indps+2)=amp[iloc]
						*(scale1*(xy)
							-ht*scale*(thetay*x[iloc]+thetax*y[iloc]));
				}
				if(ngsmod->indastig){
					P(mod, iloc, ngsmod->indastig)=amp[iloc]*(xx-yy);
					P(mod, iloc, ngsmod->indastig+1)=amp[iloc]*(xy);
				}
				if(ngsmod->indfocus){ /*remove piston in focus */
					P(mod, iloc, ngsmod->indfocus)=amp[iloc]*(xx+yy-MCC_fcp);
				}
			}
		}
		dspcell* Hat=dspcellnew(ndm, 1);
		for(int idm=0; idm<ndm; idm++){
			real hc=parms->dm[idm].ht;
			real displacex=thetax*hc;
			real displacey=thetay*hc;
			if(parms->dm[idm].isground&&HatGround){
				P(Hat,idm)=dspref(HatGround);
			} else{
			/*from DM to ploc (plocs) science beam */
				P(Hat,idm)=mkhb(P(recon->aloc,idm), loc, displacex, displacey, 1.);
				if(parms->dm[idm].isground){
					HatGround=dspref(P(Hat,idm));
				}
			}
		}
		dcellmm(&HatWHm, Hat, WHm, "nn", wt[ievl]);
		dspcellfree(Hat);
	}
	dspfree(HatGround);
	dcell* IMCC=dcellnew(1, 1);
	P(IMCC,0)=dref(ngsmod->IMCC);
	dcell* Pngs=NULL;
	dcellmm(&Pngs, IMCC, HatWHm, "nt", 1);
	dcellfree(IMCC);
	dcellfree(WHm);
	if(use_ploc){
		free(amp);
	}
	dcellfree(HatWHm);
	return Pngs;
}

/**
   DM modes for all the low order modes, defined on DM grid. It uses ngsmod2dm
   to define the modes*/

static dcell* ngsmod_dm(const parms_t* parms, recon_t* recon){
	ngsmod_t* ngsmod=recon->ngsmod;
	int ndm=parms->ndm;
	int nmod=ngsmod->nmod;
	dcell* M=dcellnew(1, 1);
	P(M,0)=dnew(nmod, 1);
	dcell* mod=dcellnew(ndm, 1);
	dcell* dmt=dcellnew(ndm, 1);
	loc_t** aloc=recon->aloc->p;
	for(int idm=0; idm<ndm; idm++){
		P(dmt,idm)=dnew(aloc[idm]->nloc, 1);
		P(mod,idm)=dnew(aloc[idm]->nloc, nmod);
	}
	for(int imod=0; imod<nmod; imod++){
		dcellzero(dmt);
		P(P(M,0),imod)=1;
		ngsmod2dm(&dmt, recon, M, 1.);
		P(P(M,0),imod)=0;
		for(int idm=0; idm<ndm; idm++){
			real piston=dsum(P(dmt,idm))/aloc[idm]->nloc;
			for(long iact=0; iact<aloc[idm]->nloc; iact++){
				P(P(mod,idm), iact, imod)=P(P(dmt,idm), iact)-piston;
			}
		}
	}
	dcellfree(M);
	dcellfree(dmt);
	return mod;
}

/**
   AHST parameters that are related to the geometry only, and
   will not be updated when estimated WFS measurement noise changes.
*/
void setup_ngsmod_prep(const parms_t* parms, recon_t* recon,
	const aper_t* aper, const powfs_t* powfs){
	if(recon->ngsmod){
		warning("Should only be called once\n");
		return;
	}
	ngsmod_t* ngsmod=recon->ngsmod=mycalloc(1, ngsmod_t);
	ngsmod->ahstfocus=parms->tomo.ahst_focus;
	const int ndm=parms->ndm;
	ngsmod->aper_fcp=aper->fcp;
	if(ndm>1&&fabs(parms->dm[0].ht)>1.e-10){
		error("Error configuration. First DM is not on ground\n");
	}
	real hs=NAN;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].lo||parms->powfs[ipowfs].skip) continue;
		if(isnan(hs)){
			hs=parms->powfs[ipowfs].hs;
		} else{
			if(isfinite(hs)||isfinite(parms->powfs[ipowfs].hs)){
				if(fabs(hs-parms->powfs[ipowfs].hs)>1000){
					//High order GS at different altitude.
					hs=INFINITY; break;
				}
			}
		}
	}
	ngsmod->nmod=2; //Basic tip/tilt mode.
	if(isfinite(hs)){
	//LGS WFS.
		if(ndm>1){//Plate scale mode
			if(parms->evl.nevl>1){//cannot determine plate scale with only 1 direction.
				ngsmod->indps=ngsmod->nmod;
				ngsmod->nmod+=3;
			}
		} else if(parms->nhiwfs>1){//Astigmatism for LTAO
			ngsmod->indastig=ngsmod->nmod;
			ngsmod->nmod+=2;
		}
		{//Always enable focus when LGS WFS is present
			ngsmod->indfocus=ngsmod->nmod;
			ngsmod->nmod+=1;
		}
	}
	info("ahst: nmod=%d, mffocus=%d, ahst_focus=%d\n", ngsmod->nmod, parms->sim.mffocus, parms->tomo.ahst_focus);
	ngsmod->hs=hs;
	if(ndm>1){
		ngsmod->ht=parms->dm[ndm-1].ht;//last DM.
	} else{
		ngsmod->ht=0;
	}
	ngsmod->scale=pow(1.-ngsmod->ht/ngsmod->hs, -2);
	/*modal cross coupling matrix along science for performance evaluation. */
	ngsmod->MCCP=ngsmod_mcc(parms, recon, aper, parms->evl.wt->p);
	if(ngsmod->MCCP->nx==1){
		ngsmod->MCC=dref(P(ngsmod->MCCP,0));
	} else{
		ngsmod->MCC=NULL;
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			dadd(&ngsmod->MCC, 1, P(ngsmod->MCCP,ievl), P(parms->evl.wt,ievl));
		}
	}
	dmat* MCCl=dchol(ngsmod->MCC);
	ngsmod->MCCu=dtrans(MCCl); dfree(MCCl);
	if(parms->save.setup){
		writebin(recon->ngsmod->MCC, "ahst_MCC");
	}
	ngsmod->IMCC=dinvspd(ngsmod->MCC);
	dmat* MCC=ngsmod->MCC;
	ngsmod->IMCC_TT=dnew(2, 2);
	P(ngsmod->IMCC_TT,0)=P(MCC, 0, 0);
	P(ngsmod->IMCC_TT,1)=P(MCC, 1, 0);
	P(ngsmod->IMCC_TT,2)=P(MCC, 0, 1);
	P(ngsmod->IMCC_TT,3)=P(MCC, 1, 1);
	dinvspd_inplace(ngsmod->IMCC_TT);
	if(ngsmod->indfocus){
		ngsmod->IMCC_F=dnew(ngsmod->nmod, ngsmod->nmod);
		const int indfocus=ngsmod->indfocus;
		P(ngsmod->IMCC_F, indfocus, indfocus)=1./P(ngsmod->MCC, indfocus, indfocus);
	}
	/*the ngsmodes defined on the DM.*/
	ngsmod->Modes=ngsmod_dm(parms, recon);
	if(recon->actstuck&&!parms->recon.modal&&parms->dbg.recon_stuck){
		warning("Apply stuck actuators to ngs modes\n");
		act_zero(recon->aloc, recon->ngsmod->Modes, recon->actstuck);
	}
   /*if(recon->actfloat){
	  We do extrapolation to float actuators, so no need to modify Pngs/Ptt.
	  warning("Apply float actuators to Pngs, Ptt\n");
	  act_zero(recon->aloc, recon->ngsmod->Modes, recon->actfloat);
	  }*/

	if(parms->recon.modal){//convert Modes to space of amod
		for(int idm=0; idm<parms->ndm; idm++){
			dmat* proj=dpinv(P(recon->amod,idm), NULL);
			dmat* tmp=0;
			dmm(&tmp, 0, proj, P(ngsmod->Modes,idm), "nn", 1);
			dfree(P(ngsmod->Modes,idm));
			P(ngsmod->Modes,idm)=tmp;
			dfree(proj);
		}
	}
	/*
	  ngsmod to NGS gradient interaction matrix. Defined in modal space for modal control
	*/
	if(parms->recon.split==1&&!parms->tomo.ahst_idealngs&&parms->ntipowfs){
		ngsmod->GM=dcellnew(parms->nwfsr, 1);
		info2("Low order control includes WFS");
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->powfs[ipowfs].skip!=3 &&
				(parms->powfs[ipowfs].lo
				||(parms->recon.split && !parms->powfs[ipowfs].trs))){
				info2(" %d", iwfs);
				for(int idm=0; idm<parms->ndm; idm++){
					if(parms->powfs[ipowfs].type==WFS_SH || parms->recon.modal){//shwfs or modal control
						dcellmm(&P(ngsmod->GM, iwfs), P(recon->GAlo, iwfs, idm), P(ngsmod->Modes, idm), "nn", 1);
					} else{//pwfs in zonal control.
						real  ht=parms->dm[idm].ht-parms->powfs[ipowfs].hc;
						real  dispx=0, dispy=0;
						dispx=parms->wfsr[iwfs].thetax*ht;
						dispy=parms->wfsr[iwfs].thetay*ht;
						dmat* tmp=pywfs_mkg(powfs[ipowfs].pywfs, P(recon->aloc,idm),
							parms->misreg.dm2wfs[iwfs+idm*parms->nwfs],
							P(ngsmod->Modes, idm), 0, dispx, dispy);
						dadd(&P(ngsmod->GM, iwfs), 1, tmp, 1);//accumulate
					}
				}
			}
		}
		info2("\n");
	}
	
	/**

	   Next, decouple LGS reconstruction from ngsmodes.

	   parms->tomo.ahst_wt control NGS modes removal from LGS DM commands
	   if ahst_wt==1
	   Rngs*GA*dmerr is zero
	   if ahst_wt==2
	   Doesn't perturb NGS modes in science direction.
	   if ahst_wt==3
	   Identity weighting.

	*/
	if(parms->tomo.ahst_wt==1){
	//Do it in setup_ngsmod_recon();
	} else if(parms->tomo.ahst_wt==2){
	/*Use science based weighting to isolate active meta pupil. */
		if(parms->dbg.wamethod==0){
			info("Wa using DM mode\n");

			tic;
			ngsmod->Wa=ngsmod_Wa(parms, recon, aper, 1);
			/*
			  Add tikhonov regularization. H is from aloc to some other loc.
			  the eigen value of H'*amp*H is about 4/aloc->nloc.
			*/
			int nact=0;
			for(int idm=0; idm<parms->ndm; idm++){
				nact+=P(recon->aloc,idm)->nloc;
			}
			real maxeig=4./nact;
			dcelladdI(ngsmod->Wa, 1e-9*maxeig);

			toc2("Wa");
			ngsmod->Pngs=dcellpinv(ngsmod->Modes, ngsmod->Wa);
			if(parms->save.setup){
				writebin(ngsmod->Wa, "ahst_Wa");
			}
			cellfree(ngsmod->Wa);
			toc2("Pngs");
		} else{
			dbg("Wa using science mode\n");
			tic;
			ngsmod->Pngs=ngsmod_Pngs_Wa(parms, recon, aper, 0);
			toc2("Pngs_Wa");
		}
	} else if(parms->tomo.ahst_wt==3){/*Identity weighting. */
		ngsmod->Pngs=dcellpinv(ngsmod->Modes, NULL);
	} else{
		error("Invalid parms->tomo.ahst_wt=%d\n", parms->tomo.ahst_wt);
	}
	if(parms->save.setup){
		writebin(ngsmod->Modes, "ahst_Modes");
		writebin(ngsmod->GM, "ahst_GM");
		if(ngsmod->Pngs) writebin(ngsmod->Pngs, "ahst_Pngs");
	}
}
/**
   Invert GM while handling rank deficiency.

   It ignores those WFS whoes mask is 0, and assumes the modes are ordered as
   below: T/T x, T/T y, PS1, PS2, PS3, Focus.
*/
static dcell* inv_gm(const dcell* GM, const dspcell* saneai, const lmat* mask, lmat** pmodvalid){
	if(GM->ny!=1){
		error("To be implemented\n");
	}
	dbg("Rngs is using wfs ");
	dcell* GM2=dcellnew(GM->nx, GM->ny);
	int nmod=0, ntt=0, nttf=0;
	for(int iwfs=0; iwfs<GM->nx; iwfs++){
		if((!mask||P(mask, iwfs))&&P(GM, iwfs)&&P(saneai, iwfs, iwfs)->px[0]>0){
			dbg0(" %d", iwfs);
			P(GM2, iwfs)=ddup(P(GM, iwfs));
			nmod=P(GM2, iwfs)->ny;
			int ng=P(GM2, iwfs)->nx;
			if(ng>=8){//TTF OIWFS
				nttf++;
			} else if(ng==2){
				ntt++;
			} else{
				error("Unknown WFS type: ng=%d\n", ng);
			}
		}
	}
	dbg("\n");
	lmat* modvalid=0;
	if(pmodvalid){
		if(!*pmodvalid){
			*pmodvalid=lnew(nmod, 1);
		} else{
			lset(*pmodvalid, 0);
		}
		modvalid=*pmodvalid;
	} else{
		modvalid=lnew(nmod, 1);
	}
	if(nttf>0||ntt>0){
		lset(modvalid, 1);
	}
	if(nttf>0){
		if(nttf==1&&ntt==0&&nmod==6){
			P(modvalid,2)=0;//disable magnofication mode
		}//else: all mode is valid
	} else{//nttf==0;
		if(nmod>=6){
			P(modvalid,5)=0;//no focus control.
		}
		if(nmod>=5){
			if(ntt<3){//1 or 2 TT OIWFS can not control all 3 PS modes
				P(modvalid,2)=0;
			}
			if(ntt<2){//1 TT OIWFS can not control any PS mode.
				P(modvalid,3)=0;
				P(modvalid,4)=0;
			}
			if(ntt==0){
				P(modvalid,0)=0;
				P(modvalid,1)=0;
			}
		}
	}

	for(int iwfs=0; iwfs<GM->nx; iwfs++){
		if(P(GM2, iwfs)){
			for(int imod=0; imod<nmod; imod++){
				if(!P(modvalid, imod)){
					int ng=P(GM2, iwfs)->nx;
					memset(PCOL(P(GM2, iwfs), imod), 0, ng*sizeof(real));
				}
			}
		}
	}
	dcell* RM=dcellpinv(GM2, saneai);
	dcellfree(GM2);
	if(!pmodvalid){
		lfree(modvalid);
	}
	return RM;
}
/**
   setup NGS modes reconstructor in ahst mode.
 */
void setup_ngsmod_recon(const parms_t* parms, recon_t* recon){
	ngsmod_t* ngsmod=recon->ngsmod;
	if(parms->recon.split==1&&!parms->tomo.ahst_idealngs&&parms->ntipowfs){
		cellfree(ngsmod->Rngs);
		ngsmod->Rngs=dccellnew(2, 1);
		P(ngsmod->Rngs,0)=inv_gm(ngsmod->GM, recon->saneai, 0, 0);
		/*
		  W is recon->saneai;
		  Rngs=(M'*G'*W*G*M)^-1*M'*G'*W
		  Pngs=Rngs*GA
		*/
		if(parms->sim.dtrat_lo!=parms->sim.dtrat_lo2){
			//Multi-rate control
			int nwfsr=parms->nwfsr;
			lmat* mask=lnew(nwfsr, 1);
			for(int iwfsr=0; iwfsr<nwfsr; iwfsr++){
				int ipowfs=parms->wfsr[iwfsr].powfs;
				if(parms->powfs[ipowfs].dtrat==parms->sim.dtrat_lo2&&P(ngsmod->GM, iwfsr)){
					P(mask, iwfsr)=1;
				}
			}
			P(ngsmod->Rngs,1)=inv_gm(ngsmod->GM, recon->saneai, mask, &ngsmod->modvalid);
			lfree(mask);
			switch(parms->dbg.lo_blend){
			case 0://Just use the two together
				ngsmod->lp2=0;
				break;
			case 1://Slower loop generates offset for faster loop
				ngsmod->lp2=-1;
				break;
			case 2://HPF on faster loop @ 1/20 of lower loop.
				ngsmod->lp2=fc2lp(0.05/parms->sim.dtrat_lo, parms->sim.dtrat_lo2);
				break;
			default:
				error("Invalid dbg.lo_blend=%d\n", parms->dbg.lo_blend);
			}
			warning("ngsmod->lp2=%g\n", ngsmod->lp2);
		}
	}

	if(parms->tomo.ahst_wt==1){
	/*Use gradient weighting. */
		dcellzero(ngsmod->Pngs);
		dcellmm(&ngsmod->Pngs, P(ngsmod->Rngs,0), recon->GAlo, "nn", 1);
		if(parms->save.setup){
			writebin(ngsmod->Pngs, "ahst_Pngs");
		}
	}

	if(parms->save.setup){
		writebin(ngsmod->Rngs, "ahst_Rngs");
	}
}
/**
   used in performance evaluation on science opds. accumulate to out*/
void calc_ngsmod_dot(real* pttr_out, real* pttrcoeff_out,
	real* ngsmod_out,
	const parms_t* parms, const ngsmod_t* ngsmod,
	const aper_t* aper, const real* opd, int ievl){
	const real* restrict amp=aper->amp->p;
	const real* restrict locx=aper->locs->locx;
	const real* restrict locy=aper->locs->locy;
	real coeff[6]={0,0,0,0,0,0};
	double tot=0; //use double for accumulation.
	const int nmod=ngsmod->nmod;
	if(nmod==2){
		for(int iloc=0; iloc<aper->locs->nloc; iloc++){
			const real junk=amp[iloc]*opd[iloc];
			tot+=junk*opd[iloc];
			const real junkx=locx[iloc]*junk;
			const real junky=locy[iloc]*junk;
			coeff[0]+=junk;
			coeff[1]+=junkx;
			coeff[2]+=junky;
		}
	} else{
		for(int iloc=0; iloc<aper->locs->nloc; iloc++){
			const real junk=amp[iloc]*opd[iloc];
			tot+=junk*opd[iloc];
			coeff[0]+=junk;//piston
			const real junkx=locx[iloc]*junk;
			const real junky=locy[iloc]*junk;
			coeff[1]+=junkx;//tip
			coeff[2]+=junky;//tilt
			coeff[3]+=locx[iloc]*junkx;//xx
			coeff[4]+=locy[iloc]*junky;//yy
			coeff[5]+=locx[iloc]*junky;//xy
		}
	}
	const real thetax=P(parms->evl.thetax,ievl);
	const real thetay=P(parms->evl.thetay,ievl);
	calc_ngsmod_post(pttr_out, pttrcoeff_out, ngsmod_out, tot, coeff, ngsmod, aper, thetax, thetay);
}
/**
   Separate post processing part so that GPU code can call it. Return non zero if error happens.
*/
int calc_ngsmod_post(real* pttr_out, real* pttrcoeff_out, real* ngsmod_out,
	real tot, const real* coeff, const ngsmod_t* ngsmod,
	const aper_t* aper, real thetax, real thetay){
	const real MCC_fcp=ngsmod->aper_fcp;
	const real ht=ngsmod->ht;
	const real scale=ngsmod->scale;
	const real scale1=1.-scale;
	int ans=0;
	if(pttrcoeff_out){
		memset(pttrcoeff_out, 0, sizeof(real)*3);
		dmulvec(pttrcoeff_out, aper->imcc, coeff, 1);
	}
	if(pttr_out){
	/*compute TT removed wavefront variance as a side product */
		real pis=aper->ipcc*coeff[0]*coeff[0];
		real ptt=dwdot(coeff, aper->imcc, coeff);
		pttr_out[0]=tot-pis;/*PR */
		pttr_out[1]=ptt-pis;/*TT */
		pttr_out[2]=tot-ptt;/*PTTR */
		if(tot+1e-18<pis||tot+1e-18<ptt){//sanity check. allow round off error
			warning("tot=%g, pis=%g, ptt=%g\n", tot, pis, ptt);
			warning("coeff=%g,%g,%g,%g,%g,%g\n",
			coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]);
			ans=1;
		}
	}
	/*don't use +=. need locking */
	ngsmod_out[0]=coeff[1];
	ngsmod_out[1]=coeff[2];

	if(ngsmod->indps){
		if(ngsmod->ahstfocus){
			ngsmod_out[ngsmod->indps]=(-2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
		} else{
			ngsmod_out[ngsmod->indps]=(scale1*(coeff[3]+coeff[4]-coeff[0]*MCC_fcp)
				-2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
		}
		ngsmod_out[ngsmod->indps+1]=(scale1*(coeff[3]-coeff[4])
			-2*scale*ht*(thetax*coeff[1]-thetay*coeff[2]));
		ngsmod_out[ngsmod->indps+2]=(scale1*(coeff[5])
			-scale*ht*(thetay*coeff[1]+thetax*coeff[2]));
	}
	if(ngsmod->indastig){
		ngsmod_out[ngsmod->indastig]=(coeff[3]-coeff[4]);
		ngsmod_out[ngsmod->indastig+1]=(coeff[5]);
	}
	if(ngsmod->indfocus){
		ngsmod_out[ngsmod->indfocus]=(coeff[3]+coeff[4]-coeff[0]*MCC_fcp);
	}
	return ans;
}
/**
   Convert NGS modes to DM actuator commands using analytical expression. For >2
   DMs, we only put NGS modes on ground and top-most DM.
*/
static void ngsmod2dm(dcell** dmc, const recon_t* recon, const dcell* M, real gain){
	if(!M||!P(M,0)) return;
	const ngsmod_t* ngsmod=recon->ngsmod;
	//const int nmod=ngsmod->nmod;
	assert(M->nx==1&&M->ny==1&&P(M,0)->nx==ngsmod->nmod);
	real scale=ngsmod->scale;
	/*The MCC_fcp depends weakly on the aperture sampling. */
	real MCC_fcp=ngsmod->aper_fcp;
	loc_t** aloc=recon->aloc->p;
	/*convert mode vector and add to dm commands */
	const int ndm=recon->aloc->nx;
	if(!*dmc){
		*dmc=dcellnew(ndm, 1);
	}
	for(int idm=0; idm<ndm; idm++){
		if(!P(*dmc,idm)){
			P(*dmc,idm)=dnew(P(recon->aloc,idm)->nloc, 1);
		}
	}

	/*first dm */
	real* pm=P(M,0)->p;

	for(int idm=0; idm<ndm; idm++){
		real* p=P(*dmc,idm)->p;
		long nloc=aloc[idm]->nloc;
		real* xloc=aloc[idm]->locx;
		real* yloc=aloc[idm]->locy;
		if(idm==0){
			real focus=0, astigx=0, astigy=0;
			if(ngsmod->indfocus){//focus is a mode
				focus+=pm[ngsmod->indfocus];
			}
			if(ngsmod->indps){//ps mode
				if(ngsmod->ahstfocus){
					focus+=pm[ngsmod->indps]*scale;//scaled to avoid focus mode in science.
				} else{
					focus+=pm[ngsmod->indps];
				}
				astigx+=pm[ngsmod->indps+1];
				astigy+=pm[ngsmod->indps+2];
			}
			if(ngsmod->indastig){
				astigx+=pm[ngsmod->indastig];
				astigy+=pm[ngsmod->indastig+1];
			}
			for(long iloc=0; iloc<nloc; iloc++){
				real xx=xloc[iloc]*xloc[iloc];
				real xy=xloc[iloc]*yloc[iloc];
				real yy=yloc[iloc]*yloc[iloc];
				p[iloc]+=gain*(pm[0]*xloc[iloc]
					+pm[1]*yloc[iloc]
					+focus*(xx+yy-MCC_fcp)
					+astigx*(xx-yy)
					+astigy*(xy));
			}
		} else if(idm+1==ndm&&ngsmod->indps){
			real scale2=-scale*gain;
			for(long iloc=0; iloc<nloc; iloc++){
				real xx=xloc[iloc]*xloc[iloc];
				real xy=xloc[iloc]*yloc[iloc];
				real yy=yloc[iloc]*yloc[iloc];
				p[iloc]+=scale2*(pm[ngsmod->indps]*(xx+yy-MCC_fcp)
					+pm[ngsmod->indps+1]*(xx-yy)
					+pm[ngsmod->indps+2]*(xy));
			}
		}
	}
}
/**
   Convert NGS mode vector to aperture grid for science directions.

   2017-09-11: Deprecated. This routine does not take into account DM 2 science misregistration.
*/

void ngsmod2science(dmat* iopd, const loc_t* loc, const ngsmod_t* ngsmod,
	real thetax, real thetay,
	const real* mod, real alpha){
	const real* locx=loc->locx;
	const real* locy=loc->locy;
	const int nmod=ngsmod->nmod;
	if(nmod==2){//tip/tilt only
		for(int iloc=0; iloc<loc->nloc; iloc++){
			real tmp=locx[iloc]*mod[0]+locy[iloc]*mod[1];
			P(iopd,iloc)+=tmp*alpha;
		}
	} else{
		const real ht=ngsmod->ht;
		const real scale=ngsmod->scale;
		const real scale1=1.-scale;
		real focus=0, ps1=0, ps2=0, ps3=0, astigx=0, astigy=0;
		if(ngsmod->indfocus){
			focus+=mod[ngsmod->indfocus];
		}
		if(ngsmod->indps){
			if(!ngsmod->ahstfocus){
				focus+=mod[ngsmod->indps]*scale1;
			}
			ps1=mod[ngsmod->indps];
			ps2=mod[ngsmod->indps+1];
			ps3=mod[ngsmod->indps+2];
		}
		if(ngsmod->indastig){
			astigx=mod[ngsmod->indastig];
			astigy=mod[ngsmod->indastig+1];
		}
		for(int iloc=0; iloc<loc->nloc; iloc++){
			real x=locx[iloc];
			real y=locy[iloc];
			real xy=x*y;
			real x2=x*x;
			real y2=y*y;
			real tmp=locx[iloc]*mod[0]
				+locy[iloc]*mod[1]
				+focus*(x2+y2)
				+ps1*(-2*scale*ht*(thetax*x+thetay*y))
				+ps2*((x2-y2)*scale1-2*scale*ht*(thetax*x-thetay*y))
				+ps3*(xy*scale1-scale*ht*(thetay*x+thetax*y))
				+astigx*(x2-y2)
				+astigy*(xy);
			P(iopd,iloc)+=tmp*alpha;
		}
	}
}
/**
   frees ngsmod_t
*/
void ngsmod_free(ngsmod_t* ngsmod){
	if(!ngsmod) return;
	cellfree(ngsmod->GM);
	cellfree(ngsmod->Rngs);
	cellfree(ngsmod->Pngs);
	cellfree(ngsmod->Modes);
	dfree(ngsmod->MCC);
	dfree(ngsmod->MCCu);
	cellfree(ngsmod->MCCP);
	dfree(ngsmod->IMCC);
	dfree(ngsmod->IMCC_TT);
	dfree(ngsmod->IMCC_F);
	lfree(ngsmod->modvalid);
	free(ngsmod);
}

/**
   remove NGS modes from LGS DM commands
   if nmod==6: make sure the global focus mode is not removed from LGS result.
*/
void remove_dm_ngsmod(sim_t* simu, dcell* dmerr){
	if(!dmerr) return;
	const recon_t* recon=simu->recon;
	const parms_t* parms=simu->parms;
	const ngsmod_t* ngsmod=recon->ngsmod;
	dcellzero(simu->Mngs);
	dcellmm(&simu->Mngs, ngsmod->Pngs, dmerr, "nn", 1);
	real* mngs=P(simu->Mngs,0)->p;
	if(ngsmod->indastig){//LTAO
	//LTAO is unable to tell where focus/astigmatism occures. NGS WFS needs to control this.
	//Testing: remove LPF'ed focus/astigmatism from LGS DM command. 
		if(!simu->ngsmodlpf){
			simu->ngsmodlpf=dnew(3, 1);
		}
		const real lpfocus=parms->sim.lpfocushi;
		info_once("HPF focus/astig from DM error signal. lpfocus=%g\n", lpfocus);
		P(simu->ngsmodlpf,0)=P(simu->ngsmodlpf,0)*(1-lpfocus)+mngs[ngsmod->indfocus]*lpfocus;
		P(simu->ngsmodlpf,1)=P(simu->ngsmodlpf,1)*(1-lpfocus)+mngs[ngsmod->indastig]*lpfocus;
		P(simu->ngsmodlpf,2)=P(simu->ngsmodlpf,2)*(1-lpfocus)+mngs[ngsmod->indastig+1]*lpfocus;
		mngs[ngsmod->indfocus]=P(simu->ngsmodlpf,0);
		mngs[ngsmod->indastig]=P(simu->ngsmodlpf,1);
		mngs[ngsmod->indastig+1]=P(simu->ngsmodlpf,2);
	}
	if(ngsmod->indfocus&&parms->sim.mffocus&&parms->dbg.ahst_keepfocus){
	/*preserve LGS focus.
	  2019-01-15: Let all focus mode be removed
	  * actually improves performance.*/
		const real scale=ngsmod->scale;
		if(ngsmod->indps&&ngsmod->ahstfocus){
			/* When ahstfocus is true, the first PS mode contains focus mode in
			 * LGS WFS. Relocate this focus mode to the global focus mode to
			 * preserve LGS measurements.*/
			mngs[ngsmod->indfocus]=-mngs[ngsmod->indps]*(scale-1);
		} else{
			/* When ahstfocus is false, the first PS mode does not create focus
			 * mode in LGS WFS.*/
			mngs[ngsmod->indfocus]=0;
		}
	}
	dcellmm(&dmerr, ngsmod->Modes, simu->Mngs, "nn", -1);

}
