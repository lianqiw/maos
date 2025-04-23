/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/*
   2017-09-11: When there is misregistration/distortion between the DM and
   science pupil, the assumed NGS mode on DM remain intact, but the influence on
   Science OPD needs to use ray tracing.
*/
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
	real* amp=P(aper->amp);
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
		const real ht=ngsmod->hdm;
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
   Convert NGS modes to DM actuator commands using analytical expression. For >2
   DMs, we only put NGS modes on ground and top-most DM. 

   Used to create ngsmod->Modes.
*/

static void ngsmod_dm(dcell** dmc, const recon_t* recon, const dcell* M, real gain){
	if(!M||!P(M,0)) return;
	const ngsmod_t* ngsmod=recon->ngsmod;
	//const int nmod=ngsmod->nmod;
	assert(NX(M)==1&&NY(M)==1&&P(M,0)->nx==ngsmod->nmod);
	real scale=ngsmod->scale;
	//The MCC_fcp depends weakly on the aperture sampling. 
	real MCC_fcp=ngsmod->aper_fcp;
	loc_t** aloc=P(recon->aloc);
	//convert mode vector and add to dm commands 
	const int ndm=NX(recon->aloc);
	if(!*dmc){
		*dmc=dcellnew(ndm, 1);
	}
	for(int idm=0; idm<ndm; idm++){
		if(!P(*dmc,idm)){
			P(*dmc,idm)=dnew(P(recon->aloc,idm)->nloc, 1);
		}
	}

	//first dm 
	real* pm=P(P(M,0));

	for(int idm=0; idm<ndm; idm++){
		real* p=P(P(*dmc,idm));
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
   DM modes for all the low order modes, defined on DM grid. It uses ngsmod_dm
   to define the modes*/

static dcell* ngsmod_mode(const parms_t* parms, recon_t* recon){
	ngsmod_t* ngsmod=recon->ngsmod;
	int ndm=parms->ndm;
	int nmod=ngsmod->nmod;
	dcell* M=dcellnew(1, 1);
	P(M,0)=dnew(nmod, 1);
	dcell* mod=dcellnew(ndm, 1);
	dcell* dmt=dcellnew(ndm, 1);
	loc_t** aloc=P(recon->aloc);
	for(int idm=0; idm<ndm; idm++){
		P(dmt,idm)=dnew(aloc[idm]->nloc, 1);
		P(mod,idm)=dnew(aloc[idm]->nloc, nmod);
	}
	for(int imod=0; imod<nmod; imod++){
		dcellzero(dmt);
		P(P(M,0),imod)=1;
		ngsmod_dm(&dmt, recon, M, 1.);
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
void ngsmod_prep(const parms_t* parms, recon_t* recon, const aper_t* aper){
	if(recon->ngsmod){
		warning("Should only be called once\n");
		return;
	}
	ngsmod_t* ngsmod=recon->ngsmod=mycalloc(1, ngsmod_t);
	ngsmod->ahstfocus=parms->tomo.ahst_focus;
	const int ndm=parms->ndm;
	ngsmod->aper_fcp=aper->fcp;
	if(ndm>1&&fabs(parms->dm[0].ht)>1000){
		error("Unsupported configuration. First DM is not on ground\n");
	}
	real hs=NAN;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].lo||parms->powfs[ipowfs].skip) continue;
		if(isnan(hs)){
			hs=parms->powfs[ipowfs].hs;
		} else{
			if(!isinf(hs)||!isinf(parms->powfs[ipowfs].hs)){
				if(fabs(hs-parms->powfs[ipowfs].hs)>1000){
					//High order GS at different altitude.
					hs=INFINITY; break;
				}
			}
		}
	}
	if(isnan(hs)){
		hs=90e3;
	}
	ngsmod->nmod=2; //Basic tip/tilt mode.
	if(!isinf(hs)){	//LGS WFS.
		if(ndm>1&&parms->evl.nevl>1){//Plate scale mode for multi-dm with fov
			ngsmod->indps=ngsmod->nmod;
			ngsmod->nmod+=3;
		} 
		//Always enable focus when LGS WFS is present
		ngsmod->indfocus=ngsmod->nmod;
		ngsmod->nmod+=1;
		if(ngsmod->nmod==3 && parms->nhiwfs>1){//Astigmatism for LTAO
			ngsmod->indastig=ngsmod->nmod;
			ngsmod->nmod+=2;
		}
	}
	info("Low order modes: nmod=%d, mffocus=%d, ahst_focus=%d, indps=%d, indfocus=%d, indastig=%d\n", 
		ngsmod->nmod, parms->sim.mffocus, parms->tomo.ahst_focus, ngsmod->indps, ngsmod->indfocus, ngsmod->indastig);
	ngsmod->hs=hs;
	if(ndm>1){
		ngsmod->hdm=parms->dm[ndm-1].ht;//last DM.
	} else{
		ngsmod->hdm=0;
	}
	ngsmod->scale=pow(1.-ngsmod->hdm/ngsmod->hs, -2);
	/*modal cross coupling matrix along science for performance evaluation. */
	ngsmod->MCCP=ngsmod_mcc(parms, recon, aper, P(parms->evl.wt));
	if(NX(ngsmod->MCCP)==1){
		ngsmod->MCC=dref(P(ngsmod->MCCP,0));
	} else{
		ngsmod->MCC=NULL;
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			dadd(&ngsmod->MCC, 1, P(ngsmod->MCCP,ievl), P(parms->evl.wt,ievl));
		}
	}
	dmat* MCCl=dchol(ngsmod->MCC);
	ngsmod->MCCu=dtrans(MCCl); dfree(MCCl);//convert mode to meter
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
	ngsmod->Modes=ngsmod_mode(parms, recon);
	if(parms->recon.modal){
		dcell *tmp=ngsmod->Modes; ngsmod->Modes=NULL;
		dcellmm(&ngsmod->Modes, recon->amodpinv, tmp, "nn", 1);
		dcellfree(tmp);
	}
	if(recon->actstuck&&!parms->recon.modal&&parms->dbg.recon_stuck){
		warning("Apply stuck actuators to ngs modes\n");
		act_zero(recon->aloc, recon->ngsmod->Modes, recon->actstuck);
	}
	/*
	  ngsmod to NGS gradient interaction matrix. In modal control, the Modes are
	  defined in actuator space, but GA is in modal space for skip==0 WFS.
	*/
	if(parms->recon.split==1&&!parms->tomo.ahst_idealngs&&parms->ntipowfs){
		ngsmod->GM=dcellnew(parms->nwfsr, 1);
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			int ipowfs=parms->wfsr[iwfs].powfs;
			if(parms->powfs[ipowfs].skip!=3 
				&&(parms->powfs[ipowfs].lo||(parms->recon.split && !parms->powfs[ipowfs].trs))){
				for(int idm=0; idm<parms->ndm; idm++){
					dcellmm(&P(ngsmod->GM, iwfs), P(recon->GAlo, iwfs, idm), P(ngsmod->Modes, idm), "nn", 1);
				}
			}
		}
		info2("Low order control includes WFS");
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			if(P(ngsmod->GM, iwfs)) info2(" %d", iwfs);
		}
		info2("\n");
	}
	/**
	   To Decouple LGS reconstruction from ngsmodes.

	   parms->tomo.ahst_wt control NGS modes removal from LGS DM commands
	   if ahst_wt==1: Rngs*GA*dmerr is zero
	   if ahst_wt==2: Doesn't perturb NGS modes in science direction. Retired. use ahst_wt==3 instead.
	   if ahst_wt==3: Identity weighting.
	   if ahst_wt==4: if ndm=2. Move all 2nd order modes from upper DM to ground DM. if ndm==1, remove 2nd order modes.

	   ahst_wt==1 is different from others as it considers aliasing errors of the NGS WFS.
	   ahst_wt==2 or 3 is virtually identical. ahst_wt=3 is retired as of 12/11/2024 since it is similar to 3 and never used.
	*/
	dcell *wts=NULL;
	if(!parms->recon.modal && parms->tomo.ahst_wt>1){
		wts=dcellnew(parms->ndm, 1);
		for(int idm=0; idm<parms->ndm; idm++){
			loc_t *aloc=P(recon->aloc, idm);
			long nloc=aloc->nloc;
			dmat *wt=P(wts,idm)=dnew(nloc, 1);
			real r2max=pow(parms->aper.d*0.5, 2);
			for(long i=0; i<PN(wt); i++){
				real r2=aloc->locx[i]*aloc->locx[i]+aloc->locy[i]*aloc->locy[i];
				P(wt, i)=r2<r2max?1:0;
			}
		}
	}
	if(parms->tomo.ahst_wt>1){/*Identity weighting. */
		ngsmod->Pngs=dcellpinv(ngsmod->Modes, wts);
	}
	if(parms->tomo.ahst_wt==4){
		/*Remove 2nd order modes from the upper DM and place on ground DM with a scaling factor
			which is handled by scaling the aperture D. This method works for more than 2 DMs.*/
		ngsmod->Modes2=dcellnew(parms->ndm, parms->ndm);
		for(int idm=0; idm<parms->ndm; idm++){
			real ht=parms->dm[idm].ht;
			real sc=1.-ht/ngsmod->hs;
			P(ngsmod->Modes2, idm, idm)=zernike(P(recon->aloc, idm), -parms->aper.d*sc, 0, ngsmod->nmod==2?1:2, 0);
			if(NY(P(ngsmod->Modes2, idm, idm))>ngsmod->nmod+1){//Modes2 include piston as an additional term
				dresize(P(ngsmod->Modes2, idm, idm), NX(P(ngsmod->Modes2, idm, idm)), ngsmod->nmod+1);
			}
		}
		if(parms->recon.modal){
			dcell *tmp=ngsmod->Modes2; ngsmod->Modes2=NULL;
			dcellmm(&ngsmod->Modes2, recon->amodpinv, tmp, "nn", 1);
			dcellfree(tmp);
		}
		ngsmod->Pngs2=dcellpinv(ngsmod->Modes2, wts);
	}
	dcellfree(wts);
	if(parms->save.setup){
		writebin(ngsmod->Modes, "ahst_Modes");
		writebin(ngsmod->GM, "ahst_GM");
		if(ngsmod->Pngs) writebin(ngsmod->Pngs, "ahst_Pngs_wt%d", parms->tomo.ahst_wt);
		if(ngsmod->Pngs2) writebin(ngsmod->Pngs2, "ahst_Pngs2_wt%d", parms->tomo.ahst_wt);
		if(ngsmod->Modes2) writebin(ngsmod->Modes2, "ahst_Modes2");
	}
}
/**
   Invert GM while handling rank deficiency.

   It ignores those WFS whoes mask is 0, and assumes the modes are ordered as
   below: T/T x, T/T y, PS1, PS2, PS3, Focus.
*/
static dcell* inv_gm(const dcell* GM, const dspcell* saneai, const lmat* mask, lmat** pmodvalid){
	if(NY(GM)!=1){
		error("To be implemented\n");
	}
	info("Low order reconstruction uses wfs ");
	dcell* GM2=dcellnew(NX(GM), NY(GM));
	int nmod=0, ntt=0, nttf=0;
	for(int iwfs=0; iwfs<NX(GM); iwfs++){
		if((!mask||P(mask, iwfs))&&P(GM, iwfs)&&P(saneai, iwfs, iwfs)->px[0]>0){
			dbg(" %d", iwfs);
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

	for(int iwfs=0; iwfs<NX(GM); iwfs++){
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
   
    When tomo.ahst_wt=1:

	Rngs=GM^\dagger=(GM'*W_N*GM)^{-1}*GM'*W_N  #W_N is recon->saneai
	M^\dagger=Rngs*GA #reconstruct mode from NGS measurements equals to original mode
	Pngs=M*M^\dagger=M*Rngs*GA satisfies Rngs*GA*Pngs=Rngs*GA
	In other words, M^\dagger=(M^T W_G M)^{-1}*M^T W_G with W_G=(GA^T*W_N*GA)

	Projector has the properties: M*Pngs*(1-M*Pngs)=0 
 */
void ngsmod_setup(const parms_t* parms, recon_t* recon){
	ngsmod_t* ngsmod=recon->ngsmod;
	if(parms->recon.split==1&&!parms->tomo.ahst_idealngs&&parms->ntipowfs){
		cellfree(ngsmod->Rngs);
		ngsmod->Rngs=dccellnew(2, 1);
		P(ngsmod->Rngs,0)=inv_gm(ngsmod->GM, recon->saneai, 0, 0);
		
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
		dcellzero(ngsmod->Pngs);
		dcellmm((cell **)&ngsmod->Pngs, P(ngsmod->Rngs, 0), recon->GAlo, "nn", 1); 
	}
	if(parms->save.setup){
		writebin(ngsmod->Rngs, "ahst_Rngs");
		if(ngsmod->Pngs) writebin(ngsmod->Pngs, "ahst_Pngs_wt%d", parms->tomo.ahst_wt);
	}
}
/**
   used in performance evaluation on science opds. accumulate to out*/
void ngsmod_dot(real* pttr_out, real* pttrcoeff_out,
	real* ngsmod_out,
	const parms_t* parms, const ngsmod_t* ngsmod,
	const aper_t* aper, const real* opd, int ievl){
	const real* restrict amp=P(aper->amp);
	const real* restrict locx=aper->locs->locx;
	const real* restrict locy=aper->locs->locy;
	double coeff[6]={0,0,0,0,0,0};//always use double to enhance precision of accumulation
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
#if CPU_SINGLE	
	real coeff2[6]={coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]};
#else
	#define coeff2 coeff
#endif	
	ngsmod_dot_post(pttr_out, pttrcoeff_out, ngsmod_out, tot, coeff2, ngsmod, aper, thetax, thetay);
}
/**
   Separate post processing part so that GPU code can call it. Return non zero if error happens.
*/
int ngsmod_dot_post(real* pttr_out, real* pttrcoeff_out, real* ngsmod_out,
	real tot, const real* coeff, const ngsmod_t* ngsmod,
	const aper_t* aper, real thetax, real thetay){
	const real MCC_fcp=ngsmod->aper_fcp;
	const real hdm=ngsmod->hdm;
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
		if(pttr_out[0]<0 || pttr_out[1]<0 || pttr_out[2]<0){
			if(pttr_out[0]<-1e-6 || pttr_out[1]<-1e-6 || pttr_out[2]<-1e-6){
				//allow round off error
				warning("tot=%g, pis=%g, ptt=%g\n", tot, pis, ptt);
				warning("coeff=%g,%g,%g,%g,%g,%g\n",
				coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]);
				ans=1;
			}
			if(pttr_out[0]<0) pttr_out[0]=0;
			if(pttr_out[1]<0) pttr_out[1]=0;
			if(pttr_out[2]<0) pttr_out[2]=0;
		}
	}
	/*don't use +=. need locking */
	ngsmod_out[0]=coeff[1];
	ngsmod_out[1]=coeff[2];

	if(ngsmod->indps){
		if(ngsmod->ahstfocus){
			ngsmod_out[ngsmod->indps]=(-2*scale*hdm*(thetax*coeff[1]+thetay*coeff[2]));
		} else{
			ngsmod_out[ngsmod->indps]=(scale1*(coeff[3]+coeff[4]-coeff[0]*MCC_fcp)
				-2*scale*hdm*(thetax*coeff[1]+thetay*coeff[2]));
		}
		ngsmod_out[ngsmod->indps+1]=(scale1*(coeff[3]-coeff[4])
			-2*scale*hdm*(thetax*coeff[1]-thetay*coeff[2]));
		ngsmod_out[ngsmod->indps+2]=(scale1*(coeff[5])
			-scale*hdm*(thetay*coeff[1]+thetax*coeff[2]));
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
   Obtain OPD on loc for NGS mode vector.
*/

void ngsmod_opd(dmat* iopd, const loc_t* loc, const ngsmod_t* ngsmod,
	real thetax, real thetay, const real* mode, real alpha){
	const real* locx=loc->locx;
	const real* locy=loc->locy;
	const int nmod=ngsmod->nmod;
	if(nmod==2){//tip/tilt only
		for(int iloc=0; iloc<loc->nloc; iloc++){
			real tmp=locx[iloc]*mode[0]+locy[iloc]*mode[1];
			P(iopd,iloc)+=tmp*alpha;
		}
	} else{
		const real hdm=ngsmod->hdm;
		const real scale=ngsmod->scale;
		const real scale1=1.-scale;
		real focus=0, ps1=0, ps2=0, ps3=0, astigx=0, astigy=0;
		if(ngsmod->indfocus){
			focus+=mode[ngsmod->indfocus];
		}
		if(ngsmod->indps){
			if(!ngsmod->ahstfocus){
				focus+=mode[ngsmod->indps]*scale1;
			}
			ps1=mode[ngsmod->indps];
			ps2=mode[ngsmod->indps+1];
			ps3=mode[ngsmod->indps+2];
		}
		if(ngsmod->indastig){
			astigx=mode[ngsmod->indastig];
			astigy=mode[ngsmod->indastig+1];
		}
		for(int iloc=0; iloc<loc->nloc; iloc++){
			real x=locx[iloc];
			real y=locy[iloc];
			real xy=x*y;
			real x2=x*x;
			real y2=y*y;
			real tmp=locx[iloc]*mode[0]
				+locy[iloc]*mode[1]
				+focus*(x2+y2)
				+ps1*(-2*scale*hdm*(thetax*x+thetay*y))
				+ps2*((x2-y2)*scale1-2*scale*hdm*(thetax*x-thetay*y))
				+ps3*(xy*scale1-scale*hdm*(thetay*x+thetax*y))
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
	cellfree(ngsmod->Pngs2);
	cellfree(ngsmod->Modes);
	cellfree(ngsmod->Modes2);
	dfree(ngsmod->MCC);
	dfree(ngsmod->MCCu);
	cellfree(ngsmod->MCCP);
	dfree(ngsmod->IMCC);
	dfree(ngsmod->IMCC_TT);
	dfree(ngsmod->IMCC_F);
	lfree(ngsmod->modvalid);
	free(ngsmod);
}
/*
	Split NGS mode error from dmerr 
*/
void ngsmod_split(dcell **Merr, sim_t *simu, dcell *dmerr){
	const ngsmod_t *ngsmod=simu->recon->ngsmod;
	dcellzero(*Merr);
	dcellmm(Merr, simu->recon->ngsmod->Pngs, dmerr, "nn", 1);
	dcellmm(&dmerr, ngsmod->Modes, *Merr, "nn", -1);
}
/**
   remove NGS modes from LGS DM commands. \todo handle 3 DM case.
*/
void ngsmod_remove(sim_t* simu, dcell* dmerr){
	const recon_t* recon=simu->recon;
	const parms_t* parms=simu->parms;
	const ngsmod_t* ngsmod=recon->ngsmod;
	if(!dmerr||!ngsmod->Pngs) return;
	dcellzero(simu->Mngs);
	if(parms->tomo.ahst_wt==4){//this is equivalent to ahst_wt=1 for 2 dms.
		dcellmm(&simu->Mngs, ngsmod->Pngs2, dmerr, "nn", 1);
		//if ndm==1, remove all 2nd order modes from ground DM since LGS cannot sense them properly
		//if ndm>1 and ahst_keepfocus=0, remove focus from all DMs. LGS focus is not reliable
		//merge 2nd order modes from upper DMs to ground DM.
		const int im0=parms->ndm==1?ngsmod->nmod:(parms->dbg.ahst_keepfocus?3:4);
		for(int im=im0; im<ngsmod->nmod; im++){
			P(P(simu->Mngs, 0), im)=0;
			for(int idm=1; idm<parms->ndm; idm++){
				P(P(simu->Mngs, 0), im)-=P(P(simu->Mngs, idm), im);
			}
		}
		dcellmm(&dmerr, ngsmod->Modes2, simu->Mngs, "nn", -1);//t/t are removed from all dms.
	}else{
		dcellmm(&simu->Mngs, ngsmod->Pngs, dmerr, "nn", 1);
		real *mngs=P(P(simu->Mngs, 0));
		if(ngsmod->indastig&&parms->sim.mffocus){//LTAO
			//LTAO is unable to tell where focus/astigmatism occures. NGS WFS needs to control this.
			//Solution: remove LPF'ed focus/astigmatism from LGS DM command. 
			const real lpfocus=parms->sim.lpfocushi;
			if(lpfocus<1){
				if(!simu->ngsmodlpf){
					simu->ngsmodlpf=dnew(3, 1);
				}
				info_once("HPF focus/astig from DM error signal. lpfocus=%g\n", lpfocus);
				P(simu->ngsmodlpf,0)=P(simu->ngsmodlpf,0)*(1-lpfocus)+mngs[ngsmod->indfocus]*lpfocus;
				P(simu->ngsmodlpf,1)=P(simu->ngsmodlpf,1)*(1-lpfocus)+mngs[ngsmod->indastig]*lpfocus;
				P(simu->ngsmodlpf,2)=P(simu->ngsmodlpf,2)*(1-lpfocus)+mngs[ngsmod->indastig+1]*lpfocus;
				mngs[ngsmod->indfocus]=P(simu->ngsmodlpf,0);
				mngs[ngsmod->indastig]=P(simu->ngsmodlpf,1);
				mngs[ngsmod->indastig+1]=P(simu->ngsmodlpf,2);
			}else{
				info_once("Remove focus/astig from DM error signal.\n");
			}
		}else if(parms->dbg.ahst_keepfocus&&ngsmod->indfocus&&parms->sim.mffocus){
			/*The LGS focus mode is filtered from gradient input to the reconstructor, so we do not remove it again here*/
			if(ngsmod->indps&&ngsmod->ahstfocus){
				/* When ahstfocus is true, the first PS mode contains focus mode in
				* LGS WFS. Relocate this focus mode to the global focus mode to
				* preserve LGS focus measurements.*/
				mngs[ngsmod->indfocus]=-mngs[ngsmod->indps]*(ngsmod->scale-1);
			} else{
				mngs[ngsmod->indfocus]=0;//preserve focus
			}
		}
		dcellmm(&dmerr, ngsmod->Modes, simu->Mngs, "nn", -1);
	}
	//dshow(P(simu->Mngs, 0), "Mngs");
}
