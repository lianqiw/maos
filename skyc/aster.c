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
#include "skyc.h"
#include "aster.h"
#include "photon.h"
#include "physim.h"
#include "mtch.h"
#include "utils.h"
#include "star.h"
/**
   \file skyc/setup_aster.c
   Routines to handle asterisms.
 */
/**
   Create combination of stars to form asterism. It has the option to put TTF
   always on the brightest for testing purpose.  */
aster_s* setup_aster_comb(int* naster, const star_s* star, int nstar, const parms_s* parms){
	if(!star || !nstar){
		*naster=0;
		return NULL;
	}
	int flag=parms->skyc.keeporder?0:(parms->skyc.ttfbrightest?2:1);
	int npowfs=parms->maos.npowfs;
	lmat *starvalid=lnew(nstar, npowfs);
	for(int istar=0; istar<nstar; istar++){
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			P(starvalid, istar, ipowfs)=star[istar].use[ipowfs];
		}
	}
	lmat *comb=comb_stars(parms->skyc.nwfsmax, starvalid, flag);
	*naster=NY(comb);
	int nwfstot=NX(comb);
	aster_s *aster=mycalloc(*naster, aster_s);
	for(int iaster=0; iaster<*naster; iaster++){
		aster[iaster].nwfs=nwfstot;
		aster[iaster].wfs=mycalloc(nwfstot, wfs_s);
		aster[iaster].iaster=iaster;
		int iwfs=0;
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			for(int jwfs=0; jwfs<P(parms->skyc.nwfsmax, ipowfs); jwfs++){
				if(iwfs<nwfstot){
					int istar=P(comb, iwfs, iaster);
					aster[iaster].wfs[iwfs].ipowfs=ipowfs;
					aster[iaster].wfs[iwfs].istar=istar;
					aster[iaster].wfs[iwfs].thetax=star[istar].thetax;
					aster[iaster].wfs[iwfs].thetay=star[istar].thetay;
					/*Magnitude */
					aster[iaster].wfs[iwfs].mags=star[istar].mags;//do not free
					/*Signal Level */
					aster[iaster].wfs[iwfs].siglev=P(star[istar].siglev,ipowfs);//do not free
					aster[iaster].wfs[iwfs].siglevtot=P(star[istar].siglevtot,ipowfs);
					aster[iaster].wfs[iwfs].bkgrnd=P(star[istar].bkgrnd,ipowfs);

					/*Pixel intensity statistics. */
					aster[iaster].wfs[iwfs].pistat=&star[istar].pistat[ipowfs];
				}
				iwfs++;
			}
		}
	}
	lfree(comb);
	lfree(starvalid);
	if(parms->skyc.verbose){
		info("Number of stars: %d, number of asterisms: %d\n", nstar, *naster);
	}
	return aster;
}
/**
   Compute Modal to gradient operator by copying from the stars. Using average
gradients. Similar to Z tilt since the mode is low order */
void setup_aster_gm(aster_s* aster, star_s* star, const parms_s* parms){
	aster->g=dcellnew(aster->nwfs, 1);
	aster->ngs=mycalloc(aster->nwfs, long);
	aster->tsa=0;
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		const int istar=aster->wfs[iwfs].istar;
		const int ipowfs=aster->wfs[iwfs].ipowfs;
		const long nsa=parms->maos.nsa[ipowfs];
		P(aster->g,iwfs)=ddup(P(star[istar].g,ipowfs));
		aster->tsa+=nsa;
		aster->ngs[iwfs]=nsa*2;
	}
	/*
	  aster->g is also used for simulation. Do not zero columns here.
	*/
	aster->gm=dcell2m(aster->g);
}

/**
   Copy time history of complex pupil function from star_s to aster_s.
 */
void setup_aster_wvf(aster_s* aster, star_s* star){
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		const int ipowfs=aster->wfs[iwfs].ipowfs;
		const int istar=aster->wfs[iwfs].istar;
		aster->wfs[iwfs].wvfout=star[istar].wvfout[ipowfs];
	}
}
/**
   Copy time history of complex pupil function from star_s to aster_s.
 */
void setup_aster_ztilt(aster_s* aster, star_s* star){
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		const int ipowfs=aster->wfs[iwfs].ipowfs;
		const int istar=aster->wfs[iwfs].istar;
		aster->wfs[iwfs].ztiltout=P(star[istar].ztiltout,ipowfs);
		if(star[istar].goff){
			aster->wfs[iwfs].goff=P(star[istar].goff,ipowfs);
		}
	}
}
/**
  Estimate wavefront error propagated from measurement error. pgm is the reconstructor. ineam is the
  error inverse.
  trace(Mcc*(pgm*neam*pgm'))
*/
static dmat* calc_recon_error(const dmat* pgm,   /**<[in] the reconstructor*/
	const dmat* neam,/**<[in] the inverse of error covariance matrix*/
	const dmat* mcc   /**<[in] NGS mode covariance matrix.*/
){
	dmat* psp=NULL;
	dmat* tmp=NULL;
	dmat* var=NULL;
	dcp(&tmp, pgm);
	dmuldiag(tmp, neam);
	dmm(&psp, 0, tmp, pgm, "nt", 1);
	dfree(tmp);
	dmm(&var, 0, mcc, psp, "nn", 1);
	dfree(psp);
	dmat* res=dnew(mcc->nx+1, 1);
	/*It is right for both ix, iy to stop at ib.*/
	for(int ib=0; ib<mcc->ny; ib++){
		P(res,ib)=P(var, ib, ib);
		if(P(res,ib)<0){
			if(P(res,ib)<-1e-30){
				dbg_once("Negative noise (%g) is set to 0.\n", P(res, ib));
			}
			P(res,ib)=0;
		}
		P(res, mcc->ny)+=P(res, ib);//trace
	}
	dfree(var);
	return res;
}

/**
   Interpolate simu->gain based on noise.
*/
static void interp_gain(real* out, const dcell* gain, const dmat* gainx,
	real sigma2){
	const long nx=gainx->nx;
	const real xsep=(log(P(gainx,nx-1))-log(P(gainx,0)))/(nx-1);
	const real xx=(log(sigma2)-log(P(gainx,0)))/xsep;
	int ig;
	if(xx<0){/*too small. */
		ig=0;
	} else if(xx>=nx-1){/*too big */
		ig=nx-1;
	} else{/*within the range */
		ig=ifloor(xx);
		/*2013-12-06: use one of the set, not interpolate*/
	}
	memcpy(out, P(P(gain,ig)), sizeof(real)*P(gain,ig)->nx);
}
/**
 * Select GM using WFS mask. Also compute direct mode mask (mdirect) if pmdirect is set.
 * direct mode are modes that are only sensed by the slower loop and therefore directly fed into the integrator.
 * */
static dmat* setup_aster_mask_gm(const dcell* gm_in, const lmat* mask, lmat** pmdirect){
	if(!gm_in||gm_in->nx==0) return NULL;
	dcell* gm=dcelldup(gm_in);
	int nwfs=gm->nx; assert(gm->ny==1);
	int nmod=P(gm,0)->ny;
	int ntt=0, nttf=0;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		if(mask&&!P(mask,iwfs)){
			dzero(P(gm, iwfs));
		} else{
			int ng=P(gm, iwfs)->nx;
			if(ng>=8){
				nttf++;
			} else if(ng==2){
				ntt++;
			} else{
				error("Unknown WFS type: ng=%d\n", ng);
			}
		}
	}
	dmat* ggm=dcell2m(gm); dcellfree(gm);
	lmat *mdirect=(pmdirect && *pmdirect)?*pmdirect:lnew(nmod,1);
	
	if(nttf>0){//there is TTF
		if(nttf+ntt==1&&nmod==6){//only 1 TTF
			//1 ttf cannot control both focus and magnification. 
			P(mdirect, 2)=1;
		}
		//P(mdirect, 3)=1;
		//P(mdirect, 4)=1;
	} else{//Only TT WFS
		if(nmod>5){//no focus control.
			P(mdirect, 5)=1;
		}
		if(ntt<3&&nmod>=5){//1 or 2 TT OIWFS can not control all 3 PS modes
			P(mdirect, 2)=1;
		}
		if(ntt<2&&nmod>=5){//1 TT OIWFS can not control any PS mode.
			P(mdirect, 3)=1;
			P(mdirect, 4)=1;
		}
		if(ntt<1){
			lzero(mdirect);
			warning("there is no WFS\n");
		}
	}
	for(int imod=0; imod<nmod; imod++){
		if(P(mdirect, imod)){
			dzerocol(ggm, imod);
		}
	}
	if(!pmdirect){
		lfree(mdirect);
	}else if(!*pmdirect){
		*pmdirect=mdirect;
	}
	return ggm;
}
/**
  	For the multirate case, setup the dtrat of each WFS.
	- First, the fastest rate of each wfs is computed based on snr threshold.
	- If skyc.ttffastest is set, the fast rate is forced to be the TTF OIWFS rate. Otherwise, it is set to the fastest OIWFS.
*/
static void
setup_aster_dtrat_multirate(aster_s* aster, const parms_s* parms){
	const int ndtrat=parms->skyc.ndtrat;
	lfree(aster->idtrats);
	lfree(aster->dtrats);
	aster->idtrats=lnew(aster->nwfs, 1);
	aster->dtrats=lnew(aster->nwfs, 1);
	int idtrat_fast=0;
	//real snr_fast=0;
	int idtrat_slow=parms->skyc.ndtrat;
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		int idtrat;
		for(idtrat=ndtrat-1; idtrat>=0; idtrat--){
			//select highest idtrat that meets snrmin for each WFS.
			if(P(aster->wfs[iwfs].pistat->snr,idtrat)>=P(parms->skyc.snrmin_mr,idtrat)){
				P(aster->idtrats,iwfs)=idtrat;
				break;
			}
		}
		if(idtrat==-1){
			warning("Invalid asterism: wfs%d idtrat=%d\n", iwfs, idtrat);
			idtrat=0;
		}
		int isttf=parms->maos.nsa[aster->wfs[iwfs].ipowfs]>1;
		if((parms->skyc.ttffastest&&isttf)||(!parms->skyc.ttffastest &&idtrat>idtrat_fast)){
			idtrat_fast=idtrat;
			//snr_fast=P(aster->wfs[iwfs].pistat->snr, idtrat);
			if(idtrat_slow==parms->skyc.ndtrat) idtrat_slow=idtrat;
		}
		if(idtrat<idtrat_slow){
			//real snr2=P(aster->wfs[iwfs].pistat->snr, idtrat_fast);
			//if(!isttf && idtrat+2 >= idtrat_fast && snr2>3){
				//We speed up TT wfs that is slower than fast WFS to the fast rate under the following conditions:
				//The difference in frequency is no more than 4 times.
				//snr is better than 3
				//snr is better than snr_fast-5
				//We don't speed up TTF wfs as it will degrade focus correction too much.
			//	idtrat=P(aster->idtrats, iwfs)=idtrat_fast;
			//}
			//info("aster %d snr=%g %g\n", aster->iaster, snr_fast, snr);
			idtrat_slow=idtrat;
		}
	}
	int idtrat_fast2=MAX(3, idtrat_fast);//force fast loop to be faster than dtrat=8
	//We use only two rates: a faster rate for tip/tilt(/focus) control and a slower rate for focus/ps control.
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		//force fast loop to be faster than dtrat=16
		if(P(aster->idtrats, iwfs)>=idtrat_fast){
			P(aster->idtrats, iwfs)=idtrat_fast2;
		}
		if(P(aster->idtrats, iwfs)<idtrat_fast){
			P(aster->idtrats, iwfs)=idtrat_slow;
		}
		P(aster->dtrats, iwfs)=P(parms->skyc.dtrats, P(aster->idtrats, iwfs));
	}
	aster->idtratest=idtrat_fast2;
	aster->idtratmin=idtrat_fast2;
	aster->idtratmax=idtrat_fast2+1;
}

/**

   Setup the least squares reconstructor and controller.

   It first computers the reconstructor, and noise propagation, and then
   optimize servo gains to reduce the combined noise propagation and residual
   error.

   We try to minimize

   \f[
   \sigma^2=\int \textrm{PSD}_{ngs,ws}H_{rej}\textrm{d}\nu + \int_0^{\nu_{nyquist}} \textrm{PSF}\textrm{d}\nu
   \f]
*/
static void setup_aster_servo(sim_s* simu, aster_s* aster, const parms_s* parms){
	const int multirate=parms->skyc.multirate;
	int ncase=0;
	if(aster->gain){
		dcellfree(aster->pgm);
		dcellfree(aster->sigman);
		dcellfree(aster->gain);
		dfree(aster->res_ws);
		dfree(aster->res_ngs);
	}
	lmat* case_valid=0;
	int max_idtrat=0;

	if(multirate){
		ncase=(1<<aster->nwfs)-1;
		setup_aster_dtrat_multirate(aster, parms);
		int min_idtrat=parms->skyc.ndtrat;
		//Find the dtrat of the slowest WFS
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
			int idtrat=P(aster->idtrats,iwfs);
			if(min_idtrat>idtrat){
				min_idtrat=idtrat;
			}
			if(max_idtrat<idtrat){
				max_idtrat=idtrat;
			}
		}
		case_valid=lnew(ncase, 1);
		//Loop over to find all possible combinations
		for(int isim=0; isim<P(parms->skyc.dtrats,min_idtrat); isim++){
			int indk=0; int count=0;
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				if((isim%P(aster->dtrats,iwfs))==0){
					indk|=(1<<iwfs);
					count++;
				}
			}
			if(count){
				P(case_valid,indk-1)=indk;
			}
		}
	} else{
		ncase=parms->skyc.ndtrat;
	}
	aster->pgm=dcellnew(ncase, 1);
	aster->sigman=dcellnew(ncase, 1);
	aster->gain=dcellnew(ncase, 1);
	aster->res_ws=dnew(ncase, 1);
	aster->res_ngs=dnew(ncase, 3);
	dmat* pres_ngs=aster->res_ngs;
	dmat* gm=multirate?0:setup_aster_mask_gm(aster->g, NULL, NULL);//no mask
	real resngs_multi=0;//multirate error
	const long nmod=parms->maos.nmod;
	int icase_fast=-1;//fast loop index into case_valid, pgm, gain, etc.
	const int icase_slow=ncase-1;//slow loop is always the last entry.
	for(int icase=0; icase<ncase; icase++){
		//Assemble the measurement error covariance matrix
		if(multirate){
		 	if(!P(case_valid,icase)){
				dfree(P(aster->pgm, icase));
				continue;
			}else{
				if(icase_fast==-1){
					icase_fast=icase;
				}
			}
		}
		dcell* nea=dcellnew3(aster->nwfs, 1, aster->ngs, 0);
		lmat* mask=multirate?lnew(aster->nwfs, 1):0;
		int idtrat=icase;//for non multirate
		if(multirate){//for multirate, idtrat is the slowest of all active WFS
			idtrat=-1;
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				if((icase+1)&(1<<iwfs)){
					if(idtrat==-1||idtrat>P(aster->idtrats, iwfs)){
						idtrat=P(aster->idtrats,iwfs);
					}
					P(mask, iwfs)=1;
				}
			}
		}
		//for multrate, we use sanea at the slower rate for fast WFS due to averaging.
		//info("aster %d: icase=%d, idtrat=%d\n", aster->iaster, icase, idtrat);

		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
			if(!mask || P(mask, iwfs)){
				dcp(&P(nea,iwfs), P(aster->wfs[iwfs].pistat->sanea,idtrat));
				dcwpow(P(nea,iwfs), -2);
			}
		}
		//dshow(nea->m, "nea");
		if(multirate){
			//Reconstructor
			dfree(gm);
			//modes not measured by the faster loop needs to be directly output to integrator by the slower loop
			//todo: implement for more than 2 rates
			gm=setup_aster_mask_gm(aster->g, mask, (icase+1!=ncase)?&aster->mdirect:NULL);
		}
		P(aster->pgm,icase)=dpinv(gm, nea->m);

		//Noise propagation
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
			if(!mask || P(mask,iwfs)){
				dcwpow(P(nea,iwfs), -1);//inverse again
			}
		}
		P(aster->sigman,icase)=calc_recon_error(P(aster->pgm,icase), nea->m, parms->maos.mcc);
		if(parms->skyc.dbg){
			writebin(nea, "%s/aster%d_nea%d", dirsetup, aster->iaster, icase);
		}
		dcellfree(nea);
		lfree(mask);

		/*gsplit:
		  0: All modes use the same gain.
		  1: Different mode use different gain (2017-0-24) was only tt/ps separate.
		  note that P(simu->psds,0) contains windshake PSD.
		*/
		real res_ngs=0;/*residual error due to signal after servo rejection. */
		real res_ngsn=0;/*residual error due to noise. */
		const int servotype=parms->skyc.servo;
		const int ng=parms->skyc.ngain;//number of gain parameters
		dmat* pgain=P(aster->gain,icase)=dnew(ng, nmod);
		for(int ipsd=0; ipsd<simu->psds->nx; ipsd++){
			real sigma=P(P(aster->sigman,icase),parms->skyc.gsplit?ipsd:nmod);
			real pg[ng+2];
			if(parms->skyc.interpg){
				interp_gain(pg, P(P(simu->gain_pre,idtrat),ipsd), simu->gain_x, sigma);
			} else{
				dmat* sigma2=dnew(1, 1);
				P(sigma2,0)=sigma;
				dcell *tmp=servo_optim(parms->maos.dt, P(parms->skyc.dtrats, idtrat), 0, parms->skyc.pmargin, 0, 0, servotype, P(simu->psds, ipsd), sigma2);
				memcpy(pg, P(P(tmp,0)), (ng+2)*sizeof(real));
				dcellfree(tmp);
				dfree(sigma2);
			}
			res_ngs+=pg[ng];
			res_ngsn+=pg[ng+1];
			memcpy(PCOL(pgain, ipsd), pg, sizeof(real)*ng);
			if(multirate){
				if(icase+1!=ncase){//this is fast loop
					if(!P(aster->mdirect, ipsd)){//mode driven by fast loop
						resngs_multi+=pg[ng]+pg[ng+1];
						//info("aster %d mode %d fast\n", aster->iaster, ipsd);
					}
				}else{//slow loop
					if(aster->mdirect&&!P(aster->mdirect, ipsd)){//there is mode driven by fast loop
						resngs_multi+=pg[ng+1];//only count noise propation
						//info("aster %d mode %d slow noisy only\n", aster->iaster, ipsd);
					}else{//driven by sloper loop
						resngs_multi+=pg[ng]+pg[ng+1];
						//info("aster %d mode %d slow\n", aster->iaster, ipsd);
					}
				}
			}
		}
		if(simu->psds->nx < nmod){
			warning_once("Use gain of mode 0 for modes with no PSD\n");
			for(int imod=simu->psds->nx; imod<nmod; imod++){
				memcpy(PCOL(pgain, imod), PCOL(pgain, 0), sizeof(real)*ng);
			}
		}
		P(pres_ngs, icase, 0)=res_ngs+res_ngsn;/*error due to signal and noise */
		P(pres_ngs, icase, 1)=res_ngs;/*error due to signal */
		P(pres_ngs, icase, 2)=res_ngsn;/*error due to noise propagation. */
		
		dmat* g_tt=drefcols(pgain, 0, 1);
		real gain_n;
		if(parms->skyc.psd_ws){
			P(aster->res_ws,icase)=servo_residual(&gain_n, parms->skyc.psd_ws,
				parms->maos.dt, P(parms->skyc.dtrats,idtrat), 0, g_tt, servotype);
		}
		dfree(g_tt);
	}//for dtrat
	if(multirate){//for setup_aster_select() to use.
		P(pres_ngs, 0, 0)=resngs_multi;
		
		aster->minest=P(aster->res_ngs, 0, 0);
		for(int imod=0; imod<nmod; imod++){
			if(aster->mdirect && P(aster->mdirect, imod)){//mode only measured by slow loop. Use the slow loop gain.
				P(P(aster->gain, icase_fast), imod)=P(P(aster->gain, icase_slow), imod);
			}/*else if(P(P(aster->gain, icase_slow), imod)>1.){
				P(P(aster->gain, icase_slow), imod)=1.;//avoid gain >1 to help stability
			}*/
		}
	}
	if(parms->skyc.dbg){
		writebin(gm, "%s/aster%d_gm", dirsetup, aster->iaster);
		writebin(aster->pgm, "%s/aster%d_pgm", dirsetup, aster->iaster);
		writebin(aster->sigman, "%s/aster%d_sigman", dirsetup, aster->iaster);
		writebin(aster->gain, "%s/aster%d_gain", dirsetup, aster->iaster);
		if(parms->skyc.psd_ws){
			writebin(aster->res_ws, "%s/aster%d_res_ws", dirsetup, aster->iaster);
		}
		writebin(aster->res_ngs, "%s/aster%d_res_ngs", dirsetup, aster->iaster);
	}
	dfree(gm);
	lfree(case_valid);
}
/**
   Setup the dtrat of other WFS. Not allowed to be faster than the first wfs (consider revise).
*/
static void setup_aster_kalman_multirate(aster_s* aster, const parms_s* parms, int idtrat_wfs0){
	if(parms->skyc.verbose){
		info("aster %d dtrat_wfs0=%3d, dtrat=", aster->iaster,
			(int)P(parms->skyc.dtrats,idtrat_wfs0));
	}
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		int idtrat=idtrat_wfs0;
		if(iwfs>0){
			/*don't allow snr to fall below 3.*/
			while(idtrat>0&&(P(aster->wfs[iwfs].pistat->snr,idtrat)<3
				||(int)P(parms->skyc.dtrats,idtrat)%(int)P(aster->dtrats,0)!=0)){
				idtrat--;
			}
		}
		P(aster->idtrats,iwfs)=idtrat;
		P(aster->dtrats,iwfs)=P(parms->skyc.dtrats,idtrat);
		int ng=P(aster->g,iwfs)->nx;
		if(idtrat>-1){
			for(int ig=0; ig<ng; ig++){
				P(P(P(aster->neam,0), iwfs, iwfs),ig,ig)=
					pow(P(P(aster->wfs[iwfs].pistat->sanea,idtrat),ig), 2);
			}
		} else{//no star available
			dset(P(P(aster->neam,0), iwfs, iwfs), 0);
		}
		if(parms->skyc.verbose){
			info("%3d ", (int)P(parms->skyc.dtrats,idtrat));
		}
	}//for iwfs
}

static void setup_aster_kalman(sim_s* simu, aster_s* aster, const parms_s* parms){
	int ndtrat=parms->skyc.ndtrat;
	if(parms->skyc.multirate){
		aster->res_ngs=dnew(1, 1);
		//assemble neam
		if(aster->neam) error("neam is already set?\n");
		aster->neam=dccellnew(1, 1);
		P(aster->neam,0)=dcellnew(aster->nwfs, aster->nwfs);
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
			int ng=P(aster->g,iwfs)->nx;
			P(P(aster->neam,0), iwfs, iwfs)=dnew(ng, ng);
		}
		aster->dtrats=lnew(aster->nwfs, 1);
		aster->idtrats=lnew(aster->nwfs, 1);
		int wfs0_min=0, wfs0_max=0;
		pistat_s* pistat0=aster->wfs[0].pistat;//&star[aster->wfs[0].istar].pistat[aster->wfs[0].ipowfs];
		for(int idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
			if(P(pistat0->snr,idtrat)>=P(parms->skyc.snrmin_mr,idtrat)){
				if(wfs0_min==0){
					wfs0_min=idtrat;
				}
				wfs0_max=idtrat;
			}
		}
		aster->kalman=mycalloc(1, kalman_t*);
		real resmin=INFINITY;
		kalman_t* kalman_min=0;
		int idtrat_min=0;
		//Try progressively lower sampling frequencies until performance starts to degrades
		for(int idtrat_limit=wfs0_max; idtrat_limit>wfs0_min; idtrat_limit--){
			setup_aster_kalman_multirate(aster, parms, idtrat_limit);//setup aster->dtrats and aster->neam
			aster->kalman[0]=sde_kalman(simu->sdecoeff, parms->maos.dt, aster->dtrats, aster->g, P(aster->neam,0), 0);
			dmat* rests=0;
#define USE_SIM 0 //ztiltout is not available here.
#if USE_SIM   //more accurate
			dmat* res=physim(parms->skyc.dbg?&rests:0, simu->mideal, simu->mideal_oa, simu->varol,
				aster, 0, parms, -1, 1, -1);
			real res0=res?P(res,0):simu->varol;
			dfree(res);
#else
			rests=kalman_test(aster->kalman[0], simu->mideal);
			real res0=calc_rms(rests, parms->maos.mcc, parms->skyc.evlstart);
#endif
			if(parms->skyc.dbg){
				writebin(rests, "isky%d_iaster%d_dtrat%d_rest", simu->isky, aster->iaster, idtrat_limit);
			}
			if(parms->skyc.verbose) info("res0=%g, resmin=%g\n", sqrt(res0)*1e9, sqrt(resmin)*1e9);
			dfree(rests);
			if(res0<resmin-100e-18){//better by 10 nm
				resmin=res0;
				kalman_free(kalman_min);
				kalman_min=aster->kalman[0];
				aster->kalman[0]=0;
				idtrat_min=idtrat_limit;
			} else{
				kalman_free(aster->kalman[0]);aster->kalman[0]=0;
				if(!isinf(resmin)&&res0>resmin*2){//stop trying.
					break;
				}
			}
		}
		setup_aster_kalman_multirate(aster, parms, idtrat_min);
		if(parms->skyc.verbose) info("selected\n");
		aster->minest=P(aster->res_ngs,0,0)=resmin;
		aster->kalman[0]=kalman_min;
		if(parms->skyc.dbg){
			kalman_write(aster->kalman[0], "%s/aster%d_kalman_%g", dirsetup, aster->iaster, P(parms->skyc.dtrats,idtrat_min));
		}
	} else{
		if(aster->neam) cellfree(aster->neam);
		aster->neam=dccellnew(ndtrat, 1);
		aster->res_ngs=dnew(ndtrat, 3);
		aster->kalman=mycalloc(ndtrat, kalman_t*);
		lmat* dtrats=lnew(aster->nwfs, 1);
		for(int idtrat=0; idtrat<ndtrat; idtrat++){
			//assemble neam
			//TIC;tic;
			P(aster->neam,idtrat)=dcellnew(aster->nwfs, aster->nwfs);
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				dmat* tmp=ddup(P(aster->wfs[iwfs].pistat->sanea,idtrat));/*in rad */
				dcwpow(tmp, 2);
				dsp* tmp2=dspnewdiag(tmp->nx, P(tmp), 1);
				dspfull(&P(P(aster->neam,idtrat), iwfs, iwfs), tmp2, 'n', 1);
				dfree(tmp); dspfree(tmp2);
			}

			int dtrat=P(parms->skyc.dtrats,idtrat);
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				P(dtrats,iwfs)=dtrat;
			}
			aster->kalman[idtrat]=sde_kalman(simu->sdecoeff, parms->maos.dt, dtrats, aster->g, P(aster->neam,idtrat), 0);
			//toc2("kalman");
#if USE_SIM 
			dmat* res=physim(0, simu->mideal, simu->mideal_oa, simu->varol, aster, 0, parms, idtrat, 1, -1);
			real rms=res?P(res,0):simu->varol;
			dfree(res);
#else
			dmat* res=kalman_test(aster->kalman[idtrat], simu->mideal);
			real rms=calc_rms(res, parms->maos.mcc, parms->skyc.evlstart);
			dfree(res);
#endif
		//toc2("estimate");
			P(aster->res_ngs, idtrat, 0)=rms;
			if(parms->skyc.dbg){
				kalman_write(aster->kalman[idtrat], "%s/aster%d_kalman_%d", dirsetup, aster->iaster, dtrat);
			}
		}
		lfree(dtrats);
	}
}
static void select_dtrat(aster_s* aster, int maxdtrat, real maxerror){
	const real wvfmargin=6400e-18; //Allow asterism that is worse by this much to be evaluated.
	aster->idtratest=-1;//idtrat at astermin
	aster->minest=maxerror;
	const int ndtrat=NX(aster->res_ngs);
		for(int idtrat=0; idtrat<ndtrat; idtrat++){
		real wfv=P(aster->res_ngs, idtrat, 0);
		if(wfv<aster->minest){
			aster->idtratest=idtrat;
			aster->minest=wfv;
			}
		}
	if(aster->idtratest!=-1){
		if(maxdtrat>1){
			real thres=MIN(aster->minest+wvfmargin, maxerror);//threshold at low freq end
			real thres2=MIN(aster->minest+wvfmargin, maxerror);//threshold at high freq end
					/*Find upper and skymin good dtrats. */
			for(int idtrat=aster->idtratest; idtrat<ndtrat; idtrat++){
				if(P(aster->res_ngs, idtrat, 0)<thres){
					aster->idtratmax=idtrat+1;
						} else{
							break;
						}
					}
			for(int idtrat=aster->idtratest; idtrat>=0; idtrat--){
				if(P(aster->res_ngs, idtrat, 0)<thres2){
					aster->idtratmin=idtrat;
						} else{
							break;
						}
					}
					//2018-02-28: changed to prefer high frame rate
			if(aster->idtratmax>aster->idtratmin+maxdtrat+1){
				aster->idtratmin=aster->idtratmax-(maxdtrat+1);
			}
		}else{//only evaluate at the best frame rate.
			aster->idtratmin=aster->idtratest;
			aster->idtratmax=aster->idtratmin+1;
					}
	}
}
void setup_aster_controller(sim_s* simu, aster_s* aster, const parms_s* parms){
	if(parms->skyc.servo<0){
		setup_aster_kalman(simu, aster, parms);
	} else{
		setup_aster_servo(simu, aster, parms);
	}
	if(!parms->skyc.multirate){
		select_dtrat(aster, parms->skyc.maxdtrat, parms->skyc.phytype==1?0.5*simu->varol:INFINITY);
			}
			if(parms->skyc.verbose){
		if(!parms->skyc.multirate){
			info("aster%3d stars=(", aster->iaster);
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				info(" %2d", aster->wfs[iwfs].istar);
				}
				info("), dtrats=(%2d,%2d,%2d), res= %.2f nm\n",
				(int)P(parms->skyc.dtrats,aster->idtratmin),
				(int)P(parms->skyc.dtrats,aster->idtratest),
				(int)P(parms->skyc.dtrats,aster->idtratmax-1),
				sqrt(aster->minest)*1e9);
		} else{//multirate
			info("aster%2d, dtrats=(", aster->iaster);
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				info(" %2ld", P(aster->dtrats,iwfs));
				}
				info("), stars=(");
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				info(" %2d", aster->wfs[iwfs].istar);
				}
			info("), est=%.2f nm\n", sqrt(aster->minest)*1e9);
		}
	}
}
/**
   for sort incrementally.
 */
static int sortdbl(const real* a, const real* b){
	return a[0]<b[0]?-1:1;
}
/**
   Select a few asterisms that have decent performance (less than maxerror) */
int setup_aster_select(real* result, aster_s* aster, int naster, star_s* star,
	real maxerror, const parms_s* parms){
	dmat* imin=dnew(2, naster);//record best performance of each asterism.
	int iastermin=-1;//index of asterism giving best answer
	real skymin=INFINITY;//minimum of this field.
	//select the asterism with the best performance
	for(int iaster=0; iaster<naster; iaster++){
		P(imin, 1, iaster)=iaster;
		if((parms->skyc.dbgaster>-1&&iaster!=parms->skyc.dbgaster)||aster[iaster].use==-1){
			P(imin, 0, iaster)=INFINITY;
			continue;
			}
		real astermin=P(imin, 0, iaster)=aster[iaster].minest;//minimum of this asterism.
		if(astermin<skymin){
			iastermin=iaster;
			skymin=astermin;
		}
	}
	if(parms->skyc.dbgsky>-1){
		writebin(imin, "sky%d_imin", parms->skyc.dbgsky);
	}

	result[0]=skymin;
	result[1]=iastermin;
	if(aster[iastermin].idtratest!=-1){
		result[2]=P(parms->skyc.fss, aster[iastermin].idtratest);
	}

	//Mark asterisms and stars to be evaluated further.
	int count=0;
	if(skymin<maxerror){
		int taster=naster;
		if(parms->skyc.maxaster>0&&naster>parms->skyc.maxaster){
			taster=parms->skyc.maxaster;
		}
		qsort(P(imin), naster, 2*sizeof(real), (int(*)(const void*, const void*))sortdbl);
		for(int jaster=0; jaster<taster; jaster++){
			if(P(imin, 0, jaster)>maxerror && naster>3) continue;
			int iaster=(int)P(imin, 1, jaster);
			if(aster[iaster].idtratest==-1) continue;
			count++;
			aster[iaster].use=1;/*mark as valid. */
			char temp1[1024]; temp1[0]='\0';
			for(int iwfs=0; iwfs<aster[iaster].nwfs; iwfs++){
				int istar=aster[iaster].wfs[iwfs].istar;
				int ipowfs=aster[iaster].wfs[iwfs].ipowfs;
				star[istar].use[ipowfs]=1;
				if(parms->skyc.estimate&&parms->skyc.nsky==1){
					char temp2[1023];
					snprintf(temp2, 1023, "(x=%.1f, y=%.1f, J=%.1f) ", star[istar].thetax*RAD2AS, star[istar].thetay*RAD2AS, P(star[istar].mags,0));
					strncat(temp1, temp2, 1023);
				}
			}
			if(parms->skyc.estimate){
				info("%s, %g Hz, %g nm\n", temp1, P(parms->skyc.fss, aster[iaster].idtratest), sqrt(aster[iaster].minest)*1e9);
			}
		}
	}
	if(parms->skyc.verbose){
		info("Minimum is found at aster %g at %.1f Hz: %.2f nm. Will evaluate %d asterisms.\n",
			result[1], result[2], sqrt(result[0])*1e9, count);
	}
	dfree(imin);
	return count;
}
/**
   Free the aster_s array.
 */
void free_aster(aster_s* aster, int naster, const parms_s* parms){
	for(int iaster=0; iaster<naster; iaster++){
		int ndtrat=parms->skyc.ndtrat;
		if(aster[iaster].kalman){
			if(parms->skyc.multirate){
				kalman_free(aster[iaster].kalman[0]);
			} else{
				for(int i=0; i<ndtrat; i++){
					kalman_free(aster[iaster].kalman[i]);
				}
			}
			cellfree(aster[iaster].neam);
			free(aster[iaster].kalman);
			aster[iaster].kalman=0;
		}
		dcellfree(aster[iaster].gain);
		dcellfree(aster[iaster].pgm);
		dcellfree(aster[iaster].sigman);
		dfree(aster[iaster].res_ws);
		dfree(aster[iaster].res_ngs);

		free(aster[iaster].wfs);
		dcellfree(aster[iaster].g);
		dfree(aster[iaster].gm);
		lfree(aster[iaster].dtrats);
		lfree(aster[iaster].idtrats);
		free(aster[iaster].ngs);
		dcellfree(aster[iaster].phyres);
		dcellfree(aster[iaster].phymres);
		lfree(aster[iaster].mdirect);
	}
	free(aster);
}
