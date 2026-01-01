/*
  Copyright 2009-2026 Lianqi Wang
  
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
		aster[iaster].use=1;//set to 0 to disable. upgrade to 2 to enable PO simulation
		int iwfs=0;
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			for(int jwfs=0; jwfs<P(parms->skyc.nwfsmax, ipowfs); jwfs++){
				if(iwfs<nwfstot){
					int istar=P(comb, iwfs, iaster);
					aster[iaster].wfs[iwfs].ipowfs=ipowfs;
					aster[iaster].wfs[iwfs].istar=istar;
					aster[iaster].wfs[iwfs].use=1;
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
 * @brief Mask modes that are not controlled by the loop.
 * 
 * @param parms supplies indps, indastig, indfocus information
 * @param nttf 	number of ttf wfs
 * @param ntt 	number of tt wfs
 * @return lmat* 
 */
static void mode_mask(lmat**pmdirect, const parms_s *parms, int nttf, int ntt){
	if(!pmdirect){
		warning("pmdirect shall be set\n");
		return;
	}
	const int indps=parms->maos.indps;
	const int indastig=parms->maos.indastig;
	const int indfocus=parms->maos.indfocus;
	const int nmod=2+(indps?3:0)+(indastig?2:0)+(indfocus?1:0);
	if(!*pmdirect){
		*pmdirect=lnew(nmod,1);
	}
	lmat *mdirect=*pmdirect;
	if(nttf>0){
		if(nttf+ntt==1 && indps){//single nttf can not control PS1
			P(mdirect, indps)=1;//TTF controls focus instead of magnification plate scale
		}//else: fast controls all modes
	}else{//fast includes only TT wfs. 
		if(indfocus){
			P(mdirect, indfocus)=1;//TT cannot control focus
		}
		if(ntt<3&&indps){//1 or 2 TT OIWFS cannot control ALL PS modes.
			P(mdirect, indps)=1;
		}
		if(ntt<2&&indps){//1 TT cannot control any PS mode.
			P(mdirect, indps+1)=1;
			P(mdirect, indps+2)=1;
		}
		if(ntt<1){
			lset(mdirect, 1);
			warning("there is no fast WFS\n");
		}
	}
}
/**
   Compute Modal to gradient operator by copying from the stars. Using average
gradients. Similar to Z tilt since the mode is low order */
static void setup_aster_g(aster_s* aster, star_s* star, const parms_s* parms){
	aster->g=dcellnew(aster->nwfs, 1);
	aster->ngrad=mycalloc(aster->nwfs, long);
	aster->tsa=0;
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		const int istar=aster->wfs[iwfs].istar;
		const int ipowfs=aster->wfs[iwfs].ipowfs;
		const long nsa=parms->maos.nsa[ipowfs];	
		//aster->g is also used for simulation. Do not zero columns here.
		P(aster->g,iwfs)=ddup(P(star[istar].g,ipowfs));
		aster->tsa+=nsa;
		aster->ngrad[iwfs]=nsa*2;
	}

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
  Estimate wavefront error propagated from measurement error. 
  pgm is the reconstructor. neam is the error covariance.
  trace(Mcc*(pgm*neam*pgm'))
*/
static dmat* calc_recon_error(const dmat* pgm,   /**<[in] the reconstructor*/
	const dmat* neam,/**<[in] the measurement error covariance matrix*/
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
static void mask_gm(dmat *ggm, const lmat *mblind){
	if(!mblind) return;
	for(int imod=0; imod<NY(ggm); imod++){
		if(P(mblind,imod)){
			dzerocol(ggm, imod);
		}
	}
}
/**
 * Select GM using WFS mask. Also compute direct mode mask (mblind) if mblind is set.
 * direct mode are modes that are only sensed by the slower loop and therefore directly fed into the integrator.
 * @param mblind	Modes blind to the current set of WFS
 * @param gm_in		Gradient interaction matrix for each WFS
 * @param mask		Masking out WFS if element is set to 0.
 * @return	dmat	Contatenated Gradient interaction matrix with masked WFS zeroed.
 * */
static dmat* setup_aster_gm_mask(lmat** mblind, const parms_s* parms, const dcell* gm_in, const lmat* mask){
	if(!gm_in||gm_in->nx==0) return NULL;
	dcell* gm=dcellnew(gm_in->nx, 1);
	int nwfs=gm->nx; 
	int nttf=0;
	int ntt=0;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		if(mask&&!P(mask,iwfs)){
			P(gm, iwfs)=dnew(NX(gm_in, iwfs), NY(gm_in, iwfs));
		} else{
			P(gm, iwfs)=dref(P(gm_in, iwfs));
			if(NX(gm_in, iwfs)>2){
				nttf++;
			}else{
				ntt++;
			}
		}
	}
	dmat* ggm=dcell2m(gm); dcellfree(gm);
	if(mblind){
		mode_mask(mblind, parms, nttf, ntt);
	}
	return ggm;
}
//Remove WFS that is marked unused.
static void aster_remove_wfs(aster_s *aster){
	int iwfs=0;
	for(int jwfs=0; jwfs<aster->nwfs; jwfs++){
		if(iwfs!=jwfs){
			aster->wfs[iwfs]=aster->wfs[jwfs];
		}
		if(aster->wfs[iwfs].use){
			iwfs++;
		}
	}
	if(aster->nwfs!=iwfs){
		dbg_once("nwfs reduced from %d to %d\n", aster->nwfs, iwfs);
	}
	aster->nwfs=iwfs;
}
/**
  	For the multirate case, setup the dtrat of each WFS.
	- First, the fastest rate of each wfs is computed based on snr threshold.
	- If skyc.ttffastest is set, the fast rate is forced to be the TTF OIWFS rate. Otherwise, it is set to the fastest OIWFS.
*/
static void
setup_aster_multirate(aster_s* aster, const parms_s* parms){
	const int ndtrat=parms->skyc.ndtrat;
	int idtrat_fast=0;
	//real snr_fast=0;
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		int idtrat;
		for(idtrat=ndtrat-1; idtrat>=0; idtrat--){
			//select highest idtrat that meets snrmin_fast for each WFS.
			if(P(aster->wfs[iwfs].pistat->snr,idtrat)>=P(parms->skyc.snrmin_fast,idtrat)){
				break;
			}
		}
		if(idtrat==-1){
			warning_once("iwfs=%d idtrat=%d\n", iwfs, idtrat);
			continue;
		}
		int isttf=parms->maos.nsa[aster->wfs[iwfs].ipowfs]>1;
		if((parms->skyc.ttffastest&&isttf)||(!parms->skyc.ttffastest &&idtrat>idtrat_fast)){
			idtrat_fast=idtrat;
		}
		//if(parms->skyc.dbg) info("iwfs=%d, idtrat=%d, idtrat_fast=%d\n", iwfs, idtrat, idtrat_fast);
	}
	int idtrat_slow=-1;
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		int idtrat;
		for(idtrat=ndtrat-1; idtrat>=0; idtrat--){
			//select highest idtrat that meets snrmin_slow for each WFS.
			if(P(aster->wfs[iwfs].pistat->snr,idtrat)>=P(parms->skyc.snrmin_slow,idtrat)){
				break;
			}
		}
		if(idtrat==-1){
			warning_once("iwfs=%d idtrat=%d\n", iwfs, idtrat);
			continue;
		}
		if(idtrat_slow==-1){
			idtrat_slow=idtrat;
		}else if(idtrat<idtrat_slow){
			idtrat_slow=idtrat;
		}
		//if(parms->skyc.dbg) info("iwfs=%d, idtrat=%d, idtrat_slow=%d\n", iwfs, idtrat, idtrat_slow);
	}
	int idtrat_fast2=MAX(0, idtrat_fast);//force fast loop to be faster than dtrat=8
	//make sure dtrat_slow is multiple of dtrat_fast2
	int idtrat_slow2=MAX(0, idtrat_slow);
	
	idtrat_slow=idtrat_slow2;
	if(idtrat_slow!=idtrat_fast2){
		while((int)P(parms->skyc.dtrats, idtrat_slow)%(int)P(parms->skyc.dtrats, idtrat_fast2)!=0){
			if(idtrat_slow>0){
				idtrat_slow--;
			}else if(idtrat_fast2>0){
				idtrat_fast2--;
				idtrat_slow=idtrat_slow2;
			}else{
				warning_once("idtrat_slow=%d, idtrat_fast=%d, disable asterism.\n", idtrat_slow, idtrat_fast);
				aster->use=0;
				return;
				break;
			}
		}
	}
	//if(parms->skyc.dbg) info("idtrat_fast2=%d, idtrat_slow=%d\n", idtrat_fast2, idtrat_slow);
	//We use only two rates: a faster rate for tip/tilt(/focus) control and a slower rate for focus/ps control.
	int fast_nttf=0;
	int fast_ntt=0;
	int nwfs_unused=0;
	real snrmin_fast=idtrat_fast==idtrat_slow?P(parms->skyc.snrmin_slow, idtrat_slow):P(parms->skyc.snrmin_fast, idtrat_fast);
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		if(P(aster->wfs[iwfs].pistat->snr,idtrat_fast)>=snrmin_fast){
			aster->wfs[iwfs].idtrat=idtrat_fast2;
			if(parms->maos.nsa[aster->wfs[iwfs].ipowfs]>1){
				fast_nttf++;
			}else{
				fast_ntt++;
			}
		}else if(P(aster->wfs[iwfs].pistat->snr,idtrat_slow)>=P(parms->skyc.snrmin_slow,idtrat_slow)){
			aster->wfs[iwfs].idtrat=idtrat_slow;
		}else{
			aster->wfs[iwfs].use=0;
			aster->wfs[iwfs].idtrat=-1;
			nwfs_unused++;
		}
	}
	if(nwfs_unused){
		aster_remove_wfs(aster);
	}
	aster->idtrats=lnew(aster->nwfs, 1);
	aster->dtrats=lnew(aster->nwfs, 1);
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		P(aster->idtrats, iwfs)=aster->wfs[iwfs].idtrat;
		P(aster->dtrats, iwfs)=P(parms->skyc.dtrats, P(aster->idtrats, iwfs));
	}
	/*if(fast_nttf+fast_ntt<3){
		info("fast_nttf=%d, fast_ntt=%d, nwfs_unused=%d\n", fast_nttf, fast_ntt, nwfs_unused);
	}*/
	mode_mask(&aster->mdirect, parms, fast_nttf, fast_ntt);//a mode is set to 1 if it is directly controlled by the slow mode.
	
	aster->idtratest=idtrat_fast2;
	aster->idtratmin=idtrat_fast2;
	aster->idtratmax=idtrat_fast2+1;
}
/**
 * @brief Compute reconstructor that is insensitive to mblind modes
 * 
 * @param gm 
 * @param nea 
 * @param mblind 
 */
static dmat* recon_mblind(const dmat* gm, const dmat* nea, const lmat *mblind){
	//old method:
	dmat *gmzero;
	if(mblind){
		gmzero=ddup(gm);
		mask_gm(gmzero, mblind);
	}else{
		gmzero=dref(gm);
	}
	//dshow(gmzero, "gmzero");
	dmat *res=dpinv(gmzero, nea);
	dfree(gmzero);
	return res;
}
/**

   Setup the least squares reconstructor and controller.

   It first computes the reconstructor, and noise propagation, and then optimize
   servo gains to reduce the combined noise propagation and residual error.

   For multirate controller, it sets out multiple (two for now) reconstructors
   for different dtrats based on the number of participating WFS. The slower
   loop provides complete set of reconstruction. The modes that are controlled
   by faster loop uses the slower loop integration output as an offset.

   The gain is optimized to minimize
   \f[ \sigma^2=\int \textrm{PSD}_{ngs,ws}H_{rej}\textrm{d}\nu +
   \int_0^{\nu_{nyquist}} \textrm{PSF}\textrm{d}\nu \f]
*/
static void setup_aster_servo(sim_s* simu, aster_s* aster, const parms_s* parms){
	const int multirate=parms->skyc.multirate;
	int ncase=0;
	int fast_idtrat=0;//for multirate fast loop
	int slow_idtrat=0;//for multirate slow loop
	const int icase_fast=0;//fast loop index into case_valid, pgm, gain, etc.
	if(multirate){
		slow_idtrat=parms->skyc.ndtrat;
		//Find the dtrat of the slowest WFS
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
			int idtrat=P(aster->idtrats,iwfs);
			if(idtrat!=-1){
				if(slow_idtrat>idtrat){
					slow_idtrat=idtrat;
				}
				if(fast_idtrat<idtrat){
					fast_idtrat=idtrat;
				}
			}
		}
		ncase=slow_idtrat==fast_idtrat?1:2;
	} else{
		ncase=parms->skyc.ndtrat;
	}
	const int icase_slow=ncase-1;//slow loop is always the last entry when all WFS participate.
	if(aster->gain){
		dcellfree(aster->pgm);
		dcellfree(aster->sigman);
		dcellfree(aster->gain);
		dfree(aster->res_ws);
		dfree(aster->res_ngs);
	}
	aster->pgm=dcellnew(ncase, 1);
	aster->sigman=dcellnew(ncase, 1);
	aster->gain=dcellnew(ncase, 1);
	aster->res_ws=dnew(multirate?1:ncase, 1);
	aster->res_ngs=dnew(multirate?1:ncase, 3);
	lmat* mblind=NULL;
	dmat* gmfull=setup_aster_gm_mask(&mblind, parms, aster->g, NULL);
	dmat* gmfast=NULL;
	const long nmod=parms->maos.nmod;
	
	dcell* nea=dcellnew3(aster->nwfs, 1, aster->ngrad, 0);
	lmat* fast_mask=multirate?lnew(aster->nwfs, 1):0;
	real res_ngs=0;/*residual error due to signal after servo rejection. */
	real res_ngsn=0;/*residual error due to noise. */
	for(int icase=0; icase<ncase; icase++){
		//We need to build the reconstructor for each valid case.
		//Assemble the measurement error covariance matrix
		int idtrat;
		if(multirate){//two cases for multirate. fast and slow loop
			idtrat=icase==icase_fast?fast_idtrat:slow_idtrat;
			if(icase==icase_fast){//fast loop may have subset of WFS
				for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
					if(P(aster->idtrats, iwfs)==idtrat){
						P(fast_mask, iwfs)=1;//has output for fast loop
					}else{
						P(fast_mask, iwfs)=0;//does not have output
					}
				}
			}
		}else{
			idtrat=icase;//for non multirate
		}
		//Since we masked gm, no need to mask NEA.
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
			dcp(&P(nea,iwfs), P(aster->wfs[iwfs].pistat->sanea, idtrat));
			dcwpow(P(nea,iwfs), -2);//in radian^-2
		}
		//dshow(nea->m, "nea");
		if(multirate){
			if(icase==icase_fast){//fast loop
				if(parms->skyc.dbg) lshow(aster->mdirect, "mdirect");
				//zero out GM for modes not measured by the faster loop. 
				gmfast=setup_aster_gm_mask(&aster->mdirect, parms, aster->g, fast_mask);
				P(aster->pgm,icase)=recon_mblind(gmfast, nea->m, aster->mdirect);
				if(ncase>=1){//do not control highly aliased modes in the fast loop when there is a slow loop.
					dmat *aliased=NULL;
					real factor=ncase>1?1:5;
					int nr=0;
					dmm(&aliased, 0, P(aster->pgm, icase), P(aster->pgm, icase), "nt", 1);
					if(parms->skyc.dbg) dshow(aliased, "aliased");
					if(parms->maos.indps){
						real sum1=0,sum2=0;
						for(int i=parms->maos.indps; i<parms->maos.indps+3; i++){
							if(P(aster->mdirect, i)) continue;
							for(int j=0; j<i; j++){
								sum1+=fabs(P(aliased, j, i));
							}
							sum2+=fabs(P(aliased, i, i));
						}
						if(parms->skyc.dbg) info("sum1=%g, sum2=%g, ratio=%g\n", sum1, sum2, sum1/sum2);
						if(sum1>sum2*factor){
							for(int i=parms->maos.indps; i<parms->maos.indps+3; i++){
								if(!P(aster->mdirect, i)){
									P(aster->mdirect, i)=1;
									nr++;
									//warning_once("mode %d is zeroed due to aliasing\n", i);
								}
							}
						}
					}
					dfree(aliased);
					if(nr){
						if(parms->skyc.dbg) lshow(aster->mdirect, "mdirect_aliased");
						dfree(P(aster->pgm,icase));
						P(aster->pgm,icase)=recon_mblind(gmfast, nea->m, aster->mdirect);
						//dshow(gmfast, "gmfast");
						//dshow(P(aster->pgm, icase), "pgmfast");
					}
				}
				
			}else{//slow loop
				//Setup modes controlled by slow loop that is not measured correctly by the fast loop
				//Mslow=(I-Gfast*Gfast^-1)
				dmat *Mslow=dnew(nmod, nmod); 
				daddI(Mslow, 1);
				dmm(&Mslow, 1, P(aster->pgm,icase_fast), gmfast, "nn", -1);
				if(parms->skyc.dbg){
					dshow(Mslow, "aster%d_M_slow", aster->iaster);
					//dshow(Mslow, "aster%d_M_slow", aster->iaster);
					//dshow(gmfast, "gm_fast");
					//dshow(gmfull, "gmfull");
					//dshow(pgmslow, "pgmslow");
					//dshow(P(aster->pgm, icase_fast), "pgm_fast");
				}
				dmat *gmslow=NULL;
				dmm(&gmslow, 0, gmfull, Mslow, "nn", 1);
				dmat *pgmslow=dpinv(gmslow, nea->m);
				dmm(&P(aster->pgm,icase), 0, Mslow, pgmslow, "nn", 1);
				dfree(Mslow);
				dfree(gmslow);
				dfree(pgmslow);
			}
		}else{
			P(aster->pgm,icase)=recon_mblind(gmfull, nea->m, mblind);
		}
		
		//Compute noise propagation
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
			dcwpow(P(nea,iwfs), -1);//in radian^2.
		}
		P(aster->sigman,icase)=calc_recon_error(P(aster->pgm,icase), nea->m, parms->maos.mcc);
		if(parms->skyc.dbg){
			dshow(P(aster->sigman, icase), "aster%d_sigman_%.0f", aster->iaster, P(parms->skyc.dtrats, idtrat));
		}
		if(parms->skyc.dbg>2){
			writebin(nea, "%s/aster%d_nea%d", dirsetup, aster->iaster, icase);
		}
		//look gain optimization.
		/*gsplit:
		  0: All modes use the same gain.
		  1: Different mode use different gain (2017-0-24) was only tt/ps separate.
		  note that P(simu->psds,0) contains windshake PSD.
		*/

		const int servotype=parms->skyc.servo==2?2:1;
		const int ng=parms->skyc.ngain;//number of gain parameters
		dmat* pgain=P(aster->gain,icase)=dnew(ng, nmod);
		dmat* sigma=dnew(1,1);
		for(int ipsd=0; ipsd<simu->psds->nx; ipsd++){//npsd=1 is gsplit=0
			P(sigma, 0)=P(P(aster->sigman,icase),parms->skyc.gsplit?ipsd:nmod);
			if(!P(sigma, 0)) continue;//not used
			real pg[ng+2];
			//For the slow loop, not the FULL PSD is in affect. How do we reduce the gain accordingly?
			if(parms->skyc.interpg){//LUT with noise level use pre-determined gain 
				interp_gain(pg, P(P(simu->gain_pre,idtrat),ipsd), simu->gain_x, P(sigma, 0));
			} else{
				dcell *tmp=servo_optim(parms->maos.dt, P(parms->skyc.dtrats, idtrat), 0, parms->skyc.pmargin, 0, 0, servotype, P(simu->psds, ipsd), sigma);
				memcpy(pg, P(P(tmp,0)), (ng+2)*sizeof(real));
				dcellfree(tmp);
			}
			if(multirate){
				if(icase==icase_fast){//this is faster loop
					if(!P(aster->mdirect, ipsd)){//mode driven only by fast loop
						res_ngs+=pg[ng];
						res_ngsn+=pg[ng+1];
					}
				}else{//slow loop
					if(!P(aster->mdirect, ipsd)){//mainly driven by fast loop
						//slower loop is not expriencing the full PSD. reduce the gain.
						pg[0]*=0.5;
						pg[ng+1]*=0.5;
						res_ngsn+=pg[ng+1];//only count noise propation
					}else{//mainly driven by slower loop
						res_ngs+=pg[ng];
						res_ngsn+=pg[ng+1];
					}
				}
			}else{
				res_ngs+=pg[ng];
				res_ngsn+=pg[ng+1];
			}
			memcpy(PCOL(pgain, ipsd), pg, sizeof(real)*ng);
		}//for ipsd
		if(parms->skyc.dbgaster!=-1){
			info("aster %d idtrat %d, est %6.1f signal %6.1f noise %6.1f nm\n", aster->iaster, idtrat, sqrt(res_ngs+res_ngsn)*1e9, sqrt(res_ngs)*1e9, sqrt(res_ngsn)*1e9);
		}
		dfree(sigma);
		if(simu->psds->nx < nmod){
			warning_once("Use gain of mode 0 for modes with no PSD\n");
			for(int imod=simu->psds->nx; imod<nmod; imod++){
				memcpy(PCOL(pgain, imod), PCOL(pgain, 0), sizeof(real)*ng);
			}
		}
		if(!multirate){
			P(aster->res_ngs, icase, 0)=res_ngs+res_ngsn;/*error due to signal and noise */
			P(aster->res_ngs, icase, 1)=res_ngs;/*error due to signal */
			P(aster->res_ngs, icase, 2)=res_ngsn;/*error due to noise propagation. */
			res_ngs=0;//reset for each dtrat
			res_ngsn=0;//reset for each dtrat
		}
		if(parms->skyc.psd_ws && !parms->skyc.addws && (!multirate || icase==icase_fast)){
			warning_once("Determine residual wind shake using servo analysis.\n");
			dmat* g_tt=drefcols(pgain, 0, 1);
			real gain_n;
			P(aster->res_ws,multirate?0:icase)=servo_residual(&gain_n, parms->skyc.psd_ws,
				parms->maos.dt, P(parms->skyc.dtrats,idtrat), 0, g_tt, servotype);
			dfree(g_tt);
		}
		if(parms->skyc.dbgaster!=-1){
			//dshow(gmfull, "gmfull_%d", icase);
			//dshow(P(aster->pgm,icase), "pgm_%d", icase);
			dshow(P(aster->gain, icase), "gain_%.0f", P(parms->skyc.dtrats, idtrat));
			//dshow(P(aster->sigman, icase), "sigman_%d", icase);
			//dshow(nea->m, "nea_%d", icase);
		}
	}//for icase
	dcellfree(nea);
	lfree(mblind);
	lfree(fast_mask);
	if(multirate){//for setup_aster_select() to use.
		P(aster->res_ngs, 0, 0)=res_ngs+res_ngsn;/*error due to signal and noise */
		P(aster->res_ngs, 0, 1)=res_ngs;/*error due to signal */
		P(aster->res_ngs, 0, 2)=res_ngsn;/*error due to noise propagation. */
		
		aster->minest=P(aster->res_ngs, 0, 0);

		if(icase_fast!=icase_slow){//there is slow loop. make gains equal.
			//Average the gain for different modes so that the control output is invisible to fast loop.
			int n=nmod;
			if(parms->maos.indfocus==n-1 && P(aster->mdirect, n-1)){//focus is by the slow loop and is already blind to the fast loop
				n--;
			}
			real gsum=0;
			//real gmax=0;
			//real gmin=1;
			real gct=0;
			for(int i=0; i<n; i++){
				if(P(aster->mdirect, i)){
					gsum+=P(P(aster->gain, icase_slow), i);
					//gmax=MAX(gmax, P(P(aster->gain, icase_slow), i));
					//gmin=MIN(gmin, P(P(aster->gain, icase_slow), i));
					gct+=1;
				}
			}
			if(gct>0){
				gsum/=gct;
				for(int i=0; i<n; i++){
					P(P(aster->gain, icase_slow), i)=gsum;
				}
			}
			if(parms->skyc.dbg){
				dshow(P(aster->gain, icase_slow), "gain_slow");
			}
		}
	}
	if(parms->skyc.dbg>3){
		writebin(gmfull, "%s/aster%d_gm", dirsetup, aster->iaster);
		writebin(aster->pgm, "%s/aster%d_pgm", dirsetup, aster->iaster);
		writebin(aster->sigman, "%s/aster%d_sigman", dirsetup, aster->iaster);
		writebin(aster->gain, "%s/aster%d_gain", dirsetup, aster->iaster);
		if(parms->skyc.psd_ws){
			writebin(aster->res_ws, "%s/aster%d_res_ws", dirsetup, aster->iaster);
		}
		writebin(aster->res_ngs, "%s/aster%d_res_ngs", dirsetup, aster->iaster);
	}
	dfree(gmfull);
	dfree(gmfast);
}
/**Place sanea**2 to diagonal of neam */
static void set_diag_pow2(dmat **pneam, dmat *sanea){
	int ng=PN(sanea);
	dinit(pneam, ng, ng);
	for(int ig=0; ig<ng; ig++){
		P(*pneam,ig,ig)=pow(P(sanea,ig), 2);
	}
}

static void setup_aster_kalman(sim_s* simu, aster_s* aster, const parms_s* parms){
	int ndtrat=parms->skyc.multirate?1:parms->skyc.ndtrat;
	if(!aster->res_ngs){
		aster->res_ngs=dnew(ndtrat, 3);
	}
	aster->kalman=mycalloc(ndtrat, kalman_t*);
	lmat* dtrats=parms->skyc.multirate?aster->dtrats:lnew(aster->nwfs, 1);	
	dcell *neam=NULL; 
	int skip=0;
	for(int idtrat=0; idtrat<ndtrat; idtrat++){
		if(!parms->skyc.multirate){
			if(!neam) neam=dcellnew(aster->nwfs, 1);
			lset(dtrats, P(parms->skyc.dtrats,idtrat));
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				set_diag_pow2(&P(neam, iwfs), P(aster->wfs[iwfs].pistat->sanea, idtrat));
			}
		}else{
			lmat *idtrats=lunique(aster->idtrats, 0);
			if(!neam) neam=dcellnew(aster->nwfs, PN(idtrats));
			for(int i=0; i<PN(idtrats); i++){
				for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
					set_diag_pow2(&P(neam, iwfs, i), P(aster->wfs[iwfs].pistat->sanea, P(idtrats, i)));
				}
			}
			lfree(idtrats);
		}
		if(!skip){
			aster->kalman[idtrat]=sde_kalman(simu->sdecoeff, parms->maos.dt, dtrats, aster->mdirect, aster->g, neam, 0);
		}
		if(!aster->kalman[idtrat]){//failed to converge
			P(aster->res_ngs, idtrat, 0)=simu->varol;
			aster->use=0;
		}else{//in multirate, use servo result
			dmat *tmp=0;
			dmm(&tmp, 0, P(aster->kalman[idtrat]->P, 0), parms->maos.mcc, "nn", 1);
			real rms=dtrace(tmp);
			dfree(tmp);
			P(aster->res_ngs, idtrat, 0)=rms;
		}
		if(parms->skyc.dbg){
			kalman_write(aster->kalman[idtrat], "%s/aster%d_kalman_%ld", dirsetup, aster->iaster, lmin(dtrats));
		}
	}
	dcellfree(neam);
	if(!parms->skyc.multirate){
		lfree(dtrats);
	}else{
		//aster->minest=P(aster->res_ngs,0,0);
	}
}
static void select_aster_dtrat(aster_s* aster, int maxdtrat, real maxerror){
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
		if(maxdtrat>1 && ndtrat>1){
			real thres=MIN(aster->minest*2+wvfmargin, maxerror);//threshold
			/*Find upper and skymin good dtrats. */
			for(int idtrat=aster->idtratest; idtrat<ndtrat; idtrat++){
				if(P(aster->res_ngs, idtrat, 0)<thres){
					aster->idtratmax=idtrat+1;
				} else{
					break;
				}
			}
			for(int idtrat=aster->idtratest; idtrat>=0; idtrat--){
				if(P(aster->res_ngs, idtrat, 0)<thres){
					aster->idtratmin=idtrat;
				} else{
					break;
				}
			}
			//2018-02-28: changed to prefer high frame rate
			//if(aster->idtratmax>aster->idtratmin+maxdtrat+1){
				//aster->idtratmin=aster->idtratmax-(maxdtrat+1);
			//}
		}else{//only evaluate at the best frame rate.
			aster->idtratmin=aster->idtratest;
			aster->idtratmax=aster->idtratmin+1;
		}
	}
}
void setup_aster_controller(sim_s* simu, aster_s* aster, star_s *star, const parms_s* parms){
	if(parms->skyc.multirate){
		setup_aster_multirate(aster, parms);//this may alter the number of wfs.
	}
	setup_aster_g(aster, star, parms);//after setup_aster_multirate as it may alter the number of wfs.
	if(!aster->use){
		aster->minest=simu->varol;
		return;
	}
	if(parms->skyc.servo>0 || parms->skyc.servo==-2){
		setup_aster_servo(simu, aster, parms);//may use this in LQG multirate mode.
	}
	if(parms->skyc.servo<0){
		setup_aster_kalman(simu, aster, parms);
	}
	if(!parms->skyc.multirate){
		select_aster_dtrat(aster, parms->skyc.maxdtrat, simu->varol);
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
		if((parms->skyc.dbgaster>-1&&iaster!=parms->skyc.dbgaster)||!aster[iaster].use){
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
		if(!parms->skyc.multirate && parms->skyc.maxaster>0&&naster>parms->skyc.maxaster){
			taster=parms->skyc.maxaster;
		}
		qsort(P(imin), naster, 2*sizeof(real), (int(*)(const void*, const void*))sortdbl);
		for(int jaster=0; jaster<taster; jaster++){
			if(P(imin, 0, jaster)>maxerror && naster>3) continue;
			int iaster=(int)P(imin, 1, jaster);
			if(aster[iaster].idtratest==-1) continue;
			count++;
			aster[iaster].use=2;/*mark as valid for physical optics simulation. */
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
		lfree(aster[iaster].dtrats);
		lfree(aster[iaster].idtrats);
		free(aster[iaster].ngrad);
		dcellfree(aster[iaster].phyres);
		dcellfree(aster[iaster].phymres);
		lfree(aster[iaster].mdirect);
	}
	free(aster);
}
