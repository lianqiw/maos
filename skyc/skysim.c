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
   \file skysim.c

   Contains the skysim() function that does the skycoverage simulation.

   The steps:
   *) Read in parameters
   *) Read in necessary struct saved by MAOS and/or prepared by user
   *) Read in star catalog.
   *) Read in Focus tracking error as a function of sampling frequency
   *) Read in NGS Mode time vector
   *) For each possible asterism in a star field.
   #) Compute NGS Gradient operator from NGS modes.
   #) Read in saved PSF and interpolate to the location of NGS. apply extra strehl degradation.
   #) Compute NGS signal level and measurement error from matched filter (be careful about PSF quality).
   #) For each sample frequency.
   $) Compute the reconstructor
   $) Optimize for the loop gain.
   $) Estimate residual windshake error.
   $) Do a rough open loop estimate of the total error. (atm + windshake +focus tracking ...)
   $) Guess for the best sampling frequency (compare with closed loop results to validate this model).
   *) Compare the asterisms and pick only a few that are within ~200% of the best asterism at their best sampling freq.
   *) For each remaining asterism in the star field.
   #) Do a temporal simulation to obtain the final result.
   #) save error due to focus tracking, telescope windshake, separately.
   #) Find the best sampling frequency.
   *) Find the best asterism.
   *) Repeat for another star field.


   Don't forget:

   *) Sampling frequency dependent rne. input a precomputed matrix for all
   possible sampling frequencies.
   *) Implementation error.
   *) Windshake (handle separately), focus tracking jitter
   *) On axis performance.

   2010-06-16: Single asterism debugged.
   2010-06-21: Bug found: wfs.thetax, thetay not copied from star.
   differences:
   MATLAB starts with a guess using QC model for NEA when regenerating pistat
   maos starts with loaded pistat. This is 1 difference.
   The gains are not sensitive to which PSD is loaded, or whether PSD is scaled or not.

   2011-01-17: Do not parallelize across different star fields, that takes too
   much memory/bandwidth to load PSF.  We parallelize across asterism and dtrats instead.

   2012-12-14 note about focus mode.
   When maos.nmod=5, global focus residual is estimated using PSD.
   When maos.nmod>5, global focus residual is not estimated using PSD.
   when skyc.addfocus==1, add focus to time series.
   when skyc.ahstfocus==1, use new 5+1 mode.
   if maos.mffocus==0, set skyc.addfocus to 0 because focus error is included in maos simulations.
   2013-12-20: removed skyc.addfocus option. The correct way to model focus is to use in maos presimulations
*/
#include <unistd.h>
#include "skyc.h"
#include "parms.h"
#include "skysim.h"
#include "types.h"
#include "powfs.h"
#include "aster.h"
#include "physim.h"
#include "mtch.h"
#include "genstars.h"
#include "genpistat.h"
#include "star.h"
#include "utils.h"
#include "nafocus.h"
static void res_set_all(dmat *res, int isky, real varol){
	P(res, 0, isky)=varol;
	P(res, 1, isky)=varol;
	P(res, 2, isky)=varol*0.707;
	P(res, 4, isky)=varol;
	P(res, 5, isky)=varol;
	P(res, 6, isky)=P(res, 2, isky);
	P(res, 8, isky)=-1;
}
/**
   The actual work horse that does the physical optics time domain simulation.
*/
static void* skysim_isky(sim_s* simu){
	unsigned int isky;
	powfs_s* powfs=simu->powfs;
	const parms_s* parms=simu->parms;
	const dcell* stars=simu->stars;
	const int noisy=parms->skyc.noisy;
	const int seed_maos=simu->seed_maos;
	
	int nstep=parms->maos.nstep;

	while((isky=atomic_fetch_add(&simu->isky, 1))<simu->isky_end){
		real tk_1=myclockd();
		/*Setup star parameters. */
		int nstar;
		star_s* star=setup_star(&nstar, simu, P(stars,isky), seed_maos);
		int naster;
		aster_s* asters=setup_aster_comb(&naster, star, nstar, parms);
		if(!asters||!naster){
			//info("Field%4d,  0 stars\n", isky);
			res_set_all(simu->res, isky, simu->varol);
			continue;
		}
		if(parms->skyc.dbgaster>=naster){
			warning("skyc.dbgaster=%d, but naster=%d\n", parms->skyc.dbgaster, naster);
		}
		/*
		  We first estimate the matched filter, reconstructor, and servo
		  loop optimization to determine the approximate wavefront error. Only
		  a few combinations are kept for each star field for further time
		  domain simulations.
		*/
	OMP_TASK_FOR(4)
		for(int iaster=0; iaster<naster; iaster++){
			if(parms->skyc.dbgaster<0||iaster==parms->skyc.dbgaster){
				/*Parallelizing over aster gives same random stream. */
				seed_rand(&asters[iaster].rand, parms->skyc.seed+iaster+isky+40);
				/*Compute the reconstructor, nea, sigman and optimize controller. */
				setup_aster_controller(simu, &asters[iaster], star, parms);
			}
		}
		/*Select asters that have good performance. */
		real res_geom[3]={0,0,0};
		setup_aster_select(res_geom, asters, naster, star, simu->varol, parms);
		real sky_min=simu->varol;//best result for this sky field
		int sky_iaster=-1;//selected asterism index for this sky field
		int sky_idtrat=-1;//idtrat for selected asterism for this sky field
		
		if(!parms->skyc.estimate){
			/*Read in physical optics data (wvf) */
			setup_star_read_ztilt(star, nstar, parms, seed_maos);
			if(parms->skyc.phystart>=0 && parms->skyc.phystart<parms->maos.nstep){
				nstep=setup_star_read_wvf(star, nstar, parms, seed_maos);
			}
			//Now begin time domain Physical Optics Simulations.
OMP_TASK_FOR(4)
			for(int iaster=0; iaster<naster; iaster++){
				aster_s* asteri=&asters[iaster];
				asteri->nstep=nstep;
				asteri->minphyres=simu->varol;
				if(asteri->use!=2){
					goto skip1;
				}

				/*if(parms->skyc.verbose){
					info("Aster %d, estimate is %.1fnm at %.0fhz. ", iaster,sqrt(asteri->minest)*1e9, P(parms->skyc.fss, asteri->idtratest));
					if(!parms->skyc.multirate){
						info("Try %.1f to %.1f Hz. \n",P(parms->skyc.fss, asteri->idtratmin), P(parms->skyc.fss, asteri->idtratmax-1));
					}
				}*/
				setup_aster_ztilt(asteri, star);
				if(parms->skyc.phystart>=0 && parms->skyc.phystart<parms->maos.nstep){
					setup_aster_wvf(asteri, star);//Assign wvf from star to aster 
				}
				/*Compute the reconstructor, nea, sigman and optimize controller again. no need redo?*/
				//setup_aster_controller(simu, asteri, parms);
				asteri->phyres=dcellnew(asteri->idtratmax, 1);
				asteri->phymres=dcellnew(asteri->idtratmax, 1);
				
				for(int idtrat=asteri->idtratmin; idtrat<asteri->idtratmax; idtrat++){
					P(asteri->phyres,idtrat)=physim(&P(asteri->phymres,idtrat),
						simu->mideal, simu->mideal_oa, simu->varol,
						asteri, powfs, parms, idtrat, noisy, parms->skyc.phystart);
					if(P(asteri->phyres,idtrat)){
						if(!parms->skyc.addws && parms->skyc.psd_ws){
							P(P(asteri->phyres,idtrat), 0)+=P(asteri->res_ws,idtrat);
						}
						if(P(P(asteri->phyres,idtrat), 0)<asteri->minphyres){
							asteri->minphyres=P(P(asteri->phyres,idtrat), 0);
							asteri->idtratphy=idtrat;
						}
					}
				}
				if(asteri->minphyres<simu->varol && simu->res_aster){//collect statistics of all asterisms
					int res_iaster=atomic_fetch_add(&simu->res_iaster, 1);
					asteri->iaster_all=res_iaster;
					if(res_iaster+20>NY(simu->res_aster)){
						OMP(critical) //inside if test to avoid overhead
						if(res_iaster+20>NY(simu->res_aster)){
							dresize(simu->res_aster, NX(simu->res_aster), NY(simu->res_aster)+1000);
						}
					}
					real *p=PCOL(simu->res_aster, res_iaster);
					*(p++)=asteri->minphyres;//phy result
					*(p++)=asteri->minest;//estimate
					*(p++)=isky;
					*(p++)=asteri->nwfs;
					for(int iwfs=0; iwfs<asteri->nwfs; iwfs++){
						const int istar=asteri->wfs[iwfs].istar;
						const int idtrat=parms->skyc.multirate?P(asteri->idtrats, iwfs):asteri->idtratphy;
						*(p++)=idtrat==-1?0:P(parms->skyc.fss, idtrat);
						*(p++)=idtrat==-1?0:P(asteri->wfs[iwfs].pistat->snr, idtrat);
						*(p++)=star[istar].thetax;
						*(p++)=star[istar].thetay;
						for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
							*(p++)=P(star[istar].mags, iwvl);
						}
					}
			
				}
skip1:;
			}//for iaster
			for(int iaster=0; iaster<naster; iaster++){
				aster_s* asteri=&asters[iaster];
				//find the best asterism
				if(asteri->minphyres<sky_min||parms->skyc.dbgaster==iaster){
					sky_iaster=iaster;
					sky_min=asteri->minphyres;
				}
				if(parms->skyc.verbose){//print outside of parallel region
					for(int idtrat=asteri->idtratmin; idtrat<asteri->idtratmax; idtrat++){
						int sstmp[3]={-1,-1,-1};
						//real sxtmp[3]={0,0,0};
						//real sytmp[3]={0,0,0};
						int dttmp[3]={0,0,0};
						int nwfs=MIN(3, asteri->nwfs);
						real snrtmp[3]={0,0,0};
						for(int iwfs=0; iwfs<nwfs; iwfs++){
							sstmp[iwfs]=asteri->wfs[iwfs].istar;
							//sxtmp[iwfs]=star[asteri->wfs[iwfs].istar].thetax*206265;
							//sytmp[iwfs]=star[asteri->wfs[iwfs].istar].thetay*206265;
						}
						for(int iwfs=0; iwfs<nwfs; iwfs++){
							if(!asteri->wfs[iwfs].use) continue;
							int idtmp=parms->skyc.multirate?P(asteri->idtrats, iwfs):idtrat;
							dttmp[iwfs]=P(parms->skyc.dtrats, idtmp);
							snrtmp[iwfs]=P(asteri->wfs[iwfs].pistat->snr, idtmp);
						}
						real phyps=INFINITY;
						real phytt=INFINITY;
						real phyres=INFINITY;
						const real addws=(parms->skyc.addws||!parms->skyc.psd_ws)?0:P(asteri->res_ws, idtrat);
						if(asteri->phyres&&P(asteri->phyres,idtrat)){
							phyres=P(P(asteri->phyres,idtrat), 0)+addws;
							phytt=P(P(asteri->phyres,idtrat), 1)+addws;
							phyps=P(P(asteri->phyres,idtrat), 0)-P(P(asteri->phyres,idtrat), 1);
						}
						info("Ast%3d: dtrat=(%2d %2d %2d) Stars=(%d %d %d) SNR=(%4.1f %4.1f %4.1f) Est %6.1f Phy %6.1f TT %6.1f PS%6.1fnm\n",
							asteri->iaster, dttmp[0], dttmp[1], dttmp[2], sstmp[0], sstmp[1], sstmp[2],
							snrtmp[0], snrtmp[1], snrtmp[2],
							sqrt(P(asteri->res_ngs, parms->skyc.multirate?0:idtrat))*1e9, 
							sqrt(phyres)*1e9, sqrt(phytt)*1e9, sqrt(phyps)*1e9);
					}
				}

			}
			if(sky_iaster!=-1 && asters[sky_iaster].phyres){
				aster_s* asteri=&asters[sky_iaster];
				sky_idtrat=asteri->idtratphy;
				dmat* phymin=P(asteri->phyres, sky_idtrat);
				if(phymin && P(phymin, 0)<simu->varol){
					/*Field Averaged Performance. */
					//May include ws/vib if skyc.addws is true.
					const real addws=(parms->skyc.addws||!parms->skyc.psd_ws)?0:P(asteri->res_ws, sky_idtrat);
					P(simu->res, 0, isky)=P(phymin,0)+addws;/*Total */
					P(simu->res, 1, isky)=P(phymin,0);/*ATM NGS Mode error. */
					P(simu->res, 2, isky)=P(phymin,1);/*ATM Tip/tilt Error. */
					P(simu->res, 3, isky)=addws;/*Residual wind shake TT*/
					P(simu->res, 4, isky)=asteri->minest;/*estimated error, changed on 12/17/2021*/
					P(simu->res, 5, isky)=P(phymin,2);/*ATM NGS Mode error for on axis only. */
					P(simu->res, 6, isky)=P(phymin,3);/*ATM Tip/tilt Error for on axis only. */
					P(simu->res, 7, isky)=P(parms->skyc.fss, sky_idtrat);
					P(simu->res, 8, isky)=asteri->iaster_all;//index to res_aster
					P(simu->res, 9, isky)=sky_iaster;
					if(parms->skyc.save){
						skysim_save(simu, asteri, PCOL(simu->res, isky), sky_iaster, sky_idtrat, isky);
					}
				}
				dcp(&P(simu->mres,isky), P(asteri->phymres, sky_idtrat));
			}else{
				res_set_all(simu->res, isky, simu->varol);
			}
		} else{//If(skyc.estimate)
			sky_min=res_geom[0];
			sky_iaster=res_geom[1];
			sky_idtrat=asters[sky_iaster].idtratest;
			res_set_all(simu->res, isky, sky_min);
			P(simu->res, 7, isky)=P(parms->skyc.fss, sky_idtrat);
			P(simu->res, 9, isky)=sky_iaster;
		}//If(skyc.estimate)
		if(sky_iaster!=-1){
			if(parms->skyc.servo>0&&!parms->skyc.multirate&&sky_idtrat!=-1){
				dcp(&P(simu->gain,isky), P(asters[sky_iaster].gain,sky_idtrat));
			}
		
			real tk_4=myclockd();
			LOCK(simu->mutex_status);
			simu->status->isim=isky;
			simu->status->tot=tk_4-tk_1;/*per step */
			simu->status->laps=tk_4-simu->tk_0;
			int nsky_tot=simu->isky_end-simu->isky_start;
			int nsky_left=simu->isky_end-simu->isky-1+nsky_tot*(parms->maos.nseed-simu->iseed-1);
			int nsky_laps=simu->isky-simu->isky_start+1+nsky_tot*simu->iseed;
			simu->status->rest=simu->status->laps*nsky_left/nsky_laps;
			simu->status->clerrlo=sqrt(P(simu->res, 2, isky))*1e9;
			simu->status->clerrhi=sqrt(P(simu->res, 0, isky))*1e9;
			scheduler_report(simu->status);
			//long totm=(long)floor(simu->status->tot/60.);
			//long tots=(long)simu->status->tot-totm*60;
			long laps_h=simu->status->laps/3600;
			long laps_m=simu->status->laps/60-laps_h*60;
			long rest_h=simu->status->rest/3600;
			long rest_m=simu->status->rest/60-rest_h*60;
			
			while(simu->isky_print<simu->isky_end && P(simu->res, 0, simu->isky_print)){
				real freqs[3]={0,0,0};
				real snrs[3]={0,0,0};
				if(P(simu->res, 8, simu->isky_print)!=-1){
					real *pall=PCOL(simu->res_aster, (int)P(simu->res, 8, simu->isky_print));
					for(int i=0; i<MIN(3,pall[3]); i++){
						real freq=(pall[4+(4+parms->maos.nwvl)*i]*parms->maos.dt);
						freqs[i]=freq?1/freq:0;
						snrs[i]=pall[5+(4+parms->maos.nwvl)*i];
					}
				}
				int sky_iaster2=P(simu->res, 9, simu->isky_print);
				info("Field%4d %3d: %2.0f %2.0f %2.0f SNR %5.1f %5.1f %5.1f Est %6.1f Phy %6.1f TT%4.0f PS%4.0f nm. Used%2ld:%02ld Left%2ld:%02ld\n",
					simu->isky_print, sky_iaster2, freqs[0], freqs[1], freqs[2], snrs[0], snrs[1], snrs[2],
					sqrt(P(simu->res, 4, simu->isky_print))*1e9,
					sqrt(P(simu->res, 0, simu->isky_print))*1e9,
					sqrt(P(simu->res, 2, simu->isky_print)+P(simu->res, 3, simu->isky_print))*1e9,
					sqrt(P(simu->res, 1, simu->isky_print)-P(simu->res, 2, simu->isky_print))*1e9, 
					laps_h, laps_m, rest_h, rest_m);
				simu->isky_print++;
			}
			UNLOCK(simu->mutex_status);
		}
		free_aster(asters, naster, parms);
		free_star(star, nstar, parms);
	}/*while isky*/
	return NULL;
}

/**
   Read in ideal NGS modes
*/
static void skysim_read_mideal(sim_s* simu){
	const parms_s* parms=simu->parms;
	dfree(simu->mideal);
	dfree(simu->mideal_oa);
	simu->mideal=dread("%s_%d.bin", parms->maos.fnmideal, simu->seed_maos);
	if(parms->maos.nmod!=simu->mideal->nx){
		if(parms->maos.nmod>5&&simu->mideal->nx==5){
			warning("Resize mideal (compatible mode)\n");
			dresize(simu->mideal, parms->maos.nmod, simu->mideal->ny);
		} else{
			error("Loaded mideal (%ld rows) does not agree with nmod=%d\n", simu->mideal->nx, parms->maos.nmod);
		}
	}
	dcell* midealp=dcellread("%s_%d.bin", parms->maos.fnmidealp, simu->seed_maos);
	simu->mideal_oa=dref(P(midealp,parms->maos.evlindoa));
	dcellfree(midealp);
	if(0){
		real scale=1;
		dscale(simu->mideal, scale);
		dscale(simu->mideal_oa, scale);
	}
	if(0){
		for(long i=0; i<simu->mideal->ny; i++){
			P(simu->mideal, 2, i)*=10;
			P(simu->mideal_oa, 2, i)*=10;
		}
		warning_once("PS1 mode magnified by 10 times\n");
	}
}
/**
   Update ideal NGS modes with focus or wind shake
*/
static void skysim_update_mideal(sim_s* simu){
	const parms_s* parms=simu->parms;
	if(parms->skyc.addws&&parms->skyc.psd_ws){
		/*Add ws to mideal. After genstars so we don't purturb it. */
		int im_ws=parms->skyc.addws==1?0:1;
		info("Add wind shake time series to mideal mode %d\n", im_ws);
		dmat* telws=psd2ts(parms->skyc.psd_ws, &simu->rand, parms->maos.dt, simu->mideal->ny);
		/*telws is in m. need to convert to rad since mideal is in this unit. */
		real alpha=1./sqrt(P(parms->maos.mcc, 0, 0));
		dscale(telws, alpha);//convert from wfe to radian.
		dmat* pm1=simu->mideal;
		dmat* pm2=simu->mideal_oa;
		
		for(long i=0; i<simu->mideal->ny; i++){
			P(pm1, im_ws, i)+=P(telws,i);
			P(pm2, im_ws, i)+=P(telws,i);
		}
		dfree(telws);
	}
}

/**
   Generate turbulence PSDs from time seris of mideal or from precomputed
   ones. Combine them with windshake.  */
static void skysim_calc_psd(sim_s* simu){
	const parms_s* parms=simu->parms;
	real var_all=0;
	if(parms->skyc.psdcalc){
		dmat* MCC=parms->maos.mcc;
		dmat* x=dtrans(simu->mideal);
		int npsd, iratio;
		if(parms->skyc.gsplit){
			npsd=parms->maos.nmod;
			iratio=1;
		} else{
			npsd=1;
			iratio=0;
		}
		simu->psds=dcellnew(npsd, 1);
		for(int im=0; im<x->ny; im++){
			dmat* xi=dsub(x, 20, 0, im, 1);
			dmat* psdi=psd1dt(xi, 1, parms->maos.dt);
			add_psd2(&P(simu->psds,im*iratio), psdi, P(MCC, im, im));//convert from rad^2 to m^2
			var_all+=psd_inte2(psdi);
			dfree(xi);
			dfree(psdi);
		}
		dfree(x);
	} else{
		error("Please implement interpolating modes\n");

	}

	info("PSD integrates to %.2f nm. varol=%.2f nm\n", sqrt(var_all)*1e9, sqrt(simu->varol)*1e9);
	if(parms->skyc.psd_ws){
		real var_ws=psd_inte2(parms->skyc.psd_ws);
		info("Windshake PSD integrates to %g nm\n", sqrt(var_ws)*1e9);
		simu->varol+=var_ws;//testing

		//add windshake PSD to ngs/tt. 
		int im_ws=parms->skyc.addws==1?0:1;
		add_psd2(&P(simu->psds, im_ws), parms->skyc.psd_ws, 1); //all in m^2 unit.
	}
	if(parms->skyc.dbg||1){
		writebin(simu->psds, "psds_m2.bin");
	}
}

static void skysim_prep_gain(sim_s* simu){
	const parms_s* parms=simu->parms;
	info("Precompute gains for different levels of noise.\n");
	/*dmat *sigma2=dlinspace(0.5e-16,1e-16, 400);// in m2. */
	dmat* sigma2=dlogspace(-20, -10, 400);/*in m2, logspace. */
	simu->gain_pre=(dcccell*)cellnew(parms->skyc.ndtrat, 1);
	int servotype=parms->skyc.servo;
	TIC;tic;
	for(int idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
		long dtrat=P(parms->skyc.dtrats,idtrat);
		P(simu->gain_pre,idtrat)=(dccell*)cellnew(simu->psds->nx, 1);
		for(int ip=0; ip<simu->psds->nx; ip++){
			P(P(simu->gain_pre,idtrat),ip)=servo_optim(parms->maos.dt,
				dtrat, 0, parms->skyc.pmargin, 0, 0, servotype, P(simu->psds, ip), sigma2);
		}
	}
	toc2("servo_optim");
	simu->gain_x=dref(sigma2);
	dfree(sigma2);
	writebin(simu->gain_pre, "gain_pre.bin");
	writebin(simu->gain_x, "gain_x.bin");
}

static void skysim_prep_sde(sim_s* simu){
	const parms_s* parms=simu->parms;
	if(!parms->skyc.psdcalc){
		error("Only support skyc.psdcalc=1\n");
	}
	/*
	  Do not scale the time series by MCC. We are working in "radian space" for
	  kalman. kalman.P is multiplied to mcc to compute Estimation error in nm (not WFE)
	*/
	sde_fit_auto(&simu->sdecoeff, simu->psds, parms->skyc.sdetmax);
	for(int im=0; im<parms->maos.nmod; im++){
		P(simu->sdecoeff,2,im)/=sqrt(P(parms->maos.mcc, im, im));//convert from m to rad
	}
	dshow(simu->sdecoeff, "sde_coeff");
	if(parms->skyc.dbg||1){
		writebin(simu->sdecoeff, "sde_coeff");
	}
}
/**
   Setup the stars fields and then calls skysim_isky() to handle each star field.
*/
void skysim(const parms_s* parms){
	/*if(parms->skyc.dbg){
	writebin(parms->skyc.resfocus, "%s/resfocus",dirsetup);
	}*/
	const int npowfs=parms->maos.npowfs;
	sim_s* simu=mycalloc(1, sim_s);
	simu->status=mycalloc(1, status_t);
	simu->status->info=S_RUNNING;
	simu->status->scale=1;
	simu->status->nseed=parms->maos.nseed;
	simu->status->nthread=parms->skyc.nthread;
	simu->powfs=mycalloc(npowfs, powfs_s);
	simu->parms=parms;
	setup_powfs(simu->powfs, parms);
	genpistat(parms, simu->powfs);
	simu->neaspec_dtrats=dread("%s/neaspec/neaspec_dtrats", dirstart);
	PINIT(simu->mutex_status);
	simu->tk_0=myclockd();
	extern int KALMAN_IN_SKYC;
	KALMAN_IN_SKYC=1; //handle deficient measurement case when there is a single WFS to measure 6 modes.
	for(int iseed_maos=0; iseed_maos<parms->maos.nseed; iseed_maos++){
		simu->iseed=iseed_maos;
		int seed_maos=simu->seed_maos=parms->maos.seeds[iseed_maos];/*loop over seed */
		seed_rand(&simu->rand, parms->skyc.seed+parms->maos.zadeg);
		simu->status->iseed=iseed_maos;
		prep_bspstrehl(simu);
		if(parms->skyc.neanonlin&&parms->skyc.phytype==1){
			simu->nonlin=wfs_nonlinearity(parms, simu->powfs, seed_maos);
		}
		/*generate star fields. moved to here to avoid seed being perturbed*/
		if(parms->skyc.stars){
			info("Loading stars from %s\n", parms->skyc.stars);
			if(zfexist("%s", parms->skyc.stars)){
				simu->stars=dcellread("%s", parms->skyc.stars);
			} else{
				simu->stars=dcellnew(1, 1);
				P(simu->stars,0)=readstr_dmat(0,0,"skyc.stars",parms->skyc.stars);
			}
			if(PN(simu->stars)<1 || NX(simu->stars,0)<2+parms->maos.nwvl){
				error("Loaded stars (%s) has wrong rows (%ld, expect %d)\n", 
					parms->skyc.stars, NX(simu->stars, 0), 2+parms->maos.nwvl);
			}
		} else{
			simu->stars=genstars(parms->skyc.nsky,
				parms->skyc.lat, parms->skyc.lon, parms->skyc.catscl,
				parms->skyc.patfov, parms->maos.nwvl,
				parms->maos.wvl, parms->skyc.maglimit, &simu->rand);
		}
		if(simu->stars->ny!=1){
			simu->stars->nx*=simu->stars->ny;
			simu->stars->ny=1;
		}
		if(simu->stars->nx>parms->skyc.nsky){
			dcellresize(simu->stars, parms->skyc.nsky, 1);
		}
		sortstars(simu->stars);//sort the stars with J from brightest to dimmest.
		writebin(simu->stars, "Res%d_%d_stars", simu->seed_maos, parms->skyc.seed);

		skysim_read_mideal(simu);

		simu->varol=calc_rms(simu->mideal, parms->maos.mcc, parms->skyc.evlstart);
		info("Open loop error: NGS: %.2f\n", sqrt(simu->varol)*1e9);
		
		skysim_calc_psd(simu);
		if(parms->skyc.interpg){
			skysim_prep_gain(simu);
		}
		if(parms->skyc.servo<0){//LQG
			skysim_prep_sde(simu);
		}
		if(parms->skyc.addws&&parms->skyc.psd_ws){
			//Note: it is better to add ws/vib PSD directly to the PSD rather
			// than determine the full PSD from modified time series because 
			//of noise in the PSD and the sde_fit will not be able to work well.
			skysim_update_mideal(simu);//add wind shake / vib time series
		}
		writebin(simu->mideal, "Res%d_%d_minput", seed_maos, parms->skyc.seed);
		int nsky=MIN(simu->stars->nx, parms->skyc.nsky);
		simu->res=dnew_mmap(10, nsky, NULL, "Res%d_%d%s", seed_maos, parms->skyc.seed, parms->skyc.dbg?"_dbg":"");//Total, ATM NGS, ATM TT, WS, 0
		if(!parms->skyc.estimate){
			simu->res_aster=dnew(4+(4+parms->maos.nwvl)*parms->skyc.nwfstot, nsky*(parms->skyc.maxaster?parms->skyc.maxaster:20));//result for all asterisms that have PO results
		}
		simu->res_iaster=0;
		int ng=parms->skyc.ngain;
		if(parms->skyc.servo>0&&!parms->skyc.multirate){
			simu->gain=dcellnewsame_mmap(nsky, 1, ng, parms->maos.nmod,
				NULL, "Res%d_%d_gain", seed_maos, parms->skyc.seed);
		}
		simu->mres=dcellnewsame_mmap(nsky, 1, parms->maos.nmod, parms->maos.nstep,
			NULL, "Res%d_%d_mres", seed_maos, parms->skyc.seed);
		simu->isky_start=parms->skyc.start;
		simu->isky_end=nsky;
		if(parms->skyc.dbgsky>-1){
			simu->isky_start=parms->skyc.dbgsky;
			simu->isky_end=parms->skyc.dbgsky+1;
		}
		simu->status->simstart=simu->isky_start;
		simu->status->simend=simu->isky_end;

		if(simu->isky_start<simu->isky_end){
			simu->isky=simu->isky_start;
			simu->isky_print=simu->isky_start;
			CALL(skysim_isky, simu, NCPU, 0);/*isky iteration. */
		}
		if(parms->skyc.dbgsky<0){
			touch("Res%d_%d.done", seed_maos, parms->skyc.seed);
		}
		close(parms->fdlock[simu->status->iseed]);
		remove(parms->fnlock[simu->status->iseed]); 
		FREE(parms->fnlock[simu->status->iseed]); 
		dcellfree(simu->stars);
		dfree(simu->res);
		if(simu->res_iaster>0){
			dresize(simu->res_aster, NX(simu->res_aster), simu->res_iaster);
			writebin(simu->res_aster, "Res%d_%d_aster", seed_maos, parms->skyc.seed);
		}
		dfree(simu->res_aster);
		dcellfree(simu->gain);
		dcellfree(simu->mres);
		dfree(simu->mideal);
		dfree(simu->mideal_oa);
		dcellfree(simu->psds);
		dfree(simu->sdecoeff);
		/*Free the data used to do bicubic spline. */
		for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
			dcellfreearr(simu->bspstrehl[ipowfs], parms->maos.nsa[ipowfs]*parms->maos.nwvl);
		}
		free(simu->bspstrehl);
		dfree(simu->bspstrehlxy);
		cellfree(simu->gain_pre);
		dfree(simu->gain_x);
		cellfree(simu->nonlin);
	}/*iseed_maos */
	free(simu->status);
	dfree(simu->neaspec_dtrats);
	free_powfs(simu->powfs, parms);
	free(simu->powfs);
	free(simu);
}
