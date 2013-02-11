/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
*/

#include "skyc.h"
#include "parms.h"
#include "skysim.h"
#include "types.h"
#include "setup_powfs.h"
#include "setup_aster.h"
#include "skysim_utils.h"
#include "mtch.h"
#include "genstars.h"
#include "genpistat.h"
#include "setup_star.h"
#include "utils.h"
#include "nafocus.h"
#define VERBOSE 2

/**
   The actual work horse that does the physical optics time domain simulation.
*/
static void skysim_isky(SIM_S *simu){
    int isky;
    POWFS_S *powfs=simu->powfs;
    const PARMS_S *parms=simu->parms;
    const dcell *stars=simu->stars;
    const int noisy=parms->skyc.noisy;
    const int do_demote_end=parms->skyc.demote?2:1;
    const int seed_maos=simu->seed_maos;
    int nstar;
    int nstep;

    PDMAT(simu->res, pres);
    PDMAT(simu->res_oa, pres_oa);
    PDMAT(simu->res_geom, pres_geom);
    while(LOCKADD(isky, simu->isky, 1)<simu->isky_end){
	double tk_1=myclockd();
	/*Setup star parameters. */
	STAR_S *star=setup_star(&nstar, simu,stars->p[isky],seed_maos);
	if(!star){
	    info2("Field %d, No stars available\n", isky);
	    continue;
	}
	int naster;
	ASTER_S *aster=setup_aster_comb(&naster, nstar, parms);
	if(!aster || aster==0){
	    info2("Field %d, Aster is empty. skip\n", isky);
	    continue;
	}
	/*
	  We first estimate the matched filter, reconstructor, and servo
	  loop optimization to determine the approximate wavefront error. Only
	  a few combinations are kept for each star field for further time
	  domain simulations.
	*/

	for(int iaster=0; iaster<naster; iaster++){
	    /*Parallelizing over aster gives same random stream. */
	    seed_rand(&aster[iaster].rand, parms->skyc.seed+iaster+40);
	    /*Compute signal level. */
	    setup_aster_copystar(&aster[iaster], star, parms);
	    /*setup gradient operator. */
	    setup_aster_g(&aster[iaster], star, powfs, parms);
	    /*Compute the reconstructor, nea, sigman */
	    setup_aster_recon(&aster[iaster], star,parms);
	    /*Optimize servo gains. */
	    setup_aster_servo(simu, &aster[iaster], parms);
	}
	/*Select asters that have good performance. */
	setup_aster_select(pres_geom[isky],aster, naster, star, 
			   parms->skyc.mtch?0.5*simu->rmsol->p[0]:INFINITY,parms); 
	/*Read in physical optics data (wvf) */
	nstep=setup_star_read_wvf(star,nstar,parms,seed_maos);
	double tk_3=myclockd();
	/*
	  Now begin time domain Physical Optics Simulations.
	*/
	double skymini=simu->rmsol->p[0];
	int selaster=0;
	int seldtrat=0;
	
	for(int iaster=0; iaster<naster; iaster++){
	    ASTER_S *asteri=&aster[iaster];
	    asteri->nstep=nstep;
	    if(parms->skyc.dbgaster>-1){
		if(iaster!=parms->skyc.dbgaster){
		    continue;
		}
	    }else{
		if(!asteri->use){
		    continue;
		}
	    }
	    if(parms->skyc.verbose>1){
		for(int iwfs=0; iwfs<aster[iaster].nwfs; iwfs++){
		    info2("wfs %d: istar=%d, ipowfs=%d\n",iwfs,aster[iaster].wfs[iwfs].istar,
			  aster[iaster].wfs[iwfs].ipowfs);
		    info2("wfs %d: at (%g,%g). siglev=%g\n",iwfs,
			  aster[iaster].wfs[iwfs].thetax*206265,
			  aster[iaster].wfs[iwfs].thetay*206265, aster[iaster].wfs[iwfs].siglevtot);
		}
	    }
	    if(parms->skyc.verbose){
		info2("Aster %d, Estimated minimum error is %.2fnm at %.1f Hz\n", iaster,
		      sqrt(asteri->mresol)*1e9, parms->skyc.fss[asteri->mdtrat]);
	    }
	    /*Copy wvf from star to aster */
	    setup_aster_wvf(asteri, star, parms);
	    /*Compute the reconstructor, nea, sigman */
	    setup_aster_recon(asteri, star, parms);
	    /*Optimize servo gains. */
	    setup_aster_servo(simu, asteri, parms);

	    double mini=simu->rmsol->p[0];
	    dmat *pmini=NULL;
	    dmat *min_imres=NULL;
	    int mdtrat=0;
	    int demote=-1;
	    for(int idtrat=aster->idtratmin; idtrat<=aster->idtratmax; idtrat++){
		/*focus and windshake residual; */
		double resadd=0;
		if(!parms->skyc.addws){
		    resadd+=asteri->res_ws->p[idtrat];
		}
		if(parms->maos.nmod<6){//no focus mode in time domain
		    resadd+=parms->skyc.resfocus->p[idtrat];
		}
		dmat *ires=NULL;
		dmat *imres=NULL;
		for(int do_demote=0; do_demote<do_demote_end; do_demote++){
		    ires=skysim_phy(&imres, simu->mideal, simu->mideal_oa, simu->rmsol->p[0], 
				    asteri, powfs, parms, idtrat, noisy, do_demote);
		    if(ires){
			if(parms->skyc.verbose){
			    if(do_demote){
				info2("  DemoteTTF: ");
			    }
			    info2("%5.1f Hz %7.2f +%7.2f =%7.2f", parms->skyc.fss[idtrat], 
				  sqrt(ires->p[0])*1e9, sqrt(resadd)*1e9,
				  sqrt(ires->p[0]+resadd)*1e9);
			}
			/*Add windshake contribution. */
			double tot_1=ires->p[0] + resadd;
			if(tot_1 < mini){
			    mini=tot_1;
			    mdtrat=idtrat;
			    demote=do_demote;
			    dfree(pmini); pmini=dref(ires);
			    dfree(min_imres); min_imres=dref(imres);
			}
		    }
		}
		if(ires && parms->skyc.verbose){
		    info2("\n");
		}
		dfree(ires);
		dfree(imres);
	    }
	    if(mini<skymini){
		selaster=iaster;
		skymini=mini;
		/*Field Averaged Performance. */
		pres[isky][1]=pmini->p[0];/*ATM NGS Mode error */
		pres[isky][2]=pmini->p[1];/*ATM Tip/tilt Error. */
		pres[isky][3]=parms->skyc.addws?0:asteri->res_ws->p[mdtrat];/*Residual wind shake TT*/
		pres[isky][4]=parms->maos.nmod>5?0:parms->skyc.resfocus->p[mdtrat];/*Residual focus tracking error. */
		pres[isky][0]=pres[isky][1]+pres[isky][3]+pres[isky][4];/*Total */
		/*On axis performance. */
		pres_oa[isky][1]=pmini->p[2];
		pres_oa[isky][2]=pmini->p[3];
		pres_oa[isky][3]=pres[isky][3];
		pres_oa[isky][4]=pres[isky][4];
		pres_oa[isky][0]=pres_oa[isky][1]+pres_oa[isky][3]+pres_oa[isky][4];
		seldtrat = mdtrat;
		simu->fss->p[isky]=parms->skyc.fss[mdtrat];
		simu->demote->p[isky]=demote;
		if(parms->skyc.verbose){
		    info2("%5.1f Hz: Update Tot: %6.2f nm NGS: %6.2f nm TT: %6.2f nm\n", 
			  simu->fss->p[isky],
			  sqrt(pres[isky][0])*1e9, sqrt(pres[isky][1])*1e9, sqrt(pres[isky][2])*1e9);
		}
		dcp(&simu->mres->p[isky], min_imres);
	    }
	    dfree(pmini);
	    dfree(min_imres);
	}/*iaster */

	PDMAT(simu->sel->p[isky],psel);
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	    psel[iwfs][0]=aster[selaster].wfs[iwfs].thetax;
	    psel[iwfs][1]=aster[selaster].wfs[iwfs].thetay;
	    for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
		psel[iwfs][iwvl+2]=aster[selaster].wfs[iwfs].mags->p[iwvl];
	    }
	}
	dcp(&simu->gain->p[isky], aster[selaster].gain->p[seldtrat]);
	if(parms->skyc.save){
	    skysim_save(simu, aster, pres[isky], selaster, seldtrat, isky);
	}
	free_aster(aster, naster, parms);
	free_star(star, nstar, parms);
	double tk_4=myclockd();
	LOCK(simu->mutex_status);
	simu->status->isim=isky;
	simu->status->tot=tk_4-tk_1;/*per step */
	simu->status->laps=tk_4-simu->tk_0;
	int nsky_tot=simu->isky_end-simu->isky_start;
	int nsky_left=simu->isky_end-simu->isky-1+nsky_tot*(parms->maos.nseed-simu->iseed-1);
	int nsky_laps=simu->isky-simu->isky_start+1+nsky_tot*simu->iseed;
	simu->status->rest=simu->status->laps*nsky_left/nsky_laps;
	simu->status->clerrlo=sqrt(pres[isky][1])*1e9;
	simu->status->clerrhi=sqrt(pres[isky][0])*1e9;
	scheduler_report(simu->status);
	UNLOCK(simu->mutex_status);
	long totm=(long)floor(simu->status->tot/60.);
	long tots=(long)simu->status->tot-totm*60;
	long laps_h=simu->status->laps/3600;
	long laps_m=simu->status->laps/60-laps_h*60;
	long rest_h=simu->status->rest/3600;
	long rest_m=simu->status->rest/60-rest_h*60;
	info2("Field %3d,%2d stars,%3d aster,%3.0f Hz: %6.2f %6.2f %6.2f nm "
	      "Phy %3.0fs Tot %ld:%02ld Used %ld:%02ld Left %ld:%02ld\n",
	      isky, nstar, naster, simu->fss->p[isky],
	      sqrt(pres[isky][0])*1e9, sqrt(pres[isky][1])*1e9, sqrt(pres[isky][2])*1e9,
	      tk_4-tk_3, totm, tots, laps_h, laps_m, rest_h, rest_m);
    }/*while */
}

/**
   Read in ideal NGS modes
*/
static void skysim_read_mideal(SIM_S *simu){
    const PARMS_S *parms=simu->parms;
    dfree(simu->mideal);
    dfree(simu->mideal_oa);
    simu->mideal=dread("%s_%d.bin",parms->maos.fnmideal,simu->seed_maos);
    if(parms->maos.nmod>5 && simu->mideal->nx==5){
	warning("Resize mideal (compatible mode)\n");
	dresize(simu->mideal, parms->maos.nmod, simu->mideal->ny);
    }	
    dcell *midealp=dcellread("%s_%d.bin",parms->maos.fnmidealp,simu->seed_maos);
    simu->mideal_oa=dref(midealp->p[parms->maos.evlindoa]);
    dcellfree(midealp);
}
/**
   Update ideal NGS modes with focus or wind shake
*/
static void skysim_update_mideal(SIM_S *simu){
    const PARMS_S *parms=simu->parms;
    if(parms->skyc.addfocus && parms->maos.nmod>5){
	warning("Incorrect and Deprecated; Please use sim.mffocus to simulation na focus variaiton in maos\n");
	dmat *range=NULL;
	if(parms->skyc.fnrange){
	    range=dread("%s", parms->skyc.fnrange);
	    info("Loading sodium range variation from %s\n", parms->skyc.fnrange);
	    if(range->nx<simu->mideal->ny){
		error("Time serials is not long enough. Need %ld, got %ld\n",
		      simu->mideal->ny, range->nx);
	    }
	}else{
	    range=nafocus_time(parms->skyc.na_alpha, parms->skyc.na_beta, 
			       parms->maos.dt, simu->mideal->ny, &simu->rand);
	}
	dwrite(range, "narange_%d", parms->skyc.seed);
	double scale=0.5*pow(1./parms->maos.hs, 2)*(1./(parms->maos.za));
	for(int istep=0; istep<parms->maos.nstep; istep++){
	    simu->mideal->p[5+istep*6]+=range->p[istep]*scale;
	}
	dfree(range);
    }
    if(parms->skyc.addws){
	/*Add ws to mideal. After genstars so we don't purturb it. */
	warning("Add tel wind shake time series to mideal\n");
	dmat *telws=psd2time(simu->psd_ws, &simu->rand, parms->maos.dt, simu->mideal->ny);
	/*telws is in m. need to convert to rad since mideal is in this unit. */
	dwrite(telws, "telws_%d", parms->skyc.seed);
	dscale(telws, 4./parms->maos.D);//convert from wfe to radian.
	PDMAT(simu->mideal, pm1); 
	PDMAT(simu->mideal_oa, pm2);
	for(long i=0; i<simu->mideal->ny; i++){
	    pm1[i][0]+=telws->p[i];
	    pm2[i][0]+=telws->p[i];
	}
	dfree(telws);
    }
}

/**
   Generate turbulence PSDs from time seris of mideal or from precomputed
   ones. Combine them with windshake.  */
static void skysim_calc_psd(SIM_S *simu){
    const PARMS_S *parms=simu->parms;
    simu->rmsol=calc_rmsol(simu->mideal, parms);
    simu->psd_ws=ddup(parms->skyc.psd_ws);
    if(parms->skyc.psdcalc){
	PDMAT(parms->maos.mcc, MCC);
	dmat *x=dtrans(simu->mideal);
	for(int im=0; im<x->ny && im<5; im++){
	    dmat *xi=dsub(x, 0, 0, im, 1);
	    dscale(xi, sqrt(MCC[im][im]));/*convert to unit of m.*/
	    dmat *psdi=psd1dt(xi, xi->nx, parms->maos.dt);
	    if(im<2){
		add_psd2(&simu->psd_tt, psdi);
	    }else if(im<5){
		add_psd2(&simu->psd_ps, psdi);
	    }
	    dfree(xi);
	    dfree(psdi);
	}
	dfree(x);
	simu->psd_ngs=add_psd(simu->psd_ps, simu->psd_tt);
	dwrite(simu->psd_tt, "psd_tt");
	dwrite(simu->psd_ps, "psd_ps");
	dwrite(simu->psd_ngs, "psd_ngs");
    }else{
	simu->psd_ngs=ddup(parms->skyc.psd_ngs);
	simu->psd_ps=ddup(parms->skyc.psd_ps);
	simu->psd_tt=ddup(parms->skyc.psd_tt);
	/*renormalize PSD */
	double rms_ngs=psd_inte2(simu->psd_ngs);
	if(parms->skyc.psd_scale && !parms->skyc.psdcalc){
	    info2("NGS PSD integrates to %.2f nm before scaling\n", sqrt(rms_ngs)*1e9);
	    double rms_ratio=simu->rmsol->p[0]/rms_ngs;
	    info2("Scaling PSD by %g\n", rms_ratio);
	    long nx=simu->psd_ngs->nx;
	    /*scale PSF in place. */
	    double *p_ngs=simu->psd_ngs->p+nx;
	    double *p_tt=simu->psd_tt->p+nx;
	    double *p_ps=simu->psd_ps->p+nx;
	    for(long i=0; i<nx; i++){
		p_ngs[i]*=rms_ratio;
		p_tt[i]*=rms_ratio;
		p_ps[i]*=rms_ratio;
	    }
	}
    }
    double rms_ngs=psd_inte2(simu->psd_ngs);
    info2("PSD integrates to %.2f nm.\n", sqrt(rms_ngs)*1e9);
    
    if(parms->maos.nmod>5){
	PDMAT(parms->maos.mcc, MCC);
	dmat *x=dtrans(simu->mideal);
	dmat *xi=dsub(x, 0, 0, 5, 1);/*focus PSD from mod5.*/
	dscale(xi, sqrt(MCC[5][5]));/*convert to unit of m.*/
	dmat *psdi=psd1dt(xi, xi->nx, parms->maos.dt);
	double rms_focus_atm=psd_inte2(psdi);
	info2("Atmosphere focus PSD integrates to %g nm\n", sqrt(rms_focus_atm)*1e9);
	if(!parms->maos.mffocus){
	    //add sodium focus PSD if we don't do focus tracking in maos
	    double alpha=parms->skyc.na_alpha;
	    /*convert height error to wfe*/
	    double scale=1./(16*sqrt(3))*pow((parms->maos.D/parms->maos.hs),2);
	    double beta2=parms->skyc.na_beta*scale*scale;
	    PDMAT(psdi, pp);
	    for(int i=1; i<psdi->nx; i++){//skip DC.
		pp[1][i]+=pow(pp[0][i], alpha)*beta2;
	    }
	    double rms_focus_na=psd_inte2(psdi)-rms_focus_atm;
	    info2("Sodium focus PSD integrates to %g nm\n", sqrt(rms_focus_na)*1e9);
	}
	add_psd2(&simu->psd_focus, psdi);
	dfree(xi);
	dfree(psdi);
	dfree(x);
    }
    double rms_ws=psd_inte2(simu->psd_ws);
    info2("Windshake PSD integrates to %g nm\n", sqrt(rms_ws)*1e9);
    simu->rmsol->p[0]+=rms_ws;//testing
    simu->rmsol->p[1]+=rms_ws;/*add wind shake to open loop error. */
    //add windshake PSD to ngs/tt
    add_psd2(&simu->psd_ngs, simu->psd_ws);
    add_psd2(&simu->psd_tt, simu->psd_ws);
    if(parms->skyc.gsplit!=1){
	add_psd2(&simu->psd_ngs, simu->psd_focus);
	add_psd2(&simu->psd_ps, simu->psd_focus);
    }
}

static void skysim_prep_gain(SIM_S *simu){
    const PARMS_S *parms=simu->parms;
    info2("Precompute gains for different levels of noise.\n");
    /*dmat *sigma2=dlinspace(0.5e-16,1e-16, 400);// in m2. */
    dmat *sigma2=dlogspace(-18,-10,400);/*in m2, logspace. */
    simu->gain_tt =calloc(parms->skyc.ndtrat, sizeof(dcell*));
    simu->gain_ps =calloc(parms->skyc.ndtrat, sizeof(dcell*));
    simu->gain_ngs=calloc(parms->skyc.ndtrat, sizeof(dcell*));
    if(parms->maos.nmod>5){
	simu->gain_focus=calloc(parms->skyc.ndtrat, sizeof(dcell*));
    }
    TIC;tic;
    for(int idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
	long dtrat=parms->skyc.dtrats[idtrat];
	simu->gain_tt[idtrat]=servo_optim(simu->psd_tt, parms->maos.dt,
					  dtrat, parms->skyc.pmargin, sigma2, 2);
	simu->gain_ps[idtrat]=servo_optim(simu->psd_ps, parms->maos.dt, 
					  dtrat, parms->skyc.pmargin, sigma2, 2);
	simu->gain_ngs[idtrat]=servo_optim(simu->psd_ngs, parms->maos.dt,
					   dtrat, parms->skyc.pmargin, sigma2, 2);
	if(parms->maos.nmod>5){
	    simu->gain_focus[idtrat]=servo_optim(simu->psd_focus, parms->maos.dt,
						 dtrat, parms->skyc.pmargin, sigma2, 2);
	}
	/*dcellwrite(simu->gain_tt[idtrat],  "gain_tt_%ld.bin", dtrat);
	  dcellwrite(simu->gain_ps[idtrat],  "gain_ps_%ld.bin", dtrat);
	  dcellwrite(simu->gain_ngs[idtrat], "gain_ngs_%ld.bin", dtrat);*/
    }
    toc2("servo_optim");
    simu->gain_x=dref(sigma2);
    dfree(sigma2);
}
/**
   Setup the stars fields and then calls skysim_isky() to handle each star field.
*/
void skysim(const PARMS_S *parms){
    if(parms->skyc.dbg){
	dwrite(parms->skyc.resfocus, "%s/resfocus",dirsetup);
    }
    const int npowfs=parms->maos.npowfs;
    SIM_S *simu=calloc(1, sizeof(SIM_S));
    simu->status=calloc(1, sizeof(STATUS_T));
    simu->status->info=S_RUNNING;
    simu->status->scale=1;
    simu->status->nseed=parms->maos.nseed;
    simu->status->nthread=parms->skyc.nthread;
    simu->status->timstart=myclocki();
    simu->powfs=calloc(npowfs, sizeof(POWFS_S));
    simu->parms=parms;
    setup_powfs(simu->powfs, parms);
    genpistat(parms, simu->powfs); 
 
    PINIT(simu->mutex_status);
    simu->tk_0=myclockd();
    for(int iseed_maos=0; iseed_maos<parms->maos.nseed; iseed_maos++){
	simu->iseed=iseed_maos;
	int seed_maos=simu->seed_maos=parms->maos.seeds[iseed_maos];/*loop over seed */
	seed_rand(&simu->rand, parms->skyc.seed+parms->maos.zadeg);
	simu->status->iseed=iseed_maos;
	prep_bspstrehl(simu);
	skysim_read_mideal(simu);
	skysim_calc_psd(simu);
	info2("Open loop error: NGS: %.2f TT: %.2f nm\n", 
	      sqrt(simu->rmsol->p[0])*1e9, sqrt(simu->rmsol->p[1])*1e9);
	
	if(parms->skyc.interpg){
	    skysim_prep_gain(simu);
	}
	/*generate star fields. */
	if(parms->skyc.stars){
	    info2("Loading stars from %s\n",parms->skyc.stars);
	    simu->stars=dcellread("%s",parms->skyc.stars);
	}else{
	    simu->stars=genstars(parms->skyc.nsky, 
				 parms->skyc.lat, parms->skyc.lon,parms->skyc.catscl,
				 parms->skyc.patfov,parms->maos.nwvl, 
				 parms->maos.wvl, &simu->rand);
	}
	if(simu->stars->ny!=1){
	    simu->stars->nx*=simu->stars->ny;
	    simu->stars->ny=1;
	}
	sortstars(simu->stars);//sort the stars with J from brightest to dimmest.
	dcellwrite(simu->stars, "Res%d_%d_stars",simu->seed_maos,parms->skyc.seed);
	if(parms->skyc.addfocus || parms->skyc.addws){
	    skysim_update_mideal(simu);
	}
	int nsky=MIN(simu->stars->nx, parms->skyc.nsky);
	simu->res   =dnew_mmap(5,nsky,NULL, "Res%d_%d", seed_maos, parms->skyc.seed);
	simu->res_oa=dnew_mmap(5,nsky,NULL, "Res%d_%d_oa", seed_maos, parms->skyc.seed);
	simu->res_geom=dnew_mmap(3,nsky,NULL, "Res%d_%d_geom", seed_maos, parms->skyc.seed);
	simu->fss   =dnew_mmap(nsky,1,NULL, "Res%d_%d_fss", seed_maos, parms->skyc.seed);
	simu->demote=dnew_mmap(nsky,1,NULL, "Res%d_%d_demote", seed_maos, parms->skyc.seed);
	simu->gain  =dcellnewsame_mmap(nsky, 1, 3, parms->maos.nmod,
				       NULL, "Res%d_%d_gain", seed_maos, parms->skyc.seed);
	simu->sel   =dcellnewsame_mmap(nsky, 1, 2+parms->maos.nwvl, parms->skyc.nwfstot,
				       NULL,"Res%d_%d_sel", seed_maos, parms->skyc.seed);
	simu->mres  =dcellnewsame_mmap(nsky, 1, parms->maos.nmod, parms->maos.nstep,
				       NULL,"Res%d_%d_mres",  seed_maos, parms->skyc.seed);
	dset(simu->res, simu->rmsol->p[0]);
	dset(simu->res_oa, simu->rmsol->p[0]);
	simu->isky_start=parms->skyc.start;
	simu->isky_end=nsky;
	if(parms->skyc.dbgsky>-1){
	    simu->isky_start=parms->skyc.dbgsky;
	    simu->isky_end=parms->skyc.dbgsky+1;
	}
	simu->status->simstart=simu->isky_start;
	simu->status->simend=simu->isky_end;
	if(simu->isky_start < simu->isky_end){
	    simu->isky=simu->isky_start;
	    CALL(skysim_isky, simu, parms->skyc.nthread,0);/*isky iteration. */
	}
	if(parms->skyc.dbgsky<0){
	    char fn[80];
	    char fnnew[80];
	    snprintf(fn, 80, "Res%d_%d.lock",seed_maos, parms->skyc.seed);
	    snprintf(fnnew, 80, "Res%d_%d.done",seed_maos, parms->skyc.seed);
	    (void)rename(fn, fnnew);
	    close(parms->fdlock[simu->status->iseed]);
	}
	dcellfree(simu->stars);
	dfree(simu->res);
	dfree(simu->res_oa);
	dfree(simu->res_geom);
	dcellfree(simu->mres);
	dcellfree(simu->sel);
	dfree(simu->fss);
	dfree(simu->mideal);
	dfree(simu->mideal_oa);
	dfree(simu->rmsol);
	dfree(simu->psd_ws);
	dfree(simu->psd_ngs);
	dfree(simu->psd_ps);
	dfree(simu->psd_ws);
	/*Free the data used to do bicubic spline. */
	for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	    dcellfreearr(simu->bspstrehl[ipowfs], parms->maos.nsa[ipowfs]*parms->maos.nwvl);
	}
	free(simu->bspstrehl);
	dfree(simu->bspstrehlxy);
	if(parms->skyc.interpg){
	    dcellfreearr(simu->gain_ps, parms->skyc.ndtrat);
	    dcellfreearr(simu->gain_tt, parms->skyc.ndtrat);
	    dcellfreearr(simu->gain_ngs, parms->skyc.ndtrat);
	    dcellfreearr(simu->gain_focus, parms->skyc.ndtrat);
	}
    }/*iseed_maos */
    free(simu->status);
    free_powfs(simu->powfs,parms);
    free(simu->powfs);
    free(simu);
}
