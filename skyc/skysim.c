/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

/**
   The actual work horse that does the physical optics time domain simulation.
*/
static void skysim_isky(SIM_S *simu){
    int isky;
    POWFS_S *powfs=simu->powfs;
    const PARMS_S *parms=simu->parms;
    const dcell *stars=simu->stars;
    const int noisy=parms->skyc.noisy;
    const int seed_maos=simu->seed_maos;
    int nstar;
    int nstep;

    dmat*  pres=simu->res;
    dmat*  pres_oa=simu->res_oa;
    dmat*  pres_geom=simu->res_geom;
    while(LOCKADD(isky, simu->isky, 1)<simu->isky_end){
	double tk_1=myclockd();
	/*Setup star parameters. */
	STAR_S *star=setup_star(&nstar, simu,stars->p[isky],seed_maos);
	if(!star || nstar==0){
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
	nstep=setup_star_read_ztilt(star,nstar,parms,seed_maos);
	for(int iaster=0; iaster<naster; iaster++)
#if _OPENMP >= 200805
#pragma omp task default(shared) firstprivate(iaster)
#endif
	{
	    /*Parallelizing over aster gives same random stream. */
	    seed_rand(&aster[iaster].rand, parms->skyc.seed+iaster+40);
	    /*Compute signal level. */
	    setup_aster_copystar(&aster[iaster], star, parms);
	    /*setup gradient operator. */
	    setup_aster_g(&aster[iaster], star, parms);
	    /*Compute the reconstructor, nea, sigman and optimize controller. */
	    setup_aster_ztilt(&aster[iaster], star, parms);
	    setup_aster_controller(simu, &aster[iaster], star, parms);
	}
#if _OPENMP >= 200805
#pragma omp taskwait
#endif
	/*Select asters that have good performance. */
	setup_aster_select(PCOL(pres_geom,isky),aster, naster, star, 
			   parms->skyc.mtch?0.5*simu->rmsol:INFINITY,parms); 
	double tk_2=myclockd();
	/*Read in physical optics data (wvf) */
	nstep=setup_star_read_wvf(star,nstar,parms,seed_maos);
	double tk_3=myclockd();
	/*
	  Now begin time domain Physical Optics Simulations.
	*/
	double skymini=simu->rmsol;
	int selaster=0;
	int seldtrat=0;
	
	for(int iaster=0; iaster<naster; iaster++)
#if _OPENMP >= 200805
#pragma omp task default(shared) firstprivate(iaster)
#endif
	{
	    ASTER_S *asteri=&aster[iaster];
	    double mini;
	    dmat *pmini=NULL, *min_imres=NULL;
	    int mdtrat=0;
	    asteri->nstep=nstep;
	    if(parms->skyc.dbgaster>-1){
		if(iaster!=parms->skyc.dbgaster){
		    goto skip1;
		}
	    }else{
		if(!asteri->use){
		    goto skip1;
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
		info2("Aster %d, Estimated minimum error is %.2fnm at %.1f Hz. Try %.1f to %.1f Hz\n", iaster,
		      sqrt(asteri->mresol)*1e9, parms->skyc.fss[asteri->mdtrat], 
		      parms->skyc.fss[asteri->idtratmin], parms->skyc.fss[asteri->idtratmax-1]);
	    }
	    /*Assign wvf from star to aster */
	    setup_aster_wvf(asteri, star, parms);
	    /*Compute the reconstructor, nea, sigman and optimize controller again. no need redo?*/
	    //setup_aster_controller(simu, asteri, parms);

	    mini=simu->rmsol;

	    for(int idtrat=asteri->idtratmin; idtrat<asteri->idtratmax; idtrat++)
#if _OPENMP >= 200805
#pragma omp task default(shared) firstprivate(idtrat)
#endif
	    {
		/*focus and windshake residual; */
		double resadd=0;
		if(!parms->skyc.addws){
		    resadd+=asteri->res_ws->p[idtrat];
		}
		dmat *ires=NULL;
		dmat *imres=NULL;
		if((ires=skysim_sim(&imres, simu->mideal, simu->mideal_oa, simu->rmsol,
				    asteri, powfs, parms, idtrat, noisy, parms->skyc.phystart))){
		    if(parms->skyc.verbose){
			info2("%5.1f Hz %7.2f +%7.2f =%7.2f\n", parms->skyc.fss[idtrat], 
			      sqrt(ires->p[0])*1e9, sqrt(resadd)*1e9,
			      sqrt(ires->p[0]+resadd)*1e9);
		    }
		    /*Add windshake contribution. */
		    double tot_1=ires->p[0] + resadd;	
#if _OPENMP >= 200805 
#pragma omp critical 
#endif
		    if(tot_1 < mini){
			mini=tot_1;
			mdtrat=idtrat;
			dfree(pmini); pmini=dref(ires);
			dfree(min_imres); min_imres=dref(imres);
			if(parms->skyc.verbose){
			    info2("Selected: %5.1f Hz %7.2f +%7.2f =%7.2f\n", parms->skyc.fss[idtrat], 
				  sqrt(ires->p[0])*1e9, sqrt(resadd)*1e9,
				  sqrt(ires->p[0]+resadd)*1e9);
			}
		    }
		    dfree(ires);
		}
		dfree(imres);
	    }
#if _OPENMP >= 200805
#pragma omp taskwait
#endif
#if _OPENMP >= 200805 
#pragma omp critical 
#endif	
	    if(mini<skymini){
		selaster=iaster;
		skymini=mini;
		/*Field Averaged Performance. */
		IND(pres,1,isky)=pmini->p[0];/*ATM NGS Mode error */
		IND(pres,2,isky)=pmini->p[1];/*ATM Tip/tilt Error. */
		IND(pres,3,isky)=parms->skyc.addws?0:asteri->res_ws->p[mdtrat];/*Residual wind shake TT*/
		IND(pres,4,isky)=0;/*always zero*/
		IND(pres,0,isky)=IND(pres,1,isky)+IND(pres,3,isky)+IND(pres,4,isky);/*Total */
		/*On axis performance. */
		IND(pres_oa,1,isky)=pmini->p[2];
		IND(pres_oa,2,isky)=pmini->p[3];
		IND(pres_oa,3,isky)=IND(pres,3,isky);
		IND(pres_oa,4,isky)=IND(pres,4,isky);
		IND(pres_oa,0,isky)=IND(pres_oa,1,isky)+IND(pres_oa,3,isky)+IND(pres_oa,4,isky);
		seldtrat = mdtrat;
		simu->fss->p[isky]=parms->skyc.fss[mdtrat];
		if(parms->skyc.verbose){
		    info2("%5.1f Hz: Update Tot: %6.2f nm NGS: %6.2f nm TT: %6.2f nm\n", 
			  simu->fss->p[isky],
			  sqrt(IND(pres,0,isky))*1e9, sqrt(IND(pres,1,isky))*1e9, sqrt(IND(pres,2,isky))*1e9);
		}
		dcp(&simu->mres->p[isky], min_imres);
	    }
	    dfree(pmini);
	    dfree(min_imres);
	  skip1:;
	}/*iaster */
#if _OPENMP >= 200805
#pragma omp taskwait
#endif
	simu->sel->p[isky]->ny=aster[selaster].nwfs;
	dmat* psel=simu->sel->p[isky];
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	    IND(psel,0,iwfs)=aster[selaster].wfs[iwfs].thetax;
	    IND(psel,1,iwfs)=aster[selaster].wfs[iwfs].thetay;
	    for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
		IND(psel,iwvl+2,iwfs)=aster[selaster].wfs[iwfs].mags->p[iwvl];
	    }
	}
	if(parms->skyc.servo>0){
	    dcp(&simu->gain->p[isky], aster[selaster].gain->p[seldtrat]);
	}
	if(parms->skyc.save){
	    skysim_save(simu, aster, PCOL(pres,isky), selaster, seldtrat, isky);
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
	simu->status->clerrlo=sqrt(IND(pres,1,isky))*1e9;
	simu->status->clerrhi=sqrt(IND(pres,0,isky))*1e9;
	scheduler_report(simu->status);
	UNLOCK(simu->mutex_status);
	long totm=(long)floor(simu->status->tot/60.);
	long tots=(long)simu->status->tot-totm*60;
	long laps_h=simu->status->laps/3600;
	long laps_m=simu->status->laps/60-laps_h*60;
	long rest_h=simu->status->rest/3600;
	long rest_m=simu->status->rest/60-rest_h*60;
	info2("Field%4d,%2d stars, aster%3d/%-3d,%3.0f Hz: %6.2f nm "
	      "Sel%3.0fs Load%3.0fs Phy%3.0fs Tot %ld:%02ld Used %ld:%02ld Left %ld:%02ld\n",
	      isky, nstar, selaster, naster, simu->fss->p[isky], sqrt(IND(pres,0,isky))*1e9,
	      tk_2-tk_1, tk_3-tk_2, tk_4-tk_3, totm, tots, laps_h, laps_m, rest_h, rest_m);
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
    if(parms->skyc.addws){
	/*Add ws to mideal. After genstars so we don't purturb it. */
	warning("Add tel wind shake time series to mideal\n");
	dmat *telws=psd2time(parms->skyc.psd_ws, &simu->rand, parms->maos.dt, simu->mideal->ny);
	/*telws is in m. need to convert to rad since mideal is in this unit. */
	dscale(telws, 4./parms->maos.D);//convert from wfe to radian.
	dmat*  pm1=simu->mideal; 
	dmat*  pm2=simu->mideal_oa;
	for(long i=0; i<simu->mideal->ny; i++){
	    IND(pm1,0,i)+=telws->p[i];
	    IND(pm2,0,i)+=telws->p[i];
	}
	dfree(telws);
    }
}

/**
   Generate turbulence PSDs from time seris of mideal or from precomputed
   ones. Combine them with windshake.  */
static void skysim_calc_psd(SIM_S *simu){
    const PARMS_S *parms=simu->parms;
    if(parms->skyc.psdcalc){
	dmat*  MCC=parms->maos.mcc;
	dmat *x=dtrans(simu->mideal);
	for(int im=0; im<x->ny; im++){
	    dmat *xi=dsub(x, 20, 0, im, 1);
	    dscale(xi, sqrt(IND(MCC,im,im)));/*convert to unit of m.*/
	    dmat *psdi=psd1dt(xi, 1, parms->maos.dt);
	    if(im<2){
		add_psd2(&simu->psd_tt, psdi);
	    }else if(im<5){
		add_psd2(&simu->psd_ps, psdi);
	    }else if(im==5){
		add_psd2(&simu->psd_focus, psdi);
	    }
	    dfree(xi);
	    dfree(psdi);
	}
	dfree(x);
	simu->psd_ngs=add_psd(simu->psd_ps, simu->psd_tt);
    }else{
	simu->psd_ngs=ddup(parms->skyc.psd_ngs);
	simu->psd_ps=ddup(parms->skyc.psd_ps);
	simu->psd_tt=ddup(parms->skyc.psd_tt);
	simu->psd_focus=ddup(parms->skyc.psd_focus);
	/*renormalize PSD*/
	double rms_ngs=psd_inte2(simu->psd_ngs);
	if(parms->skyc.psd_scale && !parms->skyc.psdcalc){
	    info2("NGS PSD integrates to %.2f nm before scaling\n", sqrt(rms_ngs)*1e9);
	    double rms_ratio=simu->rmsol/rms_ngs;
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
    //add focus PSD to ngs
    add_psd2(&simu->psd_ngs, simu->psd_focus);

    double rms_ngs=psd_inte2(simu->psd_ngs);
    info2("PSD integrates to %.2f nm.\n", sqrt(rms_ngs)*1e9);
    
    double rms_ws=psd_inte2(parms->skyc.psd_ws);
    info2("Windshake PSD integrates to %g nm\n", sqrt(rms_ws)*1e9);
    simu->rmsol+=rms_ws;//testing
 
    //add windshake PSD to ngs/tt
    add_psd2(&simu->psd_ngs, parms->skyc.psd_ws);
    add_psd2(&simu->psd_tt, parms->skyc.psd_ws);
    if(parms->skyc.dbg){
	writebin(simu->psd_tt, "psd_tt");
	writebin(simu->psd_ps, "psd_ps");
	writebin(simu->psd_ngs, "psd_ngs");
	writebin(simu->psd_focus, "psd_focus");
    }
}

static void skysim_prep_gain(SIM_S *simu){
    const PARMS_S *parms=simu->parms;
    info2("Precompute gains for different levels of noise.\n");
    /*dmat *sigma2=dlinspace(0.5e-16,1e-16, 400);// in m2. */
    dmat *sigma2=dlogspace(-18,-10,400);/*in m2, logspace. */
    simu->gain_tt =mycalloc(parms->skyc.ndtrat,dcell*);
    simu->gain_ps =mycalloc(parms->skyc.ndtrat,dcell*);
    simu->gain_ngs=mycalloc(parms->skyc.ndtrat,dcell*);
    int servotype=parms->skyc.servo;
    if(parms->maos.nmod>5){
	simu->gain_focus=mycalloc(parms->skyc.ndtrat,dcell*);
    }
    TIC;tic;
    for(int idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
	long dtrat=parms->skyc.dtrats->p[idtrat];
	simu->gain_tt[idtrat]=servo_optim(simu->psd_tt, parms->maos.dt,
					  dtrat, parms->skyc.pmargin, sigma2, servotype);
	simu->gain_ps[idtrat]=servo_optim(simu->psd_ps, parms->maos.dt, 
					  dtrat, parms->skyc.pmargin, sigma2, servotype);
	simu->gain_ngs[idtrat]=servo_optim(simu->psd_ngs, parms->maos.dt,
					   dtrat, parms->skyc.pmargin, sigma2, servotype);
	if(parms->maos.nmod>5){
	    simu->gain_focus[idtrat]=servo_optim(simu->psd_focus, parms->maos.dt,
						 dtrat, parms->skyc.pmargin, sigma2, servotype);
	}
	if(parms->skyc.dbg){
	    writebin(simu->gain_tt[idtrat],  "gain_tt_%ld.bin", dtrat);
	    writebin(simu->gain_ps[idtrat],  "gain_ps_%ld.bin", dtrat);
	    writebin(simu->gain_ngs[idtrat], "gain_ngs_%ld.bin", dtrat);
	}
    }
    toc2("servo_optim");
    simu->gain_x=dref(sigma2);
    dfree(sigma2);
}

static void skysim_prep_sde(SIM_S *simu){
    const PARMS_S *parms=simu->parms;
    if(!parms->skyc.psdcalc){
	error("Only support skyc.psdcalc=1\n");
    }
    /*
      Do not scale the time series by MCC. We are working in "radian space" for
      kalman. kalman.P is multiplied to mcc to compute Estimation error in nm (not WFE)
    */
    dmat *x=dtrans(simu->mideal);
    simu->psdi=dcellnew(x->ny, 1);
    simu->sdecoeff=dnew(3,x->ny);
    dmat*  pcoeff=simu->sdecoeff;
    for(int im=0; im<x->ny; im++){
	dmat *xi=dsub(x, 20, 0, im, 1);
	simu->psdi->p[im]=psd1dt(xi, 1, parms->maos.dt);
	dfree(xi);
	if(im==0 && parms->skyc.psd_ws){
	    //add windshake on first mode only
	    add_psd2(&simu->psdi->p[im], parms->skyc.psd_ws);
	}
	dmat *coeff=sde_fit(simu->psdi->p[im], NULL, parms->skyc.sdetmax, 0);
	if(IND(coeff,0,0)>100 && parms->skyc.sdetmax){
	    dfree(coeff);
	    coeff=sde_fit(simu->psdi->p[im], NULL, 0, 0);
	    if(IND(coeff,0,0)>100){
		warning("sde_fit returns unsuitable values\n");
	    }
	}
	if(coeff->ny>1){
	    error("Please handle this case\n");
	}
	memcpy(PCOL(pcoeff, im), coeff->p, coeff->nx*sizeof(double));
	dfree(coeff);
    }
    
    dfree(x);
    if(parms->skyc.dbg || 1){
	writebin(simu->sdecoeff, "coeff");
	writebin(simu->psdi, "psdi");
    }
}
/**
   Setup the stars fields and then calls skysim_isky() to handle each star field.
*/
void skysim(const PARMS_S *parms){
    /*if(parms->skyc.dbg){
	writebin(parms->skyc.resfocus, "%s/resfocus",dirsetup);
	}*/
    const int npowfs=parms->maos.npowfs;
    SIM_S *simu=mycalloc(1,SIM_S);
    simu->status=mycalloc(1,STATUS_T);
    simu->status->info=S_RUNNING;
    simu->status->scale=1;
    simu->status->nseed=parms->maos.nseed;
    simu->status->nthread=parms->skyc.nthread;
    simu->status->timstart=myclocki();
    simu->powfs=mycalloc(npowfs,POWFS_S);
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
	if(parms->skyc.neanonlin){
	    simu->nonlin=wfs_nonlinearity(parms, simu->powfs, seed_maos);
	}
	skysim_read_mideal(simu);
	simu->rmsol=calc_rms(simu->mideal, parms->maos.mcc, parms->skyc.evlstart);
    	info2("Open loop error: NGS: %.2f\n", sqrt(simu->rmsol)*1e9);
	if(parms->skyc.servo<0){//LQG
	    skysim_prep_sde(simu);
	}else{
	    skysim_calc_psd(simu);
	    if(parms->skyc.interpg){
		skysim_prep_gain(simu);
	    }
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
	writebin(simu->stars, "Res%d_%d_stars",simu->seed_maos,parms->skyc.seed);
	if(parms->skyc.addws){
	    skysim_update_mideal(simu);
	}
	int nsky=MIN(simu->stars->nx, parms->skyc.nsky);
	simu->res   =dnew_mmap(5,nsky,NULL, "Res%d_%d", seed_maos, parms->skyc.seed);
	simu->res_oa=dnew_mmap(5,nsky,NULL, "Res%d_%d_oa", seed_maos, parms->skyc.seed);
	simu->res_geom=dnew_mmap(3,nsky,NULL, "Res%d_%d_geom", seed_maos, parms->skyc.seed);
	simu->fss   =dnew_mmap(nsky,1,NULL, "Res%d_%d_fss", seed_maos, parms->skyc.seed);
	int ng=parms->skyc.ngain;
	simu->gain  =dcellnewsame_mmap(nsky, 1, ng, parms->maos.nmod,
				       NULL, "Res%d_%d_gain", seed_maos, parms->skyc.seed);
	simu->sel   =dcellnewsame_mmap(nsky, 1, 2+parms->maos.nwvl, parms->skyc.nwfstot,
				       NULL,"Res%d_%d_sel", seed_maos, parms->skyc.seed);
	simu->mres  =dcellnewsame_mmap(nsky, 1, parms->maos.nmod, parms->maos.nstep,
				       NULL,"Res%d_%d_mres",  seed_maos, parms->skyc.seed);
	dset(simu->res, simu->rmsol);
	dset(simu->res_oa, simu->rmsol);
	simu->isky_start=parms->skyc.start;
	simu->isky_end=nsky;
	if(parms->skyc.dbgsky>-1){
	    simu->isky_start=parms->skyc.dbgsky;
	    simu->isky_end=parms->skyc.dbgsky+1;
	}
	simu->status->simstart=simu->isky_start;
	simu->status->simend=simu->isky_end;
#if _OPENMP >= 200805//parallel within each sky
	int nthread=2;
#else
	int nthread=parms->skyc.nthread;
#endif
	if(simu->isky_start < simu->isky_end){
	    simu->isky=simu->isky_start;
	    CALL((thread_fun)skysim_isky, simu, nthread,0);/*isky iteration. */
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
	dfree(simu->fss);
	dcellfree(simu->gain);
	dcellfree(simu->sel);
	dcellfree(simu->mres);
	dfree(simu->mideal);
	dfree(simu->mideal_oa);
	dfree(simu->psd_ngs);
	dfree(simu->psd_ps);
	dfree(simu->psd_tt);
	dfree(simu->psd_focus);
	dcellfree(simu->psdi);
	dfree(simu->sdecoeff);
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
	dfree(simu->gain_x);
	dcellfreearr(simu->nonlin, parms->maos.npowfs);free(simu->nonlin);
    }/*iseed_maos */
    free(simu->status);
    dfree(simu->neaspec_dtrats);
    free_powfs(simu->powfs,parms);
    free(simu->powfs);
    free(simu);
}
