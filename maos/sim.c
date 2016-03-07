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
   \file maos/sim.c
   Call various functions to do the simulation and evaluate performance.
*/
#include <search.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "common.h"
#include "recon.h"
#include "recon_utils.h"
#include "setup_powfs.h"
#include "sim.h"
#include "sim_utils.h"
#include "fdpcg.h"
#include "save.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
#if HAVE_NUMA_H
#include <numa.h>
#endif
extern int PARALLEL;
extern int draw_single;
static double tk_0;
static double tk_1;
static double tk_atm=0;
SIM_T *maos_iseed(int iseed){
    if(iseed==0) tk_0=myclockd();
    tk_1=myclockd();
    const PARMS_T *parms=global->parms;
    POWFS_T *powfs=global->powfs;
    APER_T  *aper =global->aper;
    RECON_T *recon=global->recon;
    if(parms->fdlock && parms->fdlock->p[iseed]<0){
	warning("Another MAOS is already running. Skip seed %ld\n", 
		parms->sim.seeds->p[iseed]);
	return 0;
    }
    if(!parms->sim.pause){
	draw_single=1;//Only draw active frame.
    }else{
	draw_single=0;
    }
    global->iseed=iseed;
    SIM_T *simu=init_simu(parms,powfs,aper,recon,iseed);
    global->simu=simu;
    if(parms->atm.frozenflow){
	genatm(simu);/*Generating atmospheric screen(s) that frozen flows.*/
	if(parms->tomo.predict){
	    if(recon->HXWtomo){
		setup_recon_HXW_predict(simu);
	    }
	    if(parms->tomo.precond==1){
		fdpcg_free(recon->fdpcg);
		recon->fdpcg=fdpcg_prepare(parms, recon, powfs, parms->tomo.predict?simu->atm:NULL);
	    }
	}
    }
#if USE_CUDA
    if(parms->gpu.evl || parms->gpu.wfs){
	/*put here to avoid messing up timing due to transfering. */
	gpu_atm2gpu(simu->atm, simu->atmscale, parms, iseed, parms->sim.start);/*takes 0.4s for NFIRAOS. */
	if(parms->tomo.predict){
	    gpu_update_recon_cn2(parms, recon);
	}
    }
#endif
    return simu;
}

void maos_isim(int isim){
    const PARMS_T *parms=global->parms;
    RECON_T *recon=global->recon;
    SIM_T   *simu =global->simu;
    int iseed=global->iseed;
    int simstart=parms->sim.start;
    int simend=parms->sim.end;
    if(isim==simstart+1){//skip slow first step.
	tk_atm=myclockd();
    }
    if(isim+1+parms->sim.dtrat_hi>=simend){
	draw_single=0;
    }
    double ck_0=myclockd();
    simu->isim=isim;
    simu->status->isim=isim;
    sim_update_etf(simu);
    if(parms->atm.frozenflow){
#if USE_CUDA
	if(parms->gpu.evl || parms->gpu.wfs){
	    /*may need to copy another part */
	    gpu_atm2gpu(simu->atm, simu->atmscale, parms, iseed, isim);
	}
#endif
    }else{
	//Do not put this one inside parallel 
	genatm(simu);
	/*re-seed the atmosphere in case atm is loaded from shm/file */
	seed_rand(simu->atm_rand, lrand(simu->init_rand));
    }
    OMPTASK_SINGLE{
	if(parms->sim.dmproj){
	    /* teporarily disable FR.M so that Mfun is used.*/
	    cell *FRM=recon->FR.M; recon->FR.M=NULL; 
	    muv_solve(&simu->dmproj, &recon->FL, &recon->FR, NULL);
	    recon->FR.M=FRM;/*set FR.M back*/
	    if(parms->save.dm){
		cellarr_dcell(simu->save->dmproj, simu->isim, simu->dmproj);
	    }
	    if(!parms->fit.square){
		/* Embed DM commands to a square array for fast ray tracing */
		for(int idm=0; idm<parms->ndm; idm++){
		    loc_embed(simu->dmprojsq->p[idm], recon->aloc->p[idm], simu->dmproj->p[idm]->p);
		}
	    }
#if USE_CUDA
	    if(parms->gpu.evl || parms->gpu.wfs){
		gpu_dmproj2gpu(simu->dmprojsq);
	    }
#endif
	}
	save_dmreal(simu);
	extern int NO_RECON, NO_WFS, NO_EVL;
	if(PARALLEL){
	    /*
	      We do the big loop in parallel to make better use the
	      CPUs. Notice that the reconstructor is working on grad from
	      last time step so that there is no confliction in data access.
	    */
	    /*when we want to apply idealngs correction, wfsgrad need to wait for perfevl. */
	    long group=0;
	    if(parms->gpu.evl && !NO_EVL){
		//Queue tasks on GPU, no stream sync is done
		QUEUE_THREAD(group, simu->perf_evl_pre, 0);
	    }
	    if(!parms->tomo.ahst_idealngs && parms->gpu.wfs && !NO_WFS){
		//task for each wfs
		QUEUE_THREAD(group, simu->wfs_grad_pre, 0);
	    }
	    if(!NO_RECON){
		//don't put this first. It has cpu overhead in computing gradol
		QUEUE(group, reconstruct, simu, 1, 0);
	    }
	    if(!NO_EVL){
		if(parms->gpu.evl){
		    //wait for GPU tasks to be queued before calling sync
		    WAIT(group);
		}
		QUEUE(group, perfevl, simu, 1, 0);
	    }
	    if(!NO_WFS){
		if(parms->tomo.ahst_idealngs || (parms->gpu.wfs && !parms->gpu.evl)){
		    //in ahst_idealngs mode, weight for perfevl to finish.
		    //otherwise, wait for GPU tasks to be queued before calling sync
		    WAIT(group);
		}
		QUEUE(group, wfsgrad, simu, 1, 0);
	    }
	    if(!NO_RECON){
		//wait for all tasks to finish before modifying dmreal
		WAIT(group);
		shift_grad(simu);/*before filter() */
		filter_dm(simu);/*updates dmreal, so has to be after prefevl/wfsgrad is done. */
	    }
	    WAIT(group);
	}else{/*do the big loop in serial mode. */
	    if(parms->sim.closeloop){
		if(!NO_EVL) perfevl(simu);/*before wfsgrad so we can apply ideal NGS modes */
		if(!NO_WFS) wfsgrad(simu);/*output grads to gradcl, gradol */
		if(!NO_RECON) {
		    reconstruct(simu);/*uses grads from gradlast cl, gradlast ol. */
		    shift_grad(simu);
		    filter_dm(simu);
		}
	    }else{/*in OL mode,  */
		if(!NO_WFS) wfsgrad(simu);
		if(!NO_RECON) {
		    shift_grad(simu);
		    reconstruct(simu);
		    filter_dm(simu);
		}
		if(!NO_EVL) perfevl(simu);
	    }
	}
    }
    double ck_end=myclockd();
    long steps_done=iseed*(simend-simstart)+(isim+1-simstart);
    long steps_rest=parms->sim.nseed*(simend-simstart)-steps_done;
    if(isim!=simstart){
	simu->status->rest=(long)((ck_end-tk_0-(tk_atm-tk_1)*(iseed+1))/steps_done*steps_rest
				  +(tk_atm-tk_1)*(parms->sim.nseed-iseed-1));
	simu->status->mean=(ck_end-tk_atm)/(double)(isim-simstart);
    }
    simu->status->laps=(long)(ck_end-tk_0);
    simu->status->tot  =ck_end-ck_0;
    simu->status->wfs  =simu->tk_wfs;
    simu->status->recon=simu->tk_recon;
    simu->status->other=simu->tk_cache;
    simu->status->eval =simu->tk_eval;
    simu->status->scale=1;
    if(simu->timing){
	simu->timing->p[isim*simu->timing->nx]=get_job_mem();
	simu->timing->p[isim*simu->timing->nx+1]=simu->status->tot;
	simu->timing->p[isim*simu->timing->nx+2]=simu->status->wfs;
	simu->timing->p[isim*simu->timing->nx+3]=simu->status->recon;
	simu->timing->p[isim*simu->timing->nx+4]=simu->status->eval;
    }
    double this_time=myclockd();
    if(this_time>simu->last_report_time+1 || isim+1==simend || parms->sim.pause){
	/*we don't print out or report too frequently. */
	simu->last_report_time=this_time;
#if defined(__linux__) || defined(__APPLE__)
	scheduler_report(simu->status);
#endif
	print_progress(simu);
    }
}

/**
   Closed loop simulation main loop. It calls init_simu() to initialize the
   simulation struct. Then calls genatm() to generate atmospheric turbulence
   screens. Then for every time step, it calls perfevl() to evaluate
   performance, wfsgrad() to do wfs measurement, reconstruct() to do tomography
   and DM fit, filter() to do DM command filtering. In MOAO mode, it call calls
   moao_recon() for MOAO DM fitting.  \callgraph */
void maos_sim(){
    const PARMS_T *parms=global->parms;
    POWFS_T *powfs=global->powfs;
    RECON_T *recon=global->recon;
    APER_T *aper=global->aper;
    int simend=parms->sim.end;
    int simstart=parms->sim.start;
    if(parms->sim.skysim){
	save_skyc(powfs,recon,parms);
    }
    if(parms->evl.psfmean || parms->evl.psfhist){
	/*compute diffraction limited PSF. Save to output directory.*/
	dmat *iopdevl=dnew(aper->locs->nloc,1);
	ccell *psf2s=0;
	locfft_psf(&psf2s, aper->embed, iopdevl, parms->evl.psfsize, 0);
	const int nwvl=parms->evl.nwvl;
	dcell *evlpsfdl=cellnew(nwvl,1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cabs22d(&evlpsfdl->p[iwvl], 1, psf2s->p[iwvl], 1);
	    evlpsfdl->p[iwvl]->header=evl_header(parms, aper, -1, iwvl);
	}
	ccellfree(psf2s);
	writebin(evlpsfdl, "evlpsfdl.fits");
	dcellfree(evlpsfdl);
	dfree(iopdevl);
    }
    info2("PARALLEL=%d\n", PARALLEL);
    if(simstart>=simend) return;
    double restot=0; long rescount=0;
    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	SIM_T *simu=0;
	while(!(simu=maos_iseed(iseed))){
	    iseed++;
	}
#ifdef HAVE_NUMA_H
	numa_set_localalloc();
#endif
	for(int isim=simstart; isim<simend; isim++){
	    maos_isim(isim);
	    if(parms->sim.pause){
		mypause();
	    }
	}/*isim */
	{
	    /*Compute average performance*/
	    int isim0;
	    if(parms->sim.closeloop){
		if(parms->sim.end>100){
		    isim0=MAX(50,parms->sim.end/10);
		}else{
		    isim0=MIN(20, parms->sim.end/2);
		}
	    }else{
		isim0=0;
	    }
	    double sum=0;
	    for(int i=isim0; i<parms->sim.end; i++){
		sum+=simu->cle->p[i*parms->evl.nmod];
	    }
	    restot+=sum/(parms->sim.end-isim0);
	    rescount++;
	}
	free_simu(simu);
	global->simu=0;
    }/*seed */
    printf("%g\n", sqrt(restot/rescount)*1e9);
}
