/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <search.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "maos.h"
#include "recon.h"
#include "recon_utils.h"
#include "setup_powfs.h"
#include "sim.h"
#include "sim_utils.h"
#include "fdpcg.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
#define TIMING_MEAN 0
extern int PARALLEL;
/**
   Closed loop simulation main loop. It calls init_simu() to initialize the
   simulation struct. Then calls genscreen() to generate atmospheric turbulence
   screens. Then for every time step, it calls perfevl() to evaluate
   performance, wfsgrad() to do wfs measurement, reconstruct() to do tomography
   and DM fit, filter() to do DM command filtering. In MOAO mode, it call calls
   moao_recon() for MOAO DM fitting.  \callgraph */
double tk_start;
void sim(const PARMS_T *parms,  POWFS_T *powfs, APER_T *aper,  RECON_T *recon){
    int simend=parms->sim.end;
    int simstart=parms->sim.start;
    if(parms->sim.skysim){
	save_skyc(powfs,recon,parms);
    }
    info2("PARALLEL=%d\n", PARALLEL);
    if(simstart>=simend) return;
    double tk_0=myclockd();
    double ck_0, ck_end;
    double restot=0; long rescount=0;
    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	if(parms->fdlock && parms->fdlock[iseed]<0){
	    warning("Another MAOS is already running. Skip seed %d\n", 
		    parms->sim.seeds[iseed]);
	    continue;
	}
	extern int draw_single;
	if(!parms->pause){
	    draw_single=1;//Only draw active frame.
	}else{
	    draw_single=0;
	}
	global->iseed=iseed;
	double tk_1=myclockd();
	SIM_T *simu=init_simu(parms,powfs,aper,recon,iseed);
	global->simu=simu;
	if(recon) recon->simu=simu;
	if(parms->atm.frozenflow){
	    genscreen(simu);/*Generating atmospheric screen(s) that frozen flows.*/
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
	    gpu_atm2gpu(simu->atm, parms, iseed, simstart);/*takes 0.4s for NFIRAOS. */
	    if(parms->tomo.predict){
		gpu_update_recon(parms, recon);
	    }
	}
#endif
	double tk_atm=myclockd();
	const int CL=parms->sim.closeloop;
	for(int isim=simstart; isim<simend; isim++){
	    if(isim+1==simend) draw_single=0;
	    ck_0=myclockd();
	    simu->isim=isim;
	    simu->status->isim=isim;
	    sim_update_etf(simu);
	    if(parms->atm.frozenflow){
		if(parms->atm.evolve){
		    evolve_screen(simu);
		}
#if USE_CUDA
		if(parms->gpu.evl || parms->gpu.wfs){/*may need to copy another part */
		    gpu_atm2gpu(simu->atm, parms, iseed, isim);/*takes 0.4s for NFIRAOS 64 meter screen. */
		}
#endif
	    }else{
		genscreen(simu);
		/*re-seed the atmosphere in case atm is loaded from shm/file */
		seed_rand(simu->atm_rand, lrand(simu->init_rand));
	    }
	    if(parms->sim.dmproj){
		/* teporarily disable FR.M so that Mfun is used.*/
		spcell *FRM=recon->FR.M; recon->FR.M=NULL; 
		muv_solve(&simu->dmproj, &recon->FL, &recon->FR, NULL);
		recon->FR.M=FRM;/*set FR.M back*/
		if(parms->save.dm){
		    cellarr_dcell(simu->save->dmproj, simu->isim, simu->dmproj);
		}
		if(!parms->fit.square){
		    /* Embed DM commands to a square array for fast ray tracing */
		    for(int idm=0; idm<parms->ndm; idm++){
			loc_embed(simu->dmprojsq[idm], recon->aloc[idm], simu->dmproj->p[idm]->p);
		    }
		}
#if USE_CUDA
		if(parms->gpu.evl || parms->gpu.wfs){
		    gpu_dmproj2gpu(simu->dmprojsq, parms->ndm, NULL);
		}
#endif
	    }
	    extern int NO_RECON, NO_WFS, NO_EVL;
	    if(PARALLEL){
		tk_start=myclockd();
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
		    filter(simu);/*updates dmreal, so has to be after prefevl/wfsgrad is done. */
		}
		WAIT(group);
	    }else{/*do the big loop in serial mode. */
		if(CL){
		    tk_start=myclockd();
		    if(!NO_EVL) perfevl(simu);/*before wfsgrad so we can apply ideal NGS modes */
		    tk_start=myclockd();
		    if(!NO_WFS) wfsgrad(simu);/*output grads to gradcl, gradol */
		    tk_start=myclockd();
		    if(!NO_RECON) {
			reconstruct(simu);/*uses grads from gradlast cl, gradlast ol. */
			shift_grad(simu);
			filter(simu);
		    }
		}else{/*in OL mode,  */
		    tk_start=myclockd();
		    if(!NO_WFS) wfsgrad(simu);
		    tk_start=myclockd();
		    if(!NO_RECON) {
			shift_grad(simu);
			reconstruct(simu);
			filter(simu);
		    }
		    tk_start=myclockd();
		    if(!NO_EVL) perfevl(simu);
		}
	    }
	    ck_end=myclockd();
	    long steps_done=iseed*(simend-simstart)+(isim+1-simstart);
	    long steps_rest=parms->sim.nseed*(simend-simstart)-steps_done;
	    simu->status->rest=(long)((ck_end-tk_0-(tk_atm-tk_1)*(iseed+1))/steps_done*steps_rest
				      +(tk_atm-tk_1)*(parms->sim.nseed-iseed-1));
	    simu->status->laps=(long)(ck_end-tk_0);
	    simu->status->mean=(ck_end-tk_atm)/(double)(isim+1-simstart);
	    simu->status->tot  =ck_end-ck_0;
	    simu->status->wfs  =simu->tk_wfs;
	    simu->status->recon=simu->tk_recon;
	    simu->status->other=simu->tk_cache;
	    simu->status->eval =simu->tk_eval;
	    simu->status->scale=1;

	    int this_time=myclocki();
	    if(this_time>simu->last_report_time+1 || isim+1==simend){
		/*we don't print out or report too frequently. */
		simu->last_report_time=this_time;
#if defined(__linux__) || defined(__APPLE__)
		scheduler_report(simu->status);
#endif
		print_progress(simu);
	    }
	    if(parms->pause){
		mypause();
	    }
	}/*isim */
	{
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
    }/*seed */
    printf("%g\n", sqrt(restot/rescount)*1e9);
}
