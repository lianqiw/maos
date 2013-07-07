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
#include <assert.h>
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
/**
   Closed loop simulation main loop. It calls init_simu() to initialize the
   simulation struct. Then calls genscreen() to generate atmospheric turbulence
   screens. Then for every time step, it calls perfevl() to evaluate
   performance, wfsgrad() to do wfs measurement, reconstruct() to do tomography
   and DM fit, filter() to do DM command filtering. In MOAO mode, it call calls
   moao_recon() for MOAO DM fitting.  \callgraph */

void sim(const PARMS_T *parms,  POWFS_T *powfs, 
	 APER_T *aper,  RECON_T *recon){
    int simend=parms->sim.end;
    int simstart=parms->sim.start;
    if(parms->sim.skysim){
	save_skyc(powfs,recon,parms);
    }
    if(simstart>=simend) return;
    double tk_0=myclockd();
    double ck_0,ck_end;
    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	extern int draw_single;
	draw_single=1;//Only draw active frame.
	global->iseed=iseed;
	double tk_start=myclockd();
	SIM_T *simu=init_simu(parms,powfs,aper,recon,iseed);
	global->simu=simu;
	if(!simu) continue;/*skip. */
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
		gpu_setup_recon_predict(parms, recon);
		gpu_setup_recon_fdpcg(parms, recon);
	    }
	}
#endif
	double tk_atm=myclockd();
	const int CL=parms->sim.closeloop;
#if defined(__linux__)  && !defined(USE_CUDA)
	double cpu_evl=0, cpu_wfs=0, cpu_recon=0, cpu_cachedm=0, cpu_all=0;
#define READ_CPU(A) A=read_self_cpu()
#else
#define READ_CPU(A)
#endif
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
			long *embed=simu->recon->aembed[idm];
			double *pout=simu->dmprojsq[idm]->p;
			double *pin=simu->dmproj->p[idm]->p;
			for(long i=0; i<simu->dmproj->p[idm]->nx; i++){
			    pout[embed[i]]=pin[i];
			}
		    }
		}
#if USE_CUDA
		if(parms->gpu.evl || parms->gpu.wfs){
		    gpu_dmproj2gpu(simu->dmprojsq, parms->ndm, NULL);
		}
#endif
	    }
	    if(parms->sim.parallel){
		/*
		  We do the big loop in parallel to make better use the
		  CPUs. Notice that the reconstructor is working on grad from
		  last time step so that there is no confliction in data access.
		  This parallelized mode can not be used in skycoverage
		  presimulation case because perfevl updates dmreal before
		  wfsgrad
		*/
		read_self_cpu();/*initialize CPU usage counter */
		/*when we want to apply idealngs correction, wfsgrad need to wait for perfevl. */
		long group=0;
		thread_pool_queue(&group, (thread_fun)reconstruct, simu, 0);
		thread_pool_queue(&group, (thread_fun)perfevl, simu, 0);
		if(parms->tomo.ahst_idealngs)/*Need to wait until perfevl is done to start wfsgrad. */
		    thread_pool_wait(&group);
		thread_pool_queue(&group, (thread_fun)wfsgrad, simu, 0);
		thread_pool_wait(&group);
		shift_grad(simu);/*before filter() */
		filter(simu);/*updates dmreal, so has to be after prefevl/wfsgrad is done. */
		READ_CPU(cpu_all);
	    }else{/*do the big loop in serial mode. */
		read_self_cpu();/*initialize CPU usage counter */
		if(CL){
		    perfevl(simu);/*before wfsgrad so we can apply ideal NGS modes */
		    READ_CPU(cpu_evl);
		    wfsgrad(simu);/*output grads to gradcl, gradol */
		    READ_CPU(cpu_wfs);
		    reconstruct(simu);/*uses grads from gradlast cl, gradlast ol. */
		    READ_CPU(cpu_recon);
		    shift_grad(simu);
		    filter(simu);
		    READ_CPU(cpu_cachedm);
		}else{/*in OL mode,  */
		    wfsgrad(simu);
		    READ_CPU(cpu_wfs);
		    shift_grad(simu);
		    reconstruct(simu);
		    READ_CPU(cpu_recon);
		    filter(simu);
		    READ_CPU(cpu_cachedm);
		    perfevl(simu);
		    READ_CPU(cpu_evl);
		}
	    }
	    ck_end=myclockd();
	    long steps_done=iseed*(simend-simstart)+(isim+1-simstart);
	    long steps_rest=parms->sim.nseed*(simend-simstart)-steps_done;
	    simu->status->rest=(long)((ck_end-tk_0-(tk_atm-tk_start)*(iseed+1))/steps_done*steps_rest
				      +(tk_atm-tk_start)*(parms->sim.nseed-iseed-1));
	    simu->status->laps=(long)(ck_end-tk_0);
	    simu->status->mean=(ck_end-tk_atm)/(double)(isim+1-simstart);
	    simu->status->tot  =ck_end-ck_0;
	    simu->status->wfs  =simu->tk_wfs;
	    simu->status->recon=simu->tk_recon;
	    simu->status->cache=simu->tk_cache;
	    simu->status->eval =simu->tk_eval;
	    simu->status->scale=1;

	    int this_time=myclocki();
	    if(this_time>simu->last_report_time+1 || isim+1==simend){
		/*we don't print out or report too frequenctly. */
		simu->last_report_time=this_time;
#if defined(__linux__) || defined(__APPLE__)
		scheduler_report(simu->status);
#endif
		print_progress(simu);
#if defined(__linux__)  && !defined(USE_CUDA)
		    if(parms->sim.parallel){
			info2("CPU Usage: %.2f\n", cpu_all);
		    }else{
			info2("CPU Usage: WFS:%.2f Recon:%.2f CACHE: %.2f EVAL:%.2f Mean:%.2f\n",
			      cpu_wfs, cpu_recon, cpu_evl, cpu_cachedm,
			      (cpu_wfs+cpu_recon+cpu_evl+cpu_cachedm)*0.25);
		    }
#endif
	    }
	    if(parms->pause){
		PAUSE;
	    }
	}/*isim */
	free_simu(simu);
    }/*seed */
}
