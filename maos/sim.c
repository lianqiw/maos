/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <alloca.h>
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
#include "setup_powfs.h"
#include "sim.h"
#include "sim_utils.h"
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
	curiseed=iseed;
	double tk_start=myclockd();
	SIM_T *simu=init_simu(parms,powfs,aper,recon,iseed);
	if(!simu) continue;//skip.
	recon->simu=simu;
	if(parms->atm.frozenflow){
	    genscreen(simu);/*Generating atmospheric screen(s) that frozen flows.*/
	}
	double tk_atm=myclockd();
	const int CL=parms->sim.closeloop;
	for(int isim=simstart; isim<simend; isim++){
	    ck_0=myclockd();
	    simu->isim=isim;
	    simu->status->isim=isim;
	    sim_update_etf(simu);
	    if(parms->atm.frozenflow){
		if(parms->atm.evolve){
		    evolve_screen(simu);
		}
	    }else{
		disable_atm_shm=1;
		genscreen(simu);
		//re-seed the atmosphere in case atm is loaded from shm/file
		seed_rand(simu->atm_rand, lrand(simu->init));
	    }
	    if(parms->dbg.parallel){
		/*
		  We do the big loop in parallel to make better use the
		  CPUs. Notice that the reconstructor is working on grad from
		  last time step so that there is no confliction in data access.
		  This parallelized mode can not be used in skycoverage
		  presimulation case because perfevl updates dmreal before
		  wfsgrad
		*/
		read_self_cpu();//initialize CPU usage counter
		if(parms->tomo.ahst_idealngs){
		    //when we want to apply idealngs correction, wfsgrad need to wait for perfevl.
		    long group1=0, group2=0;
		    thread_pool_queue(&group1, (thread_fun)reconstruct, simu, 0);
		    thread_pool_queue(&group2, (thread_fun)perfevl, simu, 0);
		    thread_pool_wait(&group2);
		    thread_pool_queue(&group1, (thread_fun)wfsgrad, simu, 0);
		    thread_pool_wait(&group1);
		}else{
		    long group=0;
		    thread_pool_queue(&group, (thread_fun)reconstruct, simu, 0);
		    thread_pool_queue(&group, (thread_fun)perfevl, simu, 0);
		    thread_pool_queue(&group, (thread_fun)wfsgrad, simu, 0);
		    thread_pool_wait(&group);
		}
		filter(simu);//updates dmreal, so has to be after prefevl/wfsgrad is done.
		dcellcp(&simu->gradlastcl, simu->gradcl);
		dcellcp(&simu->gradlastol, simu->gradol);
#if defined(__linux__)
		if(simu->nthread>1 && !detached){
		    info2("CPU Usage: %.2f\n", read_self_cpu());
		}
#endif
	    }else{//do the big loop in serial mode.
		double cpu_evl, cpu_wfs, cpu_recon, cpu_cachedm;
		read_self_cpu();//initialize CPU usage counter
		if(CL){
		    perfevl(simu);//before wfsgrad so we can apply ideal NGS modes
		    cpu_evl=read_self_cpu();
		    wfsgrad(simu);//output grads to gradcl, gradol
		    cpu_wfs=read_self_cpu();
		    reconstruct(simu);//uses grads from gradlast cl, gradlast ol.
		    cpu_recon=read_self_cpu();
		    filter(simu);
		    cpu_cachedm=read_self_cpu();
		    dcellcp(&simu->gradlastcl, simu->gradcl);
		    dcellcp(&simu->gradlastol, simu->gradol);
		}else{//in OL mode, 
		    wfsgrad(simu);
		    cpu_wfs=read_self_cpu();
		    dcellcp(&simu->gradlastcl, simu->gradcl);
		    dcellcp(&simu->gradlastol, simu->gradol);
		    reconstruct(simu);
		    cpu_recon=read_self_cpu();
		    filter(simu);
		    cpu_cachedm=read_self_cpu();
		    perfevl(simu);
		    cpu_evl=read_self_cpu();
		}
#if defined(__linux__)
		if(simu->nthread>1 && !detached){
		    fprintf(stderr,"CPU Usage: WFS:%.2f Recon:%.2f CACHE: %.2f"
			    " EVAL:%.2f Mean:%.2f\n",
			    cpu_wfs, cpu_recon, cpu_evl, cpu_cachedm,
			    (cpu_wfs+cpu_recon+cpu_evl+cpu_cachedm)*0.25);
		}
#endif
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
		//we don't print out or report too frequenctly.
		simu->last_report_time=this_time;
#if defined(__linux__) || defined(__APPLE__)
		scheduler_report(simu->status);
#endif
	    }
	    print_progress(simu);
	    if(parms->pause){//does not work.
		info2("Press Enter to continue.\n");
		while(getchar()!=0x0a);
	    }
	}//isim
	free_simu(simu);
    }
}
