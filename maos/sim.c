/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "moao.h"
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

    for(int iseed=0; iseed<parms->sim.nseed; iseed++){
	double tk_1=myclockd();
	SIM_T *simu=init_simu(parms,powfs,aper,recon,iseed);
	if(!simu) continue;//skip.
	if(parms->sim.frozenflow){
	    /*Generating atmospheric screen(s) that frozen flows.*/
	    genscreen(simu);
	}
	double tk_2=myclockd();
	const int CL=parms->sim.closeloop;
	for(int isim=simstart; isim<simend; isim++){
	    read_self_cpu();//initialize CPU usage counter
	    simu->isim=isim;
	    simu->status->isim=isim;
	    sim_update_etf(simu);
	    double ck_0=myclockd();
	    if(!parms->sim.frozenflow){
		disable_atm_shm=1;
		genscreen(simu);
		//re-seed the atmosphere in case atm is loaded from shm
		seed_rand(simu->atm_rand, lrand(simu->init));
	    }
#define PARALLEL 1
#if PARALLEL == 1 //does not help in T410s. Need to test taurus
	    pthread_t thread_perfevl, thread_wfsgrad, thread_tomofit;
	    if(CL){
		pthread_create(&thread_perfevl, NULL, (void*(*)(void*))perfevl, (void*)simu);
	    }
	    pthread_create(&thread_wfsgrad, NULL, (void*(*)(void*))wfsgrad, (void*)simu);
	    pthread_join(thread_wfsgrad,NULL);
	    pthread_create(&thread_tomofit, NULL, (void*(*)(void*))reconstruct, (void*)simu);
	    pthread_join(thread_tomofit,NULL);
	    if(CL){
		pthread_join(thread_perfevl,NULL);
	    }
	    
	    double ck_1=myclockd();
	    double ck_2=myclockd();
	    double ck_3=myclockd();
#else  
	    if(CL){//before wfsgrad so we can apply ideal NGS modes
		perfevl(simu);
	    }
	    double ck_1=myclockd(); double cpu_1=read_self_cpu();
	    wfsgrad(simu);
	    double ck_2=myclockd(); double cpu_2=read_self_cpu();
	    reconstruct(simu);
	    double ck_3=myclockd(); double cpu_3=read_self_cpu();
#endif
	    filter(simu);//updates dmreal.
	    if(recon->moao){
		moao_recon(simu);
	    }
	    double ck_4=myclockd(); 
#if PARALLEL==0
	    double cpu_4=read_self_cpu();
#endif
	    if(!CL){//Only in Open loop
		perfevl(simu);
	    }
	    save_simu(simu);
	    double ck_5=myclockd();
	    long steps_done=iseed*(simend-simstart)+(isim+1-simstart);
	    long steps_rest=parms->sim.nseed*(simend-simstart)-steps_done;
	    simu->status->rest=(long)((ck_5-tk_0-(tk_2-tk_1)*(iseed+1))/steps_done*steps_rest
				      +(tk_2-tk_1)*(parms->sim.nseed-iseed-1));
	    simu->status->laps=(long)(ck_5-tk_0);
	    simu->status->mean=(ck_5-tk_2)/(double)(isim+1-simstart);
	    simu->status->wfs  =ck_2-ck_1;
	    simu->status->recon=ck_3-ck_2;
	    simu->status->cache=ck_4-ck_3;
	    simu->status->eval =ck_1-ck_0+ck_5-ck_4;
	    simu->status->tot  =ck_5-ck_0;
	    simu->status->scale=1;

	    int this_time=myclocki();
	    if(this_time>simu->last_report_time+1 || isim+1==simend){
		//we don't print out or report too frequency.
		simu->last_report_time=this_time;
#if defined(__linux__) || defined(__APPLE__)
		scheduler_report(simu->status);
#endif
	    }
	    print_progress(simu);
#if defined(__linux__) && PARALLEL == 0
	    if(simu->nthread>1 && !detached){
		fprintf(stderr, "CPU Usage: WFS:%.2f Recon:%.2f CACHE:%.2f EVAL:%.2f Mean:%.2f\n",
			cpu_2, cpu_3, cpu_4, cpu_1, (cpu_1+cpu_2+cpu_3+cpu_4)*0.25);
	    }
#endif
	}//isim
	free_simu(simu);
    }
}
