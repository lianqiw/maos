/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "maos.h"
#include "moao.h"
/**
   \file maos.c
   Contains main() and the entry into simulation maos()
*/
GLOBAL_T *global=NULL;//record for convenient access.
int use_cuda=0;
const char *dirskysim=NULL;

/** begin variable overridable by environment variable MAOS_ .  These are for
 debugging maos itself. Not pertinent to a particular simulation*/
double TOMOSCALE=1e-12;
int PARALLEL=1; //DO wfs, evl, and recon in parallel
int KEEP_MEM=0; //keep allocated memory during steps.
int NO_WFS=0;
int NO_EVL=0;
int NO_RECON=0;
/** end*/
static void read_env(){
    READ_ENV_DBL(TOMOSCALE,0,INFINITY);
    READ_ENV_INT(PARALLEL,0,1);
    READ_ENV_INT(NO_WFS,0,1);
    READ_ENV_INT(NO_EVL,0,1);
    READ_ENV_INT(NO_RECON,0,1);
    READ_ENV_INT(KEEP_MEM,0,1);
    info("TOMOSCALE=%g\n", TOMOSCALE);
}
void maos_setup(const PARMS_T *parms){
    TIC;tic;
    global=mycalloc(1,GLOBAL_T);
    global->parms=parms;
    APER_T  * aper=NULL;
    POWFS_T * powfs=NULL;
    RECON_T * recon=NULL;
    read_env();
    if(parms->sim.skysim){
	dirskysim="skysim";
	mymkdir("%s",dirskysim);
    }else{
	dirskysim=".";
    }
    if(parms->save.setup){
	mymkdir("setup");
	if(chdir("setup")){
	    error("Unable to save to folder setup\n");
	}
    }
 
    THREAD_POOL_INIT(NTHREAD);
    global->aper=aper=setup_aper(parms);
    info("After setup_aper:\t%.2f MiB\n",get_job_mem()/1024.);
    if(!parms->sim.evlol){
	global->powfs=powfs=setup_powfs_init(parms, aper);
	info("After setup_powfs:\t%.2f MiB\n",get_job_mem()/1024.);
	/*Setup DM fitting parameters so we can flatten the DM in setup_surf.c */
	global->recon=recon=setup_recon_prep(parms, aper, powfs);
	/*setting up M1/M2/M3, Instrument, Lenslet surface OPD. DM Calibration, WFS bias.*/
	setup_surf(parms, aper, powfs, recon);
	/*set up physical optics powfs data. It needs dmncpa and wfsadd from setup_surf()*/
	setup_powfs_phy(parms, powfs);
	/*calibrate gradient offset*/
	setup_powfs_calib(parms, powfs);
    }
    if(parms->plot.setup){
	plot_setup(parms, powfs, aper, recon);
    }
    global->setupdone=1;
    if(!parms->sim.evlol){
#if USE_CUDA
	extern int cuda_dedup;
	cuda_dedup=1;
	if(!parms->sim.evlol && (parms->gpu.wfs)){// || parms->gpu.tomo)){
	    gpu_wfsgrad_init(parms, powfs);
	}
	
#endif
	setup_recon_prep2(recon, parms, aper, powfs);
	//Don't put this inside parallel, otherwise svd will run single threaded.
	setup_recon(recon, parms, powfs);
	if(parms->recon.alg==0 || parms->sim.dmproj){
	    setup_recon_fit(recon, parms);
	}
	if(parms->recon.alg==0 && parms->nmoao){
	    setup_recon_moao(recon, parms);
	}
	info("After setup_recon:\t%.2f MiB\n",get_job_mem()/1024.);
	if(parms->dbg.wfslinearity!=-1){
	    int iwfs=parms->dbg.wfslinearity;
	    assert(iwfs>-1 || iwfs<parms->nwfs);
	    info("Studying wfslineariy for WFS %d\n", iwfs);
	    wfslinearity(parms, powfs, iwfs);
	    ((PARMS_T*)parms)->sim.end=parms->sim.start;//indicate no simulation
	}
	{
	    int LGS_SPH_PSD=-1;
	    READ_ENV_INT(LGS_SPH_PSD,-1,INFINITY);
	    if(LGS_SPH_PSD>-1){
		lgs_wfs_sph_psd(parms, powfs, recon, LGS_SPH_PSD);
		((PARMS_T*)parms)->sim.end=parms->sim.start;//indicate no simulation
	    }
	}
    }
#if USE_CUDA
    if(parms->gpu.evl){
	gpu_perfevl_init(parms, aper);
    }
    if(parms->gpu.wfs && powfs){
	gpu_wfssurf2gpu(parms, powfs);
    }
    if(!parms->sim.evlol && (parms->gpu.tomo || parms->gpu.fit || parms->gpu.lsr)){
	gpu_setup_recon(parms, powfs, recon);
    }
    extern int cuda_dedup;
    cuda_dedup=0;
#endif

    if(!parms->sim.evlol && parms->recon.mvm){
	setup_recon_mvm(parms, recon, powfs);
    }
    setup_recon_post(recon, parms, aper);
    if(parms->save.setup && chdir("..")){}
    /*
      Before entering real simulation, make sure to delete all variables that
      won't be used later on to save memory.
    */
    free_powfs_unused(parms, powfs);
    //free_recon_unused(parms, recon);
    toc2("Presimulation");
}

void maos_reset(){
    /*Free all allocated memory in setup_* functions. So that we
      keep track of all the memory allocation.*/
    if(!global) return;
    PARMS_T *parms=(PARMS_T*)global->parms;
    free_simu(global->simu);
    free_recon(parms,global->recon); 
    free_powfs(parms,global->powfs); 
    free_aper(global->aper);
    free_parms(parms);
#if USE_CUDA
    if(use_cuda){
	gpu_cleanup();
    }
#endif
    free(global); global=0;
}
