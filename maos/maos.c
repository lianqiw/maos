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
#include "maos.h"
#include "moao.h"
/**
   \file maos.c
   Contains main() and the entry into simulation maos()
*/
GLOBAL_T *global=NULL;//record for convenient access.
int use_cuda=0;
const char *dirsetup=NULL;
const char *dirskysim=NULL;

/** begin variable overridable by environment variable MAOS_ .  These are for
 debugging maos itself. Not pertinent to a particular simulation*/
double TOMOSCALE=1e-12;
int PARALLEL=1;
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
    info2("TOMOSCALE=%g\n", TOMOSCALE);
}
void maos_setup(const PARMS_T *parms){
    TIC;tic;
    global=calloc(1, sizeof(GLOBAL_T));
    global->parms=parms;
    APER_T  * aper=NULL;
    POWFS_T * powfs=NULL;
    RECON_T * recon=NULL;
    read_env();
    if(parms->save.setup){
	dirsetup="setup";
	mymkdir("%s",dirsetup);
    }else{
	dirsetup=".";
    }
    if(parms->sim.skysim){
	dirskysim="skysim";
	mymkdir("%s",dirskysim);
    }else{
	dirskysim=".";
    }
    THREAD_POOL_INIT(NTHREAD);
    global->aper=aper=setup_aper(parms);
    info2("After setup_aper:\t%.2f MiB\n",get_job_mem()/1024.);
    if(!parms->sim.evlol){
	global->powfs=powfs=setup_powfs_init(parms, aper);
	info2("After setup_powfs:\t%.2f MiB\n",get_job_mem()/1024.);
	/*Setup DM fitting parameters so we can flatten the DM in setup_surf.c */
	global->recon=recon=setup_recon_init(parms, aper);
	/*setting up M1/M2/M3, Instrument, Lenslet surface OPD. DM Calibration, WFS bias.*/
	setup_surf(parms, aper, powfs, recon);
	/*set up physical optics powfs data*/
	setup_powfs_phy(parms, powfs);
	//Don't put this inside parallel, otherwise svd will run single threaded.
	setup_recon(recon, parms, powfs, aper);
	if(parms->recon.alg==0 || parms->sim.dmproj){
	    setup_recon_fit(recon, parms);
	}
	if(parms->recon.alg==0 && parms->nmoao){
	    setup_recon_moao(recon, parms);
	}
	info2("After setup_recon:\t%.2f MiB\n",get_job_mem()/1024.);
	if(parms->dbg.wfslinearity!=-1){
	    int iwfs=parms->dbg.wfslinearity;
	    assert(iwfs>-1 || iwfs<parms->nwfs);
	    info2("Studying wfslineariy for WFS %d\n", iwfs);
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
    global->setupdone=1;
    if(parms->plot.setup){
	plot_setup(parms, powfs, aper, recon);
    }
#if USE_CUDA
    extern int cuda_dedup;
    cuda_dedup=1;
    if(!parms->sim.evlol && (parms->gpu.wfs || parms->gpu.tomo)){
	gpu_wfsgrad_init(parms, powfs);
    }
    if(parms->gpu.evl){
	gpu_perfevl_init(parms, aper);
    }
    if(parms->gpu.wfs && powfs){
	gpu_wfssurf2gpu(parms, powfs);
    }
    if(!parms->sim.evlol && (parms->gpu.tomo || parms->gpu.fit)){
	gpu_setup_recon(parms, powfs, recon);
    }
    cuda_dedup=0;
#endif

    if(!parms->sim.evlol && parms->recon.mvm){
	setup_recon_mvm(parms, recon, powfs);
    }
    /*
      Before entering real simulation, make sure to delete all variables that
      won't be used later on to save memory.
    */
    //free_powfs_unused(parms, powfs);
    //free_recon_unused(parms, recon);
    toc2("Presimulation");
}
/**
   This is the routine that calls various functions to do the simulation. maos()
   calls setup_aper(), setup_powfs(), and setup_recon() to set up the aperture
   (of type APER_T), wfs (of type POWFS_T), and reconstructor (of type RECON_T)
   structs and then hands control to sim(), which then stars the simulation.
   \callgraph */
void maos(const PARMS_T *parms){
    info2("\n*** Preparation started at %s in %s. ***\n\n",myasctime(),myhostname());
    maos_setup(parms);
    if(parms->sim.end>parms->sim.start){
	info2("\n*** Simulation  started at %s in %s. ***\n\n",myasctime(),myhostname());
	maos_sim();
    }
    maos_reset();
    info2("\n*** Simulation finished at %s in %s. ***\n\n",myasctime(),myhostname());
}

void maos_reset(){
    /*Free all allocated memory in setup_* functions. So that we
      keep track of all the memory allocation.*/
    if(!global) return;
    PARMS_T *parms=(PARMS_T*)global->parms;
    free_simu(global->simu);
    free_recon(parms,global->recon); 
    free_powfs(parms,global->powfs); 
    free_aper(parms, global->aper);
    free_parms(parms);
#if USE_CUDA
    if(use_cuda){
	gpu_cleanup();
    }
#endif
    free(global); global=0;
}
