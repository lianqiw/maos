/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include <unistd.h>
#include "maos.h"
#include "moao.h"
#include "ahst.h"
#include "sim.h"
/**
  \file setup.c
  Sets up maos simulation.
*/

global_t* global=NULL;//record for convenient access. It enables calling maos from matlab
int use_cuda=0;
const char *dirskysim=NULL;

/** begin variable overridable by environment variable MAOS_ .  These are for
 debugging maos itself. Not pertinent to a particular simulation*/
real TOMOSCALE=1e-12;//Without this, NEA inverse overflows in float point number mode.
int PARALLEL=1; //DO wfs, evl, and recon in parallel
int KEEP_MEM=0; //keep allocated memory during steps.
int NO_WFS=0;
int NO_EVL=0;
int NO_RECON=0;
/** end*/
static void read_sim_env(){
	READ_ENV_DBL(TOMOSCALE, 0, INFINITY);
	READ_ENV_INT(PARALLEL, 0, 2);
	READ_ENV_INT(NO_WFS, 0, 1);
	READ_ENV_INT(NO_EVL, 0, 1);
	READ_ENV_INT(NO_RECON, 0, 1);
	READ_ENV_INT(KEEP_MEM, 0, 1);
	dbg("TOMOSCALE=%g\n", TOMOSCALE);
}

/**
   Setup system before entering simulation.

   \callgraph
 */
void maos_setup(const parms_t* parms){
	TIC;tic;
	global=mycalloc(1, global_t);
	global->parms=parms;
	aper_t* aper=NULL;
	powfs_t* powfs=NULL;
	recon_t* recon=NULL;
	read_sim_env();
	if(PARALLEL&&(parms->sim.closeloop==0||parms->evl.tomo)){
		PARALLEL=0;	/*need to disable parallelizing the big loop. */
	}
	//dmproj/non-frozen flow does not work with PRALLEL=1 implementation
	if(PARALLEL==2&&parms->sim.dmproj){
		warning("dmproj does not work with PRALLEL=2 implementation. Set to 1.\n");
		PARALLEL=1;
	}
#ifndef _OPENMP
	if(PARALLEL==2){
		warning("PRALLEL=2 requires OpenMP. Set to 1.\n");
		PARALLEL=1;
	}
#endif

	if(parms->sim.skysim){
		dirskysim="skysim";
		mymkdir("%s", dirskysim);
	} else{
		dirskysim=".";
	}
	if(parms->save.setup){
		mymkdir("setup");
		if(chdir("setup")){
			error("Unable to save to folder setup\n");
		}
	}
#if USE_CUDA
	extern int cuda_dedup;
	cuda_dedup=1;
	
	if(use_cuda){
#if USE_CUDA>100
		gpu_ext_assign();
#endif		
		extern dmat *(*pywfs_mkg_ext)(const pywfs_t*pywfs, const loc_t*locin, const loc_t*locfft, const dmat*mod, real displacex, real displacey);
		pywfs_mkg_ext=gpu_pywfs_mkg;
	}
#endif
	global->aper=aper=setup_aper(parms);
	print_mem("After setup_aper");
#if USE_CUDA
	int has_pywfs=0;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].nwfs && parms->powfs[ipowfs].type==WFS_PY){
			has_pywfs=1;
		}
	}
#endif
	if(!parms->sim.evlol){
		global->powfs=powfs=setup_powfs_init(parms, parms->aper.amp);
		print_mem("After setup_powfs_init");
		//Setup geometry and DM fitting parameters so we can flatten the DM in setup_surf.c
		global->recon=recon=setup_recon_prep(parms, parms->aper.amp, powfs);
		print_mem("After setup_recon_prep");
		//pywfs_test(parms, powfs, recon);//as needed. needs recon->amod
		//setting up M1/M2/M3, Instrument, Lenslet surface OPD. DM Calibration, WFS bias.
		setup_surf(parms, aper, powfs, recon);
		print_mem("After setup_surf");
		//set up physical optics powfs data. It needs dmncpa and wfsadd from setup_surf()
		setup_shwfs_phy(parms, powfs);
		//setup gradient noise during simulation.
		setup_powfs_neasim(parms, powfs);
		//calibrate gradient offset of NCPA
		setup_powfs_calib(parms, powfs);
#if USE_CUDA
		if(parms->gpu.wfs&&powfs&&has_pywfs){
			gpu_wfsgrad_init(parms, powfs);//needed by pywfs_mkg
		}
#endif
		print_mem("After setup_powfs");
		//creates DM to WFS IA. needs GPU for pwfs. create amod for modal control.
		
		info_green("\nSetting up reconstructor\n\n");
		setup_recon_GA(recon, parms, powfs);//PWFS uses GPU data.
		setup_recon_GF(recon, parms);//GF depends on GA.
		setup_recon_GR(recon, parms);
		if(parms->recon.split||parms->evl.split){
			ngsmod_prep(parms, recon, aper);//needs GA
		}
		//assemble subaperture noise equivalent angle (SANEA) from powfs information
		//The following solvers needs SANEA
		setup_recon_saneai(recon, parms, powfs);
		setup_recon_dither_dm(recon, powfs, parms);//depends on saneai
		if(!NO_RECON){
			setup_recon_control(recon, parms, powfs);
			if(parms->recon.alg==RECON_MVR&&parms->nmoao){
				setup_recon_moao(recon, parms);
			}
		}
		if(parms->sim.wfsalias==2||parms->sim.idealwfs==2){
			setup_powfs_fit(powfs, recon, parms);
		}
		print_mem("After setup_recon");
		{//testing
			if(parms->dbg.wfslinearity!=-1){//testing wfs linearity
				int iwfs=parms->dbg.wfslinearity;
				assert(iwfs>-1||iwfs<parms->nwfs);
				info2("Studying wfslineariy for WFS %d\n", iwfs);
				wfslinearity(parms, powfs, iwfs);
				((parms_t*)parms)->sim.end=parms->sim.start;//indicate no simulation
			}
			int LGS_SPH_PSD=-1;
			READ_ENV_INT(LGS_SPH_PSD, -1, INFINITY);
			if(LGS_SPH_PSD>-1){
				lgs_wfs_sph_psd(parms, powfs, recon, LGS_SPH_PSD);
				((parms_t*)parms)->sim.end=parms->sim.start;//indicate no simulation
			}
		}
	}
	//Avoid initializing WFS/Perfevl data in GPU before this to preserve GPU memory
	if(!parms->sim.evlol&&parms->recon.mvm){
#if USE_CUDA
		if(parms->gpu.tomo&&parms->gpu.fit){
			gpu_setup_recon_mvm(parms, recon);
		}
#endif
		setup_recon_mvm(parms, recon, powfs);//use cpu to compute mvm or do the saving
	}

	if(parms->aper.misregu && (P(parms->aper.misregu,0) || P(parms->aper.misregu,1))){
		//Setup un-calibrated misregistration
		for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
			setup_powfs_amp(powfs, parms, parms->aper.amp, parms->aper.misregu, ipowfs);
		}
	}
#if USE_CUDA
	//setup_recon first because MVM assembly and transpose uses a lot of memory.
	if(parms->gpu.wfs&&powfs&&!has_pywfs){
		gpu_wfsgrad_init(parms, powfs);
	}
	if(parms->gpu.recon && !(parms->recon.mvm && parms->sim.mvmport)){
		gpu_setup_recon(parms, recon);
	}
	if(parms->gpu.evl){
		gpu_perfevl_init(parms, aper);
	}
	if(parms->gpu.wfs&&powfs){
		gpu_wfssurf2gpu(parms, powfs);
	}
#endif

	if(!parms->sim.evlol){
		setup_recon_misc(recon, parms, aper->locs, powfs);//needs MVM matrix
	}
	if(parms->plot.setup){
		plot_setup(parms, powfs, aper, recon);
	}
	global->setupdone=1;
	if(parms->save.setup&&chdir("..")){}
	/*
	  Before entering real simulation, make sure to delete all variables that
	  won't be used later on to save memory.
	*/
	free_powfs_unused(parms, powfs);
	//free_recon_unused(parms, recon);
	if(parms->evl.psfmean||parms->evl.psfhist){
	/*compute diffraction limited PSF. Save to output directory.*/
		dmat* iopdevl=dnew(aper->locs->nloc, 1);
		ccell* psf2s=0;
		locfft_psf(&psf2s, aper->locfft, iopdevl, parms->evl.psfsize, 0);
		const int nwvl=parms->evl.nwvl;
		dcell* evlpsfdl=dcellnew(nwvl, 1);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			cabs22d(&P(evlpsfdl,iwvl), 1, P(psf2s,iwvl), 1);
			P(evlpsfdl,iwvl)->keywords=evl_keywords(parms, aper, -1, iwvl, parms->evl.psfisim-1);
		}
		if(parms->plot.run||parms->plot.psf){
			plot_psf(psf2s, "PSFdl", 2, 0, parms->evl.wvl, parms->plot.psf!=2, parms->plot.psfmin);
		}
		ccellfree(psf2s);
		writebin(evlpsfdl, "evlpsfdl.fits");
		dcellfree(evlpsfdl);
		dfree(iopdevl);
	}
	if(parms->sim.skysim){
		save_skyc(powfs, recon, parms);
	}
#if USE_CUDA
	cuda_dedup=0;
#endif
	toc2("setup");
	print_mem("After setup");
}

/**
   Free all allocated memory in setup_* functions. So that we
   keep track of all the memory allocation.*/
void maos_reset(){
	if(!global) return;
	dbg("Deleting allocated data\n");
	parms_t* parms=(parms_t*)global->parms;
	free_recon(parms, global->recon);
	free_powfs(parms, global->powfs);
	free_aper(global->aper);
#if USE_CUDA
	if(use_cuda){
		gpu_cleanup();
	}
#endif
	free(global); global=NULL;
}
