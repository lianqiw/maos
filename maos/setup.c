/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/**
  \file maos.h
  Sets up maos simulation.

  The main() is separated  into main.c so that maos.c can be a part of libaos.la which is callable by
	MATLAB.
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
	//dmproj/frozen flow does not work with PRALLEL=1 implementation
	if(PARALLEL==2&&(parms->sim.dmproj||!parms->atm.frozenflow)){
		warning("dmproj/frozen flow does not work with PRALLEL=2 implementation. Set to 1.\n");
		PARALLEL=1;
	}

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
	#if CUDA_VERSION>10000
	extern void (*dsvd_large)(dmat**U_, dmat**S_, dmat**Vt_, dmat*A_);
	warning("set dsvd_large to gpu_dsvd\n");
	dsvd_large=gpu_dsvd;
	if(0){//test passed
		//gpu_dsvd is faster than dsvd for matrix larger than 500x500 (test is on cassiopeia)
		dmat *A=dnew(500, 500);
		rand_t rstat; seed_rand(&rstat, 1);
		drandn(A,1,&rstat);
		daddI(A,1);
		dmat *U=0, *Vt=0, *S=0;
		dmat *U2=0, *Vt2=0, *S2=0;
		tic;
		dsvd(&U, &S, &Vt, A);
		toc("dsvd");tic;
		gpu_dsvd(&U2, &S2, &Vt2, A);
		toc("gpu_dsvd");
		real diffu=ddiff(U,U2);
		real diffs=ddiff(S, S2);
		real diffv=ddiff(Vt, Vt2);
		info("dsvd and gpu_svd diff are %g %g %g\n", diffu, diffs, diffv);
		writebin(U, "U");
		writebin(S, "S");
		writebin(Vt, "Vt");
		writebin(U2, "U2");
		writebin(S2, "S2");
		writebin(Vt2, "Vt2");
		writebin(A, "A");
		exit(0);
	}
	#endif
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
		global->powfs=powfs=setup_powfs_init(parms, aper);
		print_mem("After setup_powfs");
		/*Setup DM fitting parameters so we can flatten the DM in setup_surf.c */
		global->recon=recon=setup_recon_prep(parms, aper, powfs);
		/*setting up M1/M2/M3, Instrument, Lenslet surface OPD. DM Calibration, WFS bias.*/
		setup_surf(parms, aper, powfs, recon);
		/*set up physical optics powfs data. It needs dmncpa and wfsadd from setup_surf()*/
		setup_powfs_phy(parms, powfs);
		/*setup gradient noise during simulation. */
		setup_powfs_neasim(parms, powfs);
		/*calibrate gradient offset*/
		setup_powfs_calib(parms, powfs);
#if USE_CUDA
		if(parms->gpu.wfs&&powfs&&has_pywfs){
			gpu_wfsgrad_init(parms, powfs);//needed by pywfs_mkg
		}
#endif
		setup_recon_prep2(recon, parms, aper, powfs);
		//Don't put this inside parallel, otherwise svd will run single threaded.
		setup_recon_control(recon, parms, powfs);
		if(parms->recon.split){
			/*split tomography */
			ngsmod_setup(parms, recon);
			if(!parms->sim.idealfit&&parms->recon.split==2&&parms->recon.alg==0){/*Need to be after fit */
				setup_recon_mvst(recon, parms);
			}
		}
		setup_recon_dither_dm(recon, powfs, parms);//depends on saneai
		//setup_recon_sodium_fit(recon, parms);//to restrict modes in sodium fitting correction
		if(parms->recon.alg==0||parms->sim.dmproj){
			setup_recon_fit(recon, parms);
		}
		if(recon->actcpl && !recon->actextrap){
			recon->actextrap=act_extrap(recon->aloc, recon->actcpl, parms->lsr.actthres, 1);
			if(parms->save.setup){
				writebin(recon->actextrap, "actextrap");
			}
		}
		if(parms->sim.wfsalias==2||parms->sim.idealwfs==2){
			setup_powfs_fit(powfs, recon, parms);
		}
		if(parms->recon.alg==0&&parms->nmoao){
			setup_recon_moao(recon, parms);
		}
		print_mem("After setup_recon");
		if(parms->dbg.wfslinearity!=-1){
			int iwfs=parms->dbg.wfslinearity;
			assert(iwfs>-1||iwfs<parms->nwfs);
			info2("Studying wfslineariy for WFS %d\n", iwfs);
			wfslinearity(parms, powfs, iwfs);
			((parms_t*)parms)->sim.end=parms->sim.start;//indicate no simulation
		}
		{
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
	//gpu_wfsgrad_init(parms, powfs); //moved to above
		gpu_wfssurf2gpu(parms, powfs);
	}
#endif

	if(!parms->sim.evlol){
		setup_recon_post(recon, parms, aper);
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
		locfft_psf(&psf2s, aper->embed, iopdevl, parms->evl.psfsize, 0);
		const int nwvl=parms->evl.nwvl;
		dcell* evlpsfdl=dcellnew(nwvl, 1);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			cabs22d(&P(evlpsfdl,iwvl), 1, P(psf2s,iwvl), 1);
			P(evlpsfdl,iwvl)->keywords=evl_keywords(parms, aper, -1, iwvl, parms->evl.psfisim-1);
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
	toc2("Presimulation");
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
