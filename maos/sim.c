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
#include <search.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "common.h"
#include "recon.h"
#include "recon_utils.h"
#include "powfs.h"
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
extern int KEEP_MEM;

int sim_pipe[2]={0,0};
/**
   \file sim.h

   Contains main simulation blocks.
*/
/**
   Initialize the simulation runtime data struct.
   \callgraph
 */
sim_t* maos_iseed(int iseed){
	const parms_t* parms=global->parms;
	powfs_t* powfs=global->powfs;
	aper_t* aper=global->aper;
	recon_t* recon=global->recon;
	if(parms->fdlock&&parms->fdlock[iseed]<0){
		warning("Another MAOS is already running. Skip seed %ld\n",
			P(parms->sim.seeds, iseed));
		return 0;
	}
	sim_t* simu=init_simu(parms, powfs, aper, recon, iseed);
	if(iseed==0) simu->tk_s0=myclockd();
	simu->tk_si=myclockd();
	global->simu=simu;
	global->iseed=iseed;

	if(parms->atm.frozenflow){
		genatm(simu);/*Generating atmospheric screen(s) that frozen flows.*/
		if(parms->tomo.predict){
			if(recon->HXWtomo){
				setup_recon_HXW(simu->recon, simu->parms, simu->atm);
			}
			if(parms->tomo.precond==1){
				fdpcg_free(recon->fdpcg);
				recon->fdpcg=fdpcg_prepare(parms, recon, powfs, parms->tomo.predict?simu->atm:NULL);
			}
		}
	}
#if USE_CUDA
	if(parms->gpu.evl||parms->gpu.wfs){
		/*Put here for the initial transfer to avoid messing up timing due to transfering. */
		gpu_atm2gpu(simu->atm, simu->atmscale, parms, iseed, parms->sim.start);/*takes 0.4s for NFIRAOS. */
		if(parms->tomo.predict){
			gpu_update_recon_cn2(parms, recon);
		}
	}
#endif
	return simu;
}
/**
	Prepare before each isim
*/
void prepare_isim(sim_t *simu){
	int isim=simu->perfisim;
	const parms_t *parms=simu->parms;
	recon_t *recon=simu->recon;

	if(!parms->atm.frozenflow){
		//Do not put this one inside parallel single so that FFT can use parallel for
		genatm(simu);
		/*re-seed the atmosphere in case atm is loaded from shm/file */
		seed_rand(simu->atm_rand, lrand(simu->init_rand));
	}
#if USE_CUDA
	if(parms->gpu.evl||parms->gpu.wfs){
	/*may need to copy another part */
		gpu_atm2gpu(simu->atm, simu->atmscale, parms, simu->iseed, isim);
	}
#endif
	update_wfsflags(simu);
	sim_update_zoom(simu);
	sim_update_etf(simu);
	if(parms->sim.dmproj){
		/* temporarily disable FR.M so that Mfun is used.*/
		cell *FRM=recon->fit->FR.M; recon->fit->FR.M=NULL;
		muv_solve(&simu->dmproj, &recon->fit->FL, &recon->fit->FR, NULL);
		recon->fit->FR.M=FRM;/*set FR.M back*/
		save_dmproj(simu);
		if(!parms->fit.square){
			/* Embed DM commands to a square array for fast ray tracing */
			for(int idm=0; idm<parms->ndm; idm++){
				loc_embed(P(simu->dmprojsq, idm), P(recon->aloc, idm), P(simu->dmproj, idm));
			}
		}
#if USE_CUDA
		if(parms->gpu.evl||parms->gpu.wfs){
			gpu_dmproj2gpu(simu->dmprojsq);
		}
#endif
	}
	if(!simu->pause&&parms->sim.end>10+parms->sim.start&&parms->sim.dtrat_hi+isim<parms->sim.end&&!signal_caught){
		draw_single=1;//Only draw active frame.
	} else{
		draw_single=0;
	}
}
/**
 * @brief Post process each isim
 * 
 */
void postproc_isim(sim_t *simu){
	if(simu->tomo_update){//This part causes random CUDA error in Geforce.
		if(simu->tomo_update==1){//Only update cn2 regularization term
			setup_recon_update_cn2(simu->recon, simu->parms);
			//Already updates GPU in this routine
		} else{//Also update noise in reconstructor.
			setup_recon_control(simu->recon, simu->parms, simu->powfs);
			dcellzero(simu->opdr);//zero out warm restart data. [optional in CPU, essential in GPU.]
#if USE_CUDA
			if(simu->parms->gpu.tomo){
				gpu_update_recon_control(simu->parms, simu->recon);
			}
#endif
		}
		simu->tomo_update=0;
	}
}

/**
   Simulation for each time step.

   Callable from matlab.
   \callgraph
*/
void maos_isim(int isim){
	sim_t* simu=global->simu;
	const parms_t* parms=simu->parms;
	int simstart=parms->sim.start;
	tp_counter_t group={0};
	extern int NO_RECON, NO_WFS, NO_EVL;
	if(isim==simstart+1){//skip slow first step.
		simu->tk_i1=myclockd();
	}
	simu->tk_istart=myclockd();//step start time
	simu->wfsisim=isim;
	simu->perfisim=isim;
	simu->status->isim=isim;
	if(!parms->sim.closeloop){
		simu->reconisim=isim;
	} else{//work on gradients from last time step for parallelization.
		simu->reconisim=isim-1;
	}
	prepare_isim(simu);
	if(PARALLEL==1){
		/*
		  We do the big loop in parallel to make better use the
		  CPUs. Notice that the reconstructor is working on grad from
		  last time step so that there is no confliction in data access.

		  Launch perfevl_pre and wfsgrad_pre to allow gpu code to execute in advance.
		*/

		if(parms->gpu.evl&&!NO_EVL){
			//Queue tasks on GPU, no stream sync is done
			QUEUE_THREAD(&group, simu->perfevl_pre, 0);
		}
		if(parms->tomo.ahst_idealngs!=1&&parms->gpu.wfs&&!NO_WFS){
			//task for each wfs
			QUEUE_THREAD(&group, simu->wfsgrad_pre, 0);
		}
		if(!NO_RECON){
			//don't put this first. It has cpu overhead in computing gradol
			QUEUE(&group, reconstruct, simu, 1, 0);
		}
		if(!NO_EVL){
			if(parms->gpu.evl){
				//wait for GPU tasks to be queued before calling sync
				WAIT(&group, 0);
			}
			QUEUE(&group, perfevl, simu, 1, 0);
		}
		if(!NO_WFS){
			if(parms->tomo.ahst_idealngs==1||(parms->gpu.wfs&&!parms->gpu.evl)){
				/*when we want to apply idealngs correction, wfsgrad need to wait for perfevl. */
				WAIT(&group, 0);
			}
			QUEUE(&group, wfsgrad, simu, 1, 0);
		}
		if(!NO_RECON){
			//wait for all tasks to finish before modifying grad and dmreal
			WAIT(&group, 0);
			shift_grad(simu);/*before filter() */
			filter_dm(simu);/*updates dmreal, so has to be after prefevl/wfsgrad is done. */
		}
		WAIT(&group, 0);
	} else{/*PARALLEL=0: do the big loop in serial mode. */
		if(parms->sim.closeloop){
			if(!NO_EVL) perfevl(simu);/*before wfsgrad so we can apply ideal NGS modes */
			if(!NO_WFS)	wfsgrad(simu);/*output grads to gradcl, gradol */
			if(!NO_RECON){
				reconstruct(simu);/*uses grads from gradlast cl, gradlast ol. */
				shift_grad(simu);
				filter_dm(simu);
			}
		} else{/*in OL mode,  */
			if(!NO_WFS)	wfsgrad(simu);
			if(!NO_RECON){
				shift_grad(simu);
				reconstruct(simu);
				filter_dm(simu);
			}
			if(!NO_EVL) perfevl(simu);
		}
	}
	simu->tk_iend=myclockd();
	print_progress(simu);
}
static void*perfevl_loop(sim_t *simu){
	int simstart=simu->parms->sim.start;
	int simend=simu->parms->sim.end;
	for(int isim=simstart; isim<simend&&!signal_caught; isim++){
		if(isim==simstart+1){//skip slow first step.
			simu->tk_i1=myclockd();
		}
#if USE_CUDA
		if(simu->parms->gpu.evl||simu->parms->gpu.wfs){
		/*may need to copy another part */
			gpu_atm2gpu(simu->atm, simu->atmscale, simu->parms, simu->iseed, isim);
		}
#endif
		simu->tk_istart=myclockd();//step start time
		simu->perfisim=isim;//info("call perfevl(%d)\n", isim);
		simu->status->isim=isim;
		perfevl(simu);
		simu->tk_iend=myclockd();
		print_progress(simu);
	}
	return NULL;
}
static void *wfsgrad_loop(sim_t *simu){
	int simstart=simu->parms->sim.start;
	int simend=simu->parms->sim.end;
	for(int isim=simstart; isim<simend&&!signal_caught; isim++){
		simu->wfsisim=isim;//info("call wfsgrad(%d)\n", isim);
		prepare_isim(simu);
		wfsgrad(simu);
		shift_grad(simu);
	}
	return NULL;
}
static void *reconstruct_loop(sim_t *simu){
	int simstart=simu->parms->sim.start;
	int simend=simu->parms->sim.end;
	for(int isim=simstart; isim<simend&&!signal_caught; isim++){
		simu->reconisim=isim-1;//info("call reconstruct(%d)\n", isim);
		reconstruct(simu);
		filter_dm(simu);
		postproc_isim(simu);
	}
	return NULL;
}
/**
   Closed loop simulation main loop.

   It calls init_simu() to initialize the simulation struct, and then calls
   maos_isim() for each simulation time step. Arranged this way so that
   maos_isim() can be called from matlab.
   \callgraph
*/
void maos_sim(){
	const parms_t* parms=global->parms;
	recon_t* recon=global->recon;
	int simend=parms->sim.end;
	int simstart=parms->sim.start;
	if(sim_pipe[0]==0){
		if(pipe(sim_pipe)){
			sim_pipe[0]=-1;//indicate error
			sim_pipe[1]=-1;//indicate error
			info_errno("pipe");
		}
	}
	dbg("PARALLEL=%d, NTHREAD=%d\n", PARALLEL, NTHREAD);
	long rescount=0;
	dmat* restot=dnew(parms->evl.nmod, 1);
	for(int iseed=0; iseed<parms->sim.nseed&&!signal_caught; iseed++){
		sim_t* simu=NULL;
		while(!(simu=maos_iseed(iseed))){
			iseed++;
		}
		if(recon&&recon->cn2est){//temporary. Should put runtime data in simu.
			cn2est_reset(recon->cn2est);
		}

#ifdef HAVE_NUMA_H
		numa_set_localalloc();
#endif
		if(PARALLEL==2){//event driven synchronization
			/*
				There are two pairs of synchronization.
				First pair: regarding update and consumption of dmreal, between filter and perfevl/wfsgrad.
				Second pair: regarding update and consumption of wfsgrad, between shift_grad and reconstruct.
				simulation index are used to indicate whether correct info is available.
				counters are used to indicate consumption and reset when there is an update.
			*/
			tp_counter_t group={0};
			QUEUE(&group, perfevl_loop, simu, 1, 0);
			QUEUE(&group, wfsgrad_loop, simu, 1, 0);
			QUEUE(&group, reconstruct_loop, simu, 1, 0);
			WAIT(&group, 0);
		} else{
			for(int isim=simstart; isim<simend&&!signal_caught; isim++){
				maos_isim(isim);
				if(simu->pause>0&&isim%simu->pause==0){
					simu->pause=mypause(0, sim_pipe[0]);
				}
			}/*isim */
		}
		long isim0=0;
		if(!parms->sim.evlol){//determine starting point
			double sumc=0;
			for(isim0=simend-1; isim0>1; isim0--){
				sumc+=P(simu->cle, 0, isim0);
				if((simend-isim0)>25 && P(simu->cle, 0, isim0-1)*(simend-isim0)>sumc*1.2){
					break;
				}
			}
		}
		/*Compute average performance*/
		const long nsim=simu->perfisim+1;
		//const long isim0=parms->sim.closeloop?MIN(MAX(20, nsim/5), nsim/2):0;
		const dmat* cle=parms->sim.evlol?simu->ole:simu->cle;
		for(long i=isim0; i<nsim; i++){
			for(long imod=0; imod<parms->evl.nmod; imod++){
				P(restot, imod)+=P(cle, imod, i);
			}
		}
		rescount+=(nsim-isim0);

		simu->status->clerrhi=sqrt(P(restot, 2)/rescount)*1e9;
		simu->status->clerrlo=sqrt(P(restot, 1)/rescount)*1e9;
		simu->status->tot=simu->status->mean;
#if defined(__linux__) || defined(__APPLE__)
		scheduler_report(simu->status);
#endif
		free_simu(simu);
		remove_lock(parms->fdlock, parms->fnlock, P(parms->sim.seeds), PN(parms->sim.seeds), simu->iseed, signal_caught==0);
		global->simu=NULL;
	}/*seed */

	info3("Mean: ");
	for(long imod=0; imod<parms->evl.nmod; imod++){
		info3("%.2f ", sqrt(P(restot, imod)/rescount)*1e9);
	}
	info3(" nm (%ld points).\n", rescount);
	dfree(restot);
}
