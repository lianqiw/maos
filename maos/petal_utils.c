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

#include "common.h"
#include "save.h"
#include "plot_utils.h"
/**
 * Petaling mode control
*/
void petal_setup_wrap(recon_t *recon, const parms_t *parms, const powfs_t *powfs){
	if(!parms->recon.petal) return;
	int idm=parms->idmground;
	recon->apetal=petal_mkh_loc(P(recon->aloc,idm), 6, parms->aper.rot);
	if(parms->save.setup){
		writebin(recon->apetal, "apetal");
	}
	recon->petal=mycalloc(2, petal_t*);
	for(int ir=0; ir<2; ir++){
		int ipowfs=0;
		if(ir==0&&parms->ittfpowfs!=-1){
			ipowfs=parms->ittfpowfs;
		} else if(ir==1&&parms->ittpowfs!=-1){
			ipowfs=parms->ittpowfs;
		} else{
			continue;
		}
		real nembed=2;
		real dsa=powfs[ipowfs].pts->dsa;
		real dtheta=parms->powfs[ipowfs].wvlmean/(nembed*dsa);
		real pdtheta=parms->powfs[ipowfs].pixtheta/dtheta;
		dbg("powfs[%d].pdtheta=%g\n", ipowfs, pdtheta);
		if(fabs(pdtheta-1)>0.01){
			warning("TODO: pdtheta!=1 requries resampling PSFs\n");
		}
		//only withtt only for t/t oiwfs unless petaltt>1. 
		//enable it for TTF OIWFS sometimes results in a clocking gradient pattern.
		int withtt=(parms->powfs[ipowfs].order==1||parms->recon.petaltt>1)?parms->recon.petaltt:0;
		recon->petal[ir]=petal_setup(powfs[ipowfs].pts->loc, powfs[ipowfs].loc->dx, P(powfs[ipowfs].amp, 0),
			pdtheta, parms->powfs[ipowfs].pixblur, parms->aper.rot, parms->recon.petalnpsf, withtt);
		if(parms->save.setup){
			petal_save(recon->petal[ir], "petal_%d", ir);
		}
	}
}
static unsigned int petal_bg_status=0;//if set, indicates that the thread wfsgrad_peta_bg is finished
/**
 * @brief Background profess to run the slow reconstruction process
 * 
 * @param simu 
 * @return void* 
 */
static void *wfsgrad_petal_thread(sim_t *simu){
	const int nrep=8;
	for(int ir=0; ir<2; ir++){
		dmat *phib=ir==1?P(simu->petal_m, 0):NULL;//2nd step use first step as input.
		petal_solve(NULL, &P(simu->petal_m, ir), simu->recon->petal[ir], P(simu->petal_i0, ir, 1), phib, nrep);
	}
	petal_bg_status=1;
	return NULL;
}
/**
 * Reconstruct petal modes. The function is now split into three parts:
 * Part 1: accumuate i0 at every time step and set do_petal if a reconstruction is due.
 * Part 2: launch a new thread to do the petal reconstruction.
 * Part 3: join the thread and do post-processing to update the reference vectors.
*/
void petal_recon(sim_t *simu){
	const parms_t *parms=simu->parms;

	if(!simu->petal_i0){
		simu->petal_i0=dccellnew(2,2);//first column for accumulation, second column for storage
	}
	if(!simu->petal_m){
		simu->petal_m=dcellnew(3,1);
	}
	int isim=simu->wfsisim;
	int do_petal=0;
	real rad2m=0;
	for(int ir=0; ir<2; ir++){
		int ipowfs=0;
		if(ir==0&&parms->ittfpowfs!=-1){
			ipowfs=parms->ittfpowfs;
		} else if(ir==1&&parms->ittpowfs!=-1){
			ipowfs=parms->ittpowfs;
		}else{
			continue;
		}
		rad2m=parms->powfs[ipowfs].wvlmean/TWOPI;//radian to m
		if(isim<parms->powfs[ipowfs].phystep || isim<parms->recon.petalstep) continue;
		int iwfs=P(parms->powfs[ipowfs].wfs, 0);
		dcelladd(&P(simu->petal_i0, ir, 0), 1, P(simu->ints, iwfs), 1);
		if((isim+1-parms->powfs[ipowfs].step)%parms->recon.petaldtrat==0){
			//normalization is not necessary
			if(parms->save.recon){
				writebin(P(simu->petal_i0, ir, 0), "petal_i0_%d_%d", ir, isim);
			}
			/*if(parms->plot.run){
				draw("Petal", (plot_opts){ .image=P(P(simu->petal_i0, ir), 0) }, "PSF of first subaperture", "x", "y", "powfs %d", ipowfs);
			}*/

			dcelladd(&P(simu->petal_i0, ir, 1), 0, P(simu->petal_i0, ir, 0), 1);
			dcellzero(P(simu->petal_i0, ir, 0));
			do_petal=1;
		}
	}
	static pthread_t thread=0;
	static int petal_isim=0;
	if(thread && (petal_bg_status>0 || do_petal || isim+1==parms->sim.end)){//join previous thread and finish up
		void *ans;
		pthread_join(thread, &ans); thread=0; petal_bg_status=0;
		//post-processing petalling results
		dadds(P(simu->petal_m, 1), -dmean(P(simu->petal_m, 1)));//remove piston
		dshow(P(simu->petal_m, 1), "Step%6d (from %d): petal output:", isim, petal_isim);
		dscale(P(simu->petal_m, 1), rad2m);//convert to m
		int idm=parms->idmground;
		dcellzero(simu->dmtmp);//this cannot be done in the thread as dmtmp is used by filter().
		real gain=parms->recon.petaldtrat==1?0.5:1;//gain of 1 can be used if peltadtrat>1
		dspmm(&P(simu->dmtmp, idm), simu->recon->apetal, P(simu->petal_m, 1), "nn", 1);
		if(1){
			warning_once("Remove p/t/t from dm petal offset.\n");
			loc_remove_ptt(P(simu->dmtmp, idm), P(simu->recon->aloc, idm), NULL, NULL, 0);
		}
		for(int jwfs=0; jwfs<parms->nwfs; jwfs++){
			int jpowfs=parms->wfs[jwfs].powfs;
			if(parms->powfs[jpowfs].lo) continue;
			dcellmm(&P(simu->gradoff, jwfs), P(simu->recon->GA, jwfs, idm), P(simu->dmtmp, idm), "nn", -gain);
		}
		if(parms->plot.run&&petal_isim%parms->plot.run==0){
			int draw_single_save=draw_single; draw_single=0;
			drawopd("DM", P(simu->recon->aloc, idm), P(simu->dmtmp, idm), parms->plot.opdmax, "DM Petal Error Signal (Hi)", "x (m)", "y (m)", "Petal %d", idm);
			plot_gradoff(simu, -1);
			draw_single=draw_single_save;
		}
		servo_add(simu->dmint, simu->dmtmp, gain);
	}
	if(do_petal){//launch a thread to do the task
		petal_isim=isim;
		pthread_create(&thread, NULL, (thread_fun)wfsgrad_petal_thread, simu);
	}
}
