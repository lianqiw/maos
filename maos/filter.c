/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "sim.h"
#include "ahst.h"
#include "save.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/*
   Collection of functions for servo filtering of DM commands.

   Applies hysterisis. Input dmcmd is command to the DM and output dmreal is the
   actual position the DM goes to.
*/
/**
   Add low order NGS modes to DM actuator commands for AHST and MVST
 */
void addlow2dm(dcell** dmval, const sim_t* simu, const dcell* low_val, real gain){
	switch(simu->parms->recon.split){
	case 0:
		break;/*nothing to do. */
	case 1:
		dcellmm(dmval, simu->recon->ngsmod->Modes, low_val, "nn", gain);
		break;
	case 2:
		dcellmm(dmval, simu->recon->MVModes, low_val, "nn", gain);
		break;
	default:
		error("Not implemented\n");
	}
}
static inline int limit_diff(real* x1, real* x2, real thres, long stuck1, long stuck2){
	real diff=*x2-*x1;
	if(fabs(diff)>thres){
		real ratio=signbit(diff)?-.49999:.49999;
		if(stuck1){
			*x2=*x1+thres*ratio*2;
		} else if(stuck2){
			*x1=*x2-thres*ratio*2;
		} else{
			real mean=0.5*(*x1+*x2);
			*x1=mean-thres*ratio;
			*x2=mean+thres*ratio;
		}
		return 1;
	}
	return 0;
}
/**
   Send LPF TT to TTM. Use DMTT, DMPTT to take into account possible stuck actuators.
*/
static inline void ttsplit_do(recon_t* recon, dcell* dmcmd, dmat* ttm, real lp){
#if 1
	int ndm=NX(dmcmd);
	real totaltt[2]={0,0};
	real tt1[2];
	for(int idm=0; idm<ndm; idm++){
		tt1[0]=tt1[1]=0;
		dmulvec(tt1, P(recon->DMPTT,idm), P(P(dmcmd,idm)), 1);
		dmulvec(P(P(dmcmd,idm)), P(recon->DMTT,idm), tt1, -1);
		for(int i=0; i<2; i++){
			totaltt[i]+=tt1[i];
		}
	}
	P(ttm,0)=P(ttm,0)*(1-lp)+lp*totaltt[0];
	P(ttm,1)=P(ttm,1)*(1-lp)+lp*totaltt[1];
	totaltt[0]-=P(ttm,0);
	totaltt[1]-=P(ttm,1);
	//Put HPF'ed to ground DM.
	dmulvec(P(P(dmcmd,0)), P(recon->DMTT,0), totaltt, 1);
#else
	//Only touch ground DM
	real tt1[2]={0,0};
	dmulvec(tt1, P(recon->DMPTT,0), P(P(dmcmd,0)), 1);
	P(ttm,0)=P(ttm,0)*(1-lp)+lp*tt1[0];
	P(ttm,1)=P(ttm,1)*(1-lp)+lp*tt1[1];
	dmulvec(P(P(dmcmd,0)), P(recon->DMTT,0), P(ttm), -1);
#endif
}

static inline void clipdm(sim_t* simu, dcell* dmcmd){
	const parms_t* parms=simu->parms;
	if(!dmcmd) return;
	/*
	  clip integrator. This both limits the output and
	  feeds back the clip since we are acting on the integrator directly.
	*/
	if(!parms->sim.dmclip) return;
	for(int idm=0; idm<parms->ndm; idm++){
		const int nact=P(dmcmd,idm)->nx;
		if(NX(parms->dm[idm].stroke)==1){
			if(NY(parms->dm[idm].stroke)!=1){
				error("dm.stroke is in wrong format\n");
			}
			static int nclip0=0;
			if(simu->reconisim<10) nclip0=0;
			int nclip=dclip(P(dmcmd,idm),
				-P(parms->dm[idm].stroke,0),
				P(parms->dm[idm].stroke,0));
			if(nclip>nclip0){
				nclip0=nclip;
				info2("step %d DM %d: %d actuators clipped\n", simu->reconisim, idm, nclip);
			}
		} else if(NX(parms->dm[idm].stroke)==nact){
			if(NY(parms->dm[idm].stroke)!=2){
				error("dm.stroke is in wrong format\n");
			}

			real* pcmd=P(P(dmcmd,idm));
			real* plow=P(parms->dm[idm].stroke);
			real* phigh=plow+nact;
			for(int iact=0; iact<nact; iact++){
				if(pcmd[iact]<plow[iact]){
					pcmd[iact]=plow[iact];
				} else if(pcmd[iact]>phigh[iact]){
					pcmd[iact]=phigh[iact];
				}
			}
		} else{
			error("Invalid format\n");
		}
	}
}
static inline void clipdm_dead(const sim_t* simu, dcell* dmcmd){
	if(!dmcmd) return;
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	if(!recon->actstuck) return;
	for(int idm=0; idm<parms->ndm; idm++){
		if(!P(recon->actstuck,idm)) continue;
		const int nact=P(recon->aloc,idm)->nloc;
		for(int iact=0; iact<nact; iact++){
			real val=P(P(recon->actstuck,idm),iact)*1e-9;
			if(val){
				P(P(dmcmd,idm),iact)=val;
			}
		}
	}
}

static inline void clipdm_ia(const sim_t* simu, dcell* dmcmd){
	if(!dmcmd) return;
	const parms_t* parms=simu->parms;
	if(!parms->sim.dmclipia) return;
	recon_t* recon=simu->recon;
	/*Clip interactuator stroke*/
	for(int idm=0; idm<parms->ndm; idm++){
	/* Embed DM commands to a square array (borrow dmrealsq) */
		real iastroke, iastroked;
		int nx=P(recon->anx,idm);
		dmat* dmr;
		dmat* dm;
		if(parms->dm[idm].iastrokescale){ //convert dm to voltage
			dm=dinterp1(P(parms->dm[idm].iastrokescale,0), 0, P(dmcmd,idm), NAN);
			iastroke=parms->dm[idm].iastroke;//voltage.
		} else{
			dm=P(dmcmd,idm);
			iastroke=parms->dm[idm].iastroke*2;//surface to opd
		}
		iastroked=iastroke*1.414;//for diagonal separation.
		if(!parms->fit.square){
			loc_embed(P(simu->dmrealsq,idm), P(recon->aloc,idm), P(dm));
			dmr=(dmat*)P(simu->dmrealsq,idm);
		} else{
			dmr=dm;
		}
		lcell* actstuck=recon->actstuck;
		//offset point by 1 because map is 1 based index.
		long* stuck=actstuck&&P(actstuck,idm)?(P(P(actstuck,idm))-1):0;
		int count=0, trials=0;
		do{
			count=0;
			dmat* map=(dmat*)P(recon->amap,idm)/*PDMAT*/;
			for(int iy=0; iy<P(recon->any,idm)-1; iy++){
				for(int ix=0; ix<nx-1; ix++){
					int iact1=P(map, ix, iy);
					if(iact1>0){
						int iact2=P(map, ix, iy+1);
						real* dmri=&P(dmr, ix, iy);
						if(iact2>0){
							count+=limit_diff(dmri, &P(dmr, ix, iy+1), iastroke,
								stuck?stuck[iact1]:0, stuck?stuck[iact2]:0);
						}
						iact2=P(map, ix+1, iy);
						if(iact2>0){
							count+=limit_diff(dmri, &P(dmr, ix+1, iy), iastroke,
								stuck?stuck[iact1]:0, stuck?stuck[iact2]:0);
						}
						iact2=P(map, ix+1, iy+1);
						if(iact2>0){
							count+=limit_diff(dmri, &P(dmr, ix+1, iy+1), iastroked,
								stuck?stuck[iact1]:0, stuck?stuck[iact2]:0);
						}
					}
				}
			}
			trials++;
			if(trials==1&&count>0){
				info2("Step %d, DM %d: %d actuators over ia limit. ", simu->reconisim, idm, count);
			}
		} while(count>0&&trials<100);
		if(trials>1){
			info2("trials=%d: %s\n", trials, count?"failed.":"success.");
		}
		if(!parms->fit.square){//copy data back
			loc_extract(P(simu->dmreal,idm), P(recon->aloc,idm), P(simu->dmrealsq,idm));
		}
		if(parms->dm[idm].iastrokescale){//convert back to opd
			dmat* dm2=dinterp1(P(parms->dm[idm].iastrokescale,1), 0, dm, NAN);
			dcp(&P(dmcmd,idm), dm2);
			dfree(dm); dfree(dm2);
		}
	}
}

/**
   Update DM command for next cycle using info from last cycle (two cycle delay)
in closed loop mode */
static void filter_cl(sim_t* simu){
	/*
	  2009-11-02: Moved to the end of isim loop to update
	  for next step.  only need to cache a single dmerrlast
	  now.

	  2009-12-23: Updated low fs to do lead filter/type II

	  2010-01-07: Create an option to merge the last
	  integrator in the hi/lo loop to simulation the actual
	  block diagram. removed dmreal_hi, Mreal_lo;

	  2010-01-08: Changed the filtering scheme by computing
	  dm command for next cycle instead of maintianing last
	  step error information.

	  2010-01-13: Implemented aphi.
	  a(n)=a(n-1)+ep*e(n-2) or
	  a(n)=0.5*(a(n-1)+a(n-2))+ep*e(n-2);
	*/
	const parms_t* parms=simu->parms;
	recon_t* recon=simu->recon;
	assert(parms->sim.closeloop);
	/*copy dm computed in last cycle. This is used in next cycle (already after perfevl) */
	const sim_cfg_t* simcfg=&(parms->sim);
	const int isim=simu->reconisim;
	//dbg("reconisim=%d\n", isim);
	{/*Auto adjusting ephi for testing different ephi*/
		static int ephi_is_auto=0;
		if(P(simcfg->ephi,0)<0){
			ephi_is_auto=1;
			P(simcfg->ephi,0)=0.5;
		}
		if(ephi_is_auto){
			if((isim*10)<parms->sim.end){//initial steps
				P(simcfg->ephi,0)=0.5;
			} else if((isim*10)%parms->sim.end==0){
				P(simcfg->ephi,0)=(real)isim/(real)parms->sim.end;
				dbg("ephi is set to %.1f at step %d\n", P(simcfg->ephi,0), isim);
			}
		}
	}

	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	//Record dmpsol for this time step for each powfs before updating it (z^-1).
	//Do not reference the data, even for dtrat==1
		if(!parms->powfs[ipowfs].psol||!parms->powfs[ipowfs].dtrat) continue;
		real alpha=(isim%parms->powfs[ipowfs].dtrat==0)?0:1;
		dcelladd(&P(simu->wfspsol,ipowfs), alpha, simu->dmpsol, 1./parms->powfs[ipowfs].dtrat);
	}
	/*Do the servo filtering. First simulate a drop frame*/
	int drop=0;
	if(simu->dmerr&&parms->sim.dtrat_skip){
		if(parms->sim.dtrat_skip>0){
			if((isim+1)%parms->sim.dtrat_skip==0){//evenly
				drop=1;
			}
		} else if(parms->sim.dtrat_skip<0){//use random draws
			real tmp=randu(simu->misc_rand);
			if(tmp*(-parms->sim.dtrat_skip)<1.){
				drop=1;
			}
		}
	}
	dcell* dmerr=0;
	if(drop){
		warning("Drop a frame at step %d\n", isim);
	} else if(simu->dmerr){
		dmerr=simu->dmerr;
	}
	//always run servo_filter even if dmerr is NULL.
	int hiout=servo_filter(simu->dmint, dmerr);

	if(parms->recon.split){
	/*Low order in split tomography only. fused integrator*/
		if(servo_filter(simu->Mint_lo, simu->Merr_lo)&&parms->sim.fuseint){
			/*accumulate to the main integrator. Use mpreint to properly account
			 * for type II controler. Gain is already applied.*/
			dcellzero(simu->dmtmp);
			addlow2dm(&simu->dmtmp, simu, simu->Mint_lo->mpreint, 1);
			servo_add(simu->dmint, simu->dmtmp, 1);
		}
	}
	/*The following are moved from the beginning to the end because the
	  gradients are now from last step.*/
	servo_output(simu->dmint, &simu->dmtmp);
	if(parms->recon.split&&!parms->sim.fuseint){
		dcell* Mtmp=0;
		servo_output(simu->Mint_lo, &Mtmp);
		addlow2dm(&simu->dmtmp, simu, Mtmp, 1);
		dcellfree(Mtmp);
	}

	if(parms->recon.modal){
		//convert DM command from modal to zonal space
		for(int idm=0; idm<NX(simu->dmcmd); idm++){
			dmm(&P(simu->dmcmd,idm), 0, P(simu->recon->amod,idm), P(simu->dmtmp,idm), "nn", 1);
		}
	} else if(simu->recon->actinterp&&!parms->recon.psol){
		//Extrapolate to edge actuators
		dcellzero(simu->dmcmd);
		dspcellmm(&simu->dmcmd, simu->recon->actinterp, simu->dmtmp, "nn", 1);
	} else{
		dcellcp(&simu->dmcmd, simu->dmtmp);
	}
	//dmpsol should not contain DM offset vector, which may not be recovered by reconstruction.
	dcellcp(&simu->dmpsol, simu->dmcmd);
	//The DM commands are always on zonal modes from this moment

	if(simu->ttmreal){
		ttsplit_do(recon, simu->dmcmd, simu->ttmreal, parms->sim.lpttm);
	}
	if(parms->sim.focus2tel&&hiout){//offloading DM focus mode to telescope.
		dcellcp(&simu->telfocusreal, simu->telfocusint);
		dcellmm(&simu->telfocusint, recon->RFdm, simu->dmcmd, "nn", parms->sim.epfocus2tel);
	}
	if(recon->dither_m){
	//Change phase in calc_dither_amp if phase of dithering is changed
	//this is for step isim+1
		real anglei=((isim+1)/recon->dither_dtrat)*(2*M_PI/recon->dither_npoint);
		dcelladd(&simu->dmcmd, 1, recon->dither_m, sin(anglei));
	}

	if(!parms->dbg.ncpa_preload&&recon->dm_ncpa){
		info_once("Add NCPA after integrator\n");
		dcelladd(&simu->dmcmd, 1, recon->dm_ncpa, 1);
	}

	if(recon->actstuck&&!parms->recon.modal&&parms->dbg.recon_stuck){
	//zero stuck actuators so that gradpsol still contains gradients caused
	//by stuck actuators, so that MVR can smooth it out.
		act_stuck_cmd(recon->aloc, simu->dmpsol, recon->actstuck);
	}
	if(parms->dbg.dmoff){
		info_once("Add injected DM offset vector\n");
		int icol=(isim+1)%NY(parms->dbg.dmoff);
		for(int idm=0; idm<parms->ndm; idm++){
			dadd(&P(simu->dmcmd,idm), 1, P(parms->dbg.dmoff, idm, icol), 1);
		}
	}
	//Need to clip
	if(parms->sim.dmclip||parms->sim.dmclipia||recon->actstuck){
		int feedback=(parms->sim.dmclip||parms->sim.dmclipia)||(recon->actstuck&&parms->dbg.recon_stuck);
		if(parms->sim.dmclip) clipdm(simu, simu->dmcmd);
		if(parms->sim.dmclipia) clipdm_ia(simu, simu->dmcmd);
		if(recon->actstuck&&parms->dbg.recon_stuck) clipdm_dead(simu, simu->dmcmd);
		if(feedback){
			dcelladd(&simu->dmtmp, 1, simu->dmcmd, -1); //find what is clipped
			servo_add(simu->dmint, simu->dmtmp, -1);//remove from integrator (anti wind up)
			dcelladd(&simu->dmpsol, 1, simu->dmtmp, -1);//remove from dmpsol.
		}
		if(recon->actstuck&&!parms->dbg.recon_stuck) clipdm_dead(simu, simu->dmcmd);
	}
	/*This is after the integrator output and clipping*/
	if(simu->dmhist){
		for(int idm=0; idm<parms->ndm; idm++){
			if(P(simu->dmhist,idm)){
				dhistfill(&P(simu->dmhist,idm), P(simu->dmcmd,idm), 0,
					parms->dm[idm].histbin, parms->dm[idm].histn);
			}
		}
	}

	if(recon->moao&&!parms->gpu.moao){
		if(simu->dm_wfs){
			const int nwfs=parms->nwfs;
			for(int iwfs=0; iwfs<nwfs; iwfs++){
				int ipowfs=parms->wfs[iwfs].powfs;
				int imoao=parms->powfs[ipowfs].moao;
				if(imoao<0) continue;
				real g=parms->moao[imoao].gdm;
				dadd(&P(simu->dm_wfs, iwfs, 0), 1-g, P(simu->dm_wfs, iwfs, 1) , g);
			}
		}
		if(simu->dm_evl){
			const int nevl=parms->evl.nevl;
			int imoao=parms->evl.moao;
			real g=parms->moao[imoao].gdm;
			for(int ievl=0; ievl<nevl; ievl++){
				dadd(&P(simu->dm_evl, ievl, 0), 1-g, P(simu->dm_evl, ievl, 1), g);
			}
		}
	}
	if(PARALLEL==2){
	//Wait until all clients has consumed the last dmreal
		if(simu->dmreal_isim!=-1){
			int count=parms->evl.nevl;
			for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
				if(simu->reconisim+1>=parms->powfs[ipowfs].step){
					count+=parms->powfs[ipowfs].nwfs;
				}
			}
			while(simu->dmreal_count<count){
			//dbg("waiting: dmreal_count is %d need %d\n", simu->dmreal_count, parms->evl.nevl + parms->nwfs);
				pthread_cond_wait(&simu->dmreal_condw, &simu->dmreal_mutex);
				pthread_mutex_unlock(&simu->dmreal_mutex);
			}
			if(simu->dmreal_count>parms->evl.nevl+parms->nwfs){
				warning("error: dmreal_count is %d need %d\n", simu->dmreal_count, parms->evl.nevl+parms->nwfs);
			}
		}
		//dbg("ready: dmreal_count is %d\n", simu->dmreal_count);
	}
	/*hysteresis. */
	if(simu->hyst){
		hyst_dcell(simu->hyst, simu->dmreal, simu->dmcmd);
	} else{
		dcellcp(&simu->dmreal, simu->dmcmd);
	}
}
/**
   Servo filter for FSM commands.
 */
void filter_fsm(sim_t* simu){
	const parms_t* parms=simu->parms;
	if(simu->fsmint){
	/*fsmerr is from gradients from this time step. so copy before update for correct delay*/
		if(parms->sim.f0fsm>0){//Apply SHO filter
			dcell* ftmp=0;
			servo_output(simu->fsmint, &ftmp);
			for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
				if(NE(simu->fsmreal, iwfs)){
					real* pin=P(P(ftmp,iwfs));
					P(P(simu->fsmreal,iwfs),0)=sho_step(simu->fsmsho[iwfs], pin[0], parms->sim.dt);
					P(P(simu->fsmreal,iwfs),1)=sho_step(simu->fsmsho[iwfs+parms->nwfs], pin[1], parms->sim.dt);
				}
			}
			dcellfree(ftmp);
		} else{//Copy directly
			servo_output(simu->fsmint, &simu->fsmreal);
		}
		if(parms->sim.commonfsm&&simu->fsmerr){//not good
			warning_once("Using common fsm\n");
			for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
				if(parms->powfs[ipowfs].llt){
					dmat* fsmerr=0;
					real scale=1./parms->powfs[ipowfs].nwfs;
					for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
						int iwfs=P(parms->powfs[ipowfs].wfs,jwfs);
						dadd(&fsmerr, 1, P(simu->fsmerr,iwfs), scale);
					}
					for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
						int iwfs=P(parms->powfs[ipowfs].wfs,jwfs);
						dcp(&P(simu->fsmerr,iwfs), fsmerr);
					}
					dfree(fsmerr);
				}
			}
		}
		servo_filter(simu->fsmint, simu->fsmerr);
		/*Inject dithering command, for step isim+1*/
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			const int ipowfs=parms->wfs[iwfs].powfs;
			if(parms->powfs[ipowfs].dither==1){//T/T dithering.
			    //adjust delay due to propagation, and computation delay.
				const real adjust=parms->sim.alfsm+1-parms->powfs[ipowfs].dtrat+0.5;//0.5 is for testing. the value here shouldn't matter
				//Use isim+1 because the command is for next time step.
				//minus adjust for delay
				real anglei=(2*M_PI/parms->powfs[ipowfs].dither_npoint);
				real angle=(simu->wfsisim+1-adjust)*anglei;
				P(P(simu->fsmreal,iwfs),0)-=parms->powfs[ipowfs].dither_amp*cos(angle);
				P(P(simu->fsmreal,iwfs),1)-=parms->powfs[ipowfs].dither_amp*sin(angle);
			}
		}
	}
	simu->fsmerr=0;
}
/**
   filter DM commands in open loop mode by simply copy the output
 */
static void filter_ol(sim_t* simu){
	const parms_t* parms=simu->parms;
	assert(!parms->sim.closeloop);
	if(simu->dmerr&&P(parms->sim.ephi,0)>0){
		dcellcp(&simu->dmcmd, simu->dmerr);
	} else{
		dcellzero(simu->dmcmd);
	}
	if(simu->Merr_lo&&P(parms->sim.eplo,0)>0){
		addlow2dm(&simu->dmcmd, simu, simu->Merr_lo, 1);
	}
	if(simu->recon->dm_ncpa){
		warning_once("Add NCPA after integrator\n");
		dcelladd(&simu->dmcmd, 1, simu->recon->dm_ncpa, 1);
	}
	if(parms->dbg.dmoff){
		info_once("Add injected DM offset vector\n");
		int icol=(simu->reconisim+1)%NY(parms->dbg.dmoff);
		for(int idm=0; idm<parms->ndm; idm++){
			dadd(&P(simu->dmcmd,idm), 1, P(parms->dbg.dmoff, idm, icol), -1);
		}
	}
	//Extrapolate to edge actuators
	if(simu->recon->actinterp&&!parms->recon.modal){
		dcellcp(&simu->dmtmp, simu->dmcmd);
		dcellzero(simu->dmcmd);
		dspcellmm(&simu->dmcmd, simu->recon->actinterp, simu->dmtmp, "nn", 1);
	}
	if(simu->ttmreal){
		ttsplit_do(simu->recon, simu->dmcmd, simu->ttmreal, parms->sim.lpttm);
	}

	/*hysterisis. */
	if(simu->hyst){
		hyst_dcell(simu->hyst, simu->dmreal, simu->dmcmd);
	} else{
		dcellcp(&simu->dmreal, simu->dmcmd);
	}
	/*moao DM is already taken care of automatically.*/
}
/**
   Simulate turbulence on the DM
*/
void turb_dm(sim_t* simu){
	const parms_t* parms=simu->parms;
	if(!simu->dmadd) return;
	for(int idm=0; idm<parms->ndm; idm++){
		if(!P(simu->dmadd,idm)) continue;
		real* restrict p2=P(P(simu->dmreal,idm));
		const int icol=(simu->reconisim+1)%P(simu->dmadd,idm)->ny;
		const real* p=P(P(simu->dmadd,idm))+P(simu->dmadd,idm)->nx*icol;
		if(P(simu->dmadd,idm)->nx==P(simu->dmreal,idm)->nx){//match
			for(long i=0; i<P(simu->dmadd,idm)->nx; i++){
				p2[i]+=p[i];
			}
		} else{
			loc_embed_add(P(simu->dmrealsq,idm), P(simu->recon->aloc,idm), p);
		}
	}
}
/**
   Update various quantities upon updating dmreal.
*/
void update_dm(sim_t* simu){
	const parms_t* parms=simu->parms;
	if(!parms->fit.square&&simu->dmrealsq){
	/* Embed DM commands to a square array for fast ray tracing */
		for(int idm=0; idm<parms->ndm; idm++){
			loc_embed(P(simu->dmrealsq,idm), P(simu->recon->aloc,idm), P(P(simu->dmreal,idm)));
		}
	}
#if USE_CUDA
	if(parms->gpu.wfs||parms->gpu.evl){
		gpu_dmreal2gpu(simu->dmrealsq);
	}
#endif
	calc_cachedm(simu);
}

/**
   Does the servo filtering by calling filter_cl() or filter_ol()
 */
void filter_dm(sim_t* simu){
	const parms_t* parms=simu->parms;
	if(parms->sim.evlol) return;
	if(parms->sim.closeloop){
		filter_cl(simu);
	} else{
		filter_ol(simu);
	}
	save_dmreal(simu);

#if USE_CUDA
	if(simu->recon->moao){
		if(parms->gpu.moao){
			gpu_moao_filter(simu);
		} else if(parms->gpu.wfs||parms->gpu.evl){//copy DM commands to GPU
			gpu_moao_2gpu(simu);
		}
	}
#endif
	turb_dm(simu);
	update_dm(simu);
	//mark no output to handle mutiple dtrat case.
	simu->dmerr=0;
	simu->Merr_lo=0;
	if(PARALLEL==2){
	//Signal perfevl and wfsgrad that dmreal is ready
		simu->dmreal_isim=simu->reconisim+2;
		simu->dmreal_count=0;//reset the counter
		pthread_cond_broadcast(&simu->dmreal_condr);
		//dbg("dmreal_isim is set to %d\n", simu->dmreal_isim);
	}
}
