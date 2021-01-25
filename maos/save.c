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
/**
   Collects routines that does save and plotting to clean
   up recon.c, wfsgrad.c and perfevl.c etc.  */
#include "save.h"
#include "ahst.h"
#include "sim_utils.h"
#include "sim.h"
/**
   Save pixel statistics for matched filter.
*/
void save_pistat(SIM_T* simu){
	const PARMS_T* parms=simu->parms;
	const int isim=simu->wfsisim;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		const int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].pistatout){
			//Should not use dtrat here since pistat is accumuilated at every time step.
			const int nstep=(isim+1-parms->powfs[ipowfs].pistatstart);
			if(nstep>0){
				dcell* pp=simu->pistatout->p[iwfs];
				dcellscale(pp, 1./(real)nstep);
				if(parms->sim.skysim){/*need peak in corner */
					for(long ic=0; ic<pp->nx*pp->ny; ic++){
						dfftshift(pp->p[ic]);
					}
					writebin(pp, "%s/pistat/pistat_seed%d_sa%d_x%g_y%g.bin",
						dirskysim, simu->seed,
						parms->powfs[ipowfs].order,
						parms->wfs[iwfs].thetax*206265,
						parms->wfs[iwfs].thetay*206265);
					for(long ic=0; ic<pp->nx*pp->ny; ic++){
						dfftshift(pp->p[ic]);
					}
				} else{/*need peak in center*/
					writebin(pp, "pistat_%d_wfs%d.bin", simu->seed, iwfs);
				}
				dcellscale(pp, nstep);
			}
		}
		if(parms->powfs[ipowfs].i0save){
			const int dtrat=parms->powfs[ipowfs].dtrat;
			const int nstep=(isim+1-parms->powfs[ipowfs].phystep)/dtrat;
			if(nstep>0){
				dcell* pp=simu->intsout->p[iwfs];
				dcellscale(pp, 1.f/(real)nstep);
				writebin(pp, "ints_%d_wfs%d.bin", simu->seed, iwfs);
				dcellscale(pp, nstep);
			}
		}
	}
}
/**
   Save open loop gradients to file and optionally occumulate gcov.
 */
void save_gradol(SIM_T* simu){
	const PARMS_T* parms=simu->parms;
	const POWFS_T* powfs=simu->powfs;
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		int ipowfs=parms->wfsr[iwfs].powfs;
		if(!parms->powfs[ipowfs].psol||!simu->gradlastol->p[iwfs]) continue;
		if(parms->plot.run){
			drawgrad("Gpol", powfs[ipowfs].saloc, simu->gradlastol->p[iwfs],
				parms->plot.grad2opd, parms->dbg.draw_opdmax->p,
				"WFS Pseudo Openloop Gradients ", "x (m)", "y (m)", "Gpol %d", iwfs);
		}
		if(simu->save->gradol[iwfs]&&(simu->reconisim+1)%parms->powfs[ipowfs].dtrat==0){
			zfarr_push(simu->save->gradol[iwfs], simu->reconisim, simu->gradlastol->p[iwfs]);
		}
	}
	if(parms->save.ngcov>0){
	/*Outputing psol gradient covariance. */
		for(int igcov=0; igcov<parms->save.ngcov; igcov++){
			int iwfs1=parms->save.gcov->p[igcov*2];
			int iwfs2=parms->save.gcov->p[igcov*2+1];
			//dbg("Computing covariance between wfs %d and %d\n",iwfs1,iwfs2);
			dmm(&simu->gcov->p[igcov], 1, simu->gradlastol->p[iwfs1], simu->gradlastol->p[iwfs2], "nt", 1);
		}
	}
}
/**
   Plot and save reconstruction data.
 */
void save_recon(SIM_T* simu){
	const PARMS_T* parms=simu->parms;
	const RECON_T* recon=simu->recon;
	if(simu->reconisim<0) return;
	if(parms->plot.run){
		if(simu->dm_wfs){
			for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
				int ipowfs=parms->wfs[iwfs].powfs;
				int imoao=parms->powfs[ipowfs].moao;
				if(imoao<0) continue;
				drawopd("moao", recon->moao[imoao].aloc->p[0], simu->dm_wfs->p[iwfs]->p, NULL,
					"MOAO", "x(m)", "y(m)", "WFS %d", iwfs);
			}
		}
		if(simu->dm_evl){
			int imoao=parms->evl.moao;
			for(int ievl=0; ievl<parms->evl.nevl&&imoao>=0; ievl++){
				drawopd("moao", recon->moao[imoao].aloc->p[0], simu->dm_evl->p[ievl]->p, parms->dbg.draw_opdmax->p,
					"MOAO", "x(m)", "y(m)", "Evl %d", ievl);
			}
		}
		for(int i=0; parms->recon.alg==0&&simu->dmfit&&i<parms->ndm; i++){
			if(simu->dmfit->p[i]){
				drawopd("DM", recon->aloc->p[i], simu->dmfit->p[i]->p, parms->dbg.draw_opdmax->p,
					"DM Fitting Output", "x (m)", "y (m)", "Fit %d", i);
			}
		}

		if(!parms->recon.modal){
			for(int idm=0; simu->dmerr&&idm<parms->ndm; idm++){
				if(simu->dmerr->p[idm]){
					drawopd("DM", recon->aloc->p[idm], simu->dmerr->p[idm]->p, parms->dbg.draw_opdmax->p,
						"DM Error Signal (Hi)", "x (m)", "y (m)",
						"Err Hi %d", idm);
				}
			}
		}

		if(parms->recon.alg==0&&simu->opdr){
			for(int i=0; i<simu->opdr->nx; i++){
				if(simu->opdr->p[i]){
					drawopd("opdr", recon->xloc->p[i], simu->opdr->p[i]->p, parms->dbg.draw_opdmax->p,
						"Reconstructed Atmosphere", "x (m)", "y (m)", "opdr %d", i);
				}
			}
		}

		if(simu->Merr_lo&&draw_current("DM", "Err lo")){
			dcell* dmlo=simu->dmtmp;
			dcellzero(dmlo);
			switch(simu->parms->recon.split){
			case 1:
				dcellmm(&dmlo, recon->ngsmod->Modes, simu->Merr_lo, "nn", 1);
				break;
			case 2:
				dcellmm(&dmlo, recon->MVModes, simu->Merr_lo, "nn", 1);
				break;
			}
			for(int idm=0; dmlo&&idm<parms->ndm; idm++){
				drawopd("DM", recon->aloc->p[idm], dmlo->p[idm]->p, parms->dbg.draw_opdmax->p,
					"DM Error Signal (Lo)", "x (m)", "y (m)",
					"Err Lo %d", idm);
			}
			//plot_points("DM", 1, NULL, simu->Merr_lo, NULL, NULL, "nn", NULL, NULL, "DM Error Signal (Lo)", "NGS Modes", "NGS Mode Strength", "Err lo");
		}
	}
	if(parms->recon.alg==0&&!parms->sim.idealfit&&!parms->recon.glao){
	/*minimum variance tomo/fit reconstructor */
		if(parms->save.opdr){
			zfarr_push(simu->save->opdr, simu->reconisim, simu->opdr);
		}
		if(parms->save.opdx||parms->plot.opdx){
			dcell* opdx=simu->opdx;
			if(!opdx){
				atm2xloc(&opdx, simu);
			}
			if(parms->save.opdx){
				zfarr_push(simu->save->opdx, simu->reconisim, opdx);
			}
			if(parms->plot.opdx){ /*draw opdx */
				for(int i=0; i<opdx->nx; i++){
					if(opdx->p[i]){
						drawopd("opdx", recon->xloc->p[i], opdx->p[i]->p, parms->dbg.draw_opdmax->p,
							"Atmosphere Projected to XLOC", "x (m)", "y (m)", "opdx %d", i);
					}
				}
			}
			if(!parms->sim.idealfit){
				dcellfree(opdx);
			}
		}
	}
	if(parms->save.dm&&(!parms->sim.closeloop||simu->reconisim>=0)){
		if(simu->dmfit){
			zfarr_push(simu->save->dmfit, simu->reconisim, simu->dmfit);
		}
		if(simu->dmerr){
			zfarr_push(simu->save->dmerr, simu->reconisim, simu->dmerr);
		}
		if(simu->Merr_lo){
			zfarr_push(simu->save->Merr_lo, simu->reconisim, simu->Merr_lo->p[0]);
		}
	}
	const int seed=simu->seed;
	if(parms->save.ngcov>0&&CHECK_SAVE(parms->sim.start, parms->sim.end-(parms->sim.closeloop?1:0), simu->reconisim, parms->save.gcovp)){
		real scale=1./(real)(simu->reconisim-parms->sim.start+1);
		dcellscale(simu->gcov, scale);
		for(int igcov=0; igcov<parms->save.ngcov; igcov++){
			writebin(simu->gcov->p[igcov], "gcov_%d_wfs%ld_%ld_%d.bin", seed,
				parms->save.gcov->p[igcov*2], parms->save.gcov->p[igcov*2+1],
				simu->reconisim+1);
		}
		dcellscale(simu->gcov, 1./scale); //2016-06-07: Do not reset.
	}
	if(parms->save.ecov&&CHECK_SAVE(parms->evl.psfisim, parms->sim.end-(parms->sim.closeloop?1:0), simu->reconisim, parms->save.ecov)){
		info2("Output PSF Recon Telemetry\n");
		long nstep=simu->reconisim+1-parms->evl.psfisim;
		real scale=1./nstep;
		dcellscale(simu->ecov, scale);
		if(!parms->dbg.useopdr||parms->sim.idealfit){
			writebin(simu->ecov, "ecov_%d_%d", seed, simu->reconisim);
		} else{/*deprecated */
			char strht[24];
			for(int ievl=0; ievl<parms->evl.nevl; ievl++){
				if(!simu->ecov->p[ievl]) continue;
				if(isfinite(parms->evl.hs->p[ievl])){
					snprintf(strht, 24, "_%g", parms->evl.hs->p[ievl]);
				} else{
					strht[0]='\0';
				}
				writebin(simu->ecov->p[ievl], "ecov_%d_x%g_y%g%s_%d.bin", seed,
					parms->evl.thetax->p[ievl]*206265,
					parms->evl.thetay->p[ievl]*206265, strht, simu->reconisim);
			}
		}
		dcellscale(simu->ecov, 1./scale); //2016-06-07: Do not reset. 
	}
}
/**
   Plot and save dmproj
*/
void save_dmproj(SIM_T* simu){
	const PARMS_T* parms=simu->parms;
	const RECON_T* recon=simu->recon;
	if(parms->save.dm){
		zfarr_push(simu->save->dmproj, simu->wfsisim, simu->dmproj);
	}
	if(parms->plot.run&&simu->dmproj){
		for(int idm=0; idm<parms->ndm; idm++){
			if(simu->dmproj->p[idm]){
				drawopd("DM", recon->aloc->p[idm], simu->dmproj->p[idm]->p, parms->dbg.draw_opdmax->p,
					"ATM to DM Projection (Hi)", "x (m)", "y (m)",
					"Proj Hi %d", idm);
			}
		}
	}
}
/**
   Plot and save dmreal
 */
void save_dmreal(SIM_T* simu){
	const PARMS_T* parms=simu->parms;
	const RECON_T* recon=simu->recon;
	if(parms->plot.run){
		if(parms->sim.closeloop){
			for(int idm=0; idm<parms->ndm; idm++){
				if(simu->dmint->mint->p[0]->p[idm]){
					drawopd("DM", recon->aloc->p[idm], simu->dmint->mint->p[0]->p[idm]->p, NULL,
						"DM Integrator (Hi)", "x (m)", "y (m)",
						"Int Hi %d", idm);
				}
			}
			/*if(!parms->sim.fuseint && simu->Mint_lo->mint->p[0]){
			dcell *dmlo=simu->dmtmp;
			dcellzero(dmlo);
			switch(simu->parms->recon.split){
			case 1:
				dcellmm(&dmlo, recon->ngsmod->Modes, simu->Mint_lo->mint->p[0], "nn", 1);
				break;
			case 2:
				dcellmm(&dmlo, recon->MVModes, simu->Mint_lo->mint->p[0], "nn", 1);
				break;
			}
			for(int idm=0; dmlo && idm<parms->ndm; idm++){
				drawopd("DM",recon->aloc->p[idm], dmlo->p[idm]->p,parms->dbg.draw_opdmax->p,
					"DM Integrator (Lo)","x (m)","y (m)",
					"Int Lo %d",idm);
			}
			}*/
			if(simu->dmreal){
				for(int idm=0; idm<parms->ndm; idm++){
					if(simu->dmreal->p[idm]){
						drawopd("DM", recon->aloc->p[idm], simu->dmreal->p[idm]->p, parms->dbg.draw_opdmax->p,
							"DM Actuator Stroke", "x (m)", "y (m)", "Real %d", idm);
					}
				}
				if(simu->ttmreal&&draw_current("DM", "Real TTM")){
					int idm=0;
					real ptt[3]={0,0,0};
					ptt[1]=simu->ttmreal->p[0];
					ptt[2]=simu->ttmreal->p[1];
					dzero(simu->dmtmp->p[idm]);
					loc_add_ptt(simu->dmtmp->p[idm]->p, ptt, recon->aloc->p[idm]);
					drawopd("DM", recon->aloc->p[idm], simu->dmtmp->p[idm]->p, parms->dbg.draw_opdmax->p,
						"DM Actuator Stroke", "x (m)", "y (m)", "Real TTM");
				}

				if(simu->cachedm){//use cachedm
					for(int idm=0; idm<parms->ndm; idm++){
						drawmap("DM", simu->cachedm->p[idm], NULL,
							"DM OPD", "x (m)", "y (m)", "RealOPD %d", idm);
					}
				}
				if(draw_current("Evldm", "OA")){

					dmat* opd=dnew(simu->aper->locs->nloc, 1);
					for(int idm=0; idm<parms->ndm; idm++){
						int ind=parms->evl.nevl*idm;
						simu->evl_propdata_dm[ind].phiout=opd->p;
						CALL_THREAD(simu->evl_prop_dm[ind], 0);
					}
					dscale(opd, -1);
					if(simu->ttmreal){
						real ptt[]={0,0,0};
						ptt[1]=simu->ttmreal->p[0];
						ptt[2]=simu->ttmreal->p[1];
						loc_add_ptt(opd->p, ptt, simu->aper->locs);
					}

					drawopd("Evldm", simu->aper->locs, opd->p, parms->dbg.draw_opdmax->p,
						"DM OPD", "x (m)", "y (m)", "OA");
					dfree(opd);
				}
			}


			for(int idm=0; idm<parms->ndm; idm++){
				if(simu->dmpsol&&simu->dmpsol->p[idm]){
					drawopd("DM", simu->recon->aloc->p[idm], simu->dmpsol->p[idm]->p, parms->dbg.draw_opdmax->p,
						"DM PSOL", "x (m)", "y (m)", "PSOL %d", idm);
				}
			}
		}
	}
	if(parms->save.dm){
		int isim=(parms->sim.closeloop?2:0)+simu->reconisim;
		if(isim>=0&&isim<parms->sim.end){
			zfarr_push(simu->save->dmreal, isim, simu->dmreal);
			zfarr_push(simu->save->dmcmd, isim, simu->dmcmd);
			if(simu->ttmreal){
				simu->save->ttmreal->p[isim*2]=simu->ttmreal->p[0];
				simu->save->ttmreal->p[isim*2+1]=simu->ttmreal->p[1];
			}
			if(parms->sim.closeloop){
				zfarr_push(simu->save->dmint, isim, simu->dmint->mint->p[0]);
				if(!parms->sim.fuseint){
					zfarr_push(simu->save->Mint_lo, isim, simu->Mint_lo->mint->p[0]->p[0]);
				}
			}
		}
	}
}
