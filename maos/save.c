/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "save.h"
#include "ahst.h"
#include "sim_utils.h"
#include "sim.h"
/**
 * \file save.h
   Collects routines that does save and plotting to clean
   up recon.c, wfsgrad.c and perfevl.c etc.  */
/**
   Save pixel statistics for matched filter.
*/
void save_pistat(sim_t* simu){
	const parms_t* parms=simu->parms;
	const int isim=simu->wfsisim;
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		const int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].pistatout){
			//Should not use dtrat here since pistat is accumuilated at every time step.
			const int nstep=(isim+1-parms->powfs[ipowfs].pistatstart);
			if(nstep>0){
				dcell* pp=P(simu->pistatout,iwfs);
				dcellscale(pp, 1./(real)nstep);
				if(parms->sim.skysim){/*need peak in corner */
					for(long ic=0; ic<NX(pp)*NY(pp); ic++){
						dfftshift(P(pp,ic));
					}
					writebin(pp, "%s/pistat/pistat_seed%d_sa%d_x%g_y%g.bin",
						dirskysim, simu->seed,
						parms->powfs[ipowfs].order,
						parms->wfs[iwfs].thetax*RAD2AS,
						parms->wfs[iwfs].thetay*RAD2AS);
					for(long ic=0; ic<NX(pp)*NY(pp); ic++){
						dfftshift(P(pp,ic));
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
				dcell* pp=P(simu->intsout,iwfs);
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
void save_gradol(sim_t* simu){
	const parms_t* parms=simu->parms;
	const powfs_t* powfs=simu->powfs;
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		int ipowfs=parms->wfsr[iwfs].powfs;
		if(!parms->powfs[ipowfs].psol||!P(simu->gradlastol,iwfs)) continue;
		if(parms->plot.run&&simu->wfsisim%parms->plot.run==0){
			int jwfs=P(parms->powfs[ipowfs].wfsind, iwfs);
			drawgrad("Gpol", powfs[ipowfs].saloc, PR(powfs[ipowfs].saa, jwfs), P(simu->gradlastol,iwfs),
				parms->plot.grad2opd, parms->powfs[ipowfs].trs, parms->plot.opdmax,
				"WFS Pseudo Openloop Gradients ", "x (m)", "y (m)", "WFS %2d", iwfs);
		}
		if(simu->save->gradol[iwfs]&&(simu->reconisim+1)%parms->powfs[ipowfs].dtrat==0){
			zfarr_push(simu->save->gradol[iwfs], simu->reconisim, P(simu->gradlastol,iwfs));
		}
	}
	if(parms->save.ngcov>0){
	/*Outputing psol gradient covariance. */
		for(int igcov=0; igcov<parms->save.ngcov; igcov++){
			int iwfs1=P(parms->save.gcov,igcov*2);
			int iwfs2=P(parms->save.gcov,igcov*2+1);
			//dbg("Computing covariance between wfs %d and %d\n",iwfs1,iwfs2);
			dmm(&P(simu->gcov,igcov), 1, P(simu->gradlastol,iwfs1), P(simu->gradlastol,iwfs2), "nt", 1);
		}
	}
}

void draw_dm(const parms_t *parms, const recon_t *recon, const dcell *ac, int modal, const char *title, const char *type){
	if(!ac) return;
	for(int idm=0; idm<NX(ac); idm++){
		if(!draw_current_format("DM", "%s %d", type, idm)) continue;
		dmat *dmc=NULL;
		if(recon->amod && modal){
			dmm(&dmc, 0, P(recon->amod, idm, idm), P(ac,idm), "nn", 1);
		}else{
			dmc=dref(P(ac,idm));
		}
		drawopd("DM", P(recon->aloc, idm), dmc, parms->plot.opdmax, title, "x (m)", "y (m)", "%s %d", type, idm);
		dfree(dmc);
	}
}
void draw_dm_lo(sim_t *simu, dcell *merr, const char *title, const char *type){
	if(!simu||!merr) return;
	int added=0;
	for(int idm=0; idm<NX(simu->dmtmp); idm++){
		if(draw_current_format("DM", "%s %d", type, idm)){
			if(!added){
				added=1;//add only once
				dcellzero(simu->dmtmp);
				addlow2dm(&simu->dmtmp, simu, merr, 1);
				break;
			}
		}
	}
	if(added){
		draw_dm(simu->parms, simu->recon, simu->dmtmp, 1, title, type);
	}
}
/**
   Plot and save reconstruction data.
 */
void save_recon(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	if(simu->reconisim<0) return;
	if(parms->plot.run&&simu->reconisim%parms->plot.run==0){
		if(simu->dm_wfs){
			for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
				int ipowfs=parms->wfs[iwfs].powfs;
				int imoao=parms->powfs[ipowfs].moao;
				if(imoao<0) continue;
				drawopd("DM", P(recon->moao[imoao].aloc, 0), P(simu->dm_wfs, iwfs), parms->plot.opdmax,
					"MOAO DM Command", "x(m)", "y(m)", "WFS %2d", iwfs);
			}
		}
		if(simu->dm_evl){
			int imoao=parms->evl.moao;
			for(int ievl=0; ievl<parms->evl.nevl&&imoao>=0; ievl++){
				drawopd("DM", P(recon->moao[imoao].aloc,0), P(simu->dm_evl,ievl), parms->plot.opdmax,
					"MOAO DM Command", "x(m)", "y(m)", "Evl %d", ievl);
			}
		}
		if(parms->recon.alg==RECON_MVR&&simu->dmrecon){
			draw_dm(parms, recon, simu->dmrecon, 1, "Deformable Mirror Fitting Output", "Fit");
		}
		draw_dm(parms, recon, simu->dmerr, 1, "Deformable Mirror Error Signal (Hi)", "Err Hi");

		if(parms->recon.alg==RECON_MVR&&simu->opdr){
			for(int i=0; i<NX(simu->opdr); i++){
				if(P(simu->opdr,i)){
					drawopd("Opdr", P(recon->xloc,i), P(simu->opdr,i), parms->plot.opdmax,
						"Reconstructed Atmosphere", "x (m)", "y (m)", "Layer %d", i);
				}
			}
		}
		draw_dm_lo(simu, simu->Merr_lo, "Deformable Mirror Error Signal (Lo)", "Err Lo");
	}
	if(parms->recon.alg==RECON_MVR&&!parms->sim.idealtomo&&!parms->recon.glao){
	/*minimum variance tomo/fit reconstructor */
		if(parms->save.opdr){
			zfarr_push(simu->save->opdr, simu->reconisim, simu->opdr);
		}
		if(parms->save.opdx||parms->plot.opdx){
			dcell* opdx=NULL;
			if(!opdx){
				atm2xloc(&opdx, simu);
			}
			if(parms->save.opdx){
				zfarr_push(simu->save->opdx, simu->reconisim, opdx);
			}
			if(parms->plot.opdx){ /*draw opdx */
				for(int i=0; i<NX(opdx); i++){
					if(P(opdx,i)){
						drawopd("opdx", P(recon->xloc,i), P(opdx,i), parms->plot.opdmax,
							"Atmosphere Projected to XLOC", "x (m)", "y (m)", "opdx %d", i);
					}
				}
			}
			dcellfree(opdx);
		}
	}
	if(parms->save.dm&&(!parms->sim.closeloop||simu->reconisim>=0)){
		if(simu->dmrecon){
			zfarr_push(simu->save->dmrecon, simu->reconisim, simu->dmrecon);
		}
		if(simu->dmerr){
			zfarr_push(simu->save->dmerr, simu->reconisim, simu->dmerr);
		}
		if(simu->Merr_lo){
			zfarr_push(simu->save->Merr_lo, simu->reconisim, P(simu->Merr_lo,0));
		}
	}
	const int seed=simu->seed;
	if(parms->save.ngcov>0&&CHECK_SAVE(parms->sim.start, parms->sim.end-(parms->sim.closeloop?1:0), simu->reconisim, parms->save.gcovp)){
		real scale=1./(real)(simu->reconisim-parms->sim.start+1);
		dcellscale(simu->gcov, scale);
		for(int igcov=0; igcov<parms->save.ngcov; igcov++){
			writebin(P(simu->gcov,igcov), "gcov_%d_wfs%ld_%ld_%d.bin", seed,
				P(parms->save.gcov,igcov*2), P(parms->save.gcov,igcov*2+1),
				simu->reconisim+1);
		}
		dcellscale(simu->gcov, 1./scale); //2016-06-07: Do not reset.
	}
	if(parms->save.ecov&&CHECK_SAVE(parms->evl.psfisim, parms->sim.end-(parms->sim.closeloop?1:0), simu->reconisim, parms->save.ecov)){
		info2("Output PSF Recon Telemetry\n");
		long nstep=simu->reconisim+1-parms->evl.psfisim;
		real scale=1./nstep;
		dcellscale(simu->ecov, scale);
		if(!parms->dbg.useopdr||parms->sim.idealtomo){
			writebin(simu->ecov, "ecov_%d_%d", seed, simu->reconisim);
		} else{/*deprecated */
			char strht[24];
			for(int ievl=0; ievl<parms->evl.nevl; ievl++){
				if(!P(simu->ecov,ievl)) continue;
				if(!isinf(P(parms->evl.hs,ievl))){
					snprintf(strht, 24, "_%g", P(parms->evl.hs,ievl));
				} else{
					strht[0]='\0';
				}
				writebin(P(simu->ecov,ievl), "ecov_%d_x%g_y%g%s_%d.bin", seed,
					P(parms->evl.thetax,ievl)*RAD2AS,
					P(parms->evl.thetay,ievl)*RAD2AS, strht, simu->reconisim);
			}
		}
		dcellscale(simu->ecov, 1./scale); //2016-06-07: Do not reset.
	}
}
/**
   Plot and save dmproj
*/
void save_dmproj(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	if(parms->save.dm){
		zfarr_push(simu->save->dmproj, simu->wfsisim, simu->dmproj);
	}
	if(parms->plot.run&&simu->dmproj&&simu->reconisim>0&&simu->reconisim%parms->plot.run==0){
		draw_dm(parms, recon, simu->dmproj, 0, "ATM to DM Projection (Hi)", "Proj Hi");
	}
}
/**
   Plot and save dmreal
 */
void save_dmreal(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	if(parms->plot.run&&simu->reconisim>=0 && simu->reconisim%parms->plot.run==0){
		if(parms->sim.closeloop){
			draw_dm(parms, recon, P(simu->dmint->mintc, 0), 1, "Deformable Mirror Integrator (Hi)", "Int");
			if(simu->Mint_lo&&!parms->sim.fuseint){
				draw_dm_lo(simu, P(simu->Mint_lo->mintc,0), "Deformable Mirror Integrator (Lo)", "Int Lo");
			}
		}
		if(simu->dmreal){
			//draw_dm(parms, recon, simu->dmcmd, "Deformable Mirror Command", "Cmd");
			draw_dm(parms, recon, simu->dmreal, 0, "Deformable Mirror Command", "Real");
			/*for(int idm=0; idm<parms->ndm; idm++){
				drawmap("DM", P(simu->dmrealsq, idm), parms->plot.opdmax,
					"Deformable Mirror Real Square", "x (m)", "y (m)", "dmrealsq %d", idm);
			}*/
			if(simu->ttmreal&&draw_current("DM", "Real TTM")){
				int idm=0;
				real ptt[3]={0,0,0};
				ptt[1]=P(simu->ttmreal,0);
				ptt[2]=P(simu->ttmreal,1);
				dmat* tmp=dnew(P(recon->aloc,idm)->nloc, 1);
				loc_add_ptt(tmp, ptt, P(recon->aloc,idm));
				drawopd("DM", P(recon->aloc,idm), tmp, parms->plot.opdmax,
					"TTM Command", "x (m)", "y (m)", "Real TTM");
				dfree(tmp);
			}

			/*if(simu->cachedm){//use cachedm
				for(int idm=0; idm<parms->ndm; idm++){
					drawmap("DM", P(simu->cachedm, idm), parms->plot.opdmax,
						"Deformable Mirror OPD", "x (m)", "y (m)", "Real OPD %d", idm);
				}
			}*/
			if(draw_current("DM", "Real OPD OA")){
				dmat* opd=dnew(simu->aper->locs->nloc, 1);
				for(int idm=0; idm<parms->ndm; idm++){
					int ind=parms->evl.nevl*idm;
					simu->evl_propdata_dm[ind].phiout=P(opd);
					CALL_THREAD(simu->evl_prop_dm[ind], 0);
				}
				dscale(opd, -1);
				if(simu->ttmreal){
					real ptt[]={0,0,0};
					ptt[1]=P(simu->ttmreal,0);
					ptt[2]=P(simu->ttmreal,1);
					loc_add_ptt(opd, ptt, simu->aper->locs);
				}

				drawopd("DM", simu->aper->locs, opd, parms->plot.opdmax,
					"Deformable Mirror OPD On Axis", "x (m)", "y (m)", "Real OPD OA");
				dfree(opd);
			}
		}
		/*for(int idm=0; idm<parms->ndm; idm++){
			if(parms->recon.psol && simu->dmpsol&&P(simu->dmpsol,idm)){
				drawopd("DM", P(simu->recon->aloc,idm), P(simu->dmpsol,idm), parms->plot.opdmax,
					"Deformable Mirror PSOL", "x (m)", "y (m)", "PSOL %d", idm);
			}
		}*/
	}
	if(parms->save.dm){
		int isim=(parms->sim.closeloop?2:0)+simu->reconisim;
		if(isim>=0&&isim<parms->sim.end){
			zfarr_push(simu->save->dmreal, isim, simu->dmreal);
			zfarr_push(simu->save->dmcmd, isim, simu->dmcmd);
			if(simu->ttmreal){
				P(simu->save->ttmreal,0,isim)=P(simu->ttmreal,0);
				P(simu->save->ttmreal,1,isim)=P(simu->ttmreal,1);
			}
			if(parms->sim.closeloop){
				zfarr_push(simu->save->dmint, isim, P(simu->dmint->mint,0));
				if(!parms->sim.fuseint){
					zfarr_push(simu->save->Mint_lo, isim, P(P(simu->Mint_lo->mint,0),0));
				}
			}
		}
	}
}

void plot_gradoff(sim_t *simu, int iwfs){
	const parms_t *parms=simu->parms;
	if(parms->plot.run&&simu->wfsisim%parms->plot.run==0){
		if(iwfs<0){
			for(iwfs=0; iwfs<parms->nwfs; iwfs++){
				plot_gradoff(simu, iwfs);
			}
		}else{
			int ipowfs=parms->wfs[iwfs].powfs;
			int jwfs=P(parms->powfs[ipowfs].wfsind, iwfs);
			int draw_single_save=draw_single;
			draw_single=0;
			drawgrad("Goff", simu->powfs[ipowfs].saloc, PR(simu->powfs[ipowfs].saa, jwfs), P(simu->gradoff, iwfs),
				parms->plot.grad2opd, parms->powfs[ipowfs].trs, parms->plot.gmax,
				"WFS Offset", "x (m)", "y (m)", "WFS %2d", iwfs);
			draw_single=draw_single_save;
		}
	}
}
void plot_psf(ccell* psf2s, const char* psfname, int type, int ievl, dmat* wvl, int zlog, real psfmin){
	dmat* psftemp=NULL;
	const char* title, * tab;
	for(int iwvl=0; iwvl<NX(psf2s); iwvl++){
		switch(type){
			case 2:
				title="Science Diffraction Limited PSF";
				tab="DL";
				break;
			case 1:
				title="Science Closed Loop PSF";
				tab="CL";
				break;
			case 0:
				title="Science Open Loop PSF";
				tab="OL";
				break;
			default:
				title="PSF";
				tab="PSF";
		}
		char tabname[64];
		snprintf(tabname, sizeof(tabname), "%s%2d %.2f", tab, ievl, P(wvl,iwvl)*1e6);
		if(draw_current(psfname, tabname)){
			if(psftemp&&NX(psftemp)!=P(psf2s,iwvl)->nx){
				dfree(psftemp);
			}
			cabs22d(&psftemp, 0, P(psf2s,iwvl), 1);
			draw(psfname, (plot_opts){.image=psftemp,.zlim={psfmin,1},.zlog=zlog},title, "x", "y", "%s", tabname);
		}
	}
	dfree(psftemp);
}
