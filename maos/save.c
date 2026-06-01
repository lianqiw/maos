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
#include "save.h"
#include "sim_utils.h"
#include "sim.h"
#include "plot_utils.h"
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
/**
   Propagate the atmosphere to closest xloc. skip wavefront sensing and
   reconstruction.

   2011-04-26: opdx was incorrectly computed when atm.ht and atmr.ht does not
   match in number. Fixed. Do not do scaling even if fit.ht is less.

*/
static void atm2xloc(dcell** opdx, const sim_t* simu){
	const recon_t* recon=simu->recon;
	const parms_t* parms=simu->parms;
	if(parms->recon.glao){
		return;
	}
	/*in close loop mode, opdr is from last time step. */
	int isim=simu->reconisim;
	if(!*opdx){
		*opdx=dcellnew(recon->npsr, 1);
	}
	for(int ipsr=0; ipsr<recon->npsr; ipsr++){
		if(!P(*opdx, ipsr)){
			P(*opdx, ipsr)=dnew(P(recon->xloc, ipsr)->nloc, 1);
		} else{
			dzero(P(*opdx, ipsr));
		}
	}
	if(simu->atm){
		for(int ips=0; ips<parms->atm.nps; ips++){
			real shiftx=-P(simu->atm, ips)->vx*isim*parms->sim.dt;
			real shifty=-P(simu->atm, ips)->vy*isim*parms->sim.dt;
			int ipsr=P(parms->atm.ipsr, ips);
			prop(&(propdata_t){.mapin=P(simu->atm, ips), .locout=P(recon->xloc, ipsr), .phiout=P(P(*opdx, ipsr)),
				.alpha=1, .shiftx=shiftx, .shifty=shifty, .wrap=1});
		}
	}
}
/**
   Plot and save reconstruction data.
 */
void save_recon(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	if(simu->reconisim<0) return;
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
   Plot and save dmreal
 */
void save_dmreal(sim_t* simu){
	const parms_t* parms=simu->parms;
	
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
