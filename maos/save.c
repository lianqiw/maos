/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
 \file save.c Collects routines that does save and plotting to clean
 up recon.c, wfsgrad.c and perfevl.c etc.  */
#include "save.h"
#include "ahst.h"
#include "sim_utils.h"
#include "sim.h"

void save_pistat(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	if(parms->powfs[ipowfs].pistatout){
	    //Should not use dtrat here since pistat is accumuilated at every time step.
	    const int nstep=(isim+1-parms->powfs[ipowfs].pistatstart);
	    if(nstep>0){
		dcell *pp=simu->pistatout->p[iwfs];
		dcellscale(pp,1./(double)nstep);
		if(parms->sim.skysim){/*need peak in corner */
		    for(long ic=0; ic<pp->nx*pp->ny; ic++){
			dfftshift(pp->p[ic]);
		    }
		    writebin(pp,"%s/pistat/pistat_seed%d_sa%d_x%g_y%g.bin",
			     dirskysim,simu->seed,
			     parms->powfs[ipowfs].order,
			     parms->wfs[iwfs].thetax*206265,
			     parms->wfs[iwfs].thetay*206265);
		    for(long ic=0; ic<pp->nx*pp->ny; ic++){
			dfftshift(pp->p[ic]);
		    }
		}else{/*need peak in center*/
		    writebin(pp,"pistat_%d_wfs%d.bin", simu->seed,iwfs);
		}
		dcellscale(pp, nstep);
	    }
	}
	if(parms->powfs[ipowfs].i0save){
	    const int dtrat=parms->powfs[ipowfs].dtrat;
	    const int nstep=(isim+1-parms->powfs[ipowfs].phystep)/dtrat;
	    if(nstep>0){
		dcell *pp=simu->intsout->p[iwfs];
		dcellscale(pp, 1.f/(double)nstep);
		writebin(pp, "ints_%d_wfs%d.bin", simu->seed,iwfs);
		dcellscale(pp, nstep);
	    }
	}
    }
}
void save_gradol(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const POWFS_T *powfs=simu->powfs;
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	const int nsa=powfs[ipowfs].saloc->nloc;
	if(!parms->powfs[ipowfs].psol || !simu->gradlastol->p[iwfs]) continue;
	if(parms->plot.run){
	    drawopd("Gpolx",powfs[ipowfs].saloc, simu->gradlastol->p[iwfs]->p, parms->dbg.draw_opdmax->p,
		    "WFS Pseudo Openloop Gradients (x)","x (m)", "y (m)", "x %d",  iwfs);
	    drawopd("Gpoly",powfs[ipowfs].saloc, simu->gradlastol->p[iwfs]->p+nsa, parms->dbg.draw_opdmax->p,
		    "WFS Pseudo Openloop Gradients (y)","x (m)", "y (m)", "y %d",  iwfs);
	}
	if(simu->save->gradol[iwfs] && (simu->reconisim+1) % parms->powfs[ipowfs].dtrat == 0){
	    zfarr_push(simu->save->gradol[iwfs], simu->reconisim, simu->gradlastol->p[iwfs]);
	}
    }
    if(parms->save.ngcov>0){
	/*Outputing psol gradient covariance. */
	for(int igcov=0; igcov<parms->save.ngcov; igcov++){
	    int iwfs1=parms->save.gcov->p[igcov*2];
	    int iwfs2=parms->save.gcov->p[igcov*2+1];
	    //dbg("Computing covariance between wfs %d and %d\n",iwfs1,iwfs2);
	    dmm(&simu->gcov->p[igcov], 1, simu->gradlastol->p[iwfs1], simu->gradlastol->p[iwfs2],"nt",1);
	}
    }
}

void save_recon(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(parms->plot.run){
	for(int idm=0; simu->dmerr && idm<parms->ndm; idm++){
	    if(simu->dmint->mint->p[0]->p[idm]){
		drawopd("DM",recon->aloc->p[idm], simu->dmint->mint->p[0]->p[idm]->p,NULL,
			"DM Integrator (Hi)","x (m)","y (m)",
			"Int Hi %d",idm);
	    }
	}
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
	    for(int ievl=0; ievl<parms->evl.nevl && imoao>=0; ievl++){
		drawopd("moao", recon->moao[imoao].aloc->p[0], simu->dm_evl->p[ievl]->p, parms->dbg.draw_opdmax->p,
			"MOAO", "x(m)", "y(m)", "Evl %d", ievl);
	    }
	}
    }
    if(parms->plot.run && simu->Merr_lo){
	dcell *dmlo=NULL;
	switch(simu->parms->recon.split){
	case 1:
	    dcellmm(&dmlo, recon->ngsmod->Modes, simu->Merr_lo, "nn", 1);
	    break;
	case 2:
	    dcellmm(&dmlo, recon->MVModes, simu->Merr_lo, "nn", 1);
	    break;
	}
	for(int idm=0; dmlo && idm<parms->ndm; idm++){
	    drawopd("DM",recon->aloc->p[idm], dmlo->p[idm]->p,parms->dbg.draw_opdmax->p,
		    "DM Error Signal (Lo)","x (m)","y (m)",
		    "Err Lo %d",idm);
	}
	dcellfree(dmlo);
    }
    if(parms->plot.run && simu->Mint_lo && simu->Mint_lo->mint->p[0]){
	dcell *dmlo=NULL;
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
	dcellfree(dmlo);
    }
    if(parms->recon.alg==0 && !parms->sim.idealfit && !parms->recon.glao){
	/*minimum variance tomo/fit reconstructor */
	if(parms->save.opdr){
	    zfarr_push(simu->save->opdr, simu->reconisim, simu->opdr);
	}
	if(parms->save.opdx || parms->plot.opdx){
	    dcell *opdx=simu->opdx;
	    if(!opdx){
		atm2xloc(&opdx, simu);
	    }
	    if(parms->save.opdx){
		zfarr_push(simu->save->opdx, simu->isim, opdx);
	    }
	    if(parms->plot.opdx){ /*draw opdx */
		for(int i=0; i<opdx->nx; i++){
		    if(opdx->p[i]){
			drawopd("opdx", recon->xloc->p[i], opdx->p[i]->p, parms->dbg.draw_opdmax->p,
				"Atmosphere Projected to XLOC","x (m)","y (m)","opdx %d",i);
		    }
		}
	    }
	    if(!parms->sim.idealfit){
		dcellfree(opdx);
	    }
	}
    }
    if(parms->save.dm && (!parms->sim.closeloop || simu->isim>0)){
	if(simu->dmfit){
	    zfarr_push(simu->save->dmfit, simu->reconisim, simu->dmfit);
	}
	if(simu->dmerr){
	    zfarr_push(simu->save->dmerr, simu->reconisim, simu->dmerr);
	}
	if(simu->dmint->mint->p[0]){
	    zfarr_push(simu->save->dmint, simu->reconisim, simu->dmint->mint->p[0]);
	}
	if(simu->Merr_lo){
	    zfarr_push(simu->save->Merr_lo, simu->reconisim, simu->Merr_lo->p[0]);
	    if(!parms->sim.fuseint && simu->Mint_lo->mint->p[0]){
		zfarr_push(simu->save->Mint_lo, simu->reconisim, simu->Mint_lo->mint->p[0]->p[0]);
	    }
	}
    }
    const int seed=simu->seed;
    if(parms->save.ngcov>0 && CHECK_SAVE(parms->sim.start, parms->sim.end-(parms->sim.closeloop?1:0), simu->reconisim, parms->save.gcovp)){
	double scale=1./(double)(simu->reconisim-parms->sim.start+1);
	dcellscale(simu->gcov, scale);
	for(int igcov=0; igcov<parms->save.ngcov; igcov++){
	    writebin(simu->gcov->p[igcov], "gcov_%d_wfs%ld_%ld_%d.bin", seed,
		   parms->save.gcov->p[igcov*2], parms->save.gcov->p[igcov*2+1],
		   simu->reconisim+1);
	}
	dcellscale(simu->gcov, 1./scale); //2016-06-07: Do not reset.
    }
    if(parms->save.ecov && CHECK_SAVE(parms->evl.psfisim, parms->sim.end-(parms->sim.closeloop?1:0), simu->reconisim, parms->save.ecov)){
	info("Output PSF Recon Telemetry\n");
	long nstep=simu->reconisim+1-parms->evl.psfisim;
	double scale=1./nstep;
	dcellscale(simu->ecov, scale);
	if(!parms->dbg.useopdr || parms->sim.idealfit){
	    writebin(simu->ecov, "ecov_%d_%d", seed, simu->reconisim);
	}else{/*deprecated */
	    char strht[24];
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!simu->ecov->p[ievl]) continue;
		if(isfinite(parms->evl.hs->p[ievl])){
		    snprintf(strht, 24, "_%g", parms->evl.hs->p[ievl]);
		}else{
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

void save_dmreal(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(parms->plot.run){
    	for(int idm=0; idm<parms->ndm; idm++){
	    if(simu->dmproj && simu->dmproj->p[idm]){
		drawopd("DM",recon->aloc->p[idm], simu->dmproj->p[idm]->p,parms->dbg.draw_opdmax->p,
			"ATM to DM Projection (Hi)","x (m)","y (m)",
			"Proj Hi %d",idm);
	    }
	}
	//2014-05-28: moved from filter.c to here for synchronous display with dmint.
	for(int idm=0; idm<parms->ndm; idm++){
	    if(simu->dmreal && simu->dmreal->p[idm]){
		drawopd("DM", simu->recon->aloc->p[idm], simu->dmreal->p[idm]->p,parms->dbg.draw_opdmax->p,
			"Actual DM Actuator Position","x (m)", "y (m)", "Real %d",idm);
	    }
	}
    }
    if(parms->save.dm){
	zfarr_push(simu->save->dmreal, simu->isim, simu->dmreal);
	zfarr_push(simu->save->dmcmd, simu->isim, simu->dmcmd);
	if(simu->ttmreal){
	    simu->save->ttmreal->p[simu->isim*2]=simu->ttmreal->p[0];
	    simu->save->ttmreal->p[simu->isim*2+1]=simu->ttmreal->p[1];
	}
    }
}
