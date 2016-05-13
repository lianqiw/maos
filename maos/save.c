/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

void save_gradstat(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    //Save pistat in the end of simulation
    for(int iwfs=0; iwfs<simu->parms->nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	const int dtrat=parms->powfs[ipowfs].dtrat;
	double scale;
	if(parms->powfs[ipowfs].usephy){
	    scale=(simu->isim+1-simu->parms->powfs[ipowfs].phystep)/dtrat;
	}else{
	    scale=(simu->isim+1)/dtrat;
	}
	if(scale<=0) continue;	    
	if(simu->pistatout && simu->pistatout->p[iwfs]){
	    int nstep=isim+1-parms->powfs[ipowfs].pistatstart;
	    scale=1./(double)nstep;
	    dcell *pp=simu->pistatout->p[iwfs];
	    dcellscale(pp,scale);
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
	    }else{/*need peak in center */
		writebin(pp,"pistat_seed%d_wfs%d.bin", simu->seed,iwfs);
	    }
	    dcellzero(pp);
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
	    drawopd("Gpolx",powfs[ipowfs].saloc, simu->gradlastol->p[iwfs]->p,NULL,
		    "WFS Pseudo Openloop Gradients (x)","x (m)", "y (m)", "x %d",  iwfs);
	    drawopd("Gpoly",powfs[ipowfs].saloc, simu->gradlastol->p[iwfs]->p+nsa, NULL,
		    "WFS Pseudo Openloop Gradients (y)","x (m)", "y (m)", "y %d",  iwfs);
	}
	if(simu->save->gradol[iwfs] && (simu->reconisim+1) % parms->powfs[ipowfs].dtrat == 0){
	    zfarr_dmat(simu->save->gradol[iwfs], simu->reconisim/parms->powfs[ipowfs].dtrat, simu->gradlastol->p[iwfs]);
	}
    }
    if(parms->save.ngcov>0){
	/*Outputing psol gradient covariance. */
	for(int igcov=0; igcov<parms->save.ngcov; igcov++){
	    int iwfs1=parms->save.gcov->p[igcov*2];
	    int iwfs2=parms->save.gcov->p[igcov*2+1];
	    //info("Computing covariance between wfs %d and %d\n",iwfs1,iwfs2);
	    dmm(&simu->gcov->p[igcov], 1, simu->gradlastol->p[iwfs1], simu->gradlastol->p[iwfs2],"nt",1);
	}
    }
}

void save_recon(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(parms->plot.run){
	if(parms->recon.alg==0){
	    for(int i=0; simu->opdr && i<simu->opdr->nx; i++){
		if(simu->opdr->p[i]){
		    drawopd("opdr", recon->xloc->p[i], simu->opdr->p[i]->p, NULL,
			    "Reconstructed Atmosphere","x (m)","y (m)","opdr %d",i);
		}
	    }
	    for(int i=0; simu->dmfit && i<simu->dmfit->nx; i++){
		if(simu->dmfit->p[i]){
		    drawopd("DM", recon->aloc->p[i], simu->dmfit->p[i]->p,NULL,
			    "DM Fitting Output","x (m)", "y (m)","Fit %d",i);
		}
	    }
	}
	if(!parms->recon.modal){
	    for(int idm=0; simu->dmerr && idm<parms->ndm; idm++){
		if(simu->dmerr->p[idm]){
		    drawopd("DM",recon->aloc->p[idm], simu->dmerr->p[idm]->p,NULL,
			    "DM Error Signal (Hi)","x (m)","y (m)",
			    "Err Hi %d",idm);
		}
	    }
	}
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
		drawopd("moao", recon->moao[imoao].aloc->p[0], simu->dm_evl->p[ievl]->p, NULL,
			"MOAO", "x(m)", "y(m)", "Evl %d", ievl);
	    }
	}
    }
    if(parms->plot.run && simu->Merr_lo){
	dcell *dmlo=NULL;
	switch(simu->parms->recon.split){
	case 1:
	    ngsmod2dm(&dmlo, recon, simu->Merr_lo, 1);
	    break;
	case 2:
	    dcellmm(&dmlo, simu->recon->MVModes, simu->Merr_lo, "nn", 1);
	    break;
	}
	for(int idm=0; dmlo && idm<parms->ndm; idm++){
	    drawopd("DM",recon->aloc->p[idm], dmlo->p[idm]->p,NULL,
		    "DM Error Signal (Lo)","x (m)","y (m)",
		    "Err Lo %d",idm);
	}
	dcellfree(dmlo);
    }
    if(parms->plot.run && simu->Mint_lo && simu->Mint_lo->mint->p[0]){
	dcell *dmlo=NULL;
	switch(simu->parms->recon.split){
	case 1:
	    ngsmod2dm(&dmlo, recon, simu->Mint_lo->mint->p[0], 1);
	    break;
	case 2:
	    dcellmm(&dmlo, simu->recon->MVModes, simu->Mint_lo->mint->p[0], "nn", 1);
	    break;
	}
	for(int idm=0; dmlo && idm<parms->ndm; idm++){
	    drawopd("DM",recon->aloc->p[idm], dmlo->p[idm]->p,NULL,
		    "DM Integrator (Lo)","x (m)","y (m)",
		    "Int Lo %d",idm);
	}
	dcellfree(dmlo);
    }
    if(parms->recon.alg==0){/*minimum variance tomo/fit reconstructor */
	if(parms->save.opdr){
	    zfarr_dcell(simu->save->opdr, simu->reconisim, simu->opdr);
	}
	if(parms->save.opdx || parms->plot.opdx){
	    dcell *opdx=simu->opdx;
	    if(!opdx){
		atm2xloc(&opdx, simu);
	    }
	    if(parms->save.opdx){
		zfarr_dcell(simu->save->opdx, simu->isim, opdx);
	    }
	    if(parms->plot.opdx){ /*draw opdx */
		for(int i=0; i<opdx->nx; i++){
		    if(opdx->p[i]){
			drawopd("opdx", recon->xloc->p[i], opdx->p[i]->p, NULL,
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
	    zfarr_dcell(simu->save->dmfit, simu->reconisim, simu->dmfit);
	}
	if(simu->dmerr){
	    zfarr_dcell(simu->save->dmerr, simu->reconisim, simu->dmerr);
	}
	if(simu->dmint->mint->p[0]){
	    zfarr_dcell(simu->save->dmint, simu->reconisim, simu->dmint->mint->p[0]);
	}
	if(simu->Merr_lo){
	    zfarr_dcell(simu->save->Merr_lo, simu->reconisim, simu->Merr_lo);
	    if(!parms->sim.fuseint && simu->Mint_lo->mint->p[0]){
		zfarr_dcell(simu->save->Mint_lo, simu->reconisim, simu->Mint_lo->mint->p[0]);
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
	dcellzero(simu->gcov);//2015-11-04: Do not cumulative average.
    }
    if(parms->sim.psfr && CHECK_SAVE(parms->evl.psfisim, parms->sim.end-(parms->sim.closeloop?1:0), simu->reconisim, parms->sim.psfr)){
	info2("Output PSF Recon Telemetry\n");
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
	dcellzero(simu->ecov);//2015-11-04: Do not cumulative average.
    }
}

void save_dmreal(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(parms->plot.run){
    	for(int idm=0; idm<parms->ndm; idm++){
	    if(simu->dmproj && simu->dmproj->p[idm]){
		drawopd("DM",recon->aloc->p[idm], simu->dmproj->p[idm]->p,NULL,
			"ATM to DM Projection (Hi)","x (m)","y (m)",
			"Proj Hi %d",idm);
	    }
	}
	//2014-05-28: moved from filter.c to here for synchronous display with dmint.
	for(int idm=0; idm<parms->ndm; idm++){
	    if(simu->dmreal && simu->dmreal->p[idm]){
		drawopd("DM", simu->recon->aloc->p[idm], simu->dmreal->p[idm]->p,NULL,
			"Actual DM Actuator Position","x (m)", "y (m)", "Real %d",idm);
	    }
	}
    }
    if(parms->save.dm){
	zfarr_dcell(simu->save->dmreal, simu->isim, simu->dmreal);
	zfarr_dcell(simu->save->dmcmd, simu->isim, simu->dmcmd);
	if(simu->ttmreal){
	    simu->save->ttmreal->p[simu->isim*2]=simu->ttmreal->p[0];
	    simu->save->ttmreal->p[simu->isim*2+1]=simu->ttmreal->p[1];
	}
    }
}
