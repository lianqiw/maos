/**
 \file save.c Collects routines that does save and plotting to clean
 up recon.c, wfsgrad.c and perfevl.c etc.  */
#include "save.h"
#include "ahst.h"
#include "sim_utils.h"
void save_gradol(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const POWFS_T *powfs=simu->powfs;
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	int ipowfs=parms->wfs[iwfs].powfs;
	const int nsa=powfs[ipowfs].pts->nsa;
	if(!simu->gradlastol->p[iwfs]) continue;
	if(parms->plot.run){
	    drawopd("Gpolx",(loc_t*)powfs[ipowfs].pts, simu->gradlastol->p[iwfs]->p,NULL,
		    "WFS Pseudo Openloop Gradients (x)","x (m)", "y (m)", "x %d",  iwfs);
	    drawopd("Gpoly",(loc_t*)powfs[ipowfs].pts, simu->gradlastol->p[iwfs]->p+nsa, NULL,
		    "WFS Pseudo Openloop Gradients (y)","x (m)", "y (m)", "y %d",  iwfs);
	}
	if(simu->save->gradol[iwfs] && (simu->reconisim+1) % parms->powfs[ipowfs].dtrat == 0){
	    cellarr_dmat(simu->save->gradol[iwfs], simu->gradlastol->p[iwfs]);
	}
    }
    if(parms->save.ngcov>0){
	//Outputing psol gradient covariance.
	for(int igcov=0; igcov<parms->save.ngcov; igcov++){
	    int iwfs1=parms->save.gcov[igcov*2];
	    int iwfs2=parms->save.gcov[igcov*2+1];
	    info("Computing covariance between wfs %d and %d\n",iwfs1,iwfs2);
	    dmm(&simu->gcov->p[igcov], simu->gradlastol->p[iwfs1], simu->gradlastol->p[iwfs2],"nt",1);
	}
    }
}

void save_recon(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(parms->plot.run){
	if(parms->sim.recon==0){
	    for(int i=0; simu->opdr && i<simu->opdr->nx; i++){
		drawopd("Recon", recon->xloc[i], simu->opdr->p[i]->p, NULL,
			"Reconstructed Atmosphere","x (m)","y (m)","opdr %d",i);
	    }
	    for(int i=0; simu->dmfit_hi && i<simu->dmfit_hi->nx; i++){
		drawopd("DM", recon->aloc[i], simu->dmfit_hi->p[i]->p,NULL,
			"DM Fitting Output","x (m)", "y (m)","Fit %d",i);
	    }
	}
	for(int idm=0; simu->dmerr_hi && idm<parms->ndm; idm++){
	    drawopd("DM",recon->aloc[idm], simu->dmerr_hi->p[idm]->p,NULL,
		    "DM Error Signal (Hi)","x (m)","y (m)",
		    "Err Hi %d",idm);
	}
    }
    if(parms->plot.run && simu->Merr_lo){
	dcell *dmlo=NULL;
	switch(simu->parms->tomo.split){
	case 1:
	    ngsmod2dm(&dmlo, recon, simu->Merr_lo, 1);
	    break;
	case 2:
	    dcellmm(&dmlo, simu->recon->MVModes, simu->Merr_lo, "nn", 1);
	    break;
	}
	for(int idm=0; dmlo && idm<parms->ndm; idm++){
	    drawopd("DM",recon->aloc[idm], dmlo->p[idm]->p,NULL,
		    "DM Error Signal (Lo)","x (m)","y (m)",
		    "Err Lo %d",idm);
	}
	dcellfree(dmlo);
    }
    if(parms->sim.recon==0){//minimum variance tomo/fit reconstructor
	if(parms->save.opdr){
	    cellarr_dcell(simu->save->opdr, simu->opdr);
	}
	if(parms->save.dm){
	    cellarr_dcell(simu->save->dmfit_hi, simu->dmfit_hi);
	}
	if(parms->save.opdx || parms->plot.opdx){
	    dcell *opdx;
	    if(parms->sim.fitonly){
		opdx=simu->opdr;
	    }else{
		opdx=atm2xloc(simu);
	    }
	    if(parms->save.opdx){
		cellarr_dcell(simu->save->opdx, opdx);
	    }
	    if(parms->plot.opdx){ //draw opdx
		for(int i=0; i<opdx->nx; i++){
		    drawopd("Recon", recon->xloc[i], opdx->p[i]->p, NULL,
			    "Atmosphere Projected to XLOC","x (m)","y (m)","opdx %d",i);
		}
	    }
	    if(!parms->sim.fitonly){
		dcellfree(opdx);
	    }
	}
    }
    if(parms->save.dm){
	cellarr_dcell(simu->save->dmerr_hi, simu->dmerr_hi);
	if(simu->save->Merr_lo){
	    cellarr_dcell(simu->save->Merr_lo, simu->Merr_lo);
	}
    }
}
