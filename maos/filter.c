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

#include "common.h"
#include "sim.h"
#include "ahst.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file filter.c
   Collection of functions for servo filtering of DM commands.
*/
/**
   Apply hysterisis. Input dmcmd is command to the DM and output dmreal is the
   actual position the DM goes to.  */

/**
   Add low order NGS modes to DM actuator commands for AHST and MVST
 */
void addlow2dm(dcell **dmval, const SIM_T *simu, 
	       const dcell *low_val, double gain){
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
INLINE int limit_diff(double *x1, double *x2, double thres, long stuck1, long stuck2){
    double diff=*x2-*x1;
    if(fabs(diff)>thres){
	double ratio=signbit(diff)?-.49999:.49999;
	if(stuck1){
	    *x2=*x1+thres*ratio*2;
	}else if(stuck2){
	    *x1=*x2-thres*ratio*2;
	}else{
	    double mean=0.5*(*x1+*x2);
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
INLINE void ttsplit_do(RECON_T *recon, dcell *dmcmd, dmat *ttm, double lp){
#if 1
    int ndm=dmcmd->nx;
    double totaltt[2]={0,0};
    double tt1[2];
    for(int idm=0; idm<ndm; idm++){
	tt1[0]=tt1[1]=0;
	dmulvec(tt1, recon->DMPTT->p[idm], dmcmd->p[idm]->p,1);
	dmulvec(dmcmd->p[idm]->p, recon->DMTT->p[idm], tt1, -1);
	for(int i=0; i<2; i++){
	    totaltt[i]+=tt1[i];
	}
    }
    ttm->p[0]=ttm->p[0]*(1-lp)+lp*totaltt[0];
    ttm->p[1]=ttm->p[1]*(1-lp)+lp*totaltt[1];
    totaltt[0]-=ttm->p[0];
    totaltt[1]-=ttm->p[1];
    //Put HPF'ed to ground DM.
    dmulvec(dmcmd->p[0]->p, recon->DMTT->p[0], totaltt, 1);
#else
    //Only touch ground DM
    double tt1[2]={0,0};
    dmulvec(tt1, recon->DMPTT->p[0], dmcmd->p[0]->p,1);
    ttm->p[0]=ttm->p[0]*(1-lp)+lp*tt1[0];
    ttm->p[1]=ttm->p[1]*(1-lp)+lp*tt1[1];
    dmulvec(dmcmd->p[0]->p, recon->DMTT->p[0], ttm->p, -1);
#endif
}

INLINE void clipdm(SIM_T *simu, dcell *dmcmd){
    const PARMS_T *parms=simu->parms;
    if(!dmcmd) return;
    /*
      clip integrator. This both limits the output and
      feeds back the clip since we are acting on the integrator directly.
    */
    if(parms->sim.dmclip){
	for(int idm=0; idm<parms->ndm; idm++){
	    const int nact=dmcmd->p[idm]->nx;
	    if(parms->dm[idm].stroke->nx==1){
		if(parms->dm[idm].stroke->ny!=1){
		    error("dm.stroke is in wrong format\n");
		}
		int nclip=dclip(dmcmd->p[idm], 
				-parms->dm[idm].stroke->p[0],
				parms->dm[idm].stroke->p[0]);
		if(nclip>0){
		    info2("step %d DM %d: %d actuators clipped\n", simu->isim, idm, nclip);
		}
	    }else if(parms->dm[idm].stroke->nx==nact){
		if(parms->dm[idm].stroke->ny!=2){
		    error("dm.stroke is in wrong format\n");
		}
		
		double *pcmd=dmcmd->p[idm]->p;
		double *plow=parms->dm[idm].stroke->p;
		double *phigh=parms->dm[idm].stroke->p+nact;
		for(int iact=0; iact<nact; iact++){
		    if(pcmd[iact]<plow[iact]){
			pcmd[iact]=plow[iact];
		    }else if(pcmd[iact]>phigh[iact]){
			pcmd[iact]=phigh[iact];
		    }
		}		    
	    }else{
		error("Invalid format\n");
	    }
	}
    }
    if(parms->sim.dmclipia){
	/*Clip interactuator stroke*/
	for(int idm=0; idm<parms->ndm; idm++){
	    /* Embed DM commands to a square array (borrow dmrealsq) */
	    double iastroke;
	    int nx=simu->recon->anx->p[idm];
	    dmat *dmr;
	    dmat *dm;
	    if(parms->dm[idm].iastrokescale){ //convert dm to voltage
		dm=dinterp1(parms->dm[idm].iastrokescale->p[0], 0, dmcmd->p[idm], NAN);
		iastroke=parms->dm[idm].iastroke;//voltage.
	    }else{
		dm=dmcmd->p[idm];
		iastroke=parms->dm[idm].iastroke*2;//surface to opd
	    }
	    if(!parms->fit.square){
		loc_embed(simu->dmrealsq->p[idm], simu->recon->aloc->p[idm], dm->p);
		dmr=(dmat*)simu->dmrealsq->p[idm];
	    }else{
		dmr=dm;
	    }
	    lcell *actstuck=simu->recon->actstuck;
	    long *stuck=actstuck?(actstuck->p[idm]?actstuck->p[idm]->p:0):0;
	    int count=0,trials=0;
	    do{
		count=0;
		dmat* map=(dmat*)simu->recon->amap->p[idm]/*PDMAT*/;
		for(int iy=0; iy<simu->recon->any->p[idm]-1; iy++){
		    for(int ix=0; ix<nx; ix++){
			int iact1=IND(map,ix,iy);
			int iact2=IND(map,ix,iy+1);
			if(iact1>0 && iact2>0){
			    count+=limit_diff(PIND(dmr,ix,iy), PIND(dmr,ix,iy+1), iastroke, 
					      stuck?stuck[iact1-1]:0, stuck?stuck[iact2-1]:0);
			}
		    } 
		}
		for(int iy=0; iy<simu->recon->any->p[idm]; iy++){
		    for(int ix=0; ix<nx-1; ix++){
			int iact1=IND(map,ix,iy);
			int iact2=IND(map,ix+1,iy);
			if(iact1>0 && iact2>0){
			    count+=limit_diff(PIND(dmr,ix,iy), PIND(dmr,ix+1,iy), iastroke, 
					      stuck?stuck[iact1-1]:0, stuck?stuck[iact2-1]:0);
			}
		    }
		}
		trials++;
		if(trials==1 && count>0) {
		    info2("Step %d, DM %d: %d actuators over ia limit. ", simu->isim, idm, count);
		}
	    }while(count>0 && trials<100);
	    if(trials>1){
		info2("trials=%d: %s\n", trials, count?"failed.":"success.");
	    }
	    if(!parms->fit.square){//copy data back
		loc_extract(simu->dmreal->p[idm], simu->recon->aloc->p[idm], simu->dmrealsq->p[idm]);
	    }
	    if(parms->dm[idm].iastrokescale){//convert back to opd
		dmat *dm2=dinterp1(parms->dm[idm].iastrokescale->p[1], 0, dm, NAN);
		dcp(&dmcmd->p[idm], dm2);
		dfree(dm); dfree(dm2);
	    }
	}
    }
}

/**
   Update DM command for next cycle using info from last cycle (two cycle delay)
in closed loop mode */
static void filter_cl(SIM_T *simu){
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

      2010-01-13: Implemented apdm. 
      a(n)=a(n-1)+ep*e(n-2) or 
      a(n)=0.5*(a(n-1)+a(n-2))+ep*e(n-2);
    */
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    assert(parms->sim.closeloop);
    /*copy dm computed in last cycle. This is used in next cycle (already after perfevl) */
    const SIM_CFG_T *simcfg=&(parms->sim);
    const int isim=simu->isim;
    {/*Auto adjusting epdm for testing different epdm*/
    	static int epdm_is_auto=0;
	if(simcfg->epdm->p[0]<0){
	    epdm_is_auto=1;
	    simcfg->epdm->p[0]=0.5;
	}
	if(epdm_is_auto){
	    if((isim*10)<parms->sim.end){//initial steps
		simcfg->epdm->p[0]=0.5;
	    }else if((isim*10)%parms->sim.end==0){
		simcfg->epdm->p[0]=(double)isim/(double)parms->sim.end;
		info("epdm is set to %.1f at step %d\n", simcfg->epdm->p[0], isim);
	    }
	}
    }
    
    /*Do the servo filtering. First simulate a drop frame*/
    int drop=0;
    if(simu->dmerr && parms->sim.dtrat_skip){
	if(parms->sim.dtrat_skip>0){
	    if((isim+1)%parms->sim.dtrat_skip==0){//evenly
		drop=1;
	    }
	}else if(parms->sim.dtrat_skip<0){//use random draws
	    double tmp=randu(simu->misc_rand);
	    if(tmp*(-parms->sim.dtrat_skip)<1.){
		drop=1;
	    }
	}
    }
    dcell *dmerr=0;
    if(drop){
	warning("Drop a frame at step %d\n", isim);
    }else if(simu->dmerr){
	dmerr=simu->dmerr;
    }
    //always run servo_filter even if dmerr is NULL.
    int hiout=servo_filter(simu->dmint, dmerr);
    if(parms->recon.split){ 
	/*Low order in split tomography only. fused integrator*/
	if(servo_filter(simu->Mint_lo, simu->Merr_lo) && parms->sim.fuseint){
	    /*accumulate to the main integrator.*/
	    addlow2dm(&simu->dmint->mint->p[0], simu, simu->Mint_lo->mpreint, 1);
	}
    }
    /*The following are moved from the beginning to the end because the
      gradients are now from last step.*/
    dcellcp(&simu->dmcmd0,simu->dmint->mint->p[0]);
    if(!parms->sim.fuseint){
	addlow2dm(&simu->dmcmd0,simu,simu->Mint_lo->mint->p[0], 1);
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	//Record dmpsol for this time step for each powfs before updating it (z^-1).
	//Do not reference the data, even for dtrat==1
	if(!parms->powfs[ipowfs].psol || !parms->powfs[ipowfs].dtrat) continue;
	double alpha=(isim % parms->powfs[ipowfs].dtrat == 0)?0:1;
	dcelladd(&simu->wfspsol->p[ipowfs], alpha, simu->dmpsol, 1./parms->powfs[ipowfs].dtrat);
    }
    dcellcp(&simu->dmpsol, simu->dmcmd0);
    if(parms->recon.modal){
	dcellzero(simu->dmcmd);
	dcellmm(&simu->dmcmd, simu->recon->amod, simu->dmcmd0, "nn", 1);
	//convert DM command from modal to zonal spae
    }else if(simu->recon->actinterp && !parms->recon.psol){
	//Extrapolate to edge actuators
	dcellzero(simu->dmcmd);
	dcellmm(&simu->dmcmd, simu->recon->actinterp, simu->dmcmd0, "nn", 1);
    }else{
	dcellcp(&simu->dmcmd, simu->dmcmd0);
    }
    
    //The DM commands are always on zonal modes from this moment

    if(simu->ttmreal){
	ttsplit_do(recon, simu->dmcmd, simu->ttmreal, parms->sim.lpttm);
    }
    if(parms->sim.focus2tel && hiout){
	dcellcp(&simu->telfocusreal, simu->telfocusint);
	dcellmm(&simu->telfocusint, recon->RFdm, simu->dmcmd, "nn", parms->sim.epfocus2tel);
    }
    if(recon->dither_m){
	//Change phase in calc_dither_amp if phase of dithering is changed
	//this is for step isim+1
	double anglei=((isim+1)/recon->dither_dtrat)*(2*M_PI/recon->dither_npoint);
	dcelladd(&simu->dmcmd, 1, recon->dither_m, sin(anglei));
    }

    if(!parms->dbg.ncpa_preload && recon->dm_ncpa){
	info_once("Add NCPA after integrator\n");
	dcelladd(&simu->dmcmd, 1, recon->dm_ncpa, 1);
    }
    if(parms->sim.dmclip || parms->sim.dmclipia || recon->actstuck){
	dcell *tmp=dcelldup(simu->dmcmd);
	if(recon->actstuck){//zero stuck actuators
	    act_stuck_cmd(recon->aloc, simu->dmerr, recon->actstuck);
	}
	clipdm(simu, simu->dmcmd);
	dcelladd(&tmp, 1, simu->dmcmd, -1); //find what is clipped
	dcelladd(&simu->dmint->mint->p[0], 1, tmp, -1);//remove from integrator (anti wind up)
	dcelladd(&simu->dmpsol, 1, tmp, -1);//also feed to PSOL (is this really necessary?)
	dcellfree(tmp);
    }
    /*This is after the integrator output and clipping*/
    if(simu->dmhist){
	for(int idm=0; idm<parms->ndm; idm++){
	    if(simu->dmhist->p[idm]){
		dhistfill(&simu->dmhist->p[idm], simu->dmcmd->p[idm],0,
			  parms->dm[idm].histbin, parms->dm[idm].histn);
	    }
	}
    }
   
    /*hysteresis. */
    if(simu->hyst){
	hyst_dcell(simu->hyst, simu->dmreal, simu->dmcmd);
    }
    
    if(recon->moao && !parms->gpu.moao){
	warning_once("moao filter implemented with LPF\n");
	if(simu->dm_wfs){
	    const int nwfs=parms->nwfs;
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		int imoao=parms->powfs[ipowfs].moao;
		if(imoao<0) continue;
		double g=parms->moao[imoao].gdm;
		dadd(&simu->dm_wfs->p[iwfs], 1-g, simu->dm_wfs->p[iwfs+nwfs], g);
	    }
	}
	if(simu->dm_evl){
	    const int nevl=parms->evl.nevl;
	    int imoao=parms->evl.moao;
	    double g=parms->moao[imoao].gdm;
	    for(int ievl=0; ievl<nevl; ievl++){
		dadd(&simu->dm_evl->p[ievl], 1-g, simu->dm_evl->p[ievl+nevl], g);
	    }
	}
    }
    if(simu->fsmint){
	/*fsmerr is from gradients from this time step. so copy before update for correct delay*/
	dcellcp(&simu->fsmreal, simu->fsmint->mint->p[0]);
	if(parms->sim.f0fsm>0){
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		if(simu->fsmreal->p[iwfs]){
		    double *pin=simu->fsmint->mint->p[0]->p[iwfs]->p;
		    simu->fsmreal->p[iwfs]->p[0]=sho_step(simu->fsmsho[iwfs], pin[0], parms->sim.dt);
		    simu->fsmreal->p[iwfs]->p[1]=sho_step(simu->fsmsho[iwfs+parms->nwfs], pin[1], parms->sim.dt);
		}
	    }
	}
	if(parms->sim.commonfsm && simu->fsmerr){
	    warning_once("Using common fsm\n");
	    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].llt){
		    dmat *fsmerr=0;
		    double scale=1./parms->powfs[ipowfs].nwfs;
		    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
			dadd(&fsmerr, 1, simu->fsmerr->p[iwfs], scale);
		    }
		    for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
			dcp(&simu->fsmerr->p[iwfs], fsmerr);
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
		const int adjust=parms->sim.alfsm+1-parms->powfs[ipowfs].dtrat;
		//Use isim+1 because the command is for next time step.
		//minus adjust for delay
		double anglei=(2*M_PI/parms->powfs[ipowfs].dither_npoint);
		double angle=((isim+1-adjust)/parms->powfs[ipowfs].dtrat)*anglei;
		simu->fsmreal->p[iwfs]->p[0]-=parms->powfs[ipowfs].dither_amp*cos(angle);
		simu->fsmreal->p[iwfs]->p[1]-=parms->powfs[ipowfs].dither_amp*sin(angle);
	    }
	}
    }
}
/**
   filter DM commands in open loop mode by simply copy the output
 */
static void filter_ol(SIM_T *simu){
    assert(!simu->parms->sim.closeloop);
    if(simu->dmerr){
	dcellcp(&simu->dmcmd, simu->dmerr);
    }else{
	dcellzero(simu->dmcmd);
    }
    if(simu->Merr_lo){
	addlow2dm(&simu->dmcmd, simu, simu->Merr_lo,1);
    }
    extern int DM_NCPA;
    if(DM_NCPA && simu->recon->dm_ncpa){
	warning_once("Add NCPA after integrator\n");
	dcelladd(&simu->dmcmd, 1, simu->recon->dm_ncpa, 1);
    }
    //Extrapolate to edge actuators
    if(simu->recon->actinterp && !simu->parms->recon.modal){
	dcellcp(&simu->dmcmd0, simu->dmcmd);
	dcellzero(simu->dmcmd);
	dcellmm(&simu->dmcmd, simu->recon->actinterp, simu->dmcmd0, "nn", 1);
    }
    if(simu->ttmreal){
	ttsplit_do(simu->recon, simu->dmcmd, simu->ttmreal, simu->parms->sim.lpttm);
    }
  
    /*hysterisis. */
    if(simu->hyst){
	hyst_dcell(simu->hyst, simu->dmreal, simu->dmcmd);
    }
    /*moao DM is already taken care of automatically.*/
}
/**
   Simulate turbulence on the DM
*/
void turb_dm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(!simu->dmadd) return;
    for(int idm=0; idm<parms->ndm; idm++){
	if(!simu->dmadd->p[idm]) continue;
	double *restrict p2=simu->dmreal->p[idm]->p;
	const int icol=(simu->isim+1)%simu->dmadd->p[idm]->ny;
	const double *p=simu->dmadd->p[idm]->p+simu->dmadd->p[idm]->nx*icol;
	if(simu->dmadd->p[idm]->nx==simu->dmreal->p[idm]->nx){//match
	    for(long i=0; i<simu->dmadd->p[idm]->nx; i++){
		p2[i]+=p[i];
	    }
	}else{
	    loc_embed_add(simu->dmrealsq->p[idm], simu->recon->aloc->p[idm], p);
	}	
    }
}
/**
   Update various quantities upon updating dmreal.
*/
void update_dm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(!parms->fit.square && simu->dmrealsq){
	/* Embed DM commands to a square array for fast ray tracing */
	for(int idm=0; idm<parms->ndm; idm++){
	    loc_embed(simu->dmrealsq->p[idm], simu->recon->aloc->p[idm], simu->dmreal->p[idm]->p);
	}
    }
#if USE_CUDA
    if(parms->gpu.wfs || parms->gpu.evl){
	gpu_dmreal2gpu(simu->dmrealsq);
    }
#endif
    calc_cachedm(simu);
}

/**
   Does the servo filtering by calling filter_cl() or filter_ol()
 */
void filter_dm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(parms->sim.evlol) return;
    if(parms->sim.closeloop){
	filter_cl(simu);
    }else{
	filter_ol(simu);
    }
#if USE_CUDA
    if(simu->recon->moao){
	if(parms->gpu.moao){
	    gpu_moao_filter(simu);
	}else if(parms->gpu.wfs || parms->gpu.evl){//copy DM commands to GPU
	    gpu_moao_2gpu(simu);
	}
    }
#endif
    turb_dm(simu);
    update_dm(simu);

    dcellzero(simu->dmerr);
    simu->dmerr=0;//mark no output.
    dcellzero(simu->Merr_lo);
    simu->Merr_lo=0;
    dcellzero(simu->fsmerr);
    simu->fsmerr=0;
}
