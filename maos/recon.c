/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "maos.h"
#include "recon.h"
#include "sim_utils.h"
#include "fdpcg.h"
#include "sim.h"
#include "recon_utils.h"
#include "mvm_client.h"
#include "ahst.h"
#include "moao.h"
#include "cn2est.h"
#include "save.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file recon.c Wavefront reconstruction and DM fitting routines. 

   Since these are related to reconstruction, we don't have access to dmreal,
   which is the *actual* location of DM actuators, and only available in
   simulation. dmint should be used for dm actuator commands.*/

/**
   Calls tomo() and fit() to do the tomography and DM fit. Do error signal and
   split tomography.  
   In closedloop mode, the gradients are from time step isim-1.
*/
void tomofit(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    int isim=simu->reconisim;
   
    if(parms->sim.idealfit){
	dcellfree(simu->opdr);
    }else if(parms->sim.idealtomo){
	atm2xloc(&simu->opdr, simu);
    }else{
	/*do tomography. */
	int maxit=parms->tomo.maxit;
	if(parms->dbg.ntomo_maxit){
	    if(isim<parms->dbg.ntomo_maxit){
		maxit=parms->dbg.tomo_maxit[isim];
		recon->RL.maxit=maxit;/*update maxit information */
		info2("Running tomo.maxit=%d\n",maxit);
	    }else{
		error("Out of range\n");
	    }
	}
#if USE_CUDA
	if(parms->gpu.tomo && parms->ndm!=0){
	    gpu_tomo(simu);
	}else
#endif
	    simu->cgres->p[0]->p[isim]=muv_solve(&simu->opdr, &recon->RL, &recon->RR, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
	
	if(parms->dbg.deltafocus){
	    if(simu->opdr && recon->RFdfx){
		//Compute the delta focus in open loop.
		dcellmm(&simu->deltafocus, recon->RFdfx, simu->opdr, "nn", 1);
	    }
	    if(parms->dbg.deltafocus==2){
		dcell *dmpsol=simu->dmpsol[parms->hipowfs[0]];
		//Compute the delta focus in closed loop.
		dcellmm(&simu->deltafocus, recon->RFdfa, dmpsol, "nn", -1);
	    }
	}
    }
    
    if(parms->ndm>0){
#if USE_CUDA
	if(parms->gpu.fit){
	    gpu_fit(simu);
	}else
#endif
	    simu->cgres->p[1]->p[isim]=muv_solve(&simu->dmfit, &recon->FL, &recon->FR, simu->opdr);
    }

    /*zero stuck actuators*/
    if(recon->actstuck){
	act_stuck_cmd(recon->aloc, simu->dmfit, recon->actstuck);
    }
    /*make floating actuators averag of neighbor.*/
    if(recon->actinterp){
	dcell *tmp=NULL;
	spcellmulmat(&tmp, recon->actinterp, simu->dmfit, 1);
	dcellcp(&simu->dmfit, tmp);
	dcellfree(tmp);
    }
    dcellcp(&simu->dmerr, simu->dmfit);/*keep dmfit for warm restart */
}
/**
   Compute pseudo open loop gradients for WFS that need. Only in close loop. In
   open loop, gradlastol is simply a reference to gradlastcl. Must add the full
   DM command to both LGS and NGS WFS grads in MV Inte or MVST mode. For dtrat>1
   cases, we copy over cl to ol in a end of multi-step integration, and add to
   the accumulated DM commands scaled by 1/dtrat. Tested works
*/
static void calc_gradol(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    dcell *dmpsol;
    if(parms->dbg.psol && !parms->sim.idealfit){
	dmpsol=simu->dmcmd;
    }else{
	dmpsol=simu->dmcmdlast;//2013-03-22
    }
    PDSPCELL(recon->GA, GA);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(!parms->powfs[ipowfs].psol) continue;
	double alpha=1.;
	if(simu->reconisim % parms->powfs[ipowfs].dtrat == 0){
	    alpha=0; /*reset accumulation. */
	}
	dcelladd(&simu->dmpsol[ipowfs], alpha, dmpsol, 1./parms->powfs[ipowfs].dtrat);
	if((simu->reconisim+1) % parms->powfs[ipowfs].dtrat == 0){/*Has output. */
	    int nindwfs=parms->recon.glao?1:parms->powfs[ipowfs].nwfs;
	    for(int indwfs=0; indwfs<nindwfs; indwfs++){
		int iwfs=parms->recon.glao?ipowfs:parms->powfs[ipowfs].wfs[indwfs];
		dcp(&simu->gradlastol->p[iwfs], simu->gradlastcl->p[iwfs]);
		for(int idm=0; idm<parms->ndm && simu->dmpsol[ipowfs]; idm++){
		    spmulmat(&simu->gradlastol->p[iwfs], GA[idm][iwfs], simu->dmpsol[ipowfs]->p[idm], 1);
		}
	    }
	}
    }
}
void recon_split(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int isim=simu->reconisim;
    int lo_output=(!parms->sim.closeloop || (isim+1)%parms->sim.dtrat_lo==0);
    if(parms->recon.split==2){
	if(!parms->gpu.tomo){
	    dcellmm(&simu->gngsmvst, recon->GXL, simu->opdr, "nn", 1./parms->sim.dtrat_lo);
	}
    }
    /*Low order WFS has output */
    if(lo_output){
	dcellzero(simu->Merr_lo);
	switch(parms->recon.split){
	case 1:{
	    NGSMOD_T *ngsmod=recon->ngsmod;
	    if(!parms->tomo.ahst_idealngs){/*Low order NGS recon. */
		dcellmm(&simu->Merr_lo,ngsmod->Rngs,simu->gradlastcl,"nn",1);
		if(parms->sim.mffocus && recon->ngsmod->nmod==6){
		    simu->ngsfocus=simu->Merr_lo->p[0]->p[5];
		    simu->Merr_lo->p[0]->p[5]=0;
		}
	    }/*else: there is ideal NGS correction done in perfevl. */
	}
	    break;
	case 2:{
	    /*A separate integrator for low order is required. Use it to form error signal*/
	    dcelladd(&simu->gradlastol, 1, simu->gngsmvst, -1);
	    dcellzero(simu->gngsmvst);/*reset accumulation. */
	    dcellmm(&simu->Merr_lo, recon->MVRngs, simu->gradlastol, "nn",1);
	    if(parms->tomo.psol){
		dcell *Mpsol_lo=simu->Mint_lo->mint[parms->dbg.psol?0:1];
		dcelladd(&simu->Merr_lo, 1., Mpsol_lo, -1);
	    }
	    if(parms->sim.mffocus){
		dcell *tmp=NULL;
		dcellmm(&tmp, recon->RFngsg, simu->gradlastcl, "nn", 1);
		dcellmm(&tmp, recon->MVFM, simu->Merr_lo, "nn", -1);
		simu->ngsfocus=tmp->p[0]->p[0]; 
		dcellfree(tmp);	
	    }
	}
	    break;
	default:
	    error("Invalid parms->recon.split: %d",parms->recon.split);
	}
	if(simu->Merr_lo && parms->sim.mffocus && parms->recon.split==1 && recon->ngsmod->nmod==6){
	    /*the global focus is handled separately.*/
	    simu->Merr_lo->p[0]->p[5]=0;
	}
	if(parms->sim.mffocus && parms->dbg.deltafocus){
	    simu->ngsfocus+=simu->deltafocus->p[0]->p[0];
	}
    }else{
	dcellfree(simu->Merr_lo);
    }
}

/**
   Wavefront reconstruction. call tomofit() to do tomo()/fit() or lsr() to do
   least square reconstruction. */
void reconstruct(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(parms->sim.evlol || !simu->gradlastcl) return;
    RECON_T *recon=simu->recon;
    int isim=simu->reconisim;
    const int nwfs=parms->nwfs;
    const int hi_output=(!parms->sim.closeloop || (isim+1)%parms->sim.dtrat_hi==0);

    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	/*New plate mode focus offset for LGS WFS. Not needed*/
	if(parms->powfs[ipowfs].llt && parms->sim.ahstfocus==2 && simu->Mint_lo->mint[1]){
	    /*In new ahst mode, the first plate scale mode contains focus for
	      lgs. But it turns out to be not necessary to remove it because the
	      HPF in the LGS path removed the influence of this focus mode. set
	      sim.ahstfocus=2 to enable adjust gradients.*/
	    double scale=simu->recon->ngsmod->scale;
	    double focus=-simu->Mint_lo->mint[1]->p[0]->p[2]*(scale-1);
	    dadd(&simu->gradlastcl->p[iwfs], 1, recon->GFall->p[ipowfs], focus);
	}

	/*Uplink FSM*/
	if(parms->powfs[ipowfs].usephy 
	   && isim>=parms->powfs[ipowfs].phystep
	   && parms->powfs[ipowfs].llt 
	   && parms->powfs[ipowfs].trs
	   && simu->gradlastcl->p[iwfs]){
	    if(!recon->PTT){
		error("powfs %d has llt, but recon->PTT is NULL",ipowfs);
	    }
	    dmat *PTT=recon->PTT->p[parms->recon.glao?(ipowfs+ipowfs*parms->npowfs):(iwfs+iwfs*nwfs)];
	    if(!PTT){
		error("powfs %d has llt, but TT removal is empty\n", ipowfs);
	    }
	    /* Compute LGS Uplink error. */
	    dmm(&simu->upterr->p[iwfs], 0, PTT, simu->gradlastcl->p[iwfs], "nn", 1);
	    /* copy upterr to output. */
	    PDMAT(simu->upterrs->p[iwfs], pupterrs);
	    pupterrs[isim][0]=simu->upterr->p[iwfs]->p[0];
	    pupterrs[isim][1]=simu->upterr->p[iwfs]->p[1];
	    /* PLL loop. We hand coded frame of 240, and dll gain of 0.5 */
	    const int nc=parms->powfs[ipowfs].dither_npll;
	    const int nskip=parms->powfs[ipowfs].dither_nskip;//Loop delay
	    if(parms->powfs[ipowfs].dither && isim>=nskip){
		double angle=M_PI*0.5*isim/parms->powfs[ipowfs].dtrat;
		DITHER_T *pd=simu->dither[iwfs];
		angle+=pd->deltam;
		double sd=sin(angle);
		double cd=cos(angle);
		double err=(-sd*simu->upterr->p[iwfs]->p[0]
			    +cd*simu->upterr->p[iwfs]->p[1]);
		err/=(parms->powfs[ipowfs].dither_amp*nc);
		pd->delta+=parms->powfs[ipowfs].dither_gpll*err;
		//To estimate the actual dithering amplitude.
		double *fsmcmd=simu->uptreal->p[iwfs]->p;
		pd->ipv+=(fsmcmd[0]*cd+fsmcmd[1]*sd);
		pd->qdv+=(fsmcmd[0]*sd-fsmcmd[1]*cd);
		/*Update DLL loop measurement. The delay is about 0.2 of a
		 * cycle, according to closed loop transfer function*/
		if((isim-nskip+1)%(nc*parms->powfs[ipowfs].dtrat)==0){
		    pd->deltam=pd->delta;
		    pd->a2m=sqrt(pd->ipv*pd->ipv+pd->qdv*pd->qdv)/nc;
		    info2("PLL step%d, wfs%d: deltam=%.2f cycle, a2m=%.1f mas\n",
			  isim, iwfs, pd->deltam/(0.5*M_PI), pd->a2m*206265000);
		    pd->ipv=pd->qdv=0;
		}
	    }
	}else if(simu->upterr){
	    dfree(simu->upterr->p[iwfs]);
	}/*LLT FSM*/
    }/*for iwfs*/

    /*high pass filter lgs focus to remove sodium range variation effect*/
    if(parms->sim.mffocus){
	focus_tracking_grads(simu);
    }

    /*Gradient offset due to mainly NCPA calibration*/
    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	const int ipowfs=parms->wfs[iwfs].powfs;
	if(simu->powfs[ipowfs].gradoff){
	    int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
	    dadd(&simu->gradlastcl->p[iwfs], 1, simu->powfs[ipowfs].gradoff->p[wfsind], -1);
	}
    }

    if(parms->sim.closeloop){
	calc_gradol(simu);
    }else if(!simu->gradlastol){
	simu->gradlastol=dcellref(simu->gradlastcl);
    }
    save_gradol(simu);/*must be here since gradol is only calculated in this file. */
    if(parms->cn2.pair){
	cn2est_isim(recon, parms, simu->gradlastol, simu->reconisim);
    }/*if cn2est */
	
    double tk_start=myclockd();
    if(hi_output){
	if(parms->recon.mvm){
	    if(!simu->dmerr){
		simu->dmerr=dcellnew3(parms->ndm, 1, simu->recon->anloc, NULL);
	    }
	    if(parms->sim.mvmport){
		mvm_client_recon(parms, simu->dmerr, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
	    }else
#if USE_CUDA
		if(parms->gpu.tomo && parms->gpu.fit){
		    gpu_recon_mvm(simu);
		}else
#endif		
		{
		    if(!simu->dmerr){
			simu->dmerr=dcellnew2(simu->dmcmd);
		    }else{
			dcellzero(simu->dmerr);
		    }
		    //This assumes skipped WFS are in the end. \todo: fix it if not.
		    dmulvec(simu->dmerr->m->p, recon->MVM, 
			    (parms->tomo.psol?simu->gradlastol:simu->gradlastcl)->m->p,1);
		}
	}else{
	    switch(parms->recon.alg){
	    case 0:
		tomofit(simu);/*tomography and fitting. */
		break;
	    case 1:
		muv_solve(&simu->dmerr,&(recon->LL), &(recon->LR), simu->gradlastcl);
		break;
	    default:
		error("recon.alg=%d is not recognized\n", parms->recon.alg);
	    }
	}
	if(parms->recon.alg==0 && parms->tomo.psol){//form error signal in PSOL mode
	    dcell *dmpsol;
	    if(parms->sim.idealfit){
		dmpsol=simu->dmcmdlast;
	    }else if(parms->sim.fuseint || parms->recon.split==1){
		dmpsol=simu->dmpsol[parms->hipowfs[0]];
	    }else{
		warning_once("Temporary solution\n");
		dmpsol=simu->dmint->mint[parms->dbg.psol?0:1];
	    }
	    dcelladd(&simu->dmerr, 1, dmpsol, -1);
	}
	if(!parms->sim.idealfit && parms->recon.split==1){/*ahst */
	    remove_dm_ngsmod(simu, simu->dmerr);
	}
	if(parms->tomo.ahst_ttr && parms->recon.split){
	    remove_dm_tt(simu, simu->dmerr);
	}
    }else{
	dcellfree(simu->dmerr);
    }
    /*low order reconstruction*/
    if(parms->recon.split){
	recon_split(simu);
    }
    /*For PSF reconstruction.*/
    if(hi_output && parms->sim.psfr && isim>=parms->evl.psfisim){
	psfr_calc(simu, simu->opdr, simu->dmpsol[parms->hipowfs[0]],
		  simu->dmerr,  simu->Merr_lo);
    }

    if(recon->moao){
#if USE_CUDA
	if(parms->gpu.moao)
	    gpu_moao_recon(simu);
	else
#endif
	    moao_recon(simu);
    }
    simu->tk_recon=myclockd()-tk_start;
    save_recon(simu);/*Moved to inside. */
}
