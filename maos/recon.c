/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
    int hi_output=(!parms->sim.closeloop || parms->sim.idealfit || (isim+1)%simu->dtrat_hi==0);
    int lo_output=(!parms->sim.closeloop || (isim+1)%simu->dtrat_lo==0);
    dcell *dmpsol[2];/*Hi and Lo */
    /*
      todo: The following need to be revised to use dmpsol, which is averaged over dtrat. 
    */
    if(parms->sim.fuseint){
	dmpsol[0]=dmpsol[1]=simu->dmint[parms->dbg.psol?0:1];
    }else{
	dmpsol[0]=simu->dmint_hi[parms->dbg.psol?0:1];
	dmpsol[1]=simu->Mint_lo[parms->dbg.psol?0:1];/*This can not be simu->dmpsol[lopowfs]. */
    }
    if(hi_output){
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
		muv_solve(&simu->opdr, &recon->RL, &recon->RR, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
	    /* wind shift does not work because tomography does not give exact
	       layers. wind estimation from tomography results does not work
	       either due to same argument. Now I implemented predictive control
	       within tomography directly by altering ray tracing.*/
	}
	if(parms->ndm>0){
#if USE_CUDA
	    if(parms->gpu.fit){
		gpu_fit(simu);
	    }else
#endif
		muv_solve(&simu->dmfit_hi, &recon->FL, &recon->FR, simu->opdr);
	}
	dcellcp(&simu->dmerr_hi, simu->dmfit_hi);/*keep dmfit_hi for warm restart */
    
	/*
	  Form error signal. Make sure what is subtracted here is what is added
	  to gradcl to form gramol.
	*/
	if(parms->tomo.psol){
	    dcelladd(&simu->dmerr_hi, 1, dmpsol[0], -1);
	}
	if(!parms->sim.idealfit && parms->recon.split==1){/*ahst */
	    remove_dm_ngsmod(simu, simu->dmerr_hi);
	}
	if(parms->tomo.ahst_rtt && parms->recon.split){
	    remove_dm_tt(simu, simu->dmerr_hi);
	}

    }else{/*if high order WFS has output */
	dcellfree(simu->dmerr_hi);
    }

    if(!parms->sim.idealfit && parms->recon.split){
	if(parms->recon.split==2){
	    dcelladd(&simu->opdrmvst, 1, simu->opdr, 1./simu->dtrat_lo);
	}
	/*Low order has output */
	if(lo_output){
	    dcellzero(simu->Merr_lo);
	    switch(parms->recon.split){
	    case 1:{
		NGSMOD_T *ngsmod=recon->ngsmod;
		if(!parms->tomo.ahst_idealngs){/*Low order NGS recon. */
		    dcellmm(&simu->Merr_lo,ngsmod->Rngs,simu->gradlastcl,"nn",1);
		}/*else: there is ideal NGS correction done in perfevl. */
	    }
		break;
	    case 2:{
		dcellmm(&simu->gradlastol, recon->GXL, simu->opdrmvst, "nn",-1);
		dcellmm(&simu->Merr_lo, recon->MVRngs, simu->gradlastol, "nn",1);
		if(parms->tomo.psol) dcelladd(&simu->Merr_lo, 1., dmpsol[1], -1);
		dcellzero(simu->opdrmvst);/*reset accumulation. */
	    }
		break;
	    default:
		error("Invalid parms->recon.split: %d",parms->recon.split);
	    }
	    if(parms->sim.psfr)
		dcellcp(&simu->Merr_lo_keep, simu->Merr_lo);
	}else{
	    dcellfree(simu->Merr_lo);/*don't have output. */
	}
    }
    if(hi_output && parms->sim.psfr && isim>=parms->evl.psfisim){
	/*
	  Since the DM fitting error is in the othorgonal of DM vector
	  space. The residuals estimated here have to be in DM space only. Do
	  not use opdr which covers more than DM space.
	*/
	if(!simu->opdr || !parms->dbg.useopdr || !parms->tomo.psol){/*opdr is not available. in sim.idealfit=1 mode. */
	    psfr_calc(simu, NULL, NULL, simu->dmerr_hi, simu->Merr_lo_keep);
	}else{
	    psfr_calc(simu, simu->opdr, dmpsol[0], NULL, simu->Merr_lo_keep);
	}
    }
}
/**
   least square reconstructor
*/
void lsr(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int isim=simu->reconisim;
    const int hi_output=(!parms->sim.closeloop || (isim+1)%simu->dtrat_hi==0);
    const int lo_output=parms->recon.split && (!parms->sim.closeloop || (isim+1)%simu->dtrat_lo==0);
   
    if(hi_output){
	muv_solve(&simu->dmerr_hi,&(recon->LL), &(recon->LR), simu->gradlastcl);
	if(parms->recon.split==1){/*ahst */
	    remove_dm_ngsmod(simu, simu->dmerr_hi);
	}
    }else{/*if high order does not has output */
	dcellfree(simu->dmerr_hi);
    }
    if(parms->recon.split){ /*Split reconstruction. */
	if(lo_output){/*Low order has output */
	    dcellzero(simu->Merr_lo);
	    NGSMOD_T *ngsmod=recon->ngsmod;
	    dcellmm(&simu->Merr_lo,ngsmod->Rngs,simu->gradlastcl,"nn",1);
	    if(parms->sim.psfr)
		dcellcp(&simu->Merr_lo_keep, simu->Merr_lo);
	}else{
	    dcellfree(simu->Merr_lo);/*don't have output. */
	}
    } 
    if(hi_output && parms->sim.psfr && isim>=parms->evl.psfisim){
	psfr_calc(simu, NULL, NULL, simu->dmerr_hi, simu->Merr_lo_keep);	
    }
   
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
    dcell *dmpsol=parms->dbg.psol?simu->dmcmd:simu->dmcmdlast;
    PDSPCELL(recon->GA, GA);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(!parms->powfs[ipowfs].psol) continue;
	dcelladd(&simu->dmpsol[ipowfs], 1, dmpsol, 1./parms->powfs[ipowfs].dtrat);
	if((simu->reconisim+1) % parms->powfs[ipowfs].dtrat == 0){/*Has output. */
	    int nindwfs=parms->recon.glao?1:parms->powfs[ipowfs].nwfs;
	    for(int indwfs=0; indwfs<nindwfs; indwfs++){
		int iwfs=parms->recon.glao?ipowfs:parms->powfs[ipowfs].wfs[indwfs];
		dcp(&simu->gradlastol->p[iwfs], simu->gradlastcl->p[iwfs]);
		for(int idm=0; idm<parms->ndm && simu->dmpsol[ipowfs]; idm++){
		    spmulmat(&simu->gradlastol->p[iwfs], GA[idm][iwfs], simu->dmpsol[ipowfs]->p[idm], 1);
		}
	    }
	    dcellzero(simu->dmpsol[ipowfs]);/*reset accumulation. */
	}
    }
}
/**
   Deformable mirror control. call tomofit() to do tomo()/fit() or lsr() to do
   least square reconstruction. */
void reconstruct(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(parms->sim.evlol) return;
    RECON_T *recon=simu->recon;
    if(simu->gradlastcl){
	if(!parms->sim.idealfit && !parms->sim.evlol){
	    if(parms->sim.closeloop){
		calc_gradol(simu);
	    }else if(!simu->gradlastol){
		simu->gradlastol=dcellref(simu->gradlastcl);
	    }
	    save_gradol(simu);/*must be here since gradol is only calculated in this file. */
	}
	if(parms->cn2.pair){
	    cn2est_isim(recon, parms, simu->gradlastol, simu->reconisim);
	}/*if cn2est */
	double tk_start=myclockd();
	switch(parms->recon.alg){
	case 0:
	    tomofit(simu);/*tomography and fitting. */
	    break;
	case 1:
	case 2:
	    lsr(simu);
	    break;
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
	if(parms->sim.mffocus){
	    focus_tracking(simu);
	}
	save_recon(simu);/*Moved to inside. */
    }
}
