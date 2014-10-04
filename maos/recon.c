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

#include "common.h"
#include "recon.h"
#include "sim_utils.h"
#include "fdpcg.h"
#include "sim.h"
#include "recon_utils.h"
#include "mvm_client.h"
#include "ahst.h"
#include "moao.h"
#include "save.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
#undef TIMING
#define TIMING 0
#if !TIMING
#define TIC_tm
#define tic_tm
#define toc_tm(A)
#else
#define TIC_tm TIC
#define tic_tm tic
#define toc_tm(A) toc2(A);tic
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
	if(parms->dbg.tomo_maxit->nx){
	    if(isim<parms->dbg.tomo_maxit->nx){
		maxit=parms->dbg.tomo_maxit->p[isim];
		recon->RL.maxit=maxit;/*update maxit information */
		info2("Running tomo.maxit=%d\n",maxit);
	    }else{
		error("Out of range\n");
	    }
	}
	TIC_tm; tic_tm;
#if USE_CUDA
	if(parms->gpu.tomo && parms->ndm!=0){
	    gpu_tomo(simu);
	}else
#endif
	    simu->cgres->p[0]->p[isim]=muv_solve(&simu->opdr, &recon->RL, &recon->RR, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
	toc_tm("Tomography");
	if(parms->dbg.deltafocus){
	    if(simu->opdr && recon->RFdfx){
		//Compute the delta focus in open loop.
		dcellmm(&simu->deltafocus, recon->RFdfx, simu->opdr, "nn", 1);
	    }
	    if(parms->dbg.deltafocus==2){
		dcell *dmpsol=simu->wfspsol->p[parms->hipowfs->p[0]];
		//Compute the delta focus in closed loop.
		dcellmm(&simu->deltafocus, recon->RFdfa, dmpsol, "nn", -1);
	    }
	}
    }
    if(parms->ndm>0){
	TIC_tm; tic_tm;
#if USE_CUDA
	if(parms->gpu.fit){
	    gpu_fit(simu);
	}else
#endif
	    simu->cgres->p[1]->p[isim]=muv_solve(&simu->dmfit, &recon->FL, &recon->FR, simu->opdr);
	toc_tm("Fitting");
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
    const RECON_T *recon=simu->recon;
    dcell *dmpsol;
    dmpsol=simu->dmpsol;
    if(parms->sim.idealfit) return;
    PDSPCELL(recon->GA, GA);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++)
  	if(parms->powfs[ipowfs].psol){
	    double alpha=1.;
	    if(simu->reconisim % parms->powfs[ipowfs].dtrat == 0){
		alpha=0; /*reset accumulation. */
	    }
	    if(parms->powfs[ipowfs].dtrat!=1){
		dcelladd(&simu->wfspsol->p[ipowfs], alpha, dmpsol, 1./parms->powfs[ipowfs].dtrat);
	    }else if(!simu->wfspsol->p[ipowfs]){
		simu->wfspsol->p[ipowfs]=dcellref(dmpsol);
	    }
	    if((simu->reconisim+1) % parms->powfs[ipowfs].dtrat == 0){/*Has output. */
		int nindwfs=parms->recon.glao?1:parms->powfs[ipowfs].nwfs;
		for(int indwfs=0; indwfs<nindwfs; indwfs++)
#if _OPENMP >= 200805 
#pragma omp task firstprivate(indwfs, alpha, ipowfs)
#endif
		{
		    int iwfs=parms->recon.glao?ipowfs:parms->powfs[ipowfs].wfs->p[indwfs];
		    dcp(&simu->gradlastol->p[iwfs], simu->gradlastcl->p[iwfs]);
		    for(int idm=0; idm<parms->ndm && simu->wfspsol->p[ipowfs]; idm++){
			dspmulmat(&simu->gradlastol->p[iwfs], GA[idm][iwfs], 
				 simu->wfspsol->p[ipowfs]->p[idm], 1);
		    }
		}
#if _OPENMP >= 200805 
#pragma omp taskwait
#endif
	    }
	}
}
void recon_split(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int isim=simu->reconisim;
    int lo_output=parms->ntipowfs && (!parms->sim.closeloop || (isim+1)%parms->sim.dtrat_lo==0);
    if(parms->recon.split==2){
	if(!parms->gpu.tomo){
	    dcellmm(&simu->gngsmvst, recon->GXL, simu->opdr, "nn", 1./parms->sim.dtrat_lo);
	}
    }
    /*Low order WFS has output */
    if(lo_output){
	simu->Merr_lo=simu->Merr_lo_store;
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
		dcell *Mpsol_lo=simu->Mint_lo->mint->p[0];
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
	    error("Invalid parms->recon.split: %d\n",parms->recon.split);
	}
	if(simu->Merr_lo && parms->sim.mffocus && parms->recon.split==1 && recon->ngsmod->nmod==6){
	    /*the global focus is handled separately.*/
	    simu->Merr_lo->p[0]->p[5]=0;
	}
	if(parms->sim.mffocus && parms->dbg.deltafocus){
	    simu->ngsfocus+=simu->deltafocus->p[0]->p[0];
	}
    }
}

/**
   Wavefront reconstruction. call tomofit() to do tomo()/fit() or lsr() to do
   least square reconstruction. */
void reconstruct(SIM_T *simu){
    double tk_start=myclockd();
    const PARMS_T *parms=simu->parms;
    if(parms->sim.evlol || !simu->gradlastcl) return;
    RECON_T *recon=simu->recon;
    int isim=simu->reconisim;
    const int hi_output=(!parms->sim.closeloop || (isim+1)%parms->sim.dtrat_hi==0);

    if(parms->sim.closeloop){
	calc_gradol(simu);
    }else if(!simu->gradlastol){
	simu->gradlastol=dcellref(simu->gradlastcl);
    }
    save_gradol(simu);//must be here since gradol is only calculated in this file. 
    if(parms->cn2.pair){
	cn2est_isim(recon, parms, parms->cn2.psol?simu->gradlastol:simu->gradlastcl);
    }//if cn2est 
    if(hi_output){
	simu->dmerr=simu->dmerr_store;
	if(parms->recon.mvm){
	    if(parms->sim.mvmport){
		mvm_client_recon(parms, simu->dmerr, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
	    }else
#if USE_CUDA
		if(parms->gpu.tomo && parms->gpu.fit){
		    gpu_recon_mvm(simu);
		}else
#endif		
		{
		    //This assumes skipped WFS are in the end. \todo: fix it if not.
		    dmulvec(simu->dmerr->m->p, recon->MVM, 
			    (parms->tomo.psol?simu->gradlastol:simu->gradlastcl)->m->p,1);
		}
	    if(parms->plot.run){
		dcellcp(&simu->dmfit, simu->dmerr);
	    }
	}else{
	    switch(parms->recon.alg){
	    case 0:
		tomofit(simu);//tomography and fitting. 
		break;
	    case 1:
		muv_solve(&simu->dmerr,&(recon->LL), &(recon->LR), simu->gradlastcl);
		break;
	    default:
		error("recon.alg=%d is not recognized\n", parms->recon.alg);
	    }
	}
	if(parms->tomo.psol && simu->recon->actinterp){
	    /* Extrapolate to edge actuators the fitting output*/
	    dcell *tmp=NULL;
	    dspcellmulmat(&tmp, simu->recon->actinterp, simu->dmerr, 1);
	    dcellcp(&simu->dmerr, tmp);
	    dcellfree(tmp);
	}
	if(parms->recon.alg==0 && parms->tomo.psol){//form error signal in PSOL mode
	    dcell *dmpsol;
	    if(parms->sim.idealfit){
		dmpsol=simu->dmpsol;
	    }else if(parms->sim.fuseint || parms->recon.split==1){
		dmpsol=simu->wfspsol->p[parms->hipowfs->p[0]];
	    }else{
		warning_once("Temporary solution\n");
		dmpsol=simu->dmint->mint->p[0];
	    }
	    dcelladd(&simu->dmerr, 1, dmpsol, -1);
	}

	if(parms->recon.alg!=1 && !parms->sim.idealfit && parms->recon.split==1){//ahst 
	    remove_dm_ngsmod(simu, simu->dmerr);
	}
	if(parms->tomo.ahst_ttr && parms->recon.split){
	    remove_dm_tt(simu, simu->dmerr);
	}
	/*zero stuck actuators*/
	if(recon->actstuck){
	    act_stuck_cmd(recon->aloc, simu->dmerr, recon->actstuck);
	}
    }
    //low order reconstruction
    if(parms->recon.split){
	recon_split(simu);
    }
    if(simu->dmerr){ //High order. 
	//global focus is the 6th mode in ngsmod->Modes
	if(parms->sim.mffocus){
	    dcellmm(&simu->dmerr, simu->recon->ngsmod->Modes, simu->ngsfocuslpf, "nn", 1);
	    //Do LPF on NGS focus measurement to drive global focus
	    double lpfocus=parms->sim.lpfocus;
	    simu->ngsfocuslpf->p[0]->p[5]=
		simu->ngsfocuslpf->p[0]->p[5]*(1.-lpfocus)+lpfocus*simu->ngsfocus;
	}
    }
    //For PSF reconstruction.
    if(hi_output && parms->sim.psfr && isim>=parms->evl.psfisim){
	psfr_calc(simu, simu->opdr, simu->wfspsol->p[parms->hipowfs->p[0]],
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
    save_recon(simu);
    simu->tk_recon=myclockd()-tk_start;
}
