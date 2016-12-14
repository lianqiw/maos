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
void tomofit(dcell **dmout, SIM_T *simu, dcell *gradin){
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
	    gpu_tomo(simu, gradin);
	}else
#endif
	    simu->cgres->p[0]->p[isim]=muv_solve(&simu->opdr, &recon->RL, &recon->RR, gradin);
	toc_tm("Tomography");
    }
    if(parms->ndm>0){
	TIC_tm; tic_tm;
#if USE_CUDA
	if(parms->gpu.fit){
	    gpu_fit(dmout, simu);
	}else
#endif
	{
	    simu->cgres->p[1]->p[isim]=muv_solve(dmout, &recon->FL, &recon->FR, simu->opdr);
	}
	toc_tm("Fitting");
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
    const RECON_T *recon=simu->recon;
    dspcell*  GA=recon->GA/*PDSPCELL*/;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
  	if(parms->powfs[ipowfs].psol){
	    if((simu->reconisim+1) % parms->powfs[ipowfs].dtrat == 0){/*Has output. */
		int nindwfs=parms->recon.glao?1:parms->powfs[ipowfs].nwfs;
		OMPTASK_FOR(indwfs, 0, nindwfs){
		    int iwfs=parms->recon.glao?ipowfs:parms->powfs[ipowfs].wfs->p[indwfs];
		    dcp(&simu->gradlastol->p[iwfs], simu->gradlastcl->p[iwfs]);
		    for(int idm=0; idm<parms->ndm && simu->wfspsol->p[ipowfs]; idm++){
			dspmm(&simu->gradlastol->p[iwfs], IND(GA,iwfs,idm), 
			      simu->wfspsol->p[ipowfs]->p[idm], "nn", 1);
		    }
		}
		OMPTASK_END;
	    }
	}
    }
}
void recon_split(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int isim=simu->reconisim;
    int lo_output=parms->ntipowfs && (!parms->sim.closeloop || (isim+1-parms->step_lo)%parms->sim.dtrat_lo==0);
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
		    //Do LPF on focus.
		    const double lpfocus=parms->sim.lpfocuslo;
		    double ngsfocus=simu->Merr_lo->p[0]->p[5];
		    simu->ngsfocuslpf=simu->ngsfocuslpf*(1-lpfocus)+lpfocus*ngsfocus;
		    simu->Merr_lo->p[0]->p[5]=simu->ngsfocuslpf;
		}
	    }/*else: there is ideal NGS correction done in perfevl. */
	}
	    break;
	case 2:{
	    /*A separate integrator for low order is required. Use it to form error signal*/
	    dcelladd(&simu->gradlastol, 1, simu->gngsmvst, -1);
	    dcellzero(simu->gngsmvst);/*reset accumulation. */
	    dcellmm(&simu->Merr_lo, recon->MVRngs, simu->gradlastol, "nn",1);
	    if(parms->recon.psol){
		dcell *Mpsol_lo=simu->Mint_lo->mint->p[0];
		dcelladd(&simu->Merr_lo, 1., Mpsol_lo, -1);
	    }
	    if(parms->sim.mffocus){
		dcell *tmp=NULL;
		dcellmm(&tmp, recon->RFngsg, simu->gradlastcl, "nn", 1);
		dcellmm(&tmp, recon->MVFM, simu->Merr_lo, "nn", -1);
		const double lpfocus=parms->sim.lpfocuslo;
		double ngsfocus=tmp->p[0]->p[0]; 
		simu->ngsfocuslpf=simu->ngsfocuslpf*(1-lpfocus)+lpfocus*ngsfocus;
		error("Please Implement: add ngsfocus to Merr_lo");
		dcellfree(tmp);	
	    }
	}
	    break;
	default:
	    error("Invalid parms->recon.split: %d\n",parms->recon.split);
	}
    }
}
void recon_servo_update(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    if(!parms->recon.psd) return;
    if(simu->dmerr && parms->sim.dtrat_hi>0){//compute PSD on dmerr.
	const int dtrat=parms->recon.psddtrat;
	const int iacc=(simu->reconisim/parms->sim.dtrat_hi);//reconstruction steps
	const int iframe=iacc % dtrat;
	//Accumulate data history.
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    for(int idm=0; idm<parms->ndm; idm++){
		dspmulvec(PCOL(IND(simu->dmerrts, ievl), iframe), IND(recon->Herr, ievl, idm),
			  simu->dmerr->p[idm]->p, 'n', 1);
	    }
	}

	if(iframe+1==dtrat){//ready for output.
	    dmat *psd=0;
	    double dthi=parms->sim.dt*parms->sim.dtrat_hi;
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		dmat *tmp=dtrans(IND(simu->dmerrts, ievl));
		dmat *psdi=psd1dt(tmp, parms->recon.psdnseg, dthi);
		dfree(tmp);
		dadd(&psd, 1, psdi, parms->evl.wt->p[ievl]);
		dfree(psdi);
	    }
	    dcellzero(simu->dmerrts);
	    //writebin(simu->dmerrts, "dmerrts_%d", simu->reconisim);
	    //writebin(psd, "psdcli_%d", simu->reconisim);
	    //average all the PSDs
	    psd_sum(psd, 1./(psd->ny-1));
	    //writebin(psd, "psdcl_%d", simu->reconisim);
	    if(simu->dmint->ep->nx==1 && simu->dmint->ep->ny==1){
		dmat *psdol=servo_rej2ol(psd, parms->sim.dt, parms->sim.dtrat_hi, simu->dmint->ep->p[0], 0);
		dcell *coeff=servo_optim(psdol, parms->sim.dt, parms->sim.dtrat_hi, M_PI*0.25, 0, 1);
		double g=0.5;
		simu->dmint->ep->p[0]=simu->dmint->ep->p[0]*(1-g)+coeff->p[0]->p[0]*g;
		info2("Step %d New Gain (high): %.3f\n", simu->reconisim, simu->dmint->ep->p[0]);
		writebin(psdol, "psdol_%d", simu->reconisim);		    
		dcellfree(coeff);
		dfree(psdol);
	    }else{
		error("Please implement\n");
	    }
	    dfree(psd);
	}
    }
    if(parms->recon.split && simu->Merr_lo && parms->recon.psddtrat_lo>0){//compute PSD on low order control
	const int iacc=(simu->reconisim/parms->sim.dtrat_lo);//reconstruction steps
	const int dtrat=parms->recon.psddtrat_lo;
	const int iframe=iacc % dtrat;
	dmulvec(PCOL(simu->Merrts, iframe), recon->ngsmod->MCCu, simu->Merr_lo->p[0]->p, 1);
	if(iframe+1==dtrat){
	    //writebin(simu->Merrts, "Merrts_%d", simu->reconisim);
	    dmat *ts=dtrans(simu->Merrts);
	    double dt=parms->sim.dt*parms->sim.dtrat_lo;
	    dmat *psd=psd1dt(ts, parms->recon.psdnseg, dt);
	    dzero(simu->Merrts);
	    dfree(ts);
	    //Sum all the PSDs
	    psd_sum(psd, 1);
	    //writebin(psd, "psdcl_lo_%d", simu->reconisim);

	    if(simu->Mint_lo->ep->nx==1 && simu->Mint_lo->ep->nx==1){
		dmat *psdol=servo_rej2ol(psd, parms->sim.dt, parms->sim.dtrat_lo, simu->Mint_lo->ep->p[0], 0);
		dcell *coeff=servo_optim(psdol, parms->sim.dt, parms->sim.dtrat_lo, M_PI*0.25, 0, 1);
		const double g=0.5;
		simu->Mint_lo->ep->p[0]=simu->Mint_lo->ep->p[0]*(1-g)+coeff->p[0]->p[0]*g;
		info2("Step %d New Gain (low) : %.3f\n", simu->reconisim, simu->Mint_lo->ep->p[0]);
		dfree(psdol);
		dcellfree(coeff);
	    }else{
		error("Please implement\n");
	    }
	    dfree(psd);
	}
    }
}
/**
   Wavefront reconstruction. call tomofit() to do tomo()/fit() or lsr() to do
   least square reconstruction. */
void reconstruct(SIM_T *simu){
    double tk_start=myclockd();
    const PARMS_T *parms=simu->parms;
    if(parms->sim.evlol) return;
    RECON_T *recon=simu->recon;
    int isim=simu->reconisim;
    if(isim<0) return;
    const int hi_output=(!parms->sim.closeloop || (isim+1-parms->step_hi)%parms->sim.dtrat_hi==0);
    if(simu->gradlastcl){
	if(parms->sim.closeloop){
	    calc_gradol(simu);
	    save_gradol(simu);//must be here since gradol is only calculated in this file. 
	}
    	if(parms->cn2.pair){
	    cn2est_isim(recon, parms, parms->cn2.psol?simu->gradlastol:simu->gradlastcl);
	}//if cn2est 
    }
    if(hi_output || parms->sim.idealfit || parms->sim.idealtomo){
	simu->dmerr=simu->dmerr_store;
	dcell *dmout, *gradin;
	if(parms->recon.psol){
	    dmout=simu->dmfit;
	    gradin=simu->gradlastol;
	}else{
	    dmout=simu->dmerr;
	    gradin=simu->gradlastcl;
	}
	if(!dmout) error("dmout cannot be empty\n");
	//The following takes gradin as input and computs dmfit in dmout.
	if(parms->recon.mvm){
	    if(parms->sim.mvmport){
		mvm_client_recon(parms->sim.mvmsize, dmout, gradin);
	    }else
#if USE_CUDA
		if((parms->gpu.tomo && parms->gpu.fit) || parms->gpu.lsr){
		    gpu_recon_mvm(&dmout, gradin);
		}else
#endif		
		{
		    dzero(dmout->m);
		    dmulvec(dmout->m->p, recon->MVM, gradin->m->p, 1);
		}
	}else{
	    switch(parms->recon.alg){
	    case 0://MVR
		tomofit(&dmout, simu, gradin);//tomography and fitting. 
		break;
	    case 1://LSR
		if(simu->gradlastcl){
#if USE_CUDA
		    if(parms->gpu.lsr){
			error("To implement. MVM is done in as in MV case\n");
		    }else
#endif
			muv_solve(&dmout,&(recon->LL), &(recon->LR), gradin);
		}
		break;
	    default:
		error("recon.alg=%d is not recognized\n", parms->recon.alg);
	    }
	}
	if(parms->recon.psol){
	    //form error signal in PSOL mode
	    dcellcp(&simu->dmerr, simu->dmfit);/*keep dmfit for warm restart */
	    if(0){
		warning_once("temporarily disable recon->actinterp\n");
	    }else if(simu->recon->actinterp){
		//extrapolate DM fitting result to float and edge actuators
		dcellcp(&simu->dmcmd0, simu->dmerr);
		dcellzero(simu->dmerr);
		dcellmm(&simu->dmerr, simu->recon->actinterp, simu->dmcmd0, "nn", 1);
	    }

	    dcell *dmpsol;
	    if(parms->sim.idealfit || parms->sim.idealtomo){
		dmpsol=simu->dmpsol;
	    }else if(parms->sim.fuseint || parms->recon.split==1){
		dmpsol=simu->wfspsol->p[parms->hipowfs->p[0]];
	    }else{
		warning_once("Temporary solution for MVST.\n");
		dmpsol=simu->dmint->mint->p[0];
	    }
	    dcelladd(&simu->dmerr, 1, dmpsol, -1);
	}

	if(parms->recon.split){
	    if(parms->recon.alg==0){//ahst 
		remove_dm_ngsmod(simu, simu->dmerr);
	    }
	}
	
	if(recon->actstuck && !parms->recon.modal){//zero stuck actuators
	    act_stuck_cmd(recon->aloc, simu->dmerr, recon->actstuck);
	}
    }
    
    if(parms->recon.split){//low order reconstruction
	recon_split(simu);
    }
    if(parms->recon.psd){
	recon_servo_update(simu);
    }
    if(hi_output && parms->sim.psfr && isim>=parms->evl.psfisim){
	//For PSF reconstruction.
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
