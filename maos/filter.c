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
#include "sim.h"
#include "ahst.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file filter.c
   Collection of functions for servo filtering of DM commands and uplink pointing loop.
*/
/**
   Apply hysterisis. Input dmcmd is command to the DM and output dmreal is the
   actual position the DM goes to.  */
void hysterisis(HYST_T **hyst, dcell *dmreal, const dcell *dmcmd){
    if(!hyst) return;
    assert(dmcmd->ny==1);
    for(int idm=0; idm<dmcmd->nx; idm++){
	if(!hyst[idm]) continue;
	double *restrict x=dmcmd->p[idm]->p;
	double *restrict xout=dmreal->p[idm]->p;
	double *restrict xlast=hyst[idm]->xlast->p;
	double *restrict dxlast=hyst[idm]->dxlast->p;
	double *restrict x0=hyst[idm]->x0->p;
	PDMAT(hyst[idm]->ylast, ylast);
	PDMAT(hyst[idm]->y0, yy0);
	PDMAT(hyst[idm]->coeff, coeff);
	int nmod=hyst[idm]->coeff->ny;
	int naloc=dmcmd->p[idm]->nx;
	for(int ia=0; ia<naloc; ia++){
	    double dx=x[ia]-xlast[ia];
	    if(fabs(dx)>1e-14){/*There is change in command */
		if(dx*dxlast[ia]<0){
		    /*Changes in moving direction, change the initial condition */
		    x0[ia]=xlast[ia];
		    for(int imod=0; imod<nmod; imod++){
			yy0[ia][imod]=ylast[ia][imod];
		    }
		}
		double alphasc=dx>0?1:-1;/*To revert the sign of alpha when dx<0 */
		for(int imod=0; imod<nmod; imod++){
		    const double alpha=alphasc*coeff[imod][1];
		    const double alphabeta=alpha*coeff[imod][2];
		    ylast[ia][imod]=x[ia]-alphabeta+(yy0[ia][imod]-x0[ia]+alphabeta)*exp(-(x[ia]-x0[ia])/alpha);
		}
		xlast[ia]=x[ia];
		dxlast[ia]=dx;
	    }/*else: no change in voltage, no change in output. */
	    /*update output. */
	    double y=0;
	    for(int imod=0; imod<nmod; imod++){
		y+=ylast[ia][imod]*coeff[imod][0];
	    }
	    xout[ia]=y;
	}
    }
}
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
static inline void cast_tt_do(SIM_T *simu, dcell *dmint){
    const PARMS_T *parms=simu->parms;
    if(parms->sim.dmttcast && dmint){
	/*
	  Here we are simulation a Tip/Tilt Mirror by
	  casting global tip/tilt out from DM commands, do
	  the saturation and histogram analysis. 
	  
	  We then add back the tip/tilt to the DM to
	  simulate the DM tip/tilt stage, or in another
	  word, to save a ray tracing from the tip/tilt
	  mirror.
	*/
	int ndm=parms->ndm;
	dcell *ptt=dcellnew(ndm,1);
	const RECON_T *recon=simu->recon;
	for(int idm=0; idm<ndm; idm++){
	    ptt->p[idm]=dnew(3,1);
	    double *ptt1=ptt->p[idm]->p;
	    loc_calc_ptt(NULL,ptt1, recon->aloc[idm],0,
			 recon->aimcc->p[idm],NULL,
			 dmint->p[idm]->p);
	    loc_remove_ptt(dmint->p[idm]->p, 
			   ptt1,recon->aloc[idm]);
	}
        /*
	  clip integrator. This both limits the output and
	  feeds back the clip since we are acting on the integrator directly.
	*/
	if(parms->sim.dmclip){
	    for(int idm=0; idm<parms->ndm; idm++){
		int nclip=dclip(dmint->p[idm], 
				-parms->dm[idm].stroke,
				parms->dm[idm].stroke);
		if(nclip>0){
		    info("DM %d: %d actuators clipped\n", idm, nclip);
		}
	    }
	}
	if(simu->dmhist){
	    for(int idm=0; idm<parms->ndm; idm++){
		if(simu->dmhist->p[idm]){
		    dhistfill(&simu->dmhist->p[idm], dmint->p[idm],0,
			      parms->dm[idm].histbin, parms->dm[idm].histn);
		}
	    }
	}
	if(parms->save.dmpttr){/*2 cycle delay. */
	    cellarr_dcell(simu->save->dmpttr, dmint);
	}
	double totalptt[3]={0,0,0};
	for(int idm=0; idm<ndm; idm++){
	    totalptt[1]+=ptt->p[idm]->p[1];
	    totalptt[2]+=ptt->p[idm]->p[2];
	}
	/*Add tip/tilt back to the ground DM only. */
	loc_add_ptt(dmint->p[0]->p, totalptt, recon->aloc[0]);
	dcellfree(ptt);
    }else if(parms->save.dmpttr){
	cellarr_dcell(simu->save->dmpttr, NULL);
    }
}
/**
   Update DM command for next cycle using info from last cycle (two cycle delay)
in closed loop mode */
void filter_cl(SIM_T *simu){
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
    if(parms->save.dm){
	cellarr_dcell(simu->save->dmreal, simu->dmreal);
	cellarr_dcell(simu->save->dmcmd, simu->dmcmd);
    }
    /*copy dm computed in last cycle. This is used in next cycle (already after perfevl) */
    const SIM_CFG_T *simcfg=&(parms->sim);
    servo_shift(simu->dmint, simcfg->apdm);
    if(!parms->sim.fuseint && parms->recon.split){
	servo_shift(simu->Mint_lo, simcfg->aplo);
    }
    if(parms->sim.mffocus){/*global focus is the 6th mode in ngsmod->Modes*/
	dcellmm(&simu->dmerr, simu->recon->ngsmod->Modes, simu->focuslpf, "nn", 1);
    }
    if(simu->dmerr){ /*High order. */
	servo_filter(simu->dmint, simu->dmerr, simu->dthi, simcfg->epdm);
    }
    if(parms->recon.split && simu->Merr_lo){ /*low order*/
	/*Low order, modal in split tomography only.  */
	servo_filter(simu->Mint_lo, simu->Merr_lo, simu->dtlo, simcfg->eplo);
	if(parms->sim.fuseint){/*accumulate to the main integrator.*/
	    addlow2dm(&simu->dmint->mint[0], simu, simu->Mint_lo->mpreint, 1);
	}
    }
    if(parms->sim.dmttcast){
	cast_tt_do(simu, simu->dmint->mint[0]);
    }
    /*The following are moved from the beginning to the end because the
      gradients are now from last step.*/
    dcellcp(&simu->dmcmd,simu->dmint->mint[0]);
    if(!parms->sim.fuseint){
	addlow2dm(&simu->dmcmd,simu,simu->Mint_lo->mint[0], 1);
    }   
    if(recon->dm_ncpa){
	dcelladd(&simu->dmcmd, 1, recon->dm_ncpa, 1);
    }
    /*hysteresis. */
    if(simu->hyst){
	hysterisis(simu->hyst, simu->dmreal, simu->dmcmd);
    }
    if(recon->actstuck){
	act_stuck_cmd(recon->aloc, simu->dmreal, recon->actstuck);
    }
    if(recon->actinterp){
	dcell *tmp=NULL;
	spcellmulmat(&tmp, recon->actinterp, simu->dmreal, 1);
	dcellcp(&simu->dmreal, tmp);
	dcellfree(tmp);
    }
    if(parms->sim.mffocus){/*gain was already applied on zoomerr*/
	dcelladd(&simu->zoomint, 1, simu->zoomerr, 1);
    }
    if(recon->moao){
	if(!parms->gpu.moao){ /*close loop filtering.*/
	    if(simu->dm_wfs){
		const int nwfs=parms->nwfs;
		for(int iwfs=0; iwfs<nwfs; iwfs++){
		    int ipowfs=parms->wfs[iwfs].powfs;
		    int imoao=parms->powfs[ipowfs].moao;
		    if(imoao<0) continue;
		    double g=parms->moao[imoao].gdm;
		    dadd(&simu->dm_wfs->p[iwfs], 1.-g, simu->dm_wfs->p[iwfs+nwfs], g);
		}
	    }
	    if(simu->dm_evl){
		const int nevl=parms->evl.nevl;
		int imoao=parms->evl.moao;
		double g=parms->moao[imoao].gdm;
		for(int ievl=0; ievl<nevl; ievl++){
		    dadd(&simu->dm_evl->p[ievl], 1.-g, simu->dm_evl->p[ievl+nevl], g);
		}
	    }
	}
    }
}
/**
   filter DM commands in open loop mode by simply copy the output
 */
void filter_ol(SIM_T *simu){
    RECON_T *recon=simu->recon;
    assert(!simu->parms->sim.closeloop);
    if(simu->dmerr){
	dcellcp(&simu->dmcmd, simu->dmerr);
    }else{
	dcellzero(simu->dmcmd);
    }
    if(simu->Merr_lo){
	addlow2dm(&simu->dmcmd, simu, simu->Merr_lo,1);
    }
    if(simu->parms->sim.dmttcast){
	cast_tt_do(simu, simu->dmcmd);
    }
    /*hysterisis. */
    if(simu->hyst){
	hysterisis(simu->hyst, simu->dmreal, simu->dmcmd);
    }
    if(recon->actstuck){
	act_stuck_cmd(recon->aloc, simu->dmreal, recon->actstuck);
    }
    if(recon->actinterp){
	dcell *tmp=NULL;
	spcellmulmat(&tmp, recon->actinterp, simu->dmreal, 1);
	dcellcp(&simu->dmreal, tmp);
	dcellfree(tmp);
    }
    if(simu->parms->save.dm){
	cellarr_dcell(simu->save->dmreal, simu->dmreal);
	cellarr_dcell(simu->save->dmcmd, simu->dmcmd);
    }
    /*moao DM is already taken care of automatically.*/
}

/**
   Update various quantities upon updating dmreal.
*/
void update_dm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    if(!parms->fit.square){
	/* Embed DM commands to a square array for fast ray tracing */
	for(int idm=0; idm<parms->ndm; idm++){
	    long *embed=simu->recon->aembed[idm];
	    double *pout=simu->dmrealsq[idm]->p;
	    double *pin=simu->dmreal->p[idm]->p;
	    for(long i=0; i<simu->dmreal->p[idm]->nx; i++){
		pout[embed[i]]=pin[i];
	    }
	}
    }
#if USE_CUDA
    if(parms->gpu.wfs || parms->gpu.evl){
	gpu_dmreal2gpu(simu->dmrealsq, parms->ndm,NULL);
    }
#endif
    calc_cachedm(simu);
    if(parms->plot.run){ /*Moved from recon.c to here. */
	for(int idm=0; simu->dmreal && idm<parms->ndm; idm++){
	    drawopd("DM", simu->recon->aloc[idm], simu->dmreal->p[idm]->p,NULL,
		    "Actual DM Actuator Commands","x (m)", "y (m)", "Real %d",idm);
	}
    }
}

/**
   Does the servo filtering by calling filter_cl() or filter_ol()
 */
void filter(SIM_T *simu){
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
	}else{
	    gpu_moao_2gpu(simu);
	}
    }
#endif
    update_dm(simu);

}
