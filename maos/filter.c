/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
	    if(fabs(dx)>1e-14){//There is change in command
		if(dx*dxlast[ia]<0){
		    //Changes in moving direction, change the initial condition
		    x0[ia]=xlast[ia];
		    for(int imod=0; imod<nmod; imod++){
			yy0[ia][imod]=ylast[ia][imod];
		    }
		}
		double alphasc=dx>0?1:-1;//To revert the sign of alpha when dx<0
		for(int imod=0; imod<nmod; imod++){
		    const double alpha=alphasc*coeff[imod][1];
		    const double alphabeta=alpha*coeff[imod][2];
		    ylast[ia][imod]=x[ia]-alphabeta+(yy0[ia][imod]-x0[ia]+alphabeta)*exp(-(x[ia]-x0[ia])/alpha);
		}
		xlast[ia]=x[ia];
		dxlast[ia]=dx;
	    }//else: no change in voltage, no change in output.
	    //update output.
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
    switch(simu->parms->tomo.split){
    case 0:
	break;//nothing to do.
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

/**
   Do type II servo filtering, except the last integrator.
*/
static void typeII_filter(TYPEII_T *MtypeII, dmat *gain, double dtngs, dcell *Merr){
    //lead filter
    double gg,ga,gs;
    if(!MtypeII->lead){
	MtypeII->lead=dcellnew(Merr->nx, Merr->ny);
	for(long imlo=0; imlo<Merr->nx*Merr->ny; imlo++){
	    long nmod=Merr->p[imlo]->nx;
	    MtypeII->lead->p[imlo]=dnew(nmod,1);
	}
    }
    if(!MtypeII->errlast){
	MtypeII->errlast=dcellnew(Merr->nx, Merr->ny);
	for(long imlo=0; imlo<Merr->nx*Merr->ny; imlo++){
	    long nmod=Merr->p[imlo]->nx;
	    MtypeII->errlast->p[imlo]=dnew(nmod,1);
	}
    }
    PDMAT(gain, pgain);
    for(long imlo=0; imlo<Merr->nx*Merr->ny; imlo++){
	const long nmod=Merr->p[imlo]->nx;
	long indmul=0;
	//if(nmod!=5){
	//   warning("nmod != 5\n");
	//}
	if(gain->ny==nmod){
	    indmul=1;
	}else if(gain->ny==1){
	    indmul=0;
	}else if(gain->ny>1 && gain->ny < nmod){
	    indmul=0;//temporary; use only first column to filter all modes
	}else{
	    error("Wrong format\n");
	}
	double *mlead=MtypeII->lead->p[imlo]->p;
	double *merr=Merr->p[imlo]->p;
	double *merrlast=MtypeII->errlast->p[imlo]->p;
		
	for(long imod=0; imod<nmod; imod++){
	    long indm=imod * indmul;
	    gg=pgain[indm][0];
	    ga=pgain[indm][1];
	    gs=pgain[indm][2]/dtngs;
		    
	    mlead[imod] = (gg/(2*ga*gs+1))*(mlead[imod]*(2*ga*gs-1)
					    +merr[imod]*(2*gs+1)
					    -merrlast[imod]*(2*gs-1));
	}
    }
    //record Merrlast for use next time
    dcellcp(&MtypeII->errlast, Merr);
    //first integrator
    dcelladd(&MtypeII->firstint, 1, MtypeII->lead, 1);
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
	if(parms->save.dmpttr){//2 cycle delay.
	    cellarr_dcell(simu->save->dmpttr, dmint);
	}
	double totalptt[3]={0,0,0};
	for(int idm=0; idm<ndm; idm++){
	    totalptt[1]+=ptt->p[idm]->p[1];
	    totalptt[2]+=ptt->p[idm]->p[2];
	}
	//Add tip/tilt back to the ground DM only.
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
    //copy dm computed in last cycle. This is used in next cycle (already after perfevl)
    const SIM_CFG_T *simt=&(parms->sim);
    if(!simu->dmerr_hi && !(parms->tomo.split && simu->Merr_lo)){
	return;
    }
    if(parms->sim.fuseint){
	shift_inte(simt->napdm,simt->apdm,simu->dmint);
    }else{
	shift_inte(simt->napdm,simt->apdm,simu->dmint_hi);
	if(parms->tomo.split){
		shift_inte(simt->napngs,simt->apngs,simu->Mint_lo);
	}
    }
    //High order.
    if(simu->dmerr_hi){
	switch(simt->servotype_hi){
	case 1://simple servo
	    if(parms->sim.fuseint){
		dcelladd(&simu->dmint[0], 1, simu->dmerr_hi, simt->epdm);
	    }else{
		dcelladd(&simu->dmint_hi[0], 1, simu->dmerr_hi, simt->epdm);
	    }
	    break;
	default:
	    error("Not implemented yet\n");
	}
    }
    //Low order, modal in split tomography only. 
    //Merr_lo is non-empty only if in split mode and (isim+1)%dtrat==0 as governed by tomofit
    if(parms->tomo.split && simu->Merr_lo){
	switch(simt->servotype_lo){
	case 1:{
	    if(parms->sim.fuseint){
		addlow2dm(&simu->dmint[0],simu,simu->Merr_lo, simt->epngs);
	    }else{
		dcelladd(&simu->Mint_lo[0], 1, simu->Merr_lo, simt->epngs);
	    }
	}
	    break;
	case 2:{ //type II with lead filter
	    //info("LoWFS DM output\n");
	    typeII_filter(simu->MtypeII_lo, simu->gtypeII_lo, simu->dtlo, simu->Merr_lo);
	    //second integrator, merged to LGS integrator.
	    if(parms->sim.fuseint){
		addlow2dm(&simu->dmint[0], simu, simu->MtypeII_lo->firstint, 1);
	    }else{
		dcelladd(&simu->Mint_lo[0], 1, simu->MtypeII_lo->firstint, 1);
	    }
	}
	    break;
	default:
	    error("Not implemented yet");
	}
    }
    if(parms->sim.dmttcast){
	cast_tt_do(simu, simu->dmint[0]);
    }
 
    /*The following are moved from the beginning to the end because the
      gradients are now from last step.*/
    if(parms->sim.fuseint){
	dcellcp(&simu->dmcmd,simu->dmint[0]);
    }else{
	dcellcp(&simu->dmcmd,simu->dmint_hi[0]);
	addlow2dm(&simu->dmcmd,simu,simu->Mint_lo[0], 1);
    }   
    //hysterisis.
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
}
/**
   filter DM commands in open loop mode by simply copy the output
 */
void filter_ol(SIM_T *simu){
    RECON_T *recon=simu->recon;
    assert(!simu->parms->sim.closeloop);
    if(simu->dmerr_hi){
	dcellcp(&simu->dmcmd, simu->dmerr_hi);
    }else{
	dcellzero(simu->dmcmd);
    }
    if(simu->Merr_lo){
	addlow2dm(&simu->dmcmd, simu, simu->Merr_lo,1);
    }
    if(simu->parms->sim.dmttcast){
	cast_tt_do(simu, simu->dmcmd);
    }
    //hysterisis.
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
    if(use_cuda){
	gpu_dm2gpu(simu->dmrealsq, parms->ndm,NULL);
    }
#endif
    calc_cachedm(simu);
  
    if(parms->plot.run){ //Moved from recon.c to here.
	for(int idm=0; simu->dmreal && idm<parms->ndm; idm++){
	    drawopd("DM", simu->recon->aloc[idm], simu->dmreal->p[idm]->p,NULL,
		    "Actual DM Actuator Commands","x (m)", "y (m)",
		    "Real %d",idm);
	}
    }
}
