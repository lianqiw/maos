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
/**
   \file filter.c
   Collection of functions for servo filtering of DM commands and uplink pointing loop.
*/
/**
   Add low order NGS modes to DM actuator commands for AHST and MVST
 */
void addlow2dm(dcell **dmval, const SIM_T *simu, 
	       const dcell *low_val, double gain){
    switch(simu->parms->tomo.split){
    case 0:
	break;//nothing to do.
    case 1:
	ngsmod2dm(dmval, simu->recon, low_val, gain);
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
    if(!parms->sim.closeloop) error("Invalid call\n");
  
    //copy dm computed in last cycle. This is used in next cycle (already after perfevl)
    const SIM_CFG_T *simt=&(parms->sim);
 
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
	    info("LoWFS DM output\n");
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
    static int initialized=0;
    static int cast_tt=0;
    static int do_clip=0;
    static int do_hist=0;
    if(!initialized){
	for(int idm=0; idm<parms->ndm; idm++){
	    if(parms->dm[idm].hist){
		do_hist=1;
	    }
	    if(!isinf(parms->dm[idm].stroke)){
		do_clip=1;
	    }
	}
	if(parms->save.dmpttr){
	    do_clip=1;
	}
	cast_tt=do_hist||do_clip;
	initialized=1;
    }
    if(cast_tt){
	if(!parms->sim.fuseint){
	    error("Sorry, clipping only works in fuseint=1 mode\n");
	}
    }
    if(cast_tt && simu->dmint[0]){
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
			 simu->dmint[0]->p[idm]->p);
	    loc_remove_ptt(simu->dmint[0]->p[idm]->p, 
			   ptt1,recon->aloc[idm]);
	}
        /*
	  clip integrator. This both limits the output and
	  feeds back the clip since we are acting on the integrator directly.
	*/
	if(do_clip){
	    for(int idm=0; idm<parms->ndm; idm++){
		int nclip=dclip(simu->dmint[0]->p[idm], 
				-parms->dm[idm].stroke,
				parms->dm[idm].stroke);
		if(nclip>0){
		    info("DM %d: %d actuators clipped\n", idm, nclip);
		}
	    }
	}
	if(do_hist){
	    if(!simu->dmhist){
		long nnx[parms->ndm];
		long nny[parms->ndm];
		for(int idm=0; idm<parms->ndm; idm++){
		    nnx[idm]=parms->dm[idm].histn;
		    nny[idm]=simu->dmint[0]->p[idm]->nx*simu->dmint[0]->p[idm]->ny;
		}
		simu->dmhist=dcellnew_mmap(parms->ndm, 1, nnx, nny, 
					   NULL, NULL,"dmhist_%d.bin",simu->seed);
	    }
	    for(int idm=0; idm<parms->ndm; idm++){
		dhistfill(&simu->dmhist->p[idm], simu->dmint[0]->p[idm],0,
			  parms->dm[idm].histbin, parms->dm[idm].histn);
	    }
	}
	if(parms->save.dmpttr){//2 cycle delay.
	    cellarr_dcell(simu->save->dmpttr, simu->dmint[0]);
	}
	double totalptt[3]={0,0,0};
	for(int idm=0; idm<ndm; idm++){
	    totalptt[1]+=ptt->p[idm]->p[1];
	    totalptt[2]+=ptt->p[idm]->p[2];
	}
	//Add tip/tilt back to the ground DM only.
	loc_add_ptt(simu->dmint[0]->p[0]->p, totalptt, recon->aloc[0]);
	dcellfree(ptt);
    }

    /*The following are moved from the beginning to the end because the
      gradients are now from last step.*/
    if(parms->sim.fuseint){
	dcellcp(&simu->dmreal,simu->dmint[0]);
    }else{
	dcellcp(&simu->dmreal,simu->dmint_hi[0]);
	addlow2dm(&simu->dmreal,simu,simu->Mint_lo[0], 1);
    }
    calc_cachedm(simu);
  
    //MOAO
    if(simu->moao_wfs){
	PDCELL(simu->moao_wfs, dmwfs);
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    dcp(&simu->moao_r_wfs->p[iwfs], dmwfs[iwfs][0]);
	}
    }
    if(simu->moao_evl){
	PDCELL(simu->moao_evl, dmevl);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    dcp(&simu->moao_r_evl->p[ievl], dmevl[ievl][0]);
	}
    }
}
/**
   filter DM commands in open loop mode by simply copy the output
 */
void filter_ol(SIM_T *simu){
    if(!simu->parms->sim.closeloop){
	if(simu->dmerr_hi){
	    dcellcp(&simu->dmreal, simu->dmerr_hi);
	}else{
	    dcellzero(simu->dmreal);
	}
	if(simu->Merr_lo){
	    addlow2dm(&simu->dmreal, simu, simu->Merr_lo,1);
	}
    }else{
	error("Invalid\n");
    }

    calc_cachedm(simu);
}
/**
   Does the servo filtering by calling filter_cl() or filter_ol()
 */
void filter(SIM_T *simu){
    if(simu->parms->sim.closeloop){
	if(simu->parms->save.dm){
	    cellarr_dcell(simu->save->dmreal, simu->dmreal);
	}
	filter_cl(simu);
    }else{
	filter_ol(simu);
	if(simu->parms->save.dm){
	    cellarr_dcell(simu->save->dmreal, simu->dmreal);
	}
    }
}
