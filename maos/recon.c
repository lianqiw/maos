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
#include "recon.h"
#include "sim_utils.h"
#include "fdpcg.h"
#include "sim.h"
#include "recon_utils.h"
#include "ahst.h"
#include "moao.h"

/**
   \file recon.c
   Wavefront reconstruction and DM fitting routines
*/

/**
   Calls tomo() and fit() to do the tomography and DM fit. Do error signal and
   split tomography.  */
void tomofit(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    //2010-12-16: replaced isim+1 by isim since recon is delayed from wfsgrad by 1 frame.
    if(!parms->sim.closeloop || parms->sim.fitonly || simu->dtrat_hi==1 || (simu->isim)%simu->dtrat_hi==0){
	if(parms->sim.fitonly){
	    dcellfree(simu->opdr);
	    //simu->opdr=atm2xloc(simu);
	}else{
	    int maxit=parms->tomo.maxit;
	    if(parms->dbg.ntomo_maxit){
		if(simu->isim<parms->dbg.ntomo_maxit){
		    maxit=parms->dbg.tomo_maxit[simu->isim];
		    recon->RL.maxit=maxit;//update maxit information
		    info2("Running tomo.maxit=%d\n",maxit);
		}else{
		    error("Out of range\n");
		}
	    }
	    dcell *rhs=NULL;
	    muv(&rhs, &recon->RR, simu->gradlastol, 1);
	    muv_solve(&simu->opdr, &recon->RL, rhs);
	    dcellfree(rhs);
	}
	if(parms->tomo.windest){
	    info2("Estimating wind direction and speed using FFT method\n");
	    windest(simu);
	    //Update wind, and interpolation matrix.
	}
	if(parms->tomo.windshift){
	    int factor=parms->tomo.windshift;
	    if(!simu->windshift){
		simu->windshift=spcellnew(recon->npsr, 1);
		for(int ips=0; ips<recon->npsr; ips++){
		    double dispx=simu->dt*simu->atm[ips]->vx*factor;//2 is two cycle delay.
		    double dispy=simu->dt*simu->atm[ips]->vy*factor;
		    info("ips=%d: dispx=%g, dispy=%g\n", ips, dispx, dispy);
		    simu->windshift->p[ips]=mkhb(recon->xloc[ips], recon->xloc[ips], NULL,
						 dispx,dispy,1,0,0);
		}
		spcellwrite(simu->windshift,"windshift");
	    }
	    info2("Using wind information to shift opdr by %d v*dt.\n", factor);
	    for(int ips=0; ips<recon->npsr; ips++){
		dmat *tmp=simu->opdr->p[ips];
		simu->opdr->p[ips]=NULL;
		spmulmat(&simu->opdr->p[ips], simu->windshift->p[ips], tmp, 1);
		dfree(tmp);
	    }
	}
	if(parms->ndm>0){//Do DM fitting
	    dcell *rhs=NULL;
	    muv(&rhs, &recon->FR, simu->opdr, 1);
	    muv_solve(&simu->dmfit_hi, &recon->FL, rhs);
	    dcellfree(rhs);
	}

	dcellcp(&simu->dmerr_hi, simu->dmfit_hi);//keep dmfit_hi for warm restart
    
	/*
	  Forming LGS error signal.
	  2010-01-07: changed dmreal_hi to dmreal to comply with
	  the block diagram. This is before NGS mode removal
	  keep dmfit for warm restart. 
	*/
 
	if(parms->sim.fuseint){
	    if(parms->dbg.psol){
		warning("Using dm for next step to form err signal\n");
		dcelladd(&simu->dmerr_hi, 1., simu->dmint[0], -1);
	    }else{
		dcelladd(&simu->dmerr_hi, 1., simu->dmint[1], -1);
	    }
	}else{
	    /**
	       2010-07-23: Moved remove_dm_ngsmod to after forming error signal. 
	    */
	    if(parms->dbg.psol){
		warning("Using dm for next step to form err signal\n");
		dcelladd(&simu->dmerr_hi, 1., simu->dmint_hi[0], -1);
	    }else{
		dcelladd(&simu->dmerr_hi, 1., simu->dmint_hi[1], -1);
	    }
	}

	if(!parms->sim.fitonly && parms->tomo.split==1){//ahst
	    remove_dm_ngsmod(simu, simu->dmerr_hi);
	}
	if(parms->tomo.ahst_rtt && parms->tomo.split){
	    remove_dm_tt(simu, simu->dmerr_hi);
	}

    }//if high order WFS has output

    if(!parms->sim.fitonly && parms->tomo.split){
	if(parms->tomo.split==2){
	    dcelladd(&simu->opdrmvst, 1, simu->opdr, 1./simu->dtrat_lo);
	}
	//Low order has output
	//2010-12-16: replaces isim+1 by isim.
	if(!parms->sim.closeloop || simu->dtrat_lo==1 || (simu->isim)%simu->dtrat_lo==0){
	    dcellzero(simu->Merr_lo);
	    switch(parms->tomo.split){
	    case 1:{
		NGSMOD_T *ngsmod=recon->ngsmod;
		if(!parms->tomo.ahst_idealngs){//Low order NGS recon.
		    dcellmm(&simu->Merr_lo,ngsmod->Rngs,simu->gradlastcl,"nn",1);
		}//else: there is ideal NGS correction done in perfevl.
	    }
		break;
	    case 2:{
		dcellmm(&simu->gradlastol, recon->GXL, simu->opdrmvst, "nn",-1);
		dcellmm(&simu->Merr_lo, recon->MVRngs, simu->gradlastol, "nn",1);
		if(parms->sim.fuseint){
		    dcelladd(&simu->Merr_lo, 1., simu->dmint[1], -1);
		    error("This mode is not finished\n");
		}else{//form error signal
		    dcelladd(&simu->Merr_lo, 1., simu->Mint_lo[1], -1);
		}
		dcellzero(simu->opdrmvst);
	    }
		break;
	    default:
		error("Invalid parms->tomo.split: %d",parms->tomo.split);
	    }
	}else{
	    dcellfree(simu->Merr_lo);//don't have output.
	}
    }
   
}
/**
   least square reconstructor
*/
void lsr(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int do_hi=(!parms->sim.closeloop || parms->sim.fitonly || 
		     simu->dtrat_hi==1 || (simu->isim)%simu->dtrat_hi==0);
    const int do_low=parms->tomo.split && (!parms->sim.closeloop || simu->dtrat_lo==1
					   || (simu->isim)%simu->dtrat_lo==0);
    dcell *graduse=NULL;
    if(do_hi || do_low){
	if(parms->sim.recon==2){//GLAO
	    graduse=dcellnew(parms->npowfs, 1);
	    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		double scale=1./parms->powfs[ipowfs].nwfs;
		for(int indwfs=0; indwfs<parms->powfs[ipowfs].nwfs; indwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs[indwfs];
		    dadd(&graduse->p[ipowfs], 1, simu->gradlastcl->p[iwfs], scale);
		}
	    }
	}else{
	    graduse=simu->gradlastcl;
	}
    }
    //2010-12-16: replaced isim+1 by isim since recon is delayed from wfsgrad by 1 frame.
    if(do_hi){
	dcell *rhs=NULL;
	muv(&rhs, &(recon->LR), graduse, 1);
	muv_solve(&simu->dmerr_hi,&(recon->LL), rhs);
	dcellfree(rhs);
	if(!parms->sim.fitonly && parms->tomo.split==1){//ahst
	    remove_dm_ngsmod(simu, simu->dmerr_hi);
	}
    }//if high order has output
    if(parms->tomo.split){
	//Low order has output
	//2010-12-16: replaced isim+1 by isim since recon is delayed from wfsgrad by 1 frame.
	if(do_low){
	    dcellzero(simu->Merr_lo);
	    NGSMOD_T *ngsmod=recon->ngsmod;
	    dcellmm(&simu->Merr_lo,ngsmod->Rngs,graduse,"nn",1);
	}else{
	    dcellfree(simu->Merr_lo);//don't have output.
	}
    } 
    if(graduse != simu->gradlastcl){
	dcellfree(graduse);
    }
}
/**
   Deformable mirror control. call tomofit() to do tomo()/fit() or lsr() to do
   least square reconstruction. */
void reconstruct(SIM_T *simu){
    double tk_start=myclockd();
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(simu->gradlastol || simu->gradlastcl){
	switch(parms->sim.recon){//mv
	case 0:
	    tomofit(simu);//tomography and fitting.
	    break;
	case 1:
	case 2:
	    lsr(simu);
	    break;
	}
	if(recon->moao){
	    moao_recon(simu);
	}
	if(parms->sim.mffocus){
	    focus_tracking(simu);
	}
    }
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
    if(parms->sim.recon==0){
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
	if(parms->sim.fuseint){
	    cellarr_dcell(simu->save->dmint, simu->dmint[0]);
	}else{
	    cellarr_dcell(simu->save->dmint_hi, simu->dmint_hi[0]);
	}
	if(simu->save->Merr_lo){
	    cellarr_dcell(simu->save->Merr_lo, simu->Merr_lo);
	}
	if(simu->Mint_lo){
	    cellarr_dcell(simu->save->Mint_lo, simu->Mint_lo[0]);
	}
    }
    simu->tk_recon=myclockd()-tk_start;
}
