/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#define TIMING 0
#include "maos.h"
#include "sim.h"
#include "ahst.h"
/**
   \file perfevl.c Peformance evaluation on science FoV. Notice that the science
   FoV can be different from the DM fitting FoV, which is tuned to better
   sharpen the NGS 
   
   \todo Write a standalone routine that can plot results, using techniques
   developped in drawdaemon.  */
/*
  2009-11-02
  improve threading effiency by devoting 1 thread to do accphi
  2nd thread to do another type of accphi, or doing evaluation.
  lean about futex.

  2009-11-03 the averaging part is not reentrant. the ievl=0
  thread may lag behind ievl=1 thread etc, causing the
  summation to be inaccurat. move it int a separate routine
  that are not threaded.
 */

/**
   Performance evaluation for each direction in parallel mode.  */
static void perfevl_ievl(SIM_T *simu){

    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int isim=simu->isim;
    const int nevl=parms->evl.nevl;
    const int nmod=parms->evl.nmod;
    const int nps=parms->atm.nps;
    const int npsr=parms->atmr.nps;
    const int imoao=parms->evl.moao;
    const double dt=simu->dt;
    const int do_psf=(parms->evl.psfmean || parms->evl.psfhist) && isim>=parms->evl.psfisim;
    dmat *iopdevl=dnew(aper->locs->nloc,1);
    int ievl;
#if TIMING==1
    static double tk[10]={0,0,0,0,0,0,0,0,0,0};
    static int tkcount=0;
#endif
 repeat:
    LOCK(simu->mutex_perfevl);
    ievl=simu->perfevl_ievl++;
    UNLOCK(simu->mutex_perfevl);
    if(ievl<nevl){
#if TIMING==1
	double ck00=myclockd();
#endif
	//Setup pointers for easy usage
	PDMAT(simu->olmp->p[ievl],polmp);//OL mode for each dir
	PDMAT(simu->olep->p[ievl],polep);//OL error for each dir
	PDMAT(simu->clmp->p[ievl],pclmp);
	PDMAT(simu->clep->p[ievl],pclep);

	if(simu->opdevlground){
	    memcpy(iopdevl->p,simu->opdevlground->p,
		   aper->locs->nloc*sizeof(double));
	}else{
	    dzero(iopdevl);
	}
	//Add surfaces along science path. prepared in setup_surf.c
	if(simu->surfevl && simu->surfevl->p[ievl]){
	    dadd(&iopdevl, 1, simu->surfevl->p[ievl], 1);
	}
	//atmosphere contribution.
	if(simu->atm){
	    for(int ips=0; ips<nps; ips++){
		if(ips==simu->perfevl_iground)
		    continue;//already done.
		double hl=simu->atm[ips]->h;
		double scale = 1. - hl/parms->evl.ht;
		double displacex=-simu->atm[ips]->vx*isim*dt
		    +parms->evl.thetax[ievl]*hl;
		double displacey=-simu->atm[ips]->vy*isim*dt
		    +parms->evl.thetay[ievl]*hl;
		prop_grid_stat(simu->atm[ips], aper->locs_stat, 
			       iopdevl->p, 1, displacex, displacey,
			       scale, 1, 0, 0, 1);
	    }
	}
#if TIMING==1
	tk[0]+=myclockd()-ck00; ck00=myclockd(); //Time:1.36s
#endif
	if(parms->save.evlopd){
	    cellarr_dmat(simu->save->evlopdol[ievl],iopdevl);
	}
	if(parms->plot.run){
	    drawopdamp("OL", aper->locs,iopdevl->p , aper->amp1, 
		       "Science Openloop OPD", "x (m)", "y (m)", "OL %d", ievl);
	}
	if(nmod==3){//evaluation piston/tip/tilt removed wve
	    loc_calc_ptt(polep[isim],polmp[isim],
			 aper->locs, aper->ipcc, aper->imcc, 
			 aper->amp, iopdevl->p);
	}else{//more general case
	    loc_calc_mod(polep[isim],polmp[isim],
			aper->mod,aper->amp,iopdevl->p);
	}
	//evaluate time averaged open loop PSF.
	if(parms->evl.psfmean 
	   &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
	      ||(parms->evl.psfol==2 && parms->evl.psf[ievl]))){
	    //Compute on axis OL psf.
	    dmat *opdevlcopy=NULL;
	    dcp(&opdevlcopy,iopdevl);
	    if(parms->evl.psfpttr){
		loc_remove_ptt(opdevlcopy->p,polmp[isim], aper->locs);
	    }
	    ccell *psf2s=psfcomp(opdevlcopy, aper->amp, aper->embed, aper->nembed,
				 parms->evl.psfsize, parms->evl.nwvl, parms->evl.psfwvl);
	    dfree(opdevlcopy);
	    int nwvl=parms->evl.nwvl;
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		cabs22d(&simu->evlpsfolmean->p[iwvl], 1, psf2s->p[iwvl], 1);
	    }
	    ccellfree(psf2s);
	}
	
	if(parms->dbg.evlol)
	    goto repeat;
#if TIMING==1
	tk[1]+=myclockd()-ck00; ck00=myclockd(); 
	//Time:0.706s using dot. 0.32 after improving
#endif
	/*
	  evaluate tomography performance: Apply ideal correction using
	  tomography output directly.
	 */
	if(parms->evl.tomo){
	    dmat *iopdevltomo=NULL;
	    //evaluate tomography error
	    PDMAT(simu->cleptomo->p[ievl],pcleptomo);
	    PDMAT(simu->clmptomo->p[ievl],pclmptomo);
	    dcp(&iopdevltomo, iopdevl);//make a copy
	    if(simu->opdr){
		for(int ipsr=0; ipsr<npsr; ipsr++){
		    double hl=parms->atmr.ht[ipsr];
		    double scale = 1. - hl/parms->evl.ht;
		    double displacex=parms->evl.thetax[ievl]*hl;
		    double displacey=parms->evl.thetay[ievl]*hl;
		    prop_nongrid(recon->xloc[ipsr], 
				 simu->opdr->p[ipsr]->p,
				 aper->locs, NULL, iopdevltomo->p, -1,
				 displacex, displacey, scale, 0, 0);
		}
	    }
	    if(parms->plot.run){
		drawopdamp("evlTomo", aper->locs, iopdevltomo->p, aper->amp1,
			   "Science Ideal MOAO Correction OPD","x (m)", "y (m)",
			   "Tomo %d",ievl);
	    }
	    
	    if(nmod==3){
		loc_calc_ptt(pcleptomo[isim],pclmptomo[isim],
			     aper->locs, aper->ipcc, 
			     aper->imcc, aper->amp, iopdevltomo->p);
	    }else{
		loc_calc_mod(pcleptomo[isim],pclmptomo[isim],
			    aper->mod,aper->amp,iopdevltomo->p);
	    }
	    //Evaluate closed loop PSF time history and time average
	    if(do_psf && parms->evl.psf[ievl] ){
		if(parms->evl.psfpttr){
		    if(isim==parms->sim.start && ievl==0){
			warning("Removing tip/tilt from PSF\n");
		    }
		    loc_remove_ptt(iopdevltomo->p,pclmp[isim], aper->locs);
		}
		ccell *psf2s=psfcomp(iopdevltomo, aper->amp, aper->embed, aper->nembed,
				     parms->evl.psfsize, parms->evl.nwvl, parms->evl.psfwvl);
		int nwvl=parms->evl.nwvl;
		if(parms->evl.psfmean){
		    PDCELL(simu->evlpsftomomean, pevlpsftomomean);
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			cabs22d(&pevlpsftomomean[ievl][iwvl], 1, psf2s->p[iwvl], 1);
		    }
		}
		if(parms->evl.psfhist){
		    cellarr_ccell(simu->evlpsftomohist[ievl], psf2s);
		}
		ccellfree(psf2s);
	    }
	    dfree(iopdevltomo);
	}
#if TIMING==1
	tk[2]+=myclockd()-ck00; ck00=myclockd(); //Time:0
#endif
	if(parms->evl.tomo==2){
	    info("Skip DM performance evalution\n");
	    goto repeat;
	}
	//Apply dm correction. tip/tilt command is contained in DM commands
	if(simu->dmreal){
	    int ndm=parms->ndm;
	    for(int idm=0; idm<ndm; idm++){
		double ht=parms->dm[idm].ht;
		double scale = 1. - ht/parms->evl.ht;
		double displacex=parms->evl.thetax[ievl]*ht;
		double displacey=parms->evl.thetay[ievl]*ht;
		if(simu->cachedm){
		    int iscale=parms->evl.scalegroup[idm];
		    prop_grid_stat(&simu->cachedm[idm][iscale], aper->locs_stat, iopdevl->p, -1, 
				   displacex, displacey, scale, 0, 0, 0, 1);

		}else{
		    if(parms->dm[idm].cubic){
			prop_nongrid_cubic(recon->aloc[idm], simu->dmreal->p[idm]->p,
					   aper->locs, NULL, iopdevl->p,-1,
					   displacex, displacey,scale,
					   parms->dm[idm].iac,0,0);
		    }else{
			prop_nongrid(recon->aloc[idm], simu->dmreal->p[idm]->p,
				     aper->locs, NULL, iopdevl->p,-1,
				     displacex, displacey,scale,0,0);
		    }
		}
	    }
#if TIMING==1
	    tk[3]+=myclockd()-ck00; ck00=myclockd(); 
	    //Time:3.5 if not cached, 0.34 if cached.
#endif
	   
	}
	if(imoao>-1){
	    PDCELL(simu->moao_evl, dmevl);
	    if(dmevl[ievl][simu->moao_evl->nx-1]){
		//TIC;tic;
		//prop is faster than spmulvec
		if(parms->moao[imoao].cubic){
		    prop_nongrid_cubic(recon->moao[imoao].aloc,dmevl[ievl][simu->moao_evl->nx-1]->p,
				       aper->locs, NULL, iopdevl->p, -1, 0,0,1,
				       parms->moao[imoao].iac, 
				       0,0);
		}else{
		    prop_nongrid(recon->moao[imoao].aloc,dmevl[ievl][simu->moao_evl->nx-1]->p,
				 aper->locs, NULL, iopdevl->p, -1, 0,0,1,0,0);
		}
	    }
	}
	if(parms->plot.run){
	    drawopdamp("CL", aper->locs, iopdevl->p, aper->amp1,
		       "Science Closeloop OPD", "x (m)", "y (m)",
		       "CL %d",ievl);
	}
	if(parms->save.evlopd){
	    cellarr_dmat(simu->save->evlopdcl[ievl],iopdevl);
	}
	//Evaluate closed loop performance.
	if(parms->tomo.split){//for split tomography
	    if(ievl==parms->evl.indoa || parms->dbg.clemp_all){
		dcp(&simu->opdevl->p[ievl], iopdevl);
	    }
	    if(parms->ndm<=2){
		PDMAT(simu->cleNGSmp->p[ievl], pcleNGSmp);
		//compute the dot product of wavefront with NGS mode for that direction
		if(nmod==3){
		    calc_ngsmod_dot(pclep[isim],pclmp[isim],
				    pcleNGSmp[isim],parms,recon,aper,
				    iopdevl->p,ievl);
		}else{
		    calc_ngsmod_dot(NULL,NULL,
				    pcleNGSmp[isim],parms,recon,aper,
				    iopdevl->p,ievl);
		    loc_calc_mod(pclep[isim],pclmp[isim],
				aper->mod,aper->amp,iopdevl->p);
		}
	    }else{
		error("Not implemented\n");
	    }
	}else{//for integrated tomography.
	    if(nmod==3){
		loc_calc_ptt(pclep[isim],pclmp[isim],
			 aper->locs, aper->ipcc, aper->imcc, 
			 aper->amp, iopdevl->p);
	    }else{
		loc_calc_mod(pclep[isim],pclmp[isim],
			    aper->mod,aper->amp,iopdevl->p);
	    }
	}
	//Evaluate closed loop PSF.
	if(do_psf && parms->evl.psf[ievl]){
	    //warning("Output PSF for direction %d\n",ievl);
	    if(parms->evl.psfpttr){
		if(isim==parms->sim.start && ievl==0){
		    warning("Removing tip/tilt from PSF\n");
		}
		loc_remove_ptt(iopdevl->p,pclmp[isim], aper->locs);
	    }
	    ccell *psf2s=psfcomp(iopdevl, aper->amp, aper->embed, aper->nembed,
				 parms->evl.psfsize, parms->evl.nwvl, parms->evl.psfwvl);
	    int nwvl=parms->evl.nwvl;
	    if(parms->evl.psfmean){
		PDCELL(simu->evlpsfmean, pevlpsfmean);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    cabs22d(&pevlpsfmean[ievl][iwvl], 1, psf2s->p[iwvl], 1);
		}
	    }
	    if(parms->evl.psfhist){
		cellarr_ccell(simu->evlpsfhist[ievl], psf2s);
	    }
	    ccellfree(psf2s);
	}
#if TIMING==1
	tk[4]+=myclockd()-ck00; ck00=myclockd(); //Time:0.586s
#endif
	goto repeat;
    }
    dfree(iopdevl);
#if TIMING==1
    tkcount++;
    double tkscale=1./(double)tkcount;
    info2("Timing:%.4f %.4f %.4f %.4f %.4f\n",
	  tk[0]*tkscale, tk[1]*tkscale, tk[2]*tkscale,
	  tk[3]*tkscale, tk[4]*tkscale);
#endif
}
/**
   Evaluation field averaged performance.
*/
static void perfevl_mean(SIM_T *simu){

    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    const int nmod=parms->evl.nmod;
    const int nevl=parms->evl.nevl;

    for(int imod=0; imod<nmod; imod++){
	int ind=imod+nmod*isim;
	simu->ole->p[ind]=0;
	simu->cle->p[ind]=0;
	for(int ievl=0; ievl<nevl; ievl++){
	    double wt=parms->evl.wt[ievl];
	    simu->ole->p[ind]+=wt*simu->olep->p[ievl]->p[ind];
	    simu->cle->p[ind]+=wt*simu->clep->p[ievl]->p[ind];
	}
	
    }
  
  
    if(parms->evl.tomo && simu->opdr){
	for(int imod=0; imod<nmod; imod++){
	    int ind=imod+nmod*isim;
	    simu->cletomo->p[ind]=0;
	    for(int ievl=0; ievl<nevl; ievl++){
		double wt=parms->evl.wt[ievl];
		simu->cletomo->p[ind]+=wt*simu->cleptomo->p[ievl]->p[ind];
	    }
	}
    }
    const RECON_T *recon=simu->recon;
    if(parms->tomo.split){
	if(parms->ndm<=2){
	    /*
	       convert cleNGSm into mode and put NGS mode WVE into clem.
	    */
	    int nngsmod=recon->ngsmod->nmod;
	    {
		if(simu->corrNGSm && simu->Mint_lo[0] && isim<parms->sim.end-1){
		    double *pcorrNGSm=simu->corrNGSm->p+(isim+1)*nngsmod;
		    for(int imod=0; imod<nngsmod; imod++){
			pcorrNGSm[imod]=simu->Mint_lo[0]->p[0]->p[imod];
		    }
		}
	    }
	    double *pcleNGSm=simu->cleNGSm->p+isim*nngsmod;
	    for(int imod=0; imod<nngsmod; imod++){
		pcleNGSm[imod]=0;
	    }
	    for(int ievl=0; ievl<nevl; ievl++){
		double wt=parms->evl.wt[ievl];
		double *pcleNGSmp=simu->cleNGSmp->p[ievl]->p+isim*nngsmod;
		for(int imod=0; imod<nngsmod; imod++){
		    pcleNGSm[imod]+=pcleNGSmp[imod]*wt;
		}
	    }
	    double tt=dwdot2(pcleNGSm, recon->ngsmod->IMCC_TT, pcleNGSm);
	    double ngs=dwdot(pcleNGSm, recon->ngsmod->IMCC, pcleNGSm);
	    double tot=simu->cle->p[isim*nmod];
	    double modngs[nngsmod];
	    memset(modngs,0,sizeof(double)*nngsmod);
	    //turn cleNGSm to modes
	    dmulvec(modngs,recon->ngsmod->IMCC,pcleNGSm,1);
	    memcpy(pcleNGSm,modngs,sizeof(double)*nngsmod);
	    dcell *Mngs=dcellnew(1,1);
	    Mngs->p[0]=dnew_ref(pcleNGSm,nngsmod,1);//ref the data
	    if(simu->parms->tomo.ahst_idealngs){
		/*
		  apply ideal ngs modes immediately to dmreal.
		  Don't forget to updated DM Cache.
		*/
		ngsmod2dm(&simu->dmreal,simu->recon, Mngs, 1.);
		calc_cachedm(simu);
		tot-=ngs; ngs=0; tt=0;
	    }
	    simu->clem->p[isim*3]=tot-ngs; //lgs mode
	    simu->clem->p[isim*3+1]=tt;    //tt mode
	    simu->clem->p[isim*3+2]=ngs;   //ngs mod
	    simu->status->clerrlo=sqrt(ngs)*1e9;
	    simu->status->clerrhi=sqrt(tot-ngs)*1e9;

	    const APER_T *aper=simu->aper;
	    //compute on axis error spliting to tip/tilt and high order.
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!parms->dbg.clemp_all && ievl!=parms->evl.indoa)
		    continue;
		ngsmod2science(simu->opdevl->p[ievl],parms,simu->recon,aper,Mngs,ievl,-1);
		double pttr[3]={0,0,0};
		loc_calc_ptt(pttr,NULL,aper->locs,aper->ipcc,aper->imcc,aper->amp,
			     simu->opdevl->p[ievl]->p);
		simu->clemp->p[ievl]->p[isim*3]=pttr[0];//LGS mode
		simu->clemp->p[ievl]->p[isim*3+1]=tt;//TT mode
		simu->clemp->p[ievl]->p[isim*3+2]=simu->clep->p[ievl]->p[nmod*isim]-pttr[0];//PR-LGS
		dfree(simu->opdevl->p[ievl]);
	    }
	    dcellfree(Mngs);//data is kept
	}else{
	    error("Not implemented\n");
	}
    }else{//if not split
	simu->status->clerrhi=sqrt(simu->cle->p[nmod*isim]*1e18);
	simu->status->clerrlo=sqrt(simu->cle->p[nmod*isim+1]*1e18);
    }//if split
    
    if(parms->dbg.noatm==0 && simu->cle->p[nmod*isim] > simu->ole->p[nmod*isim]*10){
	warning("The loop is diverging: OL: %g CL: %g\n",  
	      simu->ole->p[nmod*isim],  simu->cle->p[nmod*isim]);
    }
}
/**
   Evaluate performance by calling perfevl_ievl in parallel and then calls
perfevl_mean to field average.  */
void perfevl(SIM_T *simu){
  
    //Cache the ground layer.
    int ips=simu->perfevl_iground;
    if(ips!=-1 && simu->atm){
	simu->opdevlground=dnew(simu->aper->locs->nloc,1);
	const int isim=simu->isim;
	const double dt=simu->dt;
	double displacex=-simu->atm[ips]->vx*isim*dt;
	double displacey=-simu->atm[ips]->vy*isim*dt;
	prop_grid_stat(simu->atm[ips], simu->aper->locs_stat, 
		       simu->opdevlground->p,1,
		       displacex, displacey,
		       1., 1, 0, 0, 1);
    }
    simu->perfevl_ievl=0;
    CALL(perfevl_ievl,simu,simu->nthread);//no space allowd.
    dfree(simu->opdevlground);simu->opdevlground=NULL;
    perfevl_mean(simu);
}
