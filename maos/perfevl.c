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
#include "sim_utils.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
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
#define TIMING 0
#if TIMING == 1
#define TIM(A) double tk##A=myclockd()
#else
#define TIM(A)
#endif
static double *opdzlim=NULL;
static void perfevl_ideal_atm(SIM_T *simu, dmat *iopdevl, int ievl, double alpha){
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const double hs = parms->evl.hs[ievl];
  
    for(int idm=0; idm<parms->ndm; idm++){
	const double ht = parms->dm[idm].ht+parms->dm[idm].vmisreg;
	double dispx=ht*parms->evl.thetax[ievl]+parms->evl.misreg[0];
	double dispy=ht*parms->evl.thetay[ievl]+parms->evl.misreg[1];
	double scale=1.-ht/hs;
	if(parms->dm[idm].cubic){
	    prop_grid_cubic(simu->dmprojsq[idm],
			    aper->locs, aper->amp->p, iopdevl->p, 
			    alpha, dispx, dispy, scale, parms->dm[idm].iac, 
			    0, 0);
	}else{
	    prop_grid(simu->dmprojsq[idm],
		      aper->locs, aper->amp->p, iopdevl->p,
		      alpha, dispx, dispy, scale, 0,
		      0, 0);
	}
    }
}
static void perfevl_psfcl(const PARMS_T *parms, const APER_T *aper,
			  dcell *evlpsfmean, cellarr** evlpsfhist,
			  dmat *iopdevl, int ievl){
    /* the OPD after this time will be tilt removed. Don't use for performance
       evaluation. */
    ccell *psf2s=psfcomp(iopdevl, aper->amp->p, aper->embed, aper->nembed,
			 parms->evl.psfsize, parms->evl.nwvl, parms->evl.wvl);
    int nwvl=parms->evl.nwvl;
    if(parms->evl.psfmean){
	PDCELL(evlpsfmean, pevlpsfmean);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cabs22d(&pevlpsfmean[ievl][iwvl], 1, psf2s->p[iwvl], 1);
	}
    }
    if(parms->evl.psfhist){
	cellarr_ccell(evlpsfhist[ievl], psf2s);
    }
    if(parms->plot.run){
	dmat *psftemp=NULL;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cabs22d(&psftemp, 1, psf2s->p[iwvl], 1);
	    dcwlog10(psftemp);
	    double xylim[4]={-12,12,-12,12};
	    ddraw("CL PSF", psftemp, xylim, opdzlim, 
		  "Science Closed Loop PSF", 
		  "x", "y", "CL%2d PSF %.2f", ievl, parms->evl.wvl[iwvl]*1e6);
	    dfree(psftemp);
	}
    }
    ccellfree(psf2s);
}
/**
   Performance evaluation for each direction in parallel mode.  */
void perfevl_ievl(thread_t *info){
    SIM_T *simu=info->data;
    const int ievl=info->start;
    assert(info->end==info->start+1);/*only one evl. */
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int isim=simu->isim;
    const int nmod=parms->evl.nmod;
    const int nps=parms->atm.nps;
    const int npsr=parms->atmr.nps;
    const int imoao=parms->evl.moao;
    const double dt=simu->dt;
    const int nthread=parms->evl.nthread;
    const int do_psf_cov=(parms->evl.psfmean || parms->evl.psfhist || parms->evl.opdcov) && isim>=parms->evl.psfisim && parms->evl.psf[ievl];
    const int save_evlopd=parms->save.evlopd>0 && ((isim+1)%parms->save.evlopd)==0;
    dmat *iopdevl=dnew(aper->locs->nloc,1);
    TIM(0);
    /*Setup pointers for easy usage */
    PDMAT(simu->olmp->p[ievl],polmp);/*OL mode for each dir */
    PDMAT(simu->olep->p[ievl],polep);/*OL error for each dir */
    PDMAT(simu->clmp->p[ievl],pclmp);
    PDMAT(simu->clep->p[ievl],pclep);

    /*atmosphere contribution. */
    if(parms->sim.idealevl){
	perfevl_ideal_atm(simu, iopdevl, ievl, 1);
    }else if(simu->atm && !parms->sim.wfsalias){
	if(simu->opdevlground){
	    memcpy(iopdevl->p,simu->opdevlground->p, aper->locs->nloc*sizeof(double));
	}else{
	    dzero(iopdevl);
	}
	/*fix me: the ray tracing of the same part must be performed in the same thread. */
	for(int ips=0; ips<nps; ips++){
	    if(ips!=simu->perfevl_iground || !simu->opdevlground){
		int ind=ievl+parms->evl.nevl*ips;
		simu->evl_propdata_atm[ind].phiout=iopdevl->p;
		simu->evl_propdata_atm[ind].displacex1=-simu->atm[ips]->vx*isim*dt;
		simu->evl_propdata_atm[ind].displacey1=-simu->atm[ips]->vy*isim*dt;
		CALL_THREAD(simu->evl_prop_atm[ind], nthread, 0);
	    }
	}
    }
    if(simu->telws){/*Wind shake */
	double tmp=simu->telws->p[isim];
	double angle=simu->winddir?simu->winddir->p[0]:0;
	double ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
	loc_add_ptt(iopdevl->p, ptt, aper->locs);
    }
    /*Add surfaces along science path. prepared in setup_surf.c */
    if(aper->opdadd && aper->opdadd->p[ievl]){
	dadd(&iopdevl, 1, aper->opdadd->p[ievl], 1);
    }

    TIM(1);
    if(save_evlopd){
	cellarr_dmat(simu->save->evlopdol[ievl],iopdevl);
    }
    if(parms->plot.run){
	drawopdamp("OL", aper->locs,iopdevl->p , aper->amp1->p, opdzlim,
		   "Science Open Loop OPD", "x (m)", "y (m)", "OL %d", ievl);
    }
    if(nmod==3){/*evaluation piston/tip/tilt removed wve */
	loc_calc_ptt(polep[isim],polmp[isim],
		     aper->locs, aper->ipcc, aper->imcc, 
		     aper->amp->p, iopdevl->p);
    }else{/*more general case */
	loc_calc_mod(polep[isim],polmp[isim],
		     aper->mod,aper->amp->p,iopdevl->p);
    }

    /*evaluate time averaged open loop PSF. */
    if(parms->evl.psfmean 
       && isim>=parms->evl.psfisim 
       &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
	  ||(parms->evl.psfol==2 && parms->evl.psf[ievl]))){
	/*Compute on axis OL psf. */
	dmat *opdevlcopy=NULL;
	if(parms->evl.psfpttr[ievl]){
	    dcp(&opdevlcopy,iopdevl);
	    loc_remove_ptt(opdevlcopy->p,polmp[isim], aper->locs);
	}else{
	    opdevlcopy=dref(iopdevl);
	}
	ccell *psf2s=psfcomp(opdevlcopy, aper->amp->p, aper->embed, aper->nembed,
			     parms->evl.psfsize, parms->evl.nwvl, parms->evl.wvl);
	dfree(opdevlcopy);
	int nwvl=parms->evl.nwvl;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cabs22d(&simu->evlpsfolmean->p[iwvl], 1, psf2s->p[iwvl], 1);
	}
	if(parms->plot.run){
	    dmat *psftemp=NULL;
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		cabs22d(&psftemp, 1, psf2s->p[iwvl], 1);
		ddraw("OL PSF", psftemp, NULL, NULL, "Science Openloop PSF", 
		      "x", "y", "OL%2d PSF %.2f", ievl,  parms->evl.wvl[iwvl]*1e6);
		dfree(psftemp);
	    }
	}
	ccellfree(psf2s);
    }
    if(parms->sim.evlol) goto end;
    TIM(2);
    if(parms->evl.tomo){
	/*
	  evaluate tomography performance: Apply ideal correction using
	  tomography output directly.
	*/
	if(simu->opdr){
	    for(int ipsr=0; ipsr<npsr; ipsr++){
		double hl=parms->atmr.ht[ipsr];
		double scale = 1. - hl/parms->evl.hs[ievl];
		double displacex=parms->evl.thetax[ievl]*hl+parms->evl.misreg[0];
		double displacey=parms->evl.thetay[ievl]*hl+parms->evl.misreg[1];
		prop_nongrid(recon->xloc[ipsr], 
			     simu->opdr->p[ipsr]->p,
			     aper->locs, NULL, iopdevl->p, -1,
			     displacex, displacey, scale, 0, 0);
	    }
	}
    }else{
	/* Apply dm correction. tip/tilt command is contained in DM commands */
	if(simu->dmreal){
	    int ndm=parms->ndm;
	    for(int idm=0; idm<ndm; idm++){
		int ind=ievl+parms->evl.nevl*idm;
		simu->evl_propdata_dm[ind].phiout=iopdevl->p;
		CALL_THREAD(simu->evl_prop_dm[ind], nthread, 0);
	    }
	}
	
	TIM(4);
	if(imoao>-1){
	    dmat **dmevl=simu->moao_evl->p;
	    if(dmevl[ievl]){
		/**
		   prop is faster than spmulvec. \fixme check definition of misreg
		*/
		if(parms->moao[imoao].cubic){
		    prop_nongrid_cubic(recon->moao[imoao].aloc,dmevl[ievl]->p,
				       aper->locs, NULL, iopdevl->p, -1, 0,0,1,
				       parms->moao[imoao].iac, 
				       0,0);
		}else{
		    prop_nongrid(recon->moao[imoao].aloc,dmevl[ievl]->p,
				 aper->locs, NULL, iopdevl->p, -1, 0,0,1,0,0);
		}
	    }
	}
    }
    if(parms->plot.run){
	drawopdamp("CL", aper->locs, iopdevl->p, aper->amp1->p,NULL,
		   "Science Closed loop OPD", "x (m)", "y (m)",
		   "CL %d",ievl);
    }
    if(save_evlopd){
	cellarr_dmat(simu->save->evlopdcl[ievl],iopdevl);
    }

    /*Evaluate closed loop performance. */
    if(parms->recon.split){/*for split tomography */
	if(parms->ndm<=2){
	    PDMAT(simu->cleNGSmp->p[ievl], pcleNGSmp);
	    /*compute the dot product of wavefront with NGS mode for that direction */
	    if(nmod==3){
		calc_ngsmod_dot(pclep[isim],pclmp[isim],
				pcleNGSmp[isim],parms,recon,aper,
				iopdevl->p,ievl);
	    }else{/*since more modes are wanted. don't use ngsmod split. */
		calc_ngsmod_dot(NULL,NULL,
				pcleNGSmp[isim],parms,recon,aper,
				iopdevl->p,ievl);
		loc_calc_mod(pclep[isim],pclmp[isim],
			     aper->mod,aper->amp->p,iopdevl->p);
	    }
	}else{
	    error("Not implemented\n");
	}
    }else{/*for integrated tomography. */
	if(nmod==3){
	    loc_calc_ptt(pclep[isim],pclmp[isim],
			 aper->locs, aper->ipcc, aper->imcc, 
			 aper->amp->p, iopdevl->p);
	}else{
	    loc_calc_mod(pclep[isim],pclmp[isim],
			 aper->mod,aper->amp->p,iopdevl->p);
	}
    }
    if(do_psf_cov){
	if(parms->evl.psfngsr[ievl]!=0){
	    /* even if psfpttr=1, referencing is ok.  Change to copy if
	       incompatible in the future.*/
	    simu->evlopd->p[ievl]=dref(iopdevl);
	}
	if(parms->evl.psfngsr[ievl]!=2){/*ngsr==2 means only want ngsr. */
	    /** opdcov does not have p/t/t removed. do it in postproc is necessary*/
	    if(parms->evl.psfpttr[ievl]){
		if(isim==parms->evl.psfisim && ievl==0){
		    warning("Removing piston/tip/tilt from OPD.\n");
		}
		loc_remove_ptt(iopdevl->p, pclmp[isim], aper->locs);
	    }else if(parms->evl.opdcov){/*remove piston */
		dadds(iopdevl, -pclmp[isim][0]);
	    }
	    if(parms->evl.opdcov){
		dmm(&simu->evlopdcov->p[ievl], iopdevl, iopdevl, "nt", 1);
		dadd(&simu->evlopdmean->p[ievl], 1, iopdevl, 1);
	    }/*opdcov */
	    if(parms->evl.psfmean || parms->evl.psfhist){/*Evaluate closed loop PSF.	 */
		perfevl_psfcl(parms, aper, simu->evlpsfmean, simu->save->evlpsfhist, iopdevl, ievl);
	    }/*do_psf */
	}
    }
    TIM(5);
 end:
#if TIMING==1
    info2("Evl %d timing:ray atm %.4f evlol %.4f ray dm %.4f evlcl %.4f\n",
	  ievl, tk1-tk0, tk2-tk1, tk4-tk2, tk5-tk4);
#endif
    dfree(iopdevl);
}
/**
   Evaluation field averaged performance.
*/
static void perfevl_mean(SIM_T *simu){

    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    const int nmod=parms->evl.nmod;
    const int nevl=parms->evl.nevl;
    /*Field average the OL error */
    for(int imod=0; imod<nmod; imod++){
	int ind=imod+nmod*isim;
	simu->ole->p[ind]=0;
	for(int ievl=0; ievl<nevl; ievl++){
	    double wt=parms->evl.wt[ievl];
	    simu->ole->p[ind]+=wt*simu->olep->p[ievl]->p[ind];
	}
    }
    if(parms->sim.evlol)
	return;

    /*Field average the CL error */
    for(int imod=0; imod<nmod; imod++){
	int ind=imod+nmod*isim;
	simu->cle->p[ind]=0;
	for(int ievl=0; ievl<nevl; ievl++){
	    double wt=parms->evl.wt[ievl];
	    simu->cle->p[ind]+=wt*simu->clep->p[ievl]->p[ind];
	}
    }
  
    const RECON_T *recon=simu->recon;
    if(parms->recon.split){
	if(parms->ndm<=2){
	    /* convert cleNGSm into mode and put NGS mode WVE into clem. */
	    int nngsmod=recon->ngsmod->nmod;
	    
	    if(simu->corrNGSm && simu->Mint_lo[0] && isim<parms->sim.end-1){
		double *pcorrNGSm=simu->corrNGSm->p+(isim+1)*nngsmod;
		for(int imod=0; imod<nngsmod; imod++){
		    pcorrNGSm[imod]=simu->Mint_lo[0]->p[0]->p[imod];
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
	    /*turn cleNGSm to modes */
	    dmulvec(modngs,recon->ngsmod->IMCC,pcleNGSm,1);
	    memcpy(pcleNGSm,modngs,sizeof(double)*nngsmod);
	    dcell *Mngs=dcellnew(1,1);
	    Mngs->p[0]=dnew_ref(nngsmod,1,pcleNGSm);/*ref the data */
	    if(simu->parms->tomo.ahst_idealngs){
		/* apply ideal ngs modes immediately to dmreal.  Don't forget to
		  updated DM Cache. */
		ngsmod2dm(&simu->dmreal,simu->recon, Mngs, 1.);
		calc_cachedm(simu);
		tot-=ngs; ngs=0; tt=0;
	    }
	    int do_psf=(parms->evl.psfmean || parms->evl.psfhist);
	    if(isim>=parms->evl.psfisim && (do_psf || parms->evl.opdcov)){
#if USE_CUDA
		if(parms->gpu.evl){
		    gpu_perfevl_ngsr(simu, pcleNGSm);
		}else{
#endif
		    const APER_T *aper=simu->aper;
		    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			if(!parms->evl.psf[ievl] || !parms->evl.psfngsr[ievl]) continue;
			dmat *iopdevl=simu->evlopd->p[ievl];
			if(!iopdevl) continue;
			ngsmod2science(iopdevl, parms, recon, aper, pcleNGSm, ievl, -1);
			if(parms->evl.psfpttr[ievl]){
			    /*we cannot use clmp because the removed ngsmod
			      has tip/tilt component*/
			    double ptt[3];
			    loc_calc_ptt(NULL, ptt, aper->locs, aper->ipcc, aper->imcc, aper->amp->p, iopdevl->p);
			    loc_remove_ptt(iopdevl->p, ptt, aper->locs);
			}
			if(parms->evl.opdcov){
			    dmm(&simu->evlopdcov_ngsr->p[ievl], iopdevl, iopdevl, "nt", 1);
			    dadd(&simu->evlopdmean_ngsr->p[ievl], 1, iopdevl, 1);
			}
			if(do_psf){
			    perfevl_psfcl(parms, aper, simu->evlpsfmean_ngsr, simu->save->evlpsfhist_ngsr, iopdevl, ievl);
			}
			dfree(simu->evlopd->p[ievl]);
		    }
#if USE_CUDA
		}
#endif
	    }/*ideal ngs */
	    simu->clem->p[isim*3]=tot-ngs; /*lgs mode */
	    simu->clem->p[isim*3+1]=tt;    /*tt mode */
	    simu->clem->p[isim*3+2]=ngs;   /*ngs mod */
	    simu->status->clerrlo=sqrt(ngs)*1e9;
	    simu->status->clerrhi=sqrt(tot-ngs)*1e9;

	    /*compute error spliting to tip/tilt and high order for any
	      direction.*/
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		/*New calculation. Notice that the calculations are for
		  (phi-Hm*m)'*W*(phi-Hm*m) where m is cleNGSm, phi is opdevl
		  and W is aperture weighting. The calculation is equivalent
		  to phi'*W*phi-2*phi'*W*(Hm*m)+(Hm*m)'*W*(Hm*m).
		  phi'*W*phi is the clep phi'*W*Hm is the same as
		  pcleNGSmp'; (Hm*m)'*W*(Hm*m) = m'*(Hm'*W*Hm)*m. (Hm'*W*Hm)
		  is simply MCCp.
		*/
		if(!recon->ngsmod->MCCP->p[ievl]) continue;
 		double *pcleNGSmp=simu->cleNGSmp->p[ievl]->p+isim*nngsmod;
		double sum=0;
		for(int imod=0; imod<nngsmod; imod++){
		    sum+=pcleNGSmp[imod]*pcleNGSm[imod];
		}
		double sum2=dwdot(pcleNGSm, recon->ngsmod->MCCP->p[ievl], pcleNGSm);
		double tot2=simu->clep->p[ievl]->p[isim*nmod]-2.*sum+sum2;
		simu->clemp->p[ievl]->p[isim*3]=tot2;/*LGS mode */
		simu->clemp->p[ievl]->p[isim*3+1]=tt;/*TT mode */
		simu->clemp->p[ievl]->p[isim*3+2]=simu->clep->p[ievl]->p[nmod*isim]-tot2;/*PR-LGS */
	    }
	    dcellfree(Mngs);/*data is kept */
	}else{
	    error("Not implemented\n");
	}
    }else{/*if not split */
	simu->status->clerrhi=sqrt(simu->cle->p[nmod*isim]*1e18);
	simu->status->clerrlo=sqrt(simu->cle->p[nmod*isim+1]*1e18);
    }/*if split */
    
    if(parms->sim.noatm==0 && simu->cle->p[nmod*isim] > MAX(simu->ole->p[nmod*isim]*100, 1e-13)){
	warning("The loop is diverging: OL: %g CL: %g\n",  
		simu->ole->p[nmod*isim],  simu->cle->p[nmod*isim]);
    }
}
/**
   Save telemetry
*/
static void perfevl_save(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const int isim=simu->isim;
    if(parms->evl.psfmean && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.psfmean)){
	info2("Step %d: Output PSF\n", isim);
	dcell *psfmean=NULL;
	double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	dcelladd(&psfmean, 0, simu->evlpsfmean, scale);
	PDCELL(psfmean, pcl);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!simu->save->evlpsfmean[ievl]) continue;
	    for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
		pcl[ievl][iwvl]->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		cellarr_dmat(simu->save->evlpsfmean[ievl], pcl[ievl][iwvl]);
		free(pcl[ievl][iwvl]->header);pcl[ievl][iwvl]->header=NULL;
	    }
	}
	dcelladd(&psfmean, 0, simu->evlpsfmean_ngsr, scale);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!simu->save->evlpsfmean_ngsr[ievl]) continue;
	    for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
		pcl[ievl][iwvl]->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		cellarr_dmat(simu->save->evlpsfmean_ngsr[ievl], pcl[ievl][iwvl]);
		free(pcl[ievl][iwvl]->header);pcl[ievl][iwvl]->header=NULL;
	    }
	}
	if(simu->evlpsfolmean){
	    dcellfree(psfmean);/*different format. */
	    scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    if(parms->evl.psfol==2){
		scale=scale/parms->evl.npsf;
	    }
	    dcelladd(&psfmean, 0, simu->evlpsfolmean, scale);
	    for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
		psfmean->p[iwvl]->header=evl_header(simu->parms, simu->aper, -1, iwvl);
		cellarr_dmat(simu->save->evlpsfolmean, psfmean->p[iwvl]);
		free(psfmean->p[iwvl]->header); psfmean->p[iwvl]->header=NULL;
	    }
	}
	dcellfree(psfmean);
    }
    if(parms->evl.opdcov && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.opdcov)){
	info2("Step %d: Output opdcov\n", isim);
	long nstep=isim+1-parms->evl.psfisim;
	double scale=1./nstep;
	dmat *covmean=NULL;
	dmat *opdmean=NULL;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!simu->evlopdcov->p[ievl]) continue;
	    dadd(&covmean, 0, simu->evlopdcov->p[ievl], scale);
	    cellarr_dmat(simu->save->evlopdcov[ievl], covmean);
	    dadd(&opdmean, 0, simu->evlopdmean->p[ievl], scale);
	    cellarr_dmat(simu->save->evlopdmean[ievl], opdmean);
	}
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!simu->evlopdcov_ngsr->p[ievl]) continue;
	    dadd(&covmean, 0, simu->evlopdcov_ngsr->p[ievl], scale);
	    cellarr_dmat(simu->save->evlopdcov_ngsr[ievl], covmean);
	    dadd(&opdmean, 0, simu->evlopdmean_ngsr->p[ievl], scale);
	    cellarr_dmat(simu->save->evlopdmean_ngsr[ievl], opdmean);
	}
	dfree(opdmean);
	dfree(covmean);
    }
}
/**
   Evaluate performance by calling perfevl_ievl in parallel and then calls
   perfevl_mean to field average.  */
void perfevl(SIM_T *simu){
    double tk_start=myclockd();
    /*Cache the ground layer. */
    const PARMS_T *parms=simu->parms;
    if(!(parms->gpu.evl)){
	int ips=simu->perfevl_iground;
	if(ips!=-1 && simu->atm && !parms->sim.idealevl){
	    simu->opdevlground=dnew(simu->aper->locs->nloc,1);
	    const int ievl=0;/*doesn't matter for ground layer. */
	    int ind=ievl+parms->evl.nevl*ips;
	    const int isim=simu->isim;
	    const double dt=simu->dt;
	    simu->evl_propdata_atm[ind].phiout=simu->opdevlground->p;
	    simu->evl_propdata_atm[ind].displacex1=-simu->atm[ips]->vx*isim*dt;
	    simu->evl_propdata_atm[ind].displacey1=-simu->atm[ips]->vy*isim*dt;
	    CALL_THREAD(simu->evl_prop_atm[ind], parms->evl.nthread, 0);
	}
    }
    CALL_THREAD(simu->perf_evl, parms->evl.nevl, 0);
    if(!(parms->gpu.evl)){
	dfree(simu->opdevlground);simu->opdevlground=NULL;
    }
    perfevl_mean(simu);
#if USE_CUDA
    if(parms->gpu.evl && parms->gpu.psf){
	gpu_perfevl_save(simu);
    }else
#endif
	perfevl_save(simu);

    simu->tk_eval=myclockd()-tk_start;
}
