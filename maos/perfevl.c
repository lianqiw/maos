/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
    const double hs = parms->evl.hs->p[ievl];
  
    for(int idm=0; idm<parms->ndm; idm++){
	const double ht = parms->dm[idm].ht+parms->dm[idm].vmisreg;
	double dispx=ht*parms->evl.thetax->p[ievl];
	double dispy=ht*parms->evl.thetay->p[ievl];
	double scale=1.-ht/hs;
	loc_t *locs=aper->locs;
	if(aper->locs_dm){
	    locs=aper->locs_dm->p[ievl+idm*parms->evl.nevl];
	}
	prop_grid(simu->dmprojsq->p[idm], locs, iopdevl->p,
		  alpha, dispx, dispy, scale, 0,
		  0, 0);
    }
}
static void perfevl_psfcl(const PARMS_T *parms, const APER_T *aper,
			  dcell *evlpsfmean, cellarr** evlpsfhist,
			  dmat *iopdevl, int ievl){
    /* the OPD after this time will be tilt removed. Don't use for performance
       evaluation. */
    ccell *psf2s=0;
    locfft_psf(&psf2s, aper->embed, iopdevl, parms->evl.psfsize, 0);
    int nwvl=parms->evl.nwvl;
    if(parms->evl.psfmean){
	PDCELL(evlpsfmean, pevlpsfmean);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cabs22d(&pevlpsfmean[ievl][iwvl], 1, psf2s->p[iwvl], 1);
	}
    }
    if(parms->evl.psfhist){
	cellarr_ccell(evlpsfhist[ievl], -1, psf2s);
    }
    if(parms->plot.run){
	dmat *psftemp=NULL;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cabs22d(&psftemp, 1, psf2s->p[iwvl], 1);
	    dcwlog10(psftemp);
	    double xylim[4]={-12,12,-12,12};
	    ddraw("CL PSF", psftemp, xylim, opdzlim, 
		  "Science Closed Loop PSF", 
		  "x", "y", "CL%2d PSF %.2f", ievl, parms->evl.wvl->p[iwvl]*1e6);
	    dfree(psftemp);
	}
    }
    ccellfree(psf2s);
}
/*computes the wavefront mode/error, and ngsmod. Same operation for OL/CL
  evaluations.*/
#define PERFEVL_WFE(pclep, pclmp, cleNGSmp)				\
    if(parms->recon.split){/*for split tomography */			\
	PDMAT(cleNGSmp->p[ievl], pcleNGSmp);				\
	/*compute the dot product of wavefront with NGS mode for that direction */ \
	if(nmod==3){							\
	    calc_ngsmod_dot(pclep[isim],pclmp[isim], pcleNGSmp[isim],	\
			    parms,recon,aper, iopdevl->p,ievl);		\
	}else if(nmod==5 ||nmod==6){/*more modes are wanted. */		\
	    calc_ngsmod_dot(NULL,NULL, pcleNGSmp[isim],			\
			    parms,recon,aper, iopdevl->p,ievl);		\
	    loc_calc_mod(pclep[isim],pclmp[isim], aper->mod,aper->amp->p,iopdevl->p); \
	}else{								\
	    error("Not implemented\n");					\
	}								\
    }else{/*for integrated tomography. */				\
	if(nmod==3){							\
	    loc_calc_ptt(pclep[isim],pclmp[isim], aper->locs, aper->ipcc, aper->imcc, aper->amp->p, iopdevl->p); \
	}else{								\
	    loc_calc_mod(pclep[isim],pclmp[isim], aper->mod,aper->amp->p,iopdevl->p);	\
	}								\
    }									\
    
/**
   Performance evaluation for each direction in parallel mode.  */
void perfevl_ievl(thread_t *info){
    SIM_T *simu=info->data;
    const PARMS_T *parms=simu->parms;
    const APER_T *aper=simu->aper;
    const RECON_T *recon=simu->recon;
    const int isim=simu->isim;
    const double atmscale=simu->atmscale?simu->atmscale->p[isim]:1;
    const int nmod=parms->evl.nmod;
    const int nps=parms->atm.nps;
    const int npsr=parms->atmr.nps;
    const int imoao=parms->evl.moao;
    const double dt=simu->dt;
    dmat *iopdevl=0;
    for(int ievl=info->start; ievl<info->end; ievl++){
	const int do_psf_cov=(parms->evl.psfmean || parms->evl.psfhist || parms->evl.cov)
	    && isim>=parms->evl.psfisim && parms->evl.psf->p[ievl];
	const int save_evlopd=parms->save.evlopd>0 && ((isim+1)%parms->save.evlopd)==0;
	if(!iopdevl){
	    iopdevl=dnew(aper->locs->nloc,1);
	}else{
	    dzero(iopdevl);
	}
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
		dcp(&iopdevl, simu->opdevlground);
	    }else{
		dzero(iopdevl);
	    }
	    /*fix me: the ray tracing of the same part must be performed in the same thread. */
	    for(int ips=0; ips<nps; ips++){
		if(ips!=simu->perfevl_iground || !simu->opdevlground){
		    int ind=ievl+parms->evl.nevl*ips;
		    simu->evl_propdata_atm[ind].phiout=iopdevl->p;
		    simu->evl_propdata_atm[ind].displacex1=-simu->atm->p[ips]->vx*isim*dt;
		    simu->evl_propdata_atm[ind].displacey1=-simu->atm->p[ips]->vy*isim*dt;
		    simu->evl_propdata_atm[ind].alpha=atmscale;
		    CALL_THREAD(simu->evl_prop_atm[ind], 0);
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
	    cellarr_dmat(simu->save->evlopdol[ievl], simu->isim, iopdevl);
	}
	if(parms->plot.run){
	    drawopdamp("OL", aper->locs,iopdevl->p , aper->amp1->p, opdzlim,
		       "Science Open Loop OPD", "x (m)", "y (m)", "OL %d", ievl);
	}
	PERFEVL_WFE(polep, polmp, simu->oleNGSmp);
	/*evaluate time averaged open loop PSF. */
	if((parms->evl.psfmean || parms->evl.cov)
	   && isim>=parms->evl.psfisim 
	   &&((parms->evl.psfol==1 && ievl==parms->evl.indoa)
	      ||(parms->evl.psfol==2 && parms->evl.psf->p[ievl]))){
	    /*Compute on axis OL psf. */
	    dmat *opdevlcopy=NULL;
	    if(parms->evl.pttr->p[ievl]){
		dcp(&opdevlcopy,iopdevl);
		loc_remove_ptt(opdevlcopy->p,polmp[isim], aper->locs);
	    }else if(parms->evl.cov){
		dcp(&opdevlcopy,iopdevl);
		dadds(opdevlcopy, -polmp[isim][0]);
	    }else{
		opdevlcopy=dref(iopdevl);
	    }
	    if(parms->evl.cov){
		dmm(&simu->evlopdcovol, 1, opdevlcopy, opdevlcopy, "nt", 1);
		dadd(&simu->evlopdmeanol, 1, opdevlcopy, 1);
	    }/*opdcov*/
	    if(parms->evl.psfmean){
		ccell *psf2s=0;
		locfft_psf(&psf2s, aper->embed, opdevlcopy, parms->evl.psfsize, 0);
		int nwvl=parms->evl.nwvl;
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    cabs22d(&simu->evlpsfolmean->p[iwvl], 1, psf2s->p[iwvl], 1);
		}
		if(parms->plot.run){
		    dmat *psftemp=NULL;
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			cabs22d(&psftemp, 1, psf2s->p[iwvl], 1);
			ddraw("OL PSF", psftemp, NULL, NULL, "Science Openloop PSF", 
			      "x", "y", "OL%2d PSF %.2f", ievl,  parms->evl.wvl->p[iwvl]*1e6);
			dfree(psftemp);
		    }
		}
		ccellfree(psf2s);
	    }
	    dfree(opdevlcopy);
	}
	if(parms->sim.evlol) continue;
	TIM(2);
	if(parms->evl.tomo){
	    /*
	      evaluate tomography performance: Apply ideal correction using
	      tomography output directly.
	    */
	    if(simu->opdr){
		map_t xmap;
		for(int ipsr=0; ipsr<npsr; ipsr++){
		    double hl=parms->atmr.ht->p[ipsr];
		    double scale = 1. - hl/parms->evl.hs->p[ievl];
		    double displacex=parms->evl.thetax->p[ievl]*hl;
		    double displacey=parms->evl.thetay->p[ievl]*hl;
		    if(parms->tomo.square){
			memcpy(&xmap, recon->xmap->p[ipsr], sizeof(map_t));
			xmap.p=simu->opdr->p[ipsr]->p;
			prop_grid(&xmap, aper->locs, iopdevl->p, -1,
				  displacex, displacey, scale, 0, 0, 0);
		    }else{
			prop_nongrid(recon->xloc->p[ipsr], simu->opdr->p[ipsr]->p,
				     aper->locs,  iopdevl->p, -1,
				     displacex, displacey, scale, 0, 0);
		    }
		}
	    }
	}else{
	    /* Apply dm correction. tip/tilt command is contained in DM commands */
	    if(simu->dmreal){
		int ndm=parms->ndm;
		for(int idm=0; idm<ndm; idm++){
		    int ind=ievl+parms->evl.nevl*idm;
		    simu->evl_propdata_dm[ind].phiout=iopdevl->p;
		    CALL_THREAD(simu->evl_prop_dm[ind], 0);
		}
	    }
	}
	if(simu->ttmreal){
	    double ptt[3]={0, -simu->ttmreal->p[0], -simu->ttmreal->p[1]};
	    loc_add_ptt(iopdevl->p, ptt, aper->locs);
	}
	TIM(4);
	if(imoao>-1){
	    dmat **dmevl=simu->dm_evl->p;
	    if(dmevl[ievl]){
		/**
		   prop is faster than spmulvec. \fixme check definition of misreg
		*/
		prop_nongrid(recon->moao[imoao].aloc->p[0],dmevl[ievl]->p,
			     aper->locs, iopdevl->p, -1, 0,0,1,0,0);
	    }
	}
	
	if(parms->plot.run){
	    drawopdamp("CL", aper->locs, iopdevl->p, aper->amp1->p,NULL,
		       "Science Closed loop OPD", "x (m)", "y (m)", "CL %d",ievl);
	}
	if(save_evlopd){
	    cellarr_dmat(simu->save->evlopdcl[ievl], isim, iopdevl);
	}

	/*Evaluate closed loop performance. */
	PERFEVL_WFE(pclep, pclmp, simu->cleNGSmp);
	if(do_psf_cov){
	    if(parms->evl.psfngsr->p[ievl]!=0){
		/* even if psfpttr=1, referencing is ok.  Change to copy if
		   incompatible in the future.*/
		simu->evlopd->p[ievl]=dref(iopdevl);
	    }
	    if(parms->evl.psfngsr->p[ievl]!=2){/*ngsr==2 means only want ngsr. */
		/** opdcov does not have p/t/t removed. do it in postproc is necessary*/
		if(parms->evl.pttr->p[ievl]){
		    if(isim==parms->evl.psfisim && ievl==0){
			warning("Removing piston/tip/tilt from OPD.\n");
		    }
		    loc_remove_ptt(iopdevl->p, pclmp[isim], aper->locs);
		}else if(parms->evl.cov){/*remove piston */
		    dadds(iopdevl, -pclmp[isim][0]);
		}
		if(parms->evl.cov){
		    dmm(&simu->evlopdcov->p[ievl], 1, iopdevl, iopdevl, "nt", 1);
		    dadd(&simu->evlopdmean->p[ievl], 1, iopdevl, 1);
		}/*opdcov */
		if(parms->evl.psfmean || parms->evl.psfhist){/*Evaluate closed loop PSF.	 */
		    perfevl_psfcl(parms, aper, simu->evlpsfmean, simu->save->evlpsfhist, iopdevl, ievl);
		}/*do_psf */
	    }
	}
	TIM(5);
#if TIMING==1
	info2("Evl %d timing:ray atm %.4f evlol %.4f ray dm %.4f evlcl %.4f\n",
	      ievl, tk1-tk0, tk2-tk1, tk4-tk2, tk5-tk4);
#endif
    }
    dfree(iopdevl);
}
/**
   Evaluation field averaged performance.
*/
static void perfevl_mean(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    const int isim=simu->isim;
    const int nmod=parms->evl.nmod;
    const int nevl=parms->evl.nevl;
    /*Field average the OL error */
    for(int imod=0; imod<nmod; imod++){
	int ind=imod+nmod*isim;
	simu->ole->p[ind]=0;
	for(int ievl=0; ievl<nevl; ievl++){
	    double wt=parms->evl.wt->p[ievl];
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
	    double wt=parms->evl.wt->p[ievl];
	    simu->cle->p[ind]+=wt*simu->clep->p[ievl]->p[ind];
	}
    }
  

    if(parms->recon.split){
	/* convert cleNGSm into mode and put NGS mode WVE into clem. */
	int nngsmod=recon->ngsmod->nmod;
	    
	if(simu->corrNGSm && simu->Mint_lo->mint->p[0] && isim<parms->sim.end-1){
	    double *pcorrNGSm=simu->corrNGSm->p+(isim+1)*nngsmod;
	    for(int imod=0; imod<nngsmod; imod++){
		pcorrNGSm[imod]=simu->Mint_lo->mint->p[0]->p[0]->p[imod];
	    }
	}
	double pcleNGSdot[nngsmod];
	memset(pcleNGSdot, 0, sizeof(double)*nngsmod);
	double poleNGSdot[nngsmod];
	memset(poleNGSdot, 0, sizeof(double)*nngsmod);
	    
	for(int ievl=0; ievl<nevl; ievl++){
	    double wt=parms->evl.wt->p[ievl];
	    double *pcleNGSdotp=simu->cleNGSmp->p[ievl]->p+isim*nngsmod;
	    double *poleNGSdotp=simu->oleNGSmp->p[ievl]->p+isim*nngsmod;
	    for(int imod=0; imod<nngsmod; imod++){
		pcleNGSdot[imod]+=pcleNGSdotp[imod]*wt;
		poleNGSdot[imod]+=poleNGSdotp[imod]*wt;
	    }
	}
	double tt=dwdot2(pcleNGSdot, recon->ngsmod->IMCC_TT, pcleNGSdot);
	double ngs=dwdot(pcleNGSdot, recon->ngsmod->IMCC,    pcleNGSdot);
	double tot=simu->cle->p[isim*nmod];
	double lgs=tot-ngs;
	/*turn dot product to modes */
	double *pcleNGSm=simu->cleNGSm->p+isim*nngsmod;
	dmulvec(pcleNGSm,recon->ngsmod->IMCC,pcleNGSdot,1);
	double *poleNGSm=simu->oleNGSm->p+isim*nngsmod;
	dmulvec(poleNGSm,recon->ngsmod->IMCC,poleNGSdot,1);
	if(simu->parms->tomo.ahst_idealngs){
	    /*we no longer add ideal ngs modes to DM. Instead, we add them
	      to NGS wavefront only. This prevents it to capture all the
	      focus mode in the atmosphere.*/
	    ngs=0; tt=0;
	}
	simu->clem->p[isim*3]=lgs; /*lgs mode */
	simu->clem->p[isim*3+1]=tt;    /*tt mode */
	simu->clem->p[isim*3+2]=ngs;   /*ngs mod */
	simu->status->clerrlo=sqrt(ngs)*1e9;
	simu->status->clerrhi=sqrt(lgs)*1e9;
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
	int do_psf=(parms->evl.psfmean || parms->evl.psfhist);
	if(isim>=parms->evl.psfisim && (do_psf || parms->evl.cov)){
	    /*Only here if NGS mode removal flag is set (evl.psfngsr[ievl])*/
	    /*2013-01-23: Was using dot product before converting to modes. Fixed.*/
#if USE_CUDA
	    if(parms->gpu.evl){
		gpu_perfevl_ngsr(simu, pcleNGSm);
	    }else{
#endif
		const APER_T *aper=simu->aper;
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		    if(!parms->evl.psf->p[ievl] || !parms->evl.psfngsr->p[ievl]) continue;
		    dmat *iopdevl=simu->evlopd->p[ievl];
		    if(!iopdevl) continue;
		    ngsmod2science(iopdevl, aper->locs, recon->ngsmod, 
				   parms->evl.thetax->p[ievl], parms->evl.thetay->p[ievl],
				   pcleNGSm, -1);
		    if(parms->evl.pttr->p[ievl]){
			/*we cannot use clmp because the removed ngsmod
			  has tip/tilt component*/
			double ptt[3];
			loc_calc_ptt(NULL, ptt, aper->locs, aper->ipcc, aper->imcc, aper->amp->p, iopdevl->p);
			loc_remove_ptt(iopdevl->p, ptt, aper->locs);
		    }
		    if(parms->evl.cov){
			dmm(&simu->evlopdcov_ngsr->p[ievl], 1, iopdevl, iopdevl, "nt", 1);
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
    }else{/*if not split */
	simu->status->clerrhi=sqrt(simu->cle->p[nmod*isim]*1e18);
	simu->status->clerrlo=sqrt(simu->cle->p[nmod*isim+1]*1e18);
    }/*if split */
    
    if(parms->sim.noatm==0 && simu->cle->p[nmod*isim] > MAX(simu->ole->p[nmod*isim]*100, 1e-13)){
	static int ct=0; 
	ct++;
	if(ct>10){
	    error("Divergent simulation.");
	}
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
	double scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	if(!parms->sim.evlol){
	    dcellscale(simu->evlpsfmean, scale);
	    PDCELL(simu->evlpsfmean, pcl);
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!simu->save->evlpsfmean[ievl]
		    ||!simu->evlpsfmean->p[ievl]) continue;
		for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
		    if(!pcl[ievl][iwvl]->header){
			pcl[ievl][iwvl]->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		    }
		    cellarr_dmat(simu->save->evlpsfmean[ievl], isim*parms->evl.nwvl+iwvl, pcl[ievl][iwvl]);
		}
	    }
	    dcellscale(simu->evlpsfmean, 1./scale);//scale it back;
	}
	if(!parms->sim.evlol){
	    dcellscale(simu->evlpsfmean_ngsr, scale);
	    PDCELL(simu->evlpsfmean_ngsr, pcl);
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		if(!simu->save->evlpsfmean_ngsr[ievl]
		    ||!simu->evlpsfmean_ngsr->p[ievl]) continue;
		for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
		    if(!pcl[ievl][iwvl]->header){
			pcl[ievl][iwvl]->header=evl_header(simu->parms, simu->aper, ievl, iwvl);
		    }
		    cellarr_dmat(simu->save->evlpsfmean_ngsr[ievl], isim*parms->evl.nwvl+iwvl, pcl[ievl][iwvl]);
		}
	    }
	    dcellscale(simu->evlpsfmean_ngsr, 1./scale);//scale it back;
	}
	if(parms->evl.psfol){
	    scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    if(parms->evl.psfol==2){
		scale=scale/parms->evl.npsf;
	    }
	    dcellscale(simu->evlpsfolmean, scale);
	    dmat **pcl=simu->evlpsfolmean->p;
	    for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
		if(!simu->evlpsfolmean->p[iwvl]) continue;
		if(!pcl[iwvl]->header){
		    pcl[iwvl]->header=evl_header(simu->parms, simu->aper, -1, iwvl);
		}
		cellarr_dmat(simu->save->evlpsfolmean, isim*parms->evl.nwvl+iwvl, pcl[iwvl]);
	    }
	    dcellscale(simu->evlpsfolmean, 1./scale);//scale it back;
	}
    }
    if(parms->evl.cov && CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.cov)){
	info2("Step %d: Output opdcov\n", isim);
	long nstep=isim+1-parms->evl.psfisim;
	double scale=1./nstep;
	dcellscale(simu->evlopdcov, scale);
	dcellscale(simu->evlopdmean, scale);
	dcellscale(simu->evlopdcov_ngsr, scale);
	dcellscale(simu->evlopdmean_ngsr, scale);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!simu->evlopdcov->p[ievl]) continue;
	    cellarr_dmat(simu->save->evlopdcov[ievl], isim, simu->evlopdcov->p[ievl]);
	    cellarr_dmat(simu->save->evlopdmean[ievl], isim, simu->evlopdmean->p[ievl]);
	}
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    if(!simu->evlopdcov_ngsr->p[ievl]) continue;
	    cellarr_dmat(simu->save->evlopdcov_ngsr[ievl], isim, simu->evlopdcov_ngsr->p[ievl]);
	    cellarr_dmat(simu->save->evlopdmean_ngsr[ievl], isim, simu->evlopdmean_ngsr->p[ievl]);
	}
	dcellscale(simu->evlopdcov, 1./scale);//scale it back;
	dcellscale(simu->evlopdmean, 1./scale);//scale it back;
	dcellscale(simu->evlopdcov_ngsr, 1./scale);//scale it back;
	dcellscale(simu->evlopdmean_ngsr, 1./scale);//scale it back;
	if(parms->evl.psfol){
	    scale=1./(double)(simu->isim+1-parms->evl.psfisim);
	    if(parms->evl.psfol==2){
		scale=scale/parms->evl.npsf;
	    }
	    dscale(simu->evlopdcovol, scale);
	    dscale(simu->evlopdmeanol, scale);
	    cellarr_dmat(simu->save->evlopdcovol, isim, simu->evlopdcovol);
	    cellarr_dmat(simu->save->evlopdmeanol, isim, simu->evlopdmeanol);
	    dscale(simu->evlopdcovol, 1./scale);//scale it back;
	    dscale(simu->evlopdmeanol, 1./scale);//scale it back;
	}
    }
}
/**
   Evaluate performance by calling perfevl_ievl in parallel and then calls
   perfevl_mean to field average.  */
void perfevl(SIM_T *simu){
    double tk_start=myclockd();
    const PARMS_T *parms=simu->parms;
    if(!(parms->gpu.evl) && parms->evl.nevl>1){ //Cache the ground layer. 
	int ips=simu->perfevl_iground;
	if(ips!=-1 && simu->atm && !parms->sim.idealevl){
	    if(!simu->opdevlground){
		simu->opdevlground=dnew(simu->aper->locs->nloc,1);
	    }else{
		dzero(simu->opdevlground);
	    }
	    const int ievl=0;//doesn't matter for ground layer. 
	    int ind=ievl+parms->evl.nevl*ips;
	    const int isim=simu->isim;
	    const double dt=simu->dt;
	    const double atmscale=simu->atmscale?simu->atmscale->p[isim]:1;
	    simu->evl_propdata_atm[ind].phiout=simu->opdevlground->p;
	    simu->evl_propdata_atm[ind].displacex1=-simu->atm->p[ips]->vx*isim*dt;
	    simu->evl_propdata_atm[ind].displacey1=-simu->atm->p[ips]->vy*isim*dt;
	    simu->evl_propdata_atm[ind].alpha=atmscale;
	    CALL_THREAD(simu->evl_prop_atm[ind], 0);
	}
    }
    extern int PARALLEL;
    if(!PARALLEL || !parms->gpu.evl){
	CALL_THREAD(simu->perf_evl_pre, 0);
    }
    if(simu->perf_evl_post){//only in GPU mode
	CALL_THREAD(simu->perf_evl_post, 0);
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
