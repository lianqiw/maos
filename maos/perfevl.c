/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
 * \file perfevl.c
 * Performance evaluation
 * */
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
#define TIM(A) real tk##A=myclockd()
#else
#define TIM(A)
#endif
extern int KEEP_MEM;
static void perfevl_ideal_atm(sim_t* simu, dmat* iopdevl, int ievl, real alpha){
	const parms_t* parms=simu->parms;
	const aper_t* aper=simu->aper;
	const real hs=P(parms->evl.hs,ievl);

	for(int idm=0; idm<parms->ndm; idm++){
		const real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
		real dispx=ht*P(parms->evl.thetax,ievl);
		real dispy=ht*P(parms->evl.thetay,ievl);
		real scale=1.-ht/hs;
		if(scale<0) continue;
		loc_t* locs=aper->locs;
		if(aper->locs_dm){
			locs=P(aper->locs_dm, ievl, idm);
		}
		prop_grid(P(simu->dmprojsq,idm), locs, P(iopdevl),
			alpha, dispx, dispy, scale, 0,
			0, 0);
	}
}
static void plot_psf(ccell* psf2s, const char* psfname, int closeloop, int ievl, dmat* wvl, int uselog){
	dmat* psftemp=NULL;
	const char* title, * tab;
	for(int iwvl=0; iwvl<NX(psf2s); iwvl++){
		if(closeloop){
			title="Science Closed Loop PSF";
			tab="CL";
		} else{
			title="Science Open Loop PSF";
			tab="OL";
		}
		char tabname[64];
		snprintf(tabname, sizeof(tabname), "%s%2d %.2f", tab, ievl, P(wvl,iwvl)*1e6);
		if(draw_current(psfname, tabname)){
			if(psftemp&&NX(psftemp)!=P(psf2s,iwvl)->nx){
				dfree(psftemp);
			}
			cabs22d(&psftemp, 0, P(psf2s,iwvl), 1);
			if(uselog==2){
				dcwlog10(psftemp);
			}
			ddraw(psfname, psftemp, NULL, NULL, title, "x", "y", "%s", tabname);
		}
	}
	dfree(psftemp);
}

static void perfevl_psfcl(const parms_t* parms, const aper_t* aper, const char* psfname,
	dcell* evlpsfmean, zfarr** evlpsfhist,
	dmat* iopdevl, int ievl){
/* the OPD after this time will be tilt removed. Don't use for performance
   evaluation. */
	ccell* psf2s=0;
	locfft_psf(&psf2s, aper->embed, iopdevl, parms->evl.psfsize, 0);
	int nwvl=parms->evl.nwvl;
	if(parms->evl.psfmean){
		dcell* pevlpsfmean=evlpsfmean/*PDELL*/;
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			cabs22d(&P(pevlpsfmean, iwvl, ievl), 1, P(psf2s,iwvl), 1);
		}
	}
	if(parms->evl.psfhist){
		zfarr_push(evlpsfhist[ievl], -1, psf2s);
	}
	if(parms->plot.run){
		plot_psf(psf2s, psfname, 1, ievl, parms->evl.wvl, parms->plot.psf==2);
	}
	ccellfree(psf2s);
}
/*computes the wavefront mode/error, and ngsmod. Same operation for OL/CL
  evaluations.*/
#define PERFEVL_WFE(pclep, pclmp, cleNGSmp)				\
    if(parms->recon.split){/*for split tomography */			\
	dmat*  pcleNGSmp=P(cleNGSmp,ievl)/*PDMAT*/;			\
	/*compute the dot product of wavefront with NGS mode for that direction */ \
	if(nmod==3){							\
	    calc_ngsmod_dot(PCOL(pclep,isim),PCOL(pclmp,isim), PCOL(pcleNGSmp,isim), \
			    parms,recon->ngsmod,aper, P(iopdevl),ievl);	\
	}else{/*more modes are wanted. */				\
	    calc_ngsmod_dot(NULL,NULL, PCOL(pcleNGSmp,isim),		\
			    parms,recon->ngsmod,aper, P(iopdevl),ievl);	\
	    loc_calc_mod(PCOL(pclep,isim),PCOL(pclmp,isim), aper->mod,P(aper->amp),P(iopdevl)); \
	}								\
    }else{/*for integrated tomography. */				\
	if(nmod==3){							\
	    loc_calc_ptt(PCOL(pclep,isim),PCOL(pclmp,isim), aper->locs, aper->ipcc, aper->imcc, P(aper->amp), P(iopdevl)); \
	}else{								\
	    loc_calc_mod(PCOL(pclep,isim),PCOL(pclmp,isim), aper->mod,P(aper->amp),P(iopdevl));	\
	}								\
    }									\

/**
   Performance evaluation for each direction in parallel mode.  */
void perfevl_ievl(thread_t* info){
	sim_t* simu=(sim_t*)info->data;
	const parms_t* parms=simu->parms;
	const aper_t* aper=simu->aper;
	const recon_t* recon=simu->recon;
	const int isim=simu->perfisim;
	const real atmscale=simu->atmscale?P(simu->atmscale,isim):1;
	const int nmod=parms->evl.nmod;
	const int nps=parms->atm.nps;
	const int npsr=parms->atmr.nps;
	const int imoao=parms->evl.moao;
	const real dt=parms->sim.dt;
	dmat* iopdevl=0;
	if(KEEP_MEM){
		if(!info->thread_data){
			info->thread_data=dnew(aper->locs->nloc, 1);
		}
		iopdevl=(dmat*)info->thread_data;
	}
	for(int ievl=info->start; ievl<info->end; ievl++){
		const int do_psf_cov=(parms->evl.psfmean||parms->evl.psfhist||parms->evl.cov||parms->evl.opdmean)
			&&isim>=parms->evl.psfisim&&P(parms->evl.psf,ievl);
		const int save_evlopd=parms->save.evlopd>0&&((isim+1)%parms->save.evlopd)==0;
		if(!iopdevl){
			iopdevl=dnew(aper->locs->nloc, 1);
		} else{
			dzero(iopdevl);
		}
		TIM(0);
		/*Setup pointers for easy usage */
		dmat* polmp=P(simu->olmp,ievl)/*PDMAT*/;/*OL mode for each dir */
		dmat* polep=P(simu->olep,ievl)/*PDMAT*/;/*OL error for each dir */
		dmat* pclmp=P(simu->clmp,ievl)/*PDMAT*/;
		dmat* pclep=P(simu->clep,ievl)/*PDMAT*/;

		/*atmosphere contribution. */
		if(parms->sim.idealevl){
			perfevl_ideal_atm(simu, iopdevl, ievl, 1);
		} else if(simu->atm&&!parms->sim.wfsalias){
			if(simu->evlopdground){
				dcp(&iopdevl, simu->evlopdground);
			} else{
				dzero(iopdevl);
			}
			/*fix me: the ray tracing of the same part must be performed in the same thread. */
			for(int ips=0; ips<nps; ips++){
				if(ips!=simu->perfevl_iground||!simu->evlopdground){
					int ind=ievl+parms->evl.nevl*ips;
					simu->evl_propdata_atm[ind].phiout=iopdevl;
					simu->evl_propdata_atm[ind].displacex1=-P(simu->atm,ips)->vx*isim*dt;
					simu->evl_propdata_atm[ind].displacey1=-P(simu->atm,ips)->vy*isim*dt;
					simu->evl_propdata_atm[ind].alpha=atmscale;
					CALL_THREAD(simu->evl_prop_atm[ind], 1);
				}
			}
		}
		if(simu->telws){/*Wind shake */
			real tmp=P(simu->telws,isim);
			real angle=simu->winddir?P(simu->winddir,0):0;
			real ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
			loc_add_ptt(iopdevl, ptt, aper->locs);
		}
		if(simu->telfocusreal){
			loc_add_focus(iopdevl, aper->locs, -P(P(simu->telfocusreal,0),0));
		}
		/*Add surfaces along science path. prepared in setup_surf.c */
		if(aper->opdadd&&P(aper->opdadd,ievl)){
			dadd(&iopdevl, 1, P(aper->opdadd,ievl), 1);
		}

		TIM(1);
		if(save_evlopd && simu->save->evlopdol){
			zfarr_push(simu->save->evlopdol[ievl], simu->perfisim, iopdevl);
		}
		if(parms->plot.run){
			drawopdamp("Evlol", aper->locs, iopdevl, aper->amp1, P(parms->dbg.draw_opdmax),
				"Science Open Loop OPD", "x (m)", "y (m)", "OL %d", ievl);
		}
		PERFEVL_WFE(polep, polmp, simu->oleNGSmp);
		/*evaluate time averaged open loop PSF. */
		if((parms->evl.psfmean||parms->evl.cov||parms->evl.opdmean)
			&&isim>=parms->evl.psfisim
			&&((parms->evl.psfol==1&&ievl==parms->evl.indoa)
				||(parms->evl.psfol==2&&P(parms->evl.psf,ievl)))){
			  /*Compute on axis OL psf. */
			dmat* opdevlcopy=NULL;
			if(P(parms->evl.pttr,ievl)){
				dcp(&opdevlcopy, iopdevl);
				loc_remove_ptt(opdevlcopy, PCOL(polmp, isim), aper->locs);
			} else if(parms->evl.cov||parms->evl.opdmean){
				dcp(&opdevlcopy, iopdevl);
				dadds(opdevlcopy, -P(polmp, 0, isim));
			} else{
				opdevlcopy=dref(iopdevl);
			}
			if(parms->evl.cov){
				dmm(&simu->evlopdcovol, 1, opdevlcopy, opdevlcopy, "nt", 1);
			}
			if(parms->evl.cov || parms->evl.opdmean){
				dadd(&simu->evlopdmeanol, 1, opdevlcopy, 1);
			}/*opdcov*/
			if(parms->evl.psfmean){
				ccell* psf2s=0;
				locfft_psf(&psf2s, aper->embed, opdevlcopy, parms->evl.psfsize, 0);
				int nwvl=parms->evl.nwvl;
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					cabs22d(&P(simu->evlpsfolmean,iwvl), 1, P(psf2s,iwvl), 1);
				}
				if(parms->plot.run){
					plot_psf(psf2s, "PSFol", 0, ievl, parms->evl.wvl, parms->plot.psf==2);
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
					real hl=P(parms->atmr.ht,ipsr);
					real scale=1.-hl/P(parms->evl.hs,ievl);
					if(scale<0) continue;
					real displacex=P(parms->evl.thetax,ievl)*hl;
					real displacey=P(parms->evl.thetay,ievl)*hl;
					if(parms->tomo.square){
						memcpy(&xmap, P(recon->xmap,ipsr), sizeof(map_t));
						xmap.p=P(P(simu->opdr,ipsr));
						prop_grid(&xmap, aper->locs, P(iopdevl), -1,
							displacex, displacey, scale, 0, 0, 0);
					} else{
						prop_nongrid(P(recon->xloc,ipsr), P(P(simu->opdr,ipsr)),
							aper->locs, P(iopdevl), -1,
							displacex, displacey, scale, 0, 0);
					}
				}
			}
		} else{
			/* Apply dm correction. tip/tilt command is contained in DM commands */
			wait_dmreal(simu, simu->perfisim);
			if(simu->dmreal){
				int ndm=parms->ndm;
				for(int idm=0; idm<ndm; idm++){
					int ind=ievl+parms->evl.nevl*idm;
					simu->evl_propdata_dm[ind].phiout=iopdevl;
					CALL_THREAD(simu->evl_prop_dm[ind], 1);
				}
			}
		}
		if(simu->ttmreal){
			real ptt[3]={0, -P(simu->ttmreal,0), -P(simu->ttmreal,1)};
			loc_add_ptt(iopdevl, ptt, aper->locs);
		}
		TIM(4);
		if(imoao>-1){
			dmat** dmevl=P(simu->dm_evl);
			if(dmevl[ievl]){
			/**
			   prop is faster than spmulvec. \fixme check definition of misreg
			*/
				prop_nongrid(P(recon->moao[imoao].aloc,0), P(dmevl[ievl]),
					aper->locs, P(iopdevl), -1, 0, 0, 1, 0, 0);
			}
		}

		if(parms->plot.run){
			drawopdamp("Evlcl", aper->locs, iopdevl, aper->amp1, P(parms->dbg.draw_opdmax),
				"Science Closed Loop OPD", "x (m)", "y (m)", "CL %d", ievl);
		}
		if(save_evlopd){
			zfarr_push(simu->save->evlopdcl[ievl], isim, iopdevl);
		}

		/*Evaluate closed loop performance. */
		PERFEVL_WFE(pclep, pclmp, simu->cleNGSmp);
		TIM(5);
		if(do_psf_cov){
			if((P(parms->evl.psf,ievl)&2)){
			/* even if psfpttr=1, referencing is ok.  Change to copy if
			   incompatible in the future.*/
				P(simu->evlopd,ievl)=dref(iopdevl);
			}
			if((P(parms->evl.psf,ievl)&1)){
			/** opdcov does not have p/t/t removed. do it in postproc is necessary*/
				if(P(parms->evl.pttr,ievl)){
					warning_once("Removing piston/tip/tilt from OPD.\n");
					loc_remove_ptt(iopdevl, PCOL(pclmp, isim), aper->locs);
				} else if(parms->evl.cov||parms->evl.opdmean){/*remove piston */
					dadds(iopdevl, -P(pclmp, 0, isim));
				}
				if(P(parms->evl.psfr,ievl)){
					if(parms->evl.cov){
						dmm(&P(simu->evlopdcov,ievl), 1, iopdevl, iopdevl, "nt", 1);
					}
					if(parms->evl.cov||parms->evl.opdmean){
						dadd(&P(simu->evlopdmean,ievl), 1, iopdevl, 1);
					}
				}/*opdcov */
				if(parms->evl.psfmean||parms->evl.psfhist){/*Evaluate closed loop PSF.	 */
					perfevl_psfcl(parms, aper, "PSFcl", simu->evlpsfmean, simu->save->evlpsfhist, iopdevl, ievl);
				}/*do_psf */
			}
		}
		TIM(6);
#if TIMING==1
		info2("Evl %d timing:ray atm %.4f evlol %.4f ray dm %.4f evlcl %.4f PSF %.4f\n",
			ievl, tk1-tk0, tk2-tk1, tk4-tk2, tk5-tk4, tk6-tk5);
#endif
	}
	if(!KEEP_MEM){
		dfree(iopdevl);
	}
}
/**
   Evaluation field averaged performance.
*/
static void perfevl_mean(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int isim=simu->perfisim;
	const int nevlmod=parms->evl.nmod;
	const int nevl=parms->evl.nevl;
	/*Field average the OL error */
	for(int imod=0; imod<nevlmod; imod++){
		P(simu->ole,imod,isim)=0;
		for(int ievl=0; ievl<nevl; ievl++){
			real wt=P(parms->evl.wt,ievl);
			P(simu->ole,imod,isim)+=wt*P(P(simu->olep,ievl),imod,isim);
		}
	}
	if(parms->sim.evlol){
		simu->status->clerrhi=sqrt(P(simu->ole, 0, isim))*1e9;//pttr
		simu->status->clerrlo=sqrt(P(simu->ole, 1, isim))*1e9;//tt
		return;
	}
		/*Field average the CL error */
	for(int imod=0; imod<nevlmod; imod++){
		P(simu->cle,imod,isim)=0;
		for(int ievl=0; ievl<nevl; ievl++){
			real wt=P(parms->evl.wt,ievl);
			P(simu->cle,imod,isim)+=wt*P(P(simu->clep,ievl),imod,isim);
		}
	}


	if(parms->recon.split){
	/* convert cleNGSm into mode and put NGS mode WVE into clem. */
		int nngsmod=recon->ngsmod->nmod;
		//Record NGS mode correction time history
		if(simu->corrNGSm&&P(simu->Mint_lo->mint,0)&&isim<parms->sim.end-1){
			real* pcorrNGSm=PCOL(simu->corrNGSm,isim+1);
			for(int imod=0; imod<nngsmod; imod++){
				pcorrNGSm[imod]=P(P(P(simu->Mint_lo->mint,0),0),imod);
			}
		}
		//Compute dot product of NGS mode with OPD
		real pcleNGSdot[nngsmod];
		memset(pcleNGSdot, 0, sizeof(real)*nngsmod);
		real poleNGSdot[nngsmod];
		memset(poleNGSdot, 0, sizeof(real)*nngsmod);

		for(int ievl=0; ievl<nevl; ievl++){
			real wt=P(parms->evl.wt,ievl);
			real* pcleNGSdotp=PCOL(P(simu->cleNGSmp,ievl),isim);
			real* poleNGSdotp=PCOL(P(simu->oleNGSmp,ievl),isim);
			for(int imod=0; imod<nngsmod; imod++){
				pcleNGSdot[imod]+=pcleNGSdotp[imod]*wt;
				poleNGSdot[imod]+=poleNGSdotp[imod]*wt;
			}
		}
		//Determine NGS modes from dot product
		real tt=dwdot(pcleNGSdot, recon->ngsmod->IMCC_TT, pcleNGSdot);
		real ngs=dwdot(pcleNGSdot, recon->ngsmod->IMCC, pcleNGSdot);
		real tot=P(simu->cle, 0, isim);
		real focus=0;
		if(recon->ngsmod->indfocus){
			focus=dwdot(pcleNGSdot, recon->ngsmod->IMCC_F, pcleNGSdot);
		}
		real lgs=tot-ngs;
		//Turn dot product to modes 
		real* pcleNGSm=PCOL(simu->cleNGSm, isim);
		dmulvec(pcleNGSm, recon->ngsmod->IMCC, pcleNGSdot, 1);
		real* poleNGSm=PCOL(simu->oleNGSm, isim);
		dmulvec(poleNGSm, recon->ngsmod->IMCC, poleNGSdot, 1);
		if(simu->parms->tomo.ahst_idealngs==1){
			//Magically remove NGS modes 
			tt=0;
			ngs=0;
			focus=0;
		} else if(simu->parms->tomo.ahst_idealngs==2){
			//Simulate close loop correction using ideal NGS mod.
			simu->Merr_lo=simu->Merr_lo_store;
			memcpy(P(P(simu->Merr_lo,0)), pcleNGSm, nngsmod*sizeof(real));
			//Add current low order command to form PSOL measurement for skyc.
			for(int imod=0; imod<nngsmod; imod++){
				pcleNGSm[imod]+=P(P(P(simu->Mint_lo->mint,0),0),imod);
			}
		}
		//Record NGS mode RMS error time history
		P(simu->clem, 0, isim)=lgs; /*lgs mode */
		P(simu->clem, 1, isim)=tt;    /*tt mode */
		P(simu->clem, 2, isim)=ngs;   /*ngs mod */
		if(recon->ngsmod->indfocus){
			P(simu->clem, 3, isim)=focus;
		}
		simu->status->clerrlo=sqrt(ngs)*1e9;
		simu->status->clerrhi=sqrt(lgs)*1e9;
		/*compute error spliting to tip/tilt and high order for each direction.*/
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			/*New calculation. Notice that the calculations are for
			  (phi-Hm*m)'*W*(phi-Hm*m) where m is cleNGSm, phi is opdevl
			  and W is aperture weighting. The calculation is equivalent
			  to phi'*W*phi-2*phi'*W*(Hm*m)+(Hm*m)'*W*(Hm*m).
			  phi'*W*phi is the clep phi'*W*Hm is the same as
			  pcleNGSmp'; (Hm*m)'*W*(Hm*m) = m'*(Hm'*W*Hm)*m. (Hm'*W*Hm)
			  is simply MCCp.
			*/
			if(!P(recon->ngsmod->MCCP,ievl)) continue;
			real* pcleNGSmp=PCOL(P(simu->cleNGSmp,ievl), isim);
			real sum=0;
			for(int imod=0; imod<nngsmod; imod++){
				sum+=pcleNGSmp[imod]*pcleNGSm[imod];
			}
			real sum2=dwdot(pcleNGSm, P(recon->ngsmod->MCCP,ievl), pcleNGSm);
			real tot2=P(P(simu->clep,ievl),0,isim)-2.*sum+sum2;
			P(P(simu->clemp,ievl),0,isim)=tot2;/*LGS mode */
			P(P(simu->clemp,ievl),1,isim)=tt;/*TT mode */
			P(P(simu->clemp,ievl),2,isim)=P(P(simu->clep,ievl),0,isim)-tot2;/*PR-LGS */
		}
		int do_psf=(parms->evl.psfmean||parms->evl.psfhist);
		if(isim>=parms->evl.psfisim&&(do_psf||parms->evl.cov||parms->evl.opdmean)){
#if USE_CUDA
			if(parms->gpu.evl){
				gpu_perfevl_ngsr(simu, pcleNGSm);
			} else{
#endif
				const aper_t* aper=simu->aper;
				for(int ievl=0; ievl<parms->evl.nevl; ievl++)
					if((P(parms->evl.psf,ievl)&2))
#pragma omp task
					{

						dmat* iopdevl=P(simu->evlopd,ievl);
						ngsmod2science(iopdevl, aper->locs, recon->ngsmod,
							P(parms->evl.thetax,ievl), P(parms->evl.thetay,ievl),
							pcleNGSm, -1);
						if(parms->plot.run){
							drawopdamp("Evlcl", aper->locs, iopdevl, aper->amp1, P(parms->dbg.draw_opdmax),
								"Science Closed loop OPD", "x (m)", "y (m)", "ngsr %d", ievl);
						}
						if(P(parms->evl.pttr,ievl)){
							/*we cannot use clmp because the removed ngsmod
							  has tip/tilt component*/
							real ptt[3];
							loc_calc_ptt(NULL, ptt, aper->locs, aper->ipcc, aper->imcc, P(aper->amp), P(iopdevl));
							loc_remove_ptt(iopdevl, ptt, aper->locs);
						}
						if(P(parms->evl.psfr,ievl)){
							if(parms->evl.cov){
								dmm(&P(simu->evlopdcov_ngsr,ievl), 1, iopdevl, iopdevl, "nt", 1);
							}
							if(parms->evl.opdmean){
								dadd(&P(simu->evlopdmean_ngsr,ievl), 1, iopdevl, 1);
							}
						}
						if(do_psf){
							perfevl_psfcl(parms, aper, "PSFngsr", simu->evlpsfmean_ngsr, simu->save->evlpsfhist_ngsr, iopdevl, ievl);
						}
						dfree(P(simu->evlopd,ievl));
					}
#pragma omp taskwait
#if USE_CUDA
			}
#endif
		}/*ideal ngs */
	} else{/*if not split */
		simu->status->clerrhi=sqrt(P(simu->cle,0,isim)*1e18);
		simu->status->clerrlo=sqrt(P(simu->cle,1,isim)*1e18);
	}/*if split */

	if(parms->sim.noatm==0){
		static int ct=0;
		if(P(simu->cle,0,isim)>MAX(P(simu->ole,0,isim)*100, 1e-12)){
			ct++;
			warning("Step %5d: The loop is diverging: OL: %g CL: %g\n",
				isim, sqrt(P(simu->ole,0,isim))*1e9,
				sqrt(P(simu->cle,0,isim))*1e9);
			simu->last_report_time=0; //force print out.
			if(ct>10){
				warning("Exit after this time step.\n");
				signal_caught=1;
			}
		} else{
			ct=0;
		}
	}
}
/**
   Save telemetry.
   2015-08-19: Changed psfmean to non-cumulative average.
*/
static void perfevl_save(sim_t* simu){
	const parms_t* parms=simu->parms;
	const int isim=simu->perfisim;
	if(parms->evl.psfmean&&CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.psfmean)){
		info2("Step %d: Output PSF (cumulative average).\n", isim);
		int nacc=(simu->perfisim+1-parms->evl.psfisim);//total accumulated.
		const real scale=1./(real)nacc;
		if(!parms->sim.evlol){
			dcellscale(simu->evlpsfmean, scale);
			dcell* pcl=simu->evlpsfmean/*PDELL*/;
			for(int ievl=0; ievl<parms->evl.nevl; ievl++){
				if(!simu->save->evlpsfmean[ievl]
					||!P(simu->evlpsfmean,ievl)) continue;
				for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
					free(P(pcl, iwvl, ievl)->header);
					P(pcl, iwvl, ievl)->header=evl_header(simu->parms, simu->aper, ievl, iwvl, isim);
					zfarr_push(simu->save->evlpsfmean[ievl], isim*parms->evl.nwvl+iwvl, P(pcl, iwvl, ievl));
				}
			}
			dcellscale(simu->evlpsfmean, 1./scale);
		}
		if(!parms->sim.evlol){
			dcellscale(simu->evlpsfmean_ngsr, scale);
			dcell* pcl=simu->evlpsfmean_ngsr/*PDELL*/;
			for(int ievl=0; ievl<parms->evl.nevl; ievl++){
				if(!simu->save->evlpsfmean_ngsr[ievl]
					||!P(simu->evlpsfmean_ngsr,ievl)) continue;
				for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
					free(P(pcl, iwvl, ievl)->header);
					P(pcl, iwvl, ievl)->header=evl_header(simu->parms, simu->aper, ievl, iwvl, isim);
					zfarr_push(simu->save->evlpsfmean_ngsr[ievl], isim*parms->evl.nwvl+iwvl, P(pcl, iwvl, ievl));
				}
			}
			dcellscale(simu->evlpsfmean_ngsr, 1./scale);
		}
		if(parms->evl.psfol){
			real scaleol=(parms->evl.psfol==2)?(scale/parms->evl.npsf):(scale);
			dcellscale(simu->evlpsfolmean, scaleol);
			dmat** pcl=P(simu->evlpsfolmean);
			for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
				if(!P(simu->evlpsfolmean,iwvl)) continue;
				free(pcl[iwvl]->header);
				pcl[iwvl]->header=evl_header(simu->parms, simu->aper, -1, iwvl, isim);
				zfarr_push(simu->save->evlpsfolmean, isim*parms->evl.nwvl+iwvl, pcl[iwvl]);
			}
			dcellscale(simu->evlpsfolmean, 1./scaleol);
		}
	}
	if(parms->evl.cov&&CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.cov)){
		info2("Step %d: Output opdcov (cumulative average)\n", isim);
		int nacc=(simu->perfisim+1-parms->evl.psfisim);//total accumulated.
		const real scale=1./(real)nacc;
		dcellscale(simu->evlopdcov, scale);
		dcellscale(simu->evlopdcov_ngsr, scale);
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			if(P(simu->evlopdcov, ievl)){
				zfarr_push(simu->save->evlopdcov[ievl], isim, P(simu->evlopdcov, ievl));
			}
			if(P(simu->evlopdcov_ngsr, ievl)){
				zfarr_push(simu->save->evlopdcov_ngsr[ievl], isim, P(simu->evlopdcov_ngsr, ievl));
			}
		}
		dcellscale(simu->evlopdcov, 1./scale);
		dcellscale(simu->evlopdcov_ngsr, 1./scale);
		if(parms->evl.psfol){
			const real scaleol=(parms->evl.psfol==2)?(scale/parms->evl.npsf):(scale);
			dscale(simu->evlopdcovol, scaleol);
			zfarr_push(simu->save->evlopdcovol, isim, simu->evlopdcovol);
			dscale(simu->evlopdcovol, 1./scaleol);
		}
	}
	if(parms->evl.opdmean&&CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.opdmean)){
		info2("Step %d: Output opdmean (cumulative average)\n", isim);
		int nacc=(simu->perfisim+1-parms->evl.psfisim);//total accumulated.
		const real scale=1./(real)nacc;
		dcellscale(simu->evlopdmean, scale);
		dcellscale(simu->evlopdmean_ngsr, scale);
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			if(P(simu->evlopdmean, ievl)){
				zfarr_push(simu->save->evlopdmean[ievl], isim, P(simu->evlopdmean,ievl));
				if(parms->plot.run){
					drawopdamp("Evlclm", simu->aper->locs, P(simu->evlopdmean, ievl), simu->aper->amp1, P(parms->dbg.draw_opdmax),
						"Science Closed Loop OPD Mean", "x (m)", "y (m)", "CL %d", ievl);
				}
			}
			if(P(simu->evlopdmean_ngsr, ievl)){
				zfarr_push(simu->save->evlopdmean_ngsr[ievl], isim, P(simu->evlopdmean_ngsr,ievl));
				if(parms->plot.run){
					drawopdamp("Evlclm", simu->aper->locs, P(simu->evlopdmean_ngsr, ievl), simu->aper->amp1, P(parms->dbg.draw_opdmax),
						"Science Closed Loop OPD Mean", "x (m)", "y (m)", "ngsr %d", ievl);
				}
			}
		}

		dcellscale(simu->evlopdmean, 1./scale);
		dcellscale(simu->evlopdmean_ngsr, 1./scale);
		
		if(parms->evl.psfol){
			const real scaleol=(parms->evl.psfol==2)?(scale/parms->evl.npsf):(scale);
			dscale(simu->evlopdmeanol, scaleol);
			zfarr_push(simu->save->evlopdmeanol, isim, simu->evlopdmeanol);
			if(parms->plot.run){
				drawopdamp("Evlol", simu->aper->locs, simu->evlopdmeanol, simu->aper->amp1, P(parms->dbg.draw_opdmax),
					"Science Open Loop OPD Mean", "x (m)", "y (m)", "OL mean");
			}
			dscale(simu->evlopdmeanol, 1./scaleol);
		}
	}
}


/**

   Peformance evaluation on science FoV.

   Evaluate performance by calling perfevl_ievl in parallel and then calls
   perfevl_mean to field average. Notice that the science FoV can be different
   from the DM fitting FoV, which is tuned to better sharpen the NGS

   \todo Write a standalone routine that can plot results, using techniques
   developped in drawdaemon.  */
void perfevl(sim_t* simu){
	real tk_start=PARALLEL==1?simu->tk_0:myclockd();
	const parms_t* parms=simu->parms;
	if(!parms->gpu.evl&&parms->evl.nevl>1){ //Cache the ground layer. 
		int ips=simu->perfevl_iground;
		if(ips!=-1&&simu->atm&&!parms->sim.idealevl){
			if(!simu->evlopdground){
				simu->evlopdground=dnew(simu->aper->locs->nloc, 1);
			} else{
				dzero(simu->evlopdground);
			}
			const int ievl=0;//doesn't matter for ground layer. 
			int ind=ievl+parms->evl.nevl*ips;
			const int isim=simu->perfisim;
			const real dt=parms->sim.dt;
			const real atmscale=simu->atmscale?P(simu->atmscale,isim):1;
			simu->evl_propdata_atm[ind].phiout=simu->evlopdground;
			simu->evl_propdata_atm[ind].displacex1=-P(simu->atm,ips)->vx*isim*dt;
			simu->evl_propdata_atm[ind].displacey1=-P(simu->atm,ips)->vy*isim*dt;
			simu->evl_propdata_atm[ind].alpha=atmscale;
			CALL_THREAD(simu->evl_prop_atm[ind], 1);
		}
	}
	if(PARALLEL!=1||!parms->gpu.evl){
		CALL_THREAD(simu->perfevl_pre, 0);
	}
	if(simu->perfevl_post){//only in GPU mode
		CALL_THREAD(simu->perfevl_post, 0);
	}
	perfevl_mean(simu);
#if USE_CUDA
	if(parms->gpu.evl&&parms->gpu.psf){
		gpu_perfevl_save(simu);
	} else
#endif
		perfevl_save(simu);
	simu->tk_eval=myclockd()-tk_start;
}
