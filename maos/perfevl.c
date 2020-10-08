/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
static void perfevl_ideal_atm(SIM_T* simu, dmat* iopdevl, int ievl, real alpha){
	const PARMS_T* parms=simu->parms;
	const APER_T* aper=simu->aper;
	const real hs=parms->evl.hs->p[ievl];

	for(int idm=0; idm<parms->ndm; idm++){
		const real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
		real dispx=ht*parms->evl.thetax->p[ievl];
		real dispy=ht*parms->evl.thetay->p[ievl];
		real scale=1.-ht/hs;
		if(scale<0) continue;
		loc_t* locs=aper->locs;
		if(aper->locs_dm){
			locs=aper->locs_dm->p[ievl+idm*parms->evl.nevl];
		}
		prop_grid(simu->dmprojsq->p[idm], locs, iopdevl->p,
			alpha, dispx, dispy, scale, 0,
			0, 0);
	}
}
static void plot_psf(ccell* psf2s, const char* psfname, int closeloop, int ievl, dmat* wvl, int uselog){
	dmat* psftemp=NULL;
	const char* title, * tab;
	for(int iwvl=0; iwvl<psf2s->nx; iwvl++){
		if(closeloop){
			title="Science Closed Loop PSF";
			tab="CL";
		} else{
			title="Science Open Loop PSF";
			tab="OL";
		}
		char tabname[64];
		snprintf(tabname, sizeof(tabname), "%s%2d %.2f", tab, ievl, wvl->p[iwvl]*1e6);
		if(draw_current(psfname, tabname)){
			if(psftemp&&psftemp->nx!=psf2s->p[iwvl]->nx){
				dfree(psftemp);
			}
			cabs22d(&psftemp, 0, psf2s->p[iwvl], 1);
			if(uselog==2){
				dcwlog10(psftemp);
			}


			ddraw(psfname, psftemp, NULL, NULL, title, "x", "y", "%s", tabname);
		}
	}
	dfree(psftemp);
}

static void perfevl_psfcl(const PARMS_T* parms, const APER_T* aper, const char* psfname,
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
			cabs22d(PP(pevlpsfmean, iwvl, ievl), 1, psf2s->p[iwvl], 1);
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
	dmat*  pcleNGSmp=cleNGSmp->p[ievl]/*PDMAT*/;			\
	/*compute the dot product of wavefront with NGS mode for that direction */ \
	if(nmod==3){							\
	    calc_ngsmod_dot(PCOL(pclep,isim),PCOL(pclmp,isim), PCOL(pcleNGSmp,isim), \
			    parms,recon->ngsmod,aper, iopdevl->p,ievl);	\
	}else{/*more modes are wanted. */				\
	    calc_ngsmod_dot(NULL,NULL, PCOL(pcleNGSmp,isim),		\
			    parms,recon->ngsmod,aper, iopdevl->p,ievl);	\
	    loc_calc_mod(PCOL(pclep,isim),PCOL(pclmp,isim), aper->mod,aper->amp->p,iopdevl->p); \
	}								\
    }else{/*for integrated tomography. */				\
	if(nmod==3){							\
	    loc_calc_ptt(PCOL(pclep,isim),PCOL(pclmp,isim), aper->locs, aper->ipcc, aper->imcc, aper->amp->p, iopdevl->p); \
	}else{								\
	    loc_calc_mod(PCOL(pclep,isim),PCOL(pclmp,isim), aper->mod,aper->amp->p,iopdevl->p);	\
	}								\
    }									\

/**
   Performance evaluation for each direction in parallel mode.  */
void perfevl_ievl(thread_t* info){
	SIM_T* simu=(SIM_T*)info->data;
	const PARMS_T* parms=simu->parms;
	const APER_T* aper=simu->aper;
	const RECON_T* recon=simu->recon;
	const int isim=simu->perfisim;
	const real atmscale=simu->atmscale?simu->atmscale->p[isim]:1;
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
		const int do_psf_cov=(parms->evl.psfmean||parms->evl.psfhist||parms->evl.cov)
			&&isim>=parms->evl.psfisim&&parms->evl.psf->p[ievl];
		const int save_evlopd=parms->save.evlopd>0&&((isim+1)%parms->save.evlopd)==0;
		if(!iopdevl){
			iopdevl=dnew(aper->locs->nloc, 1);
		} else{
			dzero(iopdevl);
		}
		TIM(0);
		/*Setup pointers for easy usage */
		dmat* polmp=simu->olmp->p[ievl]/*PDMAT*/;/*OL mode for each dir */
		dmat* polep=simu->olep->p[ievl]/*PDMAT*/;/*OL error for each dir */
		dmat* pclmp=simu->clmp->p[ievl]/*PDMAT*/;
		dmat* pclep=simu->clep->p[ievl]/*PDMAT*/;

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
					simu->evl_propdata_atm[ind].phiout=iopdevl->p;
					simu->evl_propdata_atm[ind].displacex1=-simu->atm->p[ips]->vx*isim*dt;
					simu->evl_propdata_atm[ind].displacey1=-simu->atm->p[ips]->vy*isim*dt;
					simu->evl_propdata_atm[ind].alpha=atmscale;
					CALL_THREAD(simu->evl_prop_atm[ind], 0);
				}
			}
		}
		if(simu->telws){/*Wind shake */
			real tmp=simu->telws->p[isim];
			real angle=simu->winddir?simu->winddir->p[0]:0;
			real ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
			loc_add_ptt(iopdevl->p, ptt, aper->locs);
		}
		if(simu->telfocusreal){
			loc_add_focus(iopdevl->p, aper->locs, -simu->telfocusreal->p[0]->p[0]);
		}
		/*Add surfaces along science path. prepared in setup_surf.c */
		if(aper->opdadd&&aper->opdadd->p[ievl]){
			dadd(&iopdevl, 1, aper->opdadd->p[ievl], 1);
		}

		TIM(1);
		if(save_evlopd){
			zfarr_push(simu->save->evlopdol[ievl], simu->perfisim, iopdevl);
		}
		if(parms->plot.run){
			drawopdamp("Evlol", aper->locs, iopdevl->p, aper->amp1->p, parms->dbg.draw_opdmax->p,
				"Science Open Loop OPD", "x (m)", "y (m)", "OL %d", ievl);
		}
		PERFEVL_WFE(polep, polmp, simu->oleNGSmp);
		/*evaluate time averaged open loop PSF. */
		if((parms->evl.psfmean||parms->evl.cov)
			&&isim>=parms->evl.psfisim
			&&((parms->evl.psfol==1&&ievl==parms->evl.indoa)
				||(parms->evl.psfol==2&&parms->evl.psf->p[ievl]))){
			  /*Compute on axis OL psf. */
			dmat* opdevlcopy=NULL;
			if(parms->evl.pttr->p[ievl]){
				dcp(&opdevlcopy, iopdevl);
				loc_remove_ptt(opdevlcopy->p, PCOL(polmp, isim), aper->locs);
			} else if(parms->evl.cov){
				dcp(&opdevlcopy, iopdevl);
				dadds(opdevlcopy, -P(polmp, 0, isim));
			} else{
				opdevlcopy=dref(iopdevl);
			}
			if(parms->evl.cov){
				dmm(&simu->evlopdcovol, 1, opdevlcopy, opdevlcopy, "nt", 1);
				dadd(&simu->evlopdmeanol, 1, opdevlcopy, 1);
			}/*opdcov*/
			if(parms->evl.psfmean){
				ccell* psf2s=0;
				locfft_psf(&psf2s, aper->embed, opdevlcopy, parms->evl.psfsize, 0);
				int nwvl=parms->evl.nwvl;
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					cabs22d(&simu->evlpsfolmean->p[iwvl], 1, psf2s->p[iwvl], 1);
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
					real hl=parms->atmr.ht->p[ipsr];
					real scale=1.-hl/parms->evl.hs->p[ievl];
					if(scale<0) continue;
					real displacex=parms->evl.thetax->p[ievl]*hl;
					real displacey=parms->evl.thetay->p[ievl]*hl;
					if(parms->tomo.square){
						memcpy(&xmap, recon->xmap->p[ipsr], sizeof(map_t));
						xmap.p=simu->opdr->p[ipsr]->p;
						prop_grid(&xmap, aper->locs, iopdevl->p, -1,
							displacex, displacey, scale, 0, 0, 0);
					} else{
						prop_nongrid(recon->xloc->p[ipsr], simu->opdr->p[ipsr]->p,
							aper->locs, iopdevl->p, -1,
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
					simu->evl_propdata_dm[ind].phiout=iopdevl->p;
					CALL_THREAD(simu->evl_prop_dm[ind], 0);
				}
			}
		}
		if(simu->ttmreal){
			real ptt[3]={0, -simu->ttmreal->p[0], -simu->ttmreal->p[1]};
			loc_add_ptt(iopdevl->p, ptt, aper->locs);
		}
		TIM(4);
		if(imoao>-1){
			dmat** dmevl=simu->dm_evl->p;
			if(dmevl[ievl]){
			/**
			   prop is faster than spmulvec. \fixme check definition of misreg
			*/
				prop_nongrid(recon->moao[imoao].aloc->p[0], dmevl[ievl]->p,
					aper->locs, iopdevl->p, -1, 0, 0, 1, 0, 0);
			}
		}

		if(parms->plot.run){
			drawopdamp("Evlcl", aper->locs, iopdevl->p, aper->amp1->p, parms->dbg.draw_opdmax->p,
				"Science Closed loop OPD", "x (m)", "y (m)", "CL %d", ievl);
		}
		if(save_evlopd){
			zfarr_push(simu->save->evlopdcl[ievl], isim, iopdevl);
		}

		/*Evaluate closed loop performance. */
		PERFEVL_WFE(pclep, pclmp, simu->cleNGSmp);
		TIM(5);
		if(do_psf_cov){
			if(parms->evl.psfngsr->p[ievl]!=0){
			/* even if psfpttr=1, referencing is ok.  Change to copy if
			   incompatible in the future.*/
				simu->evlopd->p[ievl]=dref(iopdevl);
			}
			if(parms->evl.psfngsr->p[ievl]!=2){/*ngsr==2 means only want ngsr. */
			/** opdcov does not have p/t/t removed. do it in postproc is necessary*/
				if(parms->evl.pttr->p[ievl]){
					warning_once("Removing piston/tip/tilt from OPD.\n");
					loc_remove_ptt(iopdevl->p, PCOL(pclmp, isim), aper->locs);
				} else if(parms->evl.cov){/*remove piston */
					dadds(iopdevl, -P(pclmp, 0, isim));
				}
				if(parms->evl.cov&&parms->evl.psfr->p[ievl]){
					dmm(&simu->evlopdcov->p[ievl], 1, iopdevl, iopdevl, "nt", 1);
					dadd(&simu->evlopdmean->p[ievl], 1, iopdevl, 1);
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
static void perfevl_mean(SIM_T* simu){
	const PARMS_T* parms=simu->parms;
	const RECON_T* recon=simu->recon;
	const int isim=simu->perfisim;
	const int nevlmod=parms->evl.nmod;
	const int nevl=parms->evl.nevl;
	/*Field average the OL error */
	for(int imod=0; imod<nevlmod; imod++){
		int ind=imod+nevlmod*isim;
		simu->ole->p[ind]=0;
		for(int ievl=0; ievl<nevl; ievl++){
			real wt=parms->evl.wt->p[ievl];
			simu->ole->p[ind]+=wt*simu->olep->p[ievl]->p[ind];
		}
	}
	if(parms->sim.evlol)
		return;

		/*Field average the CL error */
	for(int imod=0; imod<nevlmod; imod++){
		int ind=imod+nevlmod*isim;
		simu->cle->p[ind]=0;
		for(int ievl=0; ievl<nevl; ievl++){
			real wt=parms->evl.wt->p[ievl];
			simu->cle->p[ind]+=wt*simu->clep->p[ievl]->p[ind];
		}
	}


	if(parms->recon.split){
	/* convert cleNGSm into mode and put NGS mode WVE into clem. */
		int nngsmod=recon->ngsmod->nmod;
		//Record NGS mode correction time history
		if(simu->corrNGSm&&simu->Mint_lo->mint->p[0]&&isim<parms->sim.end-1){
			real* pcorrNGSm=simu->corrNGSm->p+(isim+1)*nngsmod;
			for(int imod=0; imod<nngsmod; imod++){
				pcorrNGSm[imod]=simu->Mint_lo->mint->p[0]->p[0]->p[imod];
			}
		}
		//Compute dot product of NGS mode with OPD
		real pcleNGSdot[nngsmod];
		memset(pcleNGSdot, 0, sizeof(real)*nngsmod);
		real poleNGSdot[nngsmod];
		memset(poleNGSdot, 0, sizeof(real)*nngsmod);

		for(int ievl=0; ievl<nevl; ievl++){
			real wt=parms->evl.wt->p[ievl];
			real* pcleNGSdotp=simu->cleNGSmp->p[ievl]->p+isim*nngsmod;
			real* poleNGSdotp=simu->oleNGSmp->p[ievl]->p+isim*nngsmod;
			for(int imod=0; imod<nngsmod; imod++){
				pcleNGSdot[imod]+=pcleNGSdotp[imod]*wt;
				poleNGSdot[imod]+=poleNGSdotp[imod]*wt;
			}
		}
		//Determine NGS modes from dot product
		real tt=dwdot(pcleNGSdot, recon->ngsmod->IMCC_TT, pcleNGSdot);
		real ngs=dwdot(pcleNGSdot, recon->ngsmod->IMCC, pcleNGSdot);
		real tot=simu->cle->p[isim*nevlmod];
		real focus=0;
		if(recon->ngsmod->indfocus){
			focus=dwdot(pcleNGSdot, recon->ngsmod->IMCC_F, pcleNGSdot);
		}
		real lgs=tot-ngs;
		//Turn dot product to modes 
		real* pcleNGSm=simu->cleNGSm->p+isim*nngsmod;
		dmulvec(pcleNGSm, recon->ngsmod->IMCC, pcleNGSdot, 1);
		real* poleNGSm=simu->oleNGSm->p+isim*nngsmod;
		dmulvec(poleNGSm, recon->ngsmod->IMCC, poleNGSdot, 1);
		if(simu->parms->tomo.ahst_idealngs==1){
			//Magically remove NGS modes 
			tt=0;
			ngs=0;
			focus=0;
		} else if(simu->parms->tomo.ahst_idealngs==2){
			//Simulate close loop correction using ideal NGS mod.
			simu->Merr_lo=simu->Merr_lo_store;
			memcpy(simu->Merr_lo->p[0]->p, pcleNGSm, nngsmod*sizeof(real));
			//Add current low order command to form PSOL measurement for skyc.
			for(int imod=0; imod<nngsmod; imod++){
				pcleNGSm[imod]+=simu->Mint_lo->mint->p[0]->p[0]->p[imod];
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
			if(!recon->ngsmod->MCCP->p[ievl]) continue;
			real* pcleNGSmp=PCOL(simu->cleNGSmp->p[ievl], isim);
			real sum=0;
			for(int imod=0; imod<nngsmod; imod++){
				sum+=pcleNGSmp[imod]*pcleNGSm[imod];
			}
			real sum2=dwdot(pcleNGSm, recon->ngsmod->MCCP->p[ievl], pcleNGSm);
			real tot2=simu->clep->p[ievl]->p[isim*nevlmod]-2.*sum+sum2;
			simu->clemp->p[ievl]->p[isim*3]=tot2;/*LGS mode */
			simu->clemp->p[ievl]->p[isim*3+1]=tt;/*TT mode */
			simu->clemp->p[ievl]->p[isim*3+2]=simu->clep->p[ievl]->p[nevlmod*isim]-tot2;/*PR-LGS */
		}
		int do_psf=(parms->evl.psfmean||parms->evl.psfhist);
		if(isim>=parms->evl.psfisim&&(do_psf||parms->evl.cov)){
			/*Only here if NGS mode removal flag is set (evl.psfngsr[ievl])*/
			/*2013-01-23: Was using dot product before converting to modes. Fixed.*/
#if USE_CUDA
			if(parms->gpu.evl){
				gpu_perfevl_ngsr(simu, pcleNGSm);
			} else{
#endif
				const APER_T* aper=simu->aper;
				for(int ievl=0; ievl<parms->evl.nevl; ievl++)
					if(parms->evl.psf->p[ievl]&&parms->evl.psfngsr->p[ievl])
#pragma omp task
					{

						dmat* iopdevl=simu->evlopd->p[ievl];
						ngsmod2science(iopdevl, aper->locs, recon->ngsmod,
							parms->evl.thetax->p[ievl], parms->evl.thetay->p[ievl],
							pcleNGSm, -1);
						if(parms->plot.run){
							drawopdamp("Evlcl", aper->locs, iopdevl->p, aper->amp1->p, parms->dbg.draw_opdmax->p,
								"Science Closed loop OPD", "x (m)", "y (m)", "ngsr %d", ievl);
						}
						if(parms->evl.pttr->p[ievl]){
							/*we cannot use clmp because the removed ngsmod
							  has tip/tilt component*/
							real ptt[3];
							loc_calc_ptt(NULL, ptt, aper->locs, aper->ipcc, aper->imcc, aper->amp->p, iopdevl->p);
							loc_remove_ptt(iopdevl->p, ptt, aper->locs);
						}
						if(parms->evl.cov&&parms->evl.psfr->p[ievl]){
							dmm(&simu->evlopdcov_ngsr->p[ievl], 1, iopdevl, iopdevl, "nt", 1);
							dadd(&simu->evlopdmean_ngsr->p[ievl], 1, iopdevl, 1);
						}
						if(do_psf){
							perfevl_psfcl(parms, aper, "PSFngsr", simu->evlpsfmean_ngsr, simu->save->evlpsfhist_ngsr, iopdevl, ievl);
						}
						dfree(simu->evlopd->p[ievl]);
					}
#pragma omp taskwait
#if USE_CUDA
			}
#endif
		}/*ideal ngs */
	} else{/*if not split */
		simu->status->clerrhi=sqrt(simu->cle->p[nevlmod*isim]*1e18);
		simu->status->clerrlo=sqrt(simu->cle->p[nevlmod*isim+1]*1e18);
	}/*if split */

	if(parms->sim.noatm==0){
		static int ct=0;
		if(simu->cle->p[nevlmod*isim]>MAX(simu->ole->p[nevlmod*isim]*100, 1e-12)){
			ct++;
			if(ct>10){
				sync();
				error("Divergent simulation.");
			}
			warning("Step %5d: The loop is diverging: OL: %g CL: %g\n",
				isim, sqrt(simu->ole->p[nevlmod*isim])*1e9,
				sqrt(simu->cle->p[nevlmod*isim])*1e9);
			simu->last_report_time=0; //force print out.
		} else{
			ct=0;
		}
	}
}
/**
   Save telemetry.
   2015-08-19: Changed psfmean to non-cumulative average.
*/
static void perfevl_save(SIM_T* simu){
	const PARMS_T* parms=simu->parms;
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
					||!simu->evlpsfmean->p[ievl]) continue;
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
					||!simu->evlpsfmean_ngsr->p[ievl]) continue;
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
			dmat** pcl=simu->evlpsfolmean->p;
			for(int iwvl=0; iwvl<parms->evl.nwvl; iwvl++){
				if(!simu->evlpsfolmean->p[iwvl]) continue;
				free(pcl[iwvl]->header);
				pcl[iwvl]->header=evl_header(simu->parms, simu->aper, -1, iwvl, isim);
				zfarr_push(simu->save->evlpsfolmean, isim*parms->evl.nwvl+iwvl, pcl[iwvl]);
			}
			dcellscale(simu->evlpsfolmean, 1./scaleol);
		}
	}
	if(parms->evl.cov&&CHECK_SAVE(parms->evl.psfisim, parms->sim.end, isim, parms->evl.cov)){
		info2("Step %d: Output opdcov (non-cumulative average)\n", isim);
		int nacc=(simu->perfisim+1-parms->evl.psfisim);//total accumulated.
		const real scale=1./(real)nacc;
		dcellscale(simu->evlopdcov, scale);
		dcellscale(simu->evlopdmean, scale);
		dcellscale(simu->evlopdcov_ngsr, scale);
		dcellscale(simu->evlopdmean_ngsr, scale);
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			if(!simu->evlopdcov->p[ievl]) continue;
			zfarr_push(simu->save->evlopdcov[ievl], isim, simu->evlopdcov->p[ievl]);
			zfarr_push(simu->save->evlopdmean[ievl], isim, simu->evlopdmean->p[ievl]);
		}
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			if(!simu->evlopdcov_ngsr->p[ievl]) continue;
			zfarr_push(simu->save->evlopdcov_ngsr[ievl], isim, simu->evlopdcov_ngsr->p[ievl]);
			zfarr_push(simu->save->evlopdmean_ngsr[ievl], isim, simu->evlopdmean_ngsr->p[ievl]);
		}
		dcellscale(simu->evlopdcov, 1./scale);
		dcellscale(simu->evlopdmean, 1./scale);
		dcellscale(simu->evlopdcov_ngsr, 1./scale);
		dcellscale(simu->evlopdmean_ngsr, 1./scale);
		if(parms->evl.psfol){
			const real scaleol=(parms->evl.psfol==2)?(scale/parms->evl.npsf):(scale);
			dscale(simu->evlopdcovol, scaleol);
			dscale(simu->evlopdmeanol, scaleol);
			zfarr_push(simu->save->evlopdcovol, isim, simu->evlopdcovol);
			zfarr_push(simu->save->evlopdmeanol, isim, simu->evlopdmeanol);
			dscale(simu->evlopdcovol, 1./scaleol);
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
void perfevl(SIM_T* simu){
	real tk_start=PARALLEL==1?simu->tk_0:myclockd();
	const PARMS_T* parms=simu->parms;
	if(!(parms->gpu.evl)&&parms->evl.nevl>1){ //Cache the ground layer. 
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
			const real atmscale=simu->atmscale?simu->atmscale->p[isim]:1;
			simu->evl_propdata_atm[ind].phiout=simu->evlopdground->p;
			simu->evl_propdata_atm[ind].displacex1=-simu->atm->p[ips]->vx*isim*dt;
			simu->evl_propdata_atm[ind].displacey1=-simu->atm->p[ips]->vy*isim*dt;
			simu->evl_propdata_atm[ind].alpha=atmscale;
			CALL_THREAD(simu->evl_prop_atm[ind], 0);
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
