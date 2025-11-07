/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file maos/recon.h
   Carry out wavefront reconstruction and DM fitting.
*/

/*
   Since these are related to reconstruction, we don't have access to dmreal,
   which is the *actual* location of DM actuators, and only available in
   simulation. dmint should be used for dm actuator commands.*/

/**
   Calls tomo() and fit() to do the tomography and DM fit. Do error signal and
   split tomography.
   In closedloop mode, the gradients are from time step isim-1.
*/
void tomofit(dcell** dmout, sim_t* simu, dcell* gradin){
	const parms_t* parms=simu->parms;
	recon_t* recon=simu->recon;
	int isim=simu->reconisim;

	if(parms->sim.idealtomo){
		dcellfree(simu->opdr);
	} else{	/*do tomography. */
		int maxit=parms->tomo.maxit;
		if(NX(parms->dbg.tomo_maxit)){
			if(isim<NX(parms->dbg.tomo_maxit)){
				maxit=P(parms->dbg.tomo_maxit,isim);
				recon->RL.maxit=maxit;/*update maxit information */
				info2("Running tomo.maxit=%d\n", maxit);
			} else{
				error("Out of range\n");
			}
		}
		TIC_tm; tic_tm;
#if USE_CUDA
		if(parms->gpu.tomo){
			gpu_tomo(simu, gradin);
		} else
#endif
			P(P(simu->cgres,0),isim)=muv_solve(&simu->opdr, &recon->RL, &recon->RR, gradin);
		toc_tm("Tomography");
	}
	if(parms->ndm>0){//we still do DM fitting in evl.tomo so that WFS works in closed loop.
		TIC_tm; tic_tm;
#if USE_CUDA
		if(parms->gpu.fit){
			gpu_fit(dmout, simu);
		} else
#endif
		{
			P(P(simu->cgres,1),isim)=muv_solve(dmout, &recon->fit->FL, &recon->fit->FR, simu->opdr);
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
static void calc_gradol(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].psol){
			if((simu->reconisim+1)%parms->powfs[ipowfs].dtrat==0){/*Has output. */
				int nindwfs=parms->recon.glao?1:parms->powfs[ipowfs].nwfs;
OMP_TASK_FOR(4)
				for(long indwfs=0; indwfs<nindwfs; indwfs++){			
					int iwfs=parms->recon.glao?ipowfs:P(parms->powfs[ipowfs].wfs,indwfs);
					dcp(&P(simu->gradlastol,iwfs), P(simu->gradlastcl,iwfs));
					dcell *dmpsol=P(simu->wfspsol,ipowfs);
					for(int idm=0; idm<parms->ndm&&dmpsol; idm++){
						dcellmm(&P(simu->gradlastol,iwfs), P(recon->GA, iwfs, idm),
							P(dmpsol,idm), "nn", 1);
					}
				}
			}
		}
	}
	/*{
		dcell *junk=NULL;
		dcellmm(&junk, recon->PTTF, simu->gradlastcl, "nn", 1);
		dbg("gradcl TTF: %g %g %g\n", P(P(junk, 0, 0), 0), P(P(junk, 0, 0), 1), P(P(junk, 0, 0), 2));
		dcellzero(junk);
		dcellmm(&junk, recon->PTTF, simu->gradlastol, "nn", 1);
		dbg("gradol TTF: %g %g %g\n", P(P(junk, 0, 0), 0), P(P(junk, 0, 0), 1), P(P(junk, 0, 0), 2));
		cellfree(junk);
	}*/
	save_gradol(simu);//must be here since gradol is only calculated in this file. 
}
static void recon_split_lo(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int isim=simu->reconisim;
	if(parms->recon.split==2){
		if(!parms->gpu.tomo){
			dcellmm(&simu->gngsmvst, recon->GXL, simu->opdr, "nn", 1./parms->sim.dtrat_lo);
		}
	}
	int enRngs[2]={0,0}; int anyRngs=0;
	if(parms->ntipowfs&&isim>=parms->step_lo){
		if(!parms->sim.closeloop||(isim+1)%parms->sim.dtrat_lo==0){
			enRngs[0]=1;//Common rate
			anyRngs++;
		}
		if(parms->sim.dtrat_lo!=parms->sim.dtrat_lo2&&(isim+1)%parms->sim.dtrat_lo2==0){
			enRngs[1]=1;// Multi-rate control

			anyRngs++;
			if(parms->recon.split==2){
				error("Multi-rate control for MVR is to be implemented\n");
			}
		}
	}
	if(anyRngs){/*Low order WFS has output */
		simu->Merr_lo=simu->Merr_lo_store;
		dcellzero(simu->Merr_lo);

		switch(parms->recon.split){
		case 1:	{//Low order NGS recon.
			dcell* tmp=0;
			ngsmod_t* ngsmod=recon->ngsmod;
			for(int iRngs=0; iRngs<2; iRngs++){
				//For multi-rate control, iRngs=0 is slower loop and iRngs=1 is the faster loop
				if(!enRngs[iRngs]) continue;
				dcell** merr;//reconstruction output
				if((ngsmod->lp2>=0&&iRngs==1)||(ngsmod->lp2<0&&iRngs==0)){
					merr=&tmp; //output to separate array to handle LPF or TWFS
					dcellzero(tmp);
				} else{
					merr=&simu->Merr_lo;
				}
				dcellmm(merr, P(ngsmod->Rngs,iRngs), simu->gradlastcl, "nn", 1);
				//dshow(P(*merr,0), "merr");
				if(iRngs==0){
					dcellscale(*merr, parms->dbg.eploscale);
				}
				if(iRngs==1&&ngsmod->lp2>=0){ //Do LHF on measurements
					if(ngsmod->lp2>0){//HPF
						real* valpf=P(P(simu->Merr_lo2,0));
						real* val=P(P(tmp,0));
						for(int imod=0; imod<ngsmod->nmod; imod++){
							if(imod==ngsmod->indfocus){//there is no need to blend focus.
								continue;
							}
							valpf[imod]=valpf[imod]*(1.-ngsmod->lp2)+val[imod]*ngsmod->lp2;
							val[imod]-=valpf[imod];
						}
					}
					if(ngsmod->lp2>0||(ngsmod->lp2==0&&!enRngs[0])){
						dcelladd(&simu->Merr_lo, 1, tmp, 1);
					}
				} else if(ngsmod->lp2<0){//Use slower as Truth WFS mode by mode
					for(int imod=0; imod<ngsmod->nmod; imod++){
						if(P(ngsmod->modvalid, imod)){//Modes that has multi-rates
							if(iRngs==0){//Accumulate Truth mode offset
								P(P(simu->Merr_lo2,0),imod)+=P(P(tmp,0),imod)*(0.5/parms->dbg.eploscale);
							} else{//Apply truth mode offset.
								P(P(simu->Merr_lo,0),imod)+=P(P(simu->Merr_lo2,0),imod);
							}
						} else if(iRngs==0){//direct output. Avoid real integrator as above.
							P(P(simu->Merr_lo,0),imod)=P(P(tmp,0),imod);
						}
					}
				}
				//dshow(P(*merr, 0), "merr");
			}//for iRngs
			dcellfree(tmp);
			if(parms->sim.mffocus&&ngsmod->indfocus && parms->sim.lpfocushi<1&& P(P(simu->Merr_lo,0), ngsmod->indfocus)){
				if(parms->sim.dtrat_lo>1){//new method. use NGS measurement as LGS offset
					for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
						const int ipowfs=parms->wfs[iwfs].powfs;
						if(isfinite(parms->powfs[ipowfs].hs)||parms->powfs[ipowfs].llt){
							const real ep=P(parms->sim.eplo, 0);
							dadd(&P(simu->gradoff, iwfs), 1, P(recon->GFall, iwfs), -ep*P(P(simu->Merr_lo,0), ngsmod->indfocus));
						}
					}
					P(P(simu->Merr_lo,0), ngsmod->indfocus)=0;//remove from low order path
				}else{ //Do LPF on NGS focus.
					//This does not work when dtrat_lo>1. Need to upsample to LGS WFS or DM rate
					const real lpfocus=parms->sim.lpfocuslo;
					real *ngsfocus=&P(P(simu->Merr_lo,0), ngsmod->indfocus);
					simu->ngsfocuslpf=simu->ngsfocuslpf*(1-lpfocus)+lpfocus* *ngsfocus;
					*ngsfocus=simu->ngsfocuslpf;
				}
			}
		}//else: there is ideal NGS correction done in perfevl. 
		break;
		case 2:{
			/*A separate integrator for low order is required. Use it to form error signal*/
			dcelladd(&simu->gradlastol, 1, simu->gngsmvst, -1);
			dcellzero(simu->gngsmvst);/*reset accumulation. */
			dcellmm(&simu->Merr_lo, recon->MVRngs, simu->gradlastol, "nn", 1);
			if(parms->recon.psol){
				dcell* Mpsol_lo=P(simu->Mint_lo->mintc,0);
				dcelladd(&simu->Merr_lo, 1., Mpsol_lo, -1);
			}
			if(parms->sim.mffocus){
				dcell* tmp=NULL;
				dcellmm(&tmp, recon->RFngsg, simu->gradlastcl, "nn", 1);
				dcellmm(&tmp, recon->MVFM, simu->Merr_lo, "nn", -1);
				const real lpfocus=parms->sim.lpfocuslo;
				real ngsfocus=P(P(tmp,0),0);
				simu->ngsfocuslpf=simu->ngsfocuslpf*(1-lpfocus)+lpfocus*ngsfocus;
				error("Please Implement: add ngsfocus to Merr_lo");
				dcellfree(tmp);
			}
		}
		break;
		default:
			error("Invalid parms->recon.split: %d\n", parms->recon.split);
		}
	}
}
/**
 * @brief Combine separate 6 NGS modes to TT, PS, Focus
 * 
 * @param psdall 
 * @return void* 
 */
dmat *combine_ngs_psd(const dmat *psdall){
	int ncol=NY(psdall)-1;
	if(ncol!=6 && ncol!=5){
		return dref(psdall);
	}
	dmat *psd=dnew(NX(psdall), ncol==6?4:3);
	for(int ix=0; ix<NX(psd); ix++){
		P(psd, ix, 0)=P(psdall, ix, 0);
		P(psd, ix, 1)=P(psdall, ix, 1)+P(psdall, ix, 2);//TT
		P(psd, ix, 2)=P(psdall, ix, 3)+P(psdall, ix, 4)+P(psdall, ix, 5);//PS
		if(ncol==6) P(psd, ix, 3)=P(psdall, ix, 6);//Focus
	}
	return psd;
}
/**
 * Convert each PSD mode into open loop PSD and plot both the PSDs and reverse cumulative integral of psd open loop
 */
static void plot_psd(const dmat *psdall, double dt, int dtrat, 
	double al, const dmat *ep, const char *label, const char*modes[], int doplot, zfarr *save_cl, zfarr *save_ol){
	if(!ep || NX(ep)!=1){
		error("Please implement\n");
		return;
	}
	int nmod=NY(psdall)-1;
	dcell *psds=dcellnew(2*nmod, 1);
	dcell *cums=dcellnew(nmod, 1);//only integrate OL PSD to indicate residual with correction bandwidth
	char **legs=mycalloc(nmod*2, char*);
	int32_t *style=mycalloc(nmod*2, int32_t);
	
	for(int imod=0; imod<nmod; imod++){
		dmat *psd=psd_select(psdall, imod, 1, 0, 1);
		dmat *psdol=servo_cl2ol(psd, dt, dtrat, al, PR(ep, 0, imod), 0);
		if(save_cl) zfarr_push(save_cl, -1, psd);
		if(save_ol) zfarr_push(save_ol, -1, psdol);
		P(psds, imod)=dref(psd);//dsub(psd,  1,-1,0,-1); 
		P(psds, imod+nmod)=dref(psdol);//dsub(psdol, 1, -1, 0, -1);
		P(cums, imod)=psd_reverse_cumu(psdol, 1e9);
		dfree(psdol);
		dfree(psd);
		char temp[1024]; 
		snprintf(temp, sizeof(temp), "Closed Loop %s", modes?modes[imod]:""); legs[imod]=strdup(temp);
		snprintf(temp, sizeof(temp), "Open Loop %s", modes?modes[imod]:""); legs[imod+nmod]=strdup(temp);
		style[imod]=default_color(imod)<<8;   //solid  line
		style[imod+nmod]=default_color(imod)<<8|7; //dashed line
	}
	if(doplot){
		draw("PSD", (plot_opts){.dc=psds, .style=(nmod==1)?NULL:style, .xylog="yy", .legend=(const char *const *)legs, .always=1}, "", "Frequency (Hz)", "PSD (m<sup>2</sup>/Hz)", "PSD %s", label);
		draw("PSD", (plot_opts){.dc=cums, .style=(nmod==1)?NULL:style, .xylog="yy", .legend=(const char *const *)(legs+nmod), .always=1 }, "", "Frequency (Hz)", "PSD Reverse Cumulative Integral (nm)", "Cumu %s", label);
	}
	dcellfree(psds);
	dcellfree(cums);
	for(int imod=0; imod<2*nmod; imod++){
		free(legs[imod]);
	}
	free(legs);
	free(style);
}
/**
 * Servo gain optimization.
  * */
void recon_servo_update(sim_t* simu){
	const parms_t* parms=simu->parms;
	recon_t* recon=simu->recon;
	if(!parms->recon.psd) return;
	if(simu->dmerr&&parms->recon.psddtrat_hi>0){//compute PSD on dmerr.
		const int dtrat=parms->recon.psddtrat_hi;
		const int iacc=(simu->reconisim/parms->sim.dtrat_hi);//reconstruction steps
		const int iframe=iacc%dtrat;
		//Accumulate data history.
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
			for(int idm=0; idm<parms->ndm; idm++){
				dmat *out=drefcols(P(simu->dmerrts, ievl), iframe,1);
				dcellmm((cell**)&out, P(recon->Herr, ievl, idm), P(simu->dmerr,idm), "nn", 1);
				dfree(out);
			}
		}
		if(iframe+1==dtrat){//ready for output.
			dmat* psdall=0;
			real dthi=parms->sim.dt*parms->sim.dtrat_hi;
			for(int ievl=0; ievl<parms->evl.nevl; ievl++){
				dmat* tmp=dtrans(P(simu->dmerrts, ievl));
				dmat* psdi=psd1dt(tmp, parms->recon.psdnseg, dthi);
				dfree(tmp);
				dadd(&psdall, 1, psdi, P(parms->evl.wt,ievl));
				dfree(psdi);
			}
			dcellzero(simu->dmerrts);
			dmat *psd=psd_select(psdall, -1, 1, 0, 1./(NY(psdall)-1));//average all points
			dfree(psdall);
			//plot psd ol and cl before updating the gain
			plot_psd(psd, parms->sim.dt, parms->sim.dtrat_hi, parms->sim.alhi, simu->dmint->ep, "High", NULL,
				parms->plot.run, simu->save->psdcl, simu->save->psdol);
			
			if(NX(simu->dmint->ep)==1&&NY(simu->dmint->ep)==1){
				dmat* psdol=servo_cl2ol(psd, parms->sim.dt, parms->sim.dtrat_hi, parms->sim.alhi, P(simu->dmint->ep,0), 0);
				dcell* coeff=servo_optim(parms->sim.dt, parms->sim.dtrat_hi, parms->sim.alhi, M_PI*0.25,
					parms->sim.f0dm, parms->sim.zetadm, 1, psdol, NULL);
				const real g=parms->recon.psdservo_gain;
				const real oldep=P(simu->dmint->ep, 0);
				P(simu->dmint->ep,0)=oldep*(1-g)+P(P(coeff,0),0)*g;
				info("Step %5d updated HO loop gain: %5.3f->%5.3f (%ld points)\n", 
					simu->reconisim, oldep, P(simu->dmint->ep,0), NY(P(simu->dmerrts,0)));
				if(simu->save->gain) P(P(simu->save->gain,0),(iacc/dtrat))=P(simu->dmint->ep,0);
				
				dcellfree(coeff);
				dfree(psdol);
			} else{
				error("Please implement\n");
			}
			dfree(psd);
		}
	}
	if((simu->Merr_lo||parms->evl.split)&&parms->recon.psddtrat_lo>0){//compute PSD on low order control
		const int iacc=(simu->reconisim/parms->sim.dtrat_lo);//reconstruction steps
		const int dtrat=parms->recon.psddtrat_lo;
		const int iframe=iacc%dtrat;
		if(simu->Merr_lo){
			dmulvec(PCOL(simu->Merrts, iframe), recon->ngsmod->MCCu, P(P(simu->Merr_lo,0)), 1);
		}else if(!parms->recon.split && simu->dmerr){
			dcellzero(simu->Mngs);
			dcellmm(&simu->Mngs, recon->ngsmod->Pngs, simu->dmerr, "nn", 1);
			dmulvec(PCOL(simu->Merrts, iframe), recon->ngsmod->MCCu, P(P(simu->Mngs, 0)), 1);
		}
		if(iframe+1==dtrat){
			/*if(parms->save.dm){
				writebin(simu->Merrts, "Merrts_%d", simu->reconisim);
			}*/
			dmat* ts=dtrans(simu->Merrts);
			dzero(simu->Merrts);
			real dt=parms->sim.dt*parms->sim.dtrat_lo;
			dmat* psdall=psd1dt(ts, parms->recon.psdnseg, dt);
			{//plot psd ol and cl before updating the gain
				dmat *psdngs=combine_ngs_psd(psdall);
				const char *modes1[]={"Tip Tilt","Plate Scale","Focus"};
				const char *modes2[]={"Tip","Tilt","PS1","PS2","PS3","Focus"};
				const char **modes=NULL;
				if(NY(psdngs)==6||NY(psdngs)==7){
					modes=modes2;
				}else if(NY(psdngs)<=4){
					modes=modes1;
				}
				dmat *ep=simu->Mint_lo?simu->Mint_lo->ep:simu->dmint->ep;
				plot_psd(psdngs, parms->sim.dt, parms->sim.dtrat_lo, parms->sim.allo, ep, "Low", modes,
					parms->plot.run, simu->save->psdcl_lo, simu->save->psdol_lo);
				dfree(psdngs);
			}
			if(parms->recon.split){
				for(int icol=0; icol<NY(simu->Mint_lo->ep); icol++){
					dmat *psd=psd_select(psdall, NY(simu->Mint_lo->ep)==1?-1:icol, 1, 0, 1);
					//if(simu->save->psdcl_lo) zfarr_push(simu->save->psdcl_lo, -1, psd);
					if(NX(simu->Mint_lo->ep)==1){//integrator
						dmat* psdol=servo_cl2ol(psd, parms->sim.dt, parms->sim.dtrat_lo, parms->sim.allo, P(simu->Mint_lo->ep,0,icol), 0);
						//if(simu->save->psdol_lo) zfarr_push(simu->save->psdol_lo, -1, psdol);
						dcell *coeff=servo_optim(parms->sim.dt, parms->sim.dtrat_lo, parms->sim.allo, M_PI*0.25, 0, 0, 1, psdol, 0);
						const real g=parms->recon.psdservo_gain;
						const real oldep=P(simu->Mint_lo->ep, 0, icol);
						P(simu->Mint_lo->ep,0,icol)=oldep*(1.-g)+P(P(coeff,0),0)*g;
						info("Step %5d updated LO loop gain: %5.3f->%5.3f (%ld points)\n", simu->reconisim, oldep, P(simu->Mint_lo->ep,0,icol), NY(simu->Merrts));
						if(icol==0 && simu->save->gain) P(P(simu->save->gain,1+icol),(iacc/dtrat))=P(simu->Mint_lo->ep,0,icol);
						dfree(psdol);
						dcellfree(coeff);
					} else{
						error("Please implement\n");
					}
					dfree(psd);
				}
			}
			dfree(psdall);
			dfree(ts);
		}
	}
}
/**
   Wavefront reconstruction. call tomofit() to do tomo()/fit() or lsr() to do
   least square reconstruction. */
void* reconstruct(sim_t* simu){
	real tk_start=myclockd();
	const parms_t* parms=simu->parms;
	recon_t* recon=simu->recon;
	int isim=simu->reconisim;
	if(parms->sim.evlol || isim<0) {
		//dbg("reconstruct: evlol=%d, isim=%d. return\n",parms->sim.evlol, isim);
		return NULL;
	}
	if(PARALLEL==2){
		pthread_mutex_lock(&simu->wfsgrad_mutex);
		while(simu->wfsgrad_isim<simu->reconisim){
			//dbg("waiting wfsgrad_isim is %d need %d\n", simu->wfsgrad_isim, simu->reconisim);
			//if(simu->wfsgrad_isim+1==simu->reconisim){
			//}
			struct timespec ts;
			clock_gettime(CLOCK_REALTIME, &ts);
			ts.tv_nsec+=1e6;
			pthread_cond_timedwait(&simu->wfsgrad_condr, &simu->wfsgrad_mutex, &ts);
		}
		if(simu->wfsgrad_isim>simu->reconisim){
			error("waiting wfsgrad_isim is %d need %d\n", simu->wfsgrad_isim, simu->reconisim);
		}
		//dbg("ready: wfsgrad_isim is %d need %d\n", simu->wfsgrad_isim, simu->reconisim);
		pthread_mutex_unlock(&simu->wfsgrad_mutex);
	}
	const int hi_output=(!parms->sim.closeloop||((isim+1>parms->step_hi)&&(isim+1-parms->step_hi)%parms->sim.dtrat_hi==0));
	if(hi_output){
		if(simu->gradlastcl){
			if(parms->sim.closeloop){
				calc_gradol(simu);
			}
			if(recon->cn2est){
				cn2est_isim(simu->cn2res, recon, parms, parms->cn2.psol?simu->gradlastol:simu->gradlastcl, &simu->tomo_update);
			}//if cn2est 
		}
	
		simu->dmerr=simu->dmerr_store;
		dcell* dmout=simu->dmrecon;//always output to dmrecon to enable warm restart.
		dcell* gradin;
		if(parms->recon.psol){
			gradin=simu->gradlastol;
		} else{
			gradin=simu->gradlastcl;
		}
		if(!dmout) error("dmout cannot be empty\n");
		//The following takes gradin as input and computs dmrecon in dmout.
		if(parms->recon.mvm){
			if(parms->sim.mvmport){
				mvm_client_recon(parms->sim.mvmsize, dmout, gradin);
			} else
#if USE_CUDA
				if((parms->gpu.tomo&&parms->gpu.fit)||parms->gpu.lsr){
					gpu_recon_mvm(&dmout, gradin);
				} else
#endif		
				{
					dzero(dmout->m);
					dmulvec(P(dmout->m), recon->MVM, P(gradin->m), 1);
				}
		} else{
			switch(parms->recon.alg){
			case 0://MVR
				tomofit(&dmout, simu, gradin);//tomography and fitting. 
				break;
			case 1://LSR
				if(simu->gradlastcl){
#if USE_CUDA
					if(parms->gpu.lsr){
						warning_once("Not implemented. Use CPU instead\n");
					}
#endif
					muv_solve(&dmout, &(recon->LL), &(recon->LR), gradin);
				}
				break;
			default:
				error("recon.alg=%d is not recognized\n", parms->recon.alg);
			}
		}
		if(simu->dmrecon!=simu->dmerr){
			dcellcp(&simu->dmerr, simu->dmrecon);/*keep dmrecon for warm restart */
		}
		if(parms->recon.psol){
			//form error signal in PSOL mode
			if(simu->recon->actextrap){
				//extrapolate DM fitting result to float and edge actuators.
				//If not enabled, extrapolate integrator output.
				//Must be enabled if HA is altered by actextrap.
				dcell *dmtmp;
				if(simu->dmrecon!=simu->dmerr){
					dmtmp=simu->dmrecon;
				}else{
					dcellcp(&simu->dmtmp, simu->dmerr);
					dmtmp=simu->dmtmp;
				}
				dcellzero(simu->dmerr);
				dcellmm(&simu->dmerr, simu->recon->actextrap, dmtmp, "nn", 1);
			}

			dcell* dmpsol;
			if(parms->sim.idealtomo){
				dmpsol=simu->dmpsol;
			} else if(parms->sim.fuseint||parms->recon.split==1){
				dmpsol=P(simu->wfspsol,P(parms->hipowfs,0));
			} else{
				warning_once("Temporary solution for MVST.\n");
				dmpsol=P(simu->dmint->mintc,0);
			}
			dcelladd(&simu->dmerr, 1, dmpsol, -1);
		}else if(parms->recon.modal){//multi-mode dithering in modal reconstruction
			const int iwfs=0;
			const int ipowfs=parms->wfs[iwfs].powfs;
			if(parms->powfs[ipowfs].dither>1){
				//multi-mode dithering for PWFS
				if(parms->powfs[ipowfs].type==1&&simu->gradscale2&&P(simu->gradscale2, iwfs)){
					//const double gscale=P(P(simu->gradscale, iwfs),0);
					int print=simu->wfsflags[ipowfs].ogout?1:0;
					for(int idm=0; idm<parms->ndm; idm++){
						//long nmod=parms->powfs[ipowfs].dither_mode2>1?parms->powfs[ipowfs].dither_mode2:PN(simu->dmrecon, idm);
						const dmat *gs2=P(simu->gradscale2, iwfs);					
						const long nd=NX(gs2);//number of dithered modes.
						const real gs1=P(gs2,0);//do not perturb the gain of low order modes.
						//const real gs1=sqrt(P(gs2, 0)*P(gs2, PN(gs2)-1));//scaling baseline: don't use gradscale, it is computed differently from gradscale2.
						//const real gs1=P(gs2, PN(gs2)-1); //testing. when main dithering integrator uses the last loop.
						const int md=recon->dither_md;
						
						if(print){
							info("modal scaling (%d,%ld): ", md, nd);
						}
						for(int id=0; id<nd; id++){//index of dithered mode
							const int jm=md*id;//DM mode of dithered mode. 
							const int jm2=(id+1==nd&&recon->anmod)?P(recon->anmod, idm):(jm+md);
							const real scale0=P(gs2, id)/gs1;
							const real dscale1=scale0-P(gs2, MIN(nd-1, id+1))/gs1;
							if(print) info(" %4.2f", scale0);
							for(int imod=jm; imod<jm2; imod++){
								P(P(simu->dmrecon, idm), imod)*=(scale0-dscale1*(imod-jm)/(jm2-jm));
							}
						}
						if(print){
							info("\n");
						}
					}
				}
			}
		}
		if(parms->recon.split){
			ngsmod_remove(simu, simu->dmerr);
		}else if(parms->evl.split>1){//split ngs mode error from dmerr
			simu->Merr_lo=simu->Merr_lo_store;
			ngsmod_split(&simu->Merr_lo, simu, simu->dmerr);
		}
	}
	if(PARALLEL==2){
		//dbg("ready: wfsgrad_isim is %d need %d\n", simu->wfsgrad_isim, simu->reconisim);
		atomic_add_fetch(&simu->wfsgrad_count,1);
		pthread_cond_broadcast(&simu->wfsgrad_condw);
	}
	if(parms->recon.split && !parms->tomo.ahst_idealngs){//low order reconstruction
		recon_split_lo(simu);
	}
	if(parms->recon.psd&&parms->sim.closeloop){
		recon_servo_update(simu);
	}
	if(hi_output&&parms->save.ecov&&isim>=parms->evl.psfisim){
		if(!parms->recon.psol&&parms->sim.closeloop){
			error("Please enable PSOL\n");
		}
		//For PSF reconstruction.
		psfr_calc(simu, simu->opdr, P(simu->wfspsol,P(parms->hipowfs,0)),
			simu->dmerr, simu->Merr_lo);
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
	return NULL;
}
