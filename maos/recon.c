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

	if(parms->sim.idealfit){
		dcellfree(simu->opdr);
	} else if(parms->sim.idealtomo){
		atm2xloc(&simu->opdr, simu);
	} else{
	/*do tomography. */
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
		if(parms->gpu.tomo&&parms->ndm!=0){
			gpu_tomo(simu, gradin);
		} else
#endif
			P(P(simu->cgres,0),isim)=muv_solve(&simu->opdr, &recon->RL, &recon->RR, gradin);
		toc_tm("Tomography");
	}
	if(parms->ndm>0){
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
	dspcell* GA=recon->GA/*PDSPCELL*/;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(parms->powfs[ipowfs].psol){
			if((simu->reconisim+1)%parms->powfs[ipowfs].dtrat==0){/*Has output. */
				int nindwfs=parms->recon.glao?1:parms->powfs[ipowfs].nwfs;
OMP_TASK_FOR(4)
				for(long indwfs=0; indwfs<nindwfs; indwfs++){			
					int iwfs=parms->recon.glao?ipowfs:P(parms->powfs[ipowfs].wfs,indwfs);
					dcp(&P(simu->gradlastol,iwfs), P(simu->gradlastcl,iwfs));
					for(int idm=0; idm<parms->ndm&&P(simu->wfspsol,ipowfs); idm++){
						dspmm(&P(simu->gradlastol,iwfs), P(GA, iwfs, idm),
							P(P(simu->wfspsol,ipowfs),idm), "nn", 1);
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
}
static void recon_split(sim_t* simu){
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
	if(anyRngs){
	/*Low order WFS has output */
		simu->Merr_lo=simu->Merr_lo_store;
		dcellzero(simu->Merr_lo);

		switch(parms->recon.split){
		case 1:
			if(!parms->tomo.ahst_idealngs){//Low order NGS recon.
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
				}//for iRngs
				dcellfree(tmp);

				if(parms->sim.mffocus&&ngsmod->indfocus&&parms->sim.lpfocushi<1){ //Do LPF on focus.
					const real lpfocus=parms->sim.lpfocuslo;
					real ngsfocus=P(P(simu->Merr_lo,0),ngsmod->indfocus);
					if(ngsfocus){//there is output
						simu->ngsfocuslpf=simu->ngsfocuslpf*(1-lpfocus)+lpfocus*ngsfocus;
						P(P(simu->Merr_lo,0),ngsmod->indfocus)=simu->ngsfocuslpf;
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
				dcell* Mpsol_lo=P(simu->Mint_lo->mint,0);
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
				dspmulvec(PCOL(P(simu->dmerrts, ievl), iframe), P(recon->Herr, ievl, idm),
					P(P(simu->dmerr,idm)), 'n', 1);
			}
		}

		if(iframe+1==dtrat){//ready for output.
			dmat* psd=0;
			real dthi=parms->sim.dt*parms->sim.dtrat_hi;
			for(int ievl=0; ievl<parms->evl.nevl; ievl++){
				dmat* tmp=dtrans(P(simu->dmerrts, ievl));
				dmat* psdi=psd1dt(tmp, parms->recon.psdnseg, dthi);
				dfree(tmp);
				dadd(&psd, 1, psdi, P(parms->evl.wt,ievl));
				dfree(psdi);
			}
			dcellzero(simu->dmerrts);
			//writebin(simu->dmerrts, "dmerrts_%d", simu->reconisim);
			//writebin(psd, "psdcli_%d", simu->reconisim);
			//average all the PSDs
			psd_sum(psd, 1./(psd->ny-1));
			if(simu->save->psdcl) zfarr_push(simu->save->psdcl, -1, psd);
			//writebin(psd, "psdcl_%d_%d", simu->iseed, simu->reconisim);
			if(NX(simu->dmint->ep)==1&&NY(simu->dmint->ep)==1){
				dmat* psdol=servo_rej2ol(psd, parms->sim.dt, parms->sim.dtrat_hi, parms->sim.alhi, P(simu->dmint->ep,0), 0);
				dcell* coeff=servo_optim(parms->sim.dt, parms->sim.dtrat_hi, parms->sim.alhi, M_PI*0.25,
					parms->sim.f0dm, parms->sim.zetadm, 1, psdol, NULL);
				real g=0.5;
				P(simu->dmint->ep,0)=P(simu->dmint->ep,0)*(1-g)+P(P(coeff,0),0)*g;
				info("Step %5d New gain (high): %.3f\n", simu->reconisim, P(simu->dmint->ep,0));
				//if(parms->save.run){
				//writebin(psdol, "psdol_%d_%d", simu->iseed, simu->reconisim);
				if(simu->save->psdol) zfarr_push(simu->save->psdol, -1, psdol);
				//}
				dcellfree(coeff);
				dfree(psdol);
			} else{
				error("Please implement\n");
			}
			dfree(psd);
		}
	}
	if(parms->recon.split&&simu->Merr_lo&&parms->recon.psddtrat_lo>0){//compute PSD on low order control
		const int iacc=(simu->reconisim/parms->sim.dtrat_lo);//reconstruction steps
		const int dtrat=parms->recon.psddtrat_lo;
		const int iframe=iacc%dtrat;
		dmulvec(PCOL(simu->Merrts, iframe), recon->ngsmod->MCCu, P(P(simu->Merr_lo,0)), 1);
		if(iframe+1==dtrat){
			//writebin(simu->Merrts, "Merrts_%d", simu->reconisim);
			dmat* ts=dtrans(simu->Merrts);
			dzero(simu->Merrts);
			real dt=parms->sim.dt*parms->sim.dtrat_lo;
			for(int icol=0; icol<NY(ts); icol++){
				dmat* tsi=dsub(ts, icol, 1, 0, 0);
				dmat* psd=psd1dt(tsi, parms->recon.psdnseg, dt);
				if(simu->save->psdcl_lo) zfarr_push(simu->save->psdcl_lo, -1, psd);
				//writebin(psd, "psdlo%d_cl_%d", icol, simu->reconisim);
				if(NX(simu->Mint_lo->ep)==1){//integrator
					dmat* psdol=servo_rej2ol(psd, parms->sim.dt, parms->sim.dtrat_lo, parms->sim.allo, P(simu->Mint_lo->ep,0), 0);
					//writebin(psdol, "psdlo%d_ol_%d", icol, simu->reconisim);
					if(simu->save->psdol_lo) zfarr_push(simu->save->psdol_lo, -1, psdol);
					dcell *coeff=servo_optim(parms->sim.dt, parms->sim.dtrat_lo, parms->sim.allo, M_PI*0.25, 0, 0, 1, psdol, 0);
					const real g=parms->recon.psdservo_gain;
					P(simu->Mint_lo->ep,0)=P(simu->Mint_lo->ep,0)*(1-g)+P(P(coeff,0),0)*g;
					if(icol==0) dbg("Step %5d New gain (low) : %.3f\n", simu->reconisim, P(simu->Mint_lo->ep,0));
					dfree(psdol);
					dcellfree(coeff);
				} else{
					error("Please implement\n");
				}
				dfree(psd);
			}
			dfree(ts);
		}
	}
}
/**
   Wavefront reconstruction. call tomofit() to do tomo()/fit() or lsr() to do
   least square reconstruction. */
void reconstruct(sim_t* simu){
	real tk_start=myclockd();
	const parms_t* parms=simu->parms;
	if(parms->sim.evlol) return;
	recon_t* recon=simu->recon;
	int isim=simu->reconisim;
	if(isim<0) return;
	if(PARALLEL==2){
		while(simu->wfsgrad_isim<simu->reconisim){
			//dbg("waiting wfsgrad_isim is %d need %d\n", simu->wfsgrad_isim, simu->reconisim);
			//if(simu->wfsgrad_isim+1==simu->reconisim){
			//}
			pthread_cond_wait(&simu->wfsgrad_condr, &simu->wfsgrad_mutex);
		}
		//dbg("ready: wfsgrad_isim is %d need %d\n", simu->wfsgrad_isim, simu->reconisim);
		simu->wfsgrad_count++;
		pthread_mutex_unlock(&simu->wfsgrad_mutex);
	}
	const int hi_output=(!parms->sim.closeloop||(isim+1-parms->step_hi)%parms->sim.dtrat_hi==0);
	if(simu->gradlastcl){
		if(parms->sim.closeloop){
			calc_gradol(simu);
			save_gradol(simu);//must be here since gradol is only calculated in this file. 
		}
		if(recon->cn2est){
			cn2est_isim(simu->cn2res, recon, parms, parms->cn2.psol?simu->gradlastol:simu->gradlastcl, &simu->tomo_update);
		}//if cn2est 
	}
	if(hi_output||parms->sim.idealfit||parms->sim.idealtomo){
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
			if(parms->recon.alg==0 && parms->fit.actextrap){
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
				dspcellmm(&simu->dmerr, simu->recon->actextrap, dmtmp, "nn", 1);
			}

			dcell* dmpsol;
			if(parms->sim.idealfit||parms->sim.idealtomo){
				dmpsol=simu->dmpsol;
			} else if(parms->sim.fuseint||parms->recon.split==1){
				dmpsol=P(simu->wfspsol,P(parms->hipowfs,0));
			} else{
				warning_once("Temporary solution for MVST.\n");
				dmpsol=P(simu->dmint->mint,0);
			}
			dcelladd(&simu->dmerr, 1, dmpsol, -1);
		}
		if(parms->recon.split){
			remove_dm_ngsmod(simu, simu->dmerr);
		}

	}

	if(parms->recon.split){//low order reconstruction
		recon_split(simu);
	}
	if(parms->recon.psd&&parms->sim.closeloop){
		recon_servo_update(simu);
	}
	if(hi_output&&parms->save.ecov&&isim>=parms->evl.psfisim){
		if(!parms->recon.psol){
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
}
