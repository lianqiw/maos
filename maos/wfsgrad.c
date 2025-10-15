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
#include "sim.h"
#include "sim_utils.h"
#include "ahst.h"
#include "powfs_utils.h"
#include "save.h"
//#include "setup_recon.h"
#include "recon_utils.h"
#include "powfs.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   contains functions that computes WFS gradients in geometric or physical optics mode.
*/
/*#define TIMING 1
#if TIMING
#define TIM(A) real tk##A=myclockd()
#else
#define TIM(A)
#endif*/
#define TIMING 0
#if TIMING
#define TIM0 static real tk1=0,tk2=0,tk3=0,tk4=0;static int tkct=0;real tk=0,tk0=myclockd();tkct++;
#define TIM(A) tk=myclockd(); tk##A+=tk-tk0;tk0=tk;
#define TIM1 info2("wfsgrad timing: atm %.3f dm %.3f ints %.3f grad %.3f tot %.3f.\n", tk1/tkct, tk2/tkct, tk3/tkct, tk4/tkct, (tk1+tk2+tk3+tk4)/tkct)
#else
#define TIM0
#define TIM(A)
#define TIM1
#endif
/**
   Propagate only controllable component of turbulence to WFS.
*/
void wfs_ideal_atm(sim_t* simu, dmat* opd, int iwfs, real alpha){
	const parms_t* parms=simu->parms;
	const powfs_t* powfs=simu->powfs;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=P(parms->powfs[ipowfs].wfsind, iwfs);
	const real hs=parms->wfs[iwfs].hs;
	const real misregx=parms->powfs[ipowfs].type==WFS_SH?parms->wfs[iwfs].misregx:0;
	const real misregy=parms->powfs[ipowfs].type==WFS_SH?parms->wfs[iwfs].misregy:0;
	//hc is only useful for multi-sublayer raytracing
	if(parms->sim.wfsalias==2||parms->sim.idealwfs==2){
		//project turbulence onto lenslet array grid (per direction ideal)
		loc_t* aloc=P(powfs[ipowfs].fit[wfsind].aloc, 0);
		dcell* wfsopd=dcellnew(1, 1); P(wfsopd, 0)=dnew(aloc->nloc, 1);
		fit_t* fit=&powfs[ipowfs].fit[wfsind];
		muv_solve(&wfsopd, &fit->FL, &fit->FR, 0);
		loc_t *loc=powfs[ipowfs].loc;
		if(powfs[ipowfs].loc_tel){
			loc=P(powfs[ipowfs].loc_tel, wfsind);
		}
		prop_nongrid(aloc, P(P(wfsopd, 0)), loc, P(opd), alpha, misregx, misregy, 1, 0, 0);
		dcellfree(wfsopd);
	} else{
		//project tubulence onto DM grid (global ideal)
		for(int idm=0; idm<parms->ndm; idm++){
			loc_t* loc=powfs[ipowfs].loc_dm?P(powfs[ipowfs].loc_dm, wfsind, idm):(powfs[ipowfs].loc_tel?P(powfs[ipowfs].loc_tel, wfsind):powfs[ipowfs].loc);
			const real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
			const real scale=1.-ht/hs;
			real dispx=ht*parms->wfs[iwfs].thetax+scale*misregx;
			real dispy=ht*parms->wfs[iwfs].thetay+scale*misregy;
			//wfs is registered to pupil. wfs.hc only effects the cone effect.
			
			if(scale<0) continue;
			prop_grid(P(simu->dmprojsq, idm), loc, P(opd),
				alpha, dispx, dispy, scale, 0, 0, 0);
		}
	}
}
/**
 * Compute tip/tilt for llt ray trace.
 * */
void wfsgrad_llt_tt(real*ttx, real*tty, sim_t* simu, int iwfs, int isim){
	const parms_t *parms=simu->parms;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int wfsind=P(parms->powfs[ipowfs].wfsind, iwfs);
	const int illt=P(parms->powfs[ipowfs].llt->i, wfsind);
	if(NE(simu->fsmreal, iwfs)){
		*ttx=P(P(simu->fsmreal, iwfs), 0);
		*tty=P(P(simu->fsmreal, iwfs), 1);
	}
	if(simu->telws){
		real tmp=P(simu->telws, isim)*parms->powfs[ipowfs].llt->ttrat;
		real angle=simu->winddir?P(simu->winddir, 0):0;
		*ttx+=tmp*cos(angle);
		*tty+=tmp*sin(angle);
	}
	if(simu->llt_ws&&P(simu->llt_ws, ipowfs)){
		*ttx+=P(P(simu->llt_ws, ipowfs), isim, illt);//put all to x direction.
	}
	if(simu->ltpm_real&&P(simu->ltpm_real, ipowfs)){
		*ttx+=P(P(simu->ltpm_real, ipowfs), 0, illt);
		*tty+=P(P(simu->ltpm_real, ipowfs), 1, illt);

		if(simu->save->ltpm_real){
			P(P(simu->save->ltpm_real, ipowfs), 0,isim)=P(P(simu->ltpm_real, ipowfs), 0, illt);
			P(P(simu->save->ltpm_real, ipowfs), 1,isim)=P(P(simu->ltpm_real, ipowfs), 1, illt);
		}
	}
}
/**
   computes close loop and pseudo open loop gradidents for both gometric and
   physical optics WFS. Calls wfsints() to accumulate WFS subapertures images in
   physical optics mode.  */

void* wfsgrad_iwfs(thread_t* info){
	sim_t* simu=(sim_t*)info->data;
	const int isim=simu->wfsisim;
	const int iwfs=info->start;
	const parms_t* parms=simu->parms;
	const int ipowfs=parms->wfs[iwfs].powfs;
	if(isim<parms->powfs[ipowfs].step) return NULL;
	assert(iwfs<parms->nwfs);
	/*
	  simu->gradcl is CL grad output (also for warm-restart of maxapriori
	  simu->gradacc is internal, to accumulate geometric grads.
	  do not accumulate opd. accumate ints for phy, g for GS
	*/
	/*input */

	const mapcell* atm=simu->atm;
	const recon_t* recon=simu->recon;
	const powfs_t* powfs=simu->powfs;
	/*output */
	const int CL=parms->sim.closeloop;
	const int nps=parms->atm.nps;
	const real atmscale=(simu->atmscale&&!parms->atm.dtrat)?P(simu->atmscale, isim):1;
	const real dt=parms->sim.dt;
	/*The following are truly constants for this powfs */
	const int imoao=parms->powfs[ipowfs].moao;
	const int wfsind=P(parms->powfs[ipowfs].wfsind, iwfs);
	const int dtrat=parms->powfs[ipowfs].dtrat;
	const int save_gradgeom=P(parms->save.gradgeom, iwfs);
	const int save_opd=P(parms->save.wfsopd, iwfs);
	const int save_ints=P(parms->save.ints, iwfs);
	const int noisy=parms->powfs[ipowfs].noisy;
	/*The following depends on isim */
	const int do_phy=simu->wfsflags[ipowfs].do_phy;
	const int do_pistat=simu->wfsflags[ipowfs].do_pistat;
	const int do_geom=(!do_phy||save_gradgeom||do_pistat)&&parms->powfs[ipowfs].type==WFS_SH;
	const dmat* amp=PR(powfs[ipowfs].amp, wfsind);
	dmat *gradcalc=NULL; //calculation output (geometric mode)
	dmat** gradacc=&P(simu->gradacc, iwfs);//accumulation output (geometric mode)
	dmat** gradout=&P(simu->gradcl, iwfs); //final output
	dcell* ints=P(simu->ints, iwfs);
	dmat* opd=P(simu->wfsopd, iwfs);
	TIM0;
	if(isim%dtrat==0){
		dcellzero(ints);
		dzero(*gradacc);
	}
	const int nhs=parms->powfs[ipowfs].llt?parms->powfs[ipowfs].llt->nhs:1;//number of sublayer
	const real dhs=parms->powfs[ipowfs].llt?parms->powfs[ipowfs].llt->dhs/nhs:0;//spacing of sublayer
	for(int ihs=0; ihs<nhs; ihs++){
		dzero(opd);
		const real hs=(nhs>1?(ihs-(nhs-1)*0.5):0)*dhs+parms->wfs[iwfs].hs;
		const real hc=nhs>1?(parms->wfs[iwfs].hc*(1.-hs/parms->wfs[iwfs].hs)):0;//effective hc
		/* Now begin ray tracing. */
		if(atm&&((!parms->sim.idealwfs&&!parms->powfs[ipowfs].lo)
				||(!parms->sim.wfsalias&&parms->powfs[ipowfs].lo))){
			for(int ips=0; ips<nps; ips++){
				thread_t* wfs_prop=simu->wfs_prop_atm[iwfs+parms->nwfs*ips];
				propdata_t* wfs_propdata=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
				wfs_propdata->phiout=opd;
				if(parms->atm.dtrat>0){
					real wt;
					int iframe=atm_interp(&wt, ips, isim, parms->atm.dtrat, NX(atm), parms->atm.interp);
					/*int iframe=wrap_seq(isim/parms->atm.dtrat+ips, NX(atm));
					real wt2=0;
					if(nps>1&&parms->atm.interp){
						wt2=(real)(isim%parms->atm.dtrat)/parms->atm.dtrat;
						if(parms->atm.interp==2){
							wt2=pow(sin(wt2*M_PI/2), 2);//smoother interp with sin^2 function
						}
						wfs_propdata->alpha=ips==0?(1-wt2):wt2;
					}*/
					wfs_propdata->alpha=atmscale*wt;
					wfs_propdata->mapin=P(atm, iframe);
					//if(iwfs==0) dbg("wfs: isim=%d, atm frame=%d, wt1=%g\n", isim, iframe, wfs_propdata->alpha);
				}else{
					wfs_propdata->displacex1=-P(atm, ips)->vx*dt*isim;
					wfs_propdata->displacey1=-P(atm, ips)->vy*dt*isim;
					wfs_propdata->alpha=atmscale;
				}
				if(nhs>1){
					const real ht=P(parms->atm.ht, ips);
					wfs_propdata->scale=1.-(ht-hc)/(hs-hc);
					//if(iwfs==0 && (ips+1)==nps) info("wfs=%d ips=%d ihs=%d, nhs=%d, scale=%g\n", iwfs, ips, ihs, nhs, wfs_propdata->scale);
				}
				/* have to wait to finish before another phase screen. */
				CALL_THREAD(wfs_prop, 1);
			}
		}
		/*
		Propagate controllable component of atm (within range of DM) to wfs.
		wfsalias: atm - controllable.
		idealwfs: just controllable.
		*/
		/* timing: most expensive 0.10 per LGS for*/
		if(!parms->powfs[ipowfs].lo&&(parms->sim.wfsalias||parms->sim.idealwfs)){
			real alpha=parms->sim.idealwfs?1:-1;
			wfs_ideal_atm(simu, opd, iwfs, alpha);
		}


		if(simu->telws){/*Wind shake */
			real tmp=P(simu->telws, isim);
			real angle=simu->winddir?P(simu->winddir, 0):0;
			real ptt[3]={0, tmp*cos(angle), tmp*sin(angle)};
			loc_add_ptt(opd, ptt, powfs[ipowfs].loc);
		}

		real focus=wfsfocusadj(simu, iwfs);
		if(fabs(focus)>1e-20){
			loc_add_focus(opd, powfs[ipowfs].loc, focus);
		}

		/* Add surface error*/
		if(powfs[ipowfs].opdadd&&P(powfs[ipowfs].opdadd, wfsind)){
			dadd(&opd, 1, P(powfs[ipowfs].opdadd, wfsind), 1);
		}

		if(save_opd){
			zfarr_push(simu->save->wfsopdol[iwfs], isim, opd);
		}
		TIM(1);
		if(CL){
			if(PARALLEL==2){
				wait_dmreal(simu, simu->wfsisim);
			}
			for(int idm=0; idm<parms->ndm; idm++){
				thread_t* wfs_prop=simu->wfs_prop_dm[iwfs+parms->nwfs*idm];
				propdata_t* wfs_propdata=&simu->wfs_propdata_dm[iwfs+parms->nwfs*idm];
				wfs_propdata->phiout=opd;
				if(nhs>1){
					const real ht=parms->dm[idm].ht+parms->dm[idm].vmisreg;
					wfs_propdata->scale=1.-(ht-hc)/(hs-hc);
					//if(iwfs==0) info("wfs=%d idm=%d ihs=%d, nhs=%d, scale=%g\n", iwfs, idm, ihs, nhs, wfs_propdata->scale);
				}
				CALL_THREAD(wfs_prop, 1);
			}/*idm */
			real ptt[3]={0,0,0};
			if(simu->ttmreal){
				ptt[1]-=P(simu->ttmreal, 0);
				ptt[2]-=P(simu->ttmreal, 1);
			}
			//For dithering with downlink instead of uplink FSM
			if(simu->fsmreal&&NE(simu->fsmreal, iwfs)&&!powfs[ipowfs].llt){
				ptt[1]-=P(P(simu->fsmreal, iwfs), 0);
				ptt[2]-=P(P(simu->fsmreal, iwfs), 1);
			}
			if(ptt[1]||ptt[2]){
				loc_add_ptt(opd, ptt, powfs[ipowfs].loc);
			}
			if(PARALLEL==2){
				post_dmreal(simu);
			}
		}
		if(parms->powfs[ipowfs].skip&&parms->tomo.ahst_idealngs==1){
			//apply ideal NGS modes to NGS WFS
			ngsmod_opd(opd, powfs[ipowfs].loc, recon->ngsmod,
				parms->wfs[iwfs].thetax, parms->wfs[iwfs].thetay,
				PCOL(simu->cleNGSm, isim), -1);
		}
		if(imoao>-1){
			dmat** dmwfs=P(simu->dm_wfs);
			if(dmwfs[iwfs]){
				/* No need to do mis registration here since the MOAO DM is attached
				to close to the WFS.*/
				prop_nongrid_pts(P(recon->moao[imoao].aloc, 0), P(dmwfs[iwfs]),
					powfs[ipowfs].pts, P(opd), -1, 0, 0, 1, 0, 0);
			}
		}

		if(parms->powfs[ipowfs].fieldstop>0&&parms->powfs[ipowfs].type==WFS_SH){
			int jwfs=PN(powfs[ipowfs].amp)>1?wfsind:0;
			locfft_fieldstop(powfs[ipowfs].fieldstop[jwfs], opd, parms->wfs[iwfs].wvlwts);
		}

		if(save_opd){
			zfarr_push(simu->save->wfsopd[iwfs], isim, opd);
		}
		if(parms->plot.run&&isim%parms->plot.run==0){
			drawopdamp("Opdwfs", powfs[ipowfs].loc, opd, amp, parms->plot.opdmax,
				"WFS OPD", "x (m)", "y (m)", "WFS %2d", iwfs);
		}
		if(do_geom){
			/* Now Geometric Optics gradient calculations. if dtrat==1, we compute
			gradients directly to gradacc, which is the same as gradcalc. If
			dtrat>1, we compute gradients to gradcalc, and accumulate to
			gradacc. gradcalc is used to shift pistat. We DONOT include gradoff
			adjustment to gradref, but only do it on gradcl. This will make the
			pistat always peak in center no matter what NCPA is present.
			*/
			if(!do_pistat||parms->powfs[ipowfs].pistatstc||dtrat==1){
				//we do not need separate gradcalc.
				gradcalc=dref(*gradacc);
			}//else: calculate first to gradcalc then add to gradacc
			if(parms->powfs[ipowfs].gtype_sim==GTYPE_Z){ /*compute ztilt. */
				pts_ztilt(&gradcalc, powfs[ipowfs].pts,
					PR(powfs[ipowfs].saimcc, wfsind, 0),
					P(amp), P(opd));
			} else{/*G tilt */
				dspmm(&gradcalc, PR(powfs[ipowfs].GS0, wfsind, 0), opd, "nn", 1);
			}
			if(P(gradcalc)!=P(*gradacc)){
				dadd(gradacc, 1, gradcalc, 1);
			}
		}

		ccell* psfout=NULL;
		zfarr* psfoutzfarr=NULL;
		zfarr* ztiltoutzfarr=NULL;
		if(parms->powfs[ipowfs].psfout){
			psfout=P(simu->wfspsfout, iwfs);
			psfoutzfarr=simu->save->wfspsfout[iwfs];
			ztiltoutzfarr=simu->save->ztiltout[iwfs];
		}
		TIM(2);
		/* Now begin Physical Optics Intensity calculations */
		if(do_phy||psfout||do_pistat||abs(parms->powfs[ipowfs].dither)==1){
			if(nhs>1){
				error("Please implement\n");
			}
			dmat* lltopd=NULL;
			if(powfs[ipowfs].llt){//If there is LLT, apply FSM onto LLT
				if(powfs[ipowfs].llt->ncpa){
					lltopd=ddup(PR(powfs[ipowfs].llt->ncpa, wfsind, 0));
				} else{
					lltopd=dnew(powfs[ipowfs].llt->pts->nxsa, powfs[ipowfs].llt->pts->nysa);
				}
				const long illt=P(parms->powfs[ipowfs].llt->i, wfsind);
				if(atm){/*LLT OPD */
				real wt=1;
					for(int ips=0; ips<nps; ips++){
						const real hl=P(atm, ips)->h;
						const real scale=1.-hl/hs;
						if(scale<0) continue;
						const real ox=P(parms->powfs[ipowfs].llt->ox, illt);
						const real oy=P(parms->powfs[ipowfs].llt->oy, illt);
						const real thetax=parms->wfs[iwfs].thetax-ox/hs;
						const real thetay=parms->wfs[iwfs].thetay-oy/hs;
						real vx=0;
						real vy=0;
						map_t *atmi;
						if(parms->atm.dtrat>0){
							int iframe=atm_interp(&wt, ips, isim, parms->atm.dtrat, NX(atm), parms->atm.interp);
							/*int iframe=wrap_seq(isim/parms->atm.dtrat+ips, NX(atm));
							real wt2=0;
							if(nps>1&&parms->atm.interp){
								wt2=(real)(isim%parms->atm.dtrat)/parms->atm.dtrat;
								if(parms->atm.interp==2){
									wt2=pow(sin(wt2*M_PI/2), 2);//smoother interp with sin^2 function
								}
								//atmscale=ips==0?(1-wt2):wt2;
							}*/
							atmi=P(atm, iframe);
							//if(iwfs==0) dbg("lltopd: isim=%d, atm frame=%d, wt1=%g\n", isim, iframe, atmscale);
						}else{
							vx=P(atm, ips)->vx;
							vy=P(atm, ips)->vy;
							atmi=P(atm, ips);
						}
						const real displacex=-vx*isim*dt+thetax*hl+ox;
						const real displacey=-vy*isim*dt+thetay*hl+oy;
						
						prop_grid_pts(atmi, powfs[ipowfs].llt->pts,
							P(lltopd), atmscale*wt, displacex, displacey,
							scale, 1., 0, 0);
					}
				}
				if(do_pistat||parms->powfs[ipowfs].idealfsm){
					/* remove tip/tilt completely */
					dmat* lltg=dnew(2, 1);
					pts_ztilt(&lltg, powfs[ipowfs].llt->pts,
						powfs[ipowfs].llt->imcc,
						P(powfs[ipowfs].llt->amp),
						P(lltopd));
					P(P(simu->fsmreal, iwfs), 0)=-P(lltg, 0);
					P(P(simu->fsmreal, iwfs), 1)=-P(lltg, 1);
					dfree(lltg);
				}
				real ttx=0, tty=0;//uplink jitter and correction
				wfsgrad_llt_tt(&ttx, &tty, simu, iwfs, isim);
				if(ttx!=0||tty!=0){ /* add tip/tilt to llt opd */
					real ptt[3]={0, ttx, tty};
					loc_add_ptt(lltopd, ptt, powfs[ipowfs].llt->loc);
				}
				if(save_opd){
					zfarr_push(simu->save->wfslltopd[iwfs], isim, lltopd);
				}
			}
			if(parms->powfs[ipowfs].type==WFS_SH){//SHWFS
				wfsints_t* intsdata=simu->wfs_intsdata+iwfs;
				intsdata->ints=ints;
				intsdata->psfout=psfout;
				intsdata->pistatout=P(simu->pistatout, iwfs);
				if(parms->powfs[ipowfs].pistatout){
					intsdata->gradref=gradcalc;
				}
				intsdata->opd=opd;
				intsdata->lltopd=lltopd;
				intsdata->isim=isim;
				CALL_THREAD(simu->wfs_ints[iwfs], 1);
				dfree(lltopd);
				intsdata->lltopd=0;
				intsdata->opd=0;
				if(psfout){
					zfarr_push(psfoutzfarr, isim, psfout);
					zfarr_push(ztiltoutzfarr, isim, *gradacc);
				}
			} else{//Pywfs
				pywfs_ints(&P(ints, 0), powfs[ipowfs].pywfs, opd, parms->wfs[iwfs].sigsim);
			}
		}
		TIM(3);
		if(simu->wfsflags[ipowfs].gradout && (ihs+1)==nhs){
			if(do_phy){
				/* In Physical optics mode, do integration and compute
				gradients. The matched filter are in x/y coordinate even if
				radpix=1. */
				if(save_ints){
					zfarr_push(simu->save->intsnf[iwfs], isim, ints);
				}
				if(noisy){/*add noise */
					if(P(parms->save.gradnf, iwfs)){//save noise free gradients
						if(parms->powfs[ipowfs].type==WFS_SH){
							shwfs_grad(gradout, P(ints), parms, powfs, iwfs, parms->powfs[ipowfs].phytype_sim);
						} else{
							pywfs_grad(gradout, powfs[ipowfs].pywfs, P(ints, 0));
						}
						zfarr_push(simu->save->gradnf[iwfs], isim, *gradout);
					}
					const real rne=parms->powfs[ipowfs].rne;
					const real bkgrnd=parms->powfs[ipowfs].bkgrnd*dtrat;
					const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
					dmat** bkgrnd2=NULL;
					dmat** bkgrnd2c=NULL;
					if(powfs[ipowfs].bkgrnd){
						bkgrnd2=PCOLR(powfs[ipowfs].bkgrnd, wfsind);
					}
					if(powfs[ipowfs].bkgrndc){
						bkgrnd2c=PCOLR(powfs[ipowfs].bkgrndc, wfsind);
					}
					for(int isa=0; isa<NX(ints); isa++){
						dmat* bkgrnd2i=(bkgrnd2)?bkgrnd2[isa]:NULL;
						dmat* bkgrnd2ic=(bkgrnd2c)?bkgrnd2c[isa]:NULL;
						addnoise(P(ints, isa), &simu->wfs_rand[iwfs],
							bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, parms->powfs[ipowfs].qe, rne, 1.);
					}
					if(save_ints){
						zfarr_push(simu->save->intsny[iwfs], isim, ints);
					}
				}
				if(parms->powfs[ipowfs].i0save==2){
					dcelladd(&P(simu->ints, iwfs), 1, ints, 1);
				}
				if(abs(parms->powfs[ipowfs].dither)==1 && isim>=parms->powfs[ipowfs].dither_ogskip
					&&parms->powfs[ipowfs].type==WFS_SH
					&&(parms->powfs[ipowfs].dither==-1||parms->powfs[ipowfs].phytype_sim2==PTYPE_MF)){
					/*Collect statistics with dithering*/
					dither_t* pd=simu->dither[iwfs];
					dcelladd(&pd->imb, 1, ints, 1.);
					if(parms->powfs[ipowfs].dither==1){
						real cs, ss;
						dither_position(&cs, &ss, parms->powfs[ipowfs].alfsm, parms->powfs[ipowfs].dtrat,
							parms->powfs[ipowfs].dither_npoint, isim, pd->deltam);
						//accumulate for matched filter

						dcelladd(&pd->imx, 1, ints, cs);
						dcelladd(&pd->imy, 1, ints, ss);
					}
				}

				if(parms->powfs[ipowfs].type==WFS_SH){
					shwfs_grad(gradout, P(ints), parms, powfs, iwfs, parms->powfs[ipowfs].phytype_sim);
				} else{
					pywfs_grad(gradout, powfs[ipowfs].pywfs, P(ints, 0));
				}
			} else{
				/* geomtric optics accumulation mode. scale and copy results to output. */
				dcp(gradout, *gradacc);
				if(dtrat!=1||nhs!=1){
					dscale(*gradout, 1./(dtrat*nhs));/*average */
				}
				if(P(parms->save.gradnf, iwfs)){
					zfarr_push(simu->save->gradnf[iwfs], isim, *gradout);
				}

				if(noisy&&!parms->powfs[ipowfs].usephy){
					const dmat* neasim=PR(powfs[ipowfs].neasim, wfsind, 0);//neasim is the LL' decomposition
					addnoise_grad(*gradout, neasim, &simu->wfs_rand[iwfs]);
				}
			}
			if(save_gradgeom&&do_phy){
				dmat* gradtmp=NULL;
				dadd(&gradtmp, 1, *gradacc, 1./dtrat);
				zfarr_push(simu->save->gradgeom[iwfs], isim, gradtmp);/*noise free. */
				dfree(gradtmp);
			}
		}//dtrat_out
	}//for ihs
	dfree(gradcalc);
	TIM(4);
	TIM1;
	return NULL;
}


/*Compute global tip/tilt error for each WFS*/
static void wfsgrad_fsm(sim_t* simu, int iwfs){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int isim=simu->wfsisim;
	/*Uplink FSM*/
	const int ind=parms->recon.glao?(ipowfs+ipowfs*parms->npowfs):(iwfs+iwfs*parms->nwfs);
	const dmat* PTT=recon->PTT?(P(recon->PTT, ind)):0;
	if(!PTT){
		error("powfs %d has FSM, but PTT is empty\n", ipowfs);
	}
	/* Compute FSM error. */
	simu->fsmerr=simu->fsmerr_store;
	dmm(&P(simu->fsmerr, iwfs), 0, PTT, P(simu->gradcl, iwfs), "nn", 1);
	//2021-09-16: drift signal is treated as bias. do not zero fsmerr_drift
	dadd(&P(simu->fsmerr, iwfs), 1, P(simu->fsmerr_drift, iwfs), 1);
	//Save data
	P(P(simu->save->fsmerrs, iwfs), 0, isim)=P(P(simu->fsmerr, iwfs), 0);
	P(P(simu->save->fsmerrs, iwfs), 1, isim)=P(P(simu->fsmerr, iwfs), 1);
}
static void wfsgrad_tt_drift(dmat* grad, sim_t* simu, real gain, int iwfs, int remove){
	//gain can be set to 1 if the rate is slower than the main tip/tilt and focus control rate.
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int iwfsr=parms->recon.glao?ipowfs:iwfs;

	if(parms->powfs[ipowfs].trs){
		/* Update T/T drift signal to prevent matched filter from drifting*/
		dmat* tt=dnew(2, 1);
		const dmat* PTT=P(recon->PTT, iwfsr, iwfsr);
		dmm(&tt, 0, PTT, grad, "nn", 1);
		if(remove){
			dmat* TT=P(recon->TT, iwfsr, iwfsr);
			dmm(&grad, 1, TT, tt, "nn", -1);
		}
		//need to use integrator here and then use as offset
		if(gain){
			P(P(simu->fsmerr_drift, iwfs), 0)+=P(tt, 0)*gain;
			P(P(simu->fsmerr_drift, iwfs), 1)+=P(tt, 1)*gain;
			if(iwfs==P(parms->powfs[ipowfs].wfs, 0)){
				dbg("Step %5d: wfs %d uplink drift control error is (%.3f, %.3f) mas output is (%.3f, %.3f) mas.\n",
				simu->wfsisim, iwfs, P(tt, 0)*RAD2MAS, P(tt, 1)*RAD2MAS,
				P(P(simu->fsmerr_drift, iwfs), 0)*RAD2MAS, P(P(simu->fsmerr_drift, iwfs), 1)*RAD2MAS);
			}
		}
		dfree(tt);
	}
}
static void wfsgrad_focus_drift(dmat* grad, sim_t* simu, real gain, int iwfs, int remove){
	//gain can be set to 1 as zoomgain is applied at zoomint.
	const parms_t* parms=simu->parms;
	const int ipowfs=parms->wfs[iwfs].powfs;
	//Output focus error in ib to trombone error signal.
	if(parms->powfs[ipowfs].llt){
		//here we don't use RFlgsg which is noise weighted
		real focus=loc_remove_focus_grad(grad, simu->powfs[ipowfs].saloc, remove?1:0);
		if(gain){
			P(simu->zoomdrift, iwfs)+=gain*focus;//accumulate to the averager
			P(simu->zoomdrift_count, iwfs)++;
		}
	}
}
/**
   Accumulate dithering parameters
   - Every step: accumulate signal for phase detection.
   - At PLL output: determine input/output amplitude of dithering signal.
   - At Gain output:determine matched filter i0, gx, gy, or CoG gain.
   - Subtract t/t from gradients for non-comon-path (TT) dithering.

*/
static void wfsgrad_dither(sim_t* simu, int iwfs){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const powfs_t* powfs=simu->powfs;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const int iwfsr=parms->recon.glao?ipowfs:iwfs;
	const int isim=simu->wfsisim;
	const int pllrat=parms->powfs[ipowfs].dither_pllrat;
	if(!parms->powfs[ipowfs].dither||isim<parms->powfs[ipowfs].dither_pllskip){
		return;
	}
	real cs, ss; //Current phase of tip/tilt dithering signal
	dither_t* pd=simu->dither[iwfs];
	if(parms->powfs[ipowfs].dither==1){ //T/T dithering.
		//Current dithering signal phase
		dither_position(&cs, &ss, parms->powfs[ipowfs].alfsm, parms->powfs[ipowfs].dtrat,
			parms->powfs[ipowfs].dither_npoint, isim, pd->deltam);

		/* Use delay locked loop to determine the phase of actual
		dithering signal (position of LGS spot averaged over a WFS
		integration period) using measured signal (WFS global
		tip/tilt). In actual system, the LGS uplink propagation
		needs to be accounted.
		*/
		real err;
		err=(-ss*(P(P(simu->fsmerr, iwfs), 0))
			+cs*(P(P(simu->fsmerr, iwfs), 1)))/(parms->powfs[ipowfs].dither_amp);
		pd->delta+=parms->powfs[ipowfs].dither_gpll*(err/pllrat);

		//For SHWFS CoG gaim update.
		if(parms->powfs[ipowfs].type==WFS_SH&&parms->powfs[ipowfs].phytype_sim2!=PTYPE_MF&&isim>=parms->powfs[ipowfs].dither_ogskip){
			const int nsa=powfs[ipowfs].saloc->nloc;
			if(!pd->ggm){
				pd->ggm=dnew(nsa*2, 1);
			}
			for(int isa=0; isa<nsa; isa++){
				P(pd->ggm, isa)+=cs*P(P(simu->gradcl, iwfs), isa);
				P(pd->ggm, isa+nsa)+=ss*P(P(simu->gradcl, iwfs), isa+nsa);
			}
		}
	} else if(parms->powfs[ipowfs].dither>1){ //DM dithering.

		const int idm=parms->idmground;
		//Input dither signal
		dmat *mr_in=drefcols(P(pd->mr, 0), isim, 1);
		dmm(&mr_in, 0, P(recon->dither_ra, iwfs, idm), P(simu->dmreal, idm), "nn", 1);
		dfree(mr_in);
		/*P(P(pd->mr, 0), isim)=P(tmp, 0);
		if(PN(tmp)>1){//2 mode dithering
			P(P(pd->mr,2),isim)=P(tmp,1);
		}*/

		dmat *mr_out=drefcols(P(pd->mr, 1), isim, 1);
		//Measured dither signal from gradients
		dmm(&mr_out, 0, P(recon->dither_rg, iwfs, iwfs), P(simu->gradcl, iwfs), "nn", 1);

		/*P(P(pd->mr, 1), isim)=P(mr_out, 0);
		if(PN(mr_out)>1){
			P(P(pd->mr, 3), isim)=P(mr_out, 1);
		}*/
		dfree(mr_out);
	}
	if(simu->wfsflags[ipowfs].pllout&&parms->powfs[ipowfs].dither>0){
		//Synchronous detection of dither signal amplitude in input (DM) and output (gradients).
		//The ratio between the two is used for (optical) gain adjustment.
		const int npoint=parms->powfs[ipowfs].dither_npoint;
		const int ncol=(pllrat-1)*parms->powfs[ipowfs].dtrat+1;
		if(parms->powfs[ipowfs].dither==1){//TT
			//dbg("deltam=%g is updated to %g+%g=%g\n", pd->deltam, pd->delta, pd->deltao, pd->delta+pd->deltao);
			pd->deltam=pd->delta+(pd->deltao*parms->powfs[ipowfs].dither_gdrift);//output PLL
			dmat* tmp=0;
			const int detrend=parms->powfs[ipowfs].llt?0:1;
			tmp=drefcols(P(simu->save->fsmcmds, iwfs), simu->wfsisim-ncol+1, ncol);
			pd->a2m=calc_dither_amp(NULL, tmp, parms->powfs[ipowfs].dtrat, npoint, detrend, 1);
			dfree(tmp);
			tmp=drefcols(P(simu->save->fsmerrs, iwfs), simu->wfsisim-ncol+1, ncol);
			pd->a2me=calc_dither_amp(NULL, tmp, parms->powfs[ipowfs].dtrat, npoint, detrend, 1);
			dfree(tmp);
		} else if(parms->powfs[ipowfs].dither>1){//DM
			dmat* tmp=0;
			int detrend=1;//1: default.
			tmp=drefcols(P(pd->mr, 0), simu->wfsisim-ncol+1, ncol);//DM
			pd->a2m=calc_dither_amp(&pd->a2mv, tmp, parms->powfs[ipowfs].dtrat, npoint, detrend, 0);
			dfree(tmp);
			tmp=drefcols(P(pd->mr, 1), simu->wfsisim-ncol+1, ncol);//Grad
			pd->a2me=calc_dither_amp(&pd->a2mev, tmp, parms->powfs[ipowfs].dtrat, npoint, detrend, 0);
			dfree(tmp);
			/*if(PN(pd->a2mv)>1){//multi-mode
				warning_once("Use the last dithering mode\n");
				pd->a2m=P(pd->a2mv, PN(pd->a2mv)-1);
				pd->a2me=P(pd->a2mev, PN(pd->a2mev)-1);
			}*/
			/*if(PN(pd->mr)>3){
				tmp=drefcols(P(pd->mr, 2), simu->wfsisim-ncol+1, ncol);//DM
				pd->a2m2=calc_dither_amp(tmp, parms->powfs[ipowfs].dtrat, npoint, detrend);
				dfree(tmp);
				tmp=drefcols(P(pd->mr, 3), simu->wfsisim-ncol+1, ncol);//Grad
				pd->a2me2=calc_dither_amp(tmp, parms->powfs[ipowfs].dtrat, npoint, detrend);
				dfree(tmp);
			}*/
		}

		/*
		//Print PLL phase. Moved to OG gain print out
		if(iwfs==P(parms->powfs[ipowfs].wfs, 0)){
			const real anglei=(2*M_PI/parms->powfs[ipowfs].dither_npoint);
			const real scale=parms->powfs[ipowfs].dither==1?(1./parms->powfs[ipowfs].dither_amp):1;
			if(pd->a2m2){
				info2("Step %5d wfs %d PLL: delay=%.2f frame, dither amplitude=(%.2fx %.2fx), estimate=(%.2fx, %.2fx)\n",
					isim, iwfs, pd->deltam/anglei, pd->a2m*scale, pd->a2m2*scale, pd->a2me*scale, pd->a2me2*scale);
			}else{
				info2("Step %5d wfs %d PLL: delay=%.2f frame, dither amplitude=%.2fx, estimate=%.2fx\n",
					isim, iwfs, pd->deltam/anglei, pd->a2m*scale, pd->a2me*scale);
			}
		}*/
		if(simu->resdither){
			int ic=simu->wfsflags[ipowfs].pllout-1;
			P(P(simu->resdither, iwfs), 0, ic)=pd->deltam;
			P(P(simu->resdither, iwfs), 1, ic)=pd->a2m;
			P(P(simu->resdither, iwfs), 2, ic)=pd->a2me;
			if(PN(pd->a2mv)>1){
				int nm=PN(pd->a2mv);
				memcpy(&P(P(simu->resdither, iwfs), 4, ic), P(pd->a2mv), nm*sizeof(real));
				memcpy(&P(P(simu->resdither, iwfs), 4+nm, ic), P(pd->a2mev), nm*sizeof(real));
				/*P(P(simu->resdither, iwfs), 4, ic)=pd->a2m2;
				P(P(simu->resdither, iwfs), 5, ic)=pd->a2me2;*/
			}
		}
	}

	if(simu->wfsflags[ipowfs].ogacc){//Gain update statistics
		/*if(parms->dbg.gradoff_reset==2 && simu->gradoffisim0<=0
			&& parms->powfs[ipowfs].phytype_sim2==PTYPE_MF){//trigger accumulation of gradoff
			simu->gradoffisim0=simu->wfsisim;
			simu->gradoffisim=simu->wfsisim;
		}*/
		if(abs(parms->powfs[ipowfs].dither)==1){//TT Dither or i0 collection
			real scale1=1./pllrat;
			real amp=pd->a2m;
			real scale2=amp?(scale1*2./(amp)):0;
			if(pd->imb){//computer i0, gx, gy for matched filter
				dcellscale(pd->imb, scale1);
				//Accumulate data for matched filter
				dcelladd(&pd->i0, 1, pd->imb, 1);//imb was already scaled
				dcellzero(pd->imb);
				if(parms->powfs[ipowfs].dither==1){
					dcelladd(&pd->gx, 1, pd->imx, scale2);
					dcelladd(&pd->gy, 1, pd->imy, scale2);

					dcellzero(pd->imx);
					dcellzero(pd->imy);
				}
			} else if(pd->ggm&&parms->powfs[ipowfs].dither==1){//cog
				dadd(&pd->gg0, 1, pd->ggm, scale2);
				dzero(pd->ggm);
			}
		}
	}

	if(parms->powfs[ipowfs].dither==1){
		/* subtract estimated tip/tilt dithering signal to avoid perturbing the loop or dithering pattern.*/
		real amp=pd->a2me;
		real tt[2]={-cs*amp, -ss*amp};
		if(parms->powfs[ipowfs].trs){
			//info("fsmerr: %g %g %g %g\n", P(P(simu->fsmerr,iwfs),0), P(P(simu->fsmerr,iwfs),1), -tt[0], -tt[1]);
			if(!amp){//no estimate yet, do not close up FSM loop.
				P(P(simu->fsmerr, iwfs), 0)=0;
				P(P(simu->fsmerr, iwfs), 1)=0;
			} else{
				P(P(simu->fsmerr, iwfs), 0)+=tt[0];
				P(P(simu->fsmerr, iwfs), 1)+=tt[1];
			}
		}
		//also remove from gradient measurements.
		dmulvec(P(P(simu->gradcl, iwfs)), P(recon->TT, iwfsr, iwfsr), tt, 1);
	}
}

/**
   Accomplish two tasks:
   1) Use LPF'ed LGS focus measurement to drive the trombone.
   2) HPF LGS focus on gradients to remove sodium range variation effact.

   We trust the focus measurement of the LGS WFS at high temporal frequency
   which NGS cannot provide due to low frame rate. After the HPF on lgs
   gradients, our system is NO LONGER affected by sodium layer variation.

   if sim.mffocus==1: The HPF is applied to each LGS WFS independently. This largely
   removes the effect of differential focus. powfs.dfrs is no longer needed. (preferred).

   if sim.mffocus==2: We apply a LPF on the average focus from six LGS WFS, and
   then remove this value from all LGS WFS. The differential focus is still
   present and powfs.dfrs need to be set to 1 to handle it in tomography. This
   is the original focus tracking method, and is no longer recommended.

   2022-04-12: the cut off frequency is now set to inf to ignore LGS focus completely.
*/
static void wfsgrad_lgsfocus(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;

	dcell* LGSfocus=simu->LGSfocus;//computed in wfsgrad_post from gradcl.

	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		const int do_phy=simu->wfsflags[ipowfs].do_phy;
		if(!simu->wfsflags[ipowfs].gradout||!parms->powfs[ipowfs].llt
			||simu->wfsisim<parms->powfs[ipowfs].step||!do_phy){
			continue;
		}

		/*New plate mode focus offset for LGS WFS. Not really needed*/
		if(parms->tomo.ahst_focus==2
			&&simu->Mint_lo&&P(simu->Mint_lo->mint, 1)&&simu->wfsflags[ipowfs].gradout){
			 /*When tomo.ahst_focus>0, the first plate scale mode contains focus for
			   lgs. But it turns out to be not necessary to remove it because the
			   HPF in the LGS path removed the influence of this focus mode. set
			   tomo.ahst_focus=2 to enable adjust gradients. (testing only)*/

			real scale=simu->recon->ngsmod->scale;
			int indps=simu->recon->ngsmod->indps;
			dmat* mint=P(P(simu->Mint_lo->mintc, 0), 0);//2018-12-11: changed first p[1] to p[0]
			real focus=P(mint, indps)*(scale-1);
			for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
				dadd(&P(simu->gradcl, iwfs), 1, P(recon->GFall, iwfs), focus);
			}
		}
		real lgsfocusm=0;//LGS averaged focus
		if(parms->sim.mffocus==2){
			for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
				int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
				lgsfocusm+=P(P(LGSfocus, iwfs), 0);
			}
			lgsfocusm/=parms->powfs[ipowfs].nwfs;
		}

		/*Here we set trombone position according to focus in the first
		  measurement. And adjust the focus content of this
		  measurement. This simulates the initial focus acquisition
		  step. No need if start with pre-built matched filter.*/
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
			int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
			P(P(simu->LGSfocusts, iwfs), simu->wfsisim)=P(P(LGSfocus, iwfs), 0);//save time history
			if(parms->powfs[ipowfs].zoomgain){
				//Trombone from gradients. always enable
				P(simu->zoomavg, iwfs)+=P(P(LGSfocus, iwfs), 0);//zoom averager
				P(simu->zoomavg_count, iwfs)++;
			}
			if(parms->sim.mffocus){//Focus HPF
				real focus=parms->sim.mffocus==2?lgsfocusm:P(P(LGSfocus, iwfs), 0);
				/*2022-08-26: focus removal in tomo RHS (unless LHS also) has
				worse performance because it removes focus from PSOL grads, not
				CL grads.*/
				dadd(&P(simu->gradcl, iwfs), 1, P(recon->GFall, iwfs), -P(simu->lgsfocuslpf, iwfs));
				P(P(simu->LGSfocusts, iwfs), simu->wfsisim)+=-P(simu->lgsfocuslpf, iwfs);
				//LPF is after using the value to put it off critical path of the RTC.
				real lpfocus=parms->sim.lpfocushi;
				P(simu->lgsfocuslpf, iwfs)=P(simu->lgsfocuslpf, iwfs)*(1-lpfocus)+focus*lpfocus;
			}
		}//for jwfs
	}//for ipowfs
}

/**
   Every operation here should be in the Simulator not the Controller
*/
void* wfsgrad_post(thread_t* info){
	sim_t* simu=(sim_t*)info->data;
	const parms_t* parms=simu->parms;
	//Postprocessing gradients
	const int isim=simu->wfsisim;
	for(int iwfs=info->start; iwfs<info->end; iwfs++){
		const int ipowfs=parms->wfs[iwfs].powfs;
		if(isim<parms->powfs[ipowfs].step) continue;
#if USE_CUDA
		if(parms->gpu.wfs){
			gpu_wfsgrad_sync(simu, iwfs);
		}
#endif
		const int do_phy=simu->wfsflags[ipowfs].do_phy;
		dmat* gradcl=P(simu->gradcl, iwfs);
		/* copy fsmreal to output  */
		if(NE(simu->fsmreal, iwfs)){
			P(P(simu->save->fsmcmds, iwfs), 0, isim)=P(P(simu->fsmreal, iwfs), 0);
			P(P(simu->save->fsmcmds, iwfs), 1, isim)=P(P(simu->fsmreal, iwfs), 1);
		}
		if(simu->wfsflags[ipowfs].gradout){
			if(parms->plot.run&&isim%parms->plot.run==0){
				/*drawgrad("Gcl", simu->powfs[ipowfs].saloc, gradcl,
					parms->plot.grad2opd, parms->powfs[ipowfs].trs, P(parms->plot.gmax),
					"WFS Closeloop Gradients", "x (m)", "y (m)", "Gcl %d", iwfs);*/
				if(do_phy){
					drawints("Ints", simu->powfs[ipowfs].saloc, P(simu->ints, iwfs), 0,
						"WFS Subaperture Images", "x", "y", "WFS %2d", iwfs);
				}
			}

			if(P(simu->gradscale, iwfs)){
				if(PN(simu->gradscale, iwfs)==1){
					dscale(gradcl, P(P(simu->gradscale, iwfs),0));
				}else{
					dcwm(gradcl, P(simu->gradscale, iwfs));
				}
			} else{
				dscale(gradcl, parms->powfs[ipowfs].gradscale);
			}
			if(P(simu->gradoff, iwfs)){
				dadd(&P(simu->gradcl, iwfs), 1, P(simu->gradoff, iwfs), -parms->dbg.gradoff_scale);
			}

			if(parms->dbg.gradoff){
				info_once("wfs %d: add dbg.gradoff to gradient vector\n", iwfs);
				int icol=(simu->wfsisim+1)%NY(parms->dbg.gradoff);
				dadd(&P(simu->gradcl, iwfs), 1, P(parms->dbg.gradoff, iwfs, icol), -1);
			}

			if(do_phy){
				if(NE(simu->fsmerr_store, iwfs)){
					wfsgrad_fsm(simu, iwfs);
				}
				if(parms->powfs[ipowfs].dither){
					wfsgrad_dither(simu, iwfs);
				}
				if(!parms->powfs[ipowfs].trs&&parms->powfs[ipowfs].skip!=2&&simu->fsmerr){
					dzero(P(simu->fsmerr, iwfs));//do not close fsm loop when t/t is used for AO.
				}
			}
			if(parms->powfs[ipowfs].llt){
				dmm(&P(simu->LGSfocus, iwfs), 0, P(simu->recon->RFlgsg, iwfs, iwfs), P(simu->gradcl, iwfs), "nn", 1);
			}
			if(P(parms->save.grad, iwfs)){
				zfarr_push(simu->save->gradcl[iwfs], isim, gradcl);
			}
			/*if(parms->powfs[ipowfs].type==WFS_PY && isim%100==99){
				pywfs_gain_calibrate(simu->powfs[ipowfs].pywfs, gradcl, parms->atm.r0);
			}*/
		}
	}//for iwfs
	return NULL;
}
static void gradoff_acc(sim_t* simu, int ipowfs){
	(void)ipowfs;
	if(simu->parms->dbg.gradoff_reset==2&&simu->gradoffisim0>0){//accumulate gradoff before updating it.
		int nsim=(simu->wfsisim-simu->gradoffisim);
		if(nsim){
			dcelladd(&simu->gradoffacc, 1, simu->gradoff, nsim);
			info("step %d: gradoff is accumulated with factor %d\n", simu->wfsisim, nsim);
		}
		simu->gradoffisim=simu->wfsisim;
	}
}
/**
	Controls gradient drift of individual subapertures due to dithering for matched filter
*/
static void wfsgrad_sa_drift(sim_t* simu, int ipowfs){
	const parms_t* parms=simu->parms;
	const powfs_t* powfs=simu->powfs;
	if(parms->powfs[ipowfs].dither_gdrift==0||parms->powfs[ipowfs].phytype_sim!=PTYPE_MF) return;
	dmat* goff=0;
	intstat_t* intstat=simu->powfs[ipowfs].intstat;
	const int isim=simu->wfsisim;
	gradoff_acc(simu, ipowfs);

	for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
		const int nsa=NX(intstat->i0);
		//outer loop to prevent i0 from drifting
		//Compute CoG of i0 + goff and drive it toward gradoff with low gain (0.1)
		if(1){
			dzero(goff);
			shwfs_grad(&goff, PCOL(intstat->i0, jwfs), parms, powfs, iwfs, PTYPE_COG);//force cog
			dadd(&goff, 1, P(simu->gradoff, iwfs), 1);
			if(simu->gradoffdrift){//cog boot strapped
				dadd(&goff, 1, P(simu->gradoffdrift, iwfs), -1);
			} else if(simu->powfs[ipowfs].gradoff){//cmf boot strapped, gradoff is cog of i0
				dadd(&goff, 1, P(simu->powfs[ipowfs].gradoff, jwfs), -1);
			}
			dadd(&P(simu->gradoff, iwfs), 1, goff, -parms->powfs[ipowfs].dither_gdrift);
			if(jwfs==0){
				dbg("Step %5d: powfs %d gradient drift control\n", simu->wfsisim, ipowfs);
			}
		}
		if(1){
			//outer loop to prevent gx/gy direction from drifting.
			//It computes CoG of shifted images (i0+gx/gy) and make sure the angle stays the same.
			//may not be necessary.
			dmat* i0sx=0, * i0sy=0;
			real theta=0;
			const real gyoff=M_PI*0.5;
			const real gshift=parms->powfs[ipowfs].pixtheta*0.1;
			const real cogthres=parms->powfs[ipowfs].cogthres;
			const real cogoff=parms->powfs[ipowfs].cogoff;

			for(int isa=0; isa<nsa; isa++){
				real g0[2], gx[2], gy[2];
				dcp(&i0sx, PR(intstat->i0, isa, jwfs));
				dcog(g0, i0sx, 0., 0., cogthres, cogoff, 0, NULL);
				dcp(&i0sy, PR(intstat->i0, isa, jwfs));
				dadd(&i0sx, 1, PR(intstat->gx, isa, jwfs), gshift);
				dadd(&i0sy, 1, PR(intstat->gy, isa, jwfs), gshift);
				dcog(gx, i0sx, 0., 0., cogthres, cogoff, 0, NULL);
				dcog(gy, i0sy, 0., 0., cogthres, cogoff, 0, NULL);

				//Works in both x/y and r/a coordinate.
				theta+=(atan2(gx[1]-g0[1], gx[0]-g0[0])+atan2(gy[1]-g0[1], gy[0]-g0[0])-gyoff);
			}
			theta*=0.5/nsa;
			simu->dither[iwfs]->deltao=-theta;
			dfree(i0sx);
			dfree(i0sy);
			if(jwfs==0){
				dbg("Step %5d: wfs %d angle drift control deltao is %g.\n", isim, iwfs, simu->dither[iwfs]->deltao);
			}
		}
	}
	dfree(goff);
	if(parms->save.gradoff||parms->save.dither){
		writebin(simu->gradoff, "extra/gradoff_%d_drift", isim);
	}
}
/**
   Dither update: zoom corrector, matched filter, gain ajustment, TWFS.
*/
static void wfsgrad_dither_post(sim_t* simu){
	powfs_t* powfs=simu->powfs;//not const to set instat. //todo: move intstat to simu.
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int isim=simu->wfsisim;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(!parms->powfs[ipowfs].dither || isim<parms->powfs[ipowfs].step || (isim+1)%parms->powfs[ipowfs].dtrat!=0) continue;
		const int nwfs=parms->powfs[ipowfs].nwfs;

		if(simu->wfsflags[ipowfs].ogout){//This is matched filter or cog update
			const int nsa=powfs[ipowfs].saloc->nloc;
			const real scale1=(real)parms->powfs[ipowfs].dither_pllrat/(real)parms->powfs[ipowfs].dither_ograt;
			int ptype2=parms->powfs[ipowfs].phytype_sim2;
			if(parms->powfs[ipowfs].dither==-1||(ptype2==PTYPE_MF)){
				if(parms->powfs[ipowfs].dither==-1){
					if(parms->powfs[ipowfs].llt){
						info2("Step %5d: Update sodium fit for powfs %d\n", isim, ipowfs);
					}else{
						info2("Step %5d: Update i0 processing for powfs %d\n", isim, ipowfs);
					}
				} else{
					info2("Step %5d: Update matched filter for powfs %d\n", isim, ipowfs);
				}
				//For matched filter
				if(!powfs[ipowfs].intstat){
					powfs[ipowfs].intstat=mycalloc(1, intstat_t);
				}
				intstat_t* intstat=powfs[ipowfs].intstat;
				parms->powfs[ipowfs].radgx=0;//ensure derivate is interpreted as along x/y.
				int pixpsax=powfs[ipowfs].pixpsax;
				int pixpsay=powfs[ipowfs].pixpsay;
				if(!intstat->i0||NY(intstat->i0)!=nwfs){
					dcellfree(intstat->i0);
					intstat->i0=dcellnew_same(nsa, nwfs, pixpsax, pixpsay);
				}
				if(!intstat->gx||NY(intstat->gx)!=nwfs){
					dcellfree(intstat->gx);
					dcellfree(intstat->gy);
					intstat->gx=dcellnew_same(nsa, nwfs, pixpsax, pixpsay);
					intstat->gy=dcellnew_same(nsa, nwfs, pixpsax, pixpsay);
				}
				real g2=parms->powfs[ipowfs].dither_glpf;
				if(simu->wfsflags[ipowfs].ogout*g2<1){//not enough accumulations yet.
					g2=1./(simu->wfsflags[ipowfs].ogout);
				}
				gradoff_acc(simu, ipowfs);
				if(g2<1){
					info("Applying LPF with gain %.2f to i0/gx/gy update at update cycle %d\n", g2, simu->wfsflags[ipowfs].ogout);
				}
				for(int jwfs=0; jwfs<nwfs; jwfs++){
					int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
					dither_t* pd=simu->dither[iwfs];
					//Scale the output due to accumulation
					//TODO: remove the LPF which is not useful.
					//TODO: combine pd->i0 with intstat->i0
					//TODO: accumulate directly to intstat->i0 instead of to pd->imx
					for(int isa=0; isa<nsa; isa++){
						dadd(&P(intstat->i0, isa, jwfs), 1-g2, P(pd->i0, isa), scale1*g2);
						if(parms->powfs[ipowfs].dither==1){
							dadd(&P(intstat->gx, isa, jwfs), 1-g2, P(pd->gx, isa), scale1*g2);
							dadd(&P(intstat->gy, isa, jwfs), 1-g2, P(pd->gy, isa), scale1*g2);
						}
					}
					dcellzero(pd->i0);
					dcellzero(pd->gx);
					dcellzero(pd->gy);
					if(parms->powfs[ipowfs].dither==1){
						if(parms->dbg.gradoff_reset==0){
							if(jwfs==0) info("Step %5d: powfs%d reducing gradoff by grad of i0.\n", isim, ipowfs);
							dmat* goff=0;
							/*Compute the gradient of i0 using old gradient
							algorithm and subtract from the gradient offset to
							prevent sudden jump of gradient measurement.*/
							shwfs_grad(&goff, PCOL(intstat->i0, jwfs), parms, powfs, iwfs, parms->powfs[ipowfs].phytype_sim);
							dadd(&P(simu->gradoff, iwfs), 1, goff, -1);
							dfree(goff);
							if(parms->powfs[ipowfs].dither_glpf!=1){
								warning("when dbg.gradoff_reset is enabled, dither_glpf should be 1.\n");
							}
						} else if(parms->dbg.gradoff_reset==1){
							if(jwfs==0) info("Step %5d: powfs%d resetting gradoff to 0.\n", isim, ipowfs);
							dzero(P(simu->gradoff, iwfs));
						} else if(parms->dbg.gradoff_reset==2){
							if(jwfs==0) info("Step %5d: powfs%d reducing gradoff by its average.\n", isim, ipowfs);
							int nacc=simu->gradoffisim-simu->gradoffisim0;
							if(jwfs==0) info("Step %5d: powfs%d gradoffacc is scaled by 1/%d\n", isim, ipowfs, nacc);
							dscale(P(simu->gradoffacc, iwfs), 1./nacc);
							dadd(&P(simu->gradoff, iwfs), 1, P(simu->gradoffacc, iwfs), -1);
						}
					}
				}
				if(ptype2==PTYPE_MF&&parms->powfs[ipowfs].llt){
					//for LGS only. tip/tilt and focus drift control is needed for matched filter with either dithering or sodium fitting
					//2022-07-12: moved before the next block because sodium_fit_wrap() modifies intstat->i0 in place.
					//obsolete: use sodium fit for LGS instead. It works much better and avoids cog bootstrap problem and can work without TWFS.
					dmat *i0grad=0;
					for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
						int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
						shwfs_grad(&i0grad, PCOL(intstat->i0, jwfs), parms, powfs, iwfs, PTYPE_COG);
						if(parms->save.dither){
							writebin(i0grad, "extra/wfs%d_i0grad_%d", iwfs, isim);
						}
						wfsgrad_tt_drift(i0grad, simu, P(parms->sim.eplo, 0), iwfs, 0);
						wfsgrad_focus_drift(i0grad, simu, 1, iwfs, 0);
					}
					dfree(i0grad);
				}
				//i0 is collected
				if(parms->powfs[ipowfs].dither==-1){
					if(parms->powfs[ipowfs].llt){//LGS, require fiting sodium profile
						dmat* sodium=0;
						dcell* grad=0;
						//don't need gradient output for matched filter
						dcell** pgrad=(ptype2==PTYPE_COG)?&grad:NULL;
						dcell** pi0=(ptype2==PTYPE_MF)?&intstat->i0:NULL;
						dcell** pgx=(ptype2==PTYPE_MF)?&intstat->gx:NULL;
						dcell** pgy=(ptype2==PTYPE_MF)?&intstat->gy:NULL;
						//1 iteration of cog/mtch is necessary to get sodium profile, 3 iterations of mtche is necessary for gradient of i0
						int niter=parms->powfs[ipowfs].llt->na_fit_maxit?parms->powfs[ipowfs].llt->na_fit_maxit:(ptype2==PTYPE_MF?1:3);

						if(parms->save.dither){
							writebin(intstat->i0, "extra/powfs%d_i0i_%d", ipowfs, isim);
						}
						sodium_fit_wrap(&sodium, pgrad, pi0, pgx, pgy, intstat->i0, parms, powfs, ipowfs,
							recon->r0, recon->L0, niter, 1);
						if(parms->save.extra){
							writebin(sodium, "extra/powfs%d_fit_sodium_%d", ipowfs, isim);
							if(grad) writebin(grad, "extra/powfs%d_fit_grad_%d", ipowfs, isim);
						}
						if(parms->save.dither){
							if(pi0) writebin(intstat->i0, "extra/powfs%d_i0o_%d", ipowfs, isim);
							if(pgx) writebin(intstat->gx, "extra/powfs%d_gxo_%d", ipowfs, isim);
							if(pgy) writebin(intstat->gy, "extra/powfs%d_gyo_%d", ipowfs, isim);
						}

						if(ptype2==PTYPE_COG){//project pgrad to TWFS corrected modes
							dcelladd(pgrad, 1, powfs[ipowfs].gradoff, -1);
							dbg("project pgrad to TWFS corrected modes\n");
							dcell *zm=0;
							dcellmm(&zm, recon->RRlgs, *pgrad, "nn", 1);
							dcellzero(*pgrad);
							dcellmm(pgrad, recon->GRlgs, zm, "nn", 1);
							if(parms->save.dither){
								writebin(zm, "extra/powfs%d_zm_%d", ipowfs, isim);
								writebin(*pgrad, "extra/powfs%d_fit_grad_proj_%d", ipowfs, isim);
							}
							dcellfree(zm);
						}

						for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
							int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
							if(ptype2==PTYPE_COG){
								if(jwfs==0) dbg("in cog mode, gradoff+=(g_ncpa-grad)\n");
								dadd(&P(simu->gradoff, iwfs), 1, P(grad, jwfs), -1);
								//prevent gradoff from accumulating tip/tilt or focus mode if any. no need to do drift control.
								wfsgrad_tt_drift(P(simu->gradoff, iwfs), simu, 0, iwfs, 1);
								wfsgrad_focus_drift(P(simu->gradoff, iwfs), simu, 0, iwfs, 1);
							} else if(ptype2==PTYPE_MF){
								if(jwfs==0) dbg("in cmf mode, gradoff is reset to 0, and ncpa is used to create i0 with new sodium profile\n");
								//since we are building ideal matched filter with
								//the correct gradoff and sodium profile. no need
								//to use gradient reference vector.
								dzero(P(simu->gradoff, iwfs));
							}
						}
						dcellfree(grad);
						dfree(sodium);
					}else{//NGS, just use i0
						if(ptype2==PTYPE_MF){
							if(simu->gradoff) {
								for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
									int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
									dzero(P(simu->gradoff, iwfs));
								}
							}
						}else{
							warning("Not handled: powfs%d: dither=-1, phytype_sim2=%d\n", ipowfs, ptype2);
						}
					}
				}

				//there is no need to reset trombone error signal
				if((parms->save.gradoff||parms->save.dither)&&parms->dbg.gradoff_reset!=1){
					writebin(simu->gradoff, "extra/gradoff_%d_dither", isim);
				}
				if(parms->dbg.gradoff_reset==2){
					dcellzero(simu->gradoffacc);
					simu->gradoffisim0=isim;
				}
				if(!simu->gradoffdrift&&parms->powfs[ipowfs].dither==1){
					//Use gradoff before adjustment is not good. There are difference between i and i0.
					if(parms->powfs[ipowfs].phytype_sim==PTYPE_MF){
						info2("Step %5d: powfs%d set gradoffdrift to cog of initial i0\n", isim, ipowfs);
					} else{
						info2("Step %5d: powfs%d set gradoffdrift to cog of created i0 + gradoff\n", isim, ipowfs);
					}
					simu->gradoffdrift=dcellnew(nwfs, 1);
					for(int jwfs=0; jwfs<nwfs; jwfs++){
						int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
						shwfs_grad(&P(simu->gradoffdrift, iwfs), PCOLR(powfs[ipowfs].intstat->i0, jwfs), parms, powfs, iwfs, PTYPE_COG);
						if(parms->powfs[ipowfs].phytype_sim!=PTYPE_MF){
							dadd(&P(simu->gradoffdrift, iwfs), 1, P(simu->gradoff, iwfs), 1);//added on 2021-08-11
						}
					}
				}

				if(parms->powfs[ipowfs].phytype_sim!=ptype2){
					//the following parms changes need to be moved to simu. It affects the next seed.
					parms->powfs[ipowfs].phytype_sim=ptype2;
					parms->powfs[ipowfs].phytype_recon=parms->powfs[ipowfs].phytype_sim;
					info2("Step %5d: powfs %d changed to %s\n", isim, ipowfs,
						parms->powfs[ipowfs].phytype_sim==PTYPE_MF?"matched filter":"CoG");
				}
				//Generating matched filter
				if(parms->powfs[ipowfs].neareconfile||parms->powfs[ipowfs].phyusenea){
					dbg("Step %5d: powfs %d disable neareconfile and phyusenea\n", isim, ipowfs);
					parms->powfs[ipowfs].neareconfile=NULL;
					parms->powfs[ipowfs].phyusenea=0;
				}
				if(ptype2==PTYPE_MF){
					parms->powfs[ipowfs].phytype_recon=PTYPE_MF;//Make sure nea is used for reconstruction.
					genmtch(parms, powfs, ipowfs);
#if USE_CUDA
					if(parms->gpu.wfs){
						gpu_wfsgrad_update_ref(parms, powfs, ipowfs);
					}
#endif
				}
				if(parms->save.dither>1){
					writebin(intstat->i0, "extra/powfs%d_i0_%d", ipowfs, isim);
					if(ptype2==PTYPE_MF){
						writebin(intstat->gx, "extra/powfs%d_gx_%d", ipowfs, isim);
						writebin(intstat->gy, "extra/powfs%d_gy_%d", ipowfs, isim);
						writebin(intstat->mtche, "extra/powfs%d_mtche_%d", ipowfs, isim);
						writebin(powfs[ipowfs].sanea, "extra/powfs%d_sanea_%d", ipowfs, isim);
					}
				}

				if(!parms->powfs[ipowfs].lo&&parms->recon.alg==RECON_MVR){//no need to update LSR.
					simu->tomo_update=2;
				}
				if(parms->powfs[ipowfs].dither_gdrift>0&&parms->powfs[ipowfs].dither==1){
					wfsgrad_sa_drift(simu, ipowfs);
				}
			} else {
				if(parms->powfs[ipowfs].phytype_sim!=ptype2){
					error("Does not support switching to CoG.\n");
				}
				//For CoG gain
				for(int jwfs=0; jwfs<nwfs; jwfs++){
					int iwfs=P(parms->powfs[ipowfs].wfs, jwfs);
					dither_t* pd=simu->dither[iwfs];
					const int ogsingle=!pd->gg0||parms->powfs[ipowfs].dither_ogsingle;
					const int ng=ogsingle?1:(powfs[ipowfs].saloc->nloc*parms->powfs[ipowfs].ng);
					const int nd=PN(pd->a2mv);
					if(!P(simu->gradscale, iwfs)){
						P(simu->gradscale, iwfs)=dnew(ng, 1);
						dset(P(simu->gradscale, iwfs), parms->powfs[ipowfs].gradscale);
					}
					if(nd>1 && !P(simu->gradscale2, iwfs)){
						P(simu->gradscale2, iwfs)=dnew(nd, 1);
						dset(P(simu->gradscale2, iwfs), 1);
					}
					real mgold=dsum(P(simu->gradscale, iwfs))/ng;
					real mgnew=0;//,mgnew2=0;
					const char* ogtype=0;
					//gg0 is output/input of dither dithersig.
					if(ogsingle){//single gain for all subapertures. For Pyramid WFS
						ogtype="globally";
						real gerr=pd->a2m/pd->a2me;
#define HIA_G_UPDATA 0
#if HIA_G_UPDATE //HIA method.
						real adj=parms->powfs[ipowfs].dither_gog*mgold*(gerr-1);
						if(fabs(adj)>0.1) adj*=0.1/fabs(adj);
						while(adj+mgold<0){//prevent negative value
							adj*=0.5;
						}
						dadds(P(simu->gradscale, iwfs), adj);
						mgnew=mgold+adj;
#else
						real adj=pow(gerr, parms->powfs[ipowfs].dither_gog);
						CLIP(adj, 0.7, 1.5);//0.5, 1.5 //clip the adjustment to prevent divergence.
						dscale(P(simu->gradscale, iwfs), adj);
						mgnew=mgold*adj;
#endif
						if(nd>1){//multi-mode dithering
							dmat *gs2=P(simu->gradscale2, iwfs);
							for(int id=0; id<nd; id++){
								const real g2err=P(pd->a2mv, id)/P(pd->a2mev, id);
								/**
								 * A low pass filter is used to update gs2
								 * instead of an integrator like gradscale
								 * because the gradient and a2mev is not
								 * affected by its value.
								*/
								const real g=parms->powfs[ipowfs].dither_gog;
#if HIA_G_UPDATE //HIA method.
									P(gs2, id)+=(g2err-P(gs2, id))*g;
#else
									P(gs2, id)*=pow(g2err/P(gs2, id), g);
#endif
							}
							//mgnew2=P(gs2, nd-1);
							//dshow(gs2, "gradscale2");
							//info2("Scale ratios are High: %g Low: %g\n", P(simu->gradscale2, iwfs), P(P(simu->gradscale, iwfs), 0));
						}

					} else{//separate gain for each gradient. For shwfs.
						ogtype="on average";
						dscale(pd->gg0, scale1); //Scale value at end of accumulation
						for(long ig=0; ig<ng; ig++){
							if(P(pd->gg0, ig)>0.01){//skip weakly determined subapertures.
#if HIA_G_UPDATE //HIA method.
								real adj=parms->powfs[ipowfs].dither_gog*mgold*(1.-P(pd->gg0, ig));
								while(adj+P(P(simu->gradscale, iwfs), ig)<0){
									adj*=0.5;
								}
								P(P(simu->gradscale, iwfs), ig)+=adj;
#else
								real adj=pow(P(pd->gg0, ig), -parms->powfs[ipowfs].dither_gog);
								if(adj>1.5) adj=1.5; else if(adj<0.5) adj=0.5;
								P(P(simu->gradscale, iwfs), ig)*=adj;
#endif
							}
						}
						mgnew=dsum(P(simu->gradscale, iwfs))/ng;
						dzero(pd->gg0);
					}
					/*if(mgnew2){
						info2("Step %5d: wfs %d estimate/dither=%.2f, updated CoG gain=(%5.2f, %5.2f) %s\n",
							isim, iwfs, pd->a2me/pd->a2m, mgnew, mgnew2, ogtype);
					}else*/{
						info2("Step %5d: wfs %d estimate/dither=%.2f, updated CoG gain=%5.2f %s\n",
							isim, iwfs, pd->a2me/pd->a2m, mgnew, ogtype);
					}
					if(simu->resdither){
						int ic=simu->wfsflags[ipowfs].pllout-1;
						P(P(simu->resdither, iwfs), 3, ic)=mgnew;
					}
					//adjust WFS measurement dither dithersig by gain adjustment. used for dither t/t removal from gradients.
					pd->a2me*=(mgnew/mgold);//Adjust for updated gain
					dcellscale(powfs[ipowfs].sanea, pow(mgnew/mgold, 2));
					if(parms->save.dither){
						writebin(P(simu->gradscale, iwfs), "extra/gradscale_wfs%d_%d", iwfs, isim);
						if(P(simu->gradscale2, iwfs)){
							writebin(P(simu->gradscale2, iwfs), "extra/gradscale2_wfs%d_%d", iwfs, isim);
						}
					}
				}
			}
		}
	}
}
/**
   TWFS has output. Accumulate result to simu->gradoff. It is put in wfsgrad.c
   instead of recon.c to avoid race condition because it updates simu->gradoff.
*/
void wfsgrad_twfs_recon(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int itpowfs=parms->itpowfs;
	if(simu->wfsflags[itpowfs].gradout){
		info2("Step %5d: TWFS[%d] has output with gain %g\n", simu->wfsisim, itpowfs, simu->eptwfs);
		gradoff_acc(simu, parms->ilgspowfs);//todo: improve ipowfs index.
		const int nlayer=NY(recon->GRall);
		dcell* Rmod=0;
		//Build radial mode error using closed loop TWFS measurements from this time step.
		dcellmm(&Rmod, recon->RRtwfs, simu->gradcl, "nn", 1);
		if(simu->wfsflags[itpowfs].gradout<5&&parms->itwfssph>-1){
			dbg("Step %5d: TWFS output %d spherical mode (%d) gain is boosted from %g to %g\n",
				simu->wfsisim, simu->wfsflags[itpowfs].gradout, parms->itwfssph, parms->sim.eptwfs, parms->sim.eptsph);
			for(int ilayer=0; ilayer<nlayer; ilayer++){
				P(P(Rmod, ilayer), parms->itwfssph)*=(parms->sim.eptsph/simu->eptwfs);
			}
		}
		zfarr_push(simu->save->restwfs, simu->wfsflags[itpowfs].gradout-1, Rmod);

		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			if(parms->powfs[ipowfs].llt){
				for(int ilayer=0; ilayer<nlayer; ilayer++){
					dmm(&P(simu->gradoff, iwfs), 1, P(recon->GRall, iwfs, ilayer), P(Rmod, ilayer), "nn", -simu->eptwfs);
				}
				plot_gradoff(simu, iwfs);
			}
		}
		if(parms->save.gradoff){
			writebin(simu->gradoff, "extra/gradoff_%d_twfs", simu->wfsisim);
		}
		dcellfree(Rmod);
	}
}
static unsigned int petal_bg_status=0;//if set, indicates that the thread wfsgrad_peta_bg is finished
/**
 * @brief Background profess to run the slow reconstruction process
 * 
 * @param simu 
 * @return void* 
 */
void *wfsgrad_petal_bg(sim_t *simu){
	const int nrep=8;
	for(int ir=0; ir<2; ir++){
		dmat *phib=ir==1?P(simu->petal_m, 0):NULL;//2nd step use first step as input.
		petal_solve(NULL, &P(simu->petal_m, ir), simu->recon->petal[ir], P(simu->petal_i0, ir, 1), phib, nrep);
	}
	petal_bg_status=1;
	return NULL;
}
/**
 * Reconstruct petal modes. The function is now split into three parts:
 * Part 1: accumuate i0 at every time step and set do_petal if a reconstruction is due.
 * Part 2: launch a new thread to do the petal reconstruction.
 * Part 3: join the thread and do post-processing to update the reference vectors.
*/
void wfsgrad_petal_recon(sim_t *simu){
	const parms_t *parms=simu->parms;

	if(!simu->petal_i0){
		simu->petal_i0=dccellnew(2,2);//first column for accumulation, second column for storage
	}
	if(!simu->petal_m){
		simu->petal_m=dcellnew(3,1);
	}
	int isim=simu->wfsisim;
	int do_petal=0;
	for(int ir=0; ir<2; ir++){
		int ipowfs=0;
		if(ir==0&&parms->ittfpowfs!=-1){
			ipowfs=parms->ittfpowfs;
		} else if(ir==1&&parms->ittpowfs!=-1){
			ipowfs=parms->ittpowfs;
		}else{
			continue;
		}
		if(isim<parms->powfs[ipowfs].phystep || isim<parms->recon.petalstep) continue;
		int iwfs=P(parms->powfs[ipowfs].wfs, 0);
		dcelladd(&P(simu->petal_i0, ir, 0), 1, P(simu->ints, iwfs), 1);
		if((isim+1-parms->powfs[ipowfs].step)%parms->recon.petaldtrat==0){
			//normalization is not necessary
			if(parms->save.recon){
				writebin(P(simu->petal_i0, ir, 0), "petal_i0_%d_%d", ir, isim);
			}
			/*if(parms->plot.run){
				draw("Petal", (plot_opts){ .image=P(P(simu->petal_i0, ir), 0) }, "PSF of first subaperture", "x", "y", "powfs %d", ipowfs);
			}*/

			dcelladd(&P(simu->petal_i0, ir, 1), 0, P(simu->petal_i0, ir, 0), 1);
			dcellzero(P(simu->petal_i0, ir, 0));
			do_petal=1;
		}
	}
	static pthread_t thread=0;
	static int petal_isim=0;
	if(thread && (petal_bg_status>0 || do_petal || isim+1==parms->sim.end)){//join previous thread and finish up
		void *ans;
		pthread_join(thread, &ans); thread=0; petal_bg_status=0;
		//post-processing petalling results
		dadds(P(simu->petal_m, 1), -dmean(P(simu->petal_m, 1)));//remove piston
		dshow(P(simu->petal_m, 1), "Step%6d (from %d): petal output:", isim, petal_isim);
		real rad2m=parms->powfs[parms->ittfpowfs].wvlmean/TWOPI;//radian to m
		dscale(P(simu->petal_m, 1), rad2m);//convert to m
		int idm=parms->idmground;
		dcellzero(simu->dmtmp);//this cannot be done in the thread as dmtmp is used by filter().
		real gain=parms->recon.petaldtrat==1?0.5:1;//gain of 1 can be used if peltadtrat>1
		dspmm(&P(simu->dmtmp, idm), simu->recon->apetal, P(simu->petal_m, 1), "nn", 1);
		if(1){
			warning_once("Remove p/t/t from dm petal offset.\n");
			loc_remove_ptt(P(simu->dmtmp, idm), P(simu->recon->aloc, idm), NULL, NULL, 0);
		}
		for(int jwfs=0; jwfs<parms->nwfs; jwfs++){
			int jpowfs=parms->wfs[jwfs].powfs;
			if(parms->powfs[jpowfs].lo) continue;
			dcellmm(&P(simu->gradoff, jwfs), P(simu->recon->GA, jwfs, idm), P(simu->dmtmp, idm), "nn", -gain);
		}
		if(parms->plot.run&&petal_isim%parms->plot.run==0){
			int draw_single_save=draw_single; draw_single=0;
			drawopd("DM", P(simu->recon->aloc, idm), P(simu->dmtmp, idm), parms->plot.opdmax, "DM Petal Error Signal (Hi)", "x (m)", "y (m)", "Petal %d", idm);
			plot_gradoff(simu, -1);
			draw_single=draw_single_save;
		}
		servo_add(simu->dmint, simu->dmtmp, gain);
	}
	if(do_petal){//launch a thread to do the task
		petal_isim=isim;
		pthread_create(&thread, NULL, (thread_fun)wfsgrad_petal_bg, simu);
	}
}
/**
   Calls wfsgrad_iwfs() to computes WFS gradient in parallel.
   It also includes operations on Gradients before tomography.
*/
void* wfsgrad(sim_t* simu){
	real tk_start=PARALLEL==1?simu->tk_istart:myclockd();
	const parms_t* parms=simu->parms;
	if(parms->nwfs==0) return NULL;
	// call the task in parallel and wait for them to finish. It may be done in CPU or GPU.
	if(!(PARALLEL==1&&parms->tomo.ahst_idealngs!=1&&parms->gpu.wfs)){
		CALL_THREAD(simu->wfsgrad_pre, 0);
	}///else: already called by sim.c
	CALL_THREAD(simu->wfsgrad_post, 0);
	wfsgrad_dither_post(simu);//must be before wfsgrad_lgsfocus because wfsgrad_lgsfocus runs zoom integrator.
	if(parms->nlgspowfs){//high pass filter lgs focus to remove sodium range variation effect
		wfsgrad_lgsfocus(simu);
	}
	if(parms->itpowfs!=-1){
		wfsgrad_twfs_recon(simu);
	}
	if(parms->recon.petal){
		wfsgrad_petal_recon(simu);
	}
	if(parms->plot.run&&simu->wfsisim%parms->plot.run==0){
		for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
			int ipowfs=parms->wfs[iwfs].powfs;
			int jwfs=P(parms->powfs[ipowfs].wfsind, iwfs);
			drawgrad("Gcl", simu->powfs[ipowfs].saloc, PR(simu->powfs[ipowfs].saa,jwfs),P(simu->gradcl, iwfs),
				parms->plot.grad2opd, parms->powfs[ipowfs].trs, parms->plot.gmax,
				"WFS Closeloop Gradients Calibrated", "x (m)", "y (m)", "WFS %d", iwfs);
		}
	}
	//todo: split filter_fsm to per WFS.
	filter_fsm(simu);
	if(1+simu->wfsisim==parms->sim.end){
#if USE_CUDA
		if(parms->gpu.wfs){
			gpu_save_pistat(simu);
		} else
#endif
			save_pistat(simu);
	}
	simu->tk_wfs=myclockd()-tk_start;
	return NULL;
}
