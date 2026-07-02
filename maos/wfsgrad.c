/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "petal_utils.h"
#include "dither_utils.h"
#include "twfs.h"
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
		prop(&(propdata_t){.locin=aloc, .phiin=P(P(wfsopd, 0)), .locout=loc, .phiout=P(opd), 
			.alpha=alpha, .misregx=misregx, .misregy=misregy});
		dcellfree(wfsopd);
	} else{
		//project tubulence onto DM grid (global ideal)
		for(int idm=0; idm<parms->ndm; idm++){
			loc_t* loc=powfs[ipowfs].loc_dm?P(powfs[ipowfs].loc_dm, wfsind, idm):(powfs[ipowfs].loc_tel?P(powfs[ipowfs].loc_tel, wfsind):powfs[ipowfs].loc);
			const real ht=parms->dm[idm].ht;
			if(ht>hs) continue;
			prop(&(propdata_t){.mapin=P(simu->dmprojsq, idm), .locout=loc, .phiout=P(opd), .alpha=alpha, 
				.hs=hs, .thetax=parms->wfs[iwfs].thetax, .thetay=parms->wfs[iwfs].thetay, 
				.misregx=misregx, .misregy=misregy});
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
		//const real hc=nhs>1?(parms->wfs[iwfs].hc*(1.-hs/parms->wfs[iwfs].hs)):0;//effective hc (not used, but is this necessary or correct?)
		/* Now begin ray tracing. */
		if(atm&&((!parms->sim.idealwfs&&!parms->powfs[ipowfs].lo)
				||(!parms->sim.wfsalias&&parms->powfs[ipowfs].lo))){
			for(int ips=0; ips<nps; ips++){
				thread_t* wfs_prop=simu->wfs_prop_atm[iwfs+parms->nwfs*ips];
				propdata_t* wfs_propdata=&simu->wfs_propdata_atm[iwfs+parms->nwfs*ips];
				wfs_propdata->phiout=P(opd);
				if(parms->atm.dtrat>0){
					real wt;
					int iframe=atm_interp(&wt, ips, isim, parms->atm.dtrat, NX(atm), parms->atm.interp);
					wfs_propdata->alpha=atmscale*wt;
					wfs_propdata->mapin=P(atm, iframe);
				}else{
					wfs_propdata->shiftx=-P(atm, ips)->vx*dt*isim;
					wfs_propdata->shifty=-P(atm, ips)->vy*dt*isim;
					wfs_propdata->alpha=atmscale;
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
		if(simu->telfocusreal){/*Telescope focus correction*/
			real focus=-P(P(simu->telfocusreal,0),0);
			if(fabs(focus)>1e-20){
				loc_add_focus(opd, powfs[ipowfs].loc, focus);
			}
		}
		if(parms->powfs[ipowfs].llt){
			real focus=zoomfocusadj(simu, iwfs);
			if(fabs(focus)>1e-20){
				loc_add_focus_offset(opd, powfs[ipowfs].loc, focus, 
					PR(parms->powfs[ipowfs].llt->ox, wfsind), 
					PR(parms->powfs[ipowfs].llt->oy, wfsind));
			}
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
				wfs_propdata->phiout=P(opd);
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
				prop(&(propdata_t){.locin=P(recon->moao[imoao].aloc, 0), .phiin=P(dmwfs[iwfs]),
					.ptsout=powfs[ipowfs].pts, .phiout=P(opd), .alpha=-1});
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
						const real hl=P(atm, ips)->ht;
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
							atmi=P(atm, iframe);
						}else{
							vx=P(atm, ips)->vx;
							vy=P(atm, ips)->vy;
							atmi=P(atm, ips);
						}
					
						prop(&(propdata_t){.mapin=atmi, .ptsout=powfs[ipowfs].llt->pts,
							.phiout=P(lltopd), .alpha=atmscale*wt, .hs=hs, 
							.thetax=thetax, .thetay=thetay,
							.shiftx=ox-vx*isim*dt, .shifty=oy-vy*isim*dt,
							.wrap=1});
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
					dither_acc_i0(simu->parms, simu->dither[iwfs] , ints, iwfs, isim);
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

/**
   Accomplish two tasks:
   1) Use LPF'ed LGS focus measurement to drive the trombone.
   2) HPF LGS focus on gradients to remove sodium range variation effact.

   We trust the focus measurement of the LGS WFS at high temporal frequency
   which NGS cannot provide due to low frame rate. After the HPF on lgs
   gradients, our system is NO LONGER affected by sodium layer variation.

   if sim.mffocus==1: The HPF is applied to each LGS WFS independently. The remaining differential focus is negligible.

   if sim.mffocus==2: The focus measurement of all LGS is replaced by LGS averaged and HPFed focus.
*/
static void wfsgrad_lgsfocus(sim_t* simu){
	const parms_t* parms=simu->parms;
	const recon_t* recon=simu->recon;
	const int isim=simu->wfsisim;
	dcell* LGSfocus=simu->LGSfocus;//computed in wfsgrad_post from gradcl.

	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
		if(!simu->wfsflags[ipowfs].gradout||!parms->powfs[ipowfs].llt
			||isim<parms->powfs[ipowfs].step){
			continue;
		}
		real lgsfocusm=0;//LGS averaged focus
		if(parms->sim.mffocus>1){
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
			P(P(simu->LGSfocusts, iwfs), 0, isim)=P(P(LGSfocus, iwfs), 0);//save time history
			if(parms->powfs[ipowfs].zoomgain){
				//Trombone from gradients. always enable
				P(simu->zoomavg, iwfs)+=P(P(LGSfocus, iwfs), 0);//zoom averager
				P(simu->zoomavg_count, iwfs)++;
			}
			if(parms->sim.mffocus){
				real focus=0;
				if(parms->sim.mffocus==1){//Focus HPF of each LGS
					focus=P(P(LGSfocus, iwfs), 0);
				}else{//replace with global focus measurement and take HPF
					focus=lgsfocusm;
					dadd(&P(simu->gradcl, iwfs), 1, P(recon->GFall, iwfs), lgsfocusm-P(P(LGSfocus, iwfs), 0));
				}
				info_once("HPF LGS focus in gradcl\n");
				dadd(&P(simu->gradcl, iwfs), 1, P(recon->GFall, iwfs), -P(simu->lgsfocuslpf, iwfs));
				//LPF is after using the value to put it off critical path of the RTC.
				real lpfocus=parms->sim.lpfocus;
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
		if(NE(simu->fsmcmd, iwfs)){
			P(P(simu->save->fsmcmds, iwfs), 0, isim)=P(P(simu->fsmcmd, iwfs), 0);
			P(P(simu->save->fsmcmds, iwfs), 1, isim)=P(P(simu->fsmcmd, iwfs), 1);
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
					dither_acc(simu, iwfs);
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
	dither_post(simu);//must be before wfsgrad_lgsfocus because wfsgrad_lgsfocus runs zoom integrator.
	if(parms->nlgspowfs){//high pass filter lgs focus to remove sodium range variation effect
		wfsgrad_lgsfocus(simu);
	}
	if(parms->itpowfs!=-1){
		twfs_recon(simu);
	}
	if(parms->recon.petal){
		petal_recon(simu);
	}
	if(simu->gradoff){//moved here to preserve static focus in gradoff
		dcelladd(&simu->gradcl, 1, simu->gradoff, -parms->dbg.gradoff_scale);
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
