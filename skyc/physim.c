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
#include "skyc.h"
#include "parms.h"
#include "types.h"
#include "physim.h"
/**
   \file skyc/skysim_utils.c
   Utilities for skysim.c
*/

/**
   Compute Open loop NGS mode wavefront error from mode vectors.  */
real calc_rms(const dmat* mod, const dmat* mcc, int istep0){
	real rms=0;
	for(long istep=istep0; istep<mod->ny; istep++){
		rms+=dwdot(PCOL(mod, istep), mcc, PCOL(mod, istep));
	}
	return rms/(mod->ny-istep0);
}

/**
   convert mod vector to ngs WFS cloc and add to it's opd or complex pupil function.
   Notice the the coordinate for each subaperture is different for TTF.
*/
void ngsmod2wvf(cmat* wvf,            /**<[in/out] complex pupil function*/
	real wvl,           /**<[in] the wavelength*/
	const dmat* modm,     /**<[in] the NGS mode vector*/
	const powfs_s* powfs,       /**<[in] the powfs configuration*/
	int isa,              /**<[in] index of subaperture*/
	real thetax,        /**<[in] direction of WFS*/
	real thetay,        /**<[in] direction of WFS*/
	const parms_s* parms  /**<[in] the parms*/
){
	const real* mod=P(modm);
	const real wvk=TWOPI/wvl;
	real dx=powfs->dxwvf;
	real ox=powfs->saloc->locx[isa]+dx*0.5;
	real oy=powfs->saloc->locy[isa]+dx*0.5;
	int nx=powfs->nxwvf;
	assert(wvf->nx==wvf->ny);
	comp* p=P(wvf)+(wvf->nx-nx)/2*(1+wvf->nx);
	if(modm->nx==2){
		for(int iy=0; iy<nx; iy++){
			real ym=(oy+iy*dx)*mod[1];
			for(int ix=0; ix<nx; ix++){
				real x=ox+ix*dx;
				real tmp=x*mod[0]+ym;
				p[ix+wvf->nx*iy]*=EXPI(wvk*tmp);
			}
		}
	} else{
		const real hc=parms->maos.hc;
		const real hs=parms->maos.hs;
		const real scale=pow(1.-hc/hs, -2);
		const real scale1=1.-scale;
		real focus;
		if(modm->nx>5){
			focus=mod[5];
			if(!parms->maos.ahstfocus){
				focus+=mod[2]*scale1;
			}
		} else{
			focus=mod[2]*scale1;
		}
		for(int iy=0; iy<nx; iy++){
			real y=oy+iy*dx;
			for(int ix=0; ix<nx; ix++){
				real x=ox+ix*dx;
				real xy=x*y;
				real x2=x*x;
				real y2=y*y;
				real tmp=
					+x*mod[0]
					+y*mod[1]
					+focus*(x2+y2)
					+mod[2]*(-2*scale*hc*(thetax*x+thetay*y))
					+mod[3]*((x2-y2)*scale1-2*scale*hc*(thetax*x-thetay*y))
					+mod[4]*(xy*scale1-scale*hc*(thetay*x+thetax*y));
				p[ix+wvf->nx*iy]*=EXPI(wvk*tmp);
			}
		}
	}
}

/**
 * @brief Time domain physical simulation.
 * 
 * @param mresout 	Residual mode output
 * @param mideal 	Input mode
 * @param mideal_oa Input mode for on axis (?).
 * @param ngsol 	Open loop RMS WFE.
 * @param aster 	Asterism Parameters
 * @param powfs 	POWFS parameters
 * @param parms 	Parameters
 * @param idtratc 	multirate=0: common rate index. multirate=1: fastest rate index.
 * @param noisy 	0: noise free. 1: poisson and ron. 2: only poisson noise (testing)
 * @param phystart 	time step to start physical optics simulation
 * @return dmat* 	Residual RMS WFE.
 */
dmat* physim(dmat** mresout, const dmat* mideal, const dmat* mideal_oa, real ngsol,
	aster_s* aster, const powfs_s* powfs, const parms_s* parms, int idtratc, int noisy, int phystart){
	const int dtratc=P(parms->skyc.dtrats,idtratc);
	const int hasphy=(phystart>-1&&phystart<aster->nstep)?1:0;//whether physical optics is enabled
		
	const int nmod=mideal->nx;
	dmat* res=dnew(6, 1);/*Results. 
							1-2: NGS and TT modes.,
							3-4: On axis NGS and TT modes,
							4-6: On axis NGS and TT wihtout considering un-orthogonality.*/
	dmat* mreal=dnew(nmod,1);	/*modal correction at this step. */
	dmat* merr=dnew(nmod, 1);/*modal error */
	dmat* merrm=dnew(nmod, 1);//reconstruction error signal
	dmat* pmerrm=NULL;//set to merrm when there is output
	const int nstep=aster->nstep?aster->nstep:parms->maos.nstep;//number of simulation steps
	dmat* mres=dnew(nmod, nstep);//residual mode history
	dmat* rnefs=parms->skyc.rnefs;//read out noise at each sampling frequency
	dcell* zgradc=dcellnew3(aster->nwfs, 1, aster->ngrad, 0);//zernike gradient accumulation
	dcell* gradout=dcellnew3(aster->nwfs, 1, aster->ngrad, 0);
	dcell *gradslow=dcellnew3(aster->nwfs, 1, aster->ngrad, 0);//for multirate slow rate
	dmat* gradsave=0;
	if(parms->skyc.dbg>1){
		gradsave=dnew(aster->tsa*2, nstep);
	}

	servo_t* st_fast=0;
	kalman_t* kalman=0;
	int multirate=parms->skyc.multirate;
	dmat* mreal_slow=0;//servo output of slow loop
	dcell* merr_slow=0;
	int indk_fast=0;
	int divergence=0;
	int dtrat_slow=dtratc;//dtrat of slow loop
	servo_t* st_slow=NULL;
	if(multirate) {
		merr_slow=dcellnew_same(1,1,nmod, 1);
		mreal_slow=dnew(nmod, 1);
		for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
			if(!aster->wfs[iwfs].use) continue;
			if(P(aster->idtrats, iwfs)==idtratc){//active wfs at fastest rate
				indk_fast |= (1<<iwfs);
			}else{
				dtrat_slow=P(aster->dtrats,iwfs);
			}
		}
		if(dtrat_slow!=dtratc&&(parms->skyc.servo!=-1)){
			st_slow=servo_new(P(merr_slow,0), NULL, 0, parms->maos.dt*dtrat_slow, P(aster->gain, 1));
			//if(parms->skyc.dbg) dshow(st_slow->ep, "st_slow_ep");
			//lshow(aster->mdirect, "mdirect");
			//dshow(P(aster->gain, (1<<aster->nwfs)-2), "gain_slow");
		}
	}

	//aster->dtrats is only set in multirate case. dtrat of each wfs in the asterism
	if(parms->skyc.servo>0){
		st_fast=servo_new(merrm, NULL, 0, parms->maos.dt*dtratc, P(aster->gain, multirate?0:idtratc));
		//if(parms->skyc.dbg) dshow(st_fast->ep, "st_fast_ep");
	} else{
		kalman=aster->kalman[multirate?0:idtratc];
	}
	
	const long nwvl=parms->maos.nwvl;
	dcell** psf=0, ** mtche=0, ** ints=0, ** i0s=0;
	ccell* wvf=0, * wvfc=0, * otf=0;
	dmat* corr=0;
	if(hasphy){
		psf=mycalloc(aster->nwfs, dcell*);
		wvf=ccellnew(aster->nwfs, 1);
		wvfc=ccellnew(aster->nwfs, 1);
		mtche=mycalloc(aster->nwfs, dcell*);
		ints=mycalloc(aster->nwfs, dcell*);
		otf=ccellnew(aster->nwfs, 1);
		i0s=mycalloc(aster->nwfs, dcell*);
		for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
			if(!aster->wfs[iwfs].use) continue;
			const int ipowfs=aster->wfs[iwfs].ipowfs;
			const long ncomp=parms->maos.ncomp[ipowfs];
			const long nsa=parms->maos.nsa[ipowfs];
			P(wvf,iwfs)=cnew(ncomp, ncomp);
			P(wvfc,iwfs)=NULL;
			psf[iwfs]=dcellnew(nsa, nwvl);
			//cfft2plan(P(wvf,iwfs), -1);
			if(parms->skyc.multirate){
				mtche[iwfs]=aster->wfs[iwfs].pistat->mtche[(int)P(aster->idtrats,iwfs)];
			} else{
				mtche[iwfs]=aster->wfs[iwfs].pistat->mtche[idtratc];
			}
			i0s[iwfs]=aster->wfs[iwfs].pistat->i0s;
			P(otf,iwfs)=cnew(ncomp, ncomp);
			//cfft2plan(P(otf,iwfs),-1);
			//cfft2plan(P(otf,iwfs),1);
			ints[iwfs]=dcellnew(nsa, 1);
			int pixpsa=parms->skyc.pixpsa[ipowfs];
			for(long isa=0; isa<nsa; isa++){
				P(ints[iwfs],isa)=dnew(pixpsa, pixpsa);
			}
		}
	}
	zfarr* zfmerr=0;
	if(parms->skyc.dbg>1){
		zfmerr=zfarr_init(nstep, 1, "%s/merr_aster%d_dtrat%d", dirsetup, aster->iaster, dtratc);
	}
	for(int irep=0; irep<parms->skyc.navg; irep++){
		if(kalman){
			kalman_init(kalman);
		} else{
			servo_reset(st_fast);
		}
		if(st_slow){
			servo_reset(st_slow);
		}
		dcellzero(zgradc);
		dcellzero(gradout);
		if(ints){
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				dcellzero(ints[iwfs]);
			}
		}
		for(int istep=0; istep<nstep; istep++){
			//if(istep<1000) dshow(mreal, "mreal total");
			memcpy(P(merr), PCOL(mideal, istep), nmod*sizeof(real));
			dadd(&merr, 1, mreal, -1);/*form NGS mode error; */
			memcpy(PCOL(mres, istep), P(merr), sizeof(real)*nmod);
			pmerrm=NULL;
			if(istep>=parms->skyc.evlstart){/*performance evaluation*/
				real res_ngs=dwdot(P(merr), parms->maos.mcc, P(merr));
				if(res_ngs>ngsol*100){
					//dfree(res); res=NULL;
					//dbg("Loop is diverging at step %d\n", istep);
					divergence=1;
					break;
				}
				//field averaged performance
				P(res,0)+=res_ngs;
				P(res,1)+=dwdot(P(merr), parms->maos.mcc_tt, P(merr));
				real dot_oa=dwdot(P(merr), parms->maos.mcc_oa, P(merr));
				real dot_res_ideal=dwdot(P(merr), parms->maos.mcc_oa, PCOL(mideal, istep));
				real dot_res_oa=0;
				for(int imod=0; imod<nmod; imod++){
					dot_res_oa+=P(merr,imod)*P(mideal_oa, imod, istep);
				}
				P(res,2)+=dot_oa-2*dot_res_ideal+2*dot_res_oa;
				P(res,4)+=dot_oa;
				//on axis performance
				real dot_oa_tt=dwdot(P(merr), parms->maos.mcc_oa_tt, P(merr));
				/*Notice that mcc_oa_tt2 is 2x5 marix. */
				real dot_res_ideal_tt=dwdot(P(merr), parms->maos.mcc_oa_tt2, PCOL(mideal, istep));
				real dot_res_oa_tt=0;
				for(int imod=0; imod<2; imod++){
					dot_res_oa_tt+=P(merr,imod)*P(mideal_oa, imod, istep);
				}
				P(res,3)+=dot_oa_tt-2*dot_res_ideal_tt+2*dot_res_oa_tt;
				P(res,5)+=dot_oa_tt;
			}//if evl
			int indk=0;//mark wfs output.
			if(istep<phystart||phystart<0){//Geometric WFS. Ztilt
				
				for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
					if(!aster->wfs[iwfs].use) continue;
					const int ipowfs=aster->wfs[iwfs].ipowfs;
					const long ng=parms->maos.nsa[ipowfs]*2;
					dmm(&P(zgradc, iwfs), 1, P(aster->g, iwfs), merr, "nn", 1);/*grad due to residual NGS mode. */
					for(long ig=0; ig<ng; ig++){
						P(P(zgradc, iwfs), ig)+=P(aster->wfs[iwfs].ztiltout, ig, istep);
					}
					int dtrati=(multirate?P(aster->dtrats,iwfs):dtratc);
					if((istep+1)%dtrati==0){
						indk|=1<<iwfs;//has output
						dadd(&P(gradout,iwfs), 0, P(zgradc,iwfs), 1./dtrati);
						dzero(P(zgradc,iwfs));
						if(noisy){
							int idtrati=(multirate?P(aster->idtrats,iwfs):idtratc);
							dmat* nea=P(aster->wfs[iwfs].pistat->sanea,idtrati);
							//if((istep+1)==dtrati) dshow(nea, "nea[%d]", iwfs);
							for(int i=0; i<nea->nx; i++){
								P(P(gradout,iwfs),i)+=P(nea,i)*randn(&aster->rand);
							}
						}
						pmerrm=merrm;//need reconstructor output.
					}
				}
			} else{//Physical Optics WFS. Accumulate PSF intensities
				real igrad[2];
				for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
					if(!aster->wfs[iwfs].use) continue;
					const real thetax=aster->wfs[iwfs].thetax;
					const real thetay=aster->wfs[iwfs].thetay;
					const int ipowfs=aster->wfs[iwfs].ipowfs;
					const long nsa=parms->maos.nsa[ipowfs];
					ccell* wvfout=aster->wfs[iwfs].wvfout[istep];
					for(long iwvl=0; iwvl<nwvl; iwvl++){
						real wvl=parms->maos.wvl[iwvl];
						for(long isa=0; isa<nsa; isa++){
							ccp(&P(wvfc,iwfs), P(wvfout, isa, iwvl));
							/*Apply NGS mode error to PSF. */
							ngsmod2wvf(P(wvfc,iwfs), wvl, merr, powfs+ipowfs, isa, thetax, thetay, parms);
							cembed(P(wvf,iwfs), P(wvfc,iwfs), 0);
							cfft2(P(wvf,iwfs), -1); /*peak in corner. */
							cabs22d(&P(psf[iwfs],isa,iwvl), 1., P(wvf,iwfs), 1.);
						}/*isa */
					}/*iwvl */
					/*Form detector image from accumulated PSFs*/
					const int dtrati=multirate?P(aster->dtrats,iwfs):dtratc;
					const int idtrat=multirate?P(aster->idtrats,iwfs):idtratc;
					if((istep+1)%dtrati==0){/*WFS has output */
						indk|=1<<iwfs;//has output
						dcellzero(ints[iwfs]);
						for(long isa=0; isa<nsa; isa++){
							//Compute subaperture image
							for(long iwvl=0; iwvl<nwvl; iwvl++){
								real siglev=P(aster->wfs[iwfs].siglev,iwvl);
								ccpd(&P(otf,iwfs), P(psf[iwfs],isa,iwvl));
								cfft2i(P(otf,iwfs), 1); /*turn to OTF, peak in corner */
								ccwm(P(otf,iwfs), powfs[ipowfs].dtf[iwvl].nominal);
								cfft2(P(otf,iwfs), -1);
								dspmulcreal(P(P(ints[iwfs],isa)), powfs[ipowfs].dtf[iwvl].si,
									P(P(otf,iwfs)), siglev);
							}

							//Add noise 
							switch(noisy){
							case 0:/*no noise at all. */
								break;
							case 1:/*both poisson and read out noise. */
							{
								const real bkgrnd=aster->wfs[iwfs].bkgrnd*dtrati;
								addnoise(P(ints[iwfs],isa), &aster->rand, bkgrnd, bkgrnd, 0, 0, 0, P(rnefs, idtrat, ipowfs), parms->skyc.excess);
							}
							break;
							case 2:/*there is still poisson noise. */
								addnoise(P(ints[iwfs],isa), &aster->rand, 0, 0, 0, 0, 0, 0, parms->skyc.excess);
								break;
							default:
								error("Invalid noisy\n");
							}
							//Compute gradients
							igrad[0]=0;
							igrad[1]=0;
							const real pixtheta=parms->skyc.pixtheta[ipowfs];

							switch(parms->skyc.phytype){
							case 1:
								dmulvec(igrad, P(mtche[iwfs],isa), P(P(ints[iwfs],isa)), 1);
								break;
							case 2:
								dcog(igrad, P(ints[iwfs],isa), 0, 0, 3*P(rnefs, idtrat, ipowfs), 0, 0);
								igrad[0]*=pixtheta;
								igrad[1]*=pixtheta;
								break;
							case 3:
								dcorr(&corr, P(ints[iwfs],isa), P(i0s[iwfs],isa));
								dpara3(igrad, corr);
								igrad[0]*=pixtheta;
								igrad[1]*=pixtheta;
								break;
							default:
								error("Invalid phytype\n");
							}
							P(P(gradout,iwfs),isa)=igrad[0];
							P(P(gradout,iwfs),isa+nsa)=igrad[1];
						}/*isa */
						pmerrm=merrm; //mark need output
						dcellzero(psf[iwfs]);/*reset accumulation.*/
					}/*if iwfs has output*/
				}/*for wfs*/
			}/*if phystart */
			//output to mreal after using it to ensure two cycle delay.
			if(st_fast){//Type I or II control.
				servo_output(st_fast, &mreal);
			} else{//LQG control
				kalman_output(kalman, &mreal, 1, 0.5, 0);//kalman output with integrator (better than direct output)
			}
			if(dtrat_slow!=dtratc){//there is a slow loop
				if(st_slow){
					servo_output(st_slow, &mreal_slow);
				}else{
					kalman_output(kalman, &mreal_slow, 1, 0.5, 1);
				}
				for(int i=0; i<nmod; i++){
					P(mreal, i)+=P(mreal_slow, i);//directly add to corrector (temporary solution)
				}
			}
			if(indk){//has output
				if(dtratc!=dtrat_slow){
					for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
						if(!aster->wfs[iwfs].use) continue;
						if((istep+1)%P(aster->dtrats, iwfs)==0){
							dadd(&P(gradslow, iwfs), 1, P(gradout, iwfs), P(aster->dtrats, iwfs)/(real)dtrat_slow);
						}
					}
				}
				//fast rate always have output
				if(st_fast){
					dmm(&merrm, 0, P(aster->pgm, multirate?0:idtratc), gradout->m, "nn", 1);
				}else{
					kalman_update(kalman, gradout, 0);//it changes cl gradout to psol gradout
				}

				if(multirate && indk!=(indk_fast)){//slower loop has output. 
					if(st_slow){
						if(st_fast){//Both loops are integrator
							dmm(&P(merr_slow,0), 0, P(aster->pgm, 1), gradslow->m, "nn", 1);
						}else{//Fast loop is LQG. pgm not available.
							dcellzero(merr_slow);
							dcellmm(&merr_slow, P(kalman->Rlsq, 1), gradslow, "nn", 1);
						}
						servo_filter(st_slow, P(merr_slow, 0));
					}else{
						kalman_update(kalman, gradslow, 1);//it change gradslow to psol
					}
					dzero(gradslow->m);
				}
				if(zfmerr && irep==0){
					zfarr_push(zfmerr, istep, merrm);
				}
			}//if has output
			if(st_fast){
				servo_filter(st_fast, pmerrm);//do even if merrm is zero. to simulate additional latency
			}
			if(parms->skyc.dbg>1){
				memcpy(PCOL(gradsave, istep), P(gradout->m), sizeof(real)*gradsave->nx);
			}
		}/*istep; */
	}
	if(parms->skyc.dbg){
		if(parms->skyc.dbg>1){
			writebin(gradsave, "%s/grads_aster%d_dtrat%d", dirsetup, aster->iaster, dtratc);
		}
		writebin(mres, "%s/mres_aster%d_dtrat%d", dirsetup, aster->iaster, dtratc);
	}
	if(zfmerr){
		zfarr_close(zfmerr);
	}
	dfree(mreal);
	dfree(merr);
	dcellfree(merrm);
	dcellfree(zgradc);
	dcellfree(gradout);
	dcellfree(gradslow);
	dcellfree(merr_slow);
	cellfree(mreal_slow);
	dfree(gradsave);
	dfree(corr);
	if(hasphy){
		dcellfreearr(psf, aster->nwfs);
		dcellfreearr(ints, aster->nwfs);
		ccellfree(wvf);
		ccellfree(wvfc);
		ccellfree(otf);
		free(mtche);
		free(i0s);
	}
	servo_free(st_fast);
	servo_free(st_slow);
	/*dfree(mres); */
	if(mresout){
		*mresout=mres;
	} else{
		dfree(mres);
	}
	if(!divergence){
		dscale(res, 1./((nstep-parms->skyc.evlstart)*parms->skyc.navg));
	}else{
		dfree(res);
	}
	return res;
}

/**
   Save NGS WFS and other information for later use in MAOS simulations.*/
void skysim_save(const sim_s* simu, const aster_s* aster, const real* ipres, int selaster, int seldtrat, int isky){
	const parms_s* parms=simu->parms;
	const int nwvl=parms->maos.nwvl;
	char path[PATH_MAX-100];
	snprintf(path, sizeof(path), "Res%d_%d_maos/sky%d", simu->seed_maos, parms->skyc.seed, isky);
	mymkdir("%s", path);
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
		dcell* sepsf=dcelldup(aster[selaster].wfs[iwfs].pistat->psf);
		for(int ic=0; ic<sepsf->nx*sepsf->ny; ic++){
			dfftshift(P(sepsf,ic));/*put peak in center. required by MAOS. */
		}
		writebin(sepsf, "%s/pistat_wfs%d", path, iwfs+6);
		dcellfree(sepsf);
		writebin(P(aster[selaster].wfs[iwfs].pistat->sanea,seldtrat), "%s/nea_tot_wfs%d", path, iwfs+6);
		//writebin(P(aster[selaster].wfs[iwfs].pistat->sanea,seldtrat), "%s/nea_wfs%d", path, iwfs+6);
		//writebin(aster[selaster].wfs[iwfs].pistat->sanea, "%s/neafull_wfs%d", path, iwfs+6);
	}
	
	writebin(P(simu->mres,isky), "%s/mres", path);
	writebin(simu->psds, "%s/psds", path);
	char fnconf[PATH_MAX];
	snprintf(fnconf, sizeof(fnconf), "%s/base.conf", path);
	FILE* fp=fopen(fnconf, "w");

	fprintf(fp, "sim.seeds=[%d]\n", simu->seed_maos);
	fprintf(fp, "sim.end=%d\n", parms->maos.nstep);
	fprintf(fp, "sim.dt=%g\n", parms->maos.dt);
	fprintf(fp, "sim.zadeg=%g\n", parms->maos.zadeg);
	fprintf(fp, "sim.mffocus=%d\n", parms->maos.mffocus);
	fprintf(fp, "tomo.ahst_focus=%d\n", parms->maos.ahstfocus);
	fprintf(fp, "tomo.ahst_wt=3\n");
	if(parms->skyc.servo>0){
		const int idtrat=parms->skyc.multirate?0:seldtrat;
		if(PN(aster[selaster].gain, idtrat)==1){
			fprintf(fp, "sim.eplo=%g\n", P(P(aster[selaster].gain, idtrat),0));
		}else {
			fprintf(fp, "sim.eplo='gain.bin'\n");
			writebin(P(aster[selaster].gain, idtrat), "%s/gain", path);
		}
	}
	if(parms->maos.fnrange){
		fprintf(fp, "powfs0_llt.fnrange='%s'\n", parms->maos.fnrange);
	}
	fprintf(fp, "atm.r0z=%.4f\n", parms->maos.r0z);
	fprintf(fp, "atm.size=[128 128]\n");
	if(parms->maos.wddeg){
		fprintf(fp, "atm.wddeg=[");
		for(int ips=0; ips<parms->maos.nwddeg; ips++){
			fprintf(fp, "%.2f ", parms->maos.wddeg[ips]);
		}
		fprintf(fp, "]\n");
	}
	fprintf(fp, "wfs.thetax=[0 0  -33.287 -20.5725  20.5725 33.287");
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
		fprintf(fp, " %.4f", aster[selaster].wfs[iwfs].thetax*RAD2AS);
	}
	fprintf(fp, "]\n");
	fprintf(fp, "wfs.thetay=[0 35 10.8156 -28.3156 -28.3156 10.8156");
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
		fprintf(fp, " %.4f", aster[selaster].wfs[iwfs].thetay*RAD2AS);
	}
	fprintf(fp, "]\n");

	fprintf(fp, "wfs.siglev=[614 614 614 614 614 614");
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
		fprintf(fp, " %.2f", aster[selaster].wfs[iwfs].siglevtot);
	}
	fprintf(fp, "]\n");
	fprintf(fp, "wfs.wvlwts=[1 1 1 1 1 1");
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			fprintf(fp, " %.2f ", P(aster[selaster].wfs[iwfs].siglev,iwvl)
				/aster[selaster].wfs[iwfs].siglevtot);
		}
	}
	fprintf(fp, "]\n");

	if(parms->maos.npowfs<3){
		int nwfs[2]={0,0};
		for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
			nwfs[aster[selaster].wfs[iwfs].ipowfs]++;
		}
		fprintf(fp, "powfs.nwfs=[6");
		for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
			fprintf(fp, " %d", nwfs[ipowfs]);
		}
		fprintf(fp, "]\n");
	} else{
		error("Fill this out please\n");
	}
	dmat* rnefs=parms->skyc.rnefs;
	real rne=P(rnefs, seldtrat, 0);
	real bkgrnd=aster[selaster].wfs[0].bkgrnd;

	if(parms->maos.npowfs==1){
		fprintf(fp, "powfs.nwfs=[6 1]\n");
		fprintf(fp, "powfs.piinfile=[\"\" \"pistat\" ]\n");
		fprintf(fp, "powfs.neareconfile=[\"\" \"nea_tot\"]\n");
		fprintf(fp, "powfs.phyusenea=[0 1]\n");
		fprintf(fp, "powfs.dtrat=[1 %d]\n", (int)P(parms->skyc.dtrats,seldtrat));
		fprintf(fp, "powfs.bkgrnd=[0 %.2f]\n", bkgrnd);
		fprintf(fp, "powfs.rne=[3 %.2f]\n", rne);
		fprintf(fp, "powfs.phystep=[0 %ld]\n", 50+(long)P(parms->skyc.dtrats,seldtrat)*20);
		fprintf(fp, "powfs.noisy=[1 1 ]\n");
		fprintf(fp, "powfs.pixtheta=[0.8 %g]\n", parms->skyc.pixtheta[1]);
		fprintf(fp, "powfs.pixpsa=[6 %d]\n", parms->skyc.pixpsa[0]);
		fprintf(fp, "powfs.ncomp=[64 %d]\n", parms->maos.ncomp[0]);
		fprintf(fp, "powfs.nwvl=[1 %d]\n", nwvl);
		fprintf(fp, "powfs.wvl=[0.589e-6");
		for(int ip=0; ip<1; ip++){
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				fprintf(fp, " %.4g", parms->maos.wvl[iwvl]);
			}
		}
		fprintf(fp, "]\n");
	} else if(parms->maos.npowfs==2){
		fprintf(fp, "powfs.piinfile=[\"\" \"pistat\" \"pistat\"]\n");
		fprintf(fp, "powfs.neareconfile=[\"\" \"nea_tot\" \"nea_tot\"]\n");
		fprintf(fp, "powfs.phyusenea=[0 1 1]\n");
		fprintf(fp, "powfs.dtrat=[1 %d %d]\n", (int)P(parms->skyc.dtrats,seldtrat),
			(int)P(parms->skyc.dtrats,seldtrat));
		fprintf(fp, "powfs.bkgrnd=[0 %.2f %.2f]\n", bkgrnd, bkgrnd);
		fprintf(fp, "powfs.rne=[3 %.2f %.2f]\n", rne, rne);
		fprintf(fp, "powfs.phystep=[0 %ld %ld]\n",
			50+(long)P(parms->skyc.dtrats,seldtrat)*20,
			50+(long)P(parms->skyc.dtrats,seldtrat)*20);
		fprintf(fp, "powfs.noisy=[1 1 1]\n");
		fprintf(fp, "powfs.pixtheta=[0.8 %g %g]\n",
			parms->skyc.pixtheta[0]*RAD2AS,
			parms->skyc.pixtheta[1]*RAD2AS);
		fprintf(fp, "powfs.pixpsa=[6 %d %d]\n",
			parms->skyc.pixpsa[0],
			parms->skyc.pixpsa[1]);
		fprintf(fp, "powfs.ncomp=[64 %d %d]\n",
			parms->maos.ncomp[0], parms->maos.ncomp[1]);
		fprintf(fp, "powfs.nwvl=[1 %d %d]\n", nwvl, nwvl);
		fprintf(fp, "powfs.wvl=[0.589e-6");
		for(int ip=0; ip<2; ip++){
			for(int iwvl=0; iwvl<nwvl; iwvl++){
				fprintf(fp, " %.4g", parms->maos.wvl[iwvl]);
			}
		}
		fprintf(fp, "]\n");
		fprintf(fp, "powfs.wvlwts=[]\n");
	} else{
		error("Fill this out please\n");
	}

	fclose(fp);
	snprintf(fnconf, sizeof(fnconf), "%s/skyres.txt", path);
	fp=fopen(fnconf, "w");
	fprintf(fp, "TotAll\tNGS\tTT\n");
	fprintf(fp, "%g\t%g\t%g\n",
		sqrt(ipres[0])*1e9, sqrt(ipres[1])*1e9, sqrt(ipres[2])*1e9);
	fclose(fp);
}
