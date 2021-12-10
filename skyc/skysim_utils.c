/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "skysim_utils.h"
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
	const POWFS_S* powfs,       /**<[in] the powfs configuration*/
	int isa,              /**<[in] index of subaperture*/
	real thetax,        /**<[in] direction of WFS*/
	real thetay,        /**<[in] direction of WFS*/
	const PARMS_S* parms  /**<[in] the parms*/
){
	const real* mod=P(modm);
	const comp ik=COMPLEX(0, 2*M_PI/wvl);
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
				p[ix+wvf->nx*iy]*=cexp(ik*tmp);
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
				p[ix+wvf->nx*iy]*=cexp(ik*tmp);
			}
		}
	}
}


/**
   Time domain physical simulation.

   noisy:
   - 0: no noise at all;
   - 1: poisson and read out noise.
   - 2: only poisson noise.
*/
dmat* skysim_sim(dmat** mresout, const dmat* mideal, const dmat* mideal_oa, real ngsol,
	ASTER_S* aster, const POWFS_S* powfs,
	const PARMS_S* parms, int idtratc, int noisy, int phystart){
	const int dtratc=P(parms->skyc.dtrats,idtratc);
	int hasphy;
	if(phystart>-1&&phystart<aster->nstep){
		hasphy=1;
	} else{
		hasphy=0;
	}
	const int nmod=mideal->nx;
	dmat* res=dnew(6, 1);/*Results. 1-2: NGS and TT modes.,
			  3-4:On axis NGS and TT modes,
			  4-6: On axis NGS and TT wihtout considering un-orthogonality.*/
	dmat* mreal=NULL;/*modal correction at this step. */
	dmat* merr=dnew(nmod, 1);/*modal error */
	dcell* merrm=dcellnew(1, 1); P(merrm,0)=dnew(nmod, 1);
	dcell* pmerrm=NULL;
	const int nstep=aster->nstep?aster->nstep:parms->maos.nstep;
	dmat* mres=dnew(nmod, nstep);
	dmat* rnefs=parms->skyc.rnefs;
	dcell* zgradc=dcellnew3(aster->nwfs, 1, aster->ngs, 0);
	dcell* gradout=dcellnew3(aster->nwfs, 1, aster->ngs, 0);
	dcell *gradavg=dcellnew3(aster->nwfs, 1, aster->ngs, 0);//for multirate
	dmat* gradsave=0;
	if(parms->skyc.dbg){
		gradsave=dnew(aster->tsa*2, nstep);
	}

	servo_t* st2t=0;
	kalman_t* kalman=0;
	dcell* mpsol=0;
	int multirate=parms->skyc.multirate;
	dmat* moffsetint=0;
	dmat* moffseterr=0;
	int ind_fast=0;//index into gain and pgm for faster loop in multirate
	if(multirate) {
		moffsetint=dnew(nmod, 1);
		moffseterr=dnew(nmod, 1);
	}
	int dtrat_slow=0;//dtrat of slow loop
	//aster->dtrats is only set in multirate case. dtrat of each wfs in the asterism
	if(parms->skyc.servo>0){
		if(multirate){//only supports integrator
			//dmat* gtmp=dnew(1, 1); P(gtmp,0)=1;
			int indk=0;
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				if(P(aster->idtrats, iwfs)==idtratc){//active wfs at fastest rate
					indk |= (1<<iwfs);
				}else{
					dtrat_slow=P(aster->dtrats,iwfs);
				}
			}
			ind_fast=indk-1;
			//dbg("using gain[%d]\n", indk-1);
			st2t=servo_new(merrm, NULL, 0, parms->maos.dt*dtratc, P(aster->gain, ind_fast));
			//dfree(gtmp);
			//dbg("ind_fast=%d\n", ind_fast);
		} else{
			st2t=servo_new(merrm, NULL, 0, parms->maos.dt*dtratc, P(aster->gain,idtratc));
		}
	} else{
		if(multirate){
			kalman=aster->kalman[0];
		} else{
			kalman=aster->kalman[idtratc];
		}
	}
	if(kalman){
		kalman_init(kalman);
		mpsol=dcellnew(aster->nwfs, 1); //for psol grad.
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
	if(parms->skyc.dbg){
		zfmerr=zfarr_init(nstep, 1, "%s/skysim_merr_aster%d_dtrat%d", dirsetup, aster->iaster, dtratc);
	}
	for(int irep=0; irep<parms->skyc.navg; irep++){
		if(kalman){
			kalman_init(kalman);
		} else{
			servo_reset(st2t);
		}
		dcellzero(zgradc);
		dcellzero(gradout);
		if(ints){
			for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
				dcellzero(ints[iwfs]);
			}
		}
		for(int istep=0; istep<nstep; istep++){
			memcpy(P(merr), PCOL(mideal, istep), nmod*sizeof(real));
			dadd(&merr, 1, mreal, -1);/*form NGS mode error; */
			memcpy(PCOL(mres, istep), P(merr), sizeof(real)*nmod);
			if(mpsol){//collect averaged modes for PSOL.
				for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
					dadd(&P(mpsol,iwfs), 1, mreal, 1);
				}
			}
			pmerrm=0;
			if(istep>=parms->skyc.evlstart){/*performance evaluation*/
				real res_ngs=dwdot(P(merr), parms->maos.mcc, P(merr));
				if(res_ngs>ngsol*100){
					//dfree(res); res=NULL;
					dbg("Loop is diverging at step %d\n", istep);
					break;
				}
				{
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
				}
				{
					real dot_oa_tt=dwdot(P(merr), parms->maos.mcc_oa_tt, P(merr));
					/*Notice that mcc_oa_tt2 is 2x5 marix. */
					real dot_res_ideal_tt=dwdot(P(merr), parms->maos.mcc_oa_tt2, PCOL(mideal, istep));
					real dot_res_oa_tt=0;
					for(int imod=0; imod<2; imod++){
						dot_res_oa_tt+=P(merr,imod)*P(mideal_oa, imod, istep);
					}
					P(res,3)+=dot_oa_tt-2*dot_res_ideal_tt+2*dot_res_oa_tt;
					P(res,5)+=dot_oa_tt;
				}
			}//if evl

			if(istep<phystart||phystart<0){
				/*Ztilt, noise free simulation for acquisition. */
				dmm(&zgradc->m, 1, aster->gm, merr, "nn", 1);/*grad due to residual NGS mode. */
				for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
					const int ipowfs=aster->wfs[iwfs].ipowfs;
					const long ng=parms->maos.nsa[ipowfs]*2;
					for(long ig=0; ig<ng; ig++){
						P(P(zgradc, iwfs), ig)+=P(aster->wfs[iwfs].ztiltout, ig, istep);
					}
				}

				for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
					int dtrati=(multirate?P(aster->dtrats,iwfs):dtratc);
					if((istep+1)%dtrati==0){
						dadd(&P(gradout,iwfs), 0, P(zgradc,iwfs), 1./dtrati);
						dzero(P(zgradc,iwfs));
						if(noisy){
							int idtrati=(multirate?P(aster->idtrats,iwfs):idtratc);
							dmat* nea=P(aster->wfs[iwfs].pistat->sanea,idtrati);
							for(int i=0; i<nea->nx; i++){
								P(P(gradout,iwfs),i)+=P(nea,i)*randn(&aster->rand);
							}
						}
						pmerrm=merrm;//record output.
					}
				}
			} else{
				/*Accumulate PSF intensities*/
				for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
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
							ngsmod2wvf(P(wvfc,iwfs), wvl, merr, powfs+ipowfs, isa,
								thetax, thetay, parms);
							cembedc(P(wvf,iwfs), P(wvfc,iwfs), 0, C_FULL);
							cfft2(P(wvf,iwfs), -1);
							/*peak in corner. */
							cabs22d(&P(psf[iwfs],isa,iwvl), 1., P(wvf,iwfs), 1.);
						}/*isa */
					}/*iwvl */
				}/*iwfs */

				/*Form detector image from accumulated PSFs*/
				real igrad[2];
				for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
					int dtrati=dtratc, idtrat=idtratc;
					if(multirate){//multirate
						idtrat=P(aster->idtrats,iwfs);
						dtrati=P(aster->dtrats,iwfs);
					}
					if((istep+1)%dtrati==0){/*has output */
						dcellzero(ints[iwfs]);
						const int ipowfs=aster->wfs[iwfs].ipowfs;
						const long nsa=parms->maos.nsa[ipowfs];
						for(long isa=0; isa<nsa; isa++){
							for(long iwvl=0; iwvl<nwvl; iwvl++){
								real siglev=P(aster->wfs[iwfs].siglev,iwvl);
								ccpd(&P(otf,iwfs), P(psf[iwfs],isa,iwvl));
								cfft2i(P(otf,iwfs), 1); /*turn to OTF, peak in corner */
								ccwm(P(otf,iwfs), powfs[ipowfs].dtf[iwvl].nominal);
								cfft2(P(otf,iwfs), -1);
								dspmulcreal(P(P(ints[iwfs],isa)), powfs[ipowfs].dtf[iwvl].si,
									P(P(otf,iwfs)), siglev);
							}

							/*Add noise and apply matched filter. */
#if _OPENMP >= 200805 
#pragma omp critical 
#endif
							switch(noisy){
							case 0:/*no noise at all. */
								break;
							case 1:/*both poisson and read out noise. */
							{
								real bkgrnd=aster->wfs[iwfs].bkgrnd*dtrati;
								addnoise(P(ints[iwfs],isa), &aster->rand, bkgrnd, bkgrnd, 0, 0, 0, P(rnefs, idtrat, ipowfs), parms->skyc.excess);
							}
							break;
							case 2:/*there is still poisson noise. */
								addnoise(P(ints[iwfs],isa), &aster->rand, 0, 0, 0, 0, 0, 0, parms->skyc.excess);
								break;
							default:
								error("Invalid noisy\n");
							}

							igrad[0]=0;
							igrad[1]=0;
							real pixtheta=parms->skyc.pixtheta[ipowfs];

							switch(parms->skyc.phytype){
							case 1:
								dmulvec(igrad, P(mtche[iwfs],isa), P(P(ints[iwfs],isa)), 1);
								break;
							case 2:
								dcog(igrad, P(ints[iwfs],isa), 0, 0, 0, 3*P(rnefs, idtrat, ipowfs), 0);
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
						pmerrm=merrm;
						dcellzero(psf[iwfs]);/*reset accumulation.*/
					}/*if iwfs has output*/
				}/*for wfs*/
			}/*if phystart */
			//output to mreal after using it to ensure two cycle delay.
			if(st2t){//Type I or II control.
				if(P(st2t->mint,0)){//has output.
					dcp(&mreal, P(P(st2t->mint,0),0));
				}
			} else{//LQG control
				kalman_output(kalman, &mreal, 0, 1);
			}
			if(parms->skyc.servo<0){//LQG control
				int indk=0;
				//Form PSOL grads and obtain index to LQG M
				for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
					int dtrati=(multirate?P(aster->dtrats,iwfs):dtratc);
					if((istep+1)%dtrati==0){
						indk|=1<<iwfs;
						dmm(&P(gradout,iwfs), 1, P(aster->g,iwfs), P(mpsol,iwfs), "nn", 1./dtrati);
						dzero(P(mpsol,iwfs));
					}
				}
				if(indk){
					kalman_update(kalman, gradout->m, indk-1);
				}
			} else{
				if(pmerrm){
					if(!multirate){//single rate
						dmm(&P(merrm,0), 0, P(aster->pgm,idtratc), gradout->m, "nn", 1);
					} else{
						int indk=0;
						for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
							if((istep+1)%P(aster->dtrats, iwfs)==0){
								indk|=1<<iwfs;
								dadd(&P(gradavg, iwfs), 1, P(gradout, iwfs), 1);
								//if(istep<100) info("step %d accumulating wfs %d\n", istep, iwfs);
							}
						}
						dzero(P(merrm, 0));
						
						if(indk!=(ind_fast+1)){//slower loop when there is faster loop
							//slower loop (all wfs active) runs a cascaded integrator as offset to the faster loop
							for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
								real avgf=(real)P(aster->dtrats, iwfs)/(real)dtrat_slow;
								//if(istep<100) info("step %d iwfs %d scaled by %g\n", istep, iwfs, avgf);
								dscale(P(gradavg, iwfs), avgf);
							}
							dmm(&moffseterr, 0, P(aster->pgm,indk-1), gradavg->m, "nn", 1);
							dzero(gradavg->m);
							for(int imod=0; imod<nmod; imod++){
								if(aster->mdirect && P(aster->mdirect, imod)){//directly output modes not controlled by the faster loop
									P(P(merrm,0),imod)=P(moffseterr,imod);
									//if(istep<100) info("step %d copying slow to fast mode %d\n", istep, imod);
								}else if (imod!=5){//temporary: skip focus
									P(moffsetint, imod)+=P(moffseterr, imod)*P(P(aster->gain, indk-1), imod);
								}
							}
							//if(istep<100) dshow(P(merrm, 0), "merrm");
							//if(istep<100) dshow(moffseterr, "moffseterr");
							//warning("step %d: slow loop\n", istep);
						}
						if((indk & (ind_fast+1))){//check whether fast loop is active
							dmm(&P(merrm,0), 1, P(aster->pgm,ind_fast), gradout->m, "nn", 1);
							dadd(&P(merrm, 0), 1, moffsetint, 0.15);//todo: optimize the gain
							//if(istep<100) info("step %d: fast loop\n", istep);
						}
						
						if(zfmerr){
							zfarr_push(zfmerr, istep, P(merrm,0));
						}
					}
				}//if pmerrm
				servo_filter(st2t, pmerrm);//do even if merrm is zero. to simulate additional latency
			}
			if(parms->skyc.dbg){
				memcpy(PCOL(gradsave, istep), P(gradout->m), sizeof(real)*gradsave->nx);
			}
		}/*istep; */
	}
	if(parms->skyc.dbg){
		writebin(gradsave, "%s/skysim_grads_aster%d_dtrat%d", dirsetup, aster->iaster, dtratc);
		writebin(mres, "%s/skysim_mres_aster%d_dtrat%d", dirsetup, aster->iaster, dtratc);
	}
	if(zfmerr){
		zfarr_close(zfmerr);
	}
	dfree(mreal);
	dcellfree(mpsol);
	dfree(merr);
	dcellfree(merrm);
	dcellfree(zgradc);
	dcellfree(gradout);
	dcellfree(gradavg);
	dcellfree(moffsetint);
	dcellfree(moffseterr);
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
	servo_free(st2t);
	/*dfree(mres); */
	if(mresout){
		*mresout=mres;
	} else{
		dfree(mres);
	}
	dscale(res, 1./((nstep-parms->skyc.evlstart)*parms->skyc.navg));
	return res;
}

/**
   Save NGS WFS and other information for later use in MAOS simulations.*/
void skysim_save(const SIM_S* simu, const ASTER_S* aster, const real* ipres, int selaster, int seldtrat, int isky){
	const PARMS_S* parms=simu->parms;
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
		fprintf(fp, " %.4f", aster[selaster].wfs[iwfs].thetax*206265);
	}
	fprintf(fp, "]\n");
	fprintf(fp, "wfs.thetay=[0 35 10.8156 -28.3156 -28.3156 10.8156");
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
		fprintf(fp, " %.4f", aster[selaster].wfs[iwfs].thetay*206265);
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
			parms->skyc.pixtheta[0]*206265,
			parms->skyc.pixtheta[1]*206265);
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
