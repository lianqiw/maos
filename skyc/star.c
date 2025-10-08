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
/**
   \file skyc/setup_star.c
  setup stars.
*/
#include <alloca.h>
#include "skyc.h"
#include "parms.h"
#include "types.h"
#include "star.h"
#include "photon.h"
#include "mtch.h"
/**
   Create "star" data array from star information.
*/
static star_s* setup_star_create(const parms_s* parms, dmat* coord){
	if(!coord || coord->ny==0){
		return NULL;/*there are no stars available. */
	}
	int nstar=coord->ny;
	dmat* pc=coord;
	int nwvl=parms->maos.nwvl;
	star_s* star=mycalloc(nstar, star_s);
	real ngsgrid=parms->maos.ngsgrid*AS2RAD;
	real r2=pow(parms->skyc.patfov*AS2RAD/2., 2);
	real keepout=pow(parms->skyc.keepout*AS2RAD, 2);
	real minrad2=pow(parms->skyc.minrad*AS2RAD, 2);
	int jstar=0;
	if(nwvl+2>coord->nx){
		error("input coord has too few rows (need %d, has %ld\n", nwvl+2, coord->nx);
	}
	for(int istar=0; istar<nstar; istar++){
		int skip=0;
		if(parms->skyc.ngsalign){
			star[jstar].thetax=round(P(pc, 0, istar)/ngsgrid)*ngsgrid;
			star[jstar].thetay=round(P(pc, 1, istar)/ngsgrid)*ngsgrid;
			if(pow(star[jstar].thetax, 2)+pow(star[jstar].thetay, 2)>r2){
				star[jstar].thetax=trunc(P(pc, 0, istar)/ngsgrid)*ngsgrid;
				star[jstar].thetay=round(P(pc, 1, istar)/ngsgrid)*ngsgrid;
				if(pow(star[jstar].thetax, 2)+pow(star[jstar].thetay, 2)>r2){
					star[jstar].thetax=round(P(pc, 0, istar)/ngsgrid)*ngsgrid;
					star[jstar].thetay=trunc(P(pc, 1, istar)/ngsgrid)*ngsgrid;
					if(pow(star[jstar].thetax, 2)+pow(star[jstar].thetay, 2)>r2){
						star[jstar].thetax=trunc(P(pc, 0, istar)/ngsgrid)*ngsgrid;
						star[jstar].thetay=trunc(P(pc, 1, istar)/ngsgrid)*ngsgrid;
						if(pow(star[jstar].thetax, 2)+pow(star[jstar].thetay, 2)>r2){
							error("What?\n");
						}
					}
				}
			}
		} else{
			star[jstar].thetax=P(pc, 0, istar);
			star[jstar].thetay=P(pc, 1, istar);
		}
		for(int kstar=0; kstar<jstar; kstar++){
			if(pow(star[jstar].thetax-star[kstar].thetax, 2)
				+pow(star[jstar].thetay-star[kstar].thetay, 2)<keepout){
				if(P(pc, 2, istar)<P(star[kstar].mags,0)){
					memcpy(P(star[kstar].mags), PCOL(pc, istar)+2, sizeof(real)*nwvl);
					star[kstar].thetax=star[jstar].thetax;
					star[kstar].thetay=star[jstar].thetay;
				}
				dbg2("Star %d is too close to %d. Keep brightest.\n", jstar, kstar);
				skip=1;
				break;
			}
		}
		if(pow(star[istar].thetax, 2)+pow(star[istar].thetay, 2)<minrad2){
			dbg("Skip star at (%.1f, %.1f) because minrad=%g\n",
				star[istar].thetax*RAD2AS, star[istar].thetay*RAD2AS, parms->skyc.minrad);
			skip=1;
		}
		if(!skip){
			star[jstar].mags=dnew(nwvl, 1);
			memcpy(P(star[jstar].mags), PCOL(pc, istar)+2, sizeof(real)*nwvl);
			star[jstar].use=mycalloc(parms->maos.npowfs, int);
			jstar++;
		}
	}
	if(jstar<nstar){
	/*warning("%d stars dropped\n", nstar-jstar); */
		coord->ny=jstar;
		star=myrealloc(star, jstar, star_s);
	}
	return star;
}

/**
   Read in pistat information, used to compute matched filter, and SANEA.
*/
static void setup_star_read_pistat(sim_s* simu, star_s* star, int nstar, int seed){
	const parms_s* parms=simu->parms;
	const int npowfs=parms->maos.npowfs;
	const int nwvl=parms->maos.nwvl;
	const real ngsgrid=parms->maos.ngsgrid;
	for(int istar=0; istar<nstar; istar++){
		star_s* stari=&star[istar];
		stari->pistat=mycalloc(npowfs, pistat_s);
		const real thetax=stari->thetax*RAD2AS;/*in as */
		const real thetay=stari->thetay*RAD2AS;
		real thxnorm=thetax/ngsgrid;
		real thynorm=thetay/ngsgrid;
		long thxl=(long)floor(thxnorm);
		long thyl=(long)floor(thynorm);
		real wtx=thxnorm-thxl;
		real wty=thynorm-thyl;
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			const int msa=parms->maos.msa[ipowfs];
			const int nsa=parms->maos.nsa[ipowfs];
			dcell* avgpsf=NULL;
			dcell* neaspec=NULL;
			real wtsum=0;
			for(int ix=0; ix<2; ix++){
				real thx=(thxl+ix)*ngsgrid;
				for(int iy=0; iy<2; iy++){
					real thy=(thyl+iy)*ngsgrid;
					real wtxi=fabs(((1-ix)-wtx)*((1-iy)-wty));
					if(wtxi<0.01){
					/*dbg("skipping ix=%d,iy=%d because wt=%g\n",ix,iy,wtxi); */
						continue;
					}
					char fn[PATH_MAX];
					snprintf(fn, PATH_MAX, "%s/pistat/pistat_seed%d_sa%d_x%g_y%g",
						dirstart, seed, msa, thx, thy);
					if(!zfexist("%s",fn)){
					/*warning("%s doesn't exist\n",fn); */
					} else{
						dcell* avgpsfi=dcellread("%s", fn);
						dcelladd(&avgpsf, 1, avgpsfi, wtxi);
						dcellfree(avgpsfi);
						wtsum+=wtxi;

						snprintf(fn, PATH_MAX, "%s/neaspec/neaspec_seed%d_sa%d_x%g_y%g",
							dirstart, seed, msa, thx, thy);
						dcell* neaspeci=dcellread("%s", fn);
						dcelladd(&neaspec, 1, neaspeci, wtxi);
						dcellfree(neaspeci);
					}
				}
			}
			if(wtsum<0.01){
				warning("PISTAT is not available for (%g,%g) msa=%d\n", thetax, thetay, msa);
			}
			dcellscale(neaspec, 1./wtsum);
			dcellscale(avgpsf, 1./wtsum);
			dmat* scale=NULL;
			if(parms->skyc.bspstrehl){
				scale=dnew(nsa, nwvl);
				dmat* gx=dnew(1, 1); P(gx,0)=thxnorm;
				dmat* gy=dnew(1, 1); P(gy,0)=thynorm;
				if(nsa!=avgpsf->nx||nwvl!=avgpsf->ny){
					error("Mismatch: nsa=%d, nwvl=%d, avgpsf->nx=%ld, avgpsf->ny=%ld\n",
						nsa, nwvl, avgpsf->nx, avgpsf->ny);
				}
				for(int ic=0; ic<nsa*nwvl; ic++){
					dmat* val=dbspline_eval(simu->bspstrehl[ipowfs][ic],
						simu->bspstrehlxy, simu->bspstrehlxy,
						gx, gy);
					real ratio=P(val,0)/dsum(P(avgpsf,ic));
					/*dbg("strehl: bilinear: %g, cubic: %g\n", P(P(avgpsf,ic),0),P(val,0)); */
					if(ratio<0){
						warning("Ratio=%g is less than zero.\n", ratio);
						P(scale,ic)=1;
					} else{
						dscale(P(avgpsf,ic), ratio);
						P(scale,ic)=ratio;
					}
					dfree(val);
				}
				dfree(gx);
				dfree(gy);
			}

			stari->pistat[ipowfs].psf=avgpsf;/*PSF is in corner. */
			stari->pistat[ipowfs].neaspec=dcellnew(nsa*2, 1);
			for(int ig=0; ig<nsa*2; ig++){
				dmat* tmp=0;
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					dadd(&tmp, 0, P(neaspec,ig,iwvl), P(parms->skyc.wvlwt,iwvl));
				}
				P(stari->pistat[ipowfs].neaspec,ig)=dinterp1(simu->neaspec_dtrats, tmp, parms->skyc.dtrats, 0);
				dfree(tmp);
			}
			dcellfree(neaspec);
			stari->pistat[ipowfs].scale=scale;
			{/* skip stars with large PSF.*/
				real size=100;
				for(int ic=0; ic<avgpsf->nx*avgpsf->ny; ic++){
					real size0=dfwhm(P(avgpsf,ic));
					if(size0<size) size=size0;
				}
				if(size>3){
					stari->use[ipowfs]=-1;
				}
			}
			if(parms->skyc.dbg>2){
				writebin(avgpsf, "%s/star%d_ipowfs%d_avgpsf", dirsetup, istar, ipowfs);
				writebin(stari->pistat[ipowfs].neaspec, "%s/star%d_ipowfs%d_pistat_neaspec", dirsetup, istar, ipowfs);
			}
		}
	}
}

/**
   Compute Signal level
*/
static void setup_star_siglev(const parms_s* parms, star_s* star, int nstar){
	const long npowfs=parms->maos.npowfs;
	const long nwvl=parms->maos.nwvl;
	const real r2=pow(parms->skyc.patfov*AS2RAD/2., 2);
	dmat* rnefs=parms->skyc.rnefs;
	for(int istar=0; istar<nstar; istar++){
		star[istar].siglev=dcellnew(npowfs, 1);
		star[istar].bkgrnd=dnew(npowfs, 1);
		star[istar].siglevtot=dnew(npowfs, 1);
		/*Normalized angular distance */
		real th2=(pow(star[istar].thetax, 2)+pow(star[istar].thetay, 2))/r2;
		/*Field dependent error: nm^2=nma^2+nmb^2*theta_norm^2; */
		real imperrnm=sqrt(pow(parms->skyc.imperrnm, 2)+th2*pow(parms->skyc.imperrnmb, 2));
		for(long ipowfs=0; ipowfs<npowfs; ipowfs++){
			P(star[istar].siglev,ipowfs)=dnew(nwvl, 1);
			int iscircle=parms->maos.nsa[ipowfs]<=4?1:0;
			photon_flux(&parms->skyc.zb, P(P(star[istar].siglev,ipowfs)),
				&P(star[istar].siglevtot,ipowfs),
				&P(star[istar].bkgrnd,ipowfs),
				NULL, NULL,
				parms->maos.nwvl,
				parms->maos.wvl,
				P(star[istar].mags),
				parms->maos.dxsa[ipowfs], iscircle,
				parms->skyc.pixtheta[ipowfs],
				parms->maos.dt, parms->maos.za,
				NULL,
				imperrnm,
				P(rnefs, parms->skyc.ndtrat-1, ipowfs));
		}
	}
}

/**
   Setup matched filter for stars.
 */
static void setup_star_mtch(const parms_s* parms, powfs_s* powfs, star_s* star, int nstar, dccell* nonlin){
	const long nwvl=parms->maos.nwvl;
	const long npowfs=parms->maos.npowfs;
	dmat* rnefs=parms->skyc.rnefs;

	for(int istar=0; istar<nstar; istar++){
		if(!star[istar].minidtrat){
			star[istar].minidtrat=dnew(npowfs, 1);
		}
		dset(star[istar].minidtrat, -1);
		real radius=sqrt(pow(star[istar].thetax, 2)+pow(star[istar].thetay, 2));
		int igg=round(radius*RAD2AS/parms->maos.ngsgrid);
		//dbg("radius=%g as, igg=%d\n", radius*RAD2AS, igg);
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			const long nsa=parms->maos.nsa[ipowfs];
			const long pixpsa=parms->skyc.pixpsa[ipowfs];
			const real pixtheta=parms->skyc.pixtheta[ipowfs];
			//size of PSF
			const real sigma_theta=parms->skyc.wvlmean/parms->maos.dxsa[ipowfs];
			pistat_s* pistat=&star[istar].pistat[ipowfs];
			pistat->i0=dcellnew(nsa, nwvl);
			pistat->gx=dcellnew(nsa, nwvl);
			pistat->gy=dcellnew(nsa, nwvl);

			pistat->i0s=dcellnew(nsa, 1);
			pistat->gxs=dcellnew(nsa, 1);
			pistat->gys=dcellnew(nsa, 1);

			dcell* psf=pistat->psf;
			dcell* i0=pistat->i0;
			dcell* gx=pistat->gx;
			dcell* gy=pistat->gy;
			for(long iwvl=0; iwvl<nwvl; iwvl++){
				for(long isa=0; isa<nsa; isa++){
					real siglev=P(P(star[istar].siglev,ipowfs),iwvl);
					P(i0, isa, iwvl)=dnew(pixpsa, pixpsa);
					if(!parms->skyc.mtchfft){
						P(gx, isa, iwvl)=dnew(pixpsa, pixpsa);
						P(gy, isa, iwvl)=dnew(pixpsa, pixpsa);
					}
					psf2i0gxgy(P(i0, isa, iwvl), P(gx, isa, iwvl), P(gy, isa, iwvl),
						P(psf, isa, iwvl), powfs[ipowfs].dtf+iwvl, !parms->skyc.mtchfft);
					if(parms->skyc.mtchfft){
						info_once("Using derivative by FFT\n");
						P(gx, isa, iwvl)=derive_by_fft(P(i0, isa, iwvl), 0);
						P(gy, isa, iwvl)=derive_by_fft(P(i0, isa, iwvl), M_PI/2);
						dscale(P(gx, isa, iwvl), 1./pixtheta);
						dscale(P(gy, isa, iwvl), 1./pixtheta);
					}
					dadd(&P(pistat->i0s,isa), 1, P(i0, isa, iwvl), siglev);
					dadd(&P(pistat->gxs,isa), 1, P(gx, isa, iwvl), siglev);
					dadd(&P(pistat->gys,isa), 1, P(gy, isa, iwvl), siglev);
				}

			}
			if(parms->skyc.dbg>2){
				writebin(pistat->i0s, "%s/star%d_ipowfs%d_i0s", dirsetup, istar, ipowfs);
			}

			int ndtrat=parms->skyc.ndtrat;
			pistat->mtche=mycalloc(ndtrat, dcell*);
			pistat->sanea=dcellnew(ndtrat, 1);
			pistat->sanea0=dcellnew(ndtrat, 1);
			pistat->snr=dnew(ndtrat, 1);
			dcell* i0s=NULL; dcell* gxs=NULL; dcell* gys=NULL;

			for(int idtrat=0; idtrat<ndtrat; idtrat++){
				int dtrat=P(parms->skyc.dtrats,idtrat);
				dcelladd(&i0s, 0, pistat->i0s, dtrat);
				dcelladd(&gxs, 0, pistat->gxs, dtrat);
				dcelladd(&gys, 0, pistat->gys, dtrat);
				genmtch(&pistat->mtche[idtrat], &P(pistat->sanea,idtrat),
					i0s, gxs, gys, pixtheta, P(rnefs, idtrat, ipowfs),
					P(star[istar].bkgrnd,ipowfs)*dtrat, parms->skyc.mtchcr);
				if(parms->skyc.dbg>2){
					writebin(pistat->mtche[idtrat], "%s/star%d_ipowfs%d_mtche_dtrat%d", dirsetup, istar, ipowfs, dtrat);
					writebin(P(pistat->sanea,idtrat), "%s/star%d_ipowfs%d_sanea0_dtrat%d", dirsetup, istar, ipowfs, dtrat);
				}
				/*Add nolinearity*/
				if(nonlin){
					//add linearly not quadratically since the errors are related.
					dmat* nea_nonlin=dinterp1(P(P(nonlin,ipowfs),igg), NULL, P(pistat->sanea,idtrat), 0);
					for(int i=0; i<nsa*2; i++){
						//info("%g mas", P(P(pistat->sanea,idtrat),i)*RAD2MAS);
						P(P(pistat->sanea,idtrat),i)=sqrt(pow(P(P(pistat->sanea,idtrat),i), 2)
							+pow(P(nea_nonlin,i), 2));
						//info("-->%g mas\n", P(P(pistat->sanea,idtrat),i)*RAD2MAS);
					}
					dfree(nea_nonlin);
				}
				dcp(&P(pistat->sanea0,idtrat), P(pistat->sanea,idtrat));
				if(parms->skyc.neaaniso){
					for(int i=0; i<nsa*2; i++){
						P(P(pistat->sanea,idtrat),i)=sqrt(pow(P(P(pistat->sanea,idtrat),i), 2)
							+pow(P(P(star[istar].pistat[ipowfs].neaspec,i),idtrat), 2));
					}
				}
				if(parms->skyc.dbg>2){
					writebin(P(pistat->sanea, idtrat), "%s/star%d_ipowfs%d_sanea_dtrat%d", dirsetup, istar, ipowfs, dtrat);
				}
				//real nea=sqrt(dsumsq(P(pistat->sanea, idtrat))/2/nsa);//original version. averaged all subapertures
				//real nea=sqrt(dsumsq(P(pistat->sanea,idtrat))/2)/(nsa);// /2 to separate x/y. nsa is outside of sqrt(). averging sa to get tip/tilt
				//real nea=sqrt(dsumsq(P(pistat->sanea, idtrat)))/(nsa*2)*(nsa==1?2:5);// /2 to separate x/y. tt estimate 2 mode, ttf estimate 5 mode.
				//real nea=sqrt(dsumsq(P(pistat->sanea,idtrat))/(2*(nsa==1?1:2.5)));// /2 to separate x/y. nsa is outside of sqrt(). averging sa to get tip/tilt
				real nea=sqrt(dsumsq(P(pistat->sanea,idtrat))/(2*(nsa*nsa)));// /2 to separate x/y. nsa is outside of sqrt(). averging sa to get tip/tilt
				real snr=sigma_theta/nea;
				P(pistat->snr,idtrat)=snr;
				real snrmin=parms->skyc.multirate?P(parms->skyc.snrmin_fast, idtrat):parms->skyc.snrmin;
				if(snr>=snrmin){
					P(star[istar].minidtrat,ipowfs)=idtrat;
				}
			}//for idtrat
			if(P(star[istar].minidtrat,ipowfs)==-1){
				star[istar].use[ipowfs]=-1;
			}/*idtrat */
			dcellfree(i0s);
			dcellfree(gxs);
			dcellfree(gys);
		}/*for ipowfs */
	}/*for istar */
}
/**
   Compute Modal to gradient operator using average gradients. Similar to Z tilt
   since the mode is low order
 */
static void setup_star_gm(const parms_s* parms, powfs_s* powfs, star_s* star, int nstar){
	const long npowfs=parms->maos.npowfs;
	const real hc=parms->maos.hc;
	const real hs=parms->maos.hs;

	for(int istar=0; istar<nstar; istar++){
		star[istar].g=dcellnew(npowfs, 1);
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			const real thetax=star[istar].thetax;
			const real thetay=star[istar].thetay;
			
			P(star[istar].g,ipowfs)=mkg_ngs(powfs[ipowfs].sacent, thetax, thetay, hc, hs, 
				parms->maos.indps, parms->maos.indastig, parms->maos.indfocus, parms->maos.ahstfocus);	
		}
		if(parms->skyc.dbg>2){
			writebin(star[istar].g, "%s/star%d_g", dirsetup, istar);
		}
	}
}
long setup_star_read_ztilt(star_s* star, int nstar, const parms_s* parms, int seed){
	const real ngsgrid=parms->maos.ngsgrid;
	long nstep=0;
	TIC;tic;
	for(int istar=0; istar<nstar; istar++){
		star_s* stari=&star[istar];
		int npowfs=parms->maos.npowfs;
		stari->ztiltout=dcellnew(npowfs, 1);
		const real thetax=stari->thetax*RAD2AS;/*in as */
		const real thetay=stari->thetay*RAD2AS;

		real thxnorm=thetax/ngsgrid;
		real thynorm=thetay/ngsgrid;
		long thxl=(long)floor(thxnorm);/*Used to be real, but -0 appears. */
		long thyl=(long)floor(thynorm);
		real wtx=thxnorm-thxl;
		real wty=thynorm-thyl;
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			const int msa=parms->maos.msa[ipowfs];
			const int nsa=parms->maos.nsa[ipowfs];
			const int ng=nsa*2;
			char* fnztilt[2][2]={{NULL,NULL},{NULL,NULL}};
			char* fngoff[2][2]={{NULL, NULL}, {NULL, NULL}};
			real wtsum=0;
			for(int ix=0; ix<2; ix++){
				real thx=(thxl+ix)*ngsgrid;
				for(int iy=0; iy<2; iy++){
					real thy=(thyl+iy)*ngsgrid;
					real wtxi=fabs(((1-ix)-wtx)*((1-iy)-wty));

					if(wtxi<0.01){
					/*dbg("skipping ix=%d,iy=%d because wt=%g\n",ix,iy,wtxi); */
						continue;
					}
					fnztilt[iy][ix]=myalloca(PATH_MAX, char);
					if(parms->skyc.usephygrad){
						info_once("Using phygrad\n");
						snprintf(fnztilt[iy][ix], PATH_MAX, "%s/phygrad/phygrad_seed%d_sa%d_x%g_y%g",
							dirstart, seed, msa, thx, thy);
					} else{
						info_once("Using ztilt out\n");
						snprintf(fnztilt[iy][ix], PATH_MAX, "%s/ztiltout/ztiltout_seed%d_sa%d_x%g_y%g",
							dirstart, seed, msa, thx, thy);
					}
					fngoff[iy][ix]=myalloca(PATH_MAX, char);
					snprintf(fngoff[iy][ix], PATH_MAX, "%s/gradoff/gradoff_sa%d_x%g_y%g",
						dirstart, msa, thx, thy);
					if(!zfexist("%s",fnztilt[iy][ix])){
					//warning("%s doesnot exist\n",fnwvf[iy][ix]);
						fnztilt[iy][ix]=fngoff[iy][ix]=NULL;
					} else{
						wtsum+=wtxi;
					}
				}
			}
			if(wtsum<0.01){
				error("PSF is not available for (%g,%g). wtsum=%g\n", thetax, thetay, wtsum);
			}
			/*Now do the actual reading */
			for(int ix=0; ix<2; ix++){
				for(int iy=0; iy<2; iy++){
					real wtxi=fabs(((1-ix)-wtx)*((1-iy)-wty))/wtsum;
					if(fnztilt[iy][ix]){
						file_t* fp_ztilt=zfopen(fnztilt[iy][ix], "rb");
						header_t header={0,0,0,0};
						read_header(&header, fp_ztilt);

						if(iscell(&header.id)){
							// error("expected data type: %u, got %u\n",(uint32_t)MCC_ANY, header.id);
							nstep=header.nx;
							free(header.str);
							if(stari->nstep==0){
								stari->nstep=nstep;
							} else{
								if(stari->nstep!=nstep){
									error("Different type has different steps\n");
								}
							}
							if(!P(stari->ztiltout,ipowfs)){
								P(stari->ztiltout,ipowfs)=dnew(ng, nstep);
							}
							dmat* ztiltout=P(stari->ztiltout,ipowfs);
							for(long istep=0; istep<nstep; istep++){
								dmat* ztilti=dreaddata(fp_ztilt, 0);
								for(int ig=0; ig<ng; ig++){
									P(ztiltout,ig, istep)+=P(ztilti,ig)*wtxi;
								}
								dfree(ztilti);
							}
						} else{
							dmat* tmp=dreaddata(fp_ztilt, &header);
							dadd(&P(stari->ztiltout,ipowfs), 1, tmp, wtxi);
							dfree(tmp);
						}
						zfclose(fp_ztilt);
					}/* if(fnwvf) */
					if(fngoff[iy][ix]&&zfexist("%s",fngoff[iy][ix])){
						if(!stari->goff){
							stari->goff=dcellnew(npowfs, 1);
						}
						dmat* tmp=dread("%s", fngoff[iy][ix]);
						dadd(&P(stari->goff,ipowfs), 1, tmp, wtxi);
						dfree(tmp);
					}
				}/*iy */
			}/*ix */
		}/*ipowfs */
	}/*istar */
	if(parms->skyc.verbose){
		toc2("Reading PSF");
	}
	//close(fd);
	return nstep;
}
PNEW(mutex_wvf);
/**
   Read in asterism WFS wvf.*/
long setup_star_read_wvf(star_s* star, int nstar, const parms_s* parms, int seed){
	const real ngsgrid=parms->maos.ngsgrid;
	const int nwvl=parms->maos.nwvl;
	long nstep=0;
	LOCK(mutex_wvf);//allow only 1 thread to read file
	TIC;tic;
	for(int istar=0; istar<nstar; istar++){
		star_s* stari=&star[istar];
		int npowfs=parms->maos.npowfs;
		stari->wvfout=mycalloc(npowfs, ccell**);
		const real thetax=stari->thetax*RAD2AS;/*in as */
		const real thetay=stari->thetay*RAD2AS;

		real thxnorm=thetax/ngsgrid;
		real thynorm=thetay/ngsgrid;
		long thxl=(long)floor(thxnorm);/*Used to be real, but -0 appears. */
		long thyl=(long)floor(thynorm);
		real wtx=thxnorm-thxl;
		real wty=thynorm-thyl;
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			const int msa=parms->maos.msa[ipowfs];
			const int nsa=parms->maos.nsa[ipowfs];
			if(stari->use[ipowfs]==0){
				continue;
			}
			char* fnwvf[2][2]={{NULL,NULL},{NULL,NULL}};
			pistat_s* pistati=&stari->pistat[ipowfs];

			/*info("Reading PSF for (%5.1f, %5.1f), ipowfs=%d\n",thetax,thetay,ipowfs); */
			real wtsum=0;
			for(int ix=0; ix<2; ix++){
				real thx=(thxl+ix)*ngsgrid;
				for(int iy=0; iy<2; iy++){
					real thy=(thyl+iy)*ngsgrid;
					real wtxi=fabs(((1-ix)-wtx)*((1-iy)-wty));

					if(wtxi<0.01){
					/*dbg("skipping ix=%d,iy=%d because wt=%g\n",ix,iy,wtxi); */
						continue;
					}
					fnwvf[iy][ix]=myalloca(PATH_MAX, char);
					snprintf(fnwvf[iy][ix], PATH_MAX, "%s/wvfout/wvfout_seed%d_sa%d_x%g_y%g",
						dirstart, seed, msa, thx, thy);

					if(!zfexist("%s",fnwvf[iy][ix])){
					//warning("%s doesnot exist\n",fnwvf[iy][ix]);
						fnwvf[iy][ix]=0;
					} else{
						wtsum+=wtxi;
					}
				}
			}
			if(wtsum<0.01){
				error("PSF is not available for (%g,%g). wtsum=%g\n", thetax, thetay, wtsum);
			}
			/*Now do the actual reading */
			for(int ix=0; ix<2; ix++){
				for(int iy=0; iy<2; iy++){
					real wtxi=fabs(((1-ix)-wtx)*((1-iy)-wty))/wtsum;
					if(fnwvf[iy][ix]){
					/*dbg("Loading %.4f x %s\n", wtxi, fnwvf[iy][ix]); */
						file_t* fp_wvf=zfopen(fnwvf[iy][ix], "rb");
						header_t header={0,0,0,0};
						read_header(&header, fp_wvf);
						if(!iscell(&header.id)){
							error("expected data type: %u, got %u\n", (uint32_t)MCC_ANY, header.id);
						}
						nstep=header.nx;
						free(header.str);
						if(parms->skyc.limitnstep>0&&nstep>parms->skyc.limitnstep){
							nstep=parms->skyc.limitnstep;
							warning("Only read %ld steps\n", nstep);
						}
						if(stari->nstep==0){
							stari->nstep=nstep;
						} else{
							if(stari->nstep!=nstep){
								error("Different type has different steps\n");
							}
						}

						if(!stari->wvfout[ipowfs]){
							stari->wvfout[ipowfs]=mycalloc(nstep, ccell*);
						}
						ccell** pwvfout=stari->wvfout[ipowfs];
						for(long istep=0; istep<nstep; istep++){
							ccell* wvfi=ccellreaddata(fp_wvf, 0);
							ccelladd(&(pwvfout[istep]), 1, wvfi, wtxi);
							ccellfree(wvfi);
						}
						/*zfeof(fp_wvf); */
						zfclose(fp_wvf);
					}
				}/*iy */
			}/*ix */
			/*Don't bother to scale ztiltout since it does not participate in physical optics simulations. */
			if(parms->skyc.bspstrehl){
				dmat* scale=pistati->scale;
				ccell** pwvfout=stari->wvfout[ipowfs];
				for(int iwvl=0; iwvl<nwvl; iwvl++){
					for(int isa=0; isa<nsa; isa++){
					/*dbg("Scaling WVF isa %d iwvl %d with %g\n", isa, iwvl, P(scale,isa,iwvl)); */
						for(long istep=0; istep<stari->nstep; istep++){
							if(pwvfout[istep]){
								cscale(P(pwvfout[istep],isa,iwvl), P(scale, isa, iwvl));
							}
						}/*istep */
					}/*isa */
				}/*iwvl */
			}/* */
		}/*ipowfs */
	}/*istar */
	if(parms->skyc.verbose){
		toc2("Reading PSF");
	}
	UNLOCK(mutex_wvf);
	//close(fd);
	return nstep;
}
/*static int sortfun_snr(const star_s *p1, const star_s *p2){
	int ix=NX(p1->pistat[0].snr)-1;
	real s1=P(p1->pistat[0].snr, ix);
	real s2=P(p2->pistat[0].snr, ix);
	return s1<s2?1:-1; //-1: keep order. 1: reverse order
}*/
static void print_stars(const star_s *star, int nstar, const dmat *dtrats){
	if(!nstar) return;
	int npowfs=NX(star[0].minidtrat);
	int nwvl=NX(star[0].mags);
	for(int istar=0; istar<nstar; istar++){
		info("star %d at (%3.0f %3.0f)\" mag=(", istar, star[istar].thetax*RAD2AS, star[istar].thetay*RAD2AS);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			info("%4.1f%s", P(star[istar].mags,iwvl),iwvl+1<nwvl?" ":"), ");
		}
		for(int ipowfs=npowfs-1; ipowfs>=0; ipowfs--){
			int idtrat=(int)P(star[istar].minidtrat,ipowfs);
			if(idtrat>=0){
				pistat_s* pistat=&star[istar].pistat[ipowfs];
				int nsa=NX(pistat->i0);
				info("%s: snr=%4.1f at dtrat=%2.0f%s", nsa==1?"TT":"TTF", P(pistat->snr,idtrat), P(dtrats,idtrat),ipowfs==0?"":", ");
			}
		}
		info("\n");
	}
}
/**
   setup "star" data array from star information and read in average pixel
   intensitys.  Check for star PSF size. If it is larger than 5, we don't use
   the star because the PSF is too broad.
*/
star_s* setup_star(int* nstarout, sim_s* simu, dmat* stars, int seed){
	const parms_s* parms=simu->parms;
	powfs_s* powfs=simu->powfs;
	if(!stars){
		return NULL;
	}
	star_s* star=setup_star_create(parms, stars);
	int nstar=stars->ny;
	const int npowfs=parms->maos.npowfs;
	setup_star_read_pistat(simu, star, nstar, seed);
	setup_star_siglev(parms, star, nstar);
	//setup_star_gnea(parms, star, nstar, simu->neaspec_dtrats);
	setup_star_mtch(parms, powfs, star, nstar, simu->nonlin);
	setup_star_gm(parms, powfs, star, nstar);
	int jstar=0;
	for(int istar=0; istar<nstar; istar++){
		int skip=0;
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			if(star[istar].use[ipowfs]==-1){
				skip++;
			}
		}
		if(skip==npowfs||(parms->skyc.maxstar>0 && jstar>=parms->skyc.maxstar)){//remove the star;
			free_istar(star+istar, parms);
		} else{
			if(jstar!=istar){
				memcpy(&star[jstar], &star[istar], sizeof(star_s));
			}
			jstar++;
		}
	}
	if(jstar!=nstar){
		nstar=jstar;
		star=myrealloc(star, jstar, star_s);
	}
	//sort stars to have descending snr order.
	//qsort(star, jstar, sizeof(star_s), (int(*)(const void *, const void *))sortfun_snr);
	*nstarout=nstar;
	if(parms->skyc.verbose){
		print_stars(star, nstar, parms->skyc.dtrats);
	}
	return star;
}
void free_istar(star_s* star, const parms_s* parms){
	const int npowfs=parms->maos.npowfs;
	free_pistat(star->pistat, npowfs, parms);
	if(star->wvfout){
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			if(star->wvfout[ipowfs]){
				for(int istep=0; istep<star->nstep; istep++){
					ccellfree(star->wvfout[ipowfs][istep]);
				}
			}
			free(star->wvfout[ipowfs]);
		}
		free(star->wvfout);
	}
	dcellfree(star->ztiltout);
	dcellfree(star->goff);
	dcellfree(star->g);
	dfree(star->mags);
	free(star->use);
	dcellfree(star->siglev);
	dfree(star->siglevtot);
	dfree(star->bkgrnd);
	dfree(star->minidtrat);
}
/**
   Free array of star_s.
 */
void free_star(star_s* stars, int nstar, const parms_s* parms){
	for(int istar=0; istar<nstar; istar++){
		free_istar(stars+istar, parms);
	}
	free(stars);
}
/**
   Free pixel intensities.
 */
void free_pistat(pistat_s* pistat, int npistat, const parms_s* parms){
	if(!pistat) return;
	for(int ipistat=0; ipistat<npistat; ipistat++){
		dcellfree(pistat[ipistat].psf);
		dcellfree(pistat[ipistat].i0);
		dcellfree(pistat[ipistat].gx);
		dcellfree(pistat[ipistat].gy);
		dcellfree(pistat[ipistat].i0s);
		dcellfree(pistat[ipistat].gxs);
		dcellfree(pistat[ipistat].gys);
		dcellfree(pistat[ipistat].sanea);
		dcellfree(pistat[ipistat].sanea0);
		dfree(pistat[ipistat].scale);
		dcellfree(pistat[ipistat].neaspec);
		dfree(pistat[ipistat].snr);
		int ndtrat=parms->skyc.ndtrat;
		if(pistat[ipistat].mtche){
			for(int idtrat=0; idtrat<ndtrat; idtrat++){
				dcellfree(pistat[ipistat].mtche[idtrat]);
			}
			free(pistat[ipistat].mtche);
			pistat[ipistat].mtche=NULL;
		}

	}
	free(pistat);
}
