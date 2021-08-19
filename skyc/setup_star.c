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

/**
   \file skyc/setup_star.c
  setup stars.
*/
#include "skyc.h"
#include "parms.h"
#include "types.h"
#include "setup_star.h"
#include "photon.h"
#include "mtch.h"
/**
   Create "star" data array from star information.
*/
static STAR_S* setup_star_create(const PARMS_S* parms, dmat* coord){
	if(!coord){
		return NULL;/*there are no stars available. */
	}
	int nstar=coord->ny;
	dmat* pc=coord;
	int nwvl=parms->maos.nwvl;
	STAR_S* star=mycalloc(nstar, STAR_S);
	real ngsgrid=parms->maos.ngsgrid/206265.;
	real r2=pow(parms->skyc.patfov/206265./2., 2);
	real keepout=pow(parms->skyc.keepout/206265., 2);
	real minrad2=pow(parms->skyc.minrad/206265., 2);
	int jstar=0;
	assert(nwvl+2==coord->nx);
	for(int istar=0; istar<nstar; istar++){
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
			 /*warning("start %d is too close to %d. use J brightest.\n", jstar, kstar); */
				if(P(pc, 0, istar)<P(star[kstar].mags,0)){
					memcpy(star[kstar].mags->p, PCOL(pc, istar)+2, sizeof(real)*nwvl);
					star[kstar].thetax=star[jstar].thetax;
					star[kstar].thetay=star[jstar].thetay;
				}
				continue;
			}
		}
		if(pow(star[istar].thetax, 2)+pow(star[istar].thetay, 2)<minrad2){
			info("Skip star at (%.0f, %.0f) because minrad=%g\n",
				star[istar].thetax*206265, star[istar].thetay*206265, parms->skyc.minrad);
			continue;
		}
		star[jstar].mags=dnew(nwvl, 1);
		memcpy(star[jstar].mags->p, PCOL(pc, istar)+2, sizeof(real)*nwvl);
		star[jstar].use=mycalloc(parms->maos.npowfs, int);
		jstar++;
	}
	if(jstar<nstar){
	/*warning("%d stars dropped\n", nstar-jstar); */
		coord->ny=jstar;
		star=myrealloc(star, jstar, STAR_S);
	}
	return star;
}

/**
   Read in pistat information, used to compute matched filter, and SANEA.
*/
static void setup_star_read_pistat(SIM_S* simu, STAR_S* star, int nstar, int seed){
	const PARMS_S* parms=simu->parms;
	const int npowfs=parms->maos.npowfs;
	const int nwvl=parms->maos.nwvl;
	const real ngsgrid=parms->maos.ngsgrid;
	for(int istar=0; istar<nstar; istar++){
		STAR_S* stari=&star[istar];
		stari->pistat=mycalloc(npowfs, PISTAT_S);
		const real thetax=stari->thetax*206265;/*in as */
		const real thetay=stari->thetay*206265;
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
					real ratio=P(val,0)/P(P(avgpsf,ic),0);
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
			if(parms->skyc.dbg){
				writebin(avgpsf, "%s/avgpsf_star%d_ipowfs%d_psf", dirsetup, istar, ipowfs);
				writebin(stari->pistat[ipowfs].neaspec, "%s/pistat_star%d_ipowfs%d_neaspec", dirsetup, istar, ipowfs);
			}
		}
	}
}

/**
   Compute Signal level
*/
static void setup_star_siglev(const PARMS_S* parms, STAR_S* star, int nstar){
	const long npowfs=parms->maos.npowfs;
	const long nwvl=parms->maos.nwvl;
	const real r2=pow(parms->skyc.patfov/206265./2., 2);
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
			photon_flux(&parms->skyc.zb, P(star[istar].siglev,ipowfs)->p,
				&P(star[istar].siglevtot,ipowfs),
				&P(star[istar].bkgrnd,ipowfs),
				NULL, NULL,
				parms->maos.nwvl,
				parms->maos.wvl,
				star[istar].mags->p,
				parms->maos.dxsa[ipowfs], iscircle,
				parms->skyc.pixtheta[ipowfs],
				parms->maos.dt, parms->maos.za,
				NULL,
				imperrnm,
				parms->skyc.telthruput,
				parms->skyc.qe,
				P(rnefs, parms->skyc.ndtrat-1, ipowfs));
			if(parms->skyc.verbose&&ipowfs==npowfs-1){
				info("star %d at (%5.1f %5.1f)", istar,
					star[istar].thetax*206265, star[istar].thetay*206265);
				info(" bkgrnd=%5.2f, pixtheta=%4.1fmas mag=[",
					P(star[istar].bkgrnd,ipowfs), parms->skyc.pixtheta[ipowfs]*206265000);
				for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
					info("%5.2f ", P(star[istar].mags,iwvl));
				}
				info("] siglev=[");
				for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
					info("%6.1f ", P(P(star[istar].siglev,ipowfs),iwvl));
				}
				info("]\n");
			}
		}
	}
}

/**
   Setup matched filter for stars.
 */
static void setup_star_mtch(const PARMS_S* parms, POWFS_S* powfs, STAR_S* star, int nstar, dcell** nonlin){
	const long nwvl=parms->maos.nwvl;
	const long npowfs=parms->maos.npowfs;
	dmat* rnefs=parms->skyc.rnefs;

	for(int istar=0; istar<nstar; istar++){
		if(!star[istar].minidtrat){
			star[istar].minidtrat=dnew(npowfs, 1);
		}
		dset(star[istar].minidtrat, -1);
		real radius=sqrt(pow(star[istar].thetax, 2)+pow(star[istar].thetay, 2));
		int igg=round(radius*206265/parms->maos.ngsgrid);
		//dbg("radius=%g as, igg=%d\n", radius*206265, igg);
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			const long nsa=parms->maos.nsa[ipowfs];
			const long pixpsa=parms->skyc.pixpsa[ipowfs];
			const real pixtheta=parms->skyc.pixtheta[ipowfs];
			//size of PSF
			const real sigma_theta=parms->skyc.wvlmean/parms->maos.dxsa[ipowfs];
			PISTAT_S* pistat=&star[istar].pistat[ipowfs];
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
						warning_once("Using derivative by FFT\n");
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
			if(parms->skyc.dbg){
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
				/*Add nolinearity*/
				if(nonlin){
					//add linearly not quadratically since the errors are related.
					dmat* nea_nonlin=dinterp1(P(nonlin[ipowfs],igg), NULL, P(pistat->sanea,idtrat), 0);
					for(int i=0; i<nsa*2; i++){
					//info("%g mas", P(P(pistat->sanea,idtrat),i)*206265000);
						P(P(pistat->sanea,idtrat),i)=sqrt(pow(P(P(pistat->sanea,idtrat),i), 2)
							+pow(P(nea_nonlin,i), 2));
		//info("-->%g mas\n", P(P(pistat->sanea,idtrat),i)*206265000);
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
				if(parms->skyc.dbg){
					writebin(pistat->mtche[idtrat], "%s/star%d_ipowfs%d_mtche_dtrat%d",
						dirsetup, istar, ipowfs, dtrat);
				}
				real nea=sqrt(dsumsq(P(pistat->sanea,idtrat))/(nsa*2));
				real snr=sigma_theta/nea;
				P(pistat->snr,idtrat)=snr;
				if(snr>=parms->skyc.snrmin){
					P(star[istar].minidtrat,ipowfs)=idtrat;
				}
			}//for idtrat
			if(P(star[istar].minidtrat,ipowfs)==-1){
				star[istar].use[ipowfs]=-1;
				if(parms->skyc.dbg){
					info("star %2d, powfs %1d: skipped\n", istar, ipowfs);
				}
			} else{
				if(parms->skyc.dbg){
					int idtrat=(int)P(star[istar].minidtrat,ipowfs);
					info("star %2d, powfs %1d: min dtrat=%3d, snr=%4.1f\n", istar, ipowfs,
						(int)P(parms->skyc.dtrats,idtrat), P(pistat->snr,idtrat));
					  //writebin(pistat->sanea, "%s/star%d_ipowfs%d_sanea",
					  //dirsetup,istar,ipowfs);
				}
			}/*idtrat */
			dcellfree(i0s);
			dcellfree(gxs);
			dcellfree(gys);
		}/*for istar */
	}/*for ipowfs */
}
/**
   Compute Modal to gradient operator using average gradients. Similar to Z tilt
   since the mode is low order
 */
static void setup_star_gm(const PARMS_S* parms, POWFS_S* powfs, STAR_S* star, int nstar){
	const long npowfs=parms->maos.npowfs;
	const real hc=parms->maos.hc;
	const real hs=parms->maos.hs;
	const real scale=pow(1.-hc/hs, -2);
	const real scale1=1.-scale;
	const int nmod=parms->maos.nmod;
	for(int istar=0; istar<nstar; istar++){
		star[istar].g=dcellnew(npowfs, 1);
		for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
			const long nsa=parms->maos.nsa[ipowfs];
			const real thetax=star[istar].thetax;
			const real thetay=star[istar].thetay;
			P(star[istar].g,ipowfs)=dnew(nsa*2, nmod);
			dmat* pg=P(star[istar].g,ipowfs);
			for(long isa=0; isa<nsa; isa++){
				const real xm=powfs[ipowfs].locxamp[isa];/*dot of x with amp. */
				const real ym=powfs[ipowfs].locyamp[isa];

				P(pg, isa, 0)=1.;//tip
				P(pg, isa+nsa, 1)=1.;//tilt
				if(parms->maos.indps){
					int indps=parms->maos.indps;
					if(parms->maos.ahstfocus){/*This mode has no global focus*/
						P(pg, isa, indps)=(-2*thetax*hc*scale);
						P(pg, isa+nsa, indps)=(-2*thetay*hc*scale);
					} else{
						P(pg, isa, indps)=(scale1*2*xm-2*thetax*hc*scale);
						P(pg, isa+nsa, indps)=(scale1*2*ym-2*thetay*hc*scale);
					}
					P(pg, isa, indps+1)=(scale1*2*xm-2*thetax*hc*scale);
					P(pg, isa+nsa, indps+1)=(-scale1*2*ym+2*thetay*hc*scale);
					P(pg, isa, indps+2)=(scale1*ym-thetay*hc*scale);
					P(pg, isa+nsa, indps+2)=(scale1*xm-thetax*hc*scale);
				}
				if(parms->maos.indastig){
					const int indastig=parms->maos.indastig;
					P(pg, isa, indastig+1)=(2*xm);//d(x^2-y^2)/dx
					P(pg, isa+nsa, indastig+1)=(-2*ym);
					P(pg, isa, indastig+2)=(ym);
					P(pg, isa+nsa, indastig+2)=(xm);
				}
				if(parms->maos.indfocus){
					P(pg, isa, parms->maos.indfocus)=xm*2;
					P(pg, isa+nsa, parms->maos.indfocus)=ym*2;
				}
			}
		}
		if(parms->skyc.dbg){
			writebin(star[istar].g, "%s/star%d_g", dirsetup, istar);
		}
	}
}
long setup_star_read_ztilt(STAR_S* star, int nstar, const PARMS_S* parms, int seed){
	const real ngsgrid=parms->maos.ngsgrid;
	long nstep=0;
	TIC;tic;
	for(int istar=0; istar<nstar; istar++){
		STAR_S* stari=&star[istar];
		int npowfs=parms->maos.npowfs;
		stari->ztiltout=dcellnew(npowfs, 1);
		const real thetax=stari->thetax*206265;/*in as */
		const real thetay=stari->thetay*206265;

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
						warning_once("Using phygrad\n");
						snprintf(fnztilt[iy][ix], PATH_MAX, "%s/phygrad/phygrad_seed%d_sa%d_x%g_y%g",
							dirstart, seed, msa, thx, thy);
					} else{
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

						if(iscell(&header.magic)){
							// error("expected data type: %u, got %u\n",(uint32_t)MCC_ANY, header.magic);
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

/**
   Read in asterism WFS wvf.*/
long setup_star_read_wvf(STAR_S* star, int nstar, const PARMS_S* parms, int seed){
	const real ngsgrid=parms->maos.ngsgrid;
	const int nwvl=parms->maos.nwvl;
	long nstep=0;
	TIC;tic;
	for(int istar=0; istar<nstar; istar++){
		STAR_S* stari=&star[istar];
		int npowfs=parms->maos.npowfs;
		stari->wvfout=mycalloc(npowfs, ccell**);
		const real thetax=stari->thetax*206265;/*in as */
		const real thetay=stari->thetay*206265;

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
			PISTAT_S* pistati=&stari->pistat[ipowfs];

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
						if(!iscell(&header.magic)){
							error("expected data type: %u, got %u\n", (uint32_t)MCC_ANY, header.magic);
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
							if(pwvfout[istep]->p){
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
	//close(fd);
	return nstep;
}
/**
   setup "star" data array from star information and read in average pixel
   intensitys.  Check for star PSF size. If it is larger than 5, we don't use
   the star because the PSF is too broad.
*/

STAR_S* setup_star(int* nstarout, SIM_S* simu, dmat* stars, int seed){
	const PARMS_S* parms=simu->parms;
	POWFS_S* powfs=simu->powfs;
	if(!stars){
		return NULL;
	}
	STAR_S* star=setup_star_create(parms, stars);
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
		if(skip==npowfs||jstar>=parms->skyc.maxstar){//remove the star;
			free_istar(star+istar, parms);
		} else{
			if(jstar!=istar){
				memcpy(&star[jstar], &star[istar], sizeof(STAR_S));
			}
			jstar++;
		}
	}
	if(jstar!=nstar){
		nstar=jstar;
		star=myrealloc(star, jstar, STAR_S);
	}
	*nstarout=nstar;
	if(parms->skyc.verbose){
		info("There are %d stars usable from %d stars\n", jstar, nstar);
	}
	return star;
}
void free_istar(STAR_S* star, const PARMS_S* parms){
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
   Free array of STAR_S.
 */
void free_star(STAR_S* stars, int nstar, const PARMS_S* parms){
	for(int istar=0; istar<nstar; istar++){
		free_istar(stars+istar, parms);
	}
	free(stars);
}
/**
   Free pixel intensities.
 */
void free_pistat(PISTAT_S* pistat, int npistat, const PARMS_S* parms){
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
