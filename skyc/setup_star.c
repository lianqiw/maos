/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
static STAR_S *setup_star_create(const PARMS_S *parms, dmat *coord){
    if(!coord){
	return NULL;/*there are no stars available. */
    }
    int nstar=coord->ny;
    PDMAT(coord,pc);
    int nwvl=parms->maos.nwvl;
    STAR_S *star=calloc(nstar, sizeof(STAR_S));
    double ngsgrid=parms->maos.ngsgrid/206265;
    double r2=pow(parms->skyc.patfov/206265./2.,2);
    double keepout=pow(parms->skyc.keepout/206265,2);
    int jstar=0;
    assert(nwvl+2==coord->nx);
    for(int istar=0; istar<nstar; istar++){
	if(parms->skyc.ngsalign){
	    star[jstar].thetax=round(pc[istar][0]/ngsgrid)*ngsgrid;
	    star[jstar].thetay=round(pc[istar][1]/ngsgrid)*ngsgrid;
	    if(pow(star[jstar].thetax,2)+pow(star[jstar].thetay,2)>r2){
		star[jstar].thetax=trunc(pc[istar][0]/ngsgrid)*ngsgrid;
		star[jstar].thetay=round(pc[istar][1]/ngsgrid)*ngsgrid;
		if(pow(star[jstar].thetax,2)+pow(star[jstar].thetay,2)>r2){
		    star[jstar].thetax=round(pc[istar][0]/ngsgrid)*ngsgrid;
		    star[jstar].thetay=trunc(pc[istar][1]/ngsgrid)*ngsgrid;
		    if(pow(star[jstar].thetax,2)+pow(star[jstar].thetay,2)>r2){
			star[jstar].thetax=trunc(pc[istar][0]/ngsgrid)*ngsgrid;
			star[jstar].thetay=trunc(pc[istar][1]/ngsgrid)*ngsgrid;
			if(pow(star[jstar].thetax,2)+pow(star[jstar].thetay,2)>r2){
			    error("What?\n");
			}
		    }
		}
	    }
	}else{
	    star[jstar].thetax=pc[istar][0];
	    star[jstar].thetay=pc[istar][1];
	}
	for(int kstar=0; kstar<jstar; kstar++){
	    if(pow(star[jstar].thetax-star[kstar].thetax,2)
	       +pow(star[jstar].thetay-star[kstar].thetay,2)<keepout){
		/*warning("start %d is too close to %d. use J brightest.\n", jstar, kstar); */
		if(pc[istar][0]<star[kstar].mags->p[0]){
		    memcpy(star[kstar].mags->p, pc[istar]+2, sizeof(double)*nwvl);
		    star[kstar].thetax=star[jstar].thetax;
		    star[kstar].thetay=star[jstar].thetay;
		}
		continue;
	    }
	}
	star[jstar].mags=dnew(nwvl,1);
	memcpy(star[jstar].mags->p, pc[istar]+2, sizeof(double)*nwvl);
	star[jstar].use=calloc(parms->maos.npowfs, sizeof(int));
	jstar++;
    }
    if(jstar<nstar){
	/*warning2("%d stars dropped\n", nstar-jstar); */
	coord->ny=jstar;
	star=realloc(star, jstar*sizeof(STAR_S));
    }
    return star;
}

/**
   Read in pistat information, used to compute matched filter, and SANEA.
*/
static void setup_star_read_pistat(SIM_S *simu, STAR_S *star, int nstar, int seed){
    const PARMS_S *parms=simu->parms;
    const int npowfs=parms->maos.npowfs;
    const int nwvl=parms->maos.nwvl;
    const double ngsgrid=parms->maos.ngsgrid;
    for(int istar=0; istar<nstar; istar++){
	STAR_S *stari=&star[istar];
	stari->pistat=calloc(npowfs, sizeof(PISTAT_S));
	const double thetax=stari->thetax*206265;/*in as */
	const double thetay=stari->thetay*206265;
	

	double thxnorm=thetax/ngsgrid;
	double thynorm=thetay/ngsgrid;
	long thxl=(long)floor(thxnorm);
	long thyl=(long)floor(thynorm);
	double wtx=thxnorm-thxl;
	double wty=thynorm-thyl;
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	    const int msa=parms->maos.msa[ipowfs];
	    const int nsa=parms->maos.nsa[ipowfs];
	    dcell *avgpsf=NULL;
	    dmat *grad=NULL;
	    double wtsum=0;
	    for(int ix=0; ix<2; ix++){
		double thx=(thxl+ix)*ngsgrid;
		for(int iy=0; iy<2; iy++){
		    double thy=(thyl+iy)*ngsgrid;
		    double wtxi=fabs(((1-ix)-wtx)*((1-iy)-wty));
		    if(wtxi<0.01){
			/*info("skipping ix=%d,iy=%d because wt=%g\n",ix,iy,wtxi); */
			continue;
		    }
		    char fn[PATH_MAX];
		    snprintf(fn,PATH_MAX,"%s/pistat/pistat_seed%d_sa%d_x%g_y%g",
			     dirstart, seed,msa,thx,thy);
		    if(!zfexist(fn)){
			/*warning("%s doesn't exist\n",fn); */
		    }else{
			dcell *avgpsfi=dcellread("%s",fn);
			dcelladd(&avgpsf, 1, avgpsfi, wtxi);
			dcellfree(avgpsfi);
			wtsum+=wtxi;
			
			snprintf(fn,PATH_MAX,"%s/pistat/gstat_seed%d_sa%d_x%g_y%g",
				 dirstart, seed, msa, thx, thy);
			dmat *gradi=dread("%s",fn);
			dadd(&grad, 1, gradi, wtxi);
			dfree(gradi);
		    }
		}
	    }
	    if(wtsum<0.01){
		warning("PISTAT is not available for (%g,%g) msa=%d\n",thetax,thetay,msa);
	    }
	    dscale(grad, 1./wtsum);
	    dcellscale(avgpsf, 1./wtsum);
	    dmat *scale=NULL;
	    if(parms->skyc.bspstrehl){
		scale=dnew(nsa,nwvl);
		dmat *gx=dnew(1,1); gx->p[0]=thxnorm;
		dmat *gy=dnew(1,1); gy->p[0]=thynorm;
		if(nsa!=avgpsf->nx || nwvl!=avgpsf->ny){
		    error("Mismatch: nsa=%d, nwvl=%d, avgpsf->nx=%ld, avgpsf->ny=%ld\n",
			  nsa, nwvl, avgpsf->nx, avgpsf->ny);
		}
		for(int ic=0; ic<nsa*nwvl; ic++){
		    dmat *val=dbspline_eval(simu->bspstrehl[ipowfs][ic],
					    simu->bspstrehlxy,simu->bspstrehlxy,
					    gx, gy);
		    double ratio=val->p[0]/avgpsf->p[ic]->p[0];
		    /*info("strehl: bilinear: %g, cubic: %g\n", avgpsf->p[ic]->p[0],val->p[0]); */
		    if(ratio<0){
			warning("Ratio is less than zero.\n");
			scale->p[ic]=1;
		    }else{
			dscale(avgpsf->p[ic], ratio);
			scale->p[ic]=ratio;
			grad->p[ic]*=ratio;
		    }
		    dfree(val);
		}
		dfree(gx);
		dfree(gy);
	    }

	    stari->pistat[ipowfs].psf=avgpsf;/*PSF is in corner. */
	    stari->pistat[ipowfs].grad=grad;
	    stari->pistat[ipowfs].scale=scale;
	    if(parms->skyc.dbg){
		dcellwrite(avgpsf, "%s/avgpsf_star%d_ipowfs%d_psf",dirsetup,istar,ipowfs);
		dwrite(grad, "%s/pistat_star%d_ipowfs%d_grad",dirsetup,istar,ipowfs);
	    }
	}
    }
}

/**
   Compute Signal level
*/
static void setup_star_siglev(const PARMS_S *parms, STAR_S *star, int nstar){
    const long npowfs=parms->maos.npowfs;
    const long nwvl=parms->maos.nwvl;
    const double r2=pow(parms->skyc.patfov/206265./2.,2);
    PDMAT(parms->skyc.rnefs,rnefs);
    for(int istar=0; istar<nstar; istar++){
	star[istar].siglev=dnew(nwvl,npowfs);
	star[istar].bkgrnd=dnew(npowfs,1);
	star[istar].siglevtot=dnew(npowfs,1);
	/*Normalized angular distance */
	double th2=(pow(star[istar].thetax,2)+pow(star[istar].thetay,2))/r2;
	/*Field dependent error: nm^2=nma^2+nmb^2*theta_norm^2; */
	double imperrnm=sqrt(pow(parms->skyc.imperrnm,2)+th2*pow(parms->skyc.imperrnmb,2));
	for(long ipowfs=0; ipowfs<npowfs; ipowfs++){
	    int iscircle=parms->maos.nsa[ipowfs]<=4?1:0;
	    photon_flux(&star[istar].siglev->p[nwvl*ipowfs], 
			&star[istar].siglevtot->p[ipowfs],
			&star[istar].bkgrnd->p[ipowfs],
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
			rnefs[ipowfs][parms->skyc.ndtrat-1]);
	    if(parms->skyc.verbose && ipowfs==npowfs-1){
		info2("star %d at (%5.1f %5.1f)",istar, 
		      star[istar].thetax*206265,star[istar].thetay*206265);
		info2(" bkgrnd=%5.2f, pixtheta=%4.1fmas mag=[",
		      star[istar].bkgrnd->p[ipowfs],parms->skyc.pixtheta[ipowfs]*206265000);
		for(int iwvl=0; iwvl<parms->skyc.nwvl; iwvl++){
		    info2("%5.2f ", star[istar].mags->p[iwvl]);
		}
		info2("] siglev=[");
		for(int iwvl=0; iwvl<parms->skyc.nwvl; iwvl++){
		    info2("%6.1f ", star[istar].siglev->p[iwvl+nwvl*ipowfs]);
		}
		info2("]\n");
	    }
	}
    }
}
/**
   Compute an additional gradient measurement noise that is caused by gradients
   in ideal NGS mode corrected PSFs. This is to be treated as an additional
   measurement error.  */
static void setup_star_gnea(const PARMS_S *parms, STAR_S *star, int nstar){
    const long npowfs=parms->maos.npowfs;
    const int phystart=parms->skyc.phystart;
    for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	const long nsa=parms->maos.nsa[ipowfs];
	for(int istar=0; istar<nstar; istar++){
	    PISTAT_S *pistat=&star[istar].pistat[ipowfs];
	    int ndtrat=parms->skyc.ndtrat;
	    pistat->gnea=dcellnew(ndtrat, 1);
	    PDMAT(pistat->grad, pgrad);
	    for(int idtrat=0; idtrat<ndtrat; idtrat++){
		int dtrat=parms->skyc.dtrats[idtrat];
		dmat *mean_grad=NULL, *mean_gradsq=NULL;
		int count=0;
		dmat *grad=dnew(nsa*2,1);
		for(int istep=0; istep<pistat->grad->ny; istep++){
		    if(istep>phystart){/*make sure same alignment. */
			if(istep % dtrat == 0){
			    dzero(grad);
			}	
			for(int isa=0; isa<nsa*2; isa++){
			    grad->p[isa]+=pgrad[istep][isa];
			}
			if((istep+1) % dtrat == 0){/*has output */
			    dscale(grad, 1./dtrat);
			    count++;
			    dadd(&mean_grad, 1, grad, 1);
			    dcwpow(grad,2);
			    dadd(&mean_gradsq, 1, grad, 1);
			}
		    }
		}
		dscale(mean_grad, 1./count);
		dscale(mean_gradsq, 1./count);
		dcwpow(mean_grad, 2);
		/*variance=<gg>-<g><g> */
		dadd(&pistat->gnea->p[idtrat], 1, mean_gradsq, 1);
		dadd(&pistat->gnea->p[idtrat], 1, mean_grad, -1);
		dfree(grad);
		dfree(mean_grad);
		dfree(mean_gradsq);
	    }
	    if(parms->skyc.dbg){
		dcellwrite(pistat->gnea,"%s/star%d_ipowfs%d_gnea", dirsetup, istar, ipowfs);
	    }
	}
    }
}
/**
   Setup matched filter for stars.
 */
static void setup_star_mtch(const PARMS_S *parms, POWFS_S *powfs, STAR_S *star, int nstar){
    const long nwvl=parms->maos.nwvl;
    const long npowfs=parms->maos.npowfs;
    PDMAT(parms->skyc.rnefs,rnefs);
    for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	const long nsa=parms->maos.nsa[ipowfs];
	const long pixpsa=parms->skyc.pixpsa[ipowfs];
	for(int istar=0; istar<nstar; istar++){
	    PISTAT_S *pistat=&star[istar].pistat[ipowfs];
	    pistat->i0=dcellnew(nsa,nwvl);
	    pistat->gx=dcellnew(nsa,nwvl);
	    pistat->gy=dcellnew(nsa,nwvl);
	    
	    pistat->i0s=dcellnew(nsa,1);
	    pistat->gxs=dcellnew(nsa,1);
	    pistat->gys=dcellnew(nsa,1);

	    PDCELL(pistat->psf, psf);
	    PDCELL(pistat->i0, i0);
	    PDCELL(pistat->gx, gx);
	    PDCELL(pistat->gy, gy);

	    for(long iwvl=0; iwvl<nwvl; iwvl++){
		for(long isa=0; isa<nsa; isa++){
		    double siglev=star[istar].siglev->p[iwvl+nwvl*ipowfs];
		    i0[iwvl][isa]=dnew(pixpsa,pixpsa);
		    gx[iwvl][isa]=dnew(pixpsa,pixpsa);
		    gy[iwvl][isa]=dnew(pixpsa,pixpsa);
		    psf2i0gxgy(i0[iwvl][isa],gx[iwvl][isa],gy[iwvl][isa],
			       psf[iwvl][isa],powfs[ipowfs].dtf+iwvl);
		    dadd(&pistat->i0s->p[isa], 1, i0[iwvl][isa], siglev);
		    dadd(&pistat->gxs->p[isa], 1, gx[iwvl][isa], siglev);
		    dadd(&pistat->gys->p[isa], 1, gy[iwvl][isa], siglev);
		}
		 
	    }

	    const double pixtheta=parms->skyc.pixtheta[ipowfs];
	    int ndtrat=parms->skyc.ndtrat;
	    pistat->mtche=calloc(ndtrat, sizeof(dcell*));
	    pistat->sanea=dcellnew(ndtrat,1);
	    dcell *i0s=NULL; dcell *gxs=NULL; dcell *gys=NULL;

	    for(int idtrat=0; idtrat<ndtrat; idtrat++){
		int dtrat=parms->skyc.dtrats[idtrat];
		dcelladd(&i0s, 0, pistat->i0s, dtrat);
		dcelladd(&gxs, 0, pistat->gxs, dtrat);
		dcelladd(&gys, 0, pistat->gys, dtrat);
		mtch(&pistat->mtche[idtrat], &pistat->sanea->p[idtrat],
		     i0s, gxs, gys, pixtheta, rnefs[ipowfs][idtrat], 
		     star[istar].bkgrnd->p[ipowfs]*dtrat, parms->skyc.mtchcr);
		if(parms->skyc.dbg){
		    dcellwrite(pistat->mtche[idtrat], "%s/star%d_ipowfs%d_mtche_dtrat%d",
			       dirsetup,istar,ipowfs,dtrat);
		}
	    }
	    if(parms->skyc.dbg){
		dcellwrite(pistat->sanea, "%s/star%d_ipowfs%d_sanea",
			   dirsetup,istar,ipowfs);
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
static void setup_star_g(const PARMS_S *parms, POWFS_S *powfs, STAR_S *star, int nstar){
    const long npowfs=parms->maos.npowfs;
    
    const double hc=parms->maos.hc;
    const double hs=parms->maos.hs;
    const double scale=pow(1.-hc/hs, -2);
    const double scale1=1.-scale;

    for(int istar=0; istar<nstar; istar++){
	star[istar].g=dcellnew(npowfs, 1);
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	    const long nsa=parms->maos.nsa[ipowfs];
	    const double thetax=star[istar].thetax;
	    const double thetay=star[istar].thetay;
	    star[istar].g->p[ipowfs]=dnew(nsa*2, parms->maos.nmod);
	    PDMAT(star[istar].g->p[ipowfs], pg);
	    for(long isa=0; isa<nsa; isa++){
		const double xm=powfs[ipowfs].locxamp[isa];/*dot of x with amp. */
		const double ym=powfs[ipowfs].locyamp[isa];

		pg[0][isa]     = 1.;
		pg[1][isa+nsa] = 1.;
		pg[2][isa]     = (scale1*2*xm - 2*thetax*hc*scale);
		pg[2][isa+nsa] = (scale1*2*ym - 2*thetay*hc*scale);
		pg[3][isa]     = (scale1*2*xm - 2*thetax*hc*scale);
		pg[3][isa+nsa] = (-scale1*2*ym+ 2*thetay*hc*scale);
		pg[4][isa]     = (scale1*ym   - thetay*hc*scale);
		pg[4][isa+nsa] = (scale1*xm   - thetax*hc*scale);
	    }
	}
	if(parms->skyc.dbg){
	    dcellwrite(star[istar].g,"%s/star%d_g",dirsetup,istar);
	}
    }
}
/**
   Read in asterism WFS wvf and ztilt.*/
long setup_star_read_wvf(STAR_S *star, int nstar, const PARMS_S *parms, int seed){
    const double ngsgrid=parms->maos.ngsgrid;
    const int nwvl=parms->maos.nwvl;
    long nstep=0;
    TIC;tic;
    char fnlock[PATH_MAX];
    snprintf(fnlock, PATH_MAX, "%s/wvfout/wvfout.lock", dirstart);
    /*Obtain exclusive lock before proceeding so that no two process will read concurrently. */
    int fd=lock_file(fnlock, 1, 0);
    for(int istar=0; istar<nstar; istar++){
	STAR_S *stari=&star[istar];
	int npowfs=parms->maos.npowfs;
	stari->wvfout=calloc(npowfs, sizeof(ccell**));
	stari->ztiltout=calloc(npowfs, sizeof(dcell*));
	const double thetax=stari->thetax*206265;/*in as */
	const double thetay=stari->thetay*206265;

	double thxnorm=thetax/ngsgrid;
	double thynorm=thetay/ngsgrid;
	long thxl=(long)floor(thxnorm);/*Used to be double, but -0 appears. */
	long thyl=(long)floor(thynorm);
	double wtx=thxnorm-thxl;
	double wty=thynorm-thyl;
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	    const int msa=parms->maos.msa[ipowfs];
	    const int nsa=parms->maos.nsa[ipowfs];
	    if(stari->use[ipowfs]==0){
		continue;
	    }
	    char   *fnwvf[2][2]={{NULL,NULL},{NULL,NULL}};
	    char   *fnztilt[2][2]={{NULL,NULL},{NULL,NULL}};

	    PISTAT_S *pistati=&stari->pistat[ipowfs];
	    
	    /*info2("Reading PSF for (%5.1f, %5.1f), ipowfs=%d\n",thetax,thetay,ipowfs); */
	    double wtsum=0;
	    for(int ix=0; ix<2; ix++){
		double thx=(thxl+ix)*ngsgrid;
		for(int iy=0; iy<2; iy++){
		    double thy=(thyl+iy)*ngsgrid;
		    double wtxi=fabs(((1-ix)-wtx)*((1-iy)-wty));

		    if(wtxi<0.01){
			/*info("skipping ix=%d,iy=%d because wt=%g\n",ix,iy,wtxi); */
			continue;
		    }
		    fnwvf[iy][ix]=alloca(PATH_MAX*sizeof(char));
		    snprintf(fnwvf[iy][ix],PATH_MAX,"%s/wvfout/wvfout_seed%d_sa%d_x%g_y%g",
			     dirstart,seed,msa,thx,thy);
		    fnztilt[iy][ix]=alloca(PATH_MAX*sizeof(char));
		    snprintf(fnztilt[iy][ix],PATH_MAX,"%s/ztiltout/ztiltout_seed%d_sa%d_x%g_y%g",
			     dirstart,seed,msa,thx,thy);
		    if(!zfexist(fnwvf[iy][ix])){
			warning("%s doesnot exist\n",fnwvf[iy][ix]);
			fnwvf[iy][ix]=NULL;
			fnztilt[iy][ix]=NULL;
		    }else{
			if(!zfexist(fnztilt[iy][ix])){
			    error("%s exist, but %s does not\n", fnwvf[iy][ix],fnztilt[iy][ix]);
			}
			wtsum+=wtxi;
		    }
		}
	    }
	    if(wtsum<0.01){
		error("PSF is not available for (%g,%g). wtsum=%g\n",thetax,thetay, wtsum);
	    }
	    /*Now do the actual reading */
	    for(int ix=0; ix<2; ix++){
		for(int iy=0; iy<2; iy++){
		    double wtxi=fabs(((1-ix)-wtx)*((1-iy)-wty))/wtsum;
		    if(fnwvf[iy][ix]){
			/*info("Loading %.4f x %s\n", wtxi, fnwvf[iy][ix]); */
			file_t *fp_wvf=zfopen(fnwvf[iy][ix],"rb");
			header_t header;
			read_header(&header, fp_wvf);
			if(!iscell(header.magic)){
			    error("expected data type: %u, got %u\n",(uint32_t)MCC_ANY, header.magic);
			}
			nstep=header.nx;
			free(header.str);
			if(parms->skyc.limitnstep >0 && nstep>parms->skyc.limitnstep){
			    nstep=parms->skyc.limitnstep;
			    warning("Only read %ld steps\n",nstep);
			}
			if(stari->nstep==0){
			    stari->nstep=nstep;
			}else{
			    if(stari->nstep!=nstep){
				error("Different WFS has different steps\n");
			    }
			}
		    
			if(!stari->wvfout[ipowfs]){
			    stari->wvfout[ipowfs]=calloc(nstep,sizeof(ccell*));
			}
			ccell **pwvfout=stari->wvfout[ipowfs];
			for(long istep=0; istep<nstep; istep++){
			    ccell *wvfi=ccellreaddata(fp_wvf, 0);
			    ccelladd(&(pwvfout[istep]), 1, wvfi, wtxi);
			    ccellfree(wvfi);
			}
			/*zfeof(fp_wvf); */
			zfclose(fp_wvf);
		    }
		    if(fnztilt[iy][ix]){
			file_t *fp_ztilt=zfopen(fnztilt[iy][ix],"rb");
			header_t header;
			read_header(&header, fp_ztilt);
			if(!iscell(header.magic)){
			    error("expected data type: %u, got %u\n",(uint32_t)MCC_ANY, header.magic);
			}
			nstep=header.nx;
			free(header.str);
			if(parms->skyc.limitnstep >0 && nstep>parms->skyc.limitnstep){
			    nstep=parms->skyc.limitnstep;
			    warning("Only read %ld steps\n",nstep);
			}
			if(!stari->ztiltout[ipowfs]){
			    stari->ztiltout[ipowfs]=dcellnew(nstep,1);
			}
			dmat  **pztiltout=stari->ztiltout[ipowfs]->p;
			for(long istep=0; istep<nstep; istep++){
			    dmat *ztilti=dreaddata(fp_ztilt, 0);
			    dadd(&(pztiltout[istep]), 1, ztilti, wtxi);/*(2nsa)*nstep dmat array */
			    dfree(ztilti);
			}
			zfclose(fp_ztilt);
		    }/* if(fnwvf) */
		}/*iy */
	    }/*ix */
	    /*Don't bother to scale ztiltout since it does not participate in physical optics simulations. */
	    if(parms->skyc.bspstrehl){
		PDMAT(pistati->scale, scale);
		ccell **pwvfout=stari->wvfout[ipowfs];
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    for(int isa=0; isa<nsa; isa++){
			/*info("Scaling WVF isa %d iwvl %d with %g\n", isa, iwvl, scale[iwvl][isa]); */
			for(long istep=0; istep<stari->nstep; istep++){
			    cscale(pwvfout[istep]->p[isa+nsa*iwvl], scale[iwvl][isa]);
			}/*istep */
		    }/*isa */
		}/*iwvl */
	    }/* */
	}/*ipowfs */
    }/*istar */
    if(parms->skyc.verbose){
	toc2("Reading PSF");
    }
    close(fd);
    return nstep;
}
/**
   setup "star" data array from star information and read in average pixel
   intensitys.  Check for star PSF size. If it is larger than 5, we don't use
   the star because the PSF is too broad.
*/
 
STAR_S *setup_star(int *nstarout, SIM_S *simu, dmat *stars,int seed){
    const PARMS_S *parms=simu->parms;
    POWFS_S *powfs=simu->powfs;
    if(!stars){
	return NULL;
    }
    STAR_S *star=setup_star_create(parms, stars);
    int nstar=stars->ny;
    setup_star_read_pistat(simu, star, nstar, seed);
    int jstar=0;
    for(int istar=0; istar<nstar; istar++){
	dmat *pistat=star[istar].pistat[parms->skyc.npowfs-1].psf->p[0];
	double pmax=dmax(pistat);
	int size=0;
	for(long i=0; i<pistat->nx*pistat->ny; i++){
	    if(pistat->p[i]>pmax*0.5){
		size++;
	    }
	}
	if(size>6){
	    free_pistat(star[istar].pistat, parms->skyc.npowfs, parms);
	    warning("star %d doesn't have sharpen cores. size=%d\n",istar,size);
	}else{
	    if(istar!=jstar){
		memcpy(&star[jstar], &star[istar], sizeof(STAR_S));
	    }
	    jstar++;
	}
    }
    if(jstar==0){
	free(star);
	return NULL;
    }
    if(jstar>parms->skyc.maxstar){
	jstar=parms->skyc.maxstar;/*we only keed maxstar number of brightest stars. */
    }
    if(jstar!=nstar){
	star=realloc(star,sizeof(STAR_S)*jstar);
	nstar=jstar;
    }
    if(parms->skyc.verbose){
	info2("There are %d stars usable from %d stars\n",jstar,nstar);
    }
    *nstarout=nstar;
    setup_star_siglev(parms, star, nstar);
    setup_star_gnea(parms, star, nstar);
    setup_star_mtch(parms, powfs, star, nstar);
    setup_star_g(parms, powfs, star, nstar);
    return star;
}
/**
   Free array of STAR_S.
 */
void free_star(STAR_S *star, int nstar, const PARMS_S *parms){
    const long npowfs=parms->maos.npowfs;
    for(int istar=0; istar<nstar; istar++){
	free_pistat(star[istar].pistat,npowfs, parms);
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	    if(star[istar].wvfout[ipowfs]){
		for(int istep=0; istep<star[istar].nstep; istep++){
		    ccellfree(star[istar].wvfout[ipowfs][istep]);
		}
	    }
	    free(star[istar].wvfout[ipowfs]);
	    dcellfree(star[istar].ztiltout[ipowfs]);
	}
	free(star[istar].wvfout);
	free(star[istar].ztiltout);
	dcellfree(star[istar].g);
	dfree(star[istar].mags);
	free(star[istar].use);
	dfree(star[istar].siglev);
	dfree(star[istar].siglevtot);
	dfree(star[istar].bkgrnd);
    }
    free(star);
}
/**
   Free pixel intensities.
 */
void free_pistat(PISTAT_S *pistat, int npistat, const PARMS_S *parms){
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
	dfree(pistat[ipistat].scale);
	dfree(pistat[ipistat].grad);
	dcellfree(pistat[ipistat].gnea);
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
