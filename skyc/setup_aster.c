/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "setup_aster.h"
#include "photon.h"
#include "skysim_utils.h"
#include "mtch.h"
#include "utils.h"
#include "setup_star.h"
/**
   \file skyc/setup_aster.c
   Routines to handle asterisms.
 */
/**
   Computes the number of possibilities of selection k items from n items \f$C_n^k\f$.
*/
static long comb_select(long n, long k){
    /*number of possibilities of selecting k from n. */
    if(k>n||k<0) {
	error("Invalid argument\n");
	return -1;
    } else {
	/*compute factorial(n)/factorial(n-k) will overflow if n>=21 */
	/*we use double */
	double res=1;
	long nk=n-k;
	while(n>nk){
	    res*=n--;
	}
	double fk=1;
	while(k>1){
	    fk*=k--;
	}
	return (long)round(res/fk);
    }/* return factorial(n)/(factorial(k)*factorial(n-k)); */
}
/**
   initialize an initial combination composed a vector of non-negative numbers 0,1,2,...
*/
static int* comb_init(long k){
    int *comb=calloc(k, sizeof(int));
    for(int i=0; i<k; i++){
	comb[i]=i;
    }
    return comb;
}
/**
   Find the next combination.
 */
static int comb_next(int *comb, long n, long k){
    if(n<1 || k<1){
	return 0;
    }
    int i = k-1;
    comb[i]++;/*increment to next */
    while(comb[i]+k>= n+i+1 && i>0){/*out of range, increment previous one */
	i--;
	comb[i]++;
    }
    if(comb[0] + k > n){
	return 0;/*no more */
    }
    for(i=i+1; i<k; i++){
	comb[i]=comb[i-1]+1;
    }
    return 1;
}

/**
   Create combination of stars to form asterism. It has the option to put TTF
   always on the brightest for testing purpose.  */
ASTER_S *setup_aster_comb(int *naster, int nstar, const PARMS_S *parms){
    if(nstar==0){
	*naster=0;
	return NULL;
    }else if(parms->skyc.keeporder){
	/*Use the same order as input stars.*/
	ASTER_S *aster=calloc(1, sizeof(ASTER_S));
	*naster=1;
	int npowfs=parms->maos.npowfs;
	int nleft=nstar;
	int stars[npowfs];
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	    stars[ipowfs]=MIN(nleft, parms->skyc.nwfsmax[ipowfs]);
	    nleft-=stars[ipowfs];
	}
	if(nleft>0){
	    warning("skyc.keeporder is set, but there are more stars than needed, dropped the extra\n");
	}
	int ntot=nstar-nleft;
	aster[0].nwfs=ntot;
	aster[0].wfs=calloc(ntot, sizeof(WFS_S));
	aster[0].iaster=0;
	int count=0;
	for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	    for(int istar=0; istar<stars[ipowfs]; istar++){
		aster[0].wfs[count].ipowfs=ipowfs;
		aster[0].wfs[count].istar=count;
		count++;
	    }
	}
	return aster;
    }
    
    int ncomb=1;
    ASTER_S *aster;
    int npowfs=parms->skyc.npowfs;
    int nwfs[npowfs];
    int nleft;
    int nwfstot=0;
    nleft=nstar;
    for(int ipowfs=0; ipowfs<npowfs; ipowfs++){
	if(nleft>=parms->skyc.nwfsmax[ipowfs]){
	    nwfs[ipowfs]=parms->skyc.nwfsmax[ipowfs];
	}else{
	    nwfs[ipowfs]=nleft;
	}
	nwfstot+=nwfs[ipowfs];
	ncomb*=comb_select(nleft,nwfs[ipowfs]);
	nleft-=nwfs[ipowfs];
    }
    if(parms->skyc.ttfbrightest){
	if(parms->maos.msa[0]==2){
	    ncomb/=comb_select(nwfstot,nwfs[0]);
	}else{
	    error("Please revise\n");
	}
    }
    if(parms->skyc.verbose){
	info2("Number of stars: %d, number of asterisms: %d\n", nstar, ncomb);
    }
    aster=calloc(ncomb, sizeof(ASTER_S));
    int count=0;
    int *comb=comb_init(nwfstot);
    do{
	if(npowfs==1){
	    aster[count].nwfs=nwfs[0];
	    aster[count].wfs=calloc(nwfs[0], sizeof(WFS_S));
	    aster[count].iaster=count;
	    for(int iwfs=0; iwfs<nwfs[0]; iwfs++){
		aster[count].wfs[iwfs].istar=comb[iwfs];
		aster[count].wfs[iwfs].ipowfs=0;
	    }
	    count++;
	}else if(npowfs==2){
	    int mask[nwfstot];
	    int *comb2=comb_init(nwfs[0]);
	    do{
		memset(mask, 0, sizeof(int)*nwfstot);
		aster[count].nwfs=nwfstot;
		aster[count].wfs=calloc(nwfstot, sizeof(WFS_S));
		aster[count].iaster=count;
		for(int iwfs=0; iwfs<nwfs[0]; iwfs++){
		    aster[count].wfs[iwfs].ipowfs=0;
		    aster[count].wfs[iwfs].istar=comb[comb2[iwfs]];
		    mask[comb2[iwfs]]=1;
		}
		int jstar=0;
		for(int iwfs=0; iwfs<nwfs[1]; iwfs++){
		    aster[count].wfs[iwfs+nwfs[0]].ipowfs=1;
		    while(mask[jstar]) jstar++;
		    aster[count].wfs[iwfs+nwfs[0]].istar=comb[jstar];
		    mask[jstar]=1;
		}
		count++;
	    }while(comb_next(comb2,nwfstot,nwfs[0]) && !parms->skyc.ttfbrightest);
	    free(comb2);
	}
    }while(comb_next(comb,nstar,nwfstot));
    free(comb);
    if(count!=ncomb){
	warning("ncomb=%d, count=%d. They should equal.\n", ncomb, count);
    }
    *naster=count;
    return aster;
}
/**
   Compute Modal to gradient operator by copying from the stars. Using average
gradients. Similar to Z tilt since the mode is low order */
void setup_aster_g(ASTER_S *aster, STAR_S *star, POWFS_S *powfs, const PARMS_S *parms){
    /*2010-06-08: Tested against MATLAB skycoverage code. */
    aster->g=dcellnew(aster->nwfs,1);
    if(parms->maos.nmod<5 || parms->maos.nmod>6){
	error("Not compatible with the number of NGS modes\n");
    }
    aster->ngs=calloc(aster->nwfs, sizeof(long));
    aster->tsa=0;
    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
	const int istar=aster->wfs[iwfs].istar;
	const int ipowfs=aster->wfs[iwfs].ipowfs;
	const long nsa=parms->maos.nsa[ipowfs];
	aster->g->p[iwfs]=dref(star[istar].g->p[ipowfs]);
	aster->tsa+=nsa;
	aster->ngs[iwfs]=nsa*2;
    }    
    aster->gm=dcell2m(aster->g);
    
    if(aster->nwfs==1 && parms->maos.nmod==6 && aster->gm->nx==8){
	//there is a single ttf wfs and defocus needs to be estimated. remove degeneracy
	memset(aster->gm->p+2*aster->gm->nx, 0, sizeof(double)*aster->gm->nx);
    }
}
/**
   Copy information from star struct STAR_S to stars in asterism ASTER_S.
*/
void setup_aster_copystar(ASTER_S *aster, STAR_S *star, const PARMS_S *parms){
    (void)parms;
    int nwfs=aster->nwfs;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	const int ipowfs=aster->wfs[iwfs].ipowfs;
	const int istar=aster->wfs[iwfs].istar;
	/*Coordinate */
	aster->wfs[iwfs].thetax=star[istar].thetax;
	aster->wfs[iwfs].thetay=star[istar].thetay;
	/*Magnitude */
	aster->wfs[iwfs].mags=star[istar].mags;//do not free
	/*Signal Level */
	aster->wfs[iwfs].siglev=star[istar].siglev->p[ipowfs];//do not free
	aster->wfs[iwfs].siglevtot=star[istar].siglevtot->p[ipowfs];
	aster->wfs[iwfs].bkgrnd=star[istar].bkgrnd->p[ipowfs];
	
	/*Pixel intensity statistics. */
	aster->wfs[iwfs].pistat=&star[istar].pistat[ipowfs];
    }
}
/**
   Copy time history of complex pupil function from STAR_S to ASTER_S.
 */
void setup_aster_wvf(ASTER_S *aster, STAR_S *star, const PARMS_S *parms){
    (void) parms;
    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
	const int ipowfs=aster->wfs[iwfs].ipowfs;
	const int istar=aster->wfs[iwfs].istar;
	aster->wfs[iwfs].wvfout=star[istar].wvfout[ipowfs];
	aster->wfs[iwfs].ztiltout=star[istar].ztiltout[ipowfs];
	if(star[istar].goff){
	    aster->wfs[iwfs].goff=star[istar].goff->p[ipowfs];
	}
    }
}
/**
  Estimate wavefront error propagated from measurement error. pgm is the reconstructor. ineam is the
  error inverse.
  trace(Mcc*(pgm*neam*pgm'))
*/
static dmat *calc_recon_error(const dmat *pgm,   /**<[in] the reconstructor*/
			      const dmat *neam,/**<[in] the inverse of error covariance matrix*/
			      const dmat *mcc   /**<[in] NGS mode covariance matrix.*/
			      ){
    dmat *psp=NULL;
    dmat *tmp=NULL;
    dcp(&tmp, pgm);
    dmuldiag(tmp, neam);
    dmm(&psp, 0, tmp, pgm, "nt", 1);
    PDMAT(psp,ppsp); 
    PDMAT(mcc,pmcc);
    /*It is right for both ix, iy to stop at ib.*/
    double all[mcc->nx];
    for(int ib=0; ib<mcc->ny; ib++){
	all[ib]=0;
	for(int iy=0; iy<=ib; iy++){
	    for(int ix=0; ix<=ib; ix++){
		all[ib]+=ppsp[iy][ix]*pmcc[iy][ix];
	    }
	}
    }
    dfree(psp);
    dmat *res=dnew(3,1);
    res->p[0]=all[5];//total error
    res->p[1]=all[1];//total TT error
    if(mcc->nx>5){
	res->p[2]=all[5]-all[4];//focus alone
	if(res->p[2]<0){
	    res->p[2]=0;//due to coupling, this may be negative.
	}
    }
    return res;
}
/**
   Setup the reconstructor for each ASTER_S. 

   \f[R_{ngs}=(G_M^T C_n^{-1} G_M)^{-1})G_M^T C_n^{-1}\f]

   - \f$R_{ngs}\f$ is the reconstructor. 

   - \f$C_n^{-1}\f$ is the measurement error covariance matrix. Inclue both the
   measurement error and tilt anisoplanatism effects (computed as the variance
   of noise free gradient on PSFs that have ideal NGS compensation
   
   - \f$G_M^T\f$ is the NGS mode to gradient operator.
*/
void setup_aster_kalman_nea(ASTER_S *aster, STAR_S *star, const PARMS_S *parms){
    if(parms->skyc.multirate){
	dfree(aster->dtrats);
	if(aster->neam) dcellfreearr(aster->neam, 1);
	aster->neam=calloc(1, sizeof(dcell*));
	aster->neam[0]=dcellnew(aster->nwfs, aster->nwfs);
	aster->dtrats=dnew(aster->nwfs, 1);
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
	    const int istar=aster->wfs[iwfs].istar;
	    const int ipowfs=aster->wfs[iwfs].ipowfs;
	    const int idtrat=(int)star[istar].idtrat->p[ipowfs];
	    dmat *tmp=ddup(aster->wfs[iwfs].pistat->sanea->p[idtrat]);/*in rad*/
	    dcwpow(tmp, 2);
	    dsp *tmp2=spnewdiag(tmp->nx, tmp->p, 1);
	    spfull(&aster->neam[0]->p[iwfs+aster->nwfs*iwfs], tmp2,1);
	    dfree(tmp); spfree(tmp2);
	    aster->dtrats->p[iwfs]=parms->skyc.dtrats[idtrat];
	}
    }else{
	int ndtrat=parms->skyc.ndtrat;
	if(aster->neam) dcellfreearr(aster->neam, ndtrat);
	aster->neam=calloc(ndtrat, sizeof(dcell*));
	for(int idtrat=0; idtrat<ndtrat; idtrat++){
	    aster->neam[idtrat]=dcellnew(aster->nwfs, aster->nwfs);
	    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		dmat *tmp=ddup(aster->wfs[iwfs].pistat->sanea->p[idtrat]);/*in rad */
		dcwpow(tmp, 2);
		dsp *tmp2=spnewdiag(tmp->nx, tmp->p, 1);
		spfull(&aster->neam[idtrat]->p[iwfs+aster->nwfs*iwfs], tmp2,1);
		dfree(tmp); spfree(tmp2);
	    }
	}
    }
}
void setup_aster_lsr(ASTER_S *aster, STAR_S *star, const PARMS_S *parms){
    int ndtrat=parms->skyc.ndtrat;
    if(aster->pgm){
	dcellfree(aster->pgm);
	dcellfree(aster->sigman);
    }
    aster->pgm=dcellnew(ndtrat,1);
    aster->sigman=dcellnew(ndtrat,1);
    for(int idtrat=0; idtrat<ndtrat; idtrat++){
	dcell *nea=dcellnew(aster->nwfs, 1);
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
	    nea->p[iwfs]=ddup(aster->wfs[iwfs].pistat->sanea->p[idtrat]);
	    dcwpow(nea->p[iwfs], 2);
	}
	dmat *neam=dcell2m(nea);//measurement error covariance
	dcellfree(nea); 
	dcwpow(neam, -1);//inverse
	/*Reconstructor */
	aster->pgm->p[idtrat]=dpinv(aster->gm, neam, NULL);
	/*sigman is error due to noise. */
	dcwpow(neam, -1);//inverse again
	aster->sigman->p[idtrat]=calc_recon_error(aster->pgm->p[idtrat],neam,parms->maos.mcc);
	dfree(neam);
    }	
    if(parms->skyc.dbg){
	dcellwrite(aster->g,"%s/aster%d_g",dirsetup,aster->iaster);
	dcellwrite(aster->pgm,    "%s/aster%d_pgm", dirsetup,aster->iaster);
	dcellwrite(aster->sigman, "%s/aster%d_sigman", dirsetup,aster->iaster);
    }
}

/**
   Interpolate simu->gain.
*/
static void interp_gain(double *out, const dcell *gain, const dmat *gainx,
			double sigma2){
    const long nx=gainx->nx;
    const double xsep=(log(gainx->p[nx-1])-log(gainx->p[0]))/(nx-1);
    const double xx=(log(sigma2)-log(gainx->p[0]))/xsep;
    int ig;
    if(xx<0){/*too small. */
	ig=0;
    }else if(xx>=nx-1){/*too big */
	ig=nx-1;
    }else{/*within the range */
	ig=ifloor(xx);
	/*2013-12-06: use one of the set, not interpolate*/
    }
    memcpy(out, gain->p[ig]->p, sizeof(double)*gain->p[ig]->nx);
}
/**
   Optimize type II servo gains beased on measurement noise and signal PSD. We try to minimize
   \f[
   \sigma^2=\int \textrm{PSD}_{ngs,ws}H_{rej}\textrm{d}\nu + \int_0^{\nu_{nyquist}} \textrm{PSF}\textrm{d}\nu
   \f]
*/
static void setup_aster_servo(SIM_S *simu, ASTER_S *aster, const PARMS_S *parms){
    int ndtrat=parms->skyc.ndtrat;
    if(aster->gain){
	dcellfree(aster->gain);
	dfree(aster->res_ws);
	dfree(aster->res_ngs);
    }
    aster->gain=dcellnew(ndtrat,1);
    aster->res_ws=dnew(ndtrat,1);
    aster->res_ngs=dnew(ndtrat,3);
    PDMAT(aster->res_ngs, pres_ngs);
    for(int idtrat=0; idtrat<ndtrat; idtrat++){
	int dtrat=parms->skyc.dtrats[idtrat];
	double sigma_ngs= aster->sigman->p[idtrat]->p[0];
	double sigma_tt = aster->sigman->p[idtrat]->p[1];
	double sigma_ps = sigma_ngs-sigma_tt;
	double sigma_focus = aster->sigman->p[idtrat]->p[2];
	long nmod=parms->maos.nmod;
	/*gsplit:
	  0: All modes use the same gain.
	  1: PS, TT, focus (if nmod>5) use different gains. 
	  2: PS, TT use different gains. focus mode (if nmod>5) use PS gains.
	 */

	double res_ngs;/*residual error due to signal after servo rejection. */
	double res_ngsn;/*residual error due to noise. */
	const int servotype=parms->skyc.servo;
	const int ng=parms->skyc.ngain;
	aster->gain->p[idtrat]=dnew(ng,nmod);
	PDMAT(aster->gain->p[idtrat], pgain);
	if(parms->skyc.gsplit){
	    double pg_tt[ng+2];
	    double pg_ps[ng+2];
	    double pg_focus[ng+2];
	    if(parms->skyc.interpg){
		interp_gain(pg_tt, simu->gain_tt[idtrat], simu->gain_x, sigma_tt);
		interp_gain(pg_ps, simu->gain_ps[idtrat], simu->gain_x, sigma_ps);
		interp_gain(pg_focus, simu->gain_focus[idtrat], simu->gain_x, sigma_focus);
	    }else{
		dmat *sigma2=dnew(1,1); 
		dcell *tmp;
		sigma2->p[0]=sigma_tt;
		tmp=servo_optim(simu->psd_tt, parms->maos.dt, dtrat, parms->skyc.pmargin, sigma2, servotype);
		memcpy(pg_tt, tmp->p[0]->p, (ng+2)*sizeof(double)); dcellfree(tmp);

		sigma2->p[0]=sigma_ps;
		tmp=servo_optim(simu->psd_ps,    parms->maos.dt, dtrat, parms->skyc.pmargin, sigma2, servotype);
		memcpy(pg_ps, tmp->p[0]->p, (ng+2)*sizeof(double)); dcellfree(tmp);

		if(nmod>5){
		    sigma2->p[0]=sigma_focus;
		    tmp=servo_optim(simu->psd_focus, parms->maos.dt, dtrat, parms->skyc.pmargin, sigma2, servotype);
		    memcpy(pg_focus, tmp->p[0]->p, (ng+2)*sizeof(double)); dcellfree(tmp);
		}
		dfree(sigma2);
	    }
	    res_ngs  = pg_tt[ng] + pg_ps[ng] + pg_focus[ng];//residual mode
	    res_ngsn = pg_tt[ng+1] + pg_ps[ng+1] + pg_focus[ng+1];//error due to noise
	    for(int imod=0; imod<MIN(nmod,5); imod++){
		memcpy(pgain[imod], imod<2?pg_tt:pg_ps, sizeof(double)*ng);
	    }
	    if(nmod>5){
		memcpy(pgain[5], pg_focus, sizeof(double)*ng);
	    }
	}else{
	    double pg_ngs[ng+2];
	    if(parms->skyc.interpg){
		interp_gain(pg_ngs, simu->gain_ngs[idtrat], simu->gain_x, sigma_ngs);
	    }else{
		dmat *sigma2=dnew(1,1); sigma2->p[0]=sigma_ngs;
		dcell *tmp;
		tmp=servo_optim(simu->psd_ngs, parms->maos.dt, dtrat, parms->skyc.pmargin, sigma2, servotype);
		memcpy(pg_ngs, tmp->p[0]->p, (ng+2)*sizeof(double)); dcellfree(tmp);
	    }
	    res_ngs=pg_ngs[ng];
	    res_ngsn=pg_ngs[ng+1];
	    for(int imod=0; imod<nmod; imod++){
		memcpy(pgain[imod], pg_ngs, sizeof(double)*ng);
	    }
	}
	pres_ngs[0][idtrat]=res_ngs+res_ngsn;/*error due to signal and noise */
	pres_ngs[1][idtrat]=res_ngs;/*error due to signal */
	pres_ngs[2][idtrat]=res_ngsn;/*error due to noise propagation. */
	/*if(parms->skyc.reest){//estiamte error in time domain
	    dmat *sigma2=dnew(nmod,nmod);PDMAT(sigma2, psigma2);
	    PDMAT(parms->maos.mcc, pmcc);
	    //convert noise into mode space from WFE space.
	    psigma2[0][0]=psigma2[1][1]=sigma_tt/(2*pmcc[0][0]);
	    psigma2[2][2]=psigma2[3][3]=psigma2[4][4]=sigma_ps/(3*pmcc[2][2]);
	    if(nmod>5){
		psigma2[5][5]=sigma_focus/pmcc[5][5];
	    }
	    dmat *res=servo_test(simu->mideal, parms->maos.dt, dtrat, sigma2, aster->gain->p[idtrat]);
	    double rms=calc_rms(res,parms->maos.mcc);
	    pres_ngs[0][idtrat]=rms;
	    dfree(sigma2);
	    dfree(res);
	    }*/
	dmat *g_tt=dnew_ref(ng,1,pgain[0]);
	double gain_n;
	aster->res_ws->p[idtrat]=servo_residual(&gain_n, parms->skyc.psd_ws, 
						parms->maos.dt, dtrat, g_tt, 2);
	dfree(g_tt);
    }//for dtrat
    if(parms->skyc.dbg){
	dcellwrite(aster->gain,"%s/aster%d_gain",dirsetup,aster->iaster);
	dwrite(aster->res_ws,"%s/aster%d_res_ws",dirsetup,aster->iaster);
	dwrite(aster->res_ngs,"%s/aster%d_res_ngs",dirsetup,aster->iaster);
    }
}
static void setup_aster_kalman(SIM_S *simu, ASTER_S *aster, const PARMS_S *parms){
    int ndtrat=parms->skyc.ndtrat;
    
    aster->res_ngs=dnew(ndtrat,3);
    PDMAT(aster->res_ngs, pres_ngs);
    if(parms->skyc.multirate){
	aster->kalman=calloc(1, sizeof(kalman_t*));
	aster->kalman[0]=sde_kalman(simu->sdecoeff, parms->maos.dt, aster->dtrats, aster->g, aster->neam[0], 0);
    }else{
	aster->kalman=calloc(ndtrat, sizeof(kalman_t*));
	for(int idtrat=0; idtrat<ndtrat; idtrat++){
	    int dtrat=parms->skyc.dtrats[idtrat];
	    dmat *dtrats=dnew(aster->nwfs,1);
	    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		dtrats->p[iwfs]=dtrat;
	    }
	    aster->kalman[idtrat]=sde_kalman(simu->sdecoeff, parms->maos.dt, dtrats, aster->g, aster->neam[idtrat], 0);
	    dfree(dtrats);
	    dmat *res=kalman_test(aster->kalman[idtrat], simu->mideal);
	    pres_ngs[0][idtrat]=calc_rms(res, parms->maos.mcc);
	    dfree(res);
	}
    }
}
void setup_aster_controller(SIM_S *simu, ASTER_S *aster, STAR_S *star, const PARMS_S *parms){
    if(parms->skyc.servo<0){
	setup_aster_kalman_nea(aster, star, parms);
	setup_aster_kalman(simu, aster, parms);
    }else{
	setup_aster_lsr(aster, star, parms);
	setup_aster_servo(simu, aster, parms);
    }
}
/**
   for sort incrementally.
 */
static int sortfun(const double *a, const double *b){
    return a[0]<b[0]?-1:1;
}
/**
   Select a few asterisms that have decent performance (less than maxerror) */
int setup_aster_select(double *result, ASTER_S *aster, int naster, STAR_S *star, 
		       double maxerror, const PARMS_S *parms){
 
    int ndtrat=parms->skyc.ndtrat;
    dmat *res=dnew(ndtrat, naster);
    PDMAT(res,pres);
    dmat *imin=dnew(2,naster);
    PDMAT(imin, pimin);
    int master=-1;
    double mini;
    double minimum=INFINITY;
    //int dtrat_h=(parms->skyc.maxdtrat+1)/2;
    for(int iaster=0; iaster<naster; iaster++){
	mini=maxerror;
	aster[iaster].mdtrat=-1;
	for(int idtrat=0; idtrat<ndtrat; idtrat++){
	    /*should not add res_ws here since res_ngs already includes that.*/
	    double rms=aster[iaster].res_ngs->p[idtrat];
	    pres[iaster][idtrat]=rms;
	    if(rms<mini){
		mini=rms;
		aster[iaster].mdtrat=idtrat;
		aster[iaster].mresol=rms;
	    }
	}
	pimin[iaster][0]=mini;
	pimin[iaster][1]=iaster;
	if(mini<minimum){
	    master=iaster;
	    minimum=mini;
	}
	if(aster[iaster].mdtrat!=-1){
	    if(parms->skyc.maxdtrat>1){
		/*This is variance. allow a threshold */
		double thres=mini*3;
		double thres2=mini*3;
		if(thres>maxerror) thres=maxerror;
		if(thres2>maxerror) thres2=maxerror;
		/*Find upper and minimum good dtrats. */
		//MIN(ndtrat, aster[iaster].mdtrat+dtrat_h+1)
		for(int idtrat=aster[iaster].mdtrat; idtrat<ndtrat; idtrat++){
		    if(pres[iaster][idtrat]<thres){
			aster[iaster].idtratmax=idtrat;
		    }else{
			break;
		    }
		}
		//MAX(aster[iaster].mdtrat-dtrat_h, 0)
		for(int idtrat=aster[iaster].mdtrat; idtrat>=0; idtrat--){
		    if(pres[iaster][idtrat]<thres2){
			aster[iaster].idtratmin=idtrat;
		    }else{
			break;
		    }
		}
		
		if(aster[iaster].idtratmax>aster[iaster].idtratmin+parms->skyc.maxdtrat){
		    int mid=0.5*(aster[iaster].idtratmax+aster[iaster].idtratmin);
		    int min2=ceil(mid-parms->skyc.maxdtrat*0.5);
		    int max2=floor(mid+parms->skyc.maxdtrat*0.5);
		    /*info("aster%d, min=%d, max=%d, mid=%d; mdtrat=%d, min2=%d, max2=%d\n", 
			 iaster, aster[iaster].idtratmin, aster[iaster].idtratmax, 
			 mid, aster[iaster].mdtrat, min2, max2);*/
		    aster[iaster].idtratmin=min2;
		    aster[iaster].idtratmax=max2;
		}
	    }else{
		aster[iaster].idtratmin=aster[iaster].mdtrat;
		aster[iaster].idtratmax=aster[iaster].mdtrat+1;
	    }
	}else{
	    aster[iaster].idtratmin=0;
	    aster[iaster].idtratmax=0;
	}
	if(parms->skyc.verbose){
	    info2("aster%d, min=%d, max=%d, mdtrat=%d res=%g nm\n", 
		 iaster, aster[iaster].idtratmin, aster[iaster].idtratmax, 
		 aster[iaster].mdtrat, sqrt(mini)*1e9);
	}
    }
    qsort(imin->p, naster, 2*sizeof(double),(int(*)(const void*,const void*))sortfun);
    result[0]=minimum;
    result[1]=parms->skyc.fss[aster[master].mdtrat];
    int count=0;
    if(minimum<maxerror){
	double thres=MIN(minimum*3, maxerror);
	int taster=naster;
	if(parms->skyc.maxaster>0 && naster>parms->skyc.maxaster){
	    taster=parms->skyc.maxaster;
	}
	for(int jaster=0; jaster<taster; jaster++){
	    if(aster[jaster].mdtrat!=-1 && pimin[jaster][0]<thres){
		count++;
		int iaster=(int)pimin[jaster][1];
		aster[iaster].use=1;/*mark as valid. */
		for(int iwfs=0; iwfs<aster[iaster].nwfs; iwfs++){
		    int istar=aster[iaster].wfs[iwfs].istar;
		    int ipowfs=aster[iaster].wfs[iwfs].ipowfs;
		    star[istar].use[ipowfs]=1;
		}
	    }
	}
    }
    if(parms->skyc.verbose){
	info2("Minimum is found at aster %d at %.1f Hz: %.2f nm. Will evaluate %d asterisms.\n", 
	      master,result[1],sqrt(minimum)*1e9, count);
    }
    if(parms->skyc.dbg){
	dwrite(res, "%s/aster_resol",dirsetup);
    }
    dfree(res);
    dfree(imin);
    return count;
}
/**
   Free the ASTER_S array.
 */
void free_aster(ASTER_S *aster, int naster, const PARMS_S *parms){
    (void) parms;
    for(int iaster=0; iaster<naster; iaster++){
	int ndtrat=parms->skyc.ndtrat;
	if(aster[iaster].kalman){
	    if(parms->skyc.multirate){
		kalman_free(aster[iaster].kalman[0]);
		dcellfreearr(aster[iaster].neam, 1);
	
	    }else{
		for(int i=0; i<ndtrat; i++){
		    kalman_free(aster[iaster].kalman[i]);
		}
		dcellfreearr(aster[iaster].neam, ndtrat);
	    }
	    free(aster[iaster].kalman);
	    aster[iaster].kalman=0;
	}
	dcellfree(aster[iaster].gain);
	dcellfree(aster[iaster].pgm);
	dcellfree(aster[iaster].sigman);
	dfree(aster[iaster].res_ws);
	dfree(aster[iaster].res_ngs);

	free(aster[iaster].wfs);
	dcellfree(aster[iaster].g);
	dfree(aster[iaster].gm);
	dfree(aster[iaster].dtrats);
	free(aster[iaster].ngs);
    }
    free(aster);
}
