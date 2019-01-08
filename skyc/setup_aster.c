/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   Computes the number of possibilities of selection k items from n items \f$C_n^k\f$:
   factorial(n)/(factorial(k)*factorial(n-k)); 
*/
static long comb_select(long n, long k){
    return (long)round(factorial(n-k+1, n)/factorial(1, k));
}
/**
   initialize an initial combination composed a vector of non-negative numbers 0,1,2,...
*/
static int* comb_init(long k){
    int *comb=mycalloc(k,int);
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
ASTER_S *setup_aster_comb(int *naster, const STAR_S *star, int nstar, const PARMS_S *parms){
    if(nstar==0){
	*naster=0;
	return NULL;
    }else if(parms->skyc.keeporder){
	/*Use the same order as input stars.*/
	ASTER_S *aster=mycalloc(1,ASTER_S);
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
	aster[0].wfs=mycalloc(ntot,WFS_S);
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
	info("Number of stars: %d, number of asterisms: %d\n", nstar, ncomb);
    }
    aster=mycalloc(ncomb,ASTER_S);
    int count=0;
    int *comb=comb_init(nwfstot);//select nwfstot stars from all available stars
    do{
	if(npowfs==1){
	    aster[count].nwfs=nwfs[0];
	    aster[count].wfs=mycalloc(nwfs[0], WFS_S);
	    aster[count].iaster=count;
	    int skip=0;
	    for(int iwfs=0; iwfs<nwfs[0]; iwfs++){
		aster[count].wfs[iwfs].ipowfs=0;
		aster[count].wfs[iwfs].istar=comb[iwfs];
		if(star[aster[count].wfs[iwfs].istar].use[aster[count].wfs[iwfs].ipowfs]==-1){
		    skip=1;
		}
	    }
	    if(!skip) count++; else free(aster[count].wfs);
	}else if(npowfs==2){
	    int mask[nwfstot];//mask stars that are selected already.
	    int *comb2=comb_init(nwfs[0]);//select nwfs[0] from nwfstot for first ipowfs.
	    do{
		int skip=0;
		memset(mask, 0, sizeof(int)*nwfstot);
		aster[count].nwfs=nwfstot;
		aster[count].wfs=mycalloc(nwfstot,WFS_S);
		aster[count].iaster=count;
		for(int iwfs=0; iwfs<nwfs[0]; iwfs++){
		    aster[count].wfs[iwfs].ipowfs=0;
		    aster[count].wfs[iwfs].istar=comb[comb2[iwfs]];
		    if(star[aster[count].wfs[iwfs].istar].use[aster[count].wfs[iwfs].ipowfs]==-1){
			skip=1;
		    }
		    mask[comb2[iwfs]]=1;
		}
		int jstar=0;
		for(int iwfs=0; iwfs<nwfs[1]; iwfs++){
		    aster[count].wfs[iwfs+nwfs[0]].ipowfs=1;
		    while(mask[jstar]) jstar++;
		    aster[count].wfs[iwfs+nwfs[0]].istar=comb[jstar];
		    if(star[aster[count].wfs[iwfs].istar].use[aster[count].wfs[iwfs].ipowfs]==-1){
			skip=1;
		    }
		    mask[jstar]=1;
		}
		if(!skip) count++; else free(aster[count].wfs);
	    }while(comb_next(comb2,nwfstot,nwfs[0]) && !parms->skyc.ttfbrightest);
	    free(comb2);
	}
    }while(comb_next(comb,nstar,nwfstot));
    free(comb);
    *naster=count;
    return aster;
}
/**
   Compute Modal to gradient operator by copying from the stars. Using average
gradients. Similar to Z tilt since the mode is low order */
void setup_aster_g(ASTER_S *aster, STAR_S *star, const PARMS_S *parms){
    /*2010-06-08: Tested against MATLAB skycoverage code. */
    aster->g=dcellnew(aster->nwfs,1);
    /*if(parms->maos.nmod<5 || parms->maos.nmod>6){
	error("Not compatible with the number of NGS modes\n");
	}*/
    aster->ngs=mycalloc(aster->nwfs,long);
    aster->tsa=0;
    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
	const int istar=aster->wfs[iwfs].istar;
	const int ipowfs=aster->wfs[iwfs].ipowfs;
	const long nsa=parms->maos.nsa[ipowfs];
	aster->g->p[iwfs]=ddup(star[istar].g->p[ipowfs]);
	aster->tsa+=nsa;
	aster->ngs[iwfs]=nsa*2;
    }    
    /*
      aster->g is also used for simulation. Do not zero columns here.
    */
    aster->gm=dcell2m(aster->g);
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
    }
}
/**
   Copy time history of complex pupil function from STAR_S to ASTER_S.
 */
void setup_aster_ztilt(ASTER_S *aster, STAR_S *star, const PARMS_S *parms){
    (void) parms;
    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
	const int ipowfs=aster->wfs[iwfs].ipowfs;
	const int istar=aster->wfs[iwfs].istar;
	aster->wfs[iwfs].ztiltout=star[istar].ztiltout->p[ipowfs];
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
    dmat *var=NULL;
    dcp(&tmp, pgm);
    dmuldiag(tmp, neam);
    dmm(&psp, 0, tmp, pgm, "nt", 1);
    dfree(tmp);
    dmm(&var, 0, mcc, psp, "nn", 1);
    dfree(psp);
    dmat *res=dnew(mcc->nx+1,1);
    /*It is right for both ix, iy to stop at ib.*/
    for(int ib=0; ib<mcc->ny; ib++){
	res->p[ib]=IND(var, ib, ib);
	res->p[mcc->ny]+=res->p[ib];//trace
	if(res->p[ib]<0){
	    res->p[ib]=fabs(res->p[ib]);
	    /*
	    writebin(pgm, "pgm");
	    writebin(neam, "neam");
	    writebin(mcc, "mcc");
	    for(int ib2=0; ib2<mcc->nx; ib2++){
		warning("mod[%d]=%g nm^2\n", ib2, res->p[ib2]*1e18);
	    }
	    exit(0);*/
	}
    }
    dfree(var);
    return res;
}

void setup_aster_lsr(ASTER_S *aster, const PARMS_S *parms){
    int ndtrat=parms->skyc.ndtrat;
    if(aster->pgm){
	dcellfree(aster->pgm);
	dcellfree(aster->sigman);
    }
    aster->pgm=dcellnew(ndtrat,1);
    aster->sigman=dcellnew(ndtrat,1);
    dmat *gm=ddup(aster->gm);
    if(aster->nwfs==1 && parms->maos.nmod==6 && gm->nx==8){
	//focus and magnification mode degenerate in single GS mode.
	info("set 3rd column of gm to zero\n");
	memset(gm->p+gm->nx*2, 0, sizeof(double)*gm->nx);
    }

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
	aster->pgm->p[idtrat]=dpinv(gm, neam);
	/*sigman is error due to noise. */
	dcwpow(neam, -1);//inverse again
	aster->sigman->p[idtrat]=calc_recon_error(aster->pgm->p[idtrat],neam,parms->maos.mcc);
	dfree(neam);
    }	
    if(parms->skyc.dbg){
	writebin(gm,"%s/aster%d_gm",dirsetup,aster->iaster);
	writebin(aster->pgm,    "%s/aster%d_pgm", dirsetup,aster->iaster);
	writebin(aster->sigman, "%s/aster%d_sigman", dirsetup,aster->iaster);
    }
    dfree(gm);
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
    dmat*  pres_ngs=aster->res_ngs;
    for(int idtrat=0; idtrat<ndtrat; idtrat++){
	int dtrat=parms->skyc.dtrats->p[idtrat];
	long nmod=parms->maos.nmod;
	/*gsplit:
	  0: All modes use the same gain.
	  1: Different mode use defferent gain (2017-0-24) was only tt/ps separate.
	  note that simu->psds->p[0] contains windshake PSD.
	*/
	double res_ngs=0;/*residual error due to signal after servo rejection. */
	double res_ngsn=0;/*residual error due to noise. */
	const int servotype=parms->skyc.servo;
	const int ng=parms->skyc.ngain;//number of gain parameters
	dmat*  pgain=aster->gain->p[idtrat]=dnew(ng,nmod);
	for(int ipsd=0; ipsd<simu->psds->nx; ipsd++){
	    double sigma=aster->sigman->p[idtrat]->p[parms->skyc.gsplit?ipsd:nmod];
	    double pg[ng+2];
	    
	    if(parms->skyc.interpg){
		interp_gain(pg, simu->gain_pre->p[idtrat]->p[ipsd], simu->gain_x, sigma);
	    }else{
		dmat *sigma2=dnew(1,1); 
		sigma2->p[0]=sigma;
		dcell *tmp=servo_optim(simu->psds->p[ipsd], parms->maos.dt, dtrat, parms->skyc.pmargin, sigma2, servotype);
		memcpy(pg, tmp->p[0]->p, (ng+2)*sizeof(double)); 
		dcellfree(tmp);
		dfree(sigma2);
	    }
	    res_ngs  += pg[ng];
	    res_ngsn += pg[ng+1];
	    memcpy(PCOL(pgain,ipsd), pg, sizeof(double)*ng);
	}
	
	for(int imod=simu->psds->nx; imod<nmod; imod++){
	    memcpy(PCOL(pgain,imod), PCOL(pgain,0), sizeof(double)*ng);
	}
	IND(pres_ngs,idtrat,0)=res_ngs+res_ngsn;/*error due to signal and noise */
	IND(pres_ngs,idtrat,1)=res_ngs;/*error due to signal */
	IND(pres_ngs,idtrat,2)=res_ngsn;/*error due to noise propagation. */

	dmat *g_tt=dnew_ref(ng,1,PCOL(pgain,0));
	double gain_n;
	aster->res_ws->p[idtrat]=servo_residual(&gain_n, parms->skyc.psd_ws, 
						parms->maos.dt, dtrat, g_tt, servotype);
	dfree(g_tt);
    }//for dtrat
    if(parms->skyc.dbg){
	writebin(aster->gain,"%s/aster%d_gain",dirsetup,aster->iaster);
	writebin(aster->res_ws,"%s/aster%d_res_ws",dirsetup,aster->iaster);
	writebin(aster->res_ngs,"%s/aster%d_res_ngs",dirsetup,aster->iaster);
    }
}
static void setup_aster_kalman_dtrat(ASTER_S *aster, const PARMS_S *parms, int idtrat_wfs0){
    if(parms->skyc.verbose){
	info("aster %d dtrat_wfs0=%3d, dtrat=", aster->iaster, 
	      (int)parms->skyc.dtrats->p[idtrat_wfs0]);
    }
    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
	int idtrat=idtrat_wfs0;
	if(iwfs>0){
	    /*don't allow snr to fall below 3.*/
	    while(idtrat>0 && (aster->wfs[iwfs].pistat->snr->p[idtrat]<3
			       || (int)parms->skyc.dtrats->p[idtrat] % (int)aster->dtrats->p[0]!=0)){
		idtrat--;
	    }
	}
	aster->idtrats->p[iwfs]=idtrat;
	aster->dtrats->p[iwfs]=parms->skyc.dtrats->p[idtrat];
	int ng=aster->g->p[iwfs]->nx;
	if(idtrat>-1){
	    for(int ig=0; ig<ng; ig++){
		aster->neam[0]->p[iwfs+aster->nwfs*iwfs]->p[ig*(ng+1)]=
		    pow(aster->wfs[iwfs].pistat->sanea->p[idtrat]->p[ig], 2);
	    }
	}else{//no star available
	    dset(aster->neam[0]->p[iwfs+aster->nwfs*iwfs], 0);
	}
	if(parms->skyc.verbose){
	    info("%3d ", (int)parms->skyc.dtrats->p[idtrat]);
	}
    }//for iwfs
}

static void setup_aster_kalman(SIM_S *simu, ASTER_S *aster, STAR_S *star, const PARMS_S *parms){
    int ndtrat=parms->skyc.ndtrat;
    if(parms->skyc.multirate){
	aster->res_ngs=dnew(1,1);
	//assemble neam
	if(aster->neam) error("neam is already set?\n");
	aster->neam=mycalloc(1,dcell*);
	aster->neam[0]=dcellnew(aster->nwfs, aster->nwfs);
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
	    int ng=aster->g->p[iwfs]->nx;
	    aster->neam[0]->p[iwfs+aster->nwfs*iwfs]=dnew(ng, ng);
	}
	aster->dtrats=dnew(aster->nwfs, 1);
	aster->idtrats=dnew(aster->nwfs, 1);
	int wfs0_min=0, wfs0_max=0;
	PISTAT_S *pistat0=&star[aster->wfs[0].istar].pistat[aster->wfs[0].ipowfs];
	for(int idtrat=0; idtrat<parms->skyc.ndtrat; idtrat++){
	    if(wfs0_min==0 && pistat0->snr->p[idtrat]>3){
		wfs0_min=idtrat;
	    }
	    if(pistat0->snr->p[idtrat]>=parms->skyc.snrmin->p[idtrat]){
		wfs0_max=idtrat;
	    }
	}
	aster->kalman=mycalloc(1,kalman_t*);
	double resmin=INFINITY;
	kalman_t *kalman_min=0;
	int idtrat_min=0;
	//Try progressively lower sampling frequencies until performance starts to degrades
	for(int idtrat_limit=wfs0_max; idtrat_limit>wfs0_min; idtrat_limit--){
	    setup_aster_kalman_dtrat(aster, parms, idtrat_limit);
	    aster->kalman[0]=sde_kalman(simu->sdecoeff, parms->maos.dt, aster->dtrats, aster->g, aster->neam[0], 0);
	    dmat *rests=0;
#define USE_SIM 0 //ztiltout is not available here.
#if USE_SIM   //more accurate
	    dmat *res=skysim_sim(parms->skyc.dbg?&rests:0, simu->mideal, simu->mideal_oa, simu->varol,
				 aster, 0, parms, -1, 1, -1);
	    double res0=res?res->p[0]:simu->varol;
	    dfree(res);
#else
	    rests=kalman_test(aster->kalman[0], simu->mideal);
	    double res0=calc_rms(rests, parms->maos.mcc, parms->skyc.evlstart);
#endif
	    if(parms->skyc.dbg){
		writebin(rests, "isky%d_iaster%d_dtrat%d_rest", simu->isky, aster->iaster, idtrat_limit);
	    }
	    if(parms->skyc.verbose) info("res0=%g, resmin=%g\n", sqrt(res0)*1e9, sqrt(resmin)*1e9);
	    dfree(rests);
	    if(res0<resmin-100e-18){//better by 10 nm
		resmin=res0;
		kalman_free(kalman_min);
		kalman_min=aster->kalman[0];
		aster->kalman[0]=0;
		idtrat_min=idtrat_limit;
	    }else{
		kalman_free(aster->kalman[0]);aster->kalman[0]=0;
		if(isfinite(resmin) && res0>resmin*2){//stop trying.
		    break;
		}
	    }
	}
	setup_aster_kalman_dtrat(aster, parms, idtrat_min);
	if(parms->skyc.verbose) info("selected\n");
	aster->res_ngs->p[0]=resmin;
	aster->kalman[0]=kalman_min;
	if(parms->skyc.dbg){
	    kalman_write(aster->kalman[0],"%s/aster%d_kalman_%g",dirsetup,aster->iaster, parms->skyc.dtrats->p[idtrat_min]);
	}
    }else{
	if(aster->neam) dcellfreearr(aster->neam, ndtrat);
	aster->neam=mycalloc(ndtrat,dcell*);
	aster->res_ngs=dnew(ndtrat,3);
	dmat*  pres_ngs=aster->res_ngs;
	aster->kalman=mycalloc(ndtrat,kalman_t*);
	dmat *dtrats=dnew(aster->nwfs,1);
	for(int idtrat=0; idtrat<ndtrat; idtrat++){
	    //assemble neam
	    //TIC;tic;
	    aster->neam[idtrat]=dcellnew(aster->nwfs, aster->nwfs);
	    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		dmat *tmp=ddup(aster->wfs[iwfs].pistat->sanea->p[idtrat]);/*in rad */
		dcwpow(tmp, 2);
		dsp *tmp2=dspnewdiag(tmp->nx, tmp->p, 1);
		dspfull(&aster->neam[idtrat]->p[iwfs+aster->nwfs*iwfs], tmp2,'n',1);
		dfree(tmp); dspfree(tmp2);
	    }

	    int dtrat=parms->skyc.dtrats->p[idtrat];
	    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		dtrats->p[iwfs]=dtrat;
	    }
	    aster->kalman[idtrat]=sde_kalman(simu->sdecoeff, parms->maos.dt, dtrats, aster->g, aster->neam[idtrat], 0);
	    //toc("kalman");
#if USE_SIM 
	    dmat *res=skysim_sim(0, simu->mideal, simu->mideal_oa, simu->varol, aster, 0, parms, idtrat, 1, -1);
	    double rms=res?res->p[0]:simu->varol;
	    dfree(res);
#else
	    dmat *res=kalman_test(aster->kalman[idtrat], simu->mideal);
	    double rms=calc_rms(res, parms->maos.mcc, parms->skyc.evlstart);
	    dfree(res);
#endif
	    //toc("estimate");
	    IND(pres_ngs,idtrat,0)=rms;
	    if(parms->skyc.dbg){
		kalman_write(aster->kalman[idtrat],"%s/aster%d_kalman_%d",dirsetup,aster->iaster,dtrat);
	    }
	}
	dfree(dtrats);
    }
}
void setup_aster_controller(SIM_S *simu, ASTER_S *aster, STAR_S *star, const PARMS_S *parms){
    if(parms->skyc.servo<0){
	setup_aster_kalman(simu, aster, star, parms);
    }else{
	setup_aster_lsr(aster, parms);
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
    dmat *imin=dnew(2,naster);//record best performance of each asterism.
    dset(imin, INFINITY);
    int master=-1;
    double fieldMinRes=INFINITY;//minimum of this field.
    const double wvfMargin=6400e-18; //Allow asterism that is worse by this much to be evaluated.
	
    for(int iaster=0; iaster<naster; iaster++){
	aster[iaster].mdtrat=-1;//idtrat at asterMinRes
	if(parms->skyc.dbgaster>-1 && iaster != parms->skyc.dbgaster){
	    continue;
	}
	double asterMinRes=maxerror;//minimum of this asterism.
	for(int idtrat=0; idtrat<(parms->skyc.multirate?1:ndtrat); idtrat++){
	    /*should not add res_ws here since res_ngs already includes that.*/
	    double wfv=aster[iaster].res_ngs->p[idtrat];
	    IND(res,idtrat,iaster)=wfv;
	    if(wfv<asterMinRes){
		asterMinRes=wfv;
		aster[iaster].mdtrat=idtrat;
		aster[iaster].mresest=wfv;
	    }
	}
	IND(imin,0,iaster)=asterMinRes;
	IND(imin,1,iaster)=iaster;
	if(asterMinRes<fieldMinRes){
	    master=iaster;
	    fieldMinRes=asterMinRes;
	}
	if(!parms->skyc.multirate){
	    if(aster[iaster].mdtrat!=-1){
		if(parms->skyc.maxdtrat>1){
		    /*This is variance. allow a threshold */
		    double thres=MIN(asterMinRes+wvfMargin, maxerror);//threshold at low freq end
		    double thres2=MIN(asterMinRes+wvfMargin, maxerror);//threshold at high freq end
		    /*Find upper and fieldMinRes good dtrats. */
		    for(int idtrat=aster[iaster].mdtrat; idtrat<ndtrat; idtrat++){
			if(IND(res,idtrat,iaster)<thres){
			    aster[iaster].idtratmax=idtrat+1;
			}else{
			    break;
			}
		    }
		    for(int idtrat=aster[iaster].mdtrat; idtrat>=0; idtrat--){
			if(IND(res,idtrat,iaster)<thres2){
			    aster[iaster].idtratmin=idtrat;
			}else{
			    break;
			}
		    }
		    //2018-02-28: changed to prefer high frame rate
		    if(aster[iaster].idtratmax>aster[iaster].idtratmin+parms->skyc.maxdtrat+1){
			aster[iaster].idtratmin=aster[iaster].idtratmax-(parms->skyc.maxdtrat+1);
		    }
		}else{
		    aster[iaster].idtratmin=aster[iaster].mdtrat;
		    aster[iaster].idtratmax=aster[iaster].idtratmin+1;
		}
	    }
	}else{//multirate
	    aster[iaster].mdtrat=aster[iaster].idtrats->p[0];
	    aster[iaster].idtratmin=aster[iaster].idtrats->p[0];
	    aster[iaster].idtratmax=aster[iaster].idtratmin+1;
	}
	if(parms->skyc.verbose){
	    info("aster%2d, dtrats=[%2d,%2d,%2d], res=[%.2f %.2f %.2f]nm\n", 
		  iaster, (int)parms->skyc.dtrats->p[aster[iaster].idtratmin],
		  (int)parms->skyc.dtrats->p[aster[iaster].mdtrat],
		  (int)parms->skyc.dtrats->p[aster[iaster].idtratmax-1], 
		  sqrt(IND(res,aster[iaster].idtratmin,iaster))*1e9,
		  sqrt(IND(res,aster[iaster].mdtrat,iaster))*1e9,
		  sqrt(IND(res,aster[iaster].idtratmax-1,iaster))*1e9);
	}
    }
    if(parms->skyc.dbgsky>-1){
	writebin(imin, "sky%d_imin", parms->skyc.dbgsky);
    }
  
    result[0]=fieldMinRes;
    result[1]=master;
    if(aster[master].mdtrat!=-1){
	result[2]=parms->skyc.fss[aster[master].mdtrat];
    }
  
    //Mark asterisms and stars to be evaluated further.
    int count=0;
    if(fieldMinRes<maxerror){
	double thres=MIN(fieldMinRes+wvfMargin, maxerror);
	int taster=naster;
	if(parms->skyc.maxaster>0 && naster>parms->skyc.maxaster){
	    taster=parms->skyc.maxaster;
	}
	qsort(imin->p, naster, 2*sizeof(double),(int(*)(const void*,const void*))sortfun);
	for(int jaster=0; jaster<taster; jaster++){
	    if(IND(imin,0,jaster)>thres) continue;
	    int iaster=(int)IND(imin,1,jaster);
	    if(aster[iaster].mdtrat==-1) continue; 
	    count++;
	    aster[iaster].use=1;/*mark as valid. */
	    char temp1[1024]; temp1[0]='\0';
	    for(int iwfs=0; iwfs<aster[iaster].nwfs; iwfs++){
		int istar=aster[iaster].wfs[iwfs].istar;
		int ipowfs=aster[iaster].wfs[iwfs].ipowfs;
		star[istar].use[ipowfs]=1;
		if(parms->skyc.estimate && parms->skyc.nsky==1){
		    char temp2[1024];
		    snprintf(temp2, 1023,"(x=%.1f, y=%.1f, J=%.1f) ", star[istar].thetax*206265, star[istar].thetay*206265, star[istar].mags->p[0]);
		    strncat(temp1, temp2, 1023);
		}
	    }
	    if(parms->skyc.estimate){
		info("%s, %g Hz, %g nm\n", temp1, parms->skyc.fss[aster[iaster].mdtrat], sqrt(aster[iaster].mresest)*1e9);
	    }
	}
    }
    if(parms->skyc.verbose){
	info("FieldMinRes is found at aster %g at %.1f Hz: %.2f nm. Will evaluate %d asterisms.\n", 
	      result[1],result[2],sqrt(result[0])*1e9, count);
    }
    if(parms->skyc.dbg){
	writebin(res, "%s/aster_resol",dirsetup);
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
	dfree(aster[iaster].idtrats);
	free(aster[iaster].ngs);
	dcellfree(aster[iaster].phyRes);
	dcellfree(aster[iaster].phyMRes);
    }
    free(aster);
}
