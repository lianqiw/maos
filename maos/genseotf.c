/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

/*
  2009-11-26: changed to rotate OTF instead of psf to comply with the approach
  in wfsints.  this gives slightly larger intensity because max value of OTF is
  preserved which corresponds to the sum of psf.  */
/**
   \file genseotf.c contains routines to generate mean short exposure (tip/tilt
   removed) pixel intensities. Mostly used for LGS pixel intensity for its
   matched filter. Structure functions from kolmogorov spectrum is used. Not
   able to take into account outerscale yet.

   \todo find ways to factor in outerscale effect (use von karman spectrum
   instead of kolmogorov) */

#include "common.h"
#include "genseotf.h"

/**
   Master routine that generates short exposure OTF by calling genotf() in the
   library with p/t/t removal set.
*/
static void genseotf_do(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    /*create a grid representing the sub-aperture. */
    loc_t *loc=mksqloc(powfs[ipowfs].pts->nx,
		       powfs[ipowfs].pts->nx,
		       powfs[ipowfs].pts->dx,
		       powfs[ipowfs].pts->dy,
		       0,0);
    /*size of the OTF grid */
    //int ncompx=powfs[ipowfs].pts->nx*parms->powfs[ipowfs].embfac;
    //int ncompy=ncompx;
    //2015-04-27:replaced by powfs[ipowfs].ncompx/y, which maybe different
    
    /*The embeding factor for embedding the aperture */
    const int embfac=parms->powfs[ipowfs].embfac;
    const double dxsa=powfs[ipowfs].pts->dsa;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    int nsa=powfs[ipowfs].saloc->nloc;
 
    int notf=1;
    int has_ncpa=0;
    if(powfs[ipowfs].opdbias && parms->powfs[ipowfs].ncpa_method==2){
	//check whether opdbias is different between wfs[0] and following.
	int different=0;
	for(int iwfs=1; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
	    double diff=ddiff(powfs[ipowfs].opdbias->p[0], powfs[ipowfs].opdbias->p[iwfs]);
	    if(diff>1e-4){
		info("powfs[%d].opdbias[%d] is different from powfs[%d].opdbias[0] by %g.\n", 
		     ipowfs, iwfs, ipowfs, diff);
		different=1;
	    }
	}
	if(different){
	    notf=MAX(notf, parms->powfs[ipowfs].nwfs);
	}
	has_ncpa=1;
    }
    if(powfs[ipowfs].loc_tel){
	notf=MAX(notf,parms->powfs[ipowfs].nwfs);
    }
    info2("notf=%d\n", notf);
    if(powfs[ipowfs].intstat->otf){
	cellfree(powfs[ipowfs].intstat->otf);
    }
    powfs[ipowfs].intstat->notf=notf;
    powfs[ipowfs].intstat->otf=cccellnew(notf, 1);
    for(int iotf=0; iotf<notf; iotf++){
	powfs[ipowfs].intstat->otf->p[iotf]=ccellnew(nsa,nwvl);
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	double wvl=parms->powfs[ipowfs].wvl->p[iwvl];
	double dtheta=wvl/(dxsa*embfac);
	for(int iotf=0; iotf<notf; iotf++){
	    dmat* opdbias=has_ncpa?powfs[ipowfs].opdbias->p[iotf]:NULL;
	    double thres=opdbias?1:(1-1e-10);
	    info2("There is %s bias\n", opdbias?"NCPA":"no");
	    OMPTASK_SINGLE
		genotf(powfs[ipowfs].intstat->otf->p[iotf]->p+iwvl*nsa,
		       loc, powfs[ipowfs].realamp->p[iotf], opdbias, 
		       powfs[ipowfs].realsaa->p[iotf],
		       thres,wvl,dtheta,NULL,parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0, 
		       powfs[ipowfs].ncompx, 
		       powfs[ipowfs].ncompy, nsa, 1);
	}
    }/*iwvl */
    locfree(loc);
}

void genseotf(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    int npsfx=powfs[ipowfs].pts->nx *parms->powfs[ipowfs].embfac;
    int npsfy=npsfx;
    char fnprefix[PATH_MAX]; fnprefix[0]='\0';
    uint32_t key=0;
    strcat(fnprefix, "SEOTF");
    if(powfs[ipowfs].amp_tel){
	for(int iamp=0; iamp<parms->powfs[ipowfs].nwfs; iamp++){
	    key=dhash(powfs[ipowfs].amp_tel->p[iamp], key);
	}
    }else{
	key=dhash(powfs[ipowfs].amp, key);
    }
    info2("powfs %d: ncpa_method=%d, opdbias=%p\n",
	 ipowfs, parms->powfs[ipowfs].ncpa_method, powfs[ipowfs].opdbias);
    if(powfs[ipowfs].opdbias && parms->powfs[ipowfs].ncpa_method==2){
	for(int iwfs=0; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
	    key=dhash(powfs[ipowfs].opdbias->p[iwfs],key);
	}
    }
    if(key!=0){
	char tmp2[80];
	snprintf(tmp2,80,"_%ud",key);
	strncat(fnprefix, tmp2, PATH_MAX-strlen(fnprefix)-1);
    }
    char fnotf[PATH_MAX];
    char fnlock[PATH_MAX];
    snprintf(fnotf,PATH_MAX,"%s/SEOTF/",CACHE);
    if(!exist(fnotf)) {
	mymkdir("%s",fnotf);
    }else{
	remove_file_older(fnotf, 365*24*3600);
    }
    long nsa=powfs[ipowfs].saloc->nloc;
    snprintf(fnotf,PATH_MAX,"%s/SEOTF/%s_D%g_%g_"
	     "r0_%g_L0%g_dsa%g_nsa%ld_dx1_%g_"
	     "nwvl%d_%g_embfac%d_ncompx%d_%dx%d_v2",
	     CACHE, fnprefix,
	     parms->aper.d,parms->aper.din, 
	     parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0, 
	     powfs[ipowfs].pts->dsa,nsa,
	     1./powfs[ipowfs].pts->dx, 
	     parms->powfs[ipowfs].nwvl,
	     parms->powfs[ipowfs].wvl->p[0]*1.e6,
	     parms->powfs[ipowfs].embfac,
	     powfs[ipowfs].ncompx, npsfx,npsfy);
    snprintf(fnlock, PATH_MAX, "%s.lock", fnotf);
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    while(!intstat->otf){
	if(exist(fnlock) || !zfexist(fnotf)){/*need to create data */
	    int fd=lock_file(fnlock, 0, 0);/*nonblocking exclusive lock */
	    if(fd>=0){/*succeed */
		info2("Generating WFS OTF for %s...", fnotf);
		TIC;tic; genseotf_do(parms,powfs,ipowfs); toc2("done");
		writebin(intstat->otf, "%s", fnotf);
	    }else{
		warning("Waiting for previous lock to release ...");
		fd=lock_file(fnlock, 1, 0);
	    }
	    close(fd); remove(fnlock);
	}else{
	    info2("Reading WFS OTF from %s\n", fnotf);
	    intstat->otf=cccellread("%s",fnotf);
	    intstat->notf=intstat->otf->nx*intstat->otf->ny;
	    if(!intstat->notf) error("Invalid otf\n");
	}
    }
}
/**
   Creating short exposure OTF caused by turbulence within LLT uplink aperture
*/
void genselotf_do(const PARMS_T *parms,POWFS_T *powfs,int ipowfs){
    if(!parms->powfs[ipowfs].llt) return;
    loc_t *loc=pts2loc(powfs[ipowfs].llt->pts);
    const int notf=powfs[ipowfs].llt->pts->nx*parms->powfs[ipowfs].embfac;
    const int nwvl=parms->powfs[ipowfs].nwvl;

    int nlotf=1;
    dcell *ncpa=powfs[ipowfs].llt->ncpa;
    if(ncpa){
	nlotf=ncpa->nx*ncpa->ny;
    }
    if(powfs[ipowfs].intstat->lotf){
	ccellfree(powfs[ipowfs].intstat->lotf);
    }
    powfs[ipowfs].intstat->lotf=ccellnew(nwvl,nlotf);
    ccell*  lotf=powfs[ipowfs].intstat->lotf/*PCELL*/;
    if(nwvl!=1){
	warning("LGS has multi-color!\n");
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	double wvl=parms->powfs[ipowfs].wvl->p[iwvl];
	double dtheta=wvl/(notf*powfs[ipowfs].llt->pts->dx);
	double thres=1;
	for(int ilotf=0; ilotf<nlotf; ilotf++){
	    genotf(PIND(lotf,iwvl,ilotf), loc, powfs[ipowfs].llt->amp, ncpa?ncpa->p[ilotf]:NULL, 
		   0, thres, wvl, dtheta, NULL,parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0,
		   notf, notf, 1, 1);
	}
    }/*iwvl */
    locfree(loc);
}
void genselotf(const PARMS_T *parms,POWFS_T *powfs,int ipowfs){
    if(!parms->powfs[ipowfs].llt){
	return;
    }
    char fnprefix[80];
    uint32_t key=0;
    key=dhash(powfs[ipowfs].llt->amp, key);
    if(powfs[ipowfs].llt->ncpa){
	dcell *ncpa=powfs[ipowfs].llt->ncpa;
	long nlotf=ncpa->nx*ncpa->ny;
	for(long ilotf=0; ilotf<nlotf; ilotf++){
	    key=dhash(ncpa->p[ilotf], key);
	}
    }
    snprintf(fnprefix,80,"SELOTF_%0x",key);
    char fnlotf[PATH_MAX];
    snprintf(fnlotf,PATH_MAX,"%s/SELOTF/", CACHE);
    if(!exist(fnlotf)){
	mymkdir("%s",fnlotf);
    }else{
	remove_file_older(fnlotf, 365*24*3600);
    }
    snprintf(fnlotf,PATH_MAX,"%s/SELOTF/%s_"
	     "r0_%g_L0%g_lltd%g_dx1_%g_W%g_"
	     "nwvl%d_%g_embfac%d_v2", 
	     CACHE, fnprefix,
	     parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0, 
	     powfs[ipowfs].llt->pts->dsa,
	     1./powfs[ipowfs].llt->pts->dx,
	     parms->powfs[ipowfs].llt->widthp,
	     parms->powfs[ipowfs].nwvl,
	     parms->powfs[ipowfs].wvl->p[0]*1.e6,
	     parms->powfs[ipowfs].embfac);
    char fnlock[PATH_MAX];
    snprintf(fnlock, PATH_MAX, "%s.lock", fnlotf);
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    while(!intstat->lotf){
	if(exist(fnlock) || !zfexist(fnlotf)){/*need to create data */
	    int fd=lock_file(fnlock, 0, 0);/*nonblocking exclusive lock */
	    if(fd>=0){/*succeed */
		info2("Generating WFS LLT OTF for %s\n", fnlotf);
		genselotf_do(parms,powfs,ipowfs);
		writebin(intstat->lotf, "%s",fnlotf);
	    }else{//waiting by retry locking with blocking.
		warning("Waiting for previous lock to release ...");
		fd=lock_file(fnlock, 1, 0);
	    }
	    close(fd); remove(fnlock);
	}else{
	    info2("Reading WFS LLT OTF from %s\n", fnlotf);
	    intstat->lotf=ccellread("%s",fnlotf);
	    if(!intstat->lotf || !intstat->lotf->nx) error("Invalid lotf\n");
	}
    }
    if(parms->save.setup){//Save uplink PSF.
	int nwvl=intstat->lotf->nx;
	ccell*  lotf=intstat->lotf/*PCELL*/;
	int nlpsf=powfs[ipowfs].llt->pts->nx*parms->powfs[ipowfs].embfac;
	cmat *psfhat=cnew(nlpsf, nlpsf);
	dmat *psf=dnew(nlpsf, nlpsf);
	char header[64];
	zfarr *lltpsfsave=NULL;
	lltpsfsave=zfarr_init(nwvl, intstat->lotf->ny, "powfs%d_llt_psf", ipowfs);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    const double dx=powfs[ipowfs].llt->pts->dx;
	    const double wvl=parms->powfs[ipowfs].wvl->p[iwvl];
	    const double dpsf=wvl/(nlpsf*dx)*206265.;
	    snprintf(header, 64,"dtheta=%g #arcsecond\n", dpsf); 
	    free(psf->header); psf->header=strdup(header);	    
	    for(int illt=0; illt<intstat->lotf->ny; illt++){

		ccp(&psfhat, IND(lotf,iwvl,illt));
		cfftshift(psfhat);
		cfft2i(psfhat, 1);
		cfftshift(psfhat);
		creal2d(&psf, 0, psfhat, 1);
		info2("illt %d, iwvl %d has FWHM of %g\"\n",
		      illt, iwvl, sqrt(4.*(double)dfwhm(psf)/M_PI)*dpsf);
		zfarr_dmat(lltpsfsave, illt*nwvl+iwvl, psf);
	    }
	}
	zfarr_close(lltpsfsave);
	cfree(psfhat);
	dfree(psf);
    }
}

/**
   Createing subaperture short exposure PSF from the tip/tilt removed turbulence
   OTF and uplink OTF. Not including detector or elongation characteristics.  */
void gensepsf(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    const int nwvl=parms->powfs[ipowfs].nwvl;
    int nsa=powfs[ipowfs].saloc->nloc;
    int nllt;
    if(parms->powfs[ipowfs].llt)
	nllt=parms->powfs[ipowfs].llt->n;
    else
	nllt=0;
    int nlotf=0;
    if(nllt>0){
	nlotf=powfs[ipowfs].intstat->lotf->ny;
    }
    int notf=powfs[ipowfs].intstat->notf;
    powfs[ipowfs].intstat->nsepsf=notf>nlotf?notf:nlotf;
    assert(powfs[ipowfs].intstat->nsepsf==1 
	   || powfs[ipowfs].intstat->nsepsf==parms->powfs[ipowfs].nwfs);
    if(powfs[ipowfs].intstat->sepsf){
	cellfree(powfs[ipowfs].intstat->sepsf);
    }
    powfs[ipowfs].intstat->sepsf=dccellnew(powfs[ipowfs].intstat->nsepsf, 1);
    for(int isepsf=0; isepsf<powfs[ipowfs].intstat->nsepsf; isepsf++){
	int iotf=notf>1?isepsf:0;
	int ilotf=nlotf>1?isepsf:0;
	cmat **lotf=nlotf>0?(powfs[ipowfs].intstat->lotf->p+ilotf*nwvl):NULL;
	ccell* otf=powfs[ipowfs].intstat->otf->p[iotf]/*PCELL*/;
	powfs[ipowfs].intstat->sepsf->p[isepsf]=dcellnew(nsa,nwvl);
	dcell*  psepsf=powfs[ipowfs].intstat->sepsf->p[isepsf]/*PDELL*/;
	const double *area=powfs[ipowfs].realsaa->p[isepsf]->p;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    const int notfx=IND(otf,0,iwvl)->nx;
	    const int notfy=IND(otf,0,iwvl)->ny;
	    cmat *sepsf=cnew(notfx,notfy);
	    for(int isa=0; isa<nsa; isa++){
		double norm=area[isa]/((double)(notfx*notfy));
		ccp(&sepsf,IND(otf,isa,iwvl));/*peak in center */
		if(nllt>0){/*has laser launch */
		    if(sepsf->nx == lotf[iwvl]->nx){
			ccwm(sepsf,lotf[iwvl]);
		    }else{
			if(sepsf->nx > lotf[iwvl]->nx){
			    error("lotf(%ldx%ld) has smaller size than sepsf (%ldx%ld)\n", 
				  sepsf->nx, sepsf->nx, lotf[iwvl]->nx, lotf[iwvl]->nx);
			}
			cmat *tmp=cnew(sepsf->nx, sepsf->ny);
			cembedc(tmp, lotf[iwvl], 0,C_FULL);
			ccwm(sepsf, tmp);
			cfree(tmp);
		    }
		}
		cfftshift(sepsf); /*peak now in corner. */
		cfft2(sepsf,1);   /*turn to psf. FFT 1th */
		cfftshift(sepsf); /*psf with peak in center */
		creal2d(PIND(psepsf,isa,iwvl),0,sepsf,norm);/*copy to output. */
	    }
	    cfree(sepsf);
	}
    }
}
/**
   generate subaperture short exposure average pixel intensities sampled on
   detector from short expsoure PSF, the elongation transfer function of the
   sodium layer, and the detector transfer function. */
void gensei(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    const int ncompx=powfs[ipowfs].ncompx;
    const int ncompy=powfs[ipowfs].ncompy;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    const int nsa=powfs[ipowfs].saloc->nloc;
    const int nllt=parms->powfs[ipowfs].llt?parms->powfs[ipowfs].llt->n:0;
    /**
       ni0 may be greater than 1 in the following two cases
       1) multiple LLT
       2) different signal level or wvlwts
       3) powfs[ipowfs].bkgrnd contains rayleigh scatter bkgrnd for each wfs in this powfs.
    */
    const int nsepsf=intstat->nsepsf;
    int ni0=nllt<=1?nsepsf:nllt;
    if(ni0==1 && parms->powfs[ipowfs].nwfs>1){/*check wvlwts. */
	const int iwfs0=parms->powfs[ipowfs].wfs->p[0];
	const double siglev0=parms->wfs[iwfs0].siglev;
	const double *wvlwts0=parms->wfs[iwfs0].wvlwts->p;
	for(int jwfs=1; jwfs<parms->powfs[ipowfs].nwfs;jwfs++){
	    const int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    const double siglev=parms->wfs[iwfs].siglev;
	    const double *wvlwts=parms->wfs[iwfs].wvlwts->p;
	    if(fabs(siglev-siglev0)>EPS){
		ni0=parms->powfs[ipowfs].nwfs;
		warning("Different wfs for powfs %d has different siglev\n",ipowfs);
		goto done;
	    }
	    for(int iwvl=0; iwvl<parms->powfs[ipowfs].nwvl; iwvl++){
		if(fabs(wvlwts[iwvl]-wvlwts0[iwvl])>EPS){
		    ni0=parms->powfs[ipowfs].nwfs;
		    warning("Different wfs for powfs %d "
			    "has different wvlwts\n",ipowfs);
		    goto done;
		}
	    }
	}
    done:
	if(ni0!=1){
	    warning("Different wfs have different i0\n");
	}
    }
    if(powfs[ipowfs].bkgrnd && powfs[ipowfs].bkgrnd->ny>1){
	if(powfs[ipowfs].bkgrnd->ny != parms->powfs[ipowfs].nwfs){
	    error("powfs[%d].bkgrndfn must contain 1 or %d columns in this powfs\n",
		  ipowfs, parms->powfs[ipowfs].nwfs);
	    
	}
	ni0=parms->powfs[ipowfs].nwfs;
    }
    if(ni0>1){
	info2("number of i0 for matched filter is %d\n",ni0);
    }
    if(ni0!=1 && ni0!=parms->powfs[ipowfs].nwfs){
	error("Number of i0 must be either 1 or %d, but is %d\n",
	      parms->powfs[ipowfs].nwfs,ni0);
    }
    dcellfree(intstat->i0);
    dcellfree(intstat->gx);
    dcellfree(intstat->gy);
    cellfree(intstat->fotf);
    cellfree(intstat->potf);

    intstat->i0=dcellnew(nsa,ni0);
    intstat->gx=dcellnew(nsa,ni0);
    intstat->gy=dcellnew(nsa,ni0);
    if(parms->powfs[ipowfs].phytypesim==3 ){
	intstat->fotf=cccellnew(nsepsf, 1);
	for(int i=0; i<nsepsf; i++){
	    intstat->fotf->p[i]=ccellnew(nsa,nwvl);
	}
    }
    if(parms->dbg.wfslinearity!=-1 && parms->wfs[parms->dbg.wfslinearity].powfs==ipowfs){
	intstat->potf=cccellnew(nsepsf, 1);
	for(int i=0; i<nsepsf; i++){
	    intstat->potf->p[i]=ccellnew(nsa,nwvl);
	}
    }
    /* subaperture rotation angle. */
    dcell*  i0=intstat->i0/*PDELL*/;
    dcell*  gx=intstat->gx/*PDELL*/;
    dcell*  gy=intstat->gy/*PDELL*/;
  
    /*
      Notice, the generation of shifted i0s are not accurate
      because the PSF is not enough to cover the size.
      Disable the computation.
    */
    const int pixpsax=powfs[ipowfs].pixpsax;
    const int pixpsay=powfs[ipowfs].pixpsay;
  
    for(int ii0=0; ii0<ni0; ii0++){
	for(int isa=0; isa<nsa; isa++){
	    IND(i0,isa,ii0)=dnew(pixpsax,pixpsay);
	    IND(gx,isa,ii0)=dnew(pixpsax,pixpsay);
	    IND(gy,isa,ii0)=dnew(pixpsax,pixpsay);
	}
    }
  
    const int i0scale=parms->powfs[ipowfs].i0scale;
    if(i0scale){
	warning("i0 is scaled to match sa area\n");
    }
    const int isepsf_multiplier=nsepsf>1?1:0;
    const int irot_multiplier=nllt>1?1:0;
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	const int idtf_multiplier
	    =powfs[ipowfs].dtf[iwvl].si->ny>1?1:0;
	const int idtfisa_multiplier
	    =powfs[ipowfs].dtf[iwvl].si->nx>1?1:0;
	const dcomplex *Ux=powfs[ipowfs].dtf[iwvl].Ux->p;
	const dcomplex *Uy=powfs[ipowfs].dtf[iwvl].Uy->p;
	
	const double norm=1./(double)(ncompx*ncompy);
	const ccell *petf=NULL;
	void (*pccwm)(cmat*,const cmat*)=NULL;
	int rotpsf=0;
	int ietf_multiplier=0;
	if(nllt){
	    if(powfs[ipowfs].etfprep[iwvl].p1){
		petf=powfs[ipowfs].etfprep[iwvl].p1;
		pccwm=ccwmcol;
		rotpsf=1;
		if(powfs[ipowfs].etfprep[iwvl].p1->ny==1)
		    ietf_multiplier=0;
		else
		    ietf_multiplier=1;
	    }else{
		petf=powfs[ipowfs].etfprep[iwvl].p2;
		pccwm=ccwm;
		rotpsf=0;
		if(powfs[ipowfs].etfprep[iwvl].p2->ny==1)
		    ietf_multiplier=0;
		else
		    ietf_multiplier=1;
	    }
	}

	for(int ii0=0; ii0<ni0; ii0++){
	    double *area=powfs[ipowfs].realsaa->p[ii0]->p;
	    int isepsf=ii0*isepsf_multiplier;
	    int idtf=ii0*idtf_multiplier;
	    int irot=ii0*irot_multiplier;
	    int ietf=ii0*ietf_multiplier;
	    int iwfs=parms->powfs[ipowfs].wfs->p[ii0];
	    double wvlsig=parms->wfs[iwfs].wvlwts->p[iwvl]
		*parms->wfs[iwfs].siglev*parms->powfs[ipowfs].dtrat;
	    dcell*  psepsf=intstat->sepsf->p[isepsf]/*PDELL*/;
	    cmat **nominals=powfs[ipowfs].dtf[iwvl].fused?0:(powfs[ipowfs].dtf[iwvl].nominal->p+powfs[ipowfs].dtf[iwvl].nominal->nx*idtf);
	    dsp **sis=powfs[ipowfs].dtf[iwvl].si->p+powfs[ipowfs].dtf[iwvl].si->nx*idtf;
	    double *angles=nllt?(powfs[ipowfs].srot->p[irot]->p):0;
	    ccell *se_save=ccellnew(3, NTHREAD);
#ifdef _OPENMP
	    if(omp_in_parallel()){
		error("Already in parallel\n");
	    }
#endif
#pragma omp parallel shared(se_save) default(shared)
#pragma omp for 
	    for(int isa=0; isa<nsa; isa++){
		int ith=0;
		double angle=0;/*angle to rotate otf/psf */
		double angleg=0;/*angle to derivative of i0 to r/a from x/y */
		double anglegoff=0;

#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
#define seotfj se_save->p[3*ith]
#define seotfk se_save->p[3*ith+1]
#define sepsf se_save->p[3*ith+2]
		if(!seotfk){
		    seotfk=cnew(ncompx,ncompy);
		}
		cmat *nominal=NULL;
		dsp *si=NULL;

		int isadtf=isa*idtfisa_multiplier;
		if(nominals) nominal=nominals[isadtf];
		si=sis[isadtf];
		if(nllt && parms->powfs[ipowfs].radpix){
		    /*Polar CCD. */
		    if(rotpsf){/*OTF is along R/A direction. angleg is 0. */
			angle=angles[isa];
		    }else{
			angleg=angles[isa];
		    }
		}
		double pgrad[2];
		/*loaded psepsf. sum to 1 for full sa. peak in center */
		if(parms->powfs[ipowfs].mtchstc){
		    /*Forst psf to be centered. */
		    double pmax=dmax(IND(psepsf,isa,iwvl));
		    dcog(pgrad,IND(psepsf,isa,iwvl),0.5,0.5,0.1*pmax,0.2*pmax, 0);
		}
		ccpd(&sepsf,IND(psepsf,isa,iwvl));
		//cembedc(seotfk,sepsf,-angle,C_ABS);/*C_ABS causes sum of PSF to increase when there are negative values. Switch to literal copy.*/
		cembed(seotfk, sepsf, -angle);
		cfftshift(seotfk);/*PSF, peak in corner; */
		cfft2(seotfk,-1);/*turn to OTF, peak in corner, max is 1 */
		if(parms->powfs[ipowfs].mtchstc && fabs(pgrad[0])>EPS && fabs(pgrad[1])>EPS){
		    ctilt(seotfk,-pgrad[0],-pgrad[1],0);
		}
		if(nominal) ccwm(seotfk,nominal);
		cscale(seotfk, norm);/*normalized so that after fft, psf sum to 1*/
		if(intstat->potf){
		    ccp(&intstat->potf->p[isepsf]->p[iwvl*nsa+isa], seotfk);
		}
		if(nllt){/*elongation. */
		    (*pccwm)(seotfk,IND(petf,isa,ietf));
		}
		ccp(&seotfj,seotfk);/*backup */
		if(intstat->fotf){
		    ccp(&intstat->fotf->p[isepsf]->p[iwvl*nsa+isa], seotfk);
		}
		cfft2(seotfk,1);/*PSF with peak in center. sum to (pixtheta/dtheta)^2 due to nominal.*/
		/*no need fftshift becaose nominal is pre-treated */
		dspmulcreal(IND(i0,isa,ii0)->p,si,seotfk->p, wvlsig);
		ccp(&seotfk,seotfj);
		
		double ct=cos(angleg+anglegoff);
		double st=sin(angleg+anglegoff);
		
		for(int iy=0; iy<ncompy; iy++){
		    for(int ix=0; ix<ncompx; ix++){
			IND(seotfk,ix,iy)*=ct*Ux[ix]+st*Uy[iy];
			IND(seotfj,ix,iy)*=-st*Ux[ix]+ct*Uy[iy];
		    }
		}
		cfft2(seotfk,1);
		dspmulcreal(IND(gx,isa,ii0)->p,si,seotfk->p, wvlsig);
		cfft2(seotfj,1);
		dspmulcreal(IND(gy,isa,ii0)->p,si,seotfj->p, wvlsig);
		if(i0scale){
		    double scale=area[isa]/dsum(IND(i0,isa,ii0));
		    dscale(IND(i0,isa,ii0),scale);
		    dscale(IND(gx,isa,ii0),scale);
		    dscale(IND(gy,isa,ii0),scale);
		}
	    }/*for isa */
	    cellfree(se_save);
	}/*for ii0*/
    }/*iwvl */
}
