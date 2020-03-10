/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
    loc_t *loc=mksqloc_auto(powfs[ipowfs].pts->nx,
		       powfs[ipowfs].pts->nx,
		       powfs[ipowfs].pts->dx,
		       powfs[ipowfs].pts->dy);
    /*The embeding factor for embedding the aperture */
    const int embfac=parms->powfs[ipowfs].embfac;
    const int npsfx=powfs[ipowfs].pts->nx*embfac;    
    const int npsfy=powfs[ipowfs].pts->ny*embfac;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    int nsa=powfs[ipowfs].saloc->nloc;
 
    int notf=1;
    int has_ncpa=0;
    if(powfs[ipowfs].opdbias && parms->powfs[ipowfs].ncpa_method==2){
	//check whether opdbias is different between wfs[0] and following.
	int different=0;
	for(int iwfs=1; iwfs<parms->powfs[ipowfs].nwfs; iwfs++){
	    real diff=ddiff(powfs[ipowfs].opdbias->p[0], powfs[ipowfs].opdbias->p[iwfs]);
	    if(diff>1e-4){
		dbg("powfs[%d].opdbias[%d] is different from powfs[%d].opdbias[0] by %g.\n", 
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
    powfs[ipowfs].intstat->otf=cccellnew(notf, 1);
    for(int iotf=0; iotf<notf; iotf++){
	powfs[ipowfs].intstat->otf->p[iotf]=ccellnew(nsa,nwvl);
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	real wvl=parms->powfs[ipowfs].wvl->p[iwvl];
	for(int iotf=0; iotf<notf; iotf++){
	    dmat* opdbias=has_ncpa?powfs[ipowfs].opdbias->p[iotf]:NULL;
	    real thres=opdbias?1:(1-1e-10);
	    info2("There is %s bias\n", opdbias?"NCPA":"no");
	    //OTFs are always generated with native sampling. It is upsampled at gensepsf if necessary.
	    OMPTASK_SINGLE
		genotf(powfs[ipowfs].intstat->otf->p[iotf]->p+iwvl*nsa,
		       loc, powfs[ipowfs].realamp->p[iotf], opdbias, 
		       powfs[ipowfs].realsaa->p[iotf],
		       thres,wvl,NULL,parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0, 
		       npsfx, npsfy, nsa, 1);
	}
    }/*iwvl */
    locfree(loc);
}

void genseotf(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    char fnprefix[200]; fnprefix[0]='\0';
    uint32_t key=0;
    strcat(fnprefix, "SEOTF");
    if(powfs[ipowfs].amp_tel){
	for(int iamp=0; iamp<parms->powfs[ipowfs].nwfs; iamp++){
	    key=dhash(powfs[ipowfs].amp_tel->p[iamp], key);
	}
    }else{
	key=dhash(powfs[ipowfs].amp, key);
    }
    key=dhash(parms->powfs[ipowfs].wvl, key);
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
	strcat(fnprefix, tmp2);
    }
    char fnotf[PATH_MAX+20];
    char fnlock[PATH_MAX+40];
    snprintf(fnotf,PATH_MAX,"%s/SEOTF/",CACHE);
    if(!exist(fnotf)) {
	mymkdir("%s",fnotf);
    }
    long nsa=powfs[ipowfs].saloc->nloc;
    snprintf(fnotf,sizeof(fnotf),"%s/SEOTF/%s_D%g_%g_"
	     "r0_%g_L0%g_dsa%g_nsa%ld_dx1_%g_embfac%d_v2",
	     CACHE, fnprefix,
	     parms->aper.d,parms->aper.din, 
	     parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0, 
	     powfs[ipowfs].pts->dsa,nsa,
	     1./powfs[ipowfs].pts->dx, parms->powfs[ipowfs].embfac);
    snprintf(fnlock, sizeof(fnlock), "%s.lock", fnotf);
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    while(!intstat->otf){
	if(exist(fnlock) || !zfexist(fnotf)){/*need to create data */
	    int fd=lock_file(fnlock, 0, 0);/*nonblocking exclusive lock */
	    if(fd>=0){/*succeed */
		info2("Generating WFS OTF for %s...", fnotf);
		TIC;tic; genseotf_do(parms,powfs,ipowfs); toc("done");
		writebin(intstat->otf, "%s", fnotf);
	    }else{
		warning("Waiting for previous lock to release ...");
		fd=lock_file(fnlock, 1, 0);
	    }
	    close(fd); remove(fnlock);
	}else{
	    info2("Reading WFS OTF from %s\n", fnotf);
	    intstat->otf=cccellread("%s",fnotf);
	}
    }
}
/**
   Creating short exposure OTF caused by turbulence within LLT uplink aperture
*/
void genselotf_do(const PARMS_T *parms,POWFS_T *powfs,int ipowfs){
    if(!parms->powfs[ipowfs].llt) return;
    pts_t *lltpts=powfs[ipowfs].llt->pts;
    loc_t *loc=pts2loc(lltpts);
    const int nwvl=parms->powfs[ipowfs].nwvl;
    const int embfac=parms->powfs[ipowfs].embfac;
    const int npsfx=lltpts->nx*embfac;
    const int npsfy=lltpts->ny*embfac;
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
    //const int notfx=powfs[ipowfs].notfx;//keep LOTF and OTF same sampling
    //const int notfy=powfs[ipowfs].notfy;//keep LOTF and OTF same sampling
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	real wvl=parms->powfs[ipowfs].wvl->p[iwvl];
	real thres=1;
	for(int ilotf=0; ilotf<nlotf; ilotf++){
	    genotf(PP(lotf,iwvl,ilotf), loc, powfs[ipowfs].llt->amp, ncpa?ncpa->p[ilotf]:NULL, 
		   0, thres, wvl, NULL,parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0,
		   npsfx,npsfy, 1, 1);
	    if(!PP(lotf,iwvl,ilotf)){
		error("lotf is empty\n");
	    }
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
    key=dhash(parms->powfs[ipowfs].wvl, key);
    snprintf(fnprefix,80,"SELOTF_%0x",key);
    char fnlotf[PATH_MAX];
    snprintf(fnlotf,sizeof(fnlotf),"%s/SELOTF/", CACHE);
    if(!exist(fnlotf)){
	mymkdir("%s",fnlotf);
    }
    snprintf(fnlotf,sizeof(fnlotf),"%s/SELOTF/%s_"
	     "r0_%g_L0%g_lltd%g_dx1_%g_W%g_embfac%d_v2",
	     CACHE, fnprefix,
	     parms->powfs[ipowfs].r0, parms->powfs[ipowfs].L0, 
	     powfs[ipowfs].llt->pts->dsa,
	     1./powfs[ipowfs].llt->pts->dx,
	     parms->powfs[ipowfs].llt->widthp,parms->powfs[ipowfs].embfac);

    char fnlock[PATH_MAX+10];
    snprintf(fnlock, sizeof(fnlock), "%s.lock", fnlotf);
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
	cmat *psfhat=0;//cnew(nlpsf, nlpsf);
	dmat *psf=0;//dnew(nlpsf, nlpsf);
	char header[64];
	zfarr *lltpsfsave=NULL;
	lltpsfsave=zfarr_init(nwvl, intstat->lotf->ny, "powfs%d_llt_psf", ipowfs);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    const real dx=powfs[ipowfs].llt->pts->dx;
	    const real wvl=parms->powfs[ipowfs].wvl->p[iwvl];
	    const real dpsf=wvl/(nlpsf*dx)*206265.;
	    snprintf(header, 64,"dtheta=%g; #arcsecond\n", dpsf); 
	    for(int illt=0; illt<intstat->lotf->ny; illt++){
		ccp(&psfhat, P(lotf,iwvl,illt));
		cfftshift(psfhat);
		cfft2i(psfhat, 1);
		cfftshift(psfhat);
		creal2d(&psf, 0, psfhat, 1);
		info2("illt %d, iwvl %d has FWHM of %g\"\n",
		      illt, iwvl, sqrt(4.*(real)dfwhm(psf)/M_PI)*dpsf);
		free(psf->header); psf->header=strdup(header);
		zfarr_push(lltpsfsave, illt*nwvl+iwvl, psf);
	    }
	}
	zfarr_close(lltpsfsave);
	cfree(psfhat);
	dfree(psf);
    }
}
/**
   Upsample the otf in to out while preserving the PSF.
 */
static void upsample_otf(cmat *out, const cmat *in){
    if(in->nx==out->nx && in->ny==out->ny){
	ccp(&out, in);
    }else{
	cmat *temp=0;
	ccp(&temp, in);
	cfft2(temp, -1);
	cscale(temp, 1./(in->nx*in->ny));
	czero(out);
	ccpcorner(out, temp, C_FULL);
	cfft2(out, 1);
	cfree(temp);
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
    int notf=powfs[ipowfs].intstat->otf->nx;
    powfs[ipowfs].intstat->nsepsf=MAX(notf, nlotf);//notf>nlotf?notf:nlotf;
    assert(powfs[ipowfs].intstat->nsepsf==1 
	   || powfs[ipowfs].intstat->nsepsf==parms->powfs[ipowfs].nwfs);
    if(powfs[ipowfs].intstat->sepsf){
	cellfree(powfs[ipowfs].intstat->sepsf);
    }
    powfs[ipowfs].intstat->sepsf=dccellnew(powfs[ipowfs].intstat->nsepsf, 1);
    const int notfx=powfs[ipowfs].notfx;
    const int notfy=powfs[ipowfs].notfy;
    for(int isepsf=0; isepsf<powfs[ipowfs].intstat->nsepsf; isepsf++){
	//int iotf=notf>1?isepsf:0;
	//int ilotf=nlotf>1?isepsf:0;
	//cmat **lotf=nlotf>0?(powfs[ipowfs].intstat->lotf->p+ilotf*nwvl):NULL;
	ccell* otf=PR(powfs[ipowfs].intstat->otf, isepsf, 0);//->p[iotf];
	powfs[ipowfs].intstat->sepsf->p[isepsf]=dcellnew(nsa,nwvl);
	dcell*  psepsf=powfs[ipowfs].intstat->sepsf->p[isepsf];
	const real *area=powfs[ipowfs].realsaa->p[isepsf]->p;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cmat *sepsf=cnew(notfx, notfy);
	    cmat *lotf=0;
	    if(nllt>0){
		cmat *lotf2=PR(powfs[ipowfs].intstat->lotf, iwvl, isepsf);
		if(lotf2->nx!=notfx || lotf2->ny!=notfy){
		    lotf=cnew(notfx, notfy);
		    upsample_otf(lotf, lotf2);
		}else{
		    lotf=cref(lotf2);
		}
	    }
	    for(int isa=0; isa<nsa; isa++){
		real norm=area[isa]/((real)(notfx*notfy));
		if(P(otf,isa,iwvl)){
		    upsample_otf(sepsf, P(otf,isa,iwvl));/*peak in center */
		}else{
		    czero(sepsf);
		}
		if(lotf){
		    ccwm(sepsf, lotf);
		}
		cfftshift(sepsf); /*peak now in corner. */
		cfft2(sepsf,1);   /*turn to psf. FFT 1th */
		cfftshift(sepsf); /*psf with peak in center */
		creal2d(PP(psepsf,isa,iwvl),0,sepsf,norm);/*copy to output. */
	    }
	    cfree(sepsf);
	    cfree(lotf);
	}
    }
}
/**
   generate subaperture short exposure average pixel intensities sampled on
   detector from short expsoure PSF, the elongation transfer function of the
   sodium layer, and the detector transfer function. */
void gensei(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    const int notfx=powfs[ipowfs].notfx;
    const int notfy=powfs[ipowfs].notfy;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    const int nsa=powfs[ipowfs].saloc->nloc;
    const int nllt=parms->powfs[ipowfs].llt?parms->powfs[ipowfs].llt->n:0;
    const int radgx=parms->powfs[ipowfs].radgx;
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
	const real siglev0=parms->wfs[iwfs0].siglev;
	const real *wvlwts0=parms->wfs[iwfs0].wvlwts->p;
	for(int jwfs=1; jwfs<parms->powfs[ipowfs].nwfs;jwfs++){
	    const int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
	    const real siglev=parms->wfs[iwfs].siglev;
	    const real *wvlwts=parms->wfs[iwfs].wvlwts->p;
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
    if(parms->powfs[ipowfs].phytype_sim==3 ){
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
	    P(i0,isa,ii0)=dnew(pixpsax,pixpsay);
	    P(gx,isa,ii0)=dnew(pixpsax,pixpsay);
	    P(gy,isa,ii0)=dnew(pixpsax,pixpsay);
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
	const comp *Ux=powfs[ipowfs].dtf[iwvl].Ux->p;
	const comp *Uy=powfs[ipowfs].dtf[iwvl].Uy->p;
	
	const real norm=1./(real)(notfx*notfy);
	const ccell *petf=NULL;

	int ietf_multiplier=0;
	if(nllt){
	    petf=powfs[ipowfs].etfprep[iwvl].etf;

	    if(powfs[ipowfs].etfprep[iwvl].etf->ny==1)
		ietf_multiplier=0;
	    else
		ietf_multiplier=1;
	}

	for(int ii0=0; ii0<ni0; ii0++){
	    real *area=powfs[ipowfs].realsaa->p[ii0]->p;
	    int isepsf=ii0*isepsf_multiplier;
	    int idtf=ii0*idtf_multiplier;
	    int irot=ii0*irot_multiplier;
	    int ietf=ii0*ietf_multiplier;
	    int iwfs=parms->powfs[ipowfs].wfs->p[ii0];
	    real wvlsig=parms->wfs[iwfs].wvlwts->p[iwvl]
		*parms->wfs[iwfs].siglev*parms->powfs[ipowfs].dtrat;
	    dcell*  psepsf=intstat->sepsf->p[isepsf]/*PDELL*/;
	    cmat **nominals=powfs[ipowfs].dtf[iwvl].fused?0:(powfs[ipowfs].dtf[iwvl].nominal->p+powfs[ipowfs].dtf[iwvl].nominal->nx*idtf);
	    dsp **sis=powfs[ipowfs].dtf[iwvl].si->p+powfs[ipowfs].dtf[iwvl].si->nx*idtf;
	    real *angles=nllt?(powfs[ipowfs].srot->p[irot]->p):0;
	    ccell *se_save=ccellnew(3, NTHREAD);
#ifdef _OPENMP
	    if(omp_in_parallel()){
		warning("Already in parallel\n");
	    }
#endif
#pragma omp parallel default(shared)
#pragma omp for 
	    for(int isa=0; isa<nsa; isa++){
		int ith=0;
		real angleg=0;/*angle to derivative of i0 to r/a from x/y */

#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
#define seotfj se_save->p[3*ith]
#define seotfk se_save->p[3*ith+1]
//#define sepsf se_save->p[3*ith+2]
		if(!seotfk){
		    seotfk=cnew(notfx,notfy);
		}
		cmat *nominal=NULL;
		dsp *si=NULL;

		int isadtf=isa*idtfisa_multiplier;
		if(nominals) nominal=nominals[isadtf];
		si=sis[isadtf];
		if(nllt && parms->powfs[ipowfs].radpix){
		    if(radgx){
			angleg=angles[isa];
		    }
		}
		real pgrad[2];
		/*loaded psepsf. sum to 1 for full sa. peak in center */
		if(parms->powfs[ipowfs].mtchstc){
		    /*Forst psf to be centered. */
		    real pmax=dmax(P(psepsf,isa,iwvl));
		    dcog(pgrad,P(psepsf,isa,iwvl),0.5,0.5,0.1*pmax,0.2*pmax, 0);
		}
		if(dsum(P(psepsf,isa,iwvl))>1.1){
		    error("Short exposure PSF has wrong scaling. It should total to <=1\n");
		}
		/*C_ABS causes sum of PSF to increase when there are negative values. Switch to literal copy.*/
		cembedd(seotfk, P(psepsf,isa,iwvl), 0);
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
		    ccwm(seotfk,P(petf,isa,ietf));
		}
		ccp(&seotfj,seotfk);/*backup */
		if(intstat->fotf){
		    ccp(&intstat->fotf->p[isepsf]->p[iwvl*nsa+isa], seotfk);
		}
		cfft2(seotfk,1);/*PSF with peak in center. sum to (pixtheta/dtheta)^2 due to nominal.*/
		/*no need fftshift becaose nominal is pre-treated */
		dspmulcreal(P(i0,isa,ii0)->p,si,seotfk->p, wvlsig);
		ccp(&seotfk,seotfj);
		
		real ct=cos(angleg);
		real st=sin(angleg);
		
		for(int iy=0; iy<notfy; iy++){
		    for(int ix=0; ix<notfx; ix++){
			P(seotfk,ix,iy)*=ct*Ux[ix]+st*Uy[iy];
			P(seotfj,ix,iy)*=-st*Ux[ix]+ct*Uy[iy];
		    }
		}
		cfft2(seotfk,1);
		dspmulcreal(P(gx,isa,ii0)->p,si,seotfk->p, wvlsig);
		cfft2(seotfj,1);
		dspmulcreal(P(gy,isa,ii0)->p,si,seotfj->p, wvlsig);
		if(i0scale){
		    real scale=area[isa]/dsum(P(i0,isa,ii0));
		    dscale(P(i0,isa,ii0),scale);
		    dscale(P(gx,isa,ii0),scale);
		    dscale(P(gy,isa,ii0),scale);
		}
	    }/*for isa */
	    cellfree(se_save);
	}/*for ii0*/
    }/*iwvl */
}

/**
   The routine used to generate matched filter from WFS mean short exposure
   pixel intensities.
 */
void genmtch(const PARMS_T *parms, POWFS_T *powfs, const int ipowfs){
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    const real pixthetax=parms->powfs[ipowfs].radpixtheta;
    const real pixthetay=parms->powfs[ipowfs].pixtheta;
    const real rne=parms->powfs[ipowfs].rne;
    const real bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
    const real bkgrndc=bkgrnd*parms->powfs[ipowfs].bkgrndc;
    int ni0=intstat->i0->ny;
    if(ni0!=1 && ni0!=parms->powfs[ipowfs].nwfs){
	error("ni0 should be either 1 or %d\n", parms->powfs[ipowfs].nwfs);
    }
    const int nsa=powfs[ipowfs].saloc->nloc;
    //Prevent printing of NEA during recomputing of matched filter
    const int print_nea=intstat->mtche?0:1;
    dcellfree(intstat->mtche);
    dfree(intstat->i0sum);

    dcellfree(powfs[ipowfs].sanea);
    dcell *sanea=powfs[ipowfs].sanea=dcellnew(ni0,1);
    intstat->i0sum=dnew(nsa,ni0);
    intstat->i0sumsum=dnew(ni0, 1);

    const dcell *i0s=intstat->i0;
    const dcell* gxs=parms->powfs[ipowfs].mtchfft?0:intstat->gx/*PDELL*/;
    const dcell* gys=parms->powfs[ipowfs].mtchfft?0:intstat->gy/*PDELL*/;
    dmat *i0sum=intstat->i0sum;
    long npix=powfs[ipowfs].pixpsax*powfs[ipowfs].pixpsay;
    dcell *mtche=intstat->mtche=dcellnew_same(nsa,ni0,2,npix);
  
    //dcell *saneaxy=powfs[ipowfs].saneaxy;
    int nllt;
    if(parms->powfs[ipowfs].llt){
	nllt=parms->powfs[ipowfs].llt->n;
    }else{
	nllt=0;
    }
    int irot_multiplier=nllt>1?1:0;
    const int mtchadp=parms->powfs[ipowfs].mtchadp;
    real sigratio=parms->powfs[ipowfs].sigrecon>0?(parms->powfs[ipowfs].sigrecon/parms->powfs[ipowfs].siglev):1;
    real sigratior=1./sigratio;
    if(sigratio<0){
	error("sigratio cannot be negative\n");
    }

    for(int ii0=0; ii0<ni0; ii0++){
	int iwfs=parms->powfs[ipowfs].wfs->p[ii0];
	const real siglev=parms->powfs[ipowfs].dtrat*parms->wfs[iwfs].siglev;
	real i0thres=MAX(0.1*siglev, rne*10);
	real nea2thres=pixthetax*pixthetay*100;
	real *srot=NULL;
	if(powfs[ipowfs].srot){
	    int irot=ii0*irot_multiplier;
	    srot=powfs[ipowfs].srot->p[irot]->p;
	}
	sanea->p[ii0]=dnew(nsa,3);
	dmat*  psanea=sanea->p[ii0]/*PDMAT*/;
	real i0sumsum=0;
	int crdisable=0;/*adaptively disable mtched filter based in FWHM. */
	int ncrdisable=0;
	const int radgx=parms->powfs[ipowfs].radgx;
	dmat *nea2=0;
	for(int isa=0; isa<nsa; isa++){
	    real pixrot=0;//pixel rotation
	    if(srot && parms->powfs[ipowfs].radpix){
		pixrot=srot[isa]; 
	    }
	    if(mtchadp){
		long fwhm=dfwhm(P(i0s,isa,ii0));
		if(fwhm>4){
		    crdisable=0;
		}else{
		    crdisable=1;
		    ncrdisable++;
		}
	    }
	    dmat* bkgrnd2=NULL;
	    dmat* bkgrnd2c=NULL;
	    if(powfs[ipowfs].bkgrnd){
		bkgrnd2=powfs[ipowfs].bkgrnd->p[ii0*nsa+isa]; 
	    }
	    if(powfs[ipowfs].bkgrndc){
		bkgrnd2c=powfs[ipowfs].bkgrndc->p[ii0*nsa+isa]; 
	    }
	    mtch(PP(mtche,isa,ii0),&nea2, P(i0s,isa,ii0),
		 gxs?P(gxs,isa,ii0):0, gys?P(gys,isa,ii0):0, 
		 parms->powfs[ipowfs].qe,
		 bkgrnd2, bkgrnd2c, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
		 pixrot, radgx, crdisable?0:parms->powfs[ipowfs].mtchcr);
	    if(fabs(sigratio-1)>1e-5){
		dscale(P(i0s,isa,ii0), sigratio);
		if(gxs) dscale(P(gxs,isa,ii0), sigratio);
		if(gys) dscale(P(gys,isa,ii0), sigratio);
		mtch(NULL,&nea2, P(i0s,isa,ii0),
		     gxs?P(gxs,isa,ii0):0, gys?P(gys,isa,ii0):0, 
		     parms->powfs[ipowfs].qe,
		     bkgrnd2, bkgrnd2c, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
		     pixrot, radgx, crdisable?0:parms->powfs[ipowfs].mtchcr);	
		dscale(P(i0s,isa,ii0), sigratior);
		if(gxs) dscale(P(gxs,isa,ii0), sigratior);
		if(gys) dscale(P(gys,isa,ii0), sigratior);
	    }
	    P(i0sum,isa,ii0)=dsum(P(i0s,isa,ii0));
	    i0sumsum+=P(i0sum,isa,ii0);

	    if(P(i0sum,isa,ii0)<i0thres || nea2->p[0]>nea2thres || nea2->p[3]>nea2thres){
		//Signal level too low or error to high.
		nea2->p[0]=nea2->p[3]=nea2thres;
		nea2->p[1]=nea2->p[2]=0;
		dset(P(mtche,isa,ii0), 0);
	    }
	    if(parms->powfs[ipowfs].mtchcpl==0 
	       && (!parms->powfs[ipowfs].radpix || parms->powfs[ipowfs].radgx)){
		/*remove coupling between r/a (x/y) measurements. */
		nea2->p[1]=nea2->p[2]=0;
	    }
	    P(psanea,isa,0)=nea2->p[0];
	    P(psanea,isa,1)=nea2->p[3];
	    P(psanea,isa,2)=nea2->p[1];
	}/*isa  */
	dfree(nea2);

	if(mtchadp){
	    info2("Mtched filter contraint are disabled for %d subaps out of %d.\n",
		  ncrdisable, nsa);
	}
	intstat->i0sumsum->p[ii0]=i0sumsum;
    }/*ii0 */
    if(print_nea){
	info("Matched filter sanea:\n");
	if(powfs[ipowfs].sprint){/*print nea for select subapertures.*/
	    for(int ii0=0; ii0<ni0; ii0++){
		int illt=0;
		if(ni0==parms->powfs[ipowfs].llt->n){
		    illt=ii0;
		}else if(ni0==parms->powfs[ipowfs].nwfs && parms->powfs[ipowfs].llt->n==1){
		    illt=0;
		}else{
		    error("Invalid combination\n");
		}
		info("ii0 %d, llt %d.\n", ii0, illt);
		info("sa index   dist   noise equivalent angle\n");
		dmat*  psanea=sanea->p[ii0]/*PDMAT*/;
		for(int ksa=0; ksa<powfs[ipowfs].sprint->p[illt]->nx; ksa++){
		    int isa=(int)powfs[ipowfs].sprint->p[illt]->p[ksa];
		    if(isa>0){
			info("sa %4d: %5.1f m, (%6.2f, %6.2f) mas\n", 
			      isa, powfs[ipowfs].srsa->p[illt]->p[isa], 
			      sqrt(P(psanea,isa,0))*206265000,
			      sqrt(P(psanea,isa,1))*206265000);
		    }
		}
	    }
	}else{
	    real dsa=powfs[ipowfs].saloc->dx;
	    real llimit=-dsa/2;
	    real ulimit=dsa/2;
	    info("index: position noise equivalent angle\n");
	    for(int isa=0; isa<nsa; isa++){
		real locx=powfs[ipowfs].saloc->locx[isa];
		real locy=powfs[ipowfs].saloc->locy[isa];
		if((parms->powfs[ipowfs].llt && (nsa<10 || (locx>0&&locy>llimit&&locy<ulimit)))
		   ||(!parms->powfs[ipowfs].llt && locx>=0 && locx<dsa*0.6 && locy>=0 && locy<dsa*0.6)
		    || nsa<=4){
		    info("sa%4d:%6.1fm",isa, locx);
		    for(int ii0=0; ii0<ni0; ii0++){
			info(" (%4.1f,%4.1f)", 
			      sqrt(P(sanea->p[ii0],isa,0))*206265000,
			      sqrt(P(sanea->p[ii0],isa,1))*206265000);
		    }//for ii0
		    info(" mas\n");
		}
	    }/*isa  */
	}
    }
    if(parms->save.setup){
	writebin(sanea, "powfs%d_sanea", ipowfs);
    }
    if(parms->powfs[ipowfs].phytype_recon==1 && parms->recon.glao && ni0>0){
	info2("Averaging saneaxy of different WFS for GLAO mode\n");
	dmat *sanea2=0;
	real scale=1./ni0;
	for(int ii0=0; ii0<ni0; ii0++){
	    dadd(&sanea2, 1, sanea->p[ii0], scale);
	}
	dcellfree(powfs[ipowfs].sanea);
	powfs[ipowfs].sanea=dcellnew(1,1);
	powfs[ipowfs].sanea->p[0]=sanea2;
    }
}
