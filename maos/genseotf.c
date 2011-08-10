/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "maos.h"
#include "genseotf.h"
#define DEBUG_POWFSINTS 0
#ifndef ROT_OTF
#define ROT_OTF 0
#endif
#define SAVE_OTF 0

/**
   Master routine that generates short exposure OTF by calling genotf() in the
   library with p/t/t removal set.
*/
void genseotf(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    //create a grid representing the aperture.
    loc_t *loc=mksqloc(powfs[ipowfs].pts->nx,
		       powfs[ipowfs].pts->nx,
		       powfs[ipowfs].pts->dx,
		       0,0);
    //size of the OTF grid
    int ncompx=powfs[ipowfs].pts->nx*parms->powfs[ipowfs].embfac;
    int ncompy=ncompx;
    //The embeding factor for embedding the aperture
    const int embfac=parms->powfs[ipowfs].embfac;
    const double dxsa=powfs[ipowfs].pts->dsa;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    int nsa=powfs[ipowfs].pts->nsa;
 
    int notf=1;
    int has_ncpa=0;
    if(parms->powfs[ipowfs].ncpa && parms->powfs[ipowfs].ncpa_method==2){
	notf=parms->powfs[ipowfs].nwfs;
	has_ncpa=1;
    }else if(powfs[ipowfs].nlocm){
	notf=MAX(notf,powfs[ipowfs].nlocm);
    }else{
	notf=1;
    }
    powfs[ipowfs].intstat->notf=notf;
    powfs[ipowfs].intstat->otf=calloc(notf, sizeof(ccell*));
    for(int iotf=0; iotf<notf; iotf++){
	powfs[ipowfs].intstat->otf[iotf]=ccellnew(nsa,nwvl);
    }
 
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	double wvl=parms->powfs[ipowfs].wvl[iwvl];
	double dtheta=wvl/(dxsa*embfac);
	for(int iotf=0; iotf<notf; iotf++){
	    double *opdbias=has_ncpa?powfs[ipowfs].ncpa->p[iotf]->p:NULL;
	    double thres=opdbias?1:1-1e-10;
	    info2("There is %s bias\n", opdbias?"NCPA":"no");
	    genotf(powfs[ipowfs].intstat->otf[iotf]->p+iwvl*nsa,
		   loc, powfs[ipowfs].realamp[iotf], opdbias, 
		   powfs[ipowfs].realsaa[iotf],
		   thres,wvl,dtheta,NULL,parms->atm.r0, parms->atm.l0, 
		   ncompx, ncompy, nsa, 1);
	}
    }//iwvl
    locfree(loc);
}
/**
   Creating short exposure OTF caused by turbulence within LLT uplink aperture
*/
void genselotf(const PARMS_T *parms,POWFS_T *powfs,int ipowfs){
    if(!parms->powfs[ipowfs].llt) return;
    loc_t *loc=pts2loc(powfs[ipowfs].llt->pts);
    int ncompx=powfs[ipowfs].llt->pts->nx*parms->powfs[ipowfs].embfac;
    int ncompy=ncompx;
    const int nwvl=parms->powfs[ipowfs].nwvl;

    int nlotf=1;
    dcell *ncpa=powfs[ipowfs].llt->ncpa;
    if(ncpa){
	nlotf=ncpa->nx*ncpa->ny;
    }
    powfs[ipowfs].intstat->lotf=ccellnew(nwvl,nlotf);
    PCCELL(powfs[ipowfs].intstat->lotf, lotf);
    if(nwvl!=1){
	warning("LGS has multi-color!\n");
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	double wvl=parms->powfs[ipowfs].wvl[iwvl];
	double dtheta=powfs[ipowfs].dtheta->p[iwvl];
	double thres=1;
	double one=1;
	for(int ilotf=0; ilotf<nlotf; ilotf++){
	    genotf(&lotf[ilotf][iwvl], loc, powfs[ipowfs].llt->amp->p, ncpa?ncpa->p[ilotf]->p:NULL, 
		   &one, thres, wvl, dtheta, NULL,parms->atm.r0, parms->atm.l0, ncompx, ncompy, 1, 1);
	}
    }//iwvl
    locfree(loc);
}
/**
   Createing subaperture short exposure PSF from the tip/tilt removed turbulence
   OTF and uplink OTF. Not including detector or elongation characteristics.  */
void gensepsf(const PARMS_T *parms, POWFS_T *powfs, int ipowfs){
    const int nwvl=parms->powfs[ipowfs].nwvl;
    int nsa=powfs[ipowfs].pts->nsa;
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
    powfs[ipowfs].intstat->sepsf=calloc(powfs[ipowfs].intstat->nsepsf, sizeof(dcell*));
    for(int isepsf=0; isepsf<powfs[ipowfs].intstat->nsepsf; isepsf++){
	int iotf=notf>1?isepsf:0;
	int ilotf=nlotf>1?isepsf:0;
	cmat **lotf=nlotf>0?powfs[ipowfs].intstat->lotf->p+ilotf*nwvl:NULL;
	PCCELL(powfs[ipowfs].intstat->otf[iotf],otf);
	powfs[ipowfs].intstat->sepsf[isepsf]=dcellnew(nsa,nwvl);
	dmat *(*psepsf)[nsa]=(void*)powfs[ipowfs].intstat->sepsf[isepsf]->p;
	const double *area=powfs[ipowfs].realsaa[isepsf];
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    const int notfx=otf[iwvl][0]->nx;
	    const int notfy=otf[iwvl][0]->ny;
	    cmat *sepsf=cnew(notfx,notfy);
	    cfft2plan(sepsf,1);
	    
	    for(int isa=0; isa<nsa; isa++){
		double norm=area[isa]/((double)(notfx*notfy));
		ccp(&sepsf,otf[iwvl][isa]);//peak in center
		if(nllt>0){//has laser launch
		    if(sepsf->nx == lotf[iwvl]->nx){
			ccwm(sepsf,lotf[iwvl]);
		    }else{
			assert(sepsf->nx < lotf[iwvl]->nx);
			cmat *tmp=cnew(sepsf->nx, sepsf->ny);
			cembed(tmp, lotf[iwvl], 0,C_FULL);
			ccwm(sepsf, tmp);
			cfree(tmp);
		    }
		}
		cfftshift(sepsf); //peak now in corner.
		cfft2(sepsf,1);   //turn to psf. FFT 1th
		cfftshift(sepsf); //psf with peak in center
		creal2d(&psepsf[iwvl][isa],0,sepsf,norm);//copy to output.
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
    if(parms->powfs[ipowfs].radrot){
#if ROT_OTF == 1
    info2("Rotating OTF for Polar CCD\n");//Used mainly for off-axis launch
#elif ROT_OTF==0
    info2("Rotating PSF for Polar CCD\n");//Used mainly for on-axis launch
#else
    error("Invalid ROT_OTF\n");
#endif
    }
    const int ncompx=powfs[ipowfs].ncompx;
    const int ncompy=powfs[ipowfs].ncompy;
    const int nwvl=parms->powfs[ipowfs].nwvl;
    int nsa=powfs[ipowfs].pts->nsa;
    int nllt;
    if(parms->powfs[ipowfs].llt){
	nllt=parms->powfs[ipowfs].llt->n;
    }else{
	nllt=0;
    }
    /**
       ni0 may be greater than 1 in the following two cases
       1) multiple LLT
       2) different signal level or wvlwts
       3) powfs[ipowfs].bkgrnd contains rayleigh scatter bkgrnd for each wfs in this powfs.
    */
    int ni0;
    if(nllt<=1){ 
	ni0=powfs[ipowfs].intstat->nsepsf;
    }else{
	ni0=nllt;
    }
    if(ni0==1 && parms->powfs[ipowfs].nwfs>1){//check wvlwts.
	int iwfs0=parms->powfs[ipowfs].wfs[0];
	double siglev0=parms->wfs[iwfs0].siglev;
	double *wvlwts0=parms->wfs[iwfs0].wvlwts;
	for(int jwfs=1; jwfs<parms->powfs[ipowfs].nwfs;jwfs++){
	    int iwfs=parms->powfs[ipowfs].wfs[jwfs];
	    double siglev=parms->wfs[iwfs].siglev;
	    double *wvlwts=parms->wfs[iwfs].wvlwts;
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
    info2("number of i0 for matched filter is %d\n",ni0);
    if(ni0!=1 && ni0!=parms->powfs[ipowfs].nwfs){
	error("Number of i0 must be either 1 or %d, but is %d\n",
	      parms->powfs[ipowfs].nwfs,ni0);
    }
    powfs[ipowfs].intstat->i0=dcellnew(nsa,ni0);
    powfs[ipowfs].intstat->gx=dcellnew(nsa,ni0);
    powfs[ipowfs].intstat->gy=dcellnew(nsa,ni0);
    //subaperture rotation angle.
    dmat* (*i0)[nsa]=(dmat*(*)[nsa])powfs[ipowfs].intstat->i0->p;
    dmat* (*gx)[nsa]=(dmat*(*)[nsa])powfs[ipowfs].intstat->gx->p;
    dmat* (*gy)[nsa]=(dmat*(*)[nsa])powfs[ipowfs].intstat->gy->p;
    /*
      Notice, the generation of shifted i0s are not accurate
      because the PSF is not enough to cover the size.
      Disable the computation.
    */
    const int pixpsax=powfs[ipowfs].pixpsax;
    const int pixpsay=powfs[ipowfs].pixpsay;
  
    for(int ii0=0; ii0<ni0; ii0++){
	for(int isa=0; isa<nsa; isa++){
	    i0[ii0][isa]=dnew(pixpsax,pixpsay);
	    gx[ii0][isa]=dnew(pixpsax,pixpsay);
	    gy[ii0][isa]=dnew(pixpsax,pixpsay);
	}
    }
  
    const int i0scale=parms->powfs[ipowfs].i0scale;
    if(i0scale){
	warning("i0 is scaled to match sa area\n");
    }
#if SAVE_OTF==1
    ccell *saveotf=ccellnew(nsa,nwvl);
    cmat *(*psaveotf)[nsa]=(void*)saveotf->p;
    ccell *saveotfetf=ccellnew(nsa,nwvl);
    cmat *(*psaveotfetf)[nsa]=(void*)saveotfetf->p;
    ccell *savepsf=ccellnew(nsa,nwvl);
    cmat *(*psavepsf)[nsa]=(void*)savepsf->p;
#endif
    int isepsf_multiplier
	=powfs[ipowfs].intstat->nsepsf>1?1:0;
    int irot_multiplier=nllt>1?1:0;
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	int idtf_multiplier
	    =powfs[ipowfs].dtf[iwvl].si->ny>1?1:0;
	int idtfisa_multiplier
	    =powfs[ipowfs].dtf[iwvl].si->nx>1?1:0;
	cmat *seotfj=cnew(ncompx,ncompy);
	cmat *seotfk=cnew(ncompx,ncompy);
	cfft2plan(seotfk,-1);
	cfft2plan(seotfk,1);
	cfft2plan(seotfj,1);
	dcomplex *Ux=powfs[ipowfs].dtf[iwvl].Ux->p;
	dcomplex *Uy=powfs[ipowfs].dtf[iwvl].Uy->p;
	
	double norm=1./(double)(ncompx*ncompy);
	const int npsf=powfs[ipowfs].intstat->sepsf[0]->p[0]->nx;
#if ROT_OTF == 1
	double xscale=(double)npsf/(double)ncompx;
	double yscale=(double)npsf/(double)ncompy;
#endif
	cmat *sepsf=cnew(npsf,npsf);
	cfft2plan(sepsf,-1);

	const cmat *(*petf)[nsa]=NULL;
	void (*pccwm)(cmat*,const cmat*)=NULL;
	int rotpsf=0;
	int ietf_multiplier=0;
	if(nllt){
	    if(powfs[ipowfs].etfprep[iwvl].p1){
		petf=(void*)powfs[ipowfs].etfprep[iwvl].p1->p;
		pccwm=ccwmcol;
		rotpsf=1;
		if(powfs[ipowfs].etfprep[iwvl].p1->ny==1)
		    ietf_multiplier=0;
		else
		    ietf_multiplier=1;
	    }else{
		petf=(void*)powfs[ipowfs].etfprep[iwvl].p2->p;
		pccwm=ccwm;
		rotpsf=0;
		if(powfs[ipowfs].etfprep[iwvl].p2->ny==1)
		    ietf_multiplier=0;
		else
		    ietf_multiplier=1;
	    }
	}
	cmat *nominal=NULL;
	dsp *si=NULL;

	double angle=0;//angle to rotate otf/psf
	double angleg=0;//angle to derivative of i0 to r/a from x/y
	double anglegoff=0;
	for(int ii0=0; ii0<ni0; ii0++){
	    const double *area=powfs[ipowfs].realsaa[ii0];
	    int isepsf=ii0*isepsf_multiplier;
	    int idtf=ii0*idtf_multiplier;
	    int irot=ii0*irot_multiplier;
	    int ietf=ii0*ietf_multiplier;
	    int iwfs=parms->powfs[ipowfs].wfs[ii0];
	    double wvlsig=parms->wfs[iwfs].wvlwts[iwvl]
		*parms->wfs[iwfs].siglev*parms->powfs[ipowfs].dtrat;
	    info2("iwvl=%d, iwfs=%d, wvlsig=%g\n",iwvl,iwfs,wvlsig);
	    dmat *(*psepsf)[nsa]=(void*)powfs[ipowfs].intstat->sepsf[isepsf]->p;
	    double pgrad[2];
	    cmat **nominals=NULL;
	    if(!powfs[ipowfs].dtf[iwvl].fused){//may be null if fused to etf
		nominals=powfs[ipowfs].dtf[iwvl].nominal->p+nsa*idtf;
	    }
	    dsp **sis=powfs[ipowfs].dtf[iwvl].si->p+nsa*idtf;
	    const double *angles=NULL;
	    if(nllt)
		angles=powfs[ipowfs].srot->p[irot]->p;
	    for(int isa=0; isa<nsa; isa++){
		int isadtf=isa*idtfisa_multiplier;
		if(nominals)
		    nominal=nominals[isadtf];
		si=sis[isadtf];
		if(nllt && parms->powfs[ipowfs].radpix){
		    //Polar CCD.
		    if(rotpsf){//OTF is along R/A direction. angleg is 0.
			angle=angles[isa];
		    }else{
			angleg=angles[isa];
		    }
		}
	
		//loaded psepsf. sum to 1 for full sa. peak in center
		if(parms->powfs[ipowfs].mtchstc){
		    //Forst psf to be centered.
		    double pmax=dmax(psepsf[iwvl][isa]);
		    dcog(pgrad,psepsf[iwvl][isa],0.5,0.5,0.1*pmax,0.2*pmax);
		}
		ccpd(&sepsf,psepsf[iwvl][isa]);
#if ROT_OTF == 1//ROTATE OTF
		cfftshift(sepsf);//peak in corner
		cfft2(sepsf,-1);//turn to OTF;//max at 1.
		cfftshift(sepsf);//peak in center
		cembedscaleout(seotfk,sepsf,xscale,yscale,-angle,C_FULL);
		cfftshift(seotfk); //otf, peak in corner
#endif

#if ROT_OTF==0//rotate PSF
		cembed(seotfk,sepsf,-angle,C_ABS);//ABS to avoid small negative
#if SAVE_OTF == 1
		ccp(&psavepsf[iwvl][isa],seotfk);
#endif
		cfftshift(seotfk);//PSF, peak in corner;
		cfft2(seotfk,-1);//turn to OTF peak in corner
		if(parms->powfs[ipowfs].mtchstc && fabs(pgrad[0])>EPS && fabs(pgrad[1])>EPS){
		    ctilt(seotfk,-pgrad[0],-pgrad[1],0);
		}
#endif
#if SAVE_OTF ==1
		ccp(&psaveotf[iwvl][isa],seotfk);
#endif
		if(nllt){//elongation.
		    (*pccwm)(seotfk,petf[ietf][isa]);
		}
		//seotfk has peak in corner
#if SAVE_OTF ==1
		ccp(&psaveotfetf[iwvl][isa],seotfk);
#endif
		ccwm2(seotfk,nominal,norm);//NULL is handled correctly.
		ccp(&seotfj,seotfk);//backup
		cfft2(seotfk,1);//peak in center.
		//no need fftshift becaose nominal is pre-treated
		spmulcreal(i0[ii0][isa]->p,si,seotfk->p, wvlsig);
		ccp(&seotfk,seotfj);
		dcomplex(*X)[ncompx]
		    =(dcomplex(*)[ncompx])(seotfk->p);
		dcomplex(*Y)[ncompx]
		    =(dcomplex(*)[ncompx])(seotfj->p);
		
		double ct=cos(angleg+anglegoff);
		double st=sin(angleg+anglegoff);
		
		for(int iy=0; iy<ncompy; iy++){
		    for(int ix=0; ix<ncompx; ix++){
			X[iy][ix]*=ct*Ux[ix]+st*Uy[iy];
			Y[iy][ix]*=-st*Ux[ix]+ct*Uy[iy];
		    }
		}
		cfft2(seotfk,1);
		spmulcreal(gx[ii0][isa]->p,si,seotfk->p, wvlsig);
		cfft2(seotfj,1);
		spmulcreal(gy[ii0][isa]->p,si,seotfj->p, wvlsig);
		if(i0scale){
		    double scale=area[isa]/dsum(i0[ii0][isa]);
		    dscale(i0[ii0][isa],scale);
		    dscale(gx[ii0][isa],scale);
		    dscale(gy[ii0][isa],scale);
		}
	    }//for illt
	}//isa
	cfree(sepsf);
	cfree(seotfj);
	cfree(seotfk);
    }//iwvl
#if SAVE_OTF==1
    ccellwrite(saveotf,"powfs%d_saveotf",ipowfs);
    ccellwrite(savepsf,"powfs%d_savepsf",ipowfs);
    ccellwrite(saveotfetf,"powfs%d_saveotfetf",ipowfs);
#endif
}
