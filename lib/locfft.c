/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "../math/mathdef.h"
#include "locfft.h"

/**
   \file locfft.c

   For routines that embed OPDs defined on loc to square(rectangular) array and do fft on it.
*/


locfft_t *locfft_init(const loc_t *loc,       /**<[in] The loc*/
		      const dmat *amp,        /**<[in] The amplitude*/
		      const lmat *fftsize,    /**<[in] The suggested size for FFT*/
		      const dmat *wvl,        /**<[in] The wavelength*/
		      double fieldstop        /**<[in] Size of field stop (radian) if used*/
    ){
    const int nwvl=wvl->nx*wvl->ny;
    locfft_t *locfft=calloc(sizeof(locfft_t), 1);
    locfft->embed=cellnew(nwvl, 1);
    locfft->nembed=lnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	if(iwvl==0 || (fftsize && fftsize->p[iwvl]>0 && fftsize->p[iwvl]!=locfft->nembed->p[0])){
	    locfft->nembed->p[iwvl]=fftsize?fftsize->p[iwvl]:0;
	    locfft->embed->p[iwvl]=loc_create_embed(&locfft->nembed->p[iwvl], loc, 2, 1);
	}else{
	    locfft->embed->p[iwvl]=locfft->embed->p[0];
	    locfft->nembed->p[iwvl]=locfft->nembed->p[0];
	}
    }
    locfft->wvl=dref_reshape(wvl, nwvl, 1);
    locfft->amp=amp;
    locfft->loc=loc;
    locfft->ampsum=dsum(amp);
    locfft->ampnorm=dsumsq(amp);
    if(fieldstop){
	locfft->fieldstop=fieldstop;
	locfft->fieldmask=cellnew(nwvl, 1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    int nembed=locfft->nembed->p[iwvl];
	    locfft->fieldmask->p[iwvl]=dnew(nembed, nembed);
	    double dtheta=wvl->p[iwvl]/(loc->dx*nembed);//sampling of psf
	    double radius=fieldstop/(dtheta*2);
	    dcircle(locfft->fieldmask->p[iwvl], nembed/2+1, nembed/2+1, 1, 1, radius, 1);
	    dfftshift(locfft->fieldmask->p[iwvl]);
	}
    }
    return locfft;
}

void locfft_free(locfft_t *locfft){
    if(!locfft || !locfft->embed) return;
    for(int iwvl=locfft->embed->nx-1; iwvl>=0; iwvl--){
	if(locfft->embed->p[iwvl]==locfft->embed->p[0]){
	    locfft->embed->p[iwvl]=0;
	}
    }
    lcellfree(locfft->embed);
    lfree(locfft->nembed);
    dcellfree(locfft->fieldmask);
    dfree(locfft->wvl);
}
/**
   Computes strehl from OPD without doing FFT. The strehl is simply 
   
   \f$s=\sum(A*exp[\frac{2\pi i}{\lambda}*\textrm{OPD}]) \f$
   
   where A is the amplitude map.
*/
static dcomplex strehlcomp(const dmat *iopdevl, const dmat *amp, const double wvl){
    dcomplex i2pi=I*2*M_PI/wvl;
    dcomplex strehl=0;
    for(int iloc=0; iloc<iopdevl->nx; iloc++){
	strehl+=amp->p[iloc]*cexp(i2pi*iopdevl->p[iloc]);
    }
    return strehl;
}
/**
   Computes PSF from OPD by FFT. The PSF is computed as

   \f$\textrm{PSF}=\mathcal{F}[A\times exp[\frac{2\pi i}{\lambda}*\textrm{OPD}]]\f$

   The peak value (center) in the computed PSF is normalized by the peak value
   in the differaction limited PSF. In other words, the peak value in the
   computed PSF is the Strehl. Keep this in mind when you compute enclosed
   energy.
   
   Extract center part of psfsize.
*/
void locfft_psf(ccell **psf2s, locfft_t *locfft, dmat *opd, lmat *psfsize, int sum2one){
    long nwvl=locfft->wvl->nx;
    if(!*psf2s){
	*psf2s=cellnew(nwvl, 1);
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++)
#if _OPENMP>=200805
#pragma omp task
#endif
    {
	if(psfsize && psfsize->p[iwvl]==1){
	    if(!(*psf2s)->p[iwvl]){
		(*psf2s)->p[iwvl]=cnew(1,1);
	    }
	    (*psf2s)->p[iwvl]->p[0]=strehlcomp(opd, locfft->amp, locfft->wvl->p[iwvl]);
	}else{
	    long nembed=locfft->nembed->p[iwvl];
	    long *embed=locfft->embed->p[iwvl]->p;
	    const double *amp=locfft->amp->p;
	    cmat *psf2=cnew(nembed,nembed);

	    int use1d;
	    int use1d_enable=0;
	    if(psfsize && psfsize->p[iwvl]<nembed && use1d_enable){/*Want smaller PSF. */
		use1d=1;
		//cfft2partialplan(psf2, psfsize->p[iwvl], -1);
	    }else{
		use1d=0;
		//cfft2plan(psf2, -1);
	    }

	    dcomplex i2pi=I*2*M_PI/locfft->wvl->p[iwvl];
	    for(int iloc=0; iloc<opd->nx; iloc++){
		psf2->p[embed[iloc]]=amp[iloc]*cexp(i2pi*opd->p[iloc]);
	    }
	    if(use1d==1){
		cfft2partial(psf2, psfsize->p[iwvl], -1);
	    }else{
		cfft2(psf2,-1);
	    }
	    if(!(*psf2s)->p[iwvl] && (!psfsize || psfsize->p[iwvl]==nembed)){/*just reference */
		cfftshift(psf2);
		(*psf2s)->p[iwvl]=psf2;
	    }else{/*create a new array, smaller. */
		if(!(*psf2s)->p[iwvl]){
		    (*psf2s)->p[iwvl]=cnew(psfsize->p[iwvl], psfsize->p[iwvl]);
		}
		ccpcorner2center((*psf2s)->p[iwvl], psf2);
		cfree(psf2); 
	    }
	}
	double psfnorm;
	if(sum2one){
	    psfnorm=1./(sqrt(locfft->ampnorm)*locfft->nembed->p[iwvl]);
	}else{/**PSF max is strehl*/
	    psfnorm=1./locfft->ampsum;
	}
	if(fabs(psfnorm-1)>1.e-15) {
	    cscale((*psf2s)->p[iwvl], psfnorm);
	}
    }
#if _OPENMP>=200805
#pragma omp taskwait
#endif
}
/**
   Apply a field stop
*/
void locfft_fieldstop(locfft_t *locfft, dmat *opd, dmat *wvlwts){
    int nwvl=locfft->wvl->nx;
    ccell *wvfs=cellnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	int nembed=locfft->nembed->p[iwvl];
	lmat *embed=locfft->embed->p[iwvl];
	cmat *wvf=cnew(nembed, nembed);
	wvfs->p[iwvl]=wvf;
	//cfft2plan(wvf, -1); //cfft2plan(wvf, 1);
	double wvl=locfft->wvl->p[iwvl];
	dcomplex i2pi=2*M_PI/wvl*I;
	const double *amp=locfft->amp->p;
	for(int iloc=0; iloc<opd->nx; iloc++){
	    wvf->p[embed->p[iloc]]=amp[iloc]*cexp(i2pi*opd->p[iloc]);
	}
	cfft2(wvf, -1);
	ccwmd(wvf, locfft->fieldmask->p[iwvl], 1);
	cfft2(wvf, 1);
    }
    if(nwvl>1){
	/*Please fix this case. Need to do phase unwrapping first and average OPD
	 * for different wavelength result*/
	error("Not implemented yet\n");
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	double wvl=locfft->wvl->p[iwvl];
	double wvlh=wvl*0.5;
	double kki=wvl/(2*M_PI);
	cmat *wvf=wvfs->p[iwvl];
	lmat *embed=locfft->embed->p[iwvl];
	for(int iloc=0; iloc<opd->nx; iloc++){
	    double val=carg(wvf->p[embed->p[iloc]])*kki;
	    if(fabs(val-opd->p[iloc])>wvlh){//need phase unwrapping
		warning_once("phase unwrapping is needed\n");
		double diff=fmod(val-opd->p[iloc]+wvlh, wvl);
		if(diff<0) diff+=wvl;
		opd->p[iloc]+=diff-wvlh;
	    }else{
		opd->p[iloc]=val;
	    }
	}
    }
    ccellfree(wvfs);
}
