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
#include "../math/mathdef.h"
#include "locfft.h"
#define TIMING 0
#if TIMING == 1
#define TIM(A) real tk##A=myclockd()
#else
#define TIM(A)
#endif

/**
   For routines that embed OPDs defined on loc to square(rectangular) array and do fft on it.
*/


locfft_t* locfft_init(loc_t* loc,       /**<[in] The loc*/
	const dmat* amp,        /**<[in] The amplitude*/
	const dmat* wvl,        /**<[in] The wavelength*/
	const lmat* fftsize,    /**<[in] The suggested size for FFT*/
	const real oversize,  /**<[in] Factor of oversize. 2 fot FFT*/
	real fieldstop        /**<[in] Size of field stop (radian) if used*/
){
	const int nwvl=wvl->nx*wvl->ny;
	locfft_t* locfft=mycalloc(1, locfft_t);
	locfft->embed=lcellnew(nwvl, 1);
	locfft->nembed=lnew(nwvl, 1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(iwvl==0||(fftsize&&P(fftsize,iwvl)>0&&P(fftsize,iwvl)!=P(locfft->nembed,0))){
			P(locfft->nembed,iwvl)=fftsize?P(fftsize,iwvl):0;
			P(locfft->embed,iwvl)=loc_create_embed(PP(locfft->nembed,iwvl), loc, oversize, 1);
		} else{
			P(locfft->embed,iwvl)=lref(P(locfft->embed,0));
			P(locfft->nembed,iwvl)=P(locfft->nembed,0);
		}
	}
	locfft->wvl=dref_reshape(wvl, nwvl, 1);
	locfft->amp=amp;
	locfft->loc=loc;
	locfft->ampsum=dsum(amp);
	locfft->ampnorm=dsumsq(amp);
	if(fieldstop){
		locfft->fieldstop=fieldstop;
		locfft->fieldmask=dcellnew(nwvl, 1);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
			int nembed=P(locfft->nembed,iwvl);
			P(locfft->fieldmask,iwvl)=dnew(nembed, nembed);
			real dtheta=P(wvl,iwvl)/(loc->dx*nembed);//sampling of psf
			real radius=fieldstop/(dtheta*2);
			dcircle(P(locfft->fieldmask,iwvl), nembed/2+1, nembed/2+1, 1, 1, radius, 1);
			dfftshift(P(locfft->fieldmask,iwvl));
		}
	}
	return locfft;
}
/**
   Frees the struct
*/
void locfft_free(locfft_t* locfft){
	if(!locfft) return;
	cellfree(locfft->embed);
	cellfree(locfft->nembed);
	cellfree(locfft->fieldmask);
	cellfree(locfft->wvl);
	free(locfft);
}
/**
   Computes strehl from OPD without doing FFT. The strehl is simply

   \f$s=\sum(A*exp[\frac{2\pi i}{\lambda}*\textrm{OPD}]) \f$

   where A is the amplitude map.
*/
static comp strehlcomp(const dmat* iopdevl, const dmat* amp, const real wvl){
	comp i2pi=COMPLEX(0, 2*M_PI/wvl);
	comp strehl=0;
	for(int iloc=0; iloc<iopdevl->nx; iloc++){
		strehl+=P(amp,iloc)*cexp(i2pi*P(iopdevl,iloc));
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
void locfft_psf(ccell** psf2sp, const locfft_t* locfft, const dmat* opd, const lmat* psfsize, int sum2one){
	long nwvl=locfft->wvl->nx;
	if(!*psf2sp){
		*psf2sp=ccellnew(nwvl, 1);
	}
	ccell* psf2s=*psf2sp;
	if(opd->nx!=locfft->amp->nx){
		error("The length of opd should be %ld, but is %ld\n", locfft->amp->nx, opd->nx);
	}
	for(int iwvl=0; iwvl<nwvl; iwvl++)
#if _OPENMP>=200805
#pragma omp task
#endif
	{
		if(psfsize&&P(psfsize,iwvl)==1){
			if(!P(psf2s,iwvl)){
				P(psf2s,iwvl)=cnew(1, 1);
			}
			P(P(psf2s,iwvl),0)=strehlcomp(opd, locfft->amp, P(locfft->wvl,iwvl));
		} else{
			TIM(0);
			long nembed=P(locfft->nembed,iwvl);
			long* embed=P(locfft->embed,iwvl)->p;
			const real* amp=locfft->amp->p;
			const int ref=!psfsize||P(psfsize,iwvl)==nembed;
			cmat* psf2=0;
			if(ref){//Full PSF is returned
				if(!P(psf2s,iwvl)){
					P(psf2s,iwvl)=cnew(nembed, nembed);
				}
				psf2=P(psf2s,iwvl);
			} else{//Crop of PSF is returned.
				psf2=cnew(nembed, nembed);
			}

			int use1d=0;
#define USE1D_ENABLED 1
#if     USE1D_ENABLED
			if(psfsize&&P(psfsize,iwvl)+200<nembed){/*Want smaller PSF. */
				use1d=1;
			}
#endif

			comp i2pi=COMPLEX(0, 2*M_PI/P(locfft->wvl,iwvl));
			for(int iloc=0; iloc<opd->nx; iloc++){
				P(psf2,embed[iloc])=amp[iloc]*cexp(i2pi*P(opd,iloc));
			}
			TIM(1);
			if(use1d==1){
				cfft2partial(psf2, P(psfsize,iwvl), -1);
			} else{
				cfft2(psf2, -1);
			}
			TIM(2);
			if(ref){/*just reference */
				cfftshift(psf2);
			} else{/*create a new array, smaller. */
				if(!P(psf2s,iwvl)){
					P(psf2s,iwvl)=cnew(P(psfsize,iwvl), P(psfsize,iwvl));
				}
				ccpcorner2center(P(psf2s,iwvl), psf2);
				cfree(psf2);
			}
			TIM(3);
#if TIMING
			info2("locfft_psf(%d:%ldx%ld): exp %.4f, fft %.4f (%.2f GFLOPS), abs2 %.2f.\n", iwvl, nembed, nembed,
				tk1-tk0, tk2-tk1, 8L*(use1d?P(psfsize,iwvl):nembed)*nembed*log2(nembed)/(tk2-tk1)*1e-9, tk3-tk2);
#endif		
		}
		real psfnorm;
		if(sum2one){/**PSF sum to 1*/
			psfnorm=1./(sqrt(locfft->ampnorm)*P(locfft->nembed,iwvl));
		} else{/**PSF max is strehl*/
			psfnorm=1./locfft->ampsum;
		}
		if(fabs(psfnorm-1)>1.e-15){
			cscale(P(psf2s,iwvl), psfnorm);
		}
	}
#if _OPENMP>=200805
#pragma omp taskwait
#endif
}
/**
   Apply a field stop to the OPD.
*/
void locfft_fieldstop(const locfft_t* locfft, dmat* opd, const dmat* wvlwts){
	int nwvl=locfft->wvl->nx;
	if(nwvl>1){
		warning("Not tested for multi-wavelength case yet.\n");
	}
	ccell* wvfs=ccellnew(nwvl, 1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		int nembed=P(locfft->nembed,iwvl);
		lmat* embed=P(locfft->embed,iwvl);
		cmat* wvf=cnew(nembed, nembed);
		P(wvfs,iwvl)=wvf;
		//cfft2plan(wvf, -1); //cfft2plan(wvf, 1);
		real wvl=P(locfft->wvl,iwvl);
		comp i2pi=COMPLEX(0, 2*M_PI/wvl);
		const real* amp=locfft->amp->p;
		for(int iloc=0; iloc<opd->nx; iloc++){
			P(wvf, P(embed,iloc))=amp[iloc]*cexp(i2pi*P(opd,iloc));
		}
		cfft2(wvf, -1);
		ccwmd(wvf, P(locfft->fieldmask,iwvl), 1);
		cfft2(wvf, 1);
	}
	if(nwvl>1){
	/*Please fix this case. Need to do phase unwrapping first and average OPD
	 * for different wavelength result*/
		error("Not implemented yet\n");
	}
	dmat* opdold=ddup(opd); dzero(opd);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		real wvl=P(locfft->wvl,iwvl);
		real wvlh=wvl*0.5;
		real kki=wvl/(2*M_PI);
		cmat* wvf=P(wvfs,iwvl);
		lmat* embed=P(locfft->embed,iwvl);
		for(int iloc=0; iloc<opd->nx; iloc++){
			real val=carg(P(wvf,P(embed,iloc)))*kki;
			if(fabs(val-P(opdold,iloc))>wvlh){//need phase unwrapping
				warning_once("phase unwrapping is needed\n");
				real diff=fmod(val-P(opdold,iloc)+wvlh, wvl);
				if(diff<0) diff+=wvl;
				val=(diff-wvlh)+P(opdold,iloc);
			}
			P(opd,iloc)+=P(wvlwts,iwvl)*val;
		}
	}
	ccellfree(wvfs);
}
