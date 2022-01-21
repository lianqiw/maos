/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
	const int nwvl=NX(wvl)*NY(wvl);
	locfft_t* locfft=mycalloc(1, locfft_t);
	locfft->embed=lcellnew(nwvl, 1);
	locfft->nembed=lnew(nwvl, 1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(iwvl==0||(fftsize&&P(fftsize, iwvl)>0&&P(fftsize, iwvl)!=P(locfft->nembed, 0))){
			P(locfft->nembed, iwvl)=fftsize?P(fftsize, iwvl):0;
			P(locfft->embed, iwvl)=loc_create_embed(&P(locfft->nembed, iwvl), loc, oversize, 1);
		} else{
			P(locfft->embed, iwvl)=lref(P(locfft->embed, 0));
			P(locfft->nembed, iwvl)=P(locfft->nembed, 0);
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
			int nembed=P(locfft->nembed, iwvl);
			P(locfft->fieldmask, iwvl)=dnew(nembed, nembed);
			real dtheta=P(wvl, iwvl)/(loc->dx*nembed);//sampling of psf
			real radius=fieldstop/(dtheta*2);
			dcircle(P(locfft->fieldmask, iwvl), nembed/2+1, nembed/2+1, 1, 1, radius, 1);
			dfftshift(P(locfft->fieldmask, iwvl));
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
	for(int iloc=0; iloc<NX(iopdevl); iloc++){
		strehl+=P(amp, iloc)*cexp(i2pi*P(iopdevl, iloc));
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
	long nwvl=NX(locfft->wvl);
	if(!*psf2sp){
		*psf2sp=ccellnew(nwvl, 1);
	}
	ccell* psf2s=*psf2sp;
	if(NX(opd)!=NX(locfft->amp)){
		error("The length of opd should be %ld, but is %ld\n", NX(locfft->amp), NX(opd));
	}
OMP_TASK_FOR(4)	
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		if(psfsize&&P(psfsize, iwvl)==1){
			if(!P(psf2s, iwvl)){
				P(psf2s, iwvl)=cnew(1, 1);
			}
			P(P(psf2s, iwvl), 0)=strehlcomp(opd, locfft->amp, P(locfft->wvl, iwvl));
		} else{
			TIM(0);
			long nembed=P(locfft->nembed, iwvl);
			long* embed=P(P(locfft->embed, iwvl));
			const real* amp=P(locfft->amp);
			const int ref=!psfsize||P(psfsize, iwvl)==nembed;
			cmat* psf2=0;
			if(ref){//Full PSF is returned
				if(!P(psf2s, iwvl)){
					P(psf2s, iwvl)=cnew(nembed, nembed);
				}
				psf2=P(psf2s, iwvl);
			} else{//Crop of PSF is returned.
				psf2=cnew(nembed, nembed);
			}

			int use1d=0;
#define USE1D_ENABLED 1
#if     USE1D_ENABLED
			if(psfsize&&P(psfsize, iwvl)+200<nembed){/*Want smaller PSF. */
				use1d=1;
			}
#endif

			comp i2pi=COMPLEX(0, 2*M_PI/P(locfft->wvl, iwvl));
			for(int iloc=0; iloc<NX(opd); iloc++){
				P(psf2, embed[iloc])=amp[iloc]*cexp(i2pi*P(opd, iloc));
			}
			TIM(1);
			if(use1d==1){
				cfft2partial(psf2, P(psfsize, iwvl), -1);
			} else{
				cfft2(psf2, -1);
			}
			TIM(2);
			if(ref){/*just reference */
				cfftshift(psf2);
			} else{/*create a new array, smaller. */
				if(!P(psf2s, iwvl)){
					P(psf2s, iwvl)=cnew(P(psfsize, iwvl), P(psfsize, iwvl));
				}
				ccpcorner2center(P(psf2s, iwvl), psf2);
				cfree(psf2);
			}
			TIM(3);
#if TIMING
			info2("locfft_psf(%d:%ldx%ld): exp %.4f, fft %.4f (%.2f GFLOPS), abs2 %.2f.\n", iwvl, nembed, nembed,
				tk1-tk0, tk2-tk1, 8L*(use1d?P(psfsize, iwvl):nembed)*nembed*log2(nembed)/(tk2-tk1)*1e-9, tk3-tk2);
#endif		
		}
		real psfnorm;
		if(sum2one){/**PSF sum to 1*/
			psfnorm=1./(sqrt(locfft->ampnorm)*P(locfft->nembed, iwvl));
		} else{/**PSF max is strehl*/
			psfnorm=1./locfft->ampsum;
		}
		if(fabs(psfnorm-1)>1.e-15){
			cscale(P(psf2s, iwvl), psfnorm);
		}
	}
}
/**
   Apply a field stop to the OPD.
*/
void locfft_fieldstop(const locfft_t* locfft, dmat* opd, const dmat* wvlwts){
	int nwvl=NX(locfft->wvl);
	if(nwvl>1){
		warning("Not tested for multi-wavelength case yet.\n");
	}
	ccell* wvfs=ccellnew(nwvl, 1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		int nembed=P(locfft->nembed, iwvl);
		lmat* embed=P(locfft->embed, iwvl);
		cmat* wvf=cnew(nembed, nembed);
		P(wvfs, iwvl)=wvf;
		//cfft2plan(wvf, -1); //cfft2plan(wvf, 1);
		real wvl=P(locfft->wvl, iwvl);
		comp i2pi=COMPLEX(0, 2*M_PI/wvl);
		const real* amp=P(locfft->amp);
		for(int iloc=0; iloc<NX(opd); iloc++){
			P(wvf, P(embed, iloc))=amp[iloc]*cexp(i2pi*P(opd, iloc));
		}
		cfft2(wvf, -1);
		ccwmd(wvf, P(locfft->fieldmask, iwvl), 1);
		cfft2(wvf, 1);
	}
	if(nwvl>1){
	/*Please fix this case. Need to do phase unwrapping first and average OPD
	 * for different wavelength result*/
		error("Not implemented yet\n");
	}
	dmat* opdold=ddup(opd); dzero(opd);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		real wvl=P(locfft->wvl, iwvl);
		real wvlh=wvl*0.5;
		real kki=wvl/(2*M_PI);
		cmat* wvf=P(wvfs, iwvl);
		lmat* embed=P(locfft->embed, iwvl);
		for(int iloc=0; iloc<NX(opd); iloc++){
			real val=carg(P(wvf, P(embed, iloc)))*kki;
			if(fabs(val-P(opdold, iloc))>wvlh){//need phase unwrapping
				warning_once("phase unwrapping is needed\n");
				real diff=fmod(val-P(opdold, iloc)+wvlh, wvl);
				if(diff<0) diff+=wvl;
				val=(diff-wvlh)+P(opdold, iloc);
			}
			P(opd, iloc)+=P(wvlwts, iwvl)*val;
		}
	}
	ccellfree(wvfs);
}
/**
 * Apply the phase function out=in*exp(I*pi*(x^2+y^2)/(z*wvl)
 * */
static void apply_h(cmat* out, const cmat* in, real dx, real wvl, real z){
	real coeff=(M_PI/(z*wvl))*dx*dx;
	long nx2=NX(out)/2;
	long ny2=NY(out)/2;
	if(in){
		for(long iy=0; iy<NY(out); iy++){
			real ysq=(iy-ny2)*(iy-ny2);
			for(long ix=0; ix<NX(out); ix++){
				real ph=coeff*((ix-nx2)*(ix-nx2)+ysq);
				P(out, ix, iy)=P(in, ix, iy)*EXPI(ph);
			}
		}
	} else{
		for(long iy=0; iy<NY(out); iy++){
			real ysq=(iy-ny2)*(iy-ny2);
			for(long ix=0; ix<NX(out); ix++){
				real ph=coeff*((ix-nx2)*(ix-nx2)+ysq);
				P(out, ix, iy)=EXPI(ph);
			}
		}
	}
}
/**
 * Fresnel propagation using either FFT or angular spectrum method.
 *
 * For FFT method, the output sampling is z*wvl/(dxin*nxin). For angular
 * spectrum method, the output sampling equals to the input sampling. When the
 * output sampling equals, the two methods agree. Notice that at near field when
 * dp is smaller than dx, the output of angular spectrum method is periodic.
 *
 * The output is normalized so that its squared sum equals to
 * nx*ny*sum(abs(in)^2). Input argument scale should be supplied to cancel this
 * scaling. 
 *
 * The brutal force method (method=0) has the exact solution without using
 * fresnel approximation. It is however very slow. The results agrees with
 * method=1 or 2 when the sampling agrees and the approximation holds.
 */
void fresnel_prop(cmat** pout, /**<Output complex field. Sampling depends on method. Should not be empty if method is 0.*/
	real* pdxout,  /**<Sampling of output, should be supplied if method=0*/
	const cmat* in, /**<Input complex field, , Amp*exp(-i*k*opd)*/
	real dxin,   /**<Spatial sampling of in*/
	real wvl,    /**<Wavelength*/
	real z,      /**<Propagation distance*/
	real scale,  /**<Scaling to be applied to out to preserve energy. 0: determine from input*/
	int method   /**<Propagation method, 0: brutal force, 1: FFT, 2: angular spectrum, -1: auto between 1 and 2.*/
){
	if(method){
		if(!*pout){
			*pout=cnew(NX(in), NY(in));
		} else{
			if(NX(*pout)!=NX(in)){
				error("Output array should have the same dimension as input array.\n");
				return;
			}
		}

		real ratio=pow(fabs(z), 3)*8*wvl/pow(dxin*NX(in), 4);
		if(ratio<1e6){
			real zmin=pow(1e6/ratio,1./3)*z*1.1;
			warning_once("ratio is only %g, will use a two step propagation for %g and %g.\n", ratio, z+zmin, -zmin);
			real dxout1=0;
			fresnel_prop(pout, &dxout1, in, dxin, wvl, z+zmin, scale, method);
			fresnel_prop(pout, pdxout, *pout, dxout1, wvl, -zmin, scale, method);
			return;
		}
		
		if(method==-1){
			//when dp is smaller than dx, use angular spectrum method.
			real dp=z* wvl/(dxin*NX(in));
			if(dp<dxin){
				method=2;
			}else{
				method=1;
			}
		}
	}
	if(scale<=0){
		scale=1./sqrt(csumsq(in)*PN(in));
	}
	cmat* out=*pout;
	if(method==0){//brutal force method
		//Do not multiply dx*dy to the integrand or apply 1/(wvl*z) to the output 
		//Due to normalization in discrete domain.
		//Confirmed to agree with method=1 or 2 when sampling agrees.
		real k=2.*M_PI/wvl;
		real zz=z*z;
		real dxout=*pdxout;
		if(!dxout){
			*pdxout=dxout=z*wvl/(dxin*NX(in));
		}
		if(!out){
			out=*pout=cnew(32,32);
		}
		OMP_TASK_FOR(4)
		for(long iyo=0; iyo<NY(out); iyo++){
			const real y=(iyo-NY(out)/2)*dxout;
			for(long ixo=0; ixo<NX(out); ixo++){
				const real x=(ixo-NX(out)/2)*dxout;
				comp res=0;
				for(long iyi=0; iyi<NY(in); iyi++){
					real yp=(iyi-NY(in)/2)*dxin;
					for(long ixi=0; ixi<NX(in); ixi++){
						real xp=(ixi-NX(in)/2)*dxin;
						real r2=pow(x-xp, 2)+pow(y-yp, 2)+zz;
						real ph=k*sqrt(r2);
						res+=P(in, ixi, iyi)*EXPI(ph)*(zz/r2);
					}
				}
				P(out, ixo, iyo)=res*scale;
			}
		}
	} else if(method==1){//FFT method.
		apply_h(out, in, dxin, wvl, z);
		cfftshift(out);//FFT zero frequency is at corner.
		cfft2(out, -1);
		cfftshift(out);
		*pdxout=fabs(z)*wvl/(dxin*NX(in));
		apply_h(out, out, *pdxout, wvl, z);
		cscale(out, scale);
	} else if(method==2){//angular spectrum
		cmat* in2=cdup(in);	//cfftshift(in2);//no need
		cfft2(in2, -1);
		apply_h(out, NULL, dxin, wvl, z);//cfftshift(out);//no need
		cfft2(out, -1);
		ccwm(out, in2);
		cfree(in2);
		cfft2(out, 1);
		cfftshift(out);
		cscale(out, scale/(PN(out)));
		*pdxout=dxin;
	} else{
		error("Invalid method\n");
	}
}