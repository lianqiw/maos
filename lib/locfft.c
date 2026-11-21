/*
  Copyright 2009-2026 Lianqi Wang
  
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


/**
   Initialize the locfft_t struct that embeds OPD on loc into square arrays for FFT.

   @param[in] loc The loc (referenced, not copied).
   @param[in] amp The amplitude (referenced, not copied).
   @param[in] wvl The wavelength (m).
   @param[in] fftsize The suggested size for FFT, nwvl values. Can be NULL.
   @param[in] oversize Factor of oversize. 2 for Nyquist sampling.
   @param[in] fieldstop Size of field stop (radian) if used. 0 to disable.
   @return Newly allocated struct; free with locfft_free().
*/
locfft_t* locfft_init(loc_t* loc, const dmat* amp, const dmat* wvl, const lmat* fftsize, const real oversize, real fieldstop){
	const int nwvl=(int)(NX(wvl)*NY(wvl));
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
			long nembed=P(locfft->nembed, iwvl);
			P(locfft->fieldmask, iwvl)=dnew(nembed, nembed);
			real dtheta=P(wvl, iwvl)/(loc->dx*nembed);//sampling of psf
			real radius=fieldstop/(dtheta*2);
			dbg("field stop is %g, dtheta is %g, radius is %g, nembed=%ld.\n", fieldstop, dtheta, radius, nembed);
			dcircle(P(locfft->fieldmask, iwvl), nembed/2+1, nembed/2+1, 1, 1, radius, 1);
			dfftshift(P(locfft->fieldmask, iwvl));
		}
	}
	return locfft;
}
/**
   Frees the struct

   @param[in] locfft The struct to free. Can be NULL.
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

   @param[in] iopdevl The OPD (m).
   @param[in] amp The amplitude map.
   @param[in] wvl The wavelength (m).
   @return Complex sum (not normalized).
*/
static comp strehlcomp(const dmat* iopdevl, const dmat* amp, const real wvl){
	real wvk=TWOPI/wvl;
	comp strehl=0;
	for(int iloc=0; iloc<NX(iopdevl); iloc++){
		strehl+=P(amp, iloc)*EXPI(wvk*P(iopdevl, iloc));
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

   @param[in,out] psf2sp The complex PSF for each wavelength. Allocated if NULL, otherwise overwritten.
   @param[in] locfft The locfft struct.
   @param[in] opd The OPD defined on locfft->loc (m).
   @param[in] psfsize Size of PSF to extract for each wavelength. NULL for full size. 1 for Strehl only.
   @param[in] sum2one If set, normalize PSF to sum to 1; otherwise peak is Strehl.
*/
void locfft_psf(ccell** psf2sp, const locfft_t* locfft, const dmat* opd, const lmat* psfsize, int sum2one){
	int nwvl=(int)NX(locfft->wvl);
	if(!*psf2sp){
		*psf2sp=ccellnew(nwvl, 1);
	}
	ccell* psf2s=*psf2sp;
	if(NX(opd)!=NX(locfft->amp)){
		error("The length of opd should be %ld, but is %ld\n", NX(locfft->amp), NX(opd));
		return;
	}
OMP_TASK_FOR(4)	
	for(long iwvl=0; iwvl<nwvl; iwvl++){
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
				}else{
					czero(P(psf2s, iwvl));//only the embed points are written below
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

			real wvk=TWOPI/P(locfft->wvl, iwvl);
			for(int iloc=0; iloc<NX(opd); iloc++){
				P(psf2, embed[iloc])=amp[iloc]*EXPI(wvk*P(opd, iloc));
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
			info2("locfft_psf(%ld:%ldx%ld): exp %.4f, fft %.4f (%.2f GFLOPS), abs2 %.2f.\n", iwvl, nembed, nembed,
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

   @param[in] locfft The locfft struct with fieldmask set.
   @param[in,out] opd The OPD (m), replaced by the field stop filtered OPD.
   @param[in] wvlwts Weights of each wavelength.
*/
void locfft_fieldstop(const locfft_t* locfft, dmat* opd, const dmat* wvlwts){
	int nwvl=(int)NX(locfft->wvl);
	if(nwvl>1){
		warning("Not tested for multi-wavelength case yet.\n");
	}
	ccell* wvfs=ccellnew(nwvl, 1);
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		long nembed=P(locfft->nembed, iwvl);
		lmat* embed=P(locfft->embed, iwvl);
		cmat* wvf=cnew(nembed, nembed);
		P(wvfs, iwvl)=wvf;
		//cfft2plan(wvf, -1); //cfft2plan(wvf, 1);
		real wvl=P(locfft->wvl, iwvl);
		real wvk=TWOPI/wvl;
		const real* amp=P(locfft->amp);
		for(int iloc=0; iloc<NX(opd); iloc++){
			P(wvf, P(embed, iloc))=amp[iloc]*EXPI(wvk*P(opd, iloc));
		}
		cfft2(wvf, -1);
		ccwmd(wvf, P(locfft->fieldmask, iwvl));
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
	dfree(opdold);
	ccellfree(wvfs);
}
/**
   Apply the phase function out=in*exp(I*pi*(x^2+y^2)/(z*wvl)).

   @param[out] out The output complex field. Already allocated.
   @param[in] in The input complex field. If NULL, out is set to the phase function only. Can be the same as out.
   @param[in] dx Spatial sampling (m).
   @param[in] wvl Wavelength (m).
   @param[in] z Propagation distance (m).
*/
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
   Fresnel propagation using either FFT or angular spectrum method.

   For FFT method, the output sampling is z*wvl/(dxin*nxin). For angular
   spectrum method, the output sampling equals to the input sampling. When the
   output sampling equals, the two methods agree. The angular spectrum method is
   unitary (propagating by -z reverses it) and its output is periodic.

   The output is normalized so that its squared sum equals to
   nx*ny*sum(abs(in)^2). Input argument scale should be supplied to cancel this
   scaling. 

   The brutal force method (method=0) has the exact solution without using
   fresnel approximation. It is however very slow. The results agrees with
   method=1 or 2 when the sampling agrees and the approximation holds.

   @param[in,out] pout Output complex field. Sampling depends on method. Allocated if NULL (32x32 for method 0, same size as in otherwise).
   @param[in,out] pdxout Sampling of output. For method=0 it can be supplied, otherwise z*wvl/(dxin*nx) is used. Set on output.
   @param[in] in Input complex field, Amp*exp(-i*k*opd).
   @param[in] dxin Spatial sampling of in (m).
   @param[in] wvl Wavelength (m).
   @param[in] z Propagation distance (m).
   @param[in] scale Scaling to be applied to out to preserve energy. 0: determine from input.
   @param[in] method Propagation method, 0: brutal force, 1: FFT, 2: angular spectrum, -1: auto between 1 and 2.
*/
void fresnel_prop(cmat** pout, real* pdxout, const cmat* in, real dxin, real wvl, real z, real scale, int method){
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
		real wvk=TWOPI/wvl;
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
						real ph=wvk*sqrt(r2);
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
	} else if(method==2){//angular spectrum: apply the transfer function in the frequency domain.
		//The impulse response sampled in the spatial domain is aliased when dp<dxin, where this method is used.
		if(out!=in) ccp(&out, in);
		cfft2(out, -1);
		const long nx=NX(out), ny=NY(out);
		const real coeff=-M_PI*wvl*z/(dxin*dxin);//H=exp(-i*pi*wvl*z*(fx^2+fy^2))
		for(long iy=0; iy<ny; iy++){
			const real fy=(real)(iy<(ny+1)/2?iy:iy-ny)/ny;
			for(long ix=0; ix<nx; ix++){
				const real fx=(real)(ix<(nx+1)/2?ix:ix-nx)/nx;
				P(out, ix, iy)*=EXPI(coeff*(fx*fx+fy*fy));
			}
		}
		cfft2(out, 1);
		cscale(out, scale/sqrt((real)PN(out)));//same normalization as the FFT method
		*pdxout=dxin;
	} else{
		error("Invalid method\n");
	}
}
