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
/**
   \file nafocus.c
   Sodium layer focus tracking analytic results
*/
#include "skyc.h"
#include "nafocus.h"
/**
   Sodium layer height PSD beta2*nu^alpha (not divided by 2pi).

   @param[in] nu Temporal frequency (Hz).
   @param[in] alpha Power law exponent of the PSD.
   @param[in] beta2 Scaling of the PSD.
   @return PSD value at nu.
*/
static inline real nafocus_NaPSD(real nu, real alpha, real beta2){
	return beta2*pow(nu, alpha);/*we don't divide 2pi */
}

/**
   Compute Zoom optics focus corrector open loop transfer function

   @param[in] nu Frequency (Hz).
   @param[in] fs Sampling frequency (Hz).
   @param[in] tau Time delay (s).
   @param[in] zeta Damping of zoom corrector.
   @param[in] fzc Frequency of zoom corrector (Hz).
   @return Complex open loop transfer function at nu.
*/
static inline comp nafocus_Hol(real nu, real fs, real tau, real zeta, real fzc){
	const comp s=COMPLEX(0, 2*M_PI*nu);
	const comp fsovers=fs/s;
	const comp Hint=fsovers;
	const comp Hwfs=(1.-cexp(-1./fsovers))*fsovers;
	const comp Hlag=cexp(-s*tau);
	const comp Hdac=1;
	const real omega=fzc*2*M_PI;
	const comp Hzoom=omega*omega/(s*s+2*zeta*omega*s+omega*omega);
	const comp Hol=Hwfs*Hlag*Hint*Hdac*Hzoom;
	return Hol;
}

/**
   Sodium tracking Openloop transfer function. The cross over frequency
   (where gain is 1) have a phase margin of 45 (angle(Hol)=-135).

   First find the frequency fc where phase margin is 45. then scale Hol to have
   gain of 1 at fc.

   alpha, beta are the parameters of the sodium power spectrum obtained by
   fitting measurement data at UBC.

   Written: 2010-06-15
   Tested ok against RC's skycoverage code: 2010-06-15.

   @param[in] fs Sampling frequency of NGS (Hz).
   @param[in] tau Zoom corrector time delay (s).
   @param[in] zcf Zoom corrector frequency (Hz).
   @param[in] zeta Zoom corrector damping.
   @param[in] D Telescope diameter (m).
   @param[in] hs Guide star altitude (m).
   @param[in] alpha Parameter of sodium layer height PSD.
   @param[in] beta Parameter of sodium layer height PSD.
   @return RMS focus error (m) after sodium tracking.
*/
real nafocus_residual(real fs, real tau, real zcf, real zeta, real D, real hs, real alpha, real beta){
	int inu=1;
	real nu;
	real margin=M_PI/4;
	real aw=-(M_PI-margin);
	comp Hol;
	real angle;
	/*First we find out where NFIRAOS LGS zoom optics Hol have 45 degrees phase margin.*/
	while(1){
		nu=inu;
		Hol=nafocus_Hol(nu, fs, tau, zeta, zcf);
		angle=atan2(cimag(Hol), creal(Hol));
		if(angle>aw){
			inu++;
		} else{
			break;
		}
	}
	real nu1=nu-1;
	real nu2=nu;
	/*Then determine the true nu that lies in between nu1 and nu2. */
	comp Hol1=nafocus_Hol(nu1, fs, tau, zeta, zcf);/*>aw */
	comp Hol2=nafocus_Hol(nu2, fs, tau, zeta, zcf);/*<aw */
	real a1=atan2(cimag(Hol1), creal(Hol1));/*bigger than aw */
	real a2=atan2(cimag(Hol2), creal(Hol2));
	while(fabs(nu1-nu2)>1.e-5){
		nu=(nu1+nu2)/2;
		Hol=nafocus_Hol(nu, fs, tau, zeta, zcf);
		real am=atan2(cimag(Hol), creal(Hol));
		if(am>aw){
			nu1=nu;
		} else{
			nu2=nu;
		}
	}
	nu=nu1+(nu2-nu1)*(a1-aw)/(a1-a2);
	Hol=nafocus_Hol(nu, fs, tau, zeta, zcf);
	angle=atan2(cimag(Hol), creal(Hol));
	/*The gain is adjusted so that the cross over frequency (total gain of Hol
	  is 1) aligns with the location where we have 45 degrees phase margin*/
	real gain=1./cabs(Hol);
	/*We integrate over f, not f*2pi */
	dmat* nus=dlogspace(-3, 5, 2000);/*agrees good with skycoverage matlab code. */
	real rms=0;
	for(long i=0; i<nus->nx; i++){
		nu=P(nus,i);
		Hol=nafocus_Hol(nu, fs, tau, zeta, zcf);
		const comp Hrej=1./(1.+gain*Hol);
		const real NaPSD=nafocus_NaPSD(nu, alpha, beta);
		rms+=NaPSD*pow(cabs(Hrej), 2)*nu;/*we integratr f(nu)nu d(log(nu)) */
	}
	rms*=(log(P(nus,nus->nx-1))-log(P(nus,0)))/(nus->nx-1);
	real focus=sqrt(rms)*1./(16*sqrt(3))*pow((D/hs), 2);/*convert to focus error in meter. */
	dfree(nus);
	return focus;
}
/**
   Generate a random time series of sodium layer height with nstep points at
   time step dt, following the sodium layer height PSD defined by alpha and beta.

   @param[in] alpha Parameter of sodium layer height PSD.
   @param[in] beta Parameter of sodium layer height PSD.
   @param[in] dt Time step (s).
   @param[in] nstep Number of time steps.
   @param[in,out] rstat Random number generator state.
   @return Newly allocated nstep x 1 time series. The caller frees it.
*/
dmat* nafocus_time(real alpha, real beta, real dt, long nstep, rand_t* rstat){
	real df=1./(nstep*dt);
	cmat* psd=cnew(nstep, 1);
	//cfft2plan(psd, -1);
	for(int i=1; i<nstep; i++){
		P(psd,i)=sqrt(nafocus_NaPSD(df*i, alpha, beta)*df)*COMPLEX(randn(rstat), randn(rstat));
	}
	cfft2(psd, -1);
	dmat* out=NULL;
	creal2d(&out, 0, psd, 1);
	cfree(psd);
	return out;
}
