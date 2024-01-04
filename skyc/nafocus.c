/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
static inline real nafocus_NaPSD(real nu, real alpha, real beta2){
	return beta2*pow(nu, alpha);/*we don't divide 2pi */
}

/**
   Compute Zoom optics focus corrector open loop transfer function
*/
static inline comp nafocus_Hol(real nu,  /**<[in] frequency.*/
	real fs,  /**<[in] sampling frequency*/
	real tau, /**<[in] time delay*/
	real zeta,/**<[in] dampling of zoom corrector*/
	real fzc  /**<[in] frequency of zoom corrector*/
){
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

   Written: 2010-06-15
   Tested ok against RC's skycoverage code: 2010-06-15.
*/
real nafocus_residual(real fs,   /**<[in] sampling frequency of NGS*/
	real tau,  /**<[in] zoom corrector time delay*/
	real zcf,  /**<[in] zoom corrector frequency (Hz)*/
	real zeta, /**<[in] zoom corrector sampling*/
	real D,    /**<[in] telescope diameter */
	real hs,   /**<[in] guide star altitude*/
	real alpha,/**<[in] parameter of sodium layer height PSD.*/
	real beta  /**<[in] parameter of sodium layer height PSD.*/
){
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
dmat* nafocus_time(real alpha,/**<[in] parameter of sodium layer height PSD.*/
	real beta, /**<[in] parameter of sodium layer height PSD.*/
	real dt,   /**<[in] time step*/
	long nstep,/**<[in] number of steps*/
	rand_t* rstat /**<random number stat*/){
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
