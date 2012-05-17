/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
   Compute Zoom optics focus corrector open loop transfer function
*/
static inline dcomplex nafocus_Hol(double nu,  /**<[in] frequency.*/
				   double fs,  /**<[in] sampling frequency*/
				   double tau, /**<[in] time delay*/
				   double zeta,/**<[in] dampling of zoom corrector*/
				   double fzc  /**<[in] frequency of zoom corrector*/
				   ){
    const dcomplex s=2*M_PI*I*nu;
    const dcomplex fsovers=fs/s;
    const dcomplex Hint=fsovers;
    const dcomplex Hwfs=(1.-cexp(-1./fsovers))*fsovers;
    const dcomplex Hlag=cexp(-s*tau);
    const dcomplex Hdac=1;
    const double omega=fzc*2*M_PI;
    const dcomplex Hzoom=omega*omega/(s*s+2*zeta*omega*s+omega*omega);
    const dcomplex Hol=Hwfs*Hlag*Hint*Hdac*Hzoom;
    return Hol;
}
/**
   Compute sodium power spectrum density. alpha, beta are the parameters of the
   sodium power spectrum obtained by fitting measurement data at UBC. */
static inline double nafocus_NaPSD(double nu, double alpha, double beta){
    return pow(10,beta)*pow(nu,alpha);/*we don't divide 2pi */
}
/**
   Sodium tracking Openloop transfer function. The cross over frequency
   (where gain is 1) have a phase margin of 45 (angle(Hol)=-135).
       
   First find the frequency fc where phase margin is 45. then scale Hol to have
   gain of 1 at fc.

   Written: 2010-06-15
   Tested ok against RC's skycoverage code: 2010-06-15.
*/
double nafocus_residual(double fs,   /**<[in] sampling frequency of NGS*/
			double tau,  /**<[in] zoom corrector time delay*/
			double zcf,  /**<[in] zoom corrector frequency (Hz)*/
			double zeta, /**<[in] zoom corrector sampling*/
			double D,    /**<[in] telescope diameter */
			double hs,   /**<[in] guide star altitude*/
			double alpha,/**<[in] parameter of sodium layer height PSD.*/
			double beta  /**<[in] parameter of sodium layer height PSD.*/
			){
    int inu=1;
    double nu;
    double margin=M_PI/4;
    double aw=-(M_PI-margin);
    dcomplex Hol;
    double angle;
    /*First we find out where NFIRAOS LGS zoom optics Hol have 45 degrees phase margin.*/
    while(1){
	nu=inu;
	Hol=nafocus_Hol(nu, fs, tau, zeta, zcf);
	angle=atan2(cimag(Hol),creal(Hol));
	if(angle>aw){
	    inu++;
	}else{
	    break;
	}
    }
    double nu1=nu-1;
    double nu2=nu;
    /*Then determine the true nu that lies in between nu1 and nu2. */
    dcomplex Hol1=nafocus_Hol(nu1, fs, tau, zeta, zcf);/*>aw */
    dcomplex Hol2=nafocus_Hol(nu2, fs, tau, zeta, zcf);/*<aw */
    double a1=atan2(cimag(Hol1), creal(Hol1));/*bigger than aw */
    double a2=atan2(cimag(Hol2), creal(Hol2));
    while(fabs(nu1-nu2)>1.e-5){
	nu=(nu1+nu2)/2;
	Hol=nafocus_Hol(nu, fs, tau, zeta, zcf);
	double am=atan2(cimag(Hol), creal(Hol));
	if(am>aw){
	    nu1=nu;
	}else{
	    nu2=nu;
	}
    }
    nu=nu1+(nu2-nu1)*(a1-aw)/(a1-a2);
    Hol=nafocus_Hol(nu, fs, tau, zeta, zcf);
    angle=atan2(cimag(Hol), creal(Hol));
    /*The gain is adjusted so that the cross over frequency (total gain of Hol
      is 1) aligns with the location where we have 45 degrees phase margin*/
    double gain=1./cabs(Hol);
    /*We integrate over f, not f*2pi */
    dmat *nus=dlogspace(-3,5,2000);/*agrees good with skycoverage matlab code. */
    double rms=0;
    for(long i=0; i<nus->nx; i++){
	nu=nus->p[i];
	Hol=nafocus_Hol(nu, fs, tau, zeta, zcf);
	const dcomplex Hrej=1./(1.+gain*Hol);
	const double NaPSD=nafocus_NaPSD(nu, alpha, beta);
	rms+=NaPSD*pow(cabs(Hrej),2)*nu;/*we integratr f(nu)nu d(log(nu)) */
    }
    rms*=(log(nus->p[nus->nx-1])-log(nus->p[0]))/(nus->nx-1);
    double focus=sqrt(rms)*1./(16*sqrt(3))*pow((D/hs),2);/*convert to focus error in meter. */
    dfree(nus);
    return focus;
}
