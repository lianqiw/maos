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

/*
  Contains functions used to calculate star flux.
*/
#include <math.h>
#include "skyc.h"
#include "photon.h"
/*
  Zero magnitude and sky background relocated to parameter file skyc.conf
*/
//const double Z_J=1.97e7*(20*9.56);
/*
  The following were obtain by R Clare on Table IV in Bessel 1988. But the print was garbelled. Fixed.
  const double Z_H=9.6e6*(20*16.5);
  const double Z_K=4.5e6*(20*14.41);
*/
/*
const double Z_H=9.6e6*(20*14.17);
const double Z_K=4.5e6*(20*16.74);


const double MB_J=(16.7+15.8)/2;
const double MB_H=(15.0+13.8)/2;
const double MB_K=(13.0+12.7)/2;//kshort
*/
/**
   \file skyc/photon.c
   Routines that calculates star flux from magnitude.
*/
/**
   Compute subaperture star flux from magnitude, subaperture size, sample
   period, etc.
 */
void photon_flux(const ZB_S *zb,        /**<[in] Sky background and zero magnitude flux*/
		 double *Np,      /**<[out] number of total signal at each wvl.*/
		 double *Nptot,   /**<[out] total signal (sum of Np).*/
		 double *Nbtot,   /**<[out] number of background photon per pixel*/
		 double *QCSNR,   /**<[out] signal to noise ratio for a Quadcell*/
		 double *QCNEA,   /**<[out] noise equivalent angle for a Quadcell with Nyquist sampling*/
		 int nwvl,        /**<[in] number of wavelength.*/ 
		 double* wvls,    /**<[in] vector of wavelength*/
		 double *mags,    /**<[in] vector of magnitudes*/
		 double dxsa,     /**<[in] subaperture side length*/
		 int iscircle,    /**<[in] where the subaperture is circle/part of a circle. true for TT/TTFA sensors.*/
		 double pixtheta, /**<[in] pixel extense in radian*/
		 double dt,       /**<[in] sampling period in seconds*/
		 double za,       /**<[in] zenith angle*/
		 double *strehl,  /**<[in] Strehl of the image. set to 1 for full flux.*/
		 double imperrnm, /**<[in] Implementation error in nm*/
		 double *thruput, /**<[in] end to end optical throughput*/
		 double *qe,      /**<[in] detector quantum efficiency.*/
		 double rne       /**<[in] detector read out noise.*/
		 ){
    /*
       Written 2010-06-09;
       Tested PASS 2010-06-09;
     */
    double pixas=pixtheta*206265;/*in arcsec. */
    double Npsum=0;
    double Nbsum=0;
    double saa;
    if(iscircle){
	saa=M_PI*pow(dxsa*0.5,2);
    }else{
	saa=pow(dxsa,2);
    }
    if(!Np){
	Np=alloca(sizeof(double)*nwvl);
    }
    double Npwvl=0;
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	double wvl=wvls[iwvl];
	double Z=0, ZB=0;
	double imperr_rad2=pow(imperrnm*1e-9*(2*M_PI/wvl),2);
	double imperr_strehl=exp(-imperr_rad2);
	if(fabs(wvl-1.25e-6)<1.e-7){ /*J band */
	    Z=zb->ZJ; ZB=pow(10,-zb->BJ/2.5)*Z;
	}else if(fabs(wvl-1.65e-6)<1.e-7){/*H band */
	    Z=zb->ZH; ZB=pow(10,-zb->BH/2.5)*Z;
	}else if(fabs(wvl-2.2e-6)<1.e-7){/*K band */
	    Z=zb->ZK; ZB=pow(10,-zb->BK/2.5)*Z;
	}else{
	    error("Invalid");
	}
	double absorp=(1./cos(za)-1)*(-log(0.98));/*atmosphere absorption. */
	double strehl_iwvl=1;
	if(strehl){
	    strehl_iwvl=strehl[iwvl];
	}
	Np[iwvl]=dt*saa*pow(10,-(absorp+mags[iwvl])/2.5)*Z
	    *thruput[iwvl]*qe[iwvl]*strehl_iwvl*imperr_strehl;
	Npsum+=Np[iwvl];
	Nbsum+=dt*saa*pow(pixas,2)*ZB*thruput[iwvl]*qe[iwvl];
	Npwvl+=Np[iwvl]/wvl;
    }
    double wvlm=Npsum/Npwvl; /*Average wavelength 1/mean(1/wvl) with signal weighting */
    double deltheta=wvlm/dxsa;
    double thetaB=3.*M_PI*deltheta/16.;
    double snr=Npsum/sqrt(Npsum+4*Nbsum+4.*pow(rne,2));
    if(Nptot) *Nptot=Npsum;
    if(Nbtot) *Nbtot=Nbsum;
    if(QCSNR) *QCSNR=snr;
    if(QCNEA) *QCNEA=thetaB/snr;
}
