/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "skyc.h"
#include "photon.h"
/*
  Zero magnitude and sky background relocated to parameter file skyc.conf
*/
//const real Z_J=1.97e7*(20*9.56);
/*
  The following were obtain by R Clare on Table IV in Bessel 1988. But the print was garbelled. Fixed.
  const real Z_H=9.6e6*(20*16.5);
  const real Z_K=4.5e6*(20*14.41);
*/
/*
const real Z_H=9.6e6*(20*14.17);
const real Z_K=4.5e6*(20*16.74);

%these are for LA Silla from ESO slide
const real MB_J=(16.7+15.8)/2;
const real MB_H=(15.0+13.8)/2;
const real MB_K=(13.0+12.7)/2;//kshort
*/
/**
   \file skyc/photon.c
   Routines that calculates star flux from magnitude.
*/
/**
   Compute subaperture star flux from magnitude, subaperture size, sample
   period, etc.
 */
void photon_flux(const ZB_S* zb,        /**<[in] Sky background and zero magnitude flux*/
	real* Np,      /**<[out] number of total signal at each wvl.*/
	real* Nptot,   /**<[out] total signal (sum of Np).*/
	real* Nbtot,   /**<[out] number of background photon per pixel*/
	real* QCSNR,   /**<[out] signal to noise ratio for a Quadcell*/
	real* QCNEA,   /**<[out] noise equivalent angle for a Quadcell with Nyquist sampling*/
	int nwvl,        /**<[in] number of wavelength.*/
	real* wvls,    /**<[in] vector of wavelength*/
	real* mags,    /**<[in] vector of magnitudes*/
	real dxsa,     /**<[in] subaperture side length*/
	int iscircle,    /**<[in] where the subaperture is circle/part of a circle. true for TT/TTFA sensors.*/
	real pixtheta, /**<[in] pixel extense in radian*/
	real dt,       /**<[in] sampling period in seconds*/
	real za,       /**<[in] zenith angle*/
	real* strehl,  /**<[in] Strehl of the image. set to 1 for full flux.*/
	real imperrnm, /**<[in] Implementation error in nm*/
	real rne       /**<[in] detector read out noise.*/
){
/*
   Written 2010-06-09;
   Tested PASS 2010-06-09;
 */
	real pixas=pixtheta*RAD2AS;/*in arcsec. */
	real Npsum=0;
	real Nbsum=0;
	real saa;
	real Nptmp[nwvl];
	if(iscircle){
		saa=M_PI*pow(dxsa*0.5, 2);
	} else{
		saa=pow(dxsa, 2);
	}
	if(!Np){
		Np=Nptmp;
	}
	real Npwvl=0;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		real wvl=wvls[iwvl];
		real imperr_rad2=pow(imperrnm*1e-9*(2*M_PI/wvl), 2);
		real imperr_strehl=exp(-imperr_rad2);
		
		int jwvl;
		for(jwvl=0; jwvl<NX(zb->wvl); jwvl++){
			if(fabs(wvl-P(zb->wvl, jwvl))<0.1e-6){
				break;
			}
		}
		if(jwvl==NX(zb->wvl)){
			error("wvl=%g is not configured\n", wvl);
			continue;
		}
		real Z=P(zb->Z, jwvl);
		real ZB=pow(10, -P(zb->B, jwvl)/2.5)*Z;
		real thruputqe=P(zb->thruput, jwvl)*P(zb->qe, jwvl);
		real excessbkgrnd=P(zb->excessbkgrnd, jwvl);
		
		real absorp=(1./cos(za)-1)*(-log(0.98));/*atmosphere absorption. */
		real strehl_iwvl=1;
		if(strehl){
			strehl_iwvl=strehl[iwvl];
		}
		Np[iwvl]=dt*saa*pow(10, -(absorp+mags[iwvl])/2.5)*Z*thruputqe*strehl_iwvl*imperr_strehl;
		Npsum+=Np[iwvl];
		Nbsum+=dt*saa*pow(pixas, 2)*ZB*thruputqe*(1.+excessbkgrnd)/cos(za);//background scale with sec(za). 2012-10-31.
		Npwvl+=Np[iwvl]/wvl;
	}
	//warning_once("How does background scale with zenith angle? Surface brightness is independent of distance\n");
	real wvlm=Npsum/Npwvl; /*Average wavelength 1/mean(1/wvl) with signal weighting */
	real deltheta=wvlm/dxsa;
	real thetaB=3.*M_PI*deltheta/16.;
	real snr=Npsum/sqrt(Npsum+4*Nbsum+4.*pow(rne, 2));
	if(Nptot) *Nptot=Npsum;
	if(Nbtot) *Nbtot=Nbsum;
	if(QCSNR) *QCSNR=snr;
	if(QCNEA) *QCNEA=thetaB/snr;
}
