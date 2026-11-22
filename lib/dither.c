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
#include "dither.h"
/**
   Calculate expected averaged position of dithering signal during WFS integration. Called when (isim+1)%dtrat=0

   @param[out] cs The cosine of the dither angle, scaled to account for averaging.
   @param[out] ss The sine of the dither angle, scaled to account for averaging.
   @param[in] alfsm The additional latency due to uplink (in simulation steps).
   @param[in] dtrat The WFS sampling period in simulation steps.
   @param[in] npoint The number of WFS frames per dithering cycle.
   @param[in] isim The simulation step index.
   @param[in] deltam The dither phase offset (radian).
*/
void dither_position(real* cs, real* ss, int alfsm, int dtrat, int npoint, int isim, real deltam){
	//adjust for delay due to propagation, and computation delay if delay is not 2 frame.
	const int adjust=alfsm+1-dtrat;
	//adjust to get delay at beginning of integration
	const int adjust2=dtrat-1;
	const real anglei=(2*M_PI/npoint);
	const real angle0=M_PI*0.5;
	const real angle=angle0+((isim-adjust-adjust2)/dtrat)*anglei+deltam;
	const real angle2=angle0+((isim-adjust)/dtrat)*anglei+deltam;
	const real delay=(real)adjust/dtrat;
	const real beta=1+delay+floor(-delay);
	const real scale=1./(beta*beta+(1-beta)*(1-beta));
	//use average of two places during accumulation and scale
	*cs=(beta*cos(angle)+(1-beta)*cos(angle2))*scale;
	*ss=(beta*sin(angle)+(1-beta)*sin(angle2))*scale;
}
/**
   Demodulate the dithering signal to determine the amplitude. Remove trend (detrending) if detrend is set.

   @param[out] res The amplitude of each mode (nmod*1). Only used if not combine. Required if nmod>1. May be NULL.
   @param[in] signal The signal matrix. nmod*nsim, or a column vector for a single mode.
   @param[in] dtrat skip columns due to wfs/sim dt ratio
   @param[in] npoint The number of WFS frames per dithering cycle.
   @param[in] detrend Whether to detrend the signal.
   @param[in] combine Whether to combine the results. for tip/tilt only
   @return The calculated dither amplitude (of the first mode if not combine).
*/
real dither_amp_calc(dmat **res, const dmat *signal, long dtrat, long npoint, int detrend, int combine){
	long nmod=NY(signal)==1?1:NX(signal);
	long nframe=((NY(signal)==1?signal->nx:signal->ny)-1)/dtrat+1;//number of actual frames
	dmat *slope=dnew(nmod, 1);//for detrending
	long offset=(nframe/npoint-1)*npoint;//number of WFS frame separations between first and last cycle
	if(detrend&&offset){//detrending
		for(long imod=0; imod<nmod; imod++){
			for(long ip=0; ip<npoint; ip++){
				long i0=ip*dtrat*nmod+imod;
				long i1=(ip+offset)*dtrat*nmod+imod;
				P(slope,imod)+=P(signal, i1)-P(signal, i0);
			}
			P(slope,imod)/=(npoint*offset);
		}
		//dbg("slope=%g. npoint=%ld, nmod=%ld, nframe=%ld, offset=%ld\n", slope, npoint, nmod, nframe, offset);
	}
	real anglei=M_PI*2/npoint;
	real angle0=M_PI*0.5;//2023-05-29: bias added to work with npoint==2
	real ipv=0, qdv=0, a2m=0;
	if(combine){//tip and tilt dithering
		if(nmod!=2){
			error("combine only support nmod=2\n");
		}
		for(int iframe=0; iframe<nframe; iframe++){
			real angle=angle0+anglei*iframe;//position of dithering
			real cs=cos(angle);
			real ss=sin(angle);
			real ttx=P(signal, iframe*dtrat*nmod)-P(slope,0)*iframe;
			real tty=P(signal, iframe*dtrat*nmod+1)-P(slope,1)*iframe;
			ipv+=(ttx*cs+tty*ss);
			qdv+=(ttx*ss-tty*cs);
		}
		a2m=sqrt(ipv*ipv+qdv*qdv)/nframe;
		
	}else{//independent mode
		if(nmod>1&&!res){
			error("when there are more than 1 mode. res must be used to return the results\n");
		}
		if(res){
			dinit(res, nmod, 1);
		}
		dmat *ipq=dnew(2,nmod);
		for(int iframe=0; iframe<nframe; iframe++){
			real angle=angle0+anglei*iframe;//position of dithering
			real cs=cos(angle);
			real ss=sin(angle);
			for(int imod=0; imod<nmod; imod++){
				real mod=P(signal, iframe*dtrat*nmod+imod)-P(slope,imod)*iframe;
				P(ipq,0,imod)+=(mod*cs);//ipv
				P(ipq,1,imod)+=(mod*ss);//iqv
			}
		}
		for(int imod=0; imod<nmod; imod++){
			a2m=sqrt(P(ipq,0,imod)*P(ipq,0,imod)+P(ipq,1,imod)*P(ipq,1,imod))/nframe*2.;
			if(res) P(*res, imod)=a2m;
		}
		if(res){
			a2m=P(*res,0);//first mode
		}
		dfree(ipq);
	}
	dfree(slope);
	return a2m;
}
