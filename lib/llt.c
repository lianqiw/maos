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
#include "llt.h"
#include "zernike.h"

/**
   Create the points of the LLT pupil: one subaperture with the same sampling
   as the WFS subapertures, so that the uplink OTF has the same sampling, and
   at least as large as the WFS subaperture. It is centered on the LLT.

   @param[in] lltd Diameter of the LLT in meter.
   @param[in] dsa Size of the WFS subaperture in meter.
   @param[in] dx Sampling of the WFS subaperture in meter.
   @return The points. Newly allocated; the caller frees it.
*/
pts_t* llt_pts(real lltd, real dsa, real dx){
	real lltdsa=MAX(lltd, dsa);
	int nxsa=round(lltdsa/dx); lltdsa=dx*nxsa;
	pts_t* pts=ptsnew(1, lltdsa, lltdsa, nxsa, nxsa, dx, dx);
	pts->origx[0]=(dx-lltdsa)*0.5;
	pts->origy[0]=(dx-lltdsa)*0.5;
	return pts;
}
/**
   Create the amplitude of the laser beam on the LLT pupil: a Gaussian whose
   waist (the radius where the amplitude drops to 1/e) is widthp*lltd/2,
   clipped by the LLT aperture of diameter lltd. The misregistration shifts the
   beam relative to the aperture, which clips it. The amplitude is normalized
   so that max(otf)=1 for the unclipped beam.

   @param[in] pts The points of the LLT pupil (llt_pts()).
   @param[in] lltd Diameter of the LLT in meter.
   @param[in] widthp Waist of the beam relative to lltd/2.
   @param[in] misregx Misregistration of the beam along x in meter.
   @param[in] misregy Misregistration of the beam along y in meter.
   @return Amplitude on the points of pts. Newly allocated; the caller frees it.
*/
dmat* llt_amp(const pts_t* pts, real lltd, real widthp, real misregx, real misregy){
	const long nxsa=pts->nxsa;
	const real dx=pts->dx;
	dmat* amp=dnew(nxsa*nxsa, 1);
	real sumamp2=0;
	real* amps=P(amp);
	const real l2max=pow(lltd*0.5, 2);
	/*the waist is defined as the radius where amplitude
	  drop to 1/e or intensity to 1/e^2.*/
	real r2waist=pow(lltd*0.5*widthp, 2);
	//misreg is treated as pupil centering error. It clips the laser
	const real ox=pts->origx[0];
	const real oy=pts->origy[0];
	for(int iy=0; iy<nxsa; iy++){
		real yy=iy*dx+oy;
		for(int ix=0; ix<nxsa; ix++){
			real xx=ix*dx+ox;
			real r2=xx*xx+yy*yy;
			real r2m=pow(xx-misregx, 2)+pow(yy-misregy, 2);
			real ampi=exp(-r2m/r2waist);
			sumamp2+=pow(ampi, 2);
			if(r2<=l2max){
				amps[iy*nxsa+ix]=ampi;
			}
		}
	}
	/*normalized so that max(otf)=1 for unclipped beam; */
	sumamp2=1./(sqrt(sumamp2));
	dscale(amp, sumamp2);
	return amp;
}
/**
   Add focus with the given RMS (over the points where amp>0, piston removed)
   to the LLT OPD.

   @param[in,out] ncpa OPD of the LLT for each WFS (or one for all). Each element must be allocated.
   @param[in] loc The points of the LLT pupil.
   @param[in] amp Amplitude on loc.
   @param[in] lltd Diameter of the LLT in meter.
   @param[in] focus RMS focus in nm.
*/
void llt_add_focus(dcell* ncpa, const loc_t* loc, const dmat* amp, real lltd, real focus){
	real nm2rad=1e-9*2*sqrt(3)*pow(lltd, -2);
	dmat* opd=dnew(NX(P(ncpa, 0)), NY(P(ncpa, 0)));
	loc_add_focus(opd, loc, focus*nm2rad);
	real var=0, piston=0;
	long count=0;
	for(long i=0; i<loc->nloc; i++){
		if(P(amp, i)>0){
			count++;
			var+=P(opd, i)*P(opd, i);
			piston+=P(opd, i);
		}
	}
	var/=count;
	piston/=count;
	dadds(opd, -piston);
	var=sqrt(var-piston*piston);
	for(int ic=0; ic<NX(ncpa); ic++){
		dadd(&P(ncpa, ic), 1, opd, focus*1e-9/var);
	}
	dfree(opd);
}
/**
   Remove tip/tilt (and focus) from the LLT OPD, weighted by the amplitude.

   @param[in,out] ncpa OPD of the LLT for each WFS (or one for all).
   @param[in] loc The points of the LLT pupil.
   @param[in] amp Amplitude on loc.
   @param[in] lltd Diameter of the LLT in meter.
   @param[in] ttfr 1: remove tip/tilt. 2: remove tip/tilt and focus.
   @param[out] ppttf The modes removed (with piston). Not returned if NULL.
   @param[out] pproj The weighted projection to the modes. Not returned if NULL.
*/
void llt_remove_ttf(dcell* ncpa, const loc_t* loc, const dmat* amp, real lltd, int ttfr, dmat** ppttf, dmat** pproj){
	int nmod=ttfr==2?4:3;
	dmat* pttfa=zernike(loc, lltd, 0, 2, 0);
	dmat* pttf=drefcols(pttfa, 0, nmod);
	dmat* proj=dpinv(pttf, amp);
	dmat* res=dnew(nmod, 1);
	for(int ilotf=0; ilotf<NX(ncpa)*NY(ncpa); ilotf++){
		dzero(res);
		dmulvec(P(res), proj, P(P(ncpa, ilotf)), 1);
		dmulvec(P(P(ncpa, ilotf)), pttf, P(res), -1);
	}
	if(ppttf){
		*ppttf=pttf;
	} else{
		dfree(pttf);
	}
	if(pproj){
		*pproj=proj;
	} else{
		dfree(proj);
	}
	dfree(pttfa);
	dfree(res);
}
