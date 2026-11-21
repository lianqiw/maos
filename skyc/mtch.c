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
#include "skyc.h"
#include "mtch.h"
/**
   \file skyc/mtch.c
   Routings regarding matched filter computation.
*/
/**
   Compute pixel intensity i0, and gradient of pixel intensity gx, gy from the
   PSF. The results are accumulated into i0, gx, gy.
   \todo merge this file with maos/mtch.c and put in lib folder.

   @param[in,out] i0 Pixel intensity, accumulated.
   @param[in,out] gx Gradient of pixel intensity along x, accumulated. Used only if deriv is set.
   @param[in,out] gy Gradient of pixel intensity along y, accumulated. Used only if deriv is set.
   @param[in] psf PSF with peak in the corner.
   @param[in] dtf Detector transfer function and sampling matrix.
   @param[in] deriv Whether to compute gx and gy.
*/
void psf2i0gxgy(dmat* i0, dmat* gx, dmat* gy, dmat* psf, dtf_s* dtf, int deriv){
	cmat* otf=cnew(psf->nx, psf->ny);
	cmat* otfsave=0;
	ccpd(&otf, psf);//loaded psf has peak in corner 
	cfft2i(otf, -1);//turn to OTF, peak in corner.
	ccwm(otf, dtf->nominal);//apply pixel transfer function
	if(deriv){
		ccp(&otfsave, otf);//save for later
	}
	cfft2(otf, 1);//convert otf back to psf space
	dspmulcreal(P(i0), dtf->si, P(otf), 1);//sample psf to detectors.
	if(deriv){
		ccp(&otf, otfsave); //copy back otf
		//apply derivative.
		for(int iy=0; iy<otf->ny; iy++){
			for(int ix=0; ix<otf->nx; ix++){
				P(otf, ix, iy)*=P(dtf->U,ix);
				P(otfsave, ix, iy)*=P(dtf->U,iy);
			}
		}
		cfft2(otf, 1);
		cfft2(otfsave, 1);
		dspmulcreal(P(gx), dtf->si, P(otf), 1);
		dspmulcreal(P(gy), dtf->si, P(otfsave), 1);
		cfree(otfsave);
	}
	cfree(otf);
}

/**
   Compute matched filter mtche and noise equivalent angle sanea (rms, x for all
   subapertures followed by y) for each subaperture from i0, gx, gy with read out
   noise rne and background bkgrnd.

   @param[in,out] mtche Matched filter for each subaperture. Allocated if NULL.
   @param[in,out] sanea Noise equivalent angle (rad), x for all subapertures followed by y. Allocated if NULL.
   @param[in] i0 Pixel intensity for each subaperture.
   @param[in] gx Gradient of i0 along x for each subaperture.
   @param[in] gy Gradient of i0 along y for each subaperture.
   @param[in] pixtheta Pixel size (rad).
   @param[in] rne Read out noise.
   @param[in] bkgrnd Background per pixel.
   @param[in] cr Constraint flag passed to mtch().
*/
void genmtch(dcell** mtche, dmat** sanea,
	dcell* i0, dcell* gx, dcell* gy, real pixtheta,
	real rne, real bkgrnd, int cr){
	const long nsa=i0->nx;
	if(!*mtche){
		*mtche=dcellnew(nsa, 1);
	}
	if(!*sanea){
		*sanea=dnew(nsa*2, 1);
	}
	for(long isa=0; isa<nsa; isa++){
		dmat* nea2=0;
		mtch(&P(*mtche, isa), &nea2, P(i0, isa), P(gx, isa), P(gy, isa), 0, 0, 0,
			bkgrnd, bkgrnd, rne, pixtheta, pixtheta, 0, 0, cr);
	   /*Drop coupling in x/y gradients. */
		P(*sanea,isa)=sqrt(P(nea2,0));
		P(*sanea,isa+nsa)=sqrt(P(nea2,3));
		dfree(nea2);
	}
}
