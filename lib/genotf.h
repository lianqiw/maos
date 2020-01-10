/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_LIB_GENOTF_H
#define AOS_LIB_GENOTF_H
#include "../math/mathdef.h"
/**
   \file genotf.h
   Routines to generate short/long exposure OTFs of an aperture in present of
   atmosphere turbulence.
*/
void genotf(cmat **otf,    /**<The otf array for output*/
	    loc_t *loc,    /**<the common aperture grid*/
	    const dmat *amp,     /**<The amplitude map of all the (sub)apertures*/
	    const dmat *opdbias, /**<The static OPD bias. */
	    const dmat *area,    /**<normalized area of the (sub)apertures*/
	    real thres,  /**<The threshold to consider a (sub)aperture as full*/
	    real wvl,    /**<The wavelength. only needef if opdbias is not null*/
	    real dtheta, /**<Sampling of PSF.*/
	    const dmat *cov,/**<The covariance. If not supplied use r0 for kolmogorov spectrum.*/
	    real r0,     /**<Fried parameter*/
	    real l0,     /**<Outer scale*/
	    long notfx,   /**<Size of OTF*/
	    long notfy,   /**<Size of OTF*/
	    long nsa,      /**<Number of (sub)apertures*/
	    long pttr      /**<Remove piston/tip/tilt*/
	    );

cell *genotf2(loc_t *loc, const dmat *amp, const dmat *opdbias, const dmat *area, real thres, real wvl, real dtheta, const dmat *cov, real r0, real l0, long notfx, long notfy, long nsa, long pttr);

dmat* mk2dcov(loc_t *loc, const dmat *amp, real ampthres, const dmat *cov, int norm);
void mtch(dmat **mtche,   /**<[out] the matched filter*/
	  dmat **neaout,  /**<[out] the subaperture noise equivalent angle*/
	  const dmat *i0, /**<[in] Averaged subaperture image*/
	  const dmat *gx, /**<[in] derivative of i0 along x (r)*/
	  const dmat *gy, /**<[in] derivative of i0 along y (a)*/
	  const dmat *qe, /**<[in] non uniform quantum efficiency (optional)*/
	  const dmat *dbkgrnd2,  /**<[in] background*/
	  const dmat *dbkgrnd2c, /**<[in] background calibration*/
	  real bkgrnd,    /**<[in] global background*/
	  real bkgrndc,   /**<[in] global background calibration*/
	  real rne,       /**<[in] Detector read noise*/
	  real pixthetax, /**<[in] Size of pixel along x*/
	  real pixthetay, /**<[in] Size of pixel along y*/
	  real pixrot,    /**<[in] Rotation (CCW, radian) of pixel island 0 for cartesian*/
	  int radgx,      /**<[in] 1: gx/gy is along r/a coord.*/
	  int cr          /**<Constraint flag 0: disable, 1: both axis, 2: x only, 3: y only*/
    );
void mtch2(dmat **mtche, dmat **nea, const dmat *i0, const dmat *gx, const dmat *gy, int cr);
dmat *derive_by_fft(const dmat *i0, real theta);

#endif
