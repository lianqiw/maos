/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file mtch.h
*/

#ifndef AOS_LIB_MTCH_H
#define AOS_LIB_MTCH_H
#include "../math/mathdef.h"

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
void mtch_cell(
	dcell** pmtche,  /**<[out] matched filter. */
	dcell** psanea,  /**<[out] subaperture noise equivalent angle*/
	dmat** pi0sum,   /**<[out] sum of subaperture intensity*/
	dmat** pi0sumsum,/**<[out] sum of all subaperture intensity*/
	const dcell* i0s, /**<Subaperture intensities*/
	const dcell* gxs, /**<Subaperture intensity gradients along x (radial)*/
	const dcell* gys, /**<Subaperture intensity gradients along y (azimuthal)*/
	const dmat* qe, /**<Quantum efficiency of each pixel*/
	const dcell* bkgrnd2, /**<bkgrnd image*/
	const dcell* bkgrnd2c,/**<bkgrnd correction image*/
	real bkgrnd,  /**<bkgrnd per pixel*/
	real bkgrndc, /**<bkgrnd correction per pixel*/
	real rne,     /**<Read out noise per pixel*/
	real pixthetax, /**<Pixel size along x (radial)*/
	real pixthetay, /**<Pixel size along y (azimuthal)*/
	dcell* pixrots,   /**<subaperture rotation*/
	int radgx,     /**<Leave gradients at radial direction */
	int mtchcr,     /**<constraint. -1: auto*/
	real sigratio  /**<scale signal level to increase NEA (testing)*/
);
dmat *derive_by_fft(const dmat *i0, real theta);
#endif
