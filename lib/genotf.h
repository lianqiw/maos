/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "type.h"
#include "loc.h"
/**
   \file genotf.h
   Routines to generate short/long exposure OTFs of an aperture in present of
   atmosphere turbulence.
*/
void genotf(cmat **otf,    /**<The otf array for output*/
	    loc_t *loc,    /**<the common aperture grid*/
	    const double *amp,     /**<The amplitude map of all the (sub)apertures*/
	    const double *opdbias, /**<The static OPD bias. */
	    const double *area,    /**<normalized area of the (sub)apertures*/
	    double thres,  /**<The threshold to consider a (sub)aperture as full*/
	    double wvl,    /**<The wavelength. only needef if opdbias is not null*/
	    double dtheta, /**<Sampling of PSF.*/
	    const dmat *cov,/**<The covariance. If not supplied use r0 for kolmogorov spectrum.*/
	    double r0,     /**<Fried parameter*/
	    double l0,     /**<Outer scale*/
	    long ncompx,   /**<Size of OTF*/
	    long ncompy,   /**<Size of OTF*/
	    long nsa,      /**<Number of (sub)apertures*/
	    long pttr      /**<Remove piston/tip/tilt*/
	    );
void mk2dcov(dmat **cov2d, loc_t *loc, const double *amp, double ampthres, const dmat *cov, int norm);
#endif
