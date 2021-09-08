/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file genseotf.h

   Contains routines to generate mean short exposure (tip/tilt
   removed) pixel intensities. Mostly used for LGS pixel intensity for its
   matched filter. Structure functions from kolmogorov spectrum is used. Not
   able to take into account outerscale yet.

   \todo find ways to factor in outerscale effect (use von karman spectrum
   instead of kolmogorov) 
*/
#ifndef AOS_GENSEOTF_H
#define AOS_GENSEOTF_H
#include "common.h"
void genmtch(const parms_t* parms, powfs_t* powfs, const int ipowfs);
void cog_nea(real* nea, const dmat* ints, real cogthres, real cogoff, int ntry,
  rand_t* rstat, real bkgrnd, real bkgrndc, const dmat* bkgrnd2i, const dmat* bkgrnd2ic, real rne
);
void fit_sodium_profile(
  dmat** sodium, /**<The sodium profile determined by fit*/
  dcell** pgrad, /**<The gradients determined by fit.*/
  dcell** pi0,   /**<The output i0*/
  dcell** pgx,   /**<The output gx*/
  dcell** pgy,   /**<The output gy*/
  const dcell* i0i, /**<The input i0*/
  const dccell* sepsf,   /**<Short exposure PSF*/
  const dtf_t* dtf,     /**<Detector transfer function*/
  const void* saa,      /**<Subaperture area. dmat or dcell*/
  const dcell* srsa,    /**<Subaperture to LLT distance*/
  const dcell* srot,    /**<Subaperture to LLT clocking*/
  const dmat* siglev,  /**<Subaperture signal level*/
  const dmat* wvlwts,    /**<Wavelength weights*/
  real dh,      /**<The sodium profile sampling in meters*/
  real hs,      /**<LGS focusing height*/
  real htel,    /**<Telescope hegith*/
  real za,      /**<Telescope zenith angle*/
  real tikcr,    /**<Tikhonov regularization*/
  real svdthres, /**<SVD threshold*/
  int save  /**<Save results*/
);
#endif
