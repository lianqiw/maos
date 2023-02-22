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
real remove_focus_grad(const loc_t* saloc, dmat* grad, real factor);
void sodium_fit(
  dmat** sodium, /**<The sodium profile determined by fit*/
  dcell** pgrad, /**<The gradients determined by fit.*/
  dcell** pi0,   /**<The output i0*/
  dcell** pgx,   /**<The output gx*/
  dcell** pgy,   /**<The output gy*/
  const dcell* i0i, /**<The input i0*/
  const dccell* sepsf,   /**<Short exposure PSF*/
  const dtf_t* dtf,     /**<Detector transfer function*/
	const loc_t* saloc,   /**<Saloc*/
  const dcell* saa,      /**<Subaperture area. dmat or dcell*/
  const dcell* srsa,    /**<Subaperture to LLT distance*/
  const dcell* srot,    /**<Subaperture to LLT clocking*/
  const dmat* siglev,  /**<Subaperture signal level*/
  const dmat* wvlwts,    /**<Wavelength weights*/
  const dcell* gradncpa,/**<NCPA gradient to be used for pi0,pgx,pgy output.*/
  real dh,      /**<The sodium profile sampling in meters*/
  real hs,      /**<LGS focusing height*/
  real htel,    /**<Telescope hegith*/
  real za,      /**<Telescope zenith angle*/
  real svdthres, /**<SVD threshold*/
  int nrep,     /**<Number of iterations*/
  int save,      /**<Save results to file*/
  int use_cache  /**<Use cache*/
);
void sodium_fit_wrap(dmat** psodium, /**<[out] sodium profile*/
  dcell** pgrad, /**<[out] estimated actual gradient*/
  dcell** pi0,   /**<[out] The output i0*/
  dcell** pgx,   /**<[out] The output gx*/
  dcell** pgy,   /**<[out] The output gy*/
  const dcell* i0in, /**<[in]The input sa intensities. may equal to *pi0 */
  const parms_t* parms,/**<[in]parms*/
  powfs_t* powfs, /**<[in]powfs*/
  int ipowfs, /**<[in] ipowfs*/
  real r0,  /**<[in] Fried parameter*/
  real L0,  /**<[in] outer scale*/
  int nrep, /**<[in] Number of iterations. 1 for mtche, 3 for cog*/
  int use_cache /**<[in] cache intermediate results.*/
);
#endif
