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
#ifndef AOS_LIB_SHWFS_H
#define AOS_LIB_SHWFS_H
#include "loc.h"
/**
   \file shwfs.h
   Geometry of the subapertures of a Shack-Hartmann WFS: the subaperture grid,
   the points in each subaperture, the lenslet fill factor, the subaperture area,
   the selection of subapertures, the geometric gradients, and the center of
   gravity (CoG) mask, noise, and gradient offset.
*/
loc_t* shwfs_saloc(int order, real dsa, real din, real misreg, real misreg_aper);
pts_t* shwfs_pts(const loc_t* saloc, real dsa, real dx);
dmat* shwfs_safill(int nxsa, real safill2d);
dmat* shwfs_saa(const dmat* amp, long npsa);
dmat* shwfs_saa_max(dcell* saa, real alpha);
void shwfs_sa_isolated(pts_t* pts, dmat* saa, real saat);
dspcell* shwfs_mkg(loc_t* loc, dcell* amp, loc_t* saloc, dcell* saa, real saat);
dccell* shwfs_saimcc(const pts_t* pts, const dcell* amp);
void shwfs_grad_geom(dmat** grad, const pts_t* pts, const dmat* opd, const dmat* amp, const dsp* GS0, const dcell* imcc);
dcell* shwfs_cog_mask(const dcell* i0, real thres);
dmat* shwfs_cog_nea(dmat** i0, dmat** cogmask, dmat** bkgrnd2, dmat** bkgrnd2c, const real* srot, long nsa, real cogthres, real cogoff, int ntry, rand_t* rstat, real bkgrnd, real bkgrndc, real rne, real pixthetax, real pixthetay);
void shwfs_cog_gradoff(dmat* gradoff, dmat** i0, dmat** cogmask, long nsa, real cogthres, real cogoff, real pixthetax, real pixthetay);
#endif
