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
/**
   \file dither_utils.h

   Contains utility routines to dithering or sodium profile fitting to update LGS WFS gain.
*/
#ifndef AOS_DITHER_UTILS_H
#define AOS_DITHER_UTILS_H
#include "common.h"
void gradoff_acc(sim_t* simu, int ipowfs);
void dither_position(real* cs, real* ss, int alfsm, int dtrat, int npoint, int isim, real deltam);
real dither_amp_calc(dmat **res, const dmat *signal, long dtrat, long npoint, int detrend, int combine);
void dither_acc_i0(const parms_t *parms, dither_t* pd, dcell *ints, int iwfs, int isim);
void dither_acc(sim_t* simu, int iwfs);
void dither_post(sim_t* simu);
void sodium_fit_wrap(dmat** psodium, dcell** pgrad, dcell** pi0, dcell** pgx, dcell** pgy, const dcell* i0in, 
  const parms_t* parms, powfs_t* powfs, int ipowfs, real r0, real L0, int nrep, int use_cache);
void sodium_fit(dmat** sodium, dcell** pgrad, dcell** pi0, dcell** pgx, dcell** pgy, const dcell* i0i, const dccell* sepsf, 
  const dtf_t* dtf, const loc_t* saloc, const dcell* saa, const dcell* srsa, const dcell* srot, const dmat* siglev, 
  const dmat* wvlwts, const dcell* gradoff, real dh, real hs, real htel, real za, real svdthres, int nrep, int save, int use_cache);
  
#endif