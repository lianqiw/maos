/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file gensei.h

   Contains routines to generate mean short exposure (tip/tilt
   removed) pixel intensities. Mostly used for LGS pixel intensity for its
   matched filter. Structure functions from kolmogorov spectrum is used. Not
   able to take into account outerscale yet.

   \todo find ways to factor in outerscale effect (use von karman spectrum
   instead of kolmogorov)
*/
#ifndef AOS_LIB_GENSEI_H
#define AOS_LIB_GENSEI_H
#include "../math/mathdef.h"
#include "mkdtf.h"
void upsample_otf(cmat* out, const cmat* in);
cccell* genseotf(const pts_t* pts, const_anyarray amp, const dcell* opdbias,
  const_anyarray saa, const dmat* wvl, real r0, real L0, int embfac);
void gensepsf(dccell** psepsfs, const cccell* otfs, const cccell* lotf, 
  const_anyarray saa, const dmat* wvl, int notfx, int notfy);
void gensei(dcell** pi0, dcell** pgx, dcell** pgy, cccell** pfotf, 
  const dccell* sepsfs, const dtf_t* dtf, const etf_t* etf,
  const dcell* saa, const dcell* gxyrot, const dmat* siglev, const dmat* wvlwts, const dcell *goff,
  int i0scale, int shift2center);
#endif
