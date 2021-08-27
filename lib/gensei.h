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
cccell* genseotf(const pts_t* pts, const void* amp, const dcell* opdbias,
	const void* saa, const dmat* wvl, real r0, real L0, int embfac);
//void genseotf(const parms_t *parms, powfs_t *powfs, int ipowfs);
void gensepsf(dccell** psepsfs, const cccell* otfs, const cccell* lotf, const void* saa, dmat* wvl, int notfx, int notfy);
void gensei(dcell** pi0, dcell** pgx, dcell** pgy, cccell** pfotf, cccell** ppotf,
  dccell* sepsfs, dtf_t* dtf, etf_t* etf, dcell* saa, dcell* srot, dmat* siglev, dmat* wvlwts, int dtrat,
  int i0scale, int radgx, int shift2center
);
#endif
