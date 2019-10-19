/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_LIB_LAPLACIAN_C
#define AOS_LIB_LAPLACIAN_C
#include "../math/mathdef.h"
/**
   \file laplacian.h
   Compute laplacian operator for turbulence covariance approximation.
*/
void apply_laplacian_map(dmat *opdout, const dmat *opd, real dx, real r0, real weight);
dsp* mklaplacian_map(int nx, int ny, real dx, real r0, real weight);
dsp* mklaplacian_loc(loc_t *loc, real r0, real weight);
real laplacian_coef(real r0, real weight, real dx);
#endif
