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
#ifndef AOS_LIB_NEA_H
#define AOS_LIB_NEA_H
#include "../math/mathdef.h"
/**
   \file nea.h
   Noise equivalent angle (NEA) of WFS gradients. The NEA of the nsa
   subapertures is stored as nsa x 2 (xx, yy) or nsa x 3 (xx, yy, xy) matrix,
   or nsa x ng for ng>2 (diagonal, for raw PWFS).
*/
dsp *nea2sp(dmat *nea, int ll, int ur, int ng);
void nea_chol(dmat **pout, const dmat *in, int ng);
void nea_inv(dmat **pout, const dmat *in, int ng, real scale);
void nea_mm(dmat **pout, const dmat *in, int ng);
void nea_check(dmat *nea, long nsa, int ng);
#endif
