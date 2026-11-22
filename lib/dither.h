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
#ifndef AOS_LIB_DITHER_H
#define AOS_LIB_DITHER_H
#include "../math/mathdef.h"
/**
   \file dither.h
   Position of a circular dithering and the estimation of the dithering amplitude
   from the measured signal.
*/
void dither_position(real* cs, real* ss, int alfsm, int dtrat, int npoint, int isim, real deltam);
real dither_amp_calc(dmat **res, const dmat *signal, long dtrat, long npoint, int detrend, int combine);
#endif
