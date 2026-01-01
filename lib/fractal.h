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
#ifndef AOS_LIB_FRACTAL_H
#define AOS_LIB_FRACTAL_H
#include "../math/mathdef.h"
/**
   \file fractal.h

   Implementation of the fractal operation for atmospheric turbulence screen
   generation and reconstruction.
*/
#define ARGS dmat *p0, real dx, real r0, real L0, long ninit
void fractal_do(ARGS);
void fractal_inv(ARGS);
void fractal_trans(ARGS);
void fractal_inv_trans(ARGS);
void fractal_vkcov_free(void);
#endif
