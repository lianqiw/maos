/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "random.h"
#include "loc.h"
/**
   \file fractal.h

   Implementation of the fractal operation for atmospheric turbulence screen
   generation and reconstruction.
*/
#define ARGS double *p0, long nx, long ny, double dx, double r0, double L0, long ninit
void fractal(ARGS);
void fractal_inv(ARGS);
void fractal_trans(ARGS);
void fractal_inv_trans(ARGS);
void fractal_vkcov_free(void);
#endif
