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

#ifndef AOS_LAPLACIAN_C
#define AOS_LAPLACIAN_C
#include "loc.h"
void apply_laplacian_map(int nx, int ny, double dx, double r0, double weight, 
			 double *opd, double *opdout);
dsp* mklaplacian_map(int nx, int ny, double dx, double r0, double weight);
dsp* mklaplacian_loc(loc_t *loc, double r0, double weight);
double laplacian_coef(double r0, double weight, double dx);
#endif
