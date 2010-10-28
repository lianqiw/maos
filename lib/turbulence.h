/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_TURBULENCE_H
#define AOS_TURBULENCE_H

#include "loc.h"
#include "random.h"
#include "type.h"
extern int disable_atm_shm;
extern int genscreen_keep_unused;
dmat* vonkarman_spect(int nx, int ny, double dx, double r0, double L0);
dmat *vonkarman_invpsd(int nx, int ny, double dx, double r0, double L0);
MAP_T** genscreen_from_spect(struct_rand *rstat, dmat *spect, double dx,
			       double* wt, int nlayer, int nthread);
MAP_T** vonkarman_screen(struct_rand *rstat, int m, int n, double dx, 
			   double r0, double L0, double* wt, int nlayer, int nthread);
MAP_T** biharmonic_screen(struct_rand *rstat, int m, int n, double dx, 
			    double r0, double L0, double* wt, int nlayer, int nthread);
#endif
