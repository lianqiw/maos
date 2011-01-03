/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
dmat *vonkarman_psd(int nx, int ny, double dx, double r0, double L0, double power0);
dmat *biharmonic_psd(int nx, int ny, double dx, double r0, double L0, double power0);

void map_shm(map_t **screen, long totmem, int nlayer, int fd, int rw);
map_t **atmnew_shm(int *fd, int *inshm, rand_t *rstat, long nx, long ny, double dx, 
		   double r0, double L0, double *wt, int nlayer, int method);
map_t** genscreen_from_spect(rand_t *rstat, dmat *spect, double dx,double r0, double L0,
			     double* wt, int nlayer, int nthread);
map_t** vonkarman_screen(rand_t *rstat, int m, int n, double dx, 
			 double r0, double L0, double* wt, int nlayer, int nthread);
map_t** biharmonic_screen(rand_t *rstat, int m, int n, double dx, 
			  double r0, double L0, double* wt, int nlayer, int nthread);
map_t **fractal_screen(rand_t *rstat, long nx, long ny, double dx, double r0,
		       double L0, double* wt, int nlayer, int nthread);
#endif
