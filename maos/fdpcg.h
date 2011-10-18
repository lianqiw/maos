/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_RECON_FDPCG
#define AOS_RECON_FDPCG
#include "maos.h"
csp* fdpcg_saselect(long nx, long ny, double dx, loc_t *saloc, double *saa);
long *fdpcg_perm(const long *nx, const long *ny, long pos, int nps);
void fdpcg_g(cmat **gx, cmat **gy, long nx, long ny, double dx, double dsa);
csp *fdpcg_prop(long nps, long pos, const int *os, long nxg, double dx, 
		double *dispx, double *dispy);
FDPCG_T *fdpcg_prepare(const PARMS_T *parms, const RECON_T *recon, 
	       const POWFS_T *powfs);
void fdpcg_precond(dcell **xout, const void *A, const dcell *xin);
void fdpcg_free(FDPCG_T *fdpcg);
#endif
