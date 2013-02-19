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
#ifndef AOS_LIB_NR_H
#define AOS_LIB_NR_H
/*
  Contains routines borrowed from numerical recipes
 */
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);
void amoeba(double **p, double y[], int ndim, double ftol, double (*funk)(double [], void *data), void *data, int *nfunk);
typedef double (*minsearch_fun)(double *x, void *info);
int dminsearch(double *x, double *scale, int nmod, double ftol, minsearch_fun fun, void *info);
#endif
