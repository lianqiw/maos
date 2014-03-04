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

#ifndef AOS_LIB_MATHMISC_H
#define AOS_LIB_MATHMISC_H
#include "../sys/sys.h"
#include "blas.h"
/**
   \file mathmisc.h
   A few math routines
*/
double factorial(long n);
void normalize_sum(double *p, long nloc, double norm);
void normalize_max(double *p, long nloc, double max);
double dotdbl(const double *restrict p1, const double *restrict p2, 
	      const double *restrict p3, long n);
float dotflt(const float *restrict p1, const float *restrict p2, 
	     const float *restrict p3, long n);
dcomplex dotcmp(const dcomplex *restrict p1, const dcomplex *restrict p2, 
		const double *restrict p3, long n);
fcomplex dotzmp(const fcomplex *restrict p1, const fcomplex *restrict p2, 
		const float *restrict p3, long n);
double maxdbl(const double *p, long n);
double maxabs(const double *p, long n);
float maxabsf(const float *p, long n);

double dblsum(double *p, long nx);
float fltsum(float *p, long nx);
#define sinc(x) (fabs(x)<1.e-5?1:sin(M_PI*x)/(M_PI*x))

void adddbl(double *restrict out, double alpha, 
	    const double *in, int N, double beta, double gamma);
void invsq(long n, double *restrict A);
#define mysqrt(A) (A<0?-sqrt(-A):sqrt(A))
long *invperm(long *p, long np);
void cvecperm(dcomplex *restrict out, const dcomplex *restrict in, const long *perm, long nx);
void cvecpermi(dcomplex *restrict out, const dcomplex *restrict in, const long *perm, long nx);
void maxminlong(const long *restrict p, long N,
		long *restrict max, long *restrict min);
void remove_piston(double *p, long n);
long nextpow2(long n);
long nextfftsize(long n);
unsigned long mylog2(unsigned long n);
typedef double(*golden_section_fun)(void *param, double x);
double golden_section_search(golden_section_fun f, void *param, 
			     double x1, double x4, double tau);
#endif
