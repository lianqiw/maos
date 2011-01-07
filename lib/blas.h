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

#ifndef AOS_BLAS_H
#define AOS_BLAS_H

#include "common.h"
/*
The following are blas and lapack fortran function definitions. Notice that when
generating MKL custom libraries, use lp64 instead of ilp64 for blas/lapack
routines to be compatible when usual blas/lapack definitions. */
void dgemm_(const char*, const char*,   
	    const int*, const int*, const int*, const double*, 
	    const double*, const int*, 
	    const double*, const int*,    
	    const double*, double*, const int*);
void zgemm_(const char*, const char*,   
	    const int*, const int*, const int*, const dcomplex*, 
	    const dcomplex*, const int*, 
	    const dcomplex*, const int*,    
	    const dcomplex*, dcomplex*, const int*);

//Lapack
void dposv_(const char*,const int*,const int*, double*,const int*,
	    double*,const int*, int*);
void zposv_(const char*,const int*,const int*, dcomplex*,const int*,
	    dcomplex*,const int*, int*);
void dger_ (int*,int*,double*,double*,int*,double*,int*,double*,int*);
void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
void zgesv_(int*,int*,dcomplex*,int*,int*,dcomplex*,int*,int*);

/* DGESVD prototype */
extern void 
dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a,
	int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
	double* work, int* lwork, int* info );
extern void 
zgesvd_(char* jobu, char* jobvt, int* m, int* n, dcomplex* a,
	int* lda, double* s, dcomplex* u, int* ldu, dcomplex* vt, int* ldvt,
	dcomplex* work, int* lwork, double *rwork, int* info );
#if USE_MKL==1
void omp_set_num_threads(int *n);
#endif
#endif
