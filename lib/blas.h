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

#ifndef AOS_LIB_BLAS_H
#define AOS_LIB_BLAS_H

#include "../sys/sys.h"
#include <stddef.h>
/**
   \file blas.h
   The following are blas and lapack fortran function definitions. Notice that when
   generating MKL custom libraries, use lp64 instead of ilp64 for blas/lapack
   routines to be compatible when usual blas/lapack definitions. 

   2013-12-17: replaced int* by ptrdiff_t* to be compatible with blas/lapack that expect 64 bit integers. Still backward compatible with 32bit integers because only one number will be used and we are using little indians.
*/
void dgemm_(const char* tranA, const char* tranB,
	    const ptrdiff_t*nrow_C, const ptrdiff_t*ncol_C, const ptrdiff_t*ncol_opA, const double*alpha, 
	    const double*A, const ptrdiff_t*ldA, 
	    const double*B, const ptrdiff_t*ldB,    
	    const double*beta, double*C, const ptrdiff_t*ldC);
void sgemm_(const char*, const char*,   
	    const ptrdiff_t*, const ptrdiff_t*, const ptrdiff_t*, const float*, 
	    const float*, const ptrdiff_t*, 
	    const float*, const ptrdiff_t*,    
	    const float*, float*, const ptrdiff_t*);
void zgemm_(const char*, const char*,   
	    const ptrdiff_t*, const ptrdiff_t*, const ptrdiff_t*, const dcomplex*, 
	    const dcomplex*, const ptrdiff_t*, 
	    const dcomplex*, const ptrdiff_t*,    
	    const dcomplex*, dcomplex*, const ptrdiff_t*);

void dgemv_(const char* tranA, 
	    const ptrdiff_t*nrowA, const ptrdiff_t*ncolA, const double*alpha, 
	    const double*A, const ptrdiff_t*ldA, 
	    const double*X, const ptrdiff_t*incX, 
	    const double*Beta,
	    double*Y, const ptrdiff_t*IncY);
void sgemv_(const char* tranA, 
	    const ptrdiff_t*nrowA, const ptrdiff_t*ncolA, const float*alpha, 
	    const float*A, const ptrdiff_t*ldA, 
	    const float*X, const ptrdiff_t*incX, 
	    const float*Beta,
	    float*Y, const ptrdiff_t*IncY);
void zgemv_(const char* tranA, 
	    const ptrdiff_t*nrowA, const ptrdiff_t*ncolA, const dcomplex*alpha, 
	    const dcomplex*A, const ptrdiff_t*ldA, 
	    const dcomplex*X, const ptrdiff_t*incX, 
	    const dcomplex*Beta,
	    dcomplex*Y, const ptrdiff_t*IncY);
/*Lapack */
void dposv_(const char*,const ptrdiff_t*,const ptrdiff_t*, double*,const ptrdiff_t*,
	    double*,const ptrdiff_t*, ptrdiff_t*);
void sposv_(const char*,const ptrdiff_t*,const ptrdiff_t*, float*,const ptrdiff_t*,
	    float*,const ptrdiff_t*, ptrdiff_t*);
void zposv_(const char*,const ptrdiff_t*,const ptrdiff_t*, dcomplex*,const ptrdiff_t*,
	    dcomplex*,const ptrdiff_t*, ptrdiff_t*);
void dger_ (ptrdiff_t*,ptrdiff_t*,double*,double*,ptrdiff_t*,double*,ptrdiff_t*,double*,ptrdiff_t*);
void dgesv_(ptrdiff_t*,ptrdiff_t*,double*,ptrdiff_t*,ptrdiff_t*,double*,ptrdiff_t*,ptrdiff_t*);
void sger_ (ptrdiff_t*,ptrdiff_t*,float*,float*,ptrdiff_t*,float*,ptrdiff_t*,float*,ptrdiff_t*);
void sgesv_(ptrdiff_t*,ptrdiff_t*,float*,ptrdiff_t*,ptrdiff_t*,float*,ptrdiff_t*,ptrdiff_t*);
void zgesv_(ptrdiff_t*,ptrdiff_t*,dcomplex*,ptrdiff_t*,ptrdiff_t*,dcomplex*,ptrdiff_t*,ptrdiff_t*);
void dpotrf_(char* , ptrdiff_t* , double *, ptrdiff_t* , ptrdiff_t* );
void spotrf_(char* , ptrdiff_t* , float *, ptrdiff_t* , ptrdiff_t* );
void zpotrf_(char* , ptrdiff_t* , dcomplex *, ptrdiff_t* , ptrdiff_t* );
/* DGESVD prototype */
extern void 
dgesvd_(char* jobu, char* jobvt, ptrdiff_t* m, ptrdiff_t* n, double* a,
	ptrdiff_t* lda, double* s, double* u, ptrdiff_t* ldu, double* vt, ptrdiff_t* ldvt,
	double* work, ptrdiff_t* lwork, ptrdiff_t* info );
extern void 
sgesvd_(char* jobu, char* jobvt, ptrdiff_t* m, ptrdiff_t* n, float* a,
	ptrdiff_t* lda, float* s, float* u, ptrdiff_t* ldu, float* vt, ptrdiff_t* ldvt,
	float* work, ptrdiff_t* lwork, ptrdiff_t* info );
extern void 
zgesvd_(char* jobu, char* jobvt, ptrdiff_t* m, ptrdiff_t* n, dcomplex* a,
	ptrdiff_t* lda, double* s, dcomplex* u, ptrdiff_t* ldu, dcomplex* vt, ptrdiff_t* ldvt,
	dcomplex* work, ptrdiff_t* lwork, double *rwork, ptrdiff_t* info );

void dsyev_(char *jobz, char *uplo, ptrdiff_t* n, double *a,
	    ptrdiff_t* lda, double *w, double *work, ptrdiff_t* lwork,
	    ptrdiff_t* info);
void ssyev_(char *jobz, char *uplo, ptrdiff_t* n, float *a,
	    ptrdiff_t* lda, float *w, float *work, ptrdiff_t* lwork,
	    ptrdiff_t* info);
void zheev_(char *jobz, char *uplo, ptrdiff_t* n, dcomplex *a,
	    ptrdiff_t* lda, double *w, dcomplex *work, ptrdiff_t* lwork, double *rwork,
	    ptrdiff_t* info);

#if USE_MKL==1 && !defined(_OPENMP)
void omp_set_num_threads(int n);
#endif
#endif
