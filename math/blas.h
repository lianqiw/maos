/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifdef __cplusplus
extern "C"{
#endif
/*
   The following are blas and lapack fortran function definitions. Notice that when
   generating MKL custom libraries, use lp64 instead of ilp64 for blas/lapack
   routines to be compatible when usual blas/lapack definitions.

   2013-12-17: replaced int* by ptrdiff_t* to be compatible with blas/lapack that expect 64 bit integers. Still backward compatible with 32bit integers because only one number will be used and we are using little indians.
   2014-03-19: replaced ptrdiff_t by long.
*/
#define BLAS_DEF(X,T)							\
    void X(gemm)(const char* tranA, const char* tranB,			\
		 const ptrdiff_t*nrow_C, const ptrdiff_t*ncol_C, const ptrdiff_t*ncol_opA, \
		 const T*alpha,						\
		 const T*A, const ptrdiff_t*ldA,			\
		 const T*B, const ptrdiff_t*ldB,			\
		 const T*beta, T*C, const ptrdiff_t*ldC);		\
									\
    void X(gemv)(const char* tranA,					\
		 const ptrdiff_t*nrowA, const ptrdiff_t*ncolA, const T*alpha, \
		 const T*A, const ptrdiff_t*ldA,			\
		 const T*X, const ptrdiff_t*incX,			\
		 const T*Beta,						\
		 T*Y, const ptrdiff_t*IncY);				\
									\
    void X(posv)(const char*,const ptrdiff_t*,const ptrdiff_t*,		\
		 T*,const ptrdiff_t*,				\
		 T*,const ptrdiff_t*, ptrdiff_t*);			\
    void X(ger)(ptrdiff_t*,ptrdiff_t*,T*,T*,ptrdiff_t*,T*,ptrdiff_t*,T*,ptrdiff_t*); \
    void X(gesv)(ptrdiff_t*,ptrdiff_t*,T*,ptrdiff_t*,ptrdiff_t*,T*,ptrdiff_t*,ptrdiff_t*); \
    void X(potrf)(const char* , ptrdiff_t* , T *, ptrdiff_t* , ptrdiff_t* );	\

#define BLAS_DEF_REAL(X,T)						\
    void X(gesvd)(const char* jobu, const char* jobvt, ptrdiff_t* m, ptrdiff_t* n, T* a, \
		  ptrdiff_t* lda, T* s, T* u, ptrdiff_t* ldu, T* vt, ptrdiff_t* ldvt, \
		  T* work, ptrdiff_t* lwork, ptrdiff_t* info );		\
    void X(gesdd)(const char *jobz, ptrdiff_t* m, ptrdiff_t* n, T* a, ptrdiff_t* lda, \
		  T* s, T* u, ptrdiff_t* ldu, T* vt, ptrdiff_t* ldvt,	\
		  T* work, ptrdiff_t* lwork, ptrdiff_t* iwork, ptrdiff_t* info);
#define BLAS_DEF_COMP(X,T,R)						\
    void X(gesvd)(const char* jobu, const char* jobvt, ptrdiff_t* m, ptrdiff_t* n, T* a, \
		  ptrdiff_t* lda, R* s, T* u, ptrdiff_t* ldu, T* vt, ptrdiff_t* ldvt, \
		  T* work, ptrdiff_t* lwork, R*rwork, ptrdiff_t* info ); \
    void X(gesdd)(const char *jobz, ptrdiff_t* m, ptrdiff_t* n, T* a, ptrdiff_t* lda, \
		  R* s, T* u, ptrdiff_t* ldu, T* vt, ptrdiff_t* ldvt,	\
		  T* work, ptrdiff_t* lwork, R*rwork, ptrdiff_t* iwork, ptrdiff_t* info); 
#ifndef COMP_SINGLE
#define BLAS_D(A) d##A##_
#define BLAS_C(A) z##A##_
	BLAS_DEF(BLAS_D, real);
	BLAS_DEF(BLAS_C, comp);
	BLAS_DEF_REAL(BLAS_D, real);
	BLAS_DEF_COMP(BLAS_C, comp, real);
#endif

#define BLAS_S(A) s##A##_
#define BLAS_Z(A) c##A##_
	BLAS_DEF(BLAS_S, float);
	BLAS_DEF(BLAS_Z, fcomplex);
	BLAS_DEF_REAL(BLAS_S, float);
	BLAS_DEF_COMP(BLAS_Z, fcomplex, float);
#ifdef __cplusplus
}
#endif

#endif
