/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
/*
  To be included in mat.c, cell.c and matbin.c
*/
#if MAT_VERBOSE == 1
#define matinfo(A...) {fprintf(stdout, A);}
#else
#define matinfo(A...)
#endif
#ifndef MAT_TYPE
#define MAT_TYPE
#ifndef USE_SINGLE
#define XR(A) d##A
#ifndef USE_COMPLEX
/*Double */
#define X(A) d##A
#define Y(A) A
#define Z(A) d##A##_
#define T double
#define R double
#define M_T M_DBL
#define M_TT M_DMAT
#define M_SPT64 M_SP64
#define M_SPT32 M_SP32
#define REAL(A) (A)
#define ABS(A) fabs(A)
#define SQRT(A) sqrt(A)
#define RANDU(A) randu(A)
#define RANDN(A) randn(A)
#define PRINT(A) printf("%10.3e",A);
#define CONJ(x) (x)
#define DOT dotdbl
#define POW pow
#define LOG log
#else
/*Double Complex */
#define X(A) c##A
#define Y(A) c##A
#define Z(A) z##A##_ /*blas/lapack convention */
#define T dcomplex
#define R double
#define M_T M_CMP
#define M_TT M_CMAT
#define M_SPT64 M_CSP64
#define M_SPT32 M_CSP32
#define REAL(A) creal(A)
#define ABS(A) cabs(A)
#define SQRT(A) csqrt(A)
#define RANDU(A) (randu(A)+I*randu(A))
#define RANDN(A) (randn(A)+I*randn(A))
#define PRINT(A) printf("(%10.3e %10.3eI)",creal(A),cimag(A));
#define CONJ(x) conj(x)
#define DOT dotcmp
#define POW cpow
#define LOG clog
#endif
#else /*#define USE_SINGLE */
#define XR(A) s##A
/*Float */
#ifndef USE_COMPLEX
#define X(A) s##A
#define Y(A) s##A
#define Z(A) s##A##_
#define T float
#define R float
#define M_T M_FLT
#define M_TT M_SMAT
#define M_SPT64 M_SSP64
#define M_SPT32 M_SSP32
#define REAL(A) (A)
#define ABS(A) fabsf(A)
#define SQRT(A) sqrtf(A)
#define RANDU(A) (float)randu(A)
#define RANDN(A) (float)randn(A)
#define PRINT(A) printf("%10.3e",A);
#define CONJ(x) (x)
#define DOT dotflt
#define POW powf
#define LOG logf
#else
/*Single Complex */
#define X(A) z##A
#define Y(A) z##A
#define Z(A) c##A##_ /*blas/lapack convention */
#define T fcomplex
#define R float
#define M_T M_ZMP
#define M_TT M_ZMAT
#define M_SPT64 M_ZSP64
#define M_SPT32 M_ZSP32
#define REAL(A) crealf(A)
#define ABS(A) cabsf(A)
#define SQRT(A) csqrtf(A)
#define RANDU(A) ((float)randu(A)+I*(float)randu(A))
#define RANDN(A) ((float)randn(A)+I*(float)randn(A))
#define PRINT(A) printf("(%10.3e %10.3eI)",crealf(A),cimagf(A));
#define CONJ(x) conjf(x)
#define DOT dotzmp
#define POW cpowf
#define LOG clogf
#endif/*#define USE_COMPLEX */
#endif/*#define USE_SINGLE */
#endif/*#define MATTYPE */
#define PMAT(A,pp) T (*restrict pp)[(A)->nx]=(void *)(A)->p
#define PCELL(M,P) X(mat)* (*restrict P)[(M)->nx]=(void*)(M)->p
#define PSPCELL(M,P) X(sp)* (*restrict P)[(M)->nx]=(void *)(M)->p

#ifdef DLONG
#define M_SPT M_SPT64
#else
#define M_SPT M_SPT32
#endif
