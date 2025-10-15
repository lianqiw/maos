/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_MATH_DEFS_H
#define AOS_MATH_DEFS_H
#include "numtype.h"
#undef COMPLEX
#undef IMAG
#undef REAL
#undef MAT_TYPE
#undef XR
#undef XC
#undef R
#undef RD
#undef FFTW
#undef X
#undef Y
#undef Z
#undef T
#undef TD
#undef M_T
#undef M_SPT64
#undef M_SPT32
#undef RANDU
#undef RANDN
#undef PRINT
#undef PMAT
#undef PCELL
#undef PSPCELL
#undef M_SPT
#undef EXPI
#ifdef COMP_LONG
#define X(A) l##A
#else

/*
  when CPU_SINGLE is enabled, and COMP_SINGLE is not set, the function prefix
  are d, but the actual data is float.
 */
#ifndef COMP_SINGLE
#ifndef COMP_COMPLEX 
#define X(A) d##A
#else
#define X(A) c##A
#endif
#define XR(A) d##A
#define XC(A) c##A
#else 
#ifndef COMP_COMPLEX 
#define X(A) s##A
#else
#define X(A) z##A
#endif
#define XR(A) s##A
#define XC(A) z##A
#endif
#endif
//Long Integers
#ifdef COMP_LONG
#define T long
#define TD long
#define R long
#define RD long
#define ABS labs
#define FABS labs
#define REAL(A) A
#define REALC(A) A
#define IMAG(A) 0
#define CONJ(A) A
#define ARG(A) 0
#define M_T M_LONG
#define PRINT(A) info(" %ld",A);
#else //if COMP_LONG

#define MAT_TYPE

/*Double precision*/
#if !defined(COMP_SINGLE) && CPU_SINGLE == 0
#define R double
#define RD double
#define RC dcomplex
#define FFTW(A) fftw_##A

#ifndef COMP_COMPLEX //Double Real
#define Y(A) D##A
#define Z(A) d##A##_
#define T double
#define TD double
#define M_T M_DBL
#define M_SPT64 M_DSP64
#define M_SPT32 M_DSP32
#define REAL
#define REALC creal
#define IMAG(A) 0
#define CONJ
#define CARG(A) 0
#define ABS fabs
#define FABS fabs
#define RANDU(A) randu(A)
#define RANDN(A) randn(A)
#define PRINT(A) info(" %10.3g",A);
#else //Double Complex
#define Y(A) C##A
#define Z(A) z##A##_ /*blas/lapack convention */
#define T dcomplex
#define TD dcomplex
#define M_T M_CMP
#define M_SPT64 M_CSP64
#define M_SPT32 M_CSP32
#define COMPLEX DCOMPLEX
#define REAL creal
#define REALC creal
#define IMAG cimag
#define CONJ conj
#define CARG carg
#define ABS cabs
#define FABS fabs
#define RANDU(A) COMPLEX(randu(A),randu(A))
#define RANDN(A) COMPLEX(randn(A),randn(A))
#define PRINT(A) info("(%10.3e %10.3eI)",REAL(A),IMAG(A));
#define EXPI(A) COMPLEX(cos(A),sin(A))
#endif
#else 
#define COMP_SINGLE //Single Precision
#define R float
#define RD double
#define RC fcomplex
#define FFTW(A) fftwf_##A
#ifndef COMP_COMPLEX ///Float Real
#define Y(A) S##A
#define Z(A) s##A##_
#define T float
#define TD double
#define M_T M_FLT
#define M_SPT64 M_SSP64
#define M_SPT32 M_SSP32
#define REAL
#define REALC crealf
#define IMAG(A) 0
#define CONJ
#define CARG(A) 0
#define ABS fabsf
#define FABS fabsf
#define RANDU(A) (float)randu(A)
#define RANDN(A) (float)randn(A)
#define PRINT(A) info("%10.3e",A);
#else
#define Y(A) Z##A //Float Complex
#define Z(A) c##A##_ /*blas/lapack convention */
#define T fcomplex
#define TD dcomplex
#define M_T M_ZMP
#define M_SPT64 M_ZSP64
#define M_SPT32 M_ZSP32
#define COMPLEX FCOMPLEX
#define REAL crealf
#define REALC crealf
#define IMAG cimagf
#define CONJ conjf
#define CARG cargf
#define ABS cabsf
#define FABS fabsf
#define RANDU(A) COMPLEX((float)randu(A),(float)randu(A))
#define RANDN(A) COMPLEX((float)randn(A),(float)randn(A))
#define PRINT(A) info("(%10.3e %10.3eI)",REAL(A),IMAG(A));
#define EXPI(A) COMPLEX(cosf(A),sinf(A))
#endif/*#define COMP_COMPLEX */
#endif
/*
#ifndef COMP_COMPLEX
#undef conj
#define conj
#undef real
#define real
#define creal2 creal
#undef creal
#define creal
#endif
*/
#if defined(DLONG)
#define M_SPT M_SPT64
#else
#define M_SPT M_SPT32
#endif
static inline int issp(const void* id){
	return id?(*((const M_ID*)id)==M_SPT):0;
}

#endif //if COMP_LONG
#define DO_PRAGMA(A...) _Pragma(#A)
#ifdef COMP_COMPLEX 
#define OMP_SIMD(A...) //DO_PRAGMA(omp simd A) //simd is not available for complex numbers
#define OMP_SIMD_R(A...) //simd is not available for complex numbers for reduction
#else
#define OMP_SIMD(A...) DO_PRAGMA(omp simd A)
#define OMP_SIMD_R(A...) DO_PRAGMA(omp simd A)
#endif
#define ABS2(A) REAL((A)*CONJ(A))

static inline int ismat(const void* id){
	return id?(*(const uint16_t*)id==M_T):0;
}
//Check that A is valid and has mat type and has non zero size.
#define check_mat1(A) ((A)?(check(ismat(A))?1:(error("type id mismatch, got %x, expect %x\n", *(M_ID*)A, M_T),0)):0)
#define check_mat2(A,B) (check_mat1(A) && check_mat1(B))
#define check_mat3(A,B,C) (check_mat1(A) && check_mat1(B) && check_mat1(C))
#define check_mat(...) P_GET3(_0,__VA_ARGS__,check_mat3,check_mat2,check_mat1)(__VA_ARGS__)
#define check_match(A,B) ((check_mat(A,B) && check((A)->nx==(B)->nx) && check((A)->ny==(B)->ny))?1:0)
#define check_match_vec(A,B) ((check_mat(A,B) && check((A)->nx*(A)->ny==(B)->nx*(B)->ny))?1:0)
#endif //ifndef AOS_MATH_DEFS_H
