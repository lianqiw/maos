/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#define matinfo(A...) {fprintf(stderr, A);}
#else
#define matinfo(A...)
#endif
#ifndef AOS_MATH_DEFS_H
#define AOS_MATH_DEFS_H
#undef COMPLEX
#undef MAT_TYPE
#undef XR
#undef XC
#undef R
#undef FFTW
#undef X
#undef Y
#undef Z
#undef T
#undef M_T
#undef M_SPT64
#undef M_SPT32
#undef RANDU
#undef RANDN
#undef PRINT
#undef DOT
#undef PMAT
#undef PCELL
#undef PSPCELL
#undef M_SPT
#ifdef USE_LONG
#define T long
#define X(A) l##A
#undef fabs
#define fabs(A) abs(A)
#define M_T M_LONG
#define PRINT(A) fprintf(stderr," %ld",A);
#else //if USE_LONG

#define MAT_TYPE
/*Double precision*/
#ifndef USE_SINGLE
#define XR(A) d##A
#define XC(A) c##A
#define R double
#define RI dcomplex
#define FFTW(A) fftw_##A
#ifndef USE_COMPLEX
/*Double */
#define X(A) d##A
#define Y(A) D##A
#define Z(A) d##A##_
#define T double
#define M_T M_DBL
#define M_SPT64 M_DSP64
#define M_SPT32 M_DSP32
#define RANDU(A) randu(A)
#define RANDN(A) randn(A)
#define PRINT(A) fprintf(stderr," %10.3e",A);
#define DOT dotdbl
#else/*Double Complex */
#define X(A) c##A
#define Y(A) C##A
#define Z(A) z##A##_ /*blas/lapack convention */
#define T dcomplex
#define M_T M_CMP
#define M_SPT64 M_CSP64
#define M_SPT32 M_CSP32
#define COMPLEX DCOMPLEX
#define RANDU(A) COMPLEX(randu(A),randu(A))
#define RANDN(A) COMPLEX(randn(A),randn(A))
#define PRINT(A) fprintf(stderr,"(%10.3e %10.3eI)",creal(A),cimag(A));
#define DOT dotcmp
#define EXPI(A) COMPLEX(cos(A),sin(A))
#endif
#else 
/*Single Precision*/
#define XR(A) s##A
#define XC(A) z##A
#define R float
#define RI fcomplex
#define FFTW(A) fftwf_##A
/*Float */
#ifndef USE_COMPLEX
#define X(A) s##A
#define Y(A) S##A
#define Z(A) s##A##_
#define T float
#define M_T M_FLT
#define M_SPT64 M_SSP64
#define M_SPT32 M_SSP32
#define RANDU(A) (float)randu(A)
#define RANDN(A) (float)randn(A)
#define PRINT(A) fprintf(stderr,"%10.3e",A);
#define DOT dotflt
#else
/*Single Complex */
#define X(A) z##A
#define Y(A) Z##A
#define Z(A) c##A##_ /*blas/lapack convention */
#define T fcomplex
#define M_T M_ZMP
#define M_SPT64 M_ZSP64
#define M_SPT32 M_ZSP32
#define COMPLEX FCOMPLEX
#define RANDU(A) COMPLEX((float)randu(A),(float)randu(A))
#define RANDN(A) COMPLEX((float)randn(A),(float)randn(A))
#define PRINT(A) fprintf(stderr,"(%10.3e %10.3eI)",crealf(A),cimagf(A));
#define DOT dotzmp
#define EXPI(A) COMPLEX(cosf(A),sinf(A))
#endif/*#define USE_COMPLEX */
#endif/*#define USE_SINGLE */

#define PSPCELL(A,pp) PALL(X(sp)*,A,pp)

#ifdef DLONG
#define M_SPT M_SPT64
#else
#define M_SPT M_SPT32
#endif
INLINE int issp(const void *id){
    const uint32_t magic=*((const uint32_t*)id);
    return (magic==M_SPT);
}

#endif //if USE_LONG

#define ABS2(A) ((A)*conj(A))
#define PMAT(A,pp) PALL(T,A,pp)
#define PCELL(A,pp) PALL(X(mat)*,A,pp)
INLINE int ismat(const void *id){
    const uint32_t magic=*((const uint32_t*)id);
    return (magic==M_T);
}

#endif //ifndef AOS_MATH_DEFS_H
