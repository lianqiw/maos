/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_MATH_NUMTYPE_H
#define AOS_MATH_NUMTYPE_H
#include "../sys/sys.h"
#if defined(DLONG)
typedef long spint; /*Only optionally activated in AMD64. */
#else
typedef int spint;  /*This is always 32 bit. */
#endif

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

/**
   Wrap the index for dataset with total of n frames for continuity. The actual data becomes
   0, 1, 2, ..., n-2, n-1, n-2, ..., 0, 1
*/
static inline int wrap(long index, long n){
    long m=n*2-1;
    index=index%m;
    if(index<0) index+=m;
    if(index>=n) index=m-1-index;
    return index;
}

/**
   2010-01-03
   USE_MYFLOOR = 1 reduce timing by 2s
*/
#if USE_MYFLOOR
#define ifloor(A) (A<0?(int)(A)-1:(int)(A))
#else
#define ifloor(A) (int)floor(A)
#endif
#define iceil(A) (int)ceil(A)

#if defined(FP_FAST_FMA)
#define myfma fma
#else
#define myfma(x,y,z) (x)*(y)+(z)
#endif
#define SPLIT(A,B,C) {C=ifloor(A); B=(A)-(C);}

#define add_valid(dest, A, B) if((B)==(B)) dest+=(A)*(B)
#define invalid_val NAN

/**
   Include headers and finish definition for various complex number operations.
 */

#if defined(__cplusplus)
//When including by CUDA, all definitions are already available.
#ifndef AOS_CUDA_GPU_H
//C++ mode
#include <complex>
#include <cmath>
using std::real;
using std::conj;
using std::isinf;
using std::complex;
typedef complex<double> dcomplex;
typedef complex<float> fcomplex;
#define COMPLEX(A,B) dcomplex(A,B)
#define DCOMPLEX(A,B) dcomplex(A,B)
#define FCOMPLEX(A,B) fcomplex(A,B)
#define fabs std::abs
#define cabs fabs
#define cimag imag
#define creal real
#define cexp exp
#define cpow pow
#define csqrt sqrt
#define clog log
#define carg arg
#define cabsf fabs

#define cimagf imag
#define crealf real
#define conjf conj
#define cexpf exp
#define cpowf pow
#define csqrtf sqrt
#define clogf log
#define cargf arg
static inline fcomplex operator*(double A, const fcomplex &B){
    return B*(float)A;
}
static inline fcomplex operator*(const fcomplex &B, double A){
    return B*(float)A;
}
static inline dcomplex operator*(float A, const dcomplex &B){
    return B*(double)A;
}
static inline dcomplex operator*(const dcomplex &B, float A){
    return B*(double)A;
}
static inline fcomplex operator+(double A, const fcomplex &B){
    return B+(float)A;
}
static inline fcomplex operator+(const fcomplex &B, double A){
    return B+(float)A;
}
static inline dcomplex operator+(float A, const dcomplex &B){
    return B+(double)A;
}
static inline dcomplex operator+(const dcomplex &B, float A){
    return B+(double)A;
}
static inline fcomplex operator-(double A, const fcomplex &B){
    return (float)A-B;
}
static inline fcomplex operator-(const fcomplex &B, double A){
    return B-(float)A;
}
static inline dcomplex operator-(float A, const dcomplex &B){
    return (double)A-B;
}
static inline dcomplex operator-(const dcomplex &B, float A){
    return B-(double)A;
}
static inline double conj(double A){
    return A;
}
static inline double real(double A){
    return A;
}

#endif//#ifndef AOS_CUDA_GPU_H
#else //#if defined(__cplusplus) C99 mode
#include <tgmath.h> //never include tgmath.h in CUDA included headers.
typedef __complex__ double dcomplex;
typedef __complex__ float fcomplex;
#define COMPLEX(A,B) ((A)+I*(B))
#define DCOMPLEX(A,B) ((double)(A)+I*(double)(B))
#define FCOMPLEX(A,B) ((float)(A)+I*(float)(B))
#if defined(__CYGWIN__)
/*CYGWIN does not have complex.h. */
#ifndef _Complex_I
#define _Complex_I (__extension__ 1.0iF)
#endif
#define I _Complex_I
double cabs(dcomplex __z);
//double cimag(dcomplex __z);
//double creal(dcomplex __z);
//dcomplex conj(dcomplex __z);
#define cexp(z) exp(creal(z))*(cos(cimag(z))+I*sin(cimag(z)))
dcomplex cpow(dcomplex x, dcomplex z);
dcomplex csqrt(dcomplex);
dcomplex clog(dcomplex);
//double carg(dcomplex);
float cabsf(fcomplex __z);
float cimagf(fcomplex __z);
float crealf(fcomplex __z);
fcomplex conjf(fcomplex __z);
#define cexpf(z) expf(crealf(z))*(cosf(cimagf(z))+I*sinf(cimagf(z)))
fcomplex cpowf(fcomplex x, fcomplex z);
fcomplex csqrtf(fcomplex);
fcomplex clogf(fcomplex);
float cargf(fcomplex);
#else //#if defined(__CYGWIN__)
//C99 already has definitions we need
#endif
#if defined(__FreeBSD__) || defined(__NetBSD__)
/*BSD lacks few routines in C99 implementation*/
static inline dcomplex clog(dcomplex x){
  return log(cabs(x))+I*carg(x);
}
static inline fcomplex clogf(fcomplex x){
  return logf(cabsf(x))+I*cargf(x);
}
static inline dcomplex cpow(dcomplex x, dcomplex z){
  return cexp(clog(x)*z);
}
static inline fcomplex cpowf(fcomplex x, fcomplex z){
  return cexpf(clogf(x)*z);
}
#endif //defined(__FreeBSD__) || defined(__NetBSD__)
#endif //#if defined(__cplusplus) 
#endif //ifndef AOS_MATH_NUMTYPE_H
