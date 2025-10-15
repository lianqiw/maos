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
#ifndef AOS_MATH_NUMTYPE_H
#define AOS_MATH_NUMTYPE_H
#include "../sys/sys.h"

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

#define add_valid(dest, A, B) if(!isnan(B)) dest+=(A)*(B)

/**
   Include headers and finish definition for various complex number operations.
 */

#if defined(__cplusplus)
//When including by CUDA, all definitions are already available.
#if !defined(AOS_CUDA_H) 
//C++ mode
#include <ccomplex>
#include <cmath>
typedef std::complex<double> dcomplex;
typedef std::complex<float> fcomplex;
#ifdef I
#undef I
#endif
#define I dcomplex(0,1)
#define COMPLEX(A,B) dcomplex(A,B)
#define DCOMPLEX(A,B) dcomplex(A,B)
#define FCOMPLEX(A,B) fcomplex(A,B)

#define fabs std::abs
#define cabs std::abs
#define cabsf std::abs

#define cexp std::exp
#define cexpf std::exp
#define cpow std::pow
#define cpowf std::pow
#define csqrt std::sqrt
#define csqrtf std::sqrt
#define clog std::log
#define clogf std::log
#define carg std::arg
#define cargf std::arg

#define cimag std::imag
#define cimagf std::imag
#define creal std::real
#define crealf std::real
using std::conj;//cannot use #define here
#define conjf std::conj

#if !defined(IN_EXTERN_C)
//Prevent compiling failure when number is a literal that defaults to int or double.
static inline fcomplex operator*(double A, const fcomplex& B){
	return B*(float)A;
}
static inline fcomplex operator*(const fcomplex& B, double A){
	return B*(float)A;
}
static inline dcomplex operator*(float A, const dcomplex& B){
	return B*(double)A;
}
static inline dcomplex operator*(const dcomplex& B, float A){
	return B*(double)A;
}
static inline fcomplex operator+(double A, const fcomplex& B){
	return B+(float)A;
}
static inline fcomplex operator+(const fcomplex& B, double A){
	return B+(float)A;
}
static inline dcomplex operator+(float A, const dcomplex& B){
	return B+(double)A;
}
static inline dcomplex operator+(const dcomplex& B, float A){
	return B+(double)A;
}
static inline fcomplex operator-(double A, const fcomplex& B){
	return (float)A-B;
}
static inline fcomplex operator-(const fcomplex& B, double A){
	return B-(float)A;
}
static inline dcomplex operator-(float A, const dcomplex& B){
	return (double)A-B;
}
static inline dcomplex operator-(const dcomplex& B, float A){
	return B-(double)A;
}
#endif
#endif//#ifndef AOS_CUDA_H
#else //#if not defined(__cplusplus) C99 mode
//Do no use tgmath.h. It conflicts with CUDA includes.
//tgmath provides type generic macros for common math routines. It includes math.h and complex.h
#include <math.h>
#include <complex.h>
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
#ifndef I
#define I _Complex_I
#endif
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

#if CPU_SINGLE//run CPU code with single precision.
typedef float real;
typedef fcomplex comp;
#define M_REAL M_FLT
#define M_LOC  M_LOC32
#define M_MAP  M_MAP32
#define M_RMAP M_RMAP32
#define M_COMP M_ZMP
#define EXPI(A) COMPLEX(cosf(A),sinf(A))
#else //run CPU code with double
typedef double real;
typedef dcomplex comp;
#define M_REAL M_DBL
#define M_LOC  M_LOC64
#define M_MAP  M_MAP64
#define M_RMAP M_RMAP64
#define M_COMP M_CMP
#define EXPI(A) COMPLEX(cos(A),sin(A))
#endif
#endif //ifndef AOS_MATH_NUMTYPE_H
