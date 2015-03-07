/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_COMMON_H
#define AOS_COMMON_H
typedef void (*quitfun_t)(const char*);
extern quitfun_t quitfun;
void default_quitfun(const char *msg);
/**
   \file common.h
   Every source file in this folder should include this file
*/
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
#if defined(__cplusplus) && !defined(AOS_CUDA_GPU_H)
#include <csignal>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#define isfinite std::isfinite
using std::signbit;
using std::strerror;
#define isnan std::isnan
#define isinf std::isinf
#else
#include <signal.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#endif
#define INLINE inline __attribute__((always_inline))
#include <unistd.h>
#if !defined(__FreeBSD__) && !defined(__NetBSD__)
#include <alloca.h>
#endif
#include <time.h>
enum{
    T_INT=1,
    T_DBL=2,
    T_STR=3,
    T_LONG=4,
    T_INTARR=11,
    T_DBLARR=12,
    T_LONGARR=14,
};
#if defined(DLONG)
typedef unsigned long spint; /*Only optionally activated in AMD64. */
#define M_SPINT M_INT64
#else
typedef unsigned int spint;  /*This is always 32 bit. */
#define M_SPINT M_INT32
#endif
#undef	MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#undef	MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#if defined(__CYGWIN__)
/*CYGWIN does not have complex.h. */
#define complex __complex__
#define _Complex_I (__extension__ 1.0iF)
#define I _Complex_I
double cabs(dcomplex __z);
double cimag(dcomplex __z);
double creal(dcomplex __z);
dcomplex conj(dcomplex __z);
#define cexp(z) exp(creal(z))*(cos(cimag(z))+I*sin(cimag(z)))
dcomplex cpow(dcomplex x, dcomplex z);
dcomplex csqrt(dcomplex);
dcomplex clog(dcomplex);
double carg(dcomplex);
float cabsf(fcomplex __z);
float cimagf(fcomplex __z);
float crealf(fcomplex __z);
fcomplex conjf(fcomplex __z);
#define cexpf(z) expf(crealf(z))*(cosf(cimagf(z))+I*sinf(cimagf(z)))
fcomplex cpowf(fcomplex x, fcomplex z);
fcomplex csqrtf(fcomplex);
fcomplex clogf(fcomplex);
float cargf(fcomplex);
#elif defined(__FreeBSD__) || defined(__NetBSD__)
#include <complex.h>
/*BSD lacks cpow in C99 implementation*/
INLINE dcomplex clog(dcomplex x){
  return log(cabs(x))+I*carg(x);
}
INLINE fcomplex clogf(fcomplex x){
  return logf(cabsf(x))+I*cargf(x);
}
INLINE dcomplex cpow(dcomplex x, dcomplex z){
  return cexp(clog(x)*z);
}
INLINE fcomplex cpowf(fcomplex x, fcomplex z){
  return cexpf(clogf(x)*z);
}
#else
//C99 always has definitions we need
//if not included by cuda
#if !defined(__cplusplus) || defined(AOS_CUDA_GPU_H)
#include <complex.h>
#else
#include <complex>
using std::complex;
typedef complex<double> dcomplex;
typedef complex<float> fcomplex;
#define COMPLEX(A,B) dcomplex(A,B)
#define DCOMPLEX(A,B) dcomplex(A,B)
#define FCOMPLEX(A,B) fcomplex(A,B)

#define cabs abs
#define cimag imag
#define creal real
#define cexp exp
#define cpow pow
#define csqrt sqrt
#define clog log
#define carg arg
#define cabsf abs

#define cimagf imag
#define crealf real
#define conjf conj
#define cexpf exp
#define cpowf pow
#define csqrtf sqrt
#define clogf log
#define cargf arg
inline fcomplex operator*(double A, const fcomplex &B){
    return B*(float)A;
}
inline fcomplex operator*(const fcomplex &B, double A){
    return B*(float)A;
}
inline dcomplex operator*(float A, const dcomplex &B){
    return B*(double)A;
}
inline dcomplex operator*(const dcomplex &B, float A){
    return B*(double)A;
}
inline fcomplex operator+(double A, const fcomplex &B){
    return B+(float)A;
}
inline fcomplex operator+(const fcomplex &B, double A){
    return B+(float)A;
}
inline dcomplex operator+(float A, const dcomplex &B){
    return B+(double)A;
}
inline dcomplex operator+(const dcomplex &B, float A){
    return B+(double)A;
}
inline fcomplex operator-(double A, const fcomplex &B){
    return (float)A-B;
}
inline fcomplex operator-(const fcomplex &B, double A){
    return B-(float)A;
}
inline dcomplex operator-(float A, const dcomplex &B){
    return (double)A-B;
}
inline dcomplex operator-(const dcomplex &B, float A){
    return B-(double)A;
}
#endif
#endif

#ifndef COMPLEX
typedef double __complex__ dcomplex;
typedef float  __complex__ fcomplex;
#define COMPLEX(A,B) ((A)+(B)*I)
#define DCOMPLEX(A,B) ((double)(A)+(double)(B)*I)
#define FCOMPLEX(A,B) ((float)(A)+(float)(B)*I)
#endif

#include "mem.h"
#ifdef __linux__
#include <linux/limits.h> /*includes definition of PATH_MAX */
#else
#include <limits.h>
#endif/*__linux__ */
#ifndef restrict
#define restrict __restrict
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define SEC2RAD 4.848136811095360e-06 //arcsec in unit of radian
#define RAD2SEC 206264.8062470964 //radian in unit of arcsec
#ifndef EPS
#define EPS 1.e-15
#endif

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
#define not_nan(A) ((A)==(A))
#define is_nan(A) (!((A)==(A)))

#define add_valid(dest, A, B) if((B)==(B)) dest+=(A)*(B)
#define invalid_val NAN

#define BASEFILE (strrchr(__FILE__, '/') ?strrchr(__FILE__, '/')+1  : __FILE__)
long thread_id(void);
void print_backtrace();
/*
  use () to make the statements a single statement.
*/
#ifndef error
#define QUIT_FUN(A) quitfun?quitfun(A):default_quitfun(A);
#define info(A...) ({char fline[4096];char sect[4096];			\
	    snprintf(sect, 4096,"%s:%d",BASEFILE,__LINE__); snprintf(fline,4096, "%-20s",sect); \
	    snprintf(sect, 4096, A);strncat(fline,sect,4096-strlen(fline)-1); \
	    fprintf(stderr,"%s", fline);})

#define error(A...) ({char fline[4096];char sect[4096];			\
	    snprintf(sect, 4096,"%s:%d",BASEFILE,__LINE__);		\
	    snprintf(fline,4096, "\033[01;31m%-20s Fatal error(%ld): ",sect,thread_id()); \
	    snprintf(sect, 4096, A);strncat(fline,sect,4096-strlen(fline)-1); \
	    fprintf(stderr,"%s\033[00;00m", fline); QUIT_FUN("error");})

#define warning(A...) ({char fline[4096];char sect[4096];		\
	    snprintf(sect, 4096,"%s:%d",BASEFILE,__LINE__);		\
	    snprintf(fline,4096, "\033[01;31m%-20s",sect); \
	    snprintf(sect, 4096, A);strncat(fline,sect,4096-strlen(fline)-1); \
	    fprintf(stderr,"%s\033[00;00m", fline); })


#define info2(A...) fprintf(stderr, A)
#define error2(A...) ({ fprintf(stderr, "\033[01;31mFatal error\033[00;00m\t" A); QUIT_FUN("ERROR");})
#define warning2(A...) ({fprintf(stderr,"\033[00;31mWarning:\033[00;00m" A);})

#define info3(A...) ({char fline[4096];char sect[4096];			\
		      snprintf(sect, 4096,"%s:%d",BASEFILE,__LINE__);	\
		      snprintf(fline,4096, "[%s]%-20s",sect, myasctime()); \
		      snprintf(sect, 4096, A);strncat(fline,sect,4096-strlen(fline)-1); \
		      fprintf(stderr,"%s", fline);})
#define error3(A...) ({char fline[4096];char sect[4096];		\
	    snprintf(sect, 4096,"%s:%d",BASEFILE,__LINE__);		\
	    snprintf(fline,4096, "[%s]\033[01;31m%-20s Fatal error:(%ld) ",myasctime(),sect,thread_id()); \
	    snprintf(sect, 4096, A);strncat(fline,sect,4096-strlen(fline)-1); \
	    fprintf(stderr,"%s\033[00;00m", fline); QUIT_FUN("ERROR");})

#define warning3(A...) ({char fline[4096];char sect[4096];		\
	    snprintf(sect, 4096,"%s:%d",BASEFILE,__LINE__);		\
	    snprintf(fline,4096, "[%s]\033[01;31m%-20s ",myasctime(),sect); \
	    snprintf(sect, 4096, A);strncat(fline,sect,4096-strlen(sect)-1); \
	    fprintf(stderr,"%s\033[00;00m", fline);})

#define warning_once(A...) ({static int done=0; if(!done){done=1; warning(A);}})
#define info_once(A...) ({static int done=0; if(!done){done=1; info2(A);}})

#endif
#ifndef assert
#if DEBUG
#define assert(A) if(!(A)) error("assertion failed: %s\n", #A)
#else
#define assert(A)
#endif
#endif
#define error_write error("Write failed\n")
#define error_read error("Read failed\n")
/**
   Functions that return realtime:
   time(): resolution is in integer second. not enough resolution.
   gettimeofday(): must use struct as input argument.
   Functions that return tics:
   times(): the argument can not be NULL in Mac. obsolete by gettimeofday in Mac
   clock(): returns processor time, not real time.

   Decision: use customize function myclockd() in utils.c which returns a double for second. 
*/
#define TIC double tk
#define tic tk=myclockd();
#define toc(A...) ({char fline[4096];char sect[4096];			\
	    snprintf(sect, 4096,"%s:%d",BASEFILE,__LINE__); snprintf(fline,4096, "%-20s",sect); \
	    snprintf(sect, 4096, A);strncat(fline,sect,4096-strlen(fline)-1); \
	    fprintf(stderr,"%s takes %.6f seconds.\n", fline, myclockd()-tk);})
#define toc2(A...) ({char fline[4096]; snprintf(fline, 4096, A); fprintf(stderr, "%s takes %.6f seconds.\n",fline, myclockd()-tk);})
#define toc22(A) info2("%s takes %.6f seconds.\n", A, myclockd()-tk)
#define toc3 (myclockd()-tk)

#define format2fn					\
    char fnstore[PATH_MAX]; char *fn=fnstore;		\
    va_list ap;						\
    va_start(ap,format);				\
    if(format) vsnprintf(fnstore,PATH_MAX, format, ap);	else fnstore[0]='\0';\
    va_end(ap);						\
    if(strlen(fnstore)==0) fn=NULL


#if    ( __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ > 4) ) && !defined(MATLAB_MEX_FILE)
#define CHECK_ARG(n) __attribute__ ((format (printf, n, n+1)))
#else   
#define CHECK_ARG(n)
#endif /* __GNUC__ */

#if defined(__GNUC__) && (__GNUC__ > 2) && defined(__OPTIMIZE__) 
#define LIKELY(A)   (__builtin_expect (A,1))
#define UNLIKELY(A) (__builtin_expect (A,0))
#else
#define LIKELY(A)   A
#define UNLIKELY(A) A
#endif /* __GNUC__ */

#if    __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4) && !defined(__INTEL_COMPILER)
#define CHECK_UNUSED_RESULT __attribute__((warn_unused_result))
#else
#define CHECK_UNUSED_RESULT
#endif /* __GNUC__ */

#if     __GNUC__ >= 4
#define CHECK_NULL_TERMINATED __attribute__((__sentinel__))
#else
#define CHECK_NULL_TERMINATED
#endif

#define READ_ENV_INT(A,min,max)				\
    if(getenv("MAOS_"#A)){				\
	A=strtol(getenv("MAOS_"#A),NULL,10);		\
	info2(#A"=%d\n", A);				\
	if(A>max || A<min){				\
	    error("MAOS_%s: invalid range\n", #A);	\
	}						\
    }
#define READ_ENV_DBL(A,min,max)				\
    if(getenv("MAOS_"#A)){				\
	A=strtod(getenv("MAOS_"#A),NULL);		\
	info2(#A"=%g\n", A);				\
	if(A>max || A<min){				\
	    error("MAOS_%s: invalid range\n", #A);	\
	}						\
    }
#define DEF_ENV_FLAG(A,default_val)			\
    static int A=default_val;				\
    static __attribute__((constructor)) void init(){	\
	READ_ENV_INT(A, 0, 1);				\
    }
#endif

