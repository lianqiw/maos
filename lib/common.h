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

#ifndef AOS_COMMON_H
#define AOS_COMMON_H

/**
   Every source file in this folder should include this file
*/

#ifdef HAVE_CONFIG_H
#include "config.h" //not a good idea to include HAVE_CONFIG_H here
#endif
#include <signal.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <alloca.h>
#include <time.h>
typedef double __complex__ dcomplex;
typedef double ddouble;/*just for saving.*/
#ifndef __CYGWIN__
#include <complex.h>
#else
//CYGWIN does not have complex.h.
#define complex __complex__
#define _Complex_I (__extension__ 1.0iF)
#define I _Complex_I
double cimag(dcomplex __z);
double creal(dcomplex __z);
dcomplex conj(dcomplex __z);
#define cexp(z) exp(creal(z))*(cos(cimag(z))+I*sin(cimag(z)))
dcomplex cpow(dcomplex x, dcomplex z);
dcomplex csqrt(dcomplex);
dcomplex clog(dcomplex);
#endif
#include "sys/mem.h"
#ifdef __linux__
#include <linux/limits.h> //includes definition of PATH_MAX
#endif//__linux__
#ifndef restrict
#define restrict __restrict
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#ifndef EPS
#define EPS 1.e-14
#endif

/**
   2010-01-03
   USE_MYFLOOR = 1 reduce timing by 2s
*/
#define USE_MYFLOOR 1
#if USE_MYFLOOR
#define ifloor(A) (A<0?(int)(A)-1:(int)(A))
#else
#define ifloor(A) (int)floor(A)
#endif
#define iceil(A) (int)ceil(A)

#if defined(FP_FAST_FMA)
#define myfma fma
#error FP_FAST_FMA is defined . remove.
#else
#define myfma(x,y,z) (x)*(y)+z
#endif
#define SPLIT(A,B,C) {C=ifloor(A); B=A-C;}

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
#define toc(A...) ({info(A);fprintf(stderr," takes %.6f seconds.\n",myclockd()-tk);})
#define toc2(A...) ({fprintf(stderr,A);fprintf(stderr," takes %.6f seconds.\n",myclockd()-tk);})
#define toc3 myclockd()-tk

#define format2fn					\
    char fnstore[PATH_MAX]; char *fn=fnstore;		\
    va_list ap;						\
    va_start(ap,format);				\
    vsnprintf(fnstore,sizeof(fnstore), format, ap);	\
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


#if USE_MEM == 1
void print_backtrace(int sig);
#define PRINT_BACKTRACE print_backtrace(0);
#else//if USE_MEM
#define PRINT_BACKTRACE
#endif//if USE_MEM
#include <string.h>


#define BASEFILE (strrchr(__FILE__, '/') ?strrchr(__FILE__, '/')+1  : __FILE__)


/*
  use () to make the statements a single statement.
*/
#ifndef error
#define error(A...) ({char fline[80];					\
	    snprintf(fline, 80,"%s:%d",BASEFILE,__LINE__);		\
	    fprintf(stderr, "\033[01;31m%-20s Fatal error: ",fline);	\
	    fprintf(stderr, A); fprintf(stderr,"\033[00;00m");		\
	    PRINT_BACKTRACE;						\
	    raise(SIGTERM);})
#define info(A...) ({char fline[80];				\
	    snprintf(fline, 80,"%s:%d",BASEFILE,__LINE__);	\
	    fprintf(stderr, "%-20s",fline);			\
	    fprintf(stderr, A);})
#define warning(A...) ({char fline[80];				\
	    snprintf(fline, 80,"%s:%d",BASEFILE,__LINE__);	\
	    fprintf(stderr,"\033[01;31m%-20s", fline);		\
	    fprintf(stderr,A);fprintf(stderr,"\033[00;00m");})

#define error2(A...) ({							\
	    fprintf(stderr, "\033[01;31mFatal error\033[00;00m\t");	\
	    fprintf(stderr, A);						\
	    PRINT_BACKTRACE;						\
	    raise(SIGTERM);})
#define info2(A...) fprintf(stderr, A)
#define warning2(A...) ({					\
	    fprintf(stderr,"\033[00;31m");			\
	    fprintf(stderr,A);fprintf(stderr,"\033[00;00m"); }) 

#define error3(A...) ({char fline[80];					\
	    snprintf(fline, 80,"%s:%d",BASEFILE,__LINE__);		\
	    fprintf(stderr, "[%s]\033[01;31m%-20s Fatal error:",myasctime(),fline); \
	    fprintf(stderr, A); fprintf(stderr,"\033[00;00m");		\
	    PRINT_BACKTRACE;						\
	    raise(SIGTERM);})
#define warning3(A...) ({char fline[80];				\
	    snprintf(fline, 80,"%s:%d",BASEFILE,__LINE__);		\
	    fprintf(stderr,"[%s]\033[01;31m%-20s", myasctime(),fline);	\
	    fprintf(stderr,A);fprintf(stderr,"\033[00;00m");})
#define info3(A...) ({char fline[80];				\
	    snprintf(fline, 80,"%s:%d",BASEFILE,__LINE__);	\
	    fprintf(stderr, "[%s]%-20s",myasctime(),fline);	\
	    fprintf(stderr, A);})
#endif

#define register_signal_handler(func)	\
    signal(SIGBUS, func);		\
    signal(SIGILL, func);		\
    signal(SIGSEGV,func);		\
    signal(SIGINT, func);		\
    signal(SIGPIPE,SIG_IGN);		\
    signal(SIGTERM,func);		\
    signal(SIGABRT,func);		\
    signal(SIGUSR1,func);		\
    signal(SIGQUIT,func)

#define disable_signal_handler	\
    signal(SIGBUS, SIG_IGN);	\
    signal(SIGILL, SIG_IGN);	\
    signal(SIGSEGV,SIG_IGN);	\
    signal(SIGINT, SIG_IGN);	\
    signal(SIGPIPE,SIG_IGN);	\
    signal(SIGTERM,SIG_IGN);	\
    signal(SIGABRT,SIG_IGN);	\
    signal(SIGUSR1,SIG_IGN);	\
    signal(SIGQUIT,SIG_IGN)


#endif

