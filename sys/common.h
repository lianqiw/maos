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

/**
   \file common.h
   Every source file in this folder should include this file
*/

#ifndef AOS_COMMON_H
#define AOS_COMMON_H
//C standard headers
#include <unistd.h>
#include <time.h>

typedef void (*quitfun_t)(const char*);
extern quitfun_t quitfun;
void default_quitfun(const char *msg);

#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif

#if !defined(__FreeBSD__) && !defined(__NetBSD__)
#include <alloca.h>
#endif


#if defined(__cplusplus) && !defined(AOS_CUDA_GPU_H)
//c++ mode, not CUDA
#include <csignal>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using std::signbit;
using std::isfinite;
using std::isnan;
using std::strerror;
#else//C99 mode or CUDA.
#include <signal.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#endif //if defined(__cplusplus) && !defined(AOS_CUDA_GPU_H)

//GNU GCC changes definition of inline to C99 compatible since 4.4
#if __GNUC__ == 4 && __GNUC_MINOR__ < 5 && !defined(__clang__)
#define INLINE static inline __attribute__((gnu_inline, always_inline)) //GNU
#else
#define INLINE static inline __attribute__((always_inline)) //C99
#endif //if __GNUC__ == 4 && __GNUC_MINOR__ < 5

#undef	MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#undef	MIN
#define MIN(a,b) (((a)<(b))?(a):(b))

#include "mem.h"
#ifdef __linux__
#include <linux/limits.h> /*includes definition of PATH_MAX */
#else
#include <limits.h>
#endif/*__linux__ */

#ifndef restrict
#define restrict __restrict
#endif

#ifndef EPS
#define EPS 1.e-15
#endif

#define BASEFILE (strrchr(__FILE__, '/') ?strrchr(__FILE__, '/')+1  : __FILE__)

/*
  use () to make the statements a single statement.
*/
#ifndef error
#define QUIT_FUN(A) quitfun?quitfun(A):default_quitfun(A);
#define info(A...) ({char fline[4096]; int n__;			      \
	    snprintf(fline,4096, "INFO(%s:%d): ", BASEFILE, __LINE__); \
	    n__=strlen(fline); snprintf(fline+n__, 4096-n__-1, A);     \
	    fprintf(stderr,"%s", fline); })

#define error(A...) ({char fline[4096]; int n__;			\
	    snprintf(fline,4096, "\033[01;31mFATAL(%s:%d): ", BASEFILE, __LINE__); \
	    n__=strlen(fline); snprintf(fline+n__, 4096-n__-1, A);	\
	    n__=strlen(fline); snprintf(fline+n__, 4096-n__-1, "\033[00;00m"); \
	    QUIT_FUN(fline);})

#define warning(A...) ({char fline[4096]; int n__;			\
	    snprintf(fline,4096, "\033[01;31mWARN(%s:%d): ", BASEFILE, __LINE__); \
	    n__=strlen(fline); snprintf(fline+n__, 4096-n__-1, A);	\
	    fprintf(stderr,"%s\033[00;00m", fline); })

#define info2(A...) fprintf(stderr, A)
#define warning2(A...) ({fprintf(stderr,"\033[00;31mWarning:\033[00;00m" A);})

#define info_time(A...) ({char fline[4096]; int n__;			      \
	    snprintf(fline,4096, "INFO(%s:%d)[%s]: ", BASEFILE, __LINE__, myasctime()); \
	    n__=strlen(fline); snprintf(fline+n__, 4096-n__-1, A);     \
	    fprintf(stderr,"%s", fline); })

#define warning_time(A...) ({char fline[4096]; int n__;			\
	    snprintf(fline,4096, "\033[01;31mWARN(%s:%d)[%s]: ", BASEFILE, __LINE__, myasctime()); \
	    n__=strlen(fline); snprintf(fline+n__, 4096-n__-1, A);	\
	    fprintf(stderr,"%s\033[00;00m", fline); })

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

/*
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

