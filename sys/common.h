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
/**
   \file sys/common.h
   Every source file in this folder should include this file
*/

#ifndef AOS_COMMON_H
#define AOS_COMMON_H
//C standard headers

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if _OPENMP >= 201511
#include <omp.h>
#else
#undef _OPENMP
#endif

#include <stdarg.h> //for va_start
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h> //don't use tgmath here. it is included in math/numtype.h
/*#if defined(__cplusplus) && !defined(AOS_CUDA_H)
//c++ mode, not CUDA
using std::signbit;
using std::isfinite;
using std::isnan;
using std::strerror;
#endif*/ //C99 mode or CUDA.
#if defined(__cplusplus)
#define __auto_type auto
#endif
#undef	MAX
#define MAX(a,b) ({__typeof__(a) _M1=(a); __typeof__(b) _M2=(b); (_M1)>(_M2)?(_M1):(_M2);})
#undef	MIN
#define MIN(a,b) ({__typeof__(a) _m1=(a); __typeof__(b) _m2=(b); (_m1)<(_m2)?(_m1):(_m2);})
#define RSS(a,b) ({__typeof__(a) _M1=(a); __typeof__(b) _M2=(b); sqrtf(_M1*_M1+_M2*_M2);})
#define CLIP(x,d_,u_) ({__typeof__(d_) d=(d_); __typeof__(u_) u=(u_);x=x<d?d:(x>u?u:x);})
#ifdef __linux__
#include <linux/limits.h> /*includes definition of PATH_MAX */
#else
#include <limits.h>
#endif/*__linux__ */
#include <signal.h>
#include "mem.h"

#ifndef restrict
#define restrict __restrict
#endif

#undef EPS
#if CPU_SINGLE
#define EPS 1.e-6
#else
#define EPS 1.e-15
#endif
/**
   Obtain the basename of a file. Has the GNU basename behavior.
*/
static inline const char *mybasename(const char *fn){
  if(fn){
    const char *slash=strrchr(fn, '/');
    if(slash){
      if(slash[1]=='\0') return "";else return slash+1;
    } else return fn;
  } else return "";
}
#define BASEFILE mybasename(__FILE__) //(strrchr(__FILE__, '/') ?strrchr(__FILE__, '/')+1  : __FILE__)

/*
  use () to make the statements a single statement.
*/
extern int detached;
#define BLACK "\033[00;00m"
#define RED "\033[01;31m"
#define GREEN "\033[0;32m"
#define YELLOW "\033[0;33m"
#define BLUE "\033[0;34m"
#define MAGENTA "\033[0;35m"
#define CYAN "\033[0;36m"

extern int LOG_LEVEL;//default is 0; override with MAOS_LOG_LEVEL; higher value has more output
extern int signal_caught;
extern FILE* fplog;//The output to fplog is always without color unless user specified.
#define logstd(level, A...) 		      ({if(LOG_LEVEL>level){if(!detached)fprintf(stdout, A); if(fplog){fprintf(fplog, A);}}})
#define logerr(level, COLOR, format, ...) ({if(LOG_LEVEL>level){if(!detached)fprintf(stderr, COLOR format BLACK, ##__VA_ARGS__); if(fplog){fprintf(fplog, format, ##__VA_ARGS__);}}})
#define logdbg(level, COLOR, format, ...) ({if(LOG_LEVEL>level){if(!detached)fprintf(stdout, COLOR format BLACK, ##__VA_ARGS__); if(fplog){fprintf(fplog, format, ##__VA_ARGS__);}}})

#define error(format,...)      ({logerr(-4, RED,         "Error(%s:%d): " format, BASEFILE,__LINE__, ##__VA_ARGS__); default_signal_handler(SIGUSR2,0,0);})
#define warning(format,...)      logerr(-4, CYAN,      "Warning(%s:%d): " format, BASEFILE,__LINE__, ##__VA_ARGS__)
#define warning_time(format,...) logerr(-4, CYAN, "[%s] Warning(%s:%d): " format, myasctime(0),BASEFILE,__LINE__, ##__VA_ARGS__)
#define warning_once(A...)  ({static int done=0; if(!done){done=1; warning(A);}})
//all info are shown at default log level
#define info_line(format,...) logstd(-4,      "Info(%s:%d): " format, BASEFILE,__LINE__,##__VA_ARGS__)
#define info_time(format,...) logstd(-1, "[%s] Info(%s:%d): " format, myasctime(0), BASEFILE,__LINE__,##__VA_ARGS__)
#define info_green(A...) logdbg(-2, GREEN, A)
#define info(A...)  logstd(-1, A) //least important info
#define info2(A...) logstd(-2, A)
#define info3(A...) logstd(-3, A) //most important info
#define info_once(A...) ({static int done=0; if(!done){done=1; info(A);}})
#define info_progress(i,n) if((i)%(((n)>>4)+1)==0) fprintf(stderr,">") //;/*if((i)+1==(n)) fprintf(stderr,"\n");*/})
#define info_errno(A) if(errno) info(A " failed (%d): %s\n", errno, strerror(errno))
//dbg are not shown at default log level
#define dbg( A...) logdbg(0, YELLOW, A)//most important dbg
#define dbg2(A...) logdbg(1, YELLOW, A)
#define dbg3(A...) logdbg(2, YELLOW, A)//least important dbg
#define logdbg_time(level, format, ...) logdbg(level, YELLOW, "[%s] %s: " format, myasctime(0), __func__, ##__VA_ARGS__)
#define dbg_line(format,...) logdbg(0, YELLOW, "Debug(%s:%d): " format,BASEFILE,__LINE__,##__VA_ARGS__)
#define dbg_time( A...) logdbg_time(0, A)
#define dbg2_time(A...) logdbg_time(1, A)
#define dbg3_time(A...) logdbg_time(2, A)
#define dbg_once(A...) ({static int done=0; if(!done){done=1; dbg_line(A);}})
#ifndef assert
#if DEBUG
#define assert(A) if(!(A)) error("assertion failed: %s\n", #A)
#else
#define assert(A)
#endif
#endif
#if DEBUG
#define check(A) ((A)?1:(error("check failed: %s\n", #A),0))
#else
#define check(A) (A)
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
#ifdef tk
#undef tk
#endif
#define TIC double tk
#define tic tk=myclockd()
#define toc(format,...) logstd(-1, format " takes %.6f seconds.\n", ##__VA_ARGS__, myclockd()-tk)
#define toc2(format,...) dbg(format " takes %.6f seconds.\n", ##__VA_ARGS__, myclockd()-tk)
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
//check whether endptr only contains space. endstr marks the end of str (\0)
static inline int check_space(const char *endptr, const char *endstr){
   while(endptr<endstr && endptr[0]==' ') endptr++;
   return endptr!=endstr;
}
#define STR_TO_INT(A,B) strtol(A,B,10)
#define STR_TO_DBL(A,B) strtod(A,B)
#define READ_ENV_NUM(A,min,max,T,FUN)\
    if(getenv("MAOS_"#A)){\
      char *endptr=0; \
      const char *str=getenv("MAOS_"#A);\
      const char *endstr=str+strlen(str);\
		T B=(T)FUN(str,&endptr); \
      if(endptr==str || check_space(endptr, endstr) || (max>min && B>max) || B<min){  \
        error("MAOS_"#A"={%s} is invalid. Must be within [%g, %g].\n", str, (double) min, (double) max);\
      }else if(A!=B){dbg(#A" changed to %g\n", (double)B); A=B; } }

#define READ_ENV_INT(A,min,max)	READ_ENV_NUM(A,min,max,int,   STR_TO_INT)
#define READ_ENV_DBL(A,min,max)	READ_ENV_NUM(A,min,max,double,STR_TO_DBL)

#define DEF_ENV_FLAG(A,default_val)			\
    static int A=default_val;				\
    static __attribute__((constructor)) void init(){\
		READ_ENV_INT(A, 0, INT_MAX);\
    }

#define FREE(A) if(A){free(A); A=NULL;}
#endif
