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

/**
   \file sys/common.h
   Every source file in this folder should include this file
*/

#ifndef AOS_COMMON_H
#define AOS_COMMON_H
//C standard headers
#include <unistd.h>
#include <time.h>

typedef void (*quitfun_t)(const char*);
extern quitfun_t quitfun;
void default_quitfun(const char* msg);

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

#undef EPS
#if CPU_SINGLE
#define EPS 1.e-6
#else
#define EPS 1.e-15
#endif

#define BASEFILE (strrchr(__FILE__, '/') ?strrchr(__FILE__, '/')+1  : __FILE__)

/*
  use () to make the statements a single statement.
*/
extern int detached;
#define BLACK (detached?"":"\033[00;00m")
#define RED (detached?"":"\033[01;31m")
#define GREEN (detached?"":"\033[0;32m")
#define YELLOW (detached?"":"\033[0;33m")
#define BLUE (detached?"":"\033[0;34m")
#define MAGENTA (detached?"":"\033[0;35m")
#define CYAN (detached?"":"\033[0;36m")
#define QUIT_FUN(A) quitfun?quitfun(A):default_quitfun(A);

extern int LOG_LEVEL;//default is 0; override with MAOS_LOG_LEVEL; higher value has more output
extern FILE* fpconsole;
#define logstd(level, A...) ({if(LOG_LEVEL>level){fprintf(stdout, A);}})
#define logerr(level, A...) ({if(LOG_LEVEL>level){fprintf(stdout, A);fprintf(stderr, A);}})
#define error(format,...) ({logerr(-4, "%sError(%s:%d): " format "%s", RED,BASEFILE, __LINE__, ##__VA_ARGS__, BLACK); QUIT_FUN("Error happened");})
#define warning(format,...) logstd(-4, "%sWarning(%s:%d): " format "%s", YELLOW, BASEFILE, __LINE__,##__VA_ARGS__, BLACK)
#define warning_time(format,...) logstd(-4,"%s[%s] (%s:%d): " format "%s", YELLOW,myasctime(),BASEFILE, __LINE__,  ##__VA_ARGS__, BLACK)
#define warning_once(A...)  ({static int done=0; if(!done){done=1; warning(A);}})
//all info are shown at default log level
#define info(A...)  logstd(-1, A) //least important info
#define info2(A...) logstd(-2, A)
#define info3(A...) logstd(-3, A) //most important info
#define info_console(A...) ({if(LOG_LEVEL>-2 && fpconsole) fprintf(fpconsole, A);}) //stderr is not directed to file.
#define info_once(A...) ({static int done=0; if(!done){done=1; info(A);}})
//dbg are not shown at default log level
//dbg do not shown when detached
//use __func__ to indicate function name
#define logdbg(level, format, ...) logstd(level, "%s" format "%s", CYAN, ##__VA_ARGS__, BLACK)
#define dbg( A...) logdbg(0, A)
#define dbg2(A...) logdbg(1, A)
#define dbg3(A...) logdbg(2, A)
#define logdbg_time(level, format, ...) logstd(level, "%s[%s]: " format "%s", CYAN, myasctime(), ##__VA_ARGS__, BLACK)
#define dbg_time( A...) logdbg_time(0, A)
#define dbg2_time(A...) logdbg_time(1, A)
#define dbg3_time(A...) logdbg_time(2, A)

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
#define toc(format,...) info(format " takes %.6f seconds.\n", ##__VA_ARGS__, myclockd()-tk)
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

#define READ_ENV_INT(A,min,max)					\
    if(getenv("MAOS_"#A)){						\
		int B=strtol(getenv("MAOS_"#A),NULL,10);\
		if(A!=B){ 								\
			A=MIN(MAX(B, min), max);			\
			dbg(#A" changed to %d\n", A);		\
		}										\
    }
#define READ_ENV_DBL(A,min,max)					\
    if(getenv("MAOS_"#A)){						\
		double B=strtod(getenv("MAOS_"#A),NULL);	\
		if(A!=B){								\
			A=MIN(MAX(B, min), max);			\
			dbg(#A" changed to %g\n", A);		\
		}								        \
    }
#define DEF_ENV_FLAG(A,default_val)			\
    static int A=default_val;				\
    static __attribute__((constructor)) void init(){	\
		READ_ENV_INT(A, 0, 1);				\
    }
/*#define DEF_ENV_FLAG_LOCAL(A, default_val, min_val, max_val)	\
    static int A=-1;						\
    if(A==-1){					 			\
		A=default_val;						\
		READ_ENV_INT(A, min_val, max_val);	\
    }*/
#endif

