/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_SYS_COMMON_H
#define AOS_SYS_COMMON_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "mem.h"
#include <complex.h>
#include <signal.h>
#include <errno.h>
typedef double __complex__ dcomplex;
typedef double ddouble;//just for saving.
#include <stdarg.h>

#define format2fn					\
    char fnstore[PATH_MAX]; char *fn=fnstore;		\
    va_list ap;						\
    va_start(ap,format);				\
    vsnprintf(fnstore,sizeof(fnstore), format, ap);	\
    va_end(ap);						\
    if(strlen(fnstore)==0) fn=NULL

#ifndef USE_POSIX_SHM
#define USE_POSIX_SHM 0
#endif

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
#else
#define PRINT_BACKTRACE
#endif
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
    signal(SIGQUIT,func)
#define disable_signal_handler	\
    signal(SIGBUS, SIG_DFL);	\
    signal(SIGILL, SIG_DFL);	\
    signal(SIGSEGV,SIG_DFL);	\
    signal(SIGINT, SIG_DFL);	\
    signal(SIGPIPE,SIG_DFL);	\
    signal(SIGTERM,SIG_DFL);	\
    signal(SIGABRT,SIG_DFL);	\
    signal(SIGQUIT,SIG_DFL)


#endif
