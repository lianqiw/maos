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

#ifndef AOS_COMMON_H
#define AOS_COMMON_H

/**
   Every source file in this folder should include this file
*/

#include "sys/sys.h"

    
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <alloca.h>
#ifdef __linux__
//includes definition of PATH_MAX
#include <linux/limits.h>
#endif
#ifndef USE_DAEMON
#define USE_DAEMON 1
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
/*
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>   
#include <signal.h>
#include <string.h>
#include <time.h>*/
#include "misc.h"
#define error_write error("Write failed\n")
#define error_read error("Read failed\n")
#ifndef __thread
#define __thread
#endif
/**
   Functions that return realtime:
   time(): resolution is in integer second. not enough resolution.
   gettimeofday(): must use struct as input argument.
   Functions that return tics:
   times(): the argument can not be NULL in Mac. obsolete by gettimeofday in Mac
   clock(): returns processor time, not real time.

   Decision: use customize function myclockd() in utils.c which returns a double for second. 
*/
#define TIC __thread double tk
#define tic tk=myclockd();
#define toc(A...) ({info(A);fprintf(stderr," takes %.6f seconds.\n",myclockd()-tk);})
#define toc2(A...) ({fprintf(stderr,A);fprintf(stderr," takes %.6f seconds.\n",myclockd()-tk);})
//define toc2(A) info2("%s: takes %.6f seconds.\n",  A,myclockd()-tk)
#define toc3 myclockd()-tk
#endif

