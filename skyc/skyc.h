/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef __AOS_SKYC_H
#define __AOS_SKYC_H
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <math.h>
#include <string.h>
#include <search.h>
#include <assert.h>
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#ifndef restrict
#define restrict __restrict
#endif
#define AOS
#if defined(DEBUG) && DEBUG==1
#define DEBUG_SAVE 1
#else
#define DEBUG_SAVE 1
#endif
#define skyc_done(A) ({scheduler_finish(A);skyc_signal_handler(A);exit_success=1;if(A!=0) exit(A);})
#include "../lib/aos.h"
extern char *dirsetup;
extern char *dirstart;
#endif
