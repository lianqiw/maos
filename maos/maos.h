/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef __MAOS_H
#define __MAOS_H
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <math.h>
#include <string.h>
#include <search.h>
#include <assert.h>

#include "../lib/aos.h"
#include "parms.h"
#include "types.h"
#include "utils.h"
extern int exit_success;
extern char* dirsetup;
extern char *dirskysim;
#define EXIT raise(SIGUSR1)
void maos(const PARMS_T *parms);
extern const PARMS_T *curparms;
extern int curiseed;
#define adpind(A,i) ((A)->nx>1?(A)->p[i]:(A)->p[0])
#endif

