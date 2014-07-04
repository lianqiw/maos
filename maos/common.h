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

#ifndef __MAOS_COMMON_H
#define __MAOS_COMMON_H
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <math.h>
#include <string.h>
#include <search.h>

#include "../lib/aos.h"
#include "parms.h"
#include "types.h"
#include "utils.h"
extern double TOMOSCALE;
extern int exit_success;
extern const char* dirsetup;
extern const char *dirskysim;
#define EXIT raise(SIGTERM)
extern GLOBAL_T *global;
#define adpind(A,i) ((A)->nx>1?(A)->p[i]:(A)->p[0])
#define CALL_ONCE\
    {static int count=0; count++; if(count>1) warning("This function should only be called once\n");}
void maos(const PARMS_T *parms);
#endif

