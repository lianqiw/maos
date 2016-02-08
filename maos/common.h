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

#ifndef __MAOS_COMMON_H
#define __MAOS_COMMON_H
#include <sys/time.h>
#include <sys/times.h>
#include <search.h>
#include "../lib/aos.h"
#include "parms.h"
#include "types.h"
#include "utils.h"
extern int exit_fail;
extern const char *dirskysim;
#define EXIT raise(SIGTERM)
extern GLOBAL_T *global;
#define adpind(A,i) ((A)->nx>1?(A)->p[i]:(A)->p[0])
#endif

