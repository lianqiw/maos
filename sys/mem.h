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
#ifndef __UTILS_MEM_H
#define __UTILS_MEM_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

extern void *(*CALLOC)(size_t nmemb, size_t size);
extern void *(*MALLOC)(size_t size);
extern void *(*REALLOC)(void*p0, size_t size);
extern void  (*FREE)(void *p);

#define malloc MALLOC
#define calloc CALLOC
#define realloc REALLOC
#define free FREE

void register_deinit(void (*fun)(void), void *data);
void malloc_dbg_enable();
void malloc_dbg_disable();
extern int exit_success;
#endif
