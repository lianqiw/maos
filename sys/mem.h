/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifdef __cplusplus
extern "C"{
#endif
    void *calloc_maos(size_t, size_t);
    void *malloc_maos(size_t);
    void *realloc_maos(void *, size_t);
    void  free_maos(void *);
#ifdef __cplusplus
}
#endif
#define myalloca(nelem, type) (type*)alloca(nelem*sizeof(type))
#define mycalloc(nelem, type) (type*)calloc(nelem,sizeof(type))
#define mymalloc(nelem, type) (type*)malloc(nelem*sizeof(type))
#define myrealloc(p, nelem, type) (type*)realloc(p,nelem*sizeof(type))
#if defined(__linux__) || !defined(__cplusplus)
#ifndef IN_MEM_C
#define free free_maos
#define malloc malloc_maos
#define calloc calloc_maos
#define realloc realloc_maos
#endif
#endif
void register_deinit(void (*fun)(void), void *data);
void malloc_dbg_enable();
int malloc_dbg_disable(int print);
void print_backtrace();
#endif
