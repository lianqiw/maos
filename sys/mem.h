/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
  extern void *(*calloc_default)(size_t, size_t);
  extern void *(*malloc_default)(size_t);
  extern void *(*realloc_default)(void *, size_t);
  extern void  (*free_default)(void *);
#ifdef __cplusplus
}
#endif
#include <signal.h>
extern int MEM_VERBOSE;
extern int MEM_DEBUG;
extern int MEM_FUNTRACE;
extern __thread char funtrace[];//stores info about top level function ca
#define funtrace_len 64
#define funtrace_set (((MEM_DEBUG || MEM_VERBOSE) && !funtrace[0] && MEM_FUNTRACE)?(void*)(long)snprintf(funtrace, funtrace_len, "%s:%d (%s)", BASEFILE,__LINE__,__func__):NULL)
#define funtrace_unset funtrace[0]=0;
extern void  (*free_default)(void *);
//Do not convert alloca to function call. It is auto freed after function returns.
#define myalloca(nelem, type)     (type*)alloca(nelem*sizeof(type))
#define mycalloc(nelem, type)     (type*)(funtrace_set, calloc_maos(nelem, sizeof(type)))
#define mymalloc(nelem, type)     (type*)(funtrace_set, malloc_maos(nelem*sizeof(type)))
#define myrealloc(p, nelem, type) (type*)(funtrace_set, realloc_maos(p, nelem*sizeof(type)))
#define myfree(p) if(p) free_maos(p)

#if !defined(IN_MEM_C) && !defined(AOS_CUDA_H)
#define calloc(nelem, size) (funtrace_set, calloc_maos(nelem, size))
#define malloc(size)        (funtrace_set, malloc_maos(size))
#define realloc(p, size)    (funtrace_set, realloc_maos(p, size))
#define free(p) if(p) free_maos(p)
#endif
void register_malloc(void* (*ex_malloc)(size_t), void* (*ex_calloc)(size_t, size_t), void* (*ex_realloc)(void*, size_t), void (*ex_free)(void*));
void register_deinit(void (*fun)(void), void *data);
void read_sys_env();
///Check whether signal is crash
#define iscrash(sig) (sig==SIGABRT||sig==SIGSEGV||sig==SIGILL||sig==SIGFPE)
void register_signal_handler(int(*)(int));
int dummy_signal_handler(int sig);
void default_signal_handler(int sig, siginfo_t *siginfo, void *unused);
#endif
