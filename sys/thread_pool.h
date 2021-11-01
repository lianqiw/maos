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
#ifndef AOS_LIB_THREAD_POOL_H
#define AOS_LIB_THREAD_POOL_H
/**
   \file thread_pool.h
   Contains implementation of a thread pool. 
*/
#include "thread.h"
typedef struct thread_pool_t thread_pool_t;
void thread_pool_init(int nthread);
void thread_pool_queue(unsigned int *group, thread_wrapfun fun, void *arg, int njob, int urgent);
void thread_pool_wait(unsigned int *count, int urgent);
void thread_pool_wait_all(void);
void thread_pool_destroy(void);
#endif
