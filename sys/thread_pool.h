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
#ifndef AOS_LIB_THREAD_POOL_H
#define AOS_LIB_THREAD_POOL_H
typedef struct thread_pool_t thread_pool_t;
void thread_pool_init(int nthread);
void thread_pool_queue(long *count, void *(*fun)(void*), void *arg, int urgent);
void thread_pool_queue_many(long *group, thread_fun fun, void *arg, int njob, int urgent);
void thread_pool_wait(long *count);
void thread_pool_wait_all(void);
void thread_pool_destroy(void);
void thread_pool_do_urgent_job(void);
int thread_pool_do_job_once(void);
#endif
