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
#ifndef AOS_LIB_THREAD_H
#define AOS_LIB_THREAD_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef USE_PTHREAD
#if MATLAB_MEX_FILE
#define USE_PTHREAD 1
#else
#define USE_PTHREAD 0
#endif
#endif
/**
   Information about job to launch for each thread. start and end are the two indices.
*/
#ifndef AOS_ACCPHI_H
typedef struct thread_t thread_t;
#endif
typedef void *(*thread_fun)(void*);
typedef void (*thread_wrapfun)(thread_t*);
struct thread_t{
    long start;
    long end;
    long ithread;//which thread this is.
    thread_wrapfun fun;//the function, takes data as argument
    void *data;//the data to pass to the function.
};

#if USE_PTHREAD //Always use thread_pool
#include <pthread.h>
#include "thread_pool.h"
#define CALL(A,B,nthread)				\
    if(nthread>1){					\
	long thgroup=0;					\
	for(int ithread=0; ithread<nthread; ithread++){	\
	    thread_pool_queue(&thgroup,			\
			      (thread_fun)A,	\
			      (void*)B,1);		\
	}						\
	thread_pool_wait(&thgroup);			\
    }else{						\
	A(B);						\
    }
#define CALL_EACH(A,B,nthread)				\
    if(nthread>1){					\
	long thgroup=0;					\
	for(int ithread=nthread-1; ithread>-1; ithread--){	\
	    thread_pool_queue(&thgroup,			\
			      (thread_fun)A,	\
			      (void*)&(B[ithread]),1);	\
	}						\
	thread_pool_wait(&thgroup);			\
    }else{						\
	A(B);						\
    }
/**
   Queue jobs to group. Do not wait
*/
#define QUEUE_THREAD(group,A,nthread,urgent)			\
    if((nthread)>1){						\
	thread_pool_queue_many(&group,A,nthread,urgent);	\
    }else{							\
	(A)->fun(A);						\
    }
/**
   Wait for all jobs in group to finish.
 */
#define WAIT_THREAD(group) thread_pool_wait(&group)
#define THREAD_POOL_INIT(A) thread_pool_create(A)
#else
#define CALL(A,B,nthread) A(B) //no threading
#define CALL_DATAEACH(A,B,nthread) A(B)
#define QUEUE_THREAD(group,A,nthread) A->fun(A)
#define THREAD_POOL_INIT(A)  //Do nothing
#endif

/**
   Call and wait for them to finish.
*/
#define CALL_THREAD(A,nthread,urgent)		\
    if((nthread)>1){				\
	long group=0;				\
	QUEUE_THREAD(group,A,nthread,urgent);	\
	WAIT_THREAD(group);			\
    }else{					\
	(A)->fun(A);				\
    }


#if USE_PTHREAD > 0
#define LOCK(A) pthread_mutex_lock(&A)
#define UNLOCK(A) pthread_mutex_unlock(&A)
#define PINIT(A) pthread_mutex_init(&A,NULL)
#define PDEINIT(A) pthread_mutex_destroy(&A)
#define PNEW(A) static pthread_mutex_t A=PTHREAD_MUTEX_INITIALIZER
#else
#define LOCK(A)
#define UNLOCK(A)
#define PINIT(A)
#define PDEINIT(A)
#define PNEW(A)
#endif


void thread_prep(thread_t *info, long start, long tot, long nthread, 
		 thread_wrapfun fun, void *data);
#endif
