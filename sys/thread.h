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
#ifndef AOS_LIB_THREAD_H
#define AOS_LIB_THREAD_H
#include <stdio.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef USE_MEM
#if defined(__INTEL_COMPILER) || !defined(DEBUG) || defined(NDEBUG)
#define USE_MEM 0
#else
#define USE_MEM 1 
#endif
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
/**
   \file thread.h
   Functions regarding to threading
*/
/**
   Information about job to launch for each thread. start and end are the two indices.
*/
typedef struct thread_t thread_t;
typedef void *(*thread_fun)(void*);
typedef void (*thread_wrapfun)(thread_t*);
struct thread_t{
    long start;
    long end;
    long ithread;/*which thread this is. */
    long nthread;/*max number of threads*/
    thread_wrapfun fun;/*the function, takes data as argument */
    void *data;/*the data to pass to the function. */
};
/*
  For all the following calls, if urgent is 1, the job is queued in the front, otherwise in the end.

  CALL(fun,arg,nthread,urgent) executes "nthread" instances of fun with argument arg. 
  QUEUE is like CALL, but need an explicit WAIT
  
  The following use thread_t to manage the index
  QUEUE_THREAD(group, A, nthread, urgent) will queue arrays of thread_t (A) 
  CALL_THREAD(A, nthread, urgent) calls QUEUE_THREAD and waits for all to finish.
*/
#include <pthread.h>
#include "thread_pool.h"
#if _OPENMP >= 200805
/**
   Notice it is not effective to set the environment variables here.
*/
#define THREAD_YIELD	_Pragma("omp taskyield")
INLINE void THREAD_POOL_INIT(int nthread){
    fprintf(stderr, "Using OpenMP version %d\n", _OPENMP);
    omp_set_num_threads(nthread);
    omp_set_nested(0);//make sure nested is not enabled
}
INLINE void CALL(void*fun, void *arg, int nthread, int urgent){
    (void)urgent;
    for(int it=0; it<nthread; it++){
#pragma omp task untied
	((thread_fun)fun)(arg);
    }
#pragma omp taskwait
}

/*The following QUEUE, CALL, WAIT acts on function (fun) and argument (arg).*/
/*Don't turn the following into INLINE function becase task will be waited*/
#define QUEUE(group,fun,arg,nthread,urgent)	\
    (void) group; (void) urgent;		\
    for(int it=0; it<nthread; it++){		\
	_Pragma("omp task untied")		\
	fun(arg);				\
    }

#define WAIT(group)				\
    _Pragma("omp taskwait")


/*The following *_THREAD acts on thread_t array A*/

#define QUEUE_THREAD(group,A,urgent)		\
    (void)group;(void)urgent;			\
    for(int it=0; it<(A)[0].nthread; it++){	\
	if((A)[it].fun){			\
	    _Pragma("omp task untied")		\
		(A)[it].fun((A)+it);		\
        }					\
    }
#define WAIT_THREAD(group)			\
    _Pragma("omp taskwait")
/*Turn to inline function because nvcc concatenates _Pragma to } */
INLINE void CALL_THREAD(thread_t *A, int urgent){
    (void) urgent;
    for(int it=0; it<A[0].nthread; it++){		
	if(A[it].fun){
#pragma omp task untied
	    A[it].fun(A+it); 
	}
    }
#pragma omp taskwait
}

#else //using our thread_pool 

#define QUEUE(group,fun,arg,nthread,urgent)				\
    thread_pool_queue_many(&group, (thread_fun)fun, (void*)arg, nthread, urgent); \
    
#define CALL(fun,arg,nthread,urgent)					\
    if(nthread>1){							\
	long thgroup=0;							\
	thread_pool_queue_many						\
	    (&thgroup, (thread_fun)fun, (void*)arg, nthread, urgent);	\
	thread_pool_wait(&thgroup);					\
    }else{								\
	fun(arg);							\
    }

#define WAIT(group)\
    thread_pool_wait(&group);
/**
   Queue jobs to group. Do not wait
*/
#define QUEUE_THREAD(group,A,urgent)			\
    if((A[0].nthread)>1){						\
	thread_pool_queue_many(&group,NULL,A,A[0].nthread,urgent); \
    }else{							\
	(A)->fun(A);						\
    }

/**
   Queue jobs to a temp group. Wait for complete.
*/
#define CALL_THREAD(A,urgent)			\
    if((A[0].nthread)>1){			\
	long group=0;				\
	QUEUE_THREAD(group,A,urgent);		\
	WAIT_THREAD(group);			\
    }else{					\
	(A)->fun(A);				\
    }
/**
   Wait for all jobs in group to finish.
*/
#define WAIT_THREAD(group) thread_pool_wait(&group)

#define THREAD_POOL_INIT(A) ({thread_pool_init(A);fprintf(stderr, "Using thread pool\n");})
#define THREAD_YIELD thread_pool_do_job_once()
#endif

#define LOCK(A) pthread_mutex_lock(&A)
#define UNLOCK(A) pthread_mutex_unlock(&A)
#define PINIT(A) pthread_mutex_init(&A,NULL)
#define PDEINIT(A) pthread_mutex_destroy(&A)
#define PNEW(A) static pthread_mutex_t A=PTHREAD_MUTEX_INITIALIZER
#define PNEW2(A) pthread_mutex_t A=PTHREAD_MUTEX_INITIALIZER

extern pthread_mutex_t mutex_fftw;
#define LOCK_FFT LOCK(mutex_fftw)
#define UNLOCK_FFT UNLOCK(mutex_fftw)

void thread_prep(thread_t *info, long start, long end, long nthread, 
		 thread_wrapfun fun, void *data);

#define LOCKADD(dest,src,step)\
    (({__asm__ __volatile__ ("lock; xaddl %0,%1"	\
			   : "=r" (dest), "=m" (src)	\
			    : "0" (step), "m" (src));}),dest)

int lockadd(int *src, int step);

#define SPIN_LOCK(i) while(__sync_lock_test_and_set(&i, 1)) while(i)
#define SPIN_UNLOCK(i) __sync_lock_release(&i)
/**
   Create a new thread and let it go.
*/
INLINE void thread_new(thread_fun fun, void* arg){
    pthread_t temp;
    pthread_create(&temp, NULL, fun, arg);
}
void thread_block_signal();
INLINE int cmpxchg(int *ptr, int old, int newval){
    volatile int *__ptr = (volatile int *)(ptr);	
    int __ret;                                     
    __asm__ volatile( "lock; cmpxchg %2,%1"
		      : "=a" (__ret), "+m" (*__ptr) 
		      : "r" (newval), "0" (old)                     
		      : "memory");				 
    return __ret;
}
INLINE int atomicadd(int *ptr, int val){
    int old;
    do{
	old=*ptr;
    }while(cmpxchg(ptr, old, old+val)!=old);
    return old+val;
}
#endif
