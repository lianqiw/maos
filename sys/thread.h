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
#ifndef AOS_LIB_THREAD_H
#define AOS_LIB_THREAD_H
/**
   \file thread.h
   Functions regarding to threading.

   Openmp version:
   3.0 (200805): Introduced task
   4.0 (201307): Introduced taskgroup
   4.5 (201511): Introduced taskloop and priority
*/

#include "common.h"
#define DO_PRAGMA(A...) _Pragma(#A)

#if HAS_OPENMP == 1 && !defined(_OPENMP)
#error "-fopenmp is not properly passed"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#if _OPENMP >= 200805
#define OMPTASK_SINGLE				\
    DO_PRAGMA(omp parallel)			\
    DO_PRAGMA(omp single)			
#else
#define OMPTASK_SINGLE
#endif
#if _OPENMP >= 200805
#define OMP_IN_PARALLEL omp_in_parallel()
#else
#define OMP_IN_PARALLEL 0
#endif

#if _OPENMP >= 201307
#define OMP_TASKSYNC_START DO_PRAGMA(omp taskgroup)
#define OMP_TASKSYNC_END
#elif _OPENMP >= 200805
#define OMP_TASKSYNC_START
#define OMP_TASKSYNC_END DO_PRAGMA(omp taskwait)
#endif
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
long thread_id(void);
/*
  For all the following calls, if urgent is 1, the job is queued in the front, otherwise in the end.

  CALL(fun,arg,nthread,urgent) executes "nthread" instances of fun with argument arg. 
  QUEUE is like CALL, but need an explicit WAIT
  
  The following use thread_t to manage the index
  QUEUE_THREAD(group, A, nthread, urgent) will queue arrays of thread_t (A) 
  CALL_THREAD(A, urgent) calls QUEUE_THREAD and waits for all to finish.
*/
#include <pthread.h>
#include "thread_pool.h"
#if _OPENMP >= 200805
/**
   Notice it is not effective to set the environment variables here.
*/
#define THREAD_YIELD	_Pragma("omp taskyield")
INLINE void THREAD_POOL_INIT(int nthread){
    fprintf(stderr, "Using OpenMP version %d with %d threads\n", _OPENMP, nthread);
    omp_set_num_threads(nthread);
    omp_set_nested(0);//make sure nested is not enabled
}
INLINE void CALL(thread_fun fun, void *arg, int nthread, int urgent){
    (void)urgent;
    OMP_TASKSYNC_START
    for(int it=0; it<nthread; it++){
#pragma omp task 
	fun(arg);
    }
    OMP_TASKSYNC_END
	}

/*The following QUEUE, CALL, WAIT acts on function (fun) and argument (arg).*/
/*Don't turn the following into INLINE function becase task will be waited*/
#define QUEUE(group,fun,arg,nthread,urgent)	\
    (void) group; (void) urgent;		\
    for(int it=0; it<nthread; it++){		\
	_Pragma("omp task")			\
	    fun(arg);				\
    }

#define WAIT(group)				\
    _Pragma("omp taskwait")


/*The following *_THREAD acts on thread_t array A*/

#define QUEUE_THREAD(group,A,urgent)					\
    (void)group;(void)urgent;						\
    if(!OMP_IN_PARALLEL){						\
	warning_once("QUEUE_THREAD is not in parallel region\n");	\
    }									\
    for(int it=0; it<(A)[0].nthread; it++){				\
	if((A)[it].fun){						\
	    _Pragma("omp task")						\
		(A)[it].fun((A)+it);					\
        }								\
    }
#define WAIT_THREAD(group)			\
    _Pragma("omp taskwait")
/*Turn to inline function because nvcc concatenates _Pragma to } */
INLINE void CALL_THREAD_DO(thread_t *A, int urgent){
    (void)urgent;
    OMP_TASKSYNC_START
	for(int it=0; it<A[0].nthread; it++){		
	    if(A[it].fun){
#pragma omp task
		A[it].fun(A+it); 
	    }
	}
    OMP_TASKSYNC_END
	}
INLINE void CALL_THREAD(thread_t *A, int urgent){
    //Split CALL_THREAD_DO to a separate routing to avoid recursiving calling CALL_THREAD. This is to work about a bug in icc that always return 0 for omp_in_parallel when nthread==1
    (void) urgent;
    if(!OMP_IN_PARALLEL){
	//Wraps call in parallel region.
	OMPTASK_SINGLE
	    CALL_THREAD_DO(A, urgent);
    }else{
	CALL_THREAD_DO(A, urgent);
    }
}

#else //using our thread_pool 

#define QUEUE(group,fun,arg,nthread,urgent)				\
    if(nthread>1){							\
	thread_pool_queue_many(&group, (thread_fun)fun, (void*)arg, nthread, urgent); \
    }else{								\
	fun(arg);							\
    }
INLINE void CALL(thread_fun fun, void *arg, int nthread, int urgent){
    if(nthread>1){							
	long thgroup=0;							
	thread_pool_queue_many						
	    (&thgroup, fun, arg, nthread, urgent);			
	thread_pool_wait(&thgroup);					
    }else{								
	fun(arg);							
    }
}
#define WAIT(group)				\
    thread_pool_wait(&group);
/**
   Queue jobs to group. Do not wait
*/
#define QUEUE_THREAD(group,A,urgent)					\
    if((A[0].nthread)>1){						\
	thread_pool_queue_many(&group,NULL,A,A[0].nthread,urgent);	\
    }else{								\
	(A)->fun(A);							\
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

#define THREAD_POOL_INIT(A) ({thread_pool_init(A);fprintf(stderr, "Using thread pool with %d threads\n", A);})
#define THREAD_YIELD thread_pool_do_job_once()
#endif

#define LOCK(A) pthread_mutex_lock(&A)
#define TRYLOCK(A) pthread_mutex_trylock(&A)
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

#define LOCKADD(dest,src,step)				\
    (({__asm__ __volatile__ ("lock; xaddl %0,%1"	\
		: "=r" (dest), "=m" (src)		\
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

//Using taskloop causes memory double free error in icc.
#if _OPENMP >= 201511 && 0//OpenMP 4.5 (taskloop introduced)
#define OMPTASK_FOR(index, start, end, extra...)	\
    DO_PRAGMA(omp taskloop extra)			\
    for(long index=start; index<end; index++)
#define OMPTASK_END
#elif _OPENMP >= 200805 //Openmp 3.0 (task introduced). Emulating taskloop.
#define OMPTASK_FOR(index, start, end, extra...)	\
    long omp_sect=(end-start+NTHREAD-1)/NTHREAD;	\
    OMP_TASKSYNC_START					\
    /*DO_PRAGMA(omp parallel for)*//*ignored if nested*/	\
    for(long omp_j=0; omp_j<NTHREAD; omp_j++){		\
    long omp_start=start+omp_sect*omp_j;		\
    long omp_end=omp_start+omp_sect;			\
    if(omp_end>end) omp_end=end;			\
    DO_PRAGMA(omp task extra if(omp_start<omp_end))	\
    for(long index=omp_start; index<omp_end; index++)
#define OMPTASK_END } OMP_TASKSYNC_END
#else
#define OMPTASK_FOR(index,start,end, extra...)	\
    for(long index=start; index<end; index++)
#define OMPTASK_END
#endif

//For those that is only good for icc, use the following
#if _OPENMP >= 200805 && defined(__INTEL_COMPILER) 
#define ICCTASK_FOR(index, start,end, extra...)	\
    OMPTASK_FOR(index,start,end,extra)
#define ICCTASK_END OMPTASK_END
#else
#define ICCTASK_FOR(index,start,end, extra...)	\
    for(long index=start; index<end; index++)
#define ICCTASK_END 
#endif




#endif //ifndef AOS_LIB_THREAD_H
