/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   4.0 (201307): Introduced taskgroup.
   4.5 (201511): Introduced taskloop and priority. Taskloop has implicit taskgroup.
*/

#include "common.h"
#define DO_PRAGMA(A...) _Pragma(#A)

#if HAS_OPENMP && !defined(_OPENMP)
#define _OPENMP 200805 //Sometimes GPU code does not define this
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
    void *data;/*shared data between threads*/
    void *thread_data; /*private data for each thread*/
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
#define THREAD_YIELD	DO_PRAGMA(omp taskyield)
#if _OPENMP>=201511 && !defined(__INTEL_COMPILER)  //ICC-17 does not support task priority.
#define OMP_TASK(urgent) DO_PRAGMA(omp task priority(urgent))
#else
#define OMP_TASK(urgent) DO_PRAGMA(omp task)
#endif
static inline void THREAD_POOL_INIT(int nthread){
    info("Using OpenMP version %d with %d threads\n", _OPENMP, nthread);
    omp_set_num_threads(nthread);
    omp_set_nested(0);//make sure nested is not enabled
}
static inline void QUEUE(long *group, thread_wrapfun fun, void *arg, int nthread, int urgent){
    (void) group;
    (void) urgent;
    for(int it=0; it<nthread; it++){
	OMP_TASK(urgent)
	    fun((thread_t*)arg);
    }
}
static inline void CALL(thread_wrapfun fun, void *arg, int nthread, int urgent){
    (void)urgent;
    OMP_TASKSYNC_START
	QUEUE(NULL, fun, arg, nthread, urgent);
    OMP_TASKSYNC_END
}

/*The following QUEUE, CALL, WAIT acts on function (fun) and argument (arg).*/
/*Don't turn the following into static inline function becase task will be waited*/
/*
#define QUEUE(group,fun,arg,nthread,urgent)	\
    (void) group; (void) urgent;		\
    for(int it=0; it<nthread; it++){		\
	OMP_TASK(urgent)			\
	    fun(arg);				\
    }
*/
#define WAIT(group) DO_PRAGMA(omp taskwait)

/*Turn to static inline function because nvcc concatenates _Pragma to } */
static inline void QUEUE_THREAD(long *group, thread_t *A, int urgent){
    (void)urgent;
    (void)group;
    for(int it=0; it<A[0].nthread; it++){		
	if(A[it].fun){
	    OMP_TASK(urgent)
		A[it].fun((A+it)); 
	}
    }
}

static inline void CALL_THREAD(thread_t *A, int urgent){
    /*Split CALL_THREAD_DO to a separate routing to avoid recursiving calling
     * CALL_THREAD. This is to work about a bug in icc that always return 0 for
     * omp_in_parallel when nthread==1*/
    (void) urgent;
    if(!OMP_IN_PARALLEL){
	//Wraps call in parallel region. No need to sync here.
	OMPTASK_SINGLE
	    QUEUE_THREAD(NULL, A, urgent);
    }else{
	OMP_TASKSYNC_START
	    QUEUE_THREAD(NULL, A, urgent);
	OMP_TASKSYNC_END
    }
}

#else //using our thread_pool


/**
   Queue jobs to group. Do not wait
*/
#define QUEUE thread_pool_queue_many
static inline void  QUEUE_THREAD(long *group, thread_t *A, int urgent){
    thread_pool_queue_many(group,NULL,A,A[0].nthread,urgent);
}
#define WAIT(group) thread_pool_wait(&group);
/**
   Queue jobs to a temp group, Then wait for it to complete.
*/
static inline void CALL(thread_wrapfun fun, void *arg, int nthread, int urgent){
    if(nthread>1){							
	long group=0; 
	QUEUE(&group, fun, arg, nthread, urgent); 
	WAIT(group); 
    }else{								
	fun((thread_t*)arg);							
    }
}


static inline void  CALL_THREAD(thread_t *A, int urgent){
    if((A[0].nthread)>1){ 
	long group=0; 
	QUEUE_THREAD(&group,A,urgent);
	WAIT(group); 
    }else{ 
	(A)->fun(A); 
    }
}

#define THREAD_POOL_INIT(A) ({thread_pool_init(A);info("Using thread pool with %d threads\n", A);})
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
   Create a new thread and forget.
*/
int thread_new(thread_fun fun, void* arg);
void thread_block_signal();


static inline int cmpxchg(int *ptr, int old, int newval){
    volatile int *__ptr = (volatile int *)(ptr);	
    int __ret;                                     
    __asm__ volatile( "lock; cmpxchg %2,%1"
		      : "=a" (__ret), "+m" (*__ptr) 
		      : "r" (newval), "0" (old)                     
		      : "memory");				 
    return __ret;
}

static inline int atomicadd(int *ptr, int val){
    int old;
    do{
	old=*ptr;
    }while(cmpxchg(ptr, old, old+val)!=old);
    return old+val;
}

//Using taskloop causes memory double free error in icc.
#if _OPENMP >= 201511  && !defined(__INTEL_COMPILER)  //ICC V17 taskloop is not correct.
#define OMPTASK_FOR(index, start, end, extra...)	\
    DO_PRAGMA(omp taskloop extra grainsize(50))		\
    for(long index=start; index<end; index++)
#define OMPTASK_END 
#elif _OPENMP >= 200805 //Openmp 3.0 (task introduced). Emulating taskloop.
#define OMPTASK_FOR(index, start, end, extra...)	\
    const long nratio=(end-start+NTHREAD-1)/NTHREAD;	\
    const long nsect=(nratio+4)/5;			\
    const long omp_sect=(end-start+nsect-1)/nsect;	\
    OMP_TASKSYNC_START					\
    for(long omp_j=0; omp_j<nsect; omp_j++){		\
        long omp_start=start+omp_sect*omp_j;		\
        long omp_end=omp_start+omp_sect;		\
	if(omp_end>end) omp_end=end;			\
	DO_PRAGMA(omp task extra firstprivate(omp_start, omp_end))	\
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
