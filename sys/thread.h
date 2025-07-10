/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
//#include <stdatomic.h> //_Atomic with OPENMP is not supported in compilares <=6. We use builtin atomic functions instead.
#include <pthread.h>
#include "common.h"

typedef struct thread_t thread_t;
typedef void *(*thread_fun)(void *);
typedef void *(*thread_wrapfun)(thread_t *);
#include "thread_pool.h"
#define DO_PRAGMA(A...) _Pragma(#A)
#define PRINT_MACRO(x) DO_PRAGMA(message #x " is " STR(x))
#define STR(x) #x
//The following are common to OPENMP or PTHREADS 

#define LOCK(A) pthread_mutex_lock(&A)
#define TRYLOCK(A) pthread_mutex_trylock(&A)
#define UNLOCK(A) pthread_mutex_unlock(&A)
#define PINIT(A) pthread_mutex_init(&A,NULL)
#define PDEINIT(A) pthread_mutex_destroy(&A)
#define PNEW(A) static pthread_mutex_t A=PTHREAD_MUTEX_INITIALIZER
#define PNEW2(A) pthread_mutex_t A=PTHREAD_MUTEX_INITIALIZER
static inline void block_signal(int block){
	sigset_t signal_mask;
	sigemptyset (&signal_mask);
    sigaddset (&signal_mask, SIGINT);
    sigaddset (&signal_mask, SIGTERM);
	
    if(pthread_sigmask (block?SIG_BLOCK:SIG_UNBLOCK, &signal_mask, NULL)){
		warning("masking/unmasking signal failed\n");
	}
}
extern pthread_mutex_t mutex_fftw;
#define LOCK_FFT LOCK(mutex_fftw)
#define UNLOCK_FFT UNLOCK(mutex_fftw)

/**
   Information about job to launch for each thread. start and end are the two indices.
*/

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
///prepare thread information
thread_t *thread_prep(long start, long end, long nthread, thread_wrapfun fun, void *data);
///Create a thread and take it away
pthread_t thread_new(thread_fun fun, void *arg);
void thread_block_signal();
//Becareful about memory order. It determines the optimization constraints of other operations around the atomic operation
//For example, if the store before atomic store must be visible after atomic load, memory acquire/release or sequential has to be used
#define MEM_ORDER __ATOMIC_SEQ_CST
//The build in functions assume pointer is of type unsigned int
//We use functions instead of defines to enforce type
static inline unsigned int atomic_sub_fetch(unsigned int *ptr, unsigned int val){
    return __atomic_sub_fetch(ptr, val, MEM_ORDER);
}
static inline unsigned int atomic_add_fetch(unsigned int *ptr, unsigned int val){
    return __atomic_add_fetch(ptr, val, MEM_ORDER);
}
static inline unsigned int atomic_fetch_add(unsigned int *ptr, unsigned int val){
    return __atomic_fetch_add(ptr, val, MEM_ORDER);
}
static inline unsigned int atomic_fetch_sub(unsigned int *ptr, unsigned int val){
    return __atomic_fetch_sub(ptr, val, MEM_ORDER);
}
//weak is ok since we do while
#define atomic_compare_exchange_n(ptr, pexpected, desired) __atomic_compare_exchange_n(ptr, pexpected, desired, 0, MEM_ORDER, MEM_ORDER)
#define atomic_compare_exchange(ptr, pexpected, pdesired) __atomic_compare_exchange(ptr, pexpected, pdesired, 0, MEM_ORDER, MEM_ORDER) 
#define atomic_load(ptr) __atomic_load_n(ptr, MEM_ORDER) 
#define atomic_store(ptr, val) __atomic_store_n(ptr, val, MEM_ORDER) 

#if _OPENMP >= 201511
#define OMPTASK_SINGLE	\
    DO_PRAGMA(omp parallel)\
    DO_PRAGMA(omp single)			
#define OMP_IN_PARALLEL omp_in_parallel()
#if DEBUG
#define expect_level(n) if(omp_get_level()!=n) {dbg_once("omp_get_level=%d, want %d, omp_get_active_level=%d, omp_in_parallel=%d\n", omp_get_level(), n, omp_get_active_level(), omp_in_parallel());}
#else
#define expect_level(n)
#endif
#define OMP_FOR(nthread)    expect_level(0);DO_PRAGMA(omp parallel for default(shared) num_threads(nthread))
#define OMP_FOR_COLLAPSE(n, nthread) expect_level(0);DO_PRAGMA(omp parallel for default(shared) collapse(n) num_threads(nthread))
#define OMP(A) DO_PRAGMA(omp A)
#define OMP_TASK_FOR(ntask)      expect_level(1);DO_PRAGMA(omp taskloop default(shared) num_tasks(ntask) priority(1))
#define OMP_TASK_FOR_COLLAPSE(n, ntask) expect_level(1);DO_PRAGMA(omp taskloop default(shared) num_tasks(ntask) collapse(n) priority(1))
#else
#define OMPTASK_SINGLE
#define OMP_IN_PARALLEL 0
#define OMP_FOR(n)
#define OMP_FOR_COLLAPSE(n, nthread)
#define OMP(A)
#define OMP_TASK_FOR(ntask)
#define OMP_TASK_FOR_COLLAPSE(n, ntask)
#endif

/*
  For all the following calls, if urgent is 1, the job is queued in the front, otherwise in the end.

  CALL(fun,arg,nthread,urgent) executes "nthread" instances of fun with argument arg.
  QUEUE is like CALL, but need an explicit WAIT

  The following use thread_t to manage the index
  QUEUE_THREAD(counter, A, nthread, urgent) will queue arrays of thread_t (A)
  CALL_THREAD(A, urgent) calls QUEUE_THREAD and waits for all to finish.
*/

#if _OPENMP >= 201511 //minimum openmp version supported for taskloop
/**
   Notice it is not effective to set the environment variables here.
*/

static inline void THREAD_POOL_INIT(int nthread){
    info("Using OpenMP version %d with %d threads\n", _OPENMP, nthread);
    omp_set_num_threads(nthread);
}
//static inline void QUEUE(tp_counter_t *counter, thread_wrapfun fun, void *arg, int nthread, int urgent){
#define QUEUE(counter, fun, arg, nthread, urgent) \
    (void)counter; \
if(omp_in_parallel()){\
    DO_PRAGMA(omp taskloop default(shared) priority(urgent) nogroup)\
    for(int it=0; it<nthread; it++){\
        fun(arg);\
    }\
}else{\
    DO_PRAGMA(omp parallel for default(shared))\
    for(int it=0; it<nthread; it++){\
        fun(arg);\
    }\
}

//static inline void CALL(thread_wrapfun fun, void *arg, int nthread, int urgent){
#define CALL(fun, arg, nthread, urgent) \
if(omp_in_parallel()){\
    DO_PRAGMA(omp taskloop default(shared) priority(urgent))\
    for(int it=0; it<nthread; it++){\
        fun(arg);\
    }\
}else{\
    DO_PRAGMA(omp parallel for default(shared))\
    for(int it=0; it<nthread; it++){\
        fun(arg);\
    }\
}

#define WAIT(pcounter, urgent) DO_PRAGMA(omp taskwait)
/*Turn to static inline function because nvcc concatenates _Pragma to } */
//Define causes _Pragma to appear in sameline in gcc4.9

#define QUEUE_THREAD(pcounter, A, urgent) \
if(omp_in_parallel()){\
    DO_PRAGMA(omp taskloop default(shared) priority(urgent) nogroup)\
    for(long it=0; it<(A)->nthread; it++){\
        if((A)[it].fun) (A)[it].fun((&(A)[it]));\
    }\
}else{\
	DO_PRAGMA(omp parallel for default(shared))\
	for(long it=0; it<(A)->nthread; it++){\
		if((A)[it].fun) (A)[it].fun((&(A)[it]));\
	}\
}

#define CALL_THREAD(A, urgent)\
if(omp_in_parallel()){\
    DO_PRAGMA(omp taskloop default(shared) priority(urgent))\
    for(long it=0; it<(A)->nthread; it++){\
        if((A)[it].fun) (A)[it].fun((A)+it);\
    }\
}else{\
	DO_PRAGMA(omp parallel for default(shared))\
    for(long it=0; it<(A)->nthread; it++){\
        if((A)[it].fun) (A)[it].fun((A)+it);\
    }\
}
#define CALL_THREAD_ARR(A, nA, urgent)\
if(omp_in_parallel()){\
    DO_PRAGMA(omp taskloop collapse(2) default(shared) priority(urgent))\
	for(long iA=0; iA<nA; iA++){\
    	for(long it=0; it<A[0]->nthread; it++){\
        	if(A[iA][it].fun) (A[iA])[it].fun(A[iA]+it);\
    	}\
	}\
}else{\
	DO_PRAGMA(omp parallel for collapse(2) default(shared))\
    for(long iA=0; iA<nA; iA++){\
    	for(long it=0; it<A[0]->nthread; it++){\
        	if(A[iA][it].fun) (A[iA])[it].fun(A[iA]+it);\
    	}\
	}\
}

#else //using our thread_pool

#define THREAD_POOL_INIT(A) ({thread_pool_init(A);info("Using thread pool with %d threads\n", A);})

/// Queue jobs to counter. Do not wait
#define QUEUE(counter, fun, arg, njob, urgent) thread_pool_queue(counter, (thread_wrapfun)fun, arg, njob, urgent)
#define QUEUE_THREAD(counter, A, urgent) thread_pool_queue(counter, NULL, A, A[0].nthread, urgent)

//We use define so that we can print the calling line 
#if ENABLE_TP_TIMING //defined in thread_pool.h
#define WAIT(counter, urgent) ({thread_pool_wait(counter, urgent);\
dbg("tmin=%u, tmax=%u ms\n", (counter)->tmin, (counter)->tmax);(counter)->tmin=0;(counter)->tmax=0;})
#else
#define WAIT(counter, urgent) thread_pool_wait(counter, urgent)
#endif

///Queue jobs to a temp counter, Then wait for it to complete.
#define CALL(fun, arg, nthread, urgent) \
({\
    if(nthread>1){\
        tp_counter_t counter={0};\
        QUEUE(&counter, fun, arg, nthread, urgent);\
        WAIT(&counter, urgent);\
    } else{\
        fun((void *)(arg));\
    }\
})\

#define CALL_THREAD(A, urgent)\
    if((A[0].nthread)>1){\
        tp_counter_t counter={0};\
        QUEUE_THREAD(&counter, A, urgent);\
        WAIT(&counter, urgent);\
    } else{\
        (A)->fun(A);\
    }\

#define CALL_THREAD_ARR(A, nA, urgent)\
	for(long iA=0; iA<nA; iA++){\
		if((A[iA][0].nthread)>1){\
			tp_counter_t counter={0};\
			QUEUE_THREAD(&counter, A[iA], urgent);\
			WAIT(&counter, urgent);\
		} else{\
			(A[iA])->fun(A[iA]);\
		}\
	}
#endif //#if _OPENMP
#endif //ifndef AOS_LIB_THREAD_H
