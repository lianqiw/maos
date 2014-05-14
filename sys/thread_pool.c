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
/**
   Contains implementation of a thread pool. It has the following capabilities:

   1) Maintains pending taskes in a linked list. Wake up idle threads to handle
   them.

   2) Maintains a pool of fixed number of threads. The threads takes jobs from
   the pool and run them. Enters into idle state if no jobs are pending.

   3) When a threads entering into a waiting state, it will go ahead to do other
   jobs and return when the jobs it is waiting for is finished.

   A thread pool is a pool of threads that can be used to run functions without
   spawning/ending a new thread. The threads are created with no job designated,
   but waiting for a condition.
*/
/*
  The jobs are keep in a linked list. thread_pool_queue will queue the jobs to
  the head or end of the list while each thread will remove jobs from the beginning of
  the list when they are waked.  

  2010-12-13: Improved the thread scheduling by test pool->nidle instead of
  pool->jobstail. This reduces the tomography of 100 steps down to 0.65s from
  0.8s for 8 threads at taurus.

  2010-12-16: Do not interrupt a thread and let it do other things. If the
  thread is caught in a mutex locked environment (like malloc/free), there will
  be a deadlock.

  2010-12-25: spin lock is twice as fast as mutex lock. Good for short term locks. 
  2011-01-07: spin lock is not available in MAC. reverted back to pthread.
*/

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include "common.h"
#include "misc.h"
#include "thread.h"
#include "thread_pool.h"
#include "process.h"
#define USE_SPIN_LOCK 1
#if USE_SPIN_LOCK == 1 && !defined(__INTEL_COMPILER)
/*Important to be volatile, otherwise lockes up in Sandy Bridge CPUs */
#define LOCK_T volatile int 
#define LOCK_DO(A) SPIN_LOCK(A)
#define LOCK_UN(A) SPIN_UNLOCK(A)
#define LOCK_INIT(A) A=0
#define LOCK_DESTROY(A)
#else
#define LOCK_T pthread_mutex_t
#define LOCK_DO(A) pthread_mutex_lock(&A)
#define LOCK_UN(A) pthread_mutex_unlock(&A)
#define LOCK_INIT(A) pthread_mutex_init(&A, NULL)
#define LOCK_DESTROY(A) pthread_mutex_destroy(&A)
#endif
/**
 * struct of jobs for the linked first in first out (fifo) list.
 */
typedef struct jobs_t{
    thread_fun fun;     /**<The function*/
    void *arg;          /**<The argument*/
    long *count;        /**<address of the job group, which is also a counter.*/
    int urgent;         /**<whether this job is urgent (can not be interruped)*/
    struct jobs_t *next;/**<The pointer to the next entry*/
}jobs_t;
/**
 *  The thread pool struct.
 */
struct thread_pool_t{
    pthread_mutex_t mutex; /**<the mutex.*/
    LOCK_T mutex_pool;         /**<the mutex for jobpool.*/
    pthread_cond_t jobwait;/**<there are jobs jobwait.*/
    pthread_cond_t idle;   /**<all threads are idle*/
    pthread_cond_t exited; /**<all threads have exited*/
    pthread_cond_t jobdone;/**<A job is done.*/
    pthread_attr_t attr;   /**<The attribution of newly created threads.*/
    jobs_t *jobshead;      /**<Start of the fifo list of jobwait jobs*/
    jobs_t *jobstail;      /**<End of the fifo list of jobwait jobs*/
    jobs_t *jobspool;      /**<Save the allocated but finished job data here*/
    int icur; /**<the top of the threads stack.*/
    int nmax; /**<the maximum number of threads. constant.*/
    int ncur; /**<the maximum number of live threads, excluding the master thread.*/
    int quit; /**<1: quit the threads.*/
    int nidle;/**<Number of idle threads, excluding the master thread.*/
    int inited;/**<Whether the thread pool has been inited*/
};
static thread_pool_t pool;/*The default pool; */

/**
 * Do a job. The mutex should have already been locked when calling this
 * routine. It will take a job from the head of the job queue, release the mutex,
 * run the job, and then acquire the mutex, and check whether the job is
 * finished.
 */
static inline void do_job(void) {
    /*Take the job out of the todo list */
    jobs_t *job=pool.jobshead;
    pool.jobshead=job->next;
    pthread_mutex_unlock(&pool.mutex);
    /*run the job */
    job->fun(job->arg); 
    pthread_mutex_lock(&pool.mutex);
    /*decrease the counter. */
    (*job->count)--;
    if(!(*job->count)){
	/*job is done. need broadcast since multiple threads may be waiting for
	  their job to finish*/
	pthread_cond_broadcast(&pool.jobdone);
    }
    /*return job to the pool. */
    LOCK_DO(pool.mutex_pool);
    job->next=pool.jobspool;
    pool.jobspool=job;
    LOCK_UN(pool.mutex_pool);
}
/**
   In the middle of an active thread, we can use do_urgent_job() periodically to
check for and switch to more urgent jobs.  */
void thread_pool_do_urgent_job(void){
    pthread_mutex_lock(&pool.mutex);
    while(pool.jobshead && pool.jobshead->urgent){
	do_job();
    }
    pthread_mutex_unlock(&pool.mutex);
}
/**
   In the middle of an thread that waits an event, we can use do_job_once() to
   do other jobs while waiting*/
int thread_pool_do_job_once(void){
    int ans=1;
    pthread_mutex_lock(&pool.mutex);
    if(pool.jobshead){
	do_job();
	ans=0;
    }
    pthread_mutex_unlock(&pool.mutex);
    return ans;
}
/**
 *   The working function in each thread which never quits unless pool is being
 *   destroyed. The master thread does not run this routine. The whole function
 *   is guarded by mutex. The mutex is only released upon 1) cond_wait, 2)
 *   do_job, 3) exit.
 */
static void run_thread(){
    pthread_mutex_lock(&pool.mutex);
    pool.ncur++;
    while(!pool.quit){
	/*
	  Wait for the go signal. The mutex is release during jobwait and
	  locked immediately when cond is satisfied.
	*/
	if(!pool.jobshead){
	    /*this thread is idle. add to the idle count */
	    pool.nidle++;
	    /*at most ncur-1 threads can be idle. */
	    if(pool.nidle+1==pool.ncur){/*all threads are idle */
		pthread_cond_signal(&pool.idle);
	    }
	    if(pool.ncur>pool.nmax){
		break;
	    }
	    /*no more jobs to do, wait for the condition. */
	    pthread_cond_wait(&pool.jobwait, &pool.mutex);
	    pool.nidle--;/*no longer idle. */
	    /*jobshead should be not empty now. */
	    if(!pool.jobshead){/*some other task already did the job. */
		continue;
	    }
	}
	/*Take the jobs out of the job queue and run it. */
	while(pool.jobshead){
	    do_job();
	}
    }
    pool.ncur--;
    if(pool.ncur==1){/*all thread are exited except the master thread. */
	pthread_cond_signal(&pool.exited);
    }
    pthread_mutex_unlock(&pool.mutex);
}

/**
 *   Initialize the thread pool.
 */
void thread_pool_init(int nthread){
    static int inited=0;
    if(!inited){
	inited=1;
	memset(&pool, 0, sizeof(thread_pool_t));
	pool.inited=1;
	pthread_mutex_init(&pool.mutex,  NULL);
	LOCK_INIT(pool.mutex_pool);
	pthread_cond_init(&pool.idle,    NULL);
	pthread_cond_init(&pool.jobwait, NULL);
	pthread_cond_init(&pool.jobdone, NULL);
	pthread_cond_init(&pool.exited,  NULL);
	pthread_attr_init(&pool.attr);
	pthread_attr_setdetachstate(&pool.attr,PTHREAD_CREATE_DETACHED);
	pool.nidle=0;
	pool.jobspool=NULL;
	pool.jobshead=NULL;
	pool.jobstail=NULL;
	pool.ncur=1;/*counting the master thread. */
	pool.nmax=1;
    }
    int nthread_max=sysconf( _SC_NPROCESSORS_ONLN );
    if(nthread<=0) nthread=nthread_max;
    if(nthread>nthread_max) nthread=nthread_max;
    pthread_mutex_lock(&pool.mutex);
    int ncur=pool.nmax;
    pool.nmax=nthread;
    if(ncur<nthread){/*need to launch more threads. */
	pthread_t thread;/*we don't need the thread information. */
	for(int ith=0; ith<nthread-ncur; ith++){
	    if(pthread_create(&thread, &pool.attr, (thread_fun)run_thread, NULL)){
		error("Can not create thread\n");
	    }
	}
    }else if(ncur>nthread){
	/*need to quit some threads. */
	pthread_cond_broadcast(&pool.jobwait);/*wake up all threads. */
    }
    pthread_mutex_unlock(&pool.mutex);
}
/**
 *   Queue a job that belongs to group denoted by group. The argument count,
 *  will be incremented by 1 when queued and decreased by 1 if job is
 *  finished. Wait on it will clear when count is decreased to zero.
 */
void thread_pool_queue(long *group, thread_fun fun, void *arg, int urgent){
    /*Add the job to the head if urgent>0, otherwise to the tail. */
    jobs_t *job;
    if(pool.jobspool){/*take it from the pool. */
	LOCK_DO(pool.mutex_pool);
	job=pool.jobspool;
	pool.jobspool=pool.jobspool->next;
	LOCK_UN(pool.mutex_pool);
    }else{
	job=malloc(sizeof(jobs_t));
    }
    job->fun=fun;
    job->arg=arg;
    job->count=group;   
    (*job->count)++;
    job->urgent=urgent;
    job->next=NULL;
    pthread_mutex_lock(&pool.mutex);
    if(pool.jobshead){/*list is not empty */
	if(urgent){
	    /*add to head */
	    job->next=pool.jobshead;
	    pool.jobshead=job;
	}else{     
	    /*add to tail. */
	    pool.jobstail->next=job;
	    pool.jobstail=job;
	}
    }else{
	/*list is empty. need to signal the thresds. */
	pool.jobshead=job;
	pool.jobstail=job;
    }
    if(pool.nidle>0){
	pthread_cond_signal(&pool.jobwait);/*wake up one thread only. */
	pthread_cond_signal(&pool.jobdone);/*wake up the thread that is waiting for jobdone. */
    }
    pthread_mutex_unlock(&pool.mutex);
}
/**
 *   Queue njob jobs. If fun is NULL, arg is thread_t array otherwise, arg is argument to fun directly.
 */
void thread_pool_queue_many(long *group, thread_fun fun, void *arg, int njob, int urgent){
    /*
      First create a temporary list of all the jobs without acquiring lock.
    */
    jobs_t *head=NULL;
    jobs_t *tail=NULL;
    int njob2=0;
    for(int ijob=0; ijob<njob; ijob++){
	jobs_t *job=NULL;
	LOCK_DO(pool.mutex_pool);
	if(pool.jobspool){
	    job=pool.jobspool;
	    pool.jobspool=pool.jobspool->next;
	}
	LOCK_UN(pool.mutex_pool);
	if(!job){
	    job=malloc(sizeof(jobs_t));
	}
	if(fun){
	    job->fun=fun;
	    job->arg=arg;
	}else{
	    thread_t *arg2=arg;
	    job->fun=(thread_fun)(arg2[ijob].fun);
	    job->arg=arg2+ijob;
	    if(((thread_t*)job->arg)->end==((thread_t*)job->arg)->start){
		/*empty job. don't queue it. */
		LOCK_DO(pool.mutex_pool);
		job->next=pool.jobspool;
		pool.jobspool=job;
		LOCK_UN(pool.mutex_pool);
		continue;
	    }
	}
	if(!tail) tail=job;
	job->count=group;
	job->urgent=urgent;
	job->next=head;
	head=job;
	njob2++;
    }
    /*Add the temporary list of jobs to the todo queue */
    pthread_mutex_lock(&pool.mutex);
    *group+=njob2;/*Need to be inside the lock. */
    if(pool.jobshead){/*list is not empty */
	if(urgent){/*add to head */
	    tail->next=pool.jobshead;
	    pool.jobshead=head;
	}else{/*add to tail. */
	    pool.jobstail->next=head;
	    pool.jobstail=tail;
	}
    }else{
	/*list is empty. need to signal the thresds. */
	pool.jobshead=head;
	pool.jobstail=tail;
    }
    if(pool.nidle){
	pthread_cond_broadcast(&pool.jobwait);/*wake up all idle threads. */
	pthread_cond_broadcast(&pool.jobdone);/*wake up the thread that is waiting for jobdone. */
    }
    pthread_mutex_unlock(&pool.mutex);
}
/**
 * Wait for jobs in the count to be done.
 */
void thread_pool_wait(long *count){
    pthread_mutex_lock(&pool.mutex);
    /*Do some job while we are waiting. */
    while((*count)){
	if(pool.jobshead){
	    do_job();
	}else{
	    pool.nidle++;
	    pthread_cond_wait(&pool.jobdone, &pool.mutex);
	    pool.nidle--;
	}
    }
    pthread_mutex_unlock(&pool.mutex);
}
/**
 * Wait for all jobs to be done.
 */
void thread_pool_wait_all(void){
    /*
      We should lock mutex before compare nidle, otherwise another thread maybe
      modifying nidle, and emits pool.idle before we are ready to wait, thus
      hanging us here.
    */
    pthread_mutex_lock(&pool.mutex);
    while(pool.jobshead){
	do_job();
    }
    if(pool.nidle+1<pool.ncur){/*some job is still doing the last bit. */
       	/*
	  thread_pool_wait may not obtain the lock in pool.mutex immediately
	  after signal on pool.idle is emmited. thread_pool_queue may obtain
	  the lock. So we wait in a loop.*/
	pthread_cond_wait(&pool.idle, &pool.mutex);
    }
    pthread_mutex_unlock(&pool.mutex);
}

/**
 *   Exit all threads and free thread pool.
 */
void thread_pool_destroy(void){
    if(!pool.inited) return;
    thread_pool_wait_all();/*let all jobs finish. */
    /*tell all jobs to quit. */
    pool.quit=1;
    pthread_cond_broadcast(&pool.jobwait);
    pthread_mutex_lock(&pool.mutex);
    if(pool.ncur>1){
	pthread_cond_wait(&pool.exited, &pool.mutex);
    }
    pthread_mutex_unlock(&pool.mutex);
    LOCK_DO(pool.mutex_pool);
    for(jobs_t *job=pool.jobspool; job; job=pool.jobspool){
	pool.jobspool=job->next;
	free(job);
    }
    LOCK_UN(pool.mutex_pool);
    pthread_mutex_destroy(&pool.mutex);
    LOCK_DESTROY(pool.mutex_pool);
    pthread_cond_destroy(&pool.idle);
    pthread_cond_destroy(&pool.jobwait);
    pthread_cond_destroy(&pool.jobdone);
    pthread_cond_destroy(&pool.exited);
    pthread_attr_destroy(&pool.attr);
    pool.inited=0;
}
static __attribute__((destructor)) void deinit(){
    thread_pool_destroy();
}

static __attribute__((constructor)) void init(){
    register_deinit(thread_pool_destroy,NULL);/*register to mem.c */
}

