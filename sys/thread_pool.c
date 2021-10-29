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
/**
   A thread pool is a pool of threads that can be used to run functions without
   spawning/ending a new thread. The threads are created with no job designated,
   but waiting for a condition and wake up to do jobs.

   Contains implementation of a thread pool. It has the following capabilities:

   1) Maintains pending taskes in a linked list. Wake up idle threads to handle
   them.

   2) Maintains a pool of fixed number of threads. The threads takes jobs from
   the pool and run them. Enters into idle state if no jobs are pending.

   3) When a threads entering into a waiting state, it will go ahead to do other
   jobs and return when the jobs it is waiting for is finished.

*/
/*
  The jobs are keep in a linked list. thread_pool_queue will queue the jobs to
  the head or end of the list while each thread will remove jobs from the
  beginning of the list when they are waked.

  2010-12-13: Improved the thread scheduling by test pool->nidle instead of
  pool->jobstail. This reduces the tomography of 100 steps down to 0.65s from
  0.8s for 8 threads at taurus.

  2010-12-16: Do not interrupt a thread and let it do other things. If the
  thread is caught in a mutex locked environment (like malloc/free), there will
  be a deadlock.

  2010-12-25: spin lock is twice as fast as mutex lock. Good for short term
  locks. 2011-01-07: spin lock is not available in MAC. reverted back to
  pthread.

  2021-10-29: replaced mutex lock by lock free operation. It is hard to maintain
  atomically the list from both head and tail. So we split that into two
  separate queues.
*/


#include <unistd.h>
#include <pthread.h>

#include "common.h"
#include "misc.h"
#include "thread.h"
#include "thread_pool.h"
#include "process.h"

/**
 * struct of jobs for the linked first in first out (fifo) list.
 */
typedef struct jobs_t{
	thread_wrapfun fun;     /**<The function*/
	void* arg;          /**<The argument*/
	long* count;        /**<address of the job group, which is also a counter.*/
	struct jobs_t* next;/**<The pointer to the next entry*/
}jobs_t;
/**
 *  The thread pool struct.
 */
struct thread_pool_t{
	pthread_mutex_t mutex; /**<mutex to protect jobshead and jobstail.*/
	pthread_cond_t jobwait;/**<there are jobs jobwait.*/
	pthread_cond_t idle;   /**<all threads are idle*/
	pthread_cond_t exited; /**<all threads have exited*/
	pthread_attr_t attr;   /**<The attribution of newly created threads.*/
	
	int icur; /**<the top of the threads stack.*/
	int nmax; /**<the maximum number of threads. constant.*/
	int ncur; /**<the maximum number of live threads, excluding the master thread.*/
	int quit; /**<1: quit the threads.*/
	int nidle;/**<Number of idle threads, excluding the master thread.*/
	int inited;/**<Whether the thread pool has been inited*/
};

static thread_pool_t pool;/*The default pool; */
static jobs_t *jobsnormal=NULL;/**<Start of the fifo list of pending jobs*/
static jobs_t *jobsurgent=NULL;/**<Start of the fifo list of urgent jobs*/
static jobs_t* jobspool=NULL;  /**<saving unused jobs_t*/


//For a single job, head==tail
static void jobs_put(jobs_t **list, jobs_t *job){
	job->next=*list;
	while(!atomic_compare_exchange(list, &job->next, job));
	
}
static jobs_t* jobs_get(jobs_t **list){
	jobs_t* job=*list;
	while(job && !atomic_compare_exchange(list, &job, job->next));
	return job;
}

static void jobspool_free(){
	jobs_t *job=0;
	while((job=jobs_get(&jobspool))){
		free(job);
	}
}
/**
 * Do a job. The mutex should have already been locked when calling this
 * routine. It will take a job from the head of the job queue, release the mutex,
 * run the job, and then acquire the mutex, and check whether the job is
 * finished.
 */
static inline int do_job(void){
	/*Take the job out of the todo list */
	jobs_t *job=NULL;
	if((job=jobs_get(&jobsurgent)) || (job=jobs_get(&jobsnormal))){
		/*run the job */
		job->fun(job->arg);
		atomic_add_fetch(job->count, -1);
		/*return job to the pool instead of freeing, using lock free approach. */
		jobs_put(&jobspool, job);
		return 1;//more jobs pending
	}
	return 0;//no more job
}
/**
 *   The working function in each thread which never quits until pool is being
 *   resized or destroyed. The master thread does not run this routine. 
 */
static void* run_thread(void* data){
	(void)data;
	atomic_add_fetch(&pool.ncur, 1);
	while(!pool.quit){
		while(do_job());//do until no jobs are left
		pthread_mutex_lock(&pool.mutex);//acquire lock before modifying pool and wait. 	
		pool.nidle++;/*add to the idle thread count */
		if(pool.nidle+1==pool.ncur){/*at most ncur-1 threads can be idle as the main thread counts as 1. */
			pthread_cond_signal(&pool.idle);//all threads are idle 
		}
		if(pool.ncur>pool.nmax){//quit excessive threads
			break;
		}
		pthread_cond_wait(&pool.jobwait, &pool.mutex);//lock is released during wait and reacquired after wake
		pool.nidle--;/*no longer idle. */
		pthread_mutex_unlock(&pool.mutex);//release the lock
	}
	//mutex is locked
	pool.ncur--;
	if(pool.ncur==1){/*all thread are exited except the master thread. */
		pthread_cond_signal(&pool.exited);
	}
	pthread_mutex_unlock(&pool.mutex);
	return NULL;
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
		pthread_mutex_init(&pool.mutex, NULL);
		pthread_cond_init(&pool.idle, NULL);
		pthread_cond_init(&pool.jobwait, NULL);
		pthread_cond_init(&pool.exited, NULL);
		pthread_attr_init(&pool.attr);
		pthread_attr_setdetachstate(&pool.attr, PTHREAD_CREATE_DETACHED);
		pool.nidle=0;
		pool.ncur=1;/*counting the master thread. */
		pool.nmax=1;
	}
	int nthread_max=sysconf(_SC_NPROCESSORS_ONLN);
	if(nthread<=0) nthread=nthread_max;
	if(nthread>nthread_max) nthread=nthread_max;
	pthread_mutex_lock(&pool.mutex);
	int ncur=pool.nmax;
	pool.nmax=nthread;
	if(ncur<nthread){/*need to launch more threads. */
		pthread_t thread;/*we don't need the thread information. */
		for(int ith=0; ith<nthread-ncur; ith++){
			if(pthread_create(&thread, &pool.attr, run_thread, NULL)){
				error("Can not create thread\n");
			}
		}
	} else if(ncur>nthread){/*need to quit some threads. */
		pthread_cond_broadcast(&pool.jobwait);/*wake up all threads. */
	}
	pthread_mutex_unlock(&pool.mutex);
}
/**
 *  Queue a job that belongs to group denoted by group. The argument count,
 *  will be incremented by 1 when queued and decreased by 1 if job is
 *  finished. Wait on it will clear when count is decreased to zero.
 * 
 *  Job is added to jobsnormal or joburgent depending on urgent value
 */
void thread_pool_queue(long* group, thread_wrapfun fun, void* arg, int urgent){
	jobs_t *job=jobs_get(&jobspool);/*take idle job from the pool using lock free approach. */
	if(!job){
		job=(jobs_t*)malloc(sizeof(jobs_t));
	}
	job->fun=fun;
	job->arg=arg;
	job->count=group;
	atomic_add_fetch(group, 1);
	job->next=NULL;
	jobs_put(urgent?&jobsurgent:&jobsnormal, job);

	if(pool.nidle){//no need to hold mutex to call cond_signal
		pthread_cond_signal(&pool.jobwait);/*wake up one thread only. */
	}
}
/**
 *   Queue njob jobs. If fun is NULL, arg is thread_t array otherwise, arg is argument to fun directly.
 */
void thread_pool_queue_many(long* group, thread_wrapfun fun, void* arg, int njob, int urgent){
	for(int ijob=0; ijob<njob; ijob++){
		if(fun){
			thread_pool_queue(group, fun, arg, urgent);
		}else{
			thread_t *arg2=(thread_t *)arg;
			if(arg2[ijob].fun && (arg2[ijob].start < arg2[ijob].end)){
				thread_pool_queue(group, arg2[ijob].fun, &arg2[ijob], urgent);
			}
		}
	}
}
/**
 * Wait for jobs in the count to be done.
 */
void thread_pool_wait(long* count){
	while(atomic_load(count)){
		do_job();/*Do some job while we are waiting. */
	}
}
/**
 * Wait for all jobs to be done.
 */
void thread_pool_wait_all(void){
	/*
	  We should lock mutex before compare nidle, otherwise another thread maybe
	  modifying nidle, and emits pool.idle before we are ready to wait, thus
	  hanging us here.

	  thread_pool_wait may not obtain the lock in pool.mutex immediately after
	  signal on pool.idle is emmited. thread_pool_queue may obtain the lock.
	  So we wait in a loop.
	*/
	while(do_job());
	pthread_mutex_lock(&pool.mutex);
	if(pool.nidle+1<pool.ncur){/*some job is still doing the last bit. */
		pthread_cond_wait(&pool.idle, &pool.mutex);
	}
	pthread_mutex_unlock(&pool.mutex);
}

/**
 *   Exit all threads and free thread pool.
 */
void thread_pool_destroy(void){
	//only 1 thread can wait
	if(atomic_add_fetch(&pool.inited, -1)!=0) return;
	thread_pool_wait_all();/*let all jobs finish. */
	/*tell all jobs to quit. */
	pool.quit=1;
	pthread_cond_broadcast(&pool.jobwait);
	pthread_mutex_lock(&pool.mutex);
	if(pool.ncur>1){
		pthread_cond_wait(&pool.exited, &pool.mutex);
	}
	pthread_mutex_unlock(&pool.mutex);
	pthread_mutex_destroy(&pool.mutex);
	pthread_cond_destroy(&pool.idle);
	pthread_cond_destroy(&pool.jobwait);
	pthread_cond_destroy(&pool.exited);
	pthread_attr_destroy(&pool.attr);
	jobspool_free();
}
