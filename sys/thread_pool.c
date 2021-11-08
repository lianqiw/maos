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

  2021-11-01: The lock free approach suffers ABA problem that causes some job to
  be randomly lose. Replaced pointer based stack to integer based. Use
  jobsall[jobind].next instead of job->next. This solves the ABA problem by
  having a pop counter that is also an int. We didn't use pointer + counter as
  16 bit CAS is not always available.
*/

#include <unistd.h>
#include <pthread.h>
#include <limits.h>
#include "common.h"
#include "misc.h" //mysleep
#include "thread.h"
#include "thread_pool.h"
//#include "process.h"
_Thread_local unsigned int tid=0;//thread id

/**
 *  The thread pool struct.
 */
typedef struct thread_pool_t{
	pthread_mutex_t mutex; /**<mutex to protect jobshead and jobstail.*/
	pthread_cond_t jobwait;/**<there are jobs jobwait.*/
	pthread_cond_t idle;   /**<all threads are idle*/
	pthread_cond_t exited; /**<all threads have exited*/
	pthread_attr_t attr;   /**<The attribution of newly created threads.*/

	unsigned int icur; /**<the top of the threads stack.*/
	unsigned int nmax; /**<the maximum number of threads. constant.*/
	unsigned int ncur; /**<the maximum number of live threads, excluding the master thread.*/
	unsigned int quit; /**<1: quit the threads.*/
	unsigned int nidle;/**<Number of idle threads, excluding the master thread.*/
	unsigned int inited;/**<Whether the thread pool has been inited*/
}thread_pool_t;
static thread_pool_t pool;/*The default pool; */

/**
 * struct of jobs for the linked first in first out (fifo) list.
 * We use index here to emulate a 4 byte pointer
 */
typedef struct jobs_t{
	thread_wrapfun fun; /**<The function*/
	void* arg;          /**<The argument*/
	tp_counter_t *counter;/**<address of the job group, whos value keeps track of queued jobs.*/
	unsigned int next;  /**<4 byte index into jobsall for next job, replacing pointer*/
}jobs_t;

static jobs_t *jobsall=NULL;//all jobs are allocated here. The first element has next set to 0 to indicate end of list
static unsigned int jobsall_i=0; //available job to use
static unsigned int jobsall_count=0;//number of jobs in jobsall

/**
 * Counter and index of job.
 * */
typedef union{
	unsigned long state;//for CAS
	struct{
		unsigned int counter;//pop counter
		unsigned int head;   //index to jobsall of the first node
	};
}jobshead_t;

static jobshead_t jobsnormal={0L};/**<Start of the fifo list of pending jobs*/
static jobshead_t jobsurgent={0L};/**<Start of the fifo list of urgent jobs*/
static jobshead_t jobspool={0L};  /**<saving unused jobs_t*/


//Place job to the beginning of list.
static void jobs_push(jobshead_t *head, unsigned int jobheadind, unsigned int jobtailind){
	jobsall[jobtailind].next=head->head;
	while(!atomic_compare_exchange_n(&(head->head), &jobsall[jobtailind].next, jobheadind));
}
//Get job from the beginning of list
static unsigned int jobs_pop(jobshead_t *head){
	jobshead_t job, job2;
	job.counter=atomic_add_fetch(&(head->counter), 1);
	job.head=head->head;//get both values
	do{
		//compare both counter and pointer to make sure the node is not changed.
		job2.counter=job.counter;//increse the counter
		job2.head=jobsall[job.head].next;
	}while(!atomic_compare_exchange(head, &job, &job2));
	return job.head;
}

/**
 * Do a job. If urgent is set, only do urgent jobs.
 */
static inline int do_job(int urgent){
	unsigned int jobind=0;
	if((jobind=jobs_pop(&jobsurgent)) || (!urgent && (jobind=jobs_pop(&jobsnormal)))){
		jobs_t *job=&jobsall[jobind];
#if ENABLE_TP_TIMING
		double tk0=myclockd();
#endif		
		job->fun(job->arg);//run the job
#if ENABLE_TP_TIMING
		{
			unsigned int tk1=(unsigned int)ceil((myclockd()-tk0)*1e3);
			unsigned int old_value=job->counter->tmax;
			while(old_value<tk1 && !atomic_compare_exchange_n(&job->counter->tmax, &old_value, tk1));
			old_value=job->counter->tmin;
			unsigned int zero=0;
			if(!atomic_compare_exchange_n(&job->counter->tmin, &zero, tk1)){//not zero
				while(old_value>tk1&&!atomic_compare_exchange_n(&job->counter->tmin, &old_value, tk1));
			}
		}
#endif		
		atomic_sub_fetch(&job->counter->group, 1);
		jobs_push(&jobspool, jobind, jobind);//return job to the pool.
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
	tid=atomic_add_fetch(&pool.ncur, 1);//increment pool counter
	while(1){
		while(do_job(0));//do until no jobs are left
		pthread_mutex_lock(&pool.mutex);//acquire lock before modifying pool and wait. 	
		pool.nidle++;//increment idle counter
		if(pool.nidle+1==pool.ncur){//at most ncur-1 threads can be idle as the main thread counts as 1.
			pthread_cond_signal(&pool.idle);//all threads are idle 
		}
		pthread_cond_wait(&pool.jobwait, &pool.mutex);//lock is released during wait and reacquired after wake
		pool.nidle--;//no longer idle. 
		if(pool.quit||pool.ncur>pool.nmax){//quit this thread
			pool.ncur--;
			if(pool.ncur==1){//all thread are exited except the master thread. 
				pthread_cond_signal(&pool.exited);
			}
			pthread_mutex_unlock(&pool.mutex);
			break;
		}
		pthread_mutex_unlock(&pool.mutex);//release the lock
	}
	return NULL;
}

/**
 *  Queue jobs that belongs to group denoted by group. The argument group, will
 *  be incremented by 1 for each job queued and decreased by 1 after each job is
 *  finished. Wait on it will clear when count is decreased to zero.
 *
 *  Job is added to jobsnormal or joburgent depending on urgent value
 */
void thread_pool_queue(tp_counter_t *counter, thread_wrapfun fun, void* arg, int njob, int urgent){
	unsigned int headind=0, tailind=0;
	thread_t *arg2=(thread_t *)arg;
	for(int ijob=0; ijob<njob; ijob++){
		unsigned int jobind;
		jobs_t *job=NULL;
	retry:
		jobind=jobs_pop(&jobspool);//take idle job from the pool using lock free approach. 
		if(jobind==0){//no more jobs in the pool, allocate from jobsall
			jobind=atomic_fetch_add(&jobsall_i, 1);
			//need to check >= not = as jobs_pop may fail multiple times.
			if(jobind>=jobsall_count){//all jobs are being used
				//do not try to realloc jobsall as it is used by running threads.
				if(headind){//queue jobs already processed
					jobs_push(urgent?&jobsurgent:&jobsnormal, headind, tailind);
					headind=0;
					tailind=0;
				}
				if(pool.nidle){//wake up idle threads. no need to hold mutex to call cond_signal
					pthread_cond_broadcast(&pool.jobwait);
				}
				while(do_job(urgent));//help clear jobs
				goto retry;//retry
			}
		}
		job=&jobsall[jobind];
		if(fun){
			job->fun=fun;
			job->arg=arg;
		}else{
			if(!arg2[ijob].fun || arg2[ijob].start>=arg2[ijob].end){//invalid job
				jobs_push(&jobspool, jobind, jobind);
				continue;
			}
			job->fun=arg2[ijob].fun;
			job->arg=&arg2[ijob];
		}
		job->counter=counter;
		atomic_add_fetch(&counter->group, 1);
		job->next=0;
		if(!headind){
			headind=jobind;
		}else{
			jobsall[tailind].next=jobind;
		}
		tailind=jobind;
	}
	if(headind){
		jobs_push(urgent?&jobsurgent:&jobsnormal, headind, tailind);
	}
	if(pool.nidle){//no need to hold mutex to call cond_signal
		pthread_cond_broadcast(&pool.jobwait);
	}
}

/**
 * Wait for jobs in the count to be done. 
 */
void thread_pool_wait(tp_counter_t *counter, int urgent){
	while(counter->group){
		//if urgent is true, do not do do_job(0) to prevent slow down
		if(!do_job(urgent)){
			//do not use cond waiting. The wakeup is not gauranteed.
			mysleep(1e-6);
		}
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
	while(do_job(0));
	pthread_mutex_lock(&pool.mutex);
	if(pool.nidle+1<pool.ncur){/*some job is still doing the last bit. */
		pthread_cond_wait(&pool.idle, &pool.mutex);
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
		pthread_mutex_init(&pool.mutex, NULL);
		pthread_cond_init(&pool.idle, NULL);
		pthread_cond_init(&pool.jobwait, NULL);
		pthread_cond_init(&pool.exited, NULL);
		pthread_attr_init(&pool.attr);
		pthread_attr_setdetachstate(&pool.attr, PTHREAD_CREATE_DETACHED);
		pool.nidle=0;
		pool.ncur=1;/*counting the master thread. */
		pool.nmax=1;
		tid=1;//master thread
		jobsall_i=1;//start with 1. 0 indicate empty job.
		jobsall_count=nthread*100;
		jobsall=mycalloc(jobsall_count, jobs_t);
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
 *   Exit all threads and free thread pool.
 */
void thread_pool_destroy(void){
	//only 1 thread can wait
	if(atomic_sub_fetch(&pool.inited, 1)!=0) return;
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
	free(jobsall);
}
