/*
  Copyright 2009-2026 Lianqi Wang
  
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
  be randomly lost. Replaced pointer based stack to integer based. Use
  jobsall[jobind].next instead of job->next. This solves the ABA problem by
  having a pop counter that is also an int. We didn't use pointer + counter as
  16 bit CAS is not always available.
*/

#include <unistd.h>
#include <pthread.h>
#include <limits.h>
#include <signal.h>
#include "common.h"
#include "misc.h" //mysleep
#include "thread.h"
#include "thread_pool.h"
__thread unsigned int tid=0;
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
/**The default pool;*/
static thread_pool_t pool;
static const unsigned int jobsall_mask=0xFFF;//number of jobs in jobsall. Must be all 1's
static int do_job(int urgent);
#define USE_QUEUE 0
/**
   Run a thread_t job with its own function. Used for jobs queued by QUEUE_THREAD.

   @param[in] arg The thread_t.
   @return The return value of its function.
*/
static void *run_thread_job(void *arg){
	thread_t *t=(thread_t *)arg;
	return t->fun(t);
}
#if USE_QUEUE //New implementation use queue. 2022-02-08. Not yet working.
/**
 * struct of jobs for the linked first in first out (fifo) list.
 * We use index here to emulate a 4 byte pointer
 */
typedef struct jobs_t{
	thread_fun fun; /**<The function*/
	void *arg;          /**<The argument*/
	tp_counter_t *counter;/**<address of the job group, whos value keeps track of queued jobs.*/
}jobs_t;

/**all jobs are allocated here. The first element has next set to 0 to indicate end of list*/
static jobs_t *jobsall=NULL;
static unsigned int jobsall_head=0;
static unsigned int jobsall_tail=0;
static unsigned int jobsall_njob=0;
#define ATOMIC_ORDER __ATOMIC_SEQ_CST
/**
   Enqueue a new job by copying it into a free slot of jobsall. Does jobs
   while the queue is full.
   All jobs are currently queued at the head (queueing urgent jobs at the tail
   is disabled as it is slow for small jobs).
   The jobsall[].fun can have the following states:
   	0: empty slot
   	1: enqueue is in progress
   	else: slot is ready to be processed

   @param[in] job Job to enqueue (copied).
   @param[in] urgent Passed to do_job() while waiting for a free slot.
*/
static void jobs_enqueue(jobs_t *job, int urgent){
	unsigned int ind, ntot2;
	thread_fun empty=NULL;
	thread_fun one=(thread_fun)1L;
	//First, check and occupy an empty slot.
	__atomic_load(&jobsall_njob, &ntot2, ATOMIC_ORDER);
	do{
		if(ntot2>=jobsall_mask){//if slot is not available
			do_job(urgent);
			__atomic_load(&jobsall_njob, &ntot2, ATOMIC_ORDER);//load again
		}
	} while(!__atomic_compare_exchange_n(&jobsall_njob, &ntot2, ntot2+1, 0, ATOMIC_ORDER, ATOMIC_ORDER));
	//Then determine the slot location.
	if(urgent && 0){//queue from the tail. subtraction wraps around
		//This is not efficient and also causes error below with counter.
		ind=jobsall_mask & __atomic_sub_fetch(&jobsall_tail, 1, ATOMIC_ORDER);
	}else{//queue from the head
		ind=jobsall_mask & __atomic_fetch_add(&jobsall_head, 1, ATOMIC_ORDER);
	}
	//keep retrying until we can fill the slot. Cannot skip it. An empty slot will block dequeue.
	//It may fail when the dequeue is still working on the last slot. (queue is full)
	while(!__atomic_compare_exchange_n(&jobsall[ind].fun, &empty, one, 0, ATOMIC_ORDER, ATOMIC_ORDER)){
		if(pool.nidle){
			pthread_cond_broadcast(&pool.jobwait);
		}
		empty=NULL;
	}
	jobsall[ind].arg=job->arg;
	jobsall[ind].counter=job->counter;
	//after setting the value, dequeue can take the node. We need to enfore memory order ? Doesn't seem to prevent retry
	__atomic_store_n(&jobsall[ind].fun, job->fun, ATOMIC_ORDER);
}
/**
   Dequeue a job from the tail. Set job if job is available.

   @param[out] job Dequeued job.
   @param[in] urgent Unused.
   @return 1 if a job was dequeued, 0 if no job is available.
*/
static unsigned int jobs_dequeue(jobs_t *job, int urgent){
	(void)urgent;
	unsigned int ind, ntot2;
	//First check and claim a filled slot.
	__atomic_load(&jobsall_njob, &ntot2, ATOMIC_ORDER);
	do{	//compare_exchange loads ntot to ntot2 if compare fails
		if(!ntot2){
			return 0;//no jobs available.
		}
	} while(!__atomic_compare_exchange_n(&jobsall_njob, &ntot2, ntot2-1, 0, ATOMIC_ORDER, ATOMIC_ORDER));
	//Determine the slot location
	ind=jobsall_mask&__atomic_fetch_add(&jobsall_tail, 1, ATOMIC_ORDER);
	//Read the slot information. This may fail even through we forced memory order (why?)
	do{
		__atomic_load(&jobsall[ind].fun, &job->fun, ATOMIC_ORDER);
	} while(job->fun==NULL||job->fun==(thread_fun)1);
	//after .fun is availabel, the following are all available.
	job->arg=jobsall[ind].arg;
	job->counter=jobsall[ind].counter;
	__atomic_store_n(&jobsall[ind].fun, NULL, ATOMIC_ORDER);
	return 1;
}
/**
   Do a job. If urgent is set, only do urgent jobs.

   @param[in] urgent Whether to only do urgent jobs.
   @return 1 if a job was done, 0 if no job is available.
*/
static int do_job(int urgent){
	jobs_t job={0};
	if(jobs_dequeue(&job, urgent)){
#if ENABLE_TP_TIMING
		double tk0=myclockd();
#endif		
		job.fun(job.arg);//run the job
		if(!__atomic_fetch_sub(&job.counter->group, 1, ATOMIC_ORDER)){
			//Why this one fails?
			error("already 0\n");
		}
#if ENABLE_TP_TIMING
		{
			unsigned int tk1=(unsigned int)ceil((myclockd()-tk0)*1e3);
			unsigned int old_value=job.counter->tmax;
			while(old_value<tk1&&!atomic_compare_exchange_n(&job.counter->tmax, &old_value, tk1));
			old_value=job.counter->tmin;
			unsigned int zero=0;
			if(!atomic_compare_exchange_n(&job.counter->tmin, &zero, tk1)){//not zero
				while(old_value>tk1&&!atomic_compare_exchange_n(&job.counter->tmin, &old_value, tk1));
			}
		}
#endif		
		return 1;//more jobs pending
	} else{
		return 0;//no more job
	}
}

/**
   Queue jobs that belongs to group denoted by counter. counter->group will be
   incremented by 1 for each job queued and decreased by 1 after each job is
   finished. Wait on it will clear when count is decreased to zero.

   Job is added to jobsnormal or joburgent depending on urgent value

   @param[in,out] counter Job group counter.
   @param[in] fun Function to call with arg for each job. If NULL, arg is an array of thread_t and each valid element's own fun is called with it.
   @param[in] arg Argument of fun, or array of njob thread_t if fun is NULL.
   @param[in] njob Number of jobs.
   @param[in] urgent Whether the jobs are urgent.
*/
void thread_pool_queue(tp_counter_t *counter, thread_fun fun, void *arg, int njob, int urgent){
	thread_t *arg2=(thread_t *)arg;
	jobs_t job={0};
	if(!jobsall){
		thread_pool_init(NTHREAD);
	}
	for(int ijob=0; ijob<njob; ijob++){
		if(!fun&&(!arg2[ijob].fun||arg2[ijob].start>=arg2[ijob].end)){
			continue;
		}
		job.counter=counter;
		atomic_add_fetch(&counter->group, 1);
		if(fun){
			job.arg=arg;
			job.fun=fun;
		} else{
			job.arg=&arg2[ijob];
			job.fun=run_thread_job;
		}
		jobs_enqueue(&job, urgent);//take idle job from the pool using lock free approach. 
	
		if(ijob%100==0 && pool.nidle){//no need to hold mutex to call cond_signal
			pthread_cond_broadcast(&pool.jobwait);
		}
	}
}
#else //use array based stack
typedef struct jobs_t{
	thread_fun fun; /**<The function*/
	void *arg;          /**<The argument*/
	tp_counter_t *counter;/**<address of the job group, whos value keeps track of queued jobs.*/
	unsigned int next;  /**<4 byte index into jobsall for next job, replacing pointer*/
}jobs_t;

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

/**all jobs are allocated here. The first element has next set to 0 to indicate end of list*/
static jobs_t *jobsall=NULL;
/**Start of the fifo list of pending jobs*/
static jobshead_t jobsnormal={0L};
/**Start of the fifo list of urgent jobs*/
static jobshead_t jobsurgent={0L};
/**saving unused jobs_t*/
static jobshead_t jobspool={0L};


/**
   Place job list (jobhead to jobtail) to the beginning of list (head) by connecting jobtail.next to head and replacing
   head by jobhead.

   @param[in,out] head Head of the list.
   @param[in] jobheadind Index in jobsall of the first job of the list to push.
   @param[in] jobtailind Index in jobsall of the last job of the list to push.
*/
static void jobs_push(jobshead_t *head, unsigned int jobheadind, unsigned int jobtailind){
	do{
		jobsall[jobtailind].next=head->head;
	}while(!__atomic_compare_exchange_n(&(head->head), &jobsall[jobtailind].next, jobheadind, 
										1, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE));
}
/**
   Get a job from the beginning of the list.

   @param[in,out] head Head of the list.
   @return Index in jobsall of the job removed, or 0 if the list is empty.
*/
static unsigned int jobs_pop(jobshead_t *head){
	jobshead_t job, job2;
	job.state=head->state;
	do{
		//compare both counter and pointer to make sure the node is not changed.
		job2.counter=job.counter+1;//increse the counter
		job2.head=jobsall[job.head].next;
	} while(!__atomic_compare_exchange(&head->state, &job.state, &job2.state, 
									   1, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE));
	return job.head;
}

/**
   Do a job. If urgent is set, only do urgent jobs.

   @param[in] urgent Whether to only do urgent jobs.
   @return 1 if a job was done, 0 if no job is available.
*/
static int do_job(int urgent){
	unsigned int jobind=0;
	if((jobind=jobs_pop(&jobsurgent))||(!urgent&&(jobind=jobs_pop(&jobsnormal)))){
		jobs_t *job=&jobsall[jobind];
#if ENABLE_TP_TIMING
		double tk0=myclockd();
#endif		
		job->fun(job->arg);//run the job
#if ENABLE_TP_TIMING
		{
			unsigned int tk1=(unsigned int)ceil((myclockd()-tk0)*1e3);
			unsigned int old_value=job->counter->tmax;
			while(old_value<tk1&&!atomic_compare_exchange_n(&job->counter->tmax, &old_value, tk1));
			old_value=job->counter->tmin;
			unsigned int zero=0;
			if(!atomic_compare_exchange_n(&job->counter->tmin, &zero, tk1)){//not zero
				while(old_value>tk1&&!atomic_compare_exchange_n(&job->counter->tmin, &old_value, tk1));
			}
		}
#endif	
		if(__atomic_sub_fetch(&job->counter->group, 1, __ATOMIC_ACQUIRE)==4294967295){
			warning("group %p is now -1.\n", job->counter);
		}
		jobs_push(&jobspool, jobind, jobind);//return job to the pool.	
		return 1;//more jobs pending
	}
	return 0;//no more job
}

/**
   Queue jobs that belongs to group denoted by counter. counter->group will be
   incremented by 1 for each job queued and decreased by 1 after each job is
   finished. Wait on it will clear when count is decreased to zero.

   Job is added to jobsnormal or joburgent depending on urgent value

   @param[in,out] counter Job group counter.
   @param[in] fun Function to call with arg for each job. If NULL, arg is an array of thread_t and each valid element's own fun is called with it.
   @param[in] arg Argument of fun, or array of njob thread_t if fun is NULL.
   @param[in] njob Number of jobs.
   @param[in] urgent Whether the jobs are urgent.
*/
void thread_pool_queue(tp_counter_t *counter, thread_fun fun, void *arg, int njob, int urgent){
	unsigned int headind=0, tailind=0;
	thread_t *arg2=(thread_t *)arg;
	if(!jobsall){
		extern int NTHREAD;
		thread_pool_init(NTHREAD);
	}
	for(int ijob=0; ijob<njob; ijob++){
		unsigned int jobind=0;
		jobs_t *job=NULL;
		if(fun || (arg2[ijob].fun && arg2[ijob].start < arg2[ijob].end)){//valid job
			jobind=jobs_pop(&jobspool);//take idle job from the pool using lock free approach. 
			if(jobind==0){//no slot is available.
				ijob--;
			}else{
				job=&jobsall[jobind];
				if(fun){
					job->fun=fun;
					job->arg=arg;
				} else{
					job->fun=run_thread_job;
					job->arg=&arg2[ijob];
				}
				job->counter=counter;
				job->next=0;
				//_ATOMIC_RELEASE after store, and _ATOMIC_ACQUIRE before load. Always in pairs.
				//_ATOMIC_RELEASE ensures operations before this line is not moved to after.
				__atomic_add_fetch(&counter->group, 1, __ATOMIC_RELEASE);
				if(!headind){
					headind=jobind;
				} else{
					jobsall[tailind].next=jobind;
				}
				tailind=jobind;
			}
		}
		if(((ijob%100==99||ijob+1==njob)&&headind)||!jobind){
			if(headind){
				jobs_push(urgent?&jobsurgent:&jobsnormal, headind, tailind);
				headind=0; 
				tailind=0;
			}
			if(pool.nidle){//no need to hold mutex to call cond_signal
				pthread_cond_broadcast(&pool.jobwait);
			}
			if(!jobind){
				do_job(urgent);
			}
		}
	}
}

#endif



/**
   The working function in each thread which never quits until pool is being
   resized or destroyed. The master thread does not run this routine.

   @param[in] data Thread index cast to pointer. Stored in tid.
   @return NULL.
*/
static void *run_thread(void *data){
	(void)data;
	tid=(int)(long)data;
	while(1){
		while(do_job(0));//do until no jobs are left
		pthread_mutex_lock(&pool.mutex);//acquire lock before modifying pool and wait. 	
		pool.nidle++;//increment idle counter. no need atomic op. protected by mutex.
		if(pool.nidle+1==pool.ncur){//at most ncur-1 threads can be idle as the main thread counts as 1.
			pthread_cond_broadcast(&pool.idle);//all threads are idle 
		}
		if(pool.ncur>pool.nmax){//quit this thread
			pool.ncur--;
			pool.nidle--;
			if(pool.ncur==1){//all thread are exited except the master thread. 
				pthread_cond_broadcast(&pool.exited);
			}
			pthread_mutex_unlock(&pool.mutex);
			break;
		}
		pthread_cond_wait(&pool.jobwait, &pool.mutex);//lock is released during wait and reacquired after wake
		pool.nidle--;//no longer idle. 
		pthread_mutex_unlock(&pool.mutex);//release the lock
	}
	return NULL;
}

/**
   Wait for jobs in the count to be done. Does jobs while waiting. Has the
   potential to deadlock if it calls do_job and the job inside calls
   pthread_cond_wait.

   @param[in] counter Job group counter to wait for.
   @param[in] urgent Whether to only do urgent jobs while waiting.
*/
void thread_pool_wait(tp_counter_t *counter, int urgent){
	while(counter->group){
		//if urgent is true, only do urgent jobs.
		do_job(urgent);
		//do not use cond waiting. The wakeup is not gauranteed.
	}
}
/**
   Wait for all jobs to be done.
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
	//while(do_job(0));//causes segmentation fault.
	pthread_mutex_lock(&pool.mutex);
	while(pool.nidle+1<pool.ncur){/*some job is still doing the last bit. */
		pthread_cond_wait(&pool.idle, &pool.mutex);
	}
	pthread_mutex_unlock(&pool.mutex);
}
/**
   Initialize the thread pool. Repeated call to resize the pool

   @param[in] nthread Number of threads including the master thread. 0 or more than the number of processors uses the
   number of processors.
*/
void thread_pool_init(unsigned int nthread){
	if(!jobsall){
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
		tid=pool.ncur;
		jobsall=mycalloc(jobsall_mask+1, jobs_t);
#if USE_QUEUE
		jobsall_head=0;
		jobsall_tail=0;
		jobsall_njob=0;
#else
		jobsnormal.state=0;
		jobsurgent.state=0;
		jobspool.state=0;
		for(unsigned int i=1; i<=jobsall_mask; i++){//put all into the pool
			jobs_push(&jobspool, i, i);
		}
#endif


	}
	unsigned int nthread_max=(unsigned int)sysconf(_SC_NPROCESSORS_ONLN);
	if(nthread<=0) nthread=nthread_max;
	if(nthread>nthread_max) nthread=nthread_max;
	pthread_mutex_lock(&pool.mutex);
	pool.nmax=nthread;
	if(pool.ncur<pool.nmax){/*need to launch more threads. we change ncur here instead of in new thread for consistency. */
		pthread_t thread;/*we don't need the thread information. */
		for(; pool.ncur<pool.nmax; pool.ncur++){
			if(pthread_create(&thread, &pool.attr, run_thread, (void*)(long)pool.ncur)){
				error("Can not create thread\n");
				pool.nmax=pool.ncur;//only count threads that are running
				break;
			}
		}
	} else if(pool.ncur>pool.nmax){/*need to quit some threads. */
		pthread_cond_broadcast(&pool.jobwait);/*wake up all threads. */
	}
	pthread_mutex_unlock(&pool.mutex);
}

/**
   Exit all threads and free thread pool.
*/
void thread_pool_destroy(void){
	//only 1 thread can destroy the pool
	if(!jobsall) return;
	/*tell all jobs to quit. */
	pthread_mutex_lock(&pool.mutex);
	pool.nmax=1;
	pthread_cond_broadcast(&pool.jobwait);
	while(pool.ncur>1){
		pthread_cond_wait(&pool.exited, &pool.mutex);
	}
	pthread_mutex_unlock(&pool.mutex);
	pthread_mutex_destroy(&pool.mutex);
	pthread_cond_destroy(&pool.idle);
	pthread_cond_destroy(&pool.jobwait);
	pthread_cond_destroy(&pool.exited);
	pthread_attr_destroy(&pool.attr);
	free(jobsall);
	jobsall=0;
}
