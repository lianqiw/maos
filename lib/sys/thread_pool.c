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
/**
   \file thread_pool.c
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
 */

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include "../common.h"
#include "../misc.h"
#include "thread.h"
#include "thread_pool.h"
#include "process.h"
/**
   struct of jobs for the linked first in first out (fifo) list. (private)
*/
typedef struct jobs_t{
    thread_fun fun;     /**<The function*/
    void *arg;          /**<The argument*/
    long *count;        /**<address of the job group, which is also a counter.*/
    int urgent;         /**<whether this job is urgent (can not be interruped)*/
    struct jobs_t *next;/**<The pointer to the next entry*/
}jobs_t;
/**
   The thread pool struct. (public)
*/
struct thread_pool_t{
    pthread_mutex_t mutex; /**<the mutex.*/
    pthread_spinlock_t spin; /**<the mutex.*/
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
};
static thread_pool_t pool;//The default pool;

/**
   Do a job. The mutex should have already been locked when calling this
routine. It will take a job from the head of the job queue, release the mutex,
run the job, and then acquire the mutex, and check whether the job is
finished. */
static inline void do_job(void) {
    jobs_t *job;
    job=pool.jobshead;
    pool.jobshead=job->next;
    pthread_mutex_unlock(&pool.mutex);
    job->fun(job->arg); //run the job
    pthread_mutex_lock(&pool.mutex);
    (*job->count)--;//decrease the count.
    if(!(*job->count)){
	/*job is done. need broadcast since multiple threads may be waiting for
	  their job to finish*/
	pthread_cond_broadcast(&pool.jobdone);
    }
    pthread_spin_lock(&pool.spin);
    job->next=pool.jobspool;
    pool.jobspool=job;
    pthread_spin_unlock(&pool.spin);
}

/**
   The working function in each thread.
*/
static void run_thread(){
    while(!pool.quit){
	/*
	  Wait for the go signal. The mutex is release during jobwait and
	  locked immediately when cond is satisfied.
	*/
	pthread_mutex_lock(&pool.mutex);
	if(!pool.jobshead){
	    //this thread is idle. add to the idle count
	    pool.nidle++;
	    //at most ncur-1 threads can be idle.
	    if(pool.nidle+1==pool.ncur){//all threads are idle
		pthread_cond_signal(&pool.idle);
	    }
	    //no more jobs to do, wait for the condition.
	    pthread_cond_wait(&pool.jobwait, &pool.mutex);
	    pool.nidle--;//no longer idle.
	    //jobshead should be not empty now.
	    if(!pool.jobshead){//some other task already did the job.
		pthread_mutex_unlock(&pool.mutex);
		continue;
	    }
	}
	//Take the jobs out of the job queue and run it.
	while(pool.jobshead){
	    do_job();
	}
	pthread_mutex_unlock(&pool.mutex);
    }
    pool.ncur--;
    if(pool.ncur==1){//all thread are exited except the master thread.
	pthread_cond_signal(&pool.exited);
    }
    pthread_mutex_unlock(&pool.mutex);
}
/**
   Create a new thread
*/
static void thread_new(){//the caller is responsible to lock mutex.
    pthread_t thread;//we don't need the thread information.
    if(pthread_create(&thread, &pool.attr, (thread_fun)run_thread, NULL)){
	error("Can not create thread\n");
    }
    pool.ncur++;
}
/**
   Initialize the thread pool.
*/
void thread_pool_init(int nthread){
    memset(&pool, 0, sizeof(thread_pool_t));
    pthread_mutex_init(&pool.mutex,NULL);
    pthread_spin_init(&pool.spin, 0);
    pthread_cond_init(&pool.idle, NULL);
    pthread_cond_init(&pool.jobwait, NULL);
    pthread_cond_init(&pool.jobdone, NULL);
    pthread_cond_init(&pool.exited, NULL);
    pthread_attr_init(&pool.attr);
    pthread_attr_setdetachstate(&pool.attr,PTHREAD_CREATE_DETACHED);
    pool.nidle=0;
    pool.nmax=nthread;
    pool.jobspool=NULL;
    pool.jobshead=NULL;
    pool.jobstail=NULL;
    pool.ncur=1;//counting the master thread.
    pthread_mutex_lock(&pool.mutex);
    for(int ith=0; ith<nthread-1; ith++){
	//launch nthread-1 threads because the master thread is another thread.
	thread_new();
    }
    pthread_mutex_unlock(&pool.mutex);
}
/**
   Queue a job that belongs to group denoted by group. The argument count, will
   be incremented by 1 when queued and decreased by 1 if job is finished. Wait
   on it will clear when count is decreased to zero. */
void thread_pool_queue(long *group, thread_fun fun, void *arg, int urgent){
    //Add the job to the head if urgent>0, otherwise to the tail.
    jobs_t *job;
    if(pool.jobspool){//take it from the pool.
	pthread_spin_lock(&pool.spin);
	job=pool.jobspool;
	pool.jobspool=pool.jobspool->next;
	pthread_spin_unlock(&pool.spin);
    }else{
	job=malloc(sizeof(jobs_t));
    }
    job->fun=fun;
    job->arg=arg;
    job->count=group;
    job->urgent=urgent;
    job->next=NULL;
    pthread_mutex_lock(&pool.mutex);
    (*job->count)++;
    if(pool.jobshead){//list is not empty
	if(urgent){
	    //add to head
	    job->next=pool.jobshead;
	    pool.jobshead=job;
	}else{     
	    //add to tail.
	    pool.jobstail->next=job;
	    pool.jobstail=job;
	}
    }else{
	//list is empty. need to signal the thresds.
	pool.jobshead=job;
	pool.jobstail=job;
    }
    if(pool.nidle>0){
	pthread_cond_signal(&pool.jobwait);//wake up one thread only.
	pthread_cond_signal(&pool.jobdone);//wake up the thread that is waiting for jobdone.
    }
    pthread_mutex_unlock(&pool.mutex);
}
/**
   Queue njob jobs with same arguments.
 */
void thread_pool_queue_many_same(long *group, thread_fun fun, void *arg, int njob, int urgent){
    /*
       First create a list of all the jobs without acquiring lock.
    */
    jobs_t *head=NULL;
    jobs_t *tail=NULL;
    pthread_spin_lock(&pool.spin);
    for(int ijob=0; ijob<njob; ijob++){
	jobs_t *job;
	if(pool.jobspool){
	    job=pool.jobspool;
	    pool.jobspool=pool.jobspool->next;
	}else{
	    job=malloc(sizeof(jobs_t));
	}
	if(!ijob){
	    tail=job;
	}
	job->fun=fun;
	job->arg=arg;
	job->count=group;
	job->urgent=urgent;
	job->next=head;
	head=job;
    }
    pthread_spin_unlock(&pool.spin);
    //Add the job to queue
    pthread_mutex_lock(&pool.mutex);
    (*group)+=njob;
    if(pool.jobshead){//list is not empty
	if(urgent){//add to head
	    tail->next=pool.jobshead;
	    pool.jobshead=head;
	}else{//add to tail.
	    pool.jobstail->next=head;
	    pool.jobstail=tail;
	}
    }else{
	//list is empty. need to signal the thresds.
	pool.jobshead=head;
	pool.jobstail=tail;
    }
    pthread_mutex_unlock(&pool.mutex);
    if(pool.nidle){
	pthread_cond_broadcast(&pool.jobwait);//wake up all idle threads.
	pthread_cond_broadcast(&pool.jobdone);//wake up the thread that is waiting for jobdone.
    }
}
/**
   Queue njob jobs with arguments array.
 */
void thread_pool_queue_many(long *group, thread_t *arg, int njob, int urgent){
    /*
       First create a list of all the jobs without acquiring lock.
    */
    jobs_t *head=NULL;
    jobs_t *tail=NULL;
    pthread_spin_lock(&pool.spin);
    for(int ijob=0; ijob<njob; ijob++){
	jobs_t *job;
	if(pool.jobspool){
	    job=pool.jobspool;
	    pool.jobspool=pool.jobspool->next;
	}else{
	    job=malloc(sizeof(jobs_t));
	}
	if(!ijob){//first.
	    tail=job;
	}
	job->fun=(thread_fun)arg[ijob].fun;
	job->arg=arg+ijob;
	job->count=group;
	job->urgent=urgent;
	job->next=head;
	head=job;
    }
    pthread_spin_unlock(&pool.spin);
    //Add the job to queue
    pthread_mutex_lock(&pool.mutex);
    (*group)+=njob;
    if(pool.jobshead){//list is not empty
	if(urgent){//add to head
	    tail->next=pool.jobshead;
	    pool.jobshead=head;
	}else{//add to tail.
	    pool.jobstail->next=head;
	    pool.jobstail=tail;
	}
    }else{
	//list is empty. need to signal the thresds.
	pool.jobshead=head;
	pool.jobstail=tail;
    }
    if(pool.nidle){
	pthread_cond_broadcast(&pool.jobwait);//wake up all threads.
	pthread_cond_broadcast(&pool.jobdone);//wake up the thread that is waiting for jobdone.
    }
    pthread_mutex_unlock(&pool.mutex);
}
/**
   Wait for jobs in the count to be done.
*/
void thread_pool_wait(long *count){
    pthread_mutex_lock(&pool.mutex);
    //Do some job while we are waiting.
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
   Wait for all jobs to be done.
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
    if(pool.nidle+1<pool.ncur){//some job is still doing the last bit.
       	/*
	  thread_pool_wait may not obtain the lock in pool.mutex immediately
	  after signal on pool.idle is emmited. thread_pool_queue may obtain
	  the lock. So we wait in a loop.*/
	pthread_cond_wait(&pool.idle, &pool.mutex);
    }
    pthread_mutex_unlock(&pool.mutex);
}

/**
   Exit all threads and free thread pool.
*/
void thread_pool_destroy(void){
    thread_pool_wait_all();//let all jobs finish.
    //tell all jobs to quit.
    pool.quit=1;
    pthread_cond_broadcast(&pool.jobwait);
    pthread_mutex_lock(&pool.mutex);
    if(pool.ncur>0){
	pthread_cond_wait(&pool.exited, &pool.mutex);
    }
    pthread_mutex_unlock(&pool.mutex);
}
