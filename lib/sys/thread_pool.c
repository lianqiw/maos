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
   Contains implementation of a simple thread pool.

   A thread pool is a pool of threads that can be used to run functions without
   spawning/ending a new thread. The threads are created with no job designated,
   but waiting for a condition.

   The jobs are keep in a linked list. thread_pool_queue will queue the jobs to
   the head or end of the list while each thread will remove jobs from the beginning of
   the list when they are waked.  

   2010-12-13: Improved the thread scheduling by test pool->nidle instead of
   pool->jobstail. This reduces the tomography of 100 steps down to 0.65s from
   0.8s for 8 threads at taurus.
 */

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "../common.h"
#include "thread.h"
#include "thread_pool.h"
/**
   struct of jobs for the linked first in first out (fifo) list. (private)
*/
typedef struct jobs_t{
    thread_fun fun;     /**<The function*/
    void *arg;          /**<The argument*/
    long *count;        /**<address of the job group, which is also a counter.*/
    struct jobs_t *next;/**<The pointer to the next entry*/
}jobs_t;
/**
   The thread pool struct. (public)
*/
struct thread_pool_t{
    pthread_mutex_t mutex; /**<the mutex.*/
    pthread_cond_t jobwait;/**<there are jobs jobwait.*/
    pthread_cond_t idle;   /**<all threads are idle*/
    pthread_cond_t exited; /**<all threads have exited*/
    pthread_attr_t attr;   /**<The attribution of newly created threads.*/
    //pthread_t *threads;   /**<the pool of threads (stack)*/
    jobs_t *jobshead;      /**<Start of the fifo list of jobwait jobs*/
    jobs_t *jobstail;      /**<End of the fifo list of jobwait jobs*/
    int icur; /**<the top of the threads stack.*/
    int nmax; /**<the maximum number of threads. constant.*/
    int ncur; /**<the maximum number of live threads, excluding the master thread.*/
    int quit; /**<1: quit the threads.*/
    int nidle;/**<Number of idle threads, excluding the master thread.*/
    int nwait;/**<Number of jobs that are waiting for jobs done.*/
};

thread_pool_t *pool=NULL;//The default pool;
/**
   The working function in each thread.
*/
static void run_thread(void){
    /*
       We should keep the mutex locked so that the signal won't come when we are
       not jobwait.
    */
    jobs_t *job;
    while(1){
	/*
	   Wait for the go signal. The mutex is release during jobwait and
	   locked immediately when cond is satisfied.
	 */
	pthread_mutex_lock(&pool->mutex);
	if(pool->quit){//quit everything
	    break;
	}
	while(pool->nidle<0){//too many active threads.
	    pool->nidle++;
	    pthread_cond_broadcast(&pool->idle);//there is an idling thread.
	    pthread_cond_wait(&pool->jobwait, &pool->mutex);
	    pool->nidle--;
	}
	if(!pool->jobshead){ //this thread is idle. add to the idle count
	    pool->nidle++;
	    if(pool->nidle==pool->ncur){//all threads are idle
		pthread_cond_broadcast(&pool->idle);
		/*
		  thread_pool_wait may not obtain the lock in pool->mutex
		  immediately after we release the lock on mutex.
		  thread_pool_queue may obtain the lock. */
	    }
	    //no more jobs to do, wait for the condition.
	    pthread_cond_wait(&pool->jobwait, &pool->mutex);
	    pool->nidle--;//no longer idle.
	    if(pool->nidle<0){
		error("Please report to author: nidle=%d\n", pool->nidle);
	    }
	    //jobshead should be not empty now.
	    if(!pool->jobshead){//some other task already did the job.
		pthread_mutex_unlock(&pool->mutex);
		continue;
	    }
	}

	//Take the jobs out of the job queue and run it.
    repeat:
	job=pool->jobshead;
	pool->jobshead=job->next;
	if(!pool->jobshead){//empty.
	    pool->jobstail=NULL;
	}
	pthread_mutex_unlock(&pool->mutex);
	job->fun(job->arg);//run the job
	pthread_mutex_lock(&pool->mutex);
	(*job->count)--;//decrease the count.
	if(!(*job->count)){//job is done.
	    pthread_cond_broadcast(&pool->idle);
	}
	free(job);
	if(pool->jobshead){
	    goto repeat;//by pass the while loop and a cycle of unlock/lock.
	}
	pthread_mutex_unlock(&pool->mutex);
    }
    pool->ncur--;
    if(pool->ncur==1){//all thread are exited except the master thread.
	pthread_cond_signal(&pool->exited);
    }
    pthread_cond_broadcast(&pool->idle);
    pthread_mutex_unlock(&pool->mutex);
}
/**
   Create a new thread
*/
static void thread_new(){//the caller is responsible to lock mutex.
    pool->ncur++;
    pthread_t thread;//we don't need the thread information.
    if(pthread_create(&thread, &pool->attr, (thread_fun)run_thread, NULL)){
	error("Can not create thread\n");
    }
}
/**
   Initialize the thread pool.
*/
void thread_pool_create(int nthread){
    pool=calloc(1, sizeof(thread_pool_t));
    pthread_mutex_init(&pool->mutex,NULL);
    pthread_cond_init(&pool->idle, NULL);
    pthread_cond_init(&pool->jobwait, NULL);
    pthread_cond_init(&pool->exited, NULL);
    pthread_attr_init(&pool->attr);
    pthread_attr_setdetachstate(&pool->attr,PTHREAD_CREATE_DETACHED);
    pool->nidle=0;
    pool->nmax=nthread;
    pool->ncur=1;//counting the master thread.
    pool->jobshead=NULL;
    pool->jobstail=NULL;
    pthread_mutex_lock(&pool->mutex);
    for(int ith=0; ith<nthread; ith++){
	thread_new();
    }
    pthread_mutex_unlock(&pool->mutex);
}
/**
   Queue a job that belongs to group denoted by group. The argument count, will
   be incremented by 1 when queued and decreased by 1 if job is finished. Wait
   on it will clear when count is decreased to zero. */
void thread_pool_queue(long *group, thread_fun fun, void *arg, int athead){
    //Add the job to the tail
    jobs_t *job=malloc(sizeof(jobs_t));
    job->count=group;
    job->fun=fun;
    job->arg=arg;
    job->next=NULL;
    pthread_mutex_lock(&pool->mutex);
    (*job->count)++;
    if(pool->jobstail){//list is not empty
	if(athead){//add to head
	    job->next=pool->jobshead;
	    pool->jobshead=job;
	}else{//add to tail.
	    pool->jobstail->next=job;
	    pool->jobstail=job;
	}
    }else{
	//list is empty. need to signal the thresds.
	pool->jobshead=job;
	pool->jobstail=job;
    }
    if(pool->nidle){
	pthread_cond_signal(&pool->jobwait);//wake up one thread only.
    }
    pthread_mutex_unlock(&pool->mutex);
}
/**
   Queue njob jobs with same arguments.
 */
void thread_pool_queue_many_same(long *group, thread_fun fun, void *arg, int njob, int athead){
    /*
       First create a list of all the jobs without acquiring lock.
    */
    jobs_t *head=NULL;
    jobs_t *tail=NULL;
    for(int ijob=0; ijob<njob; ijob++){
	jobs_t *job=malloc(sizeof(jobs_t));
	if(!ijob){
	    tail=job;
	}
	job->count=group;
	job->fun=fun;
	job->arg=arg;
	job->next=head;
	head=job;
    }
    //Add the job to queue
    pthread_mutex_lock(&pool->mutex);
    (*group)+=njob;
    if(pool->jobstail){//list is not empty
	if(athead){//add to head
	    tail->next=pool->jobshead;
	    pool->jobshead=head;
	}else{//add to tail.
	    pool->jobstail->next=head;
	    pool->jobstail=tail;
	}
    }else{
	//list is empty. need to signal the thresds.
	pool->jobshead=head;
	pool->jobstail=tail;
    }
    if(pool->nidle){
	pthread_cond_broadcast(&pool->jobwait);//wake up all threads.
    }
    pthread_mutex_unlock(&pool->mutex);
}/**
   Queue njob jobs with arguments array.
 */
void thread_pool_queue_many(long *group, thread_t *arg, int njob, int athead){
    /*
       First create a list of all the jobs without acquiring lock.
    */
    jobs_t *head=NULL;
    jobs_t *tail=NULL;
    for(int ijob=0; ijob<njob; ijob++){
	jobs_t *job=malloc(sizeof(jobs_t));
	if(!ijob){
	    tail=job;
	}
	job->count=group;
	job->fun=(thread_fun)arg[ijob].fun;
	job->arg=arg+ijob;
	job->next=head;
	head=job;
    }
    //Add the job to queue
    pthread_mutex_lock(&pool->mutex);
    (*group)+=njob;
    if(pool->jobstail){//list is not empty
	if(athead){//add to head
	    tail->next=pool->jobshead;
	    pool->jobshead=head;
	}else{//add to tail.
	    pool->jobstail->next=head;
	    pool->jobstail=tail;
	}
    }else{
	//list is empty. need to signal the thresds.
	pool->jobshead=head;
	pool->jobstail=tail;
    }
    if(pool->nidle){
	pthread_cond_broadcast(&pool->jobwait);//wake up all threads.
    }
    pthread_mutex_unlock(&pool->mutex);
}
/**
   Wait for jobs in the count to be done.
*/
void thread_pool_wait(long *count){
    pthread_mutex_lock(&pool->mutex);
    pool->nwait++;
    pool->nidle++;//Mark this thread as sleeping when we are waiting.
    if(pool->jobshead && pool->nwait+pool->nmax > pool->ncur+1){
	//not enough active threads, create one.
	thread_new();
    }
    while(*count>0){//there are jobs to be done.
	pthread_cond_wait(&pool->idle, &pool->mutex);
    }
    //Now we are ready to proceed.
    pool->nidle--;//we are active
    pool->nwait--;
    while(pool->nidle<0){//too many active threads.
	//notify threads that we are waiting.
	pthread_cond_broadcast(&pool->jobwait);
	pthread_cond_wait(&pool->idle, &pool->mutex);
    }
    pthread_mutex_unlock(&pool->mutex);
}
/**
   Wait for all jobs to be done.
*/
void thread_pool_wait_all(void){
    /*
      We should lock mutex before compare nidle, otherwise another thread maybe
      modifying nidle, and emits pool->idle before we are ready to wait, thus
      hanging us here.
    */
    pthread_mutex_lock(&pool->mutex);
    pool->nidle++;//Mark this thread as sleeping when we are waiting. important.
    while(pool->jobshead || pool->nidle<pool->ncur){
       	/*
	  thread_pool_wait may not obtain the lock in pool->mutex immediately
	  after signal on pool->idle is emmited. thread_pool_queue may obtain
	  the lock. So we wait in a loop.*/
	pthread_cond_wait(&pool->idle, &pool->mutex);
    }
    pool->nidle--;//we are active
    pthread_mutex_unlock(&pool->mutex);
}

/**
   Exit all threads and free thread pool.
*/
void thread_pool_destroy(void){
    thread_pool_wait_all();
    pool->quit=1;
    pthread_cond_broadcast(&pool->jobwait);
    pthread_mutex_lock(&pool->mutex);
    if(pool->ncur>0){
	pthread_cond_wait(&pool->exited, &pool->mutex);
    }
    pthread_mutex_unlock(&pool->mutex);
}

#ifdef TEST
void* nothing2(void *a){
    (void)a;
    double b;
    for(int i=0; i<1000; i++){
    	b=i*i;
    }
    return NULL;
}
void* nothing(void* a){
    (void)a;
    int i=1;
    i=i+1;
    long group=0;
    for(int icase=0; icase<10; icase++){
	thread_pool_queue(&group, (thread_fun) nothing2, NULL, 1);
    }
    thread_pool_wait(&group);
    return NULL;
}
int main(int argc, char **argv){
    int nthread=1;
    if(argc>1){
	nthread=strtod(argv[1], NULL);
    }
    int ncase=1000;
    if(argc>2){
	ncase=strtod(argv[2], NULL);
    }
    thread_pool_create(nthread);
    long group=0;
    info("Running %d cases with %d threads\n", ncase, nthread);
    //TIC;.tic;
     thread_pool_queue_many_same(&group, (thread_fun)nothing,NULL,ncase, 1);
    //info("Queued");
    
    for(int icase=0; icase<ncase; icase++){
	//thread_pool_queue(&group,(thread_fun)nothing, NULL,0);
	}
    thread_pool_wait(&group);
    thread_pool_wait_all();
    info("done\n");
}
#endif
#ifdef TEST2
int main(int argc, char **argv){
    int nthread=100000000;
    if(argc>1){
	nthread=strtod(argv[1], NULL);
    }
    pthread_mutex_t mutex;
    pthread_mutex_init(&mutex, NULL);
    long i=0;
    for(int ithread=0; ithread<nthread; ithread++){
	//	pthread_mutex_lock(&mutex);
	//pthread_mutex_unlock(&mutex);
	i=i+1;
    }
}
#endif
