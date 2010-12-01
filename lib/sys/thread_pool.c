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
   the end of the list while each thread will remove jobs from the beginning of
   the list when they are waked.  */

#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "../common.h"
#include "thread_pool.h"
/**
   struct for each thread. (private)
*/
typedef struct threads_t{
    pthread_t id;       /**<thread id*/
    thread_fun fun;     /**<the function to run;*/
    void *arg;          /**<the data for the function;*/
    struct thread_pool_t *pool;/**<pointer to the pool*/
}threads_t;
/**
   struct of jobs for the linked first in first out (fifo) list. (private)
*/
typedef struct jobs_t{
    thread_fun fun;     /**<The function*/
    void *arg;          /**<The argument*/
    struct jobs_t *next;/**<The pointer to the next entry*/
}jobs_t;
/**
   The thread pool struct. (public)
*/
struct thread_pool_t{
    pthread_mutex_t mutex; /**<the mutex.*/
    pthread_cond_t jobwait;/**<there are jobs jobwait.*/
    pthread_cond_t idle;  /**<all threads are idle*/
    pthread_cond_t exited; /**<all threads have exited*/
    threads_t **threads;  /**<the pool of threads (stack)*/
    jobs_t *jobshead;     /**<Start of the fifo list of jobwait jobs*/
    jobs_t *jobstail;     /**<End of the fifo list of jobwait jobs*/
    int icur;/**<the top of the threads stack.*/
    int nmax;/**<the maximum number of threads.*/
    int ncur;/**<the maximum number of live threads.*/
    int quit;/**<1: quit the threads.*/
    int nidle;/**<Number of idle threads*/
};

thread_pool_t *pool=NULL;//The default pool;
/**
   The working function in each thread.
*/
static void run_thread(threads_t *thread){
    thread_pool_t *pool=thread->pool;
    /*
       We should keep the mutex locked so that the signal won't come when we are
       not jobwait.
    */
    while(1){
	/*
	   Wait for the go signal. The mutex is release during jobwait and
	   locked immediately when cond is satisfied.
	 */
	pthread_mutex_lock(&pool->mutex);
	if(!pool->jobshead){
	    //this thread is idle. add to the idle count
	    pool->nidle++;
	    if(pool->nidle==pool->ncur){//all threads are idle
		pthread_cond_signal(&pool->idle);
	    }
	    //no more jobs to do, wait for the condition.
	    pthread_cond_wait(&pool->jobwait, &pool->mutex);

	    pool->nidle--;//no longer idle.
	    if(pool->nidle<0){
		error("Please report to author: nidle=%d\n", pool->nidle);
	    }
	    //jobshead should be not empty now.
	    if(!pool->jobshead && !pool->quit){//some other task already did the job.
		pthread_mutex_unlock(&pool->mutex);
		continue;
	    }
	}
	if(pool->quit){//call it quit.
	    pool->ncur--;
	    if(pool->ncur==0){//all thread are exited.
		pthread_cond_signal(&pool->exited);
	    }
	    free(thread);
	    pthread_mutex_unlock(&pool->mutex);
	    break;
	}
	//Take the jobs out of the job queue and run it.
	jobs_t *jobshead=pool->jobshead;
	thread->fun=jobshead->fun;
	thread->arg=jobshead->arg;
	pool->jobshead=jobshead->next;
	if(!pool->jobshead){//empty.
	    pool->jobstail=NULL;
	}
	free(jobshead);
	pthread_mutex_unlock(&pool->mutex);
	thread->fun(thread->arg);//run the job.
    }
}
/**
   Initialize the thread pool.
*/
void thread_pool_create(int nthread){
    pool=calloc(1, sizeof(thread_pool_t));
    pool->threads=calloc(nthread, sizeof(threads_t*));
    pthread_mutex_init(&pool->mutex,NULL);
    pthread_cond_init(&pool->idle, NULL);
    pthread_cond_init(&pool->jobwait, NULL);
    pthread_cond_init(&pool->exited, NULL);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED );
    pool->nidle=0;
    pool->nmax=nthread;
    pool->ncur=nthread;
    pool->jobshead=NULL;
    pool->jobstail=NULL;
    for(int ith=0; ith<nthread; ith++){
	pool->threads[ith]=calloc(1, sizeof(threads_t));
	pool->threads[ith]->pool=pool;
	if(pthread_create(&pool->threads[ith]->id, &attr, 
			  (thread_fun)run_thread, pool->threads[ith])){
	    error("Can not create thread\n");
	}
    }
}
/**
   Queue a job.
*/
void  thread_pool_queue(thread_fun fun, void *arg){
    if(!pool){
	error("Please call thread_pool_init first\n");
    }
    pthread_mutex_lock(&pool->mutex);
    //Add the job to the tail
    jobs_t *job=malloc(sizeof(jobs_t));
    job->fun=fun;
    job->arg=arg;
    job->next=NULL;
    if(pool->jobstail){//list is not empty.
	pool->jobstail->next=job;
	pool->jobstail=job;
    }else{
	//list is empty. need to signal the thresds.
	pool->jobshead=job;
        pool->jobstail=job;
	pthread_cond_broadcast(&pool->jobwait);
    }
    pthread_mutex_unlock(&pool->mutex);
}

/**
   Wait for all jobs to be done.
*/
void thread_pool_wait(void){
    /*
      We should lock mutex before compare nidle, otherwise another thread maybe
      modifying nidle, and emits pool->idle before we are ready to wait, thus
      hanging us here.
    */
    pthread_mutex_lock(&pool->mutex);
    if(pool->jobshead || pool->nidle<pool->nmax){
	pthread_cond_wait(&pool->idle, &pool->mutex);
    }
    if(pool->jobshead){
	error("Job is not empty\n");
    }
    pthread_mutex_unlock(&pool->mutex);
}

/**
   Exit all threads and free thread pool.
*/
void thread_pool_destroy(void){
    thread_pool_wait();
    pool->quit=1;
    pthread_cond_broadcast(&pool->jobwait);
    pthread_mutex_lock(&pool->mutex);
    if(pool->ncur>0){
	pthread_cond_wait(&pool->exited, &pool->mutex);
    }
    pthread_mutex_unlock(&pool->mutex);
}

