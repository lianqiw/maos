/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifdef __linux__
#define _GNU_SOURCE //for SYSCALL
#endif
#include <signal.h>
#include "common.h"
#include "thread.h"
/**
   Functions regarding to threading
*/
PNEW2(mutex_fftw);

/**
   Create a new thread and forget.
*/
int thread_new(thread_fun fun, void* arg){
	int ans;
	pthread_t temp;
	pthread_attr_t attr;
	if((ans=pthread_attr_init(&attr))){
		error("pthread_attr_init failed with %d\n", ans);
	}
	if((ans=pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED))){
		error("pthread_attr_setdetachstate failed with %d\n", ans);
	}
	if((ans=pthread_create(&temp, &attr, fun, arg))){
		error("pthread_create failed with answer %d\n", ans);
	}
	return ans;
}
/**
   Break out the job to be executed by multiple threads.
*/
thread_t* thread_prep(long start, long end, long nthread,
	thread_wrapfun fun, void* data){
	if(nthread==0) return NULL;
	thread_t *thd=mycalloc(nthread, thread_t);
	long nt=(end-start+nthread-1)/nthread;
	int skip=0;
	int nthread_active=0; //active threads
	for(long ithread=0; ithread<nthread; ithread++){
		thd[ithread].ithread=ithread;
		thd[ithread].nthread=nthread;
		thd[ithread].data=data;
		thd[ithread].fun=fun;
		if(skip){//skip the rest
			thd[ithread].start=0;
			thd[ithread].end=0;
		} else{
			thd[ithread].start=start;
			start+=nt;
			if(start>=end){
				start=end;
				skip=1;
				data=0;
				fun=0;
				nthread_active=ithread+1;
			}
			thd[ithread].end=start;
		}
	}
	thd[0].nthread=nthread_active;
	/*Make sure we terminate at the right place. */
	if(thd[nthread-1].end&&thd[nthread-1].end!=end){
		error("Not correctly terminated\n");
	}
	return thd;
}

/**
   After calling this routine, this thread will nolonger receive signals
 */
void thread_block_signal(){
	sigset_t set;
	sigfillset(&set);
	pthread_sigmask(SIG_BLOCK, &set, NULL);
}
#if _OPENMP
long thread_id(){
	return (long)omp_get_thread_num();
}
#else
#ifdef __linux__
#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>
long thread_id(){
	return (long)syscall(SYS_gettid);
}
#else
long thread_id(){
	return (long)pthread_self();
}
#endif
#endif
