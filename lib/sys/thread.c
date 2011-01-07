/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "thread.h"
#include "common.h"
/**
   Functions regarding to threading
*/

/**
   Break out the job to be executed by multiple threads. 
   if interlaced==0: partition the job to consecutive segments.
   if interlaced!=0:  different thread do the job interelaced.
*/
void thread_prep(thread_t *info, long start, long tot, long nthread, 
		 thread_wrapfun fun, void *data){
    if(tot==0 && nthread>1){
	error("Need to specify tot\n");
    }
    long nt=(tot-start)/nthread;
    long ithread;
    for(ithread=0; ithread<nthread; ithread++){
	info[ithread].ithread=ithread;
	info[ithread].data=data;
	info[ithread].fun=fun;
	info[ithread].start=start;
	start+=nt;
	info[ithread].end=start;
	if(start>=tot){
	    info[ithread].end=tot;
	    ithread++;
	    break;
	}
    }
    for(;ithread<nthread;ithread++){//skip these threads.
	info[ithread].ithread=ithread;
	info[ithread].start=0;
	info[ithread].end=0;
	info[ithread].data=data;
	info[ithread].fun=fun;
    }
    if(info[nthread-1].end){
	info[nthread-1].end=tot;
    }
}
static pthread_mutex_t atomic_lock=PTHREAD_MUTEX_INITIALIZER;
int assign_increment(int *src, int step){
    int result;
    pthread_mutex_lock(&atomic_lock);
    result=*src;
    *src+=step;
    pthread_mutex_unlock(&atomic_lock);
    return result;
}
