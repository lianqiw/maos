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
#ifndef AOS_LIB_THREAD_H
#define AOS_LIB_THREAD_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef USE_PTHREAD
#if MATLAB_MEX_FILE
#define USE_PTHREAD 1
#else
#define USE_PTHREAD 0
#endif
#endif
#if USE_PTHREAD
#include <pthread.h>
#endif
#if USE_PTHREAD==1 //Conventional pthreading
#define CALL(A,B,nthread)				\
    if(nthread>1){					\
	pthread_t threads[nthread];			\
	for(int ithread=0; ithread<nthread; ithread++){	\
	    pthread_create(&threads[ithread], NULL,	\
			   (void*(*)(void*))A,		\
			   (void*)B);			\
	}						\
	for(int ithread=0; ithread<nthread; ithread++){	\
	    pthread_join(threads[ithread], NULL);	\
	}						\
    }else{						\
	A(B);						\
    }
#define CALL_EACH(A,B,nthread)					\
    if(nthread>1){						\
	pthread_t threads[nthread];				\
	for(int ithread=0; ithread<nthread; ithread++){		\
	    pthread_create(&threads[ithread], NULL,		\
			   (void*(*)(void*))A,			\
			   (void*)&(B[ithread]));		\
	}							\
	for(int ithread=0; ithread<nthread; ithread++){		\
	    pthread_join(threads[ithread], NULL);		\
	}							\
    }else{							\
	A(B);							\
    }
#define THREAD_POOL_INIT(A)  //Do nothing
#elif USE_PTHREAD==2 //Use pthread pool. 2010-11-30: updated with new, simpler implementation
#include "thread_pool.h"
#define CALL(A,B,nthread)				\
    if(nthread>1){					\
	for(int ithread=0; ithread<nthread; ithread++){	\
	    thread_pool_queue((void*(*)(void*))A,	\
			      (void*)B);		\
	}						\
	thread_pool_wait();				\
    }else{						\
	A(B);						\
    }
#define CALL_EACH(A,B,nthread)				\
    if(nthread>1){					\
	for(int ithread=0; ithread<nthread; ithread++){	\
	    thread_pool_queue((void*(*)(void*))A,	\
			      (void*)&(B[ithread]));	\
	}						\
	thread_pool_wait();				\
    }else{						\
	A(B);						\
    }
#define THREAD_POOL_INIT(A) thread_pool_create(A)
#else
#define CALL(A,B,nthread) A(B) //no threading
#define CALL_DATAEACH(A,B,nthread) A(B)
#define THREAD_POOL_INIT(A)  //Do nothing
#endif

#if USE_PTHREAD > 0
#define LOCK(A) pthread_mutex_lock(&A)
#define UNLOCK(A) pthread_mutex_unlock(&A)
#define PINIT(A) pthread_mutex_init(&A,NULL)
#define PDEINIT(A) pthread_mutex_destroy(&A)
#define PNEW(A) static pthread_mutex_t A=PTHREAD_MUTEX_INITIALIZER
#else
#define LOCK(A)
#define UNLOCK(A)
#define PINIT(A)
#define PDEINIT(A)
#define PNEW(A)
#endif
/**
   Information about job to launch for each thread. start and end are the two indices.
*/
struct thread_t{
    long start;
    long end;
    long step;
    long ithread;//which thread this is.
    void *data;
};
#ifndef AOS_ACCPHI_H
typedef struct thread_t thread_t;
#endif
void thread_prep(thread_t *info, long start, long tot, long interlaced, long nthread, void *data);
#endif
