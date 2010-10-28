/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#elif USE_PTHREAD==2 //Use pthread pool. not always available
#include "thr_pool.h"
#define CALL(A,B,nthread)				\
    if(nthread>1){					\
	for(int ithread=0; ithread<nthread; ithread++){	\
	    thr_pool_queue(default_pool,		\
			   (void*(*)(void*))A,		\
			   (void*)B);			\
	}						\
	thr_pool_wait(default_pool);			\
    }else{						\
	A(B);						\
    }
#define CALL_EACH(A,B,nthread)				\
    if(nthread>1){					\
	for(int ithread=0; ithread<nthread; ithread++){	\
	    thr_pool_queue(default_pool,		\
			   (void*(*)(void*))A,		\
			   (void*)&(B[ithread]));	\
	}						\
	thr_pool_wait(default_pool);			\
    }else{						\
	A(B);						\
    }
#else
#define CALL(A,B,nthread) A(B) //no threading
#define CALL_DATAEACH(A,B,nthread) A(B)
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

typedef struct thread_t{
    long start;
    long end;
    long step;
    long ithread;//which thread this is.
    void *data;
}thread_t;
void thread_prep(thread_t *info, long start, long tot, long interlaced, long nthread, void *data);
#endif
