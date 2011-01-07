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

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
int exit_success=0;
void (*call_freepath)(void)=NULL;

#define MEM_VERBOSE 0
#include "mem.h"
#include "thread.h"
#if USE_MEM == 1
#include "scheduler_client.h"
#include "../path.h"
#include "../fractal.h"
#include <math.h>
#include <search.h>
#include <string.h>
#include <execinfo.h>
#include <sys/stat.h>
#include <unistd.h>
/*
  Record allocated memory and their size.  Warn if some memory is not freed.
  Notice that the pointer return by tsearch and the pointer passwd to the
  function called by twalk is the addess of the tree node, which points to the
  key. so **key is required.

  don't use fdopen to open the pipes! it will block the routine. 
  use read/write instead

  read/write is very good.  when one process finish writing, they other will read
  the same thing. even binary data can be handled.

*/


#include "misc.h"
#undef malloc
#undef calloc
#undef free
#undef realloc
PNEW(mutex_mem);
#if MEM_VERBOSE == 1
#define meminfo(A...) {fprintf(stderr,A);}
#else
#define meminfo(A...)
#endif

static void *MROOT=NULL;
static long long memcnt=0;

//max depth in backtrace
#define DT 12
typedef struct T_MEMKEY{
    void *p;
    void *func[DT];
    int nfunc;
    size_t size;
}T_MEMKEY;

static void print_usage(const void *key, VISIT value, int level){
    const T_MEMKEY* key2=*((const T_MEMKEY**)key);
    (void) level;
    if(value>1){
	fprintf(stderr,"%p: size %8zu KiB", key2->p, 
		(key2->size)/1024);
	print_backtrace_symbol(key2->func, key2->nfunc-2);
    }
}
static __attribute__((constructor)) void init(){
    warning2("Memory management is in use\n");
}
static __attribute__((destructor)) void deinit(){
    thread_pool_destroy();
    fractal_vkcov_free();
    freepath();
    if(exit_success){
	if(MROOT){
	    warning("%lld allocated memory not freed!!!\n",memcnt);
	    if(memcnt<100000){
		twalk(MROOT,print_usage);
		warning("%lld allocated memory not freed!!!\n",memcnt);
	    }else{
		warning("Too many, will not print backtrace.\n");
	    }
	}else{
	    info("All allocated memory are freed.\n");
	    if(memcnt>0){
		warning("But memory count is still none zero: %lld\n",memcnt);
	    }
	}
    }else{
	info("exit_success=%d\n", exit_success);
    }
}
static int key_cmp(const void *a, const void *b){
    void *p1,*p2;
    if(!a || !b) return 1; 
    p1=((T_MEMKEY*)a)->p;
    p2=((T_MEMKEY*)b)->p;
    if(p1<p2) return -1;
    else if(p1>p2) return 1;
    else return 0;
}
static void memkey_add(void *p,size_t size){
    if(!p) return;
    LOCK(mutex_mem);
    T_MEMKEY *key=calloc(1,sizeof(T_MEMKEY));
    key->p=p;
    key->size=size;
    key->nfunc=backtrace(key->func,DT);
    if(!tsearch(key, &MROOT, key_cmp)){
	error("Error inserting to tree\n");
    }
    memcnt++;
    UNLOCK(mutex_mem);
}
static void memkey_update(void*p, size_t size){
    if(!p) return;
    LOCK(mutex_mem);
    T_MEMKEY key;
    key.p=p;
    void **found;
    if(!(found=tfind(&key, &MROOT, key_cmp))){
	error("Record not found for memkey_update %p\n",p);
    }
    T_MEMKEY*key0=*found;
    key0->size=size;
    key0->nfunc=backtrace(key0->func,DT);
    meminfo("%p realloced with %lu bytes\n", p, size);
    UNLOCK(mutex_mem);
}
static void memkey_del(void*p){
    if(!p) return;
    LOCK(mutex_mem);
    void **found;
    T_MEMKEY key;
    key.p=p;
    if((found=tfind(&key, &MROOT, key_cmp))){//found.
	T_MEMKEY* key1=*found;//the address of allocated T_MEMKEY.
	if(!tdelete(&key, &MROOT, key_cmp)){//return parent.
	    error("Error deleting old record\n");
	}
	free(key1);
	memcnt--;
    }
    UNLOCK(mutex_mem);
}
void *CALLOC(size_t nmemb, size_t size){
    void *p=calloc(nmemb,size);
    if(p){
	memkey_add(p,size);
    }else if(size!=0){
	error("calloc for %zu, %zu failed\n",nmemb, size);
    }
    meminfo("%p calloced with %zu bytes\n",p,nmemb*size);
    return p;
}
void *MALLOC(size_t size){
    void *p=malloc(size);
    if(p){
	memkey_add(p,size);
    }else if(size!=0){
	error("malloc for %zu failed\n",size);
    }
    meminfo("%p malloced with %lu bytes\n",p, size);
    return p;
}
void *REALLOC(void*p0, size_t size){
    if(!p0) return MALLOC(size);
    void *p=realloc(p0,size);
    if(p){
	//return p;
	if(p==p0){
	    memkey_update(p,size);
	}else{
	    memkey_del(p0);
	    memkey_add(p,size);
	}
    }else{
	memkey_del(p0);
	if(size!=0){
	    error("realloc for size %zu failed\n",size);
	}
    }
    return p;
}
void FREE(void *p){
    if(!p) return;
    memkey_del(p);
    free(p);
}
#endif

#if USE_MEM == 1
void mem_usage(void){
    info2("There are %lld allocated memory not freed!!!\n",memcnt);
}
#endif
#if USE_MEM == 1
size_t memsize(void *p){
    size_t tot=0;
    if(p){
	T_MEMKEY key,**key2;
	key.p=p;
	key2=tfind(&key,&MROOT,key_cmp);
	if(key2){
	    tot=(*key2)->size;
	}else{
	    tot=0;
	}
    }
    return tot;
}
#endif


#ifdef TEST
#include <signal.h>
int main(){
    int N=1024;
    double *a,*b,*c;
    signal(SIGSEGV, print_backtrace);
    a=MALLOC(sizeof(double)*N);
    b=MALLOC(sizeof(double)*N);
    b=REALLOC(b,sizeof(double)*N/2);
    c=CALLOC(N, sizeof(double));
    FREE(a);
    FREE(b);
    FREE(c);
}
#endif
