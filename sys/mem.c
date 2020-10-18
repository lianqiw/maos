/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE //For RTLD_NEXT in linux
#endif
#include <search.h>
#ifndef __CYGWIN__
#include <execinfo.h>
#endif
#include <sys/stat.h>

#include <dlfcn.h>
int exit_fail=0;
#define IN_MEM_C 1
#include "mem.h"
#include "thread.h"
#include "scheduler_client.h"
#include "process.h"
#include "misc.h"
/*
  Record allocated memory and their size.  Warn if some memory is not freed.
  Notice that the pointer return by tsearch and the pointer passwd to the
  function called by twalk is the addess of the tree node, which points to the
  key. so **key is required.

  We use the backtrace to see the chain of calling when memory is
  allocated. This is better than using BASEFILE, __LINE__ which does not tell
  where the call from eventually from ans is not useful.

  2016-06-28: The memory debugging is switched on by environment various instead
  of statically compiled. This has another advantage is that when launched from
  matlab, the memory allocations are replaced by matlab equivalents to have
  matlab manage the memory during referencing..
  */


void* (*calloc_default)(size_t, size_t);
void* (*malloc_default)(size_t);
void* (*realloc_default)(void*, size_t);
void  (*free_default)(void*);

/*void* (*calloc_custom)(size_t, size_t);
void* (*malloc_custom)(size_t);
void* (*realloc_custom)(void*, size_t);
void  (*free_custom)(void*);
*/
static int MEM_VERBOSE=0;
static int MEM_DEBUG=0;
int LOG_LEVEL=0;
PNEW(mutex_mem);
static void* MROOT=NULL;
static long  memcnt=0;
static double memalloc=0, memfree=0;
static void* MSTATROOT=NULL;
/*max depth in backtrace */
#define DT 16
typedef struct{
	void* p;
	void* func[DT];
	int nfunc;
	size_t size;
	size_t nbyte;
}T_MEMKEY;
typedef struct{
	void* p;
	void** func;
	int nfunc;
	size_t size;
	size_t nbyte;
}T_STATKEY;
typedef int(*compar)(const void*, const void*);
static int stat_cmp(const T_STATKEY* pa, const T_STATKEY* pb){
	if(pa->nfunc<pb->nfunc){
		return -1;
	} else if(pa->nfunc>pb->nfunc){
		return 1;
	}
	int nfunc=MIN(pa->nfunc, pb->nfunc);
	for(int i=0; i<nfunc; i++){
		if(pa->func[i]<pb->func[i]){
			return -1;
		} else if(pa->func[i]>pb->func[i]){
			return 1;
		}
	}
	return 0;
}
static void stat_usage(const void* key, VISIT which, int level){
	const T_MEMKEY* key2=*((const T_MEMKEY**)key);
	(void)level;
	if(which==leaf||which==postorder){
	/*printf("%p: size %8zu B", key2->p, (key2->size));
	  print_backtrace_symbol(key2->func, key2->nfunc-2);*/
		void* found;
		T_STATKEY key3;
		key3.func=(void**)key2->func;
		key3.nfunc=key2->nfunc-2;
		if(!(found=tfind(&key3, &MSTATROOT, (compar)stat_cmp))){
			T_STATKEY* keynew=(T_STATKEY*)calloc_default(1, sizeof(T_STATKEY));
			keynew->p=key2->p;
			keynew->func=(void**)malloc_default(key3.nfunc*sizeof(void*));
			memcpy(keynew->func, key3.func, sizeof(void*)*key3.nfunc);
			keynew->nfunc=key3.nfunc;
			keynew->size=key2->size;
			keynew->nbyte=key2->nbyte;
			if(!tsearch(keynew, &MSTATROOT, (compar)stat_cmp)){
				error("Error inserting to tree\n");
			}
		} else{
			T_STATKEY* keynew=*(T_STATKEY**)found;
			keynew->size+=key2->size;
		}
	}
}
static void print_usage(const void* key, VISIT which, int level){
	const T_STATKEY* key2=*((const T_STATKEY**)key);
	(void)level;
	if(which==leaf||which==postorder){
		info("size %4zu(%2zu) B@%p", key2->size, key2->nbyte, key2->p);
		print_backtrace_symbol(key2->func, key2->nfunc);
	}
}
typedef struct T_DEINIT{/*contains either fun or data that need to be freed. */
	void (*fun)(void);
	void* data;
	struct T_DEINIT* next;
}T_DEINIT;
T_DEINIT* DEINIT=NULL;

static int key_cmp(const void* a, const void* b){
	void* p1, * p2;
	if(!a||!b) return 1;
	p1=((T_MEMKEY*)a)->p;
	p2=((T_MEMKEY*)b)->p;
	if(p1<p2) return -1;
	else if(p1>p2) return 1;
	else return 0;
}
static void memkey_add(void* p, size_t nbyte, size_t size){
	if(!p){
		if(size){
			error("memory allocation for %zu failed\n", size);
		}
		return;
	}

	T_MEMKEY* key=(T_MEMKEY*)calloc_default(1, sizeof(T_MEMKEY));
	key->p=p;
	key->nbyte=nbyte;
	key->size=size;
#ifndef __CYGWIN__
	key->nfunc=backtrace(key->func, DT);
#endif
	LOCK(mutex_mem);
	if(tfind(key, &MROOT, key_cmp)){
		warning("%p already exists\n", p);
	}
	if(!tsearch(key, &MROOT, key_cmp)){
		warning("Error inserting to tree\n");
	}
	memcnt++;
	memalloc+=size*nbyte;
	UNLOCK(mutex_mem);
	if(MEM_VERBOSE==1){
		dbg("%p malloced with %zu bytes\n", p, size);
		print_backtrace();
	} else if(MEM_VERBOSE==2&&size>1024){
		info("Alloc:%.3f MB mem used\n", (memalloc-memfree)/1024./1024.);
	}

}
//Return 0 if success. 1 if failed.
static int memkey_del(void* p){
	if(!p) return 0;
	void* found=0;
	T_MEMKEY key;
	key.p=p;
	LOCK(mutex_mem);
	found=tfind(&key, &MROOT, key_cmp);
	if(found){
		T_MEMKEY* key1=*(T_MEMKEY**)found;/*the address of allocated T_MEMKEY. */
		memfree+=key1->size*key1->nbyte;
		memcnt--;
		if(MEM_VERBOSE==1){
			info("Free: %p freed with %zu (%2zu) bytes\n", p, key1->size, key1->nbyte);
		} else if(MEM_VERBOSE==2&&key1->size>1024){
			info("Free: %.3f MB mem used\n", (memalloc-memfree)/1024./1024.);
		}
		if(!tdelete(&key, &MROOT, key_cmp)){/*return parent. */
			warning("Free: Error deleting old record\n");
		}
		UNLOCK(mutex_mem);
		free_default(key1);
		return 0;
	} else{
		UNLOCK(mutex_mem);
		//info("Free: %p not found. Allocated from outside.\n", p);
		//print_backtrace();
		return 1;
	}

}
static void* calloc_dbg(size_t nmemb, size_t size){
	void* p=calloc_default(nmemb, size);
	memkey_add(p, size, nmemb);
	return p;
}
static void* malloc_dbg(size_t size){
	void* p=malloc_default(size);
	memkey_add(p, 1, size);
	return p;
}
static void* realloc_dbg(void* p0, size_t size){
	if(p0) memkey_del(p0);
	void* p=realloc_default(p0, size);
	memkey_add(p, 1, size);
	return p;
}
static void free_dbg(void* p){
	if(!p) return;
	memkey_del(p);
	free_default(p);
}

/**
   Register a function or data to call or free upon exit
*/
void register_deinit(void (*fun)(void), void* data){
	if(MEM_DEBUG){
		T_DEINIT* node=(T_DEINIT*)calloc_default(1, sizeof(T_DEINIT));
		node->fun=fun;
		node->data=data;
		LOCK(mutex_mem);
		node->next=DEINIT;
		DEINIT=node;
		UNLOCK(mutex_mem);
	}
}
#ifdef __cpluspluc
unamespace std{
#endif
	void* malloc_maos(size_t size){
		return MEM_DEBUG?malloc_dbg(size):malloc_default(size);
	}
	void* calloc_maos(size_t nmemb, size_t size){
		return MEM_DEBUG?calloc_dbg(nmemb, size):calloc_default(nmemb, size);
	}
	void* realloc_maos(void* p, size_t size){
		return MEM_DEBUG?realloc_dbg(p, size):realloc_default(p,size);
	}
	void free_maos(void* p){
		if(MEM_DEBUG) free_dbg(p); else free_default(p);
	}
#ifdef __cpluspluc
}
#endif
static void print_mem_debug(){
	if(MROOT){
		warning("%ld (%.3f MB) allocated memory not freed!!!\n",
			memcnt, (memalloc-memfree)/1024./1024.);
		twalk(MROOT, stat_usage);//walk over the recording tree and combine records with the same backtrace
		twalk(MSTATROOT, print_usage);//print results.
	} else{
		info("All allocated memory are freed.\n");
		if(memcnt>0){
			warning("But memory count is still none zero: %ld\n", memcnt);
		}
	}
	info("Total allocated memory is %.3f MB\n", memalloc/1024./1024.);
	info("Total freed     memory is %.3f MB\n", memfree/1024./1024.);
}
void read_sys_env(){
	READ_ENV_INT(MEM_DEBUG, 0, 1);
	READ_ENV_INT(MEM_VERBOSE, 0, 2);
	READ_ENV_INT(LOG_LEVEL, -5, 5);
	/*if(MEM_DEBUG){
		malloc_custom=malloc_dbg;
		realloc_custom=realloc_dbg;
		calloc_custom=calloc_dbg;
		free_custom=free_dbg;
	}*/
}
static __attribute__((constructor)) void init(){
#define RTLD_MINE RTLD_DEFAULT
	calloc_default=(void* (*)(size_t, size_t))dlsym(RTLD_MINE, "calloc");
	malloc_default=(void* (*)(size_t))dlsym(RTLD_MINE, "malloc");
	realloc_default=(void* (*)(void*, size_t))dlsym(RTLD_MINE, "realloc");
	free_default=(void(*)(void*))dlsym(RTLD_MINE, "free");
	read_sys_env();
	void init_process(void);
	init_process();
	init_hosts();
}
/**
   Register routines to be called with mem.c is unloading (deinit).
 */
static __attribute__((destructor)) void deinit(){
	void freepath();
	void thread_pool_destroy();
	//remove files that are 365 days old.
	remove_file_older(CACHE, 1, 30*24*3600);//1 month
	freepath();
	free_hosts();
	thread_pool_destroy();
	for(T_DEINIT* p1=DEINIT;p1;p1=DEINIT){
		DEINIT=p1->next;
		//if(p1->fun) p1->fun();
		//if(p1->data) myfree(p1->data);
		free_default(p1);
	}
	if(MEM_DEBUG){
		if(!exit_fail){
			if(!MEM_DEBUG) return;
			print_mem_debug();
		} else{
			dbg("exit_fail=%d\n", exit_fail);
		}
	}
}
