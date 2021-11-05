/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#include <signal.h>
#include <unistd.h>
#include <search.h>
#ifndef __CYGWIN__
#include <execinfo.h>
#endif
#include <sys/stat.h>
#include <errno.h>
#include <dlfcn.h>
int signal_caught=0;//indicate that signal is caught.
#define IN_MEM_C 1
#include "mem.h"
#include "thread.h"
#include "scheduler_client.h"
#include "process.h"
#include "misc.h"
/*
  Record allocated memory and their nelem.  Warn if some memory is not freed.
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

__thread char funtrace[funtrace_len]={0};
#undef strdup
#undef strndup
void *(*calloc_default)(size_t, size_t)=NULL;
void *(*malloc_default)(size_t)=NULL;
void *(*realloc_default)(void *, size_t)=NULL;
void  (*free_default)(void *)=NULL;

/*void* (*calloc_custom)(size_t, size_t);
void* (*malloc_custom)(size_t);
void* (*realloc_custom)(void*, size_t);
void  (*free_custom)(void*);
*/
int MEM_VERBOSE=0;
int MEM_DEBUG=DEBUG+0;
int LOG_LEVEL=0;
PNEW(mutex_mem);
PNEW(mutex_deinit);

static long  memcnt=0;
static size_t memalloc=0, memfree=0;
static void *MROOT=NULL;
static void *MSTATROOT=NULL;//reuses keys from MROOT.
/*max depth in backtrace */
#define DT 16
typedef struct{
	void *p;
	char *funtrace;//new way of recording initial call location (not yet used)
	void *func[DT];
	int nfunc;
	size_t nelem;
	size_t nbyte;
}T_MEMKEY;

typedef int(*compar)(const void *, const void *);
static int stat_cmp(const T_MEMKEY *pa, const T_MEMKEY *pb){
	if(!pa->nfunc&&!pb->nfunc){
		if(pa->funtrace&&pb->funtrace){
			return strcmp(pa->funtrace, pb->funtrace);
		} else{
			dbg("nfunc not set but also funtrace not set\n");
			return 1;
		}
	}
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
static void stat_usage(const void *key, VISIT which, int level){
	const T_MEMKEY *key2=*((const T_MEMKEY **)key);
	(void)level;
	//merge calls from the same location
	if(which==leaf||which==postorder){
		void *found;
		if(!(found=tfind(key2, &MSTATROOT, (compar)stat_cmp))){//not found
			if(!tsearch(key2, &MSTATROOT, (compar)stat_cmp)){
				warning("Error inserting to tree\n");
			}
		} else{
			(*(T_MEMKEY **)found)->nelem+=key2->nelem;
		}
	}
}
static int print_count=0;
static void print_usage(const void *key, VISIT which, int level){
	const T_MEMKEY *key2=*((const T_MEMKEY **)key);
	(void)level;
	if(which==leaf||which==postorder){
		print_count++;
		if(print_count<100){
			info3("%9zu(%4zu)", key2->nelem, key2->nbyte);
			if(key2->funtrace[0]){
				info3(" %s\n", key2->funtrace);
			}else if(key2->nfunc){
				int offset=key2->nfunc>3?1:0;
				if(print_backtrace_symbol(key2->func, key2->nfunc-offset)){
					info3("print_backtrace_symbol failed, stop.\n");
					print_count=101;
				}
			} 
		} else if(print_count==100){
			info3("Stop after %d records\n", print_count);
		}
	}
}
typedef struct deinit_t{/*contains either fun or data that need to be freed. */
	void (*fun)(void);
	void *data;
	struct deinit_t *next;
}deinit_t;
deinit_t *deinit_head=NULL;

static int key_cmp(const void *a, const void *b){
	void *p1, *p2;
	if(!a||!b) return 1;
	p1=((T_MEMKEY *)a)->p;
	p2=((T_MEMKEY *)b)->p;
	if(p1<p2) return -1;
	else if(p1>p2) return 1;
	else return 0;
}
static void memkey_add(void *p, size_t nbyte, size_t nelem){
	if(!p){
		if(nelem){
			warning("memory allocation for %zu failed\n", nelem);
		}
		return;
	}

	T_MEMKEY *key=(T_MEMKEY *)calloc_default(1, sizeof(T_MEMKEY));
	key->p=p;
	key->nbyte=nbyte;
	key->nelem=nelem;
	if(funtrace[0]){//do not use backtrace
		key->funtrace=strdup(funtrace);
	} else{
#ifndef __CYGWIN__
		dbg("Using backtrace\n");
		key->nfunc=backtrace(key->func, DT);//this step requires locking and severely limits multi-thread performance
#endif
	}
	LOCK(mutex_mem);
	if(tfind(key, &MROOT, key_cmp)){
		//print_backtrace();
		warning("%p already exists\n", p);
		free(key);
	} else{
		if(!tsearch(key, &MROOT, key_cmp)){
			warning("Error inserting to tree\n");
		} else{
			memcnt++;
			memalloc+=nelem*nbyte;
		}
	}
	UNLOCK(mutex_mem);
}
//Return 0 if success. 1 if failed.
static int memkey_del(void *p){
	if(!p) return 0;
	int ans=1;
	void *found=0;
	T_MEMKEY key;
	key.p=p;
	LOCK(mutex_mem);
	found=tfind(&key, &MROOT, key_cmp);
	if(found){
		T_MEMKEY *key1=*(T_MEMKEY **)found;/*the address of allocated T_MEMKEY. */
		size_t nelem=key1->nelem;
		size_t nbyte=key1->nbyte;
		memfree+=nelem*nbyte;
		memcnt--;
		if(!tdelete(&key, &MROOT, key_cmp)){/*return parent. */
			warning("%p free: Record deleting error\n", p);
		}
		UNLOCK(mutex_mem);
		free_default(key1);
		ans=0;
	} else{
		UNLOCK(mutex_mem);
		warning("%p free: Record not found\n", p);
		ans=1;
	}
	return ans;
}
static void memkey_update(void *p, size_t nbyte, size_t nelem){
	T_MEMKEY key;
	key.p=p;
	LOCK(mutex_mem);
	void* found=tfind(&key, &MROOT, key_cmp);
	if(found){
		T_MEMKEY *key1=*(T_MEMKEY **)found;
		memalloc+=(nbyte*nelem-key1->nbyte*key1->nelem);
		key1->nbyte=nbyte;
		key1->nelem=nelem;
	}else{
		warning("%p realloc: Record not found\n", p);
	}
	UNLOCK(mutex_mem);
}
void *calloc_maos(size_t nbyte, size_t nelem){
	void *p=calloc_default(nbyte, nelem);
	if(MEM_DEBUG){
		memkey_add(p, nbyte, nelem);
	}
	if(MEM_VERBOSE){
		info("%s calloc: %p with %zu (%2zu) bytes\n", funtrace, p, nelem, nbyte);
	}
	funtrace_unset;
	return p;
}
void *malloc_maos(size_t size){
	void *p=malloc_default(size);
	if(MEM_DEBUG){
		memkey_add(p, 1, size);
	}
	if(MEM_VERBOSE){
		info("%s malloc: %p with %zu bytes\n", funtrace, p, size);
	}
	funtrace_unset;
	return p;
}
void *realloc_maos(void *p0, size_t size){
	void *p=realloc_default(p0, size);
	if(MEM_DEBUG){
		if(p!=p0){
			if(p0) memkey_del(p0);
			memkey_add(p, 1, size);
		}else{
			memkey_update(p, 1, size);
		}
	}
	if(MEM_VERBOSE){
		info("%s realloc: %p with %zu bytes\n", funtrace, p, size);
	}
	funtrace_unset;
	return p;
}
void free_maos(void *p){
	if(MEM_DEBUG){
		memkey_del(p);
	}
	free_default(p);//must be after memkey_del to avoid race condition (add before del)
	if(MEM_VERBOSE){
		info("%s free: %p\n", funtrace, p);
	}
	funtrace_unset;
}

static void print_mem_debug(){
	if(MROOT){
		info3("%ld (%lu MB) allocated memory not freed!!!\n", memcnt, (memalloc-memfree)>>20);
		if(!signal_caught){
			print_count=0;
			twalk(MROOT, stat_usage);//walk over the recording tree and combine records with the same backtrace
			twalk(MSTATROOT, print_usage);//print results.
		} else{
			info3("Printing of not freed memory is enabled only when signal_caught=0.\n");
		}
	} else{
		info3("All allocated memory are freed.\n");
		if(memcnt>0){
			info3("Memory counter is still not zero: %ld\n", memcnt);
		}
	}
	info3("Memory allocation is %lu kB, still in use %lu KB.\n", memalloc>>10, (memalloc-memfree)>>10);
}
void read_sys_env(){
	READ_ENV_INT(MEM_DEBUG, 0, 1);
	READ_ENV_INT(MEM_VERBOSE, 0, 2);
	READ_ENV_INT(LOG_LEVEL, -5, 5);
}
FILE *fpconsole=NULL;
int err2out=1;
int std2out=1;
static void init_mem(){
	if(!calloc_default){
		calloc_default=(void *(*)(size_t, size_t))dlsym(RTLD_DEFAULT, "calloc");
		malloc_default=(void *(*)(size_t))dlsym(RTLD_DEFAULT, "malloc");
		realloc_default=(void *(*)(void *, size_t))dlsym(RTLD_DEFAULT, "realloc");
		free_default=(void(*)(void *))dlsym(RTLD_DEFAULT, "free");
		read_sys_env();
		void init_process(void);
		init_process();
		init_hosts();
		fpconsole=fdopen(dup(fileno(stdout)), "a");
	}
}
static __attribute__((constructor)) void init(){
	init_mem();
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
	for(deinit_t *p1=deinit_head;p1;p1=deinit_head){
		deinit_head=p1->next;
		if(p1->fun) p1->fun();
		if(p1->data) myfree(p1->data);
		free_default(p1);
	}
	if(MEM_DEBUG){
		print_mem_debug();
	}
}

/**
   Register a function or data to call or free upon exit
*/
void register_deinit(void (*fun)(void), void *data){
	if(MEM_DEBUG){
		LOCK(mutex_deinit);
		init_mem();
		deinit_t *node=NULL;
		for(node=deinit_head; node; node=node->next){
			if(fun&&fun==node->fun){
				fun=NULL;
			}
			if(data&&data==node->data){
				data=NULL;
			}
		}
		if(fun||data){
			node=(deinit_t *)calloc_default(1, sizeof(deinit_t));
			node->fun=fun;
			node->data=data;

			node->next=deinit_head;
			deinit_head=node;
		}
		UNLOCK(mutex_deinit);
	}
}

quitfun_t quitfun=NULL;
void default_quitfun(const char *msg){
	info("%s", msg);
	//sync();
	if(strncmp(msg, "FATAL", 5)){
		print_backtrace();
		//sync();
		raise(1);
	}
	exit(0);
}
static int (*signal_handler)(int)=0;
static volatile sig_atomic_t fatal_error_in_progress=0;
void default_signal_handler(int sig, siginfo_t *siginfo, void *unused){
	(void)unused;
	signal_caught=sig;
	info("\ndefault_signal_handler: %s (%d)\n", strsignal(sig), sig);
	if(sig==SIGTERM){
		char sender[PATH_MAX]={0};
		get_job_progname(sender, PATH_MAX, siginfo->si_pid);
		info("Code is %d, send by %d (uid=%d, %s).\n",
			siginfo->si_code, siginfo->si_pid, siginfo->si_uid, sender);
	}
	//sync();
	int cancel_action=0;
	/*
	{
	struct sigaction act={0};
	act.sa_handler=SIG_DFL;
	sigaction(sig, &act, 0);//prevent recursive call of handler
	sync();
	}
	*/
	if(sig==0){
		dbg_time("Signal 0 caught. do nothing\n");
		return;
	}
	if(fatal_error_in_progress){
		info("Signal handler is already in progress. force abort.\n");
		_Exit(1);
		return;
	}
	fatal_error_in_progress++;
	if(iscrash(sig)){
		if(siginfo&&siginfo->si_addr){
			info("Memory location: %p\n", siginfo->si_addr);
		}
		//It is not safe to call backtrace in SIGSEGV, so may hang.
		print_backtrace();
	}
	if(signal_handler){
		dbg_time("Call signal_handler %p\n", signal_handler);
	}
	if(signal_handler&&signal_handler(sig)){
		cancel_action=1;
	}
	//sync();
	if(!cancel_action){//Propagate signal to default handler.
		struct sigaction act={0};
		act.sa_handler=SIG_DFL;
		sigaction(sig, &act, 0);
		raise(sig);
	}/*else{//cancel signal, keep going
		struct sigaction act={0};
		act.sa_handler=NULL;
		act.sa_sigaction=default_signal_handler;
		act.sa_flags=SA_SIGINFO;
		sigaction(sig, &act, 0);
	}*/
}

/**
   Register signal handler
*/
void register_signal_handler(int (*func)(int)){
	struct sigaction act={0};
	act.sa_sigaction=default_signal_handler;
	act.sa_flags=SA_SIGINFO;
	sigaction(SIGBUS, &act, 0);
	sigaction(SIGILL, &act, 0);
	sigaction(SIGSEGV, &act, 0);
	sigaction(SIGINT, &act, 0);
	sigaction(SIGTERM, &act, 0);
	sigaction(SIGABRT, &act, 0);
	//sigaction(SIGHUP, &act, 0);
	signal(SIGHUP, SIG_IGN);//ignore the signal.
	sigaction(SIGUSR1, &act, 0);
	sigaction(SIGQUIT, &act, 0);
	signal_handler=func;
}