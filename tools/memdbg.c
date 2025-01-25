/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#if defined(__linux__)
#define _GNU_SOURCE


#include <stddef.h>
#include <pthread.h>

#include <search.h>
#include <dlfcn.h>

#include <limits.h>
#include <sys/syscall.h>
/**
   Notice: When destructor or at exit is called on this file, using
   malloc/calloc etc will cause this file to be reloaded (constructor called)
   and therefore stuck in infinite loop.
 */
static volatile int disable_memdbg=0;
static int memdbg_toreport=0;
static void* (*malloc_default)(size_t size)=0;
static void* (*calloc_default)(size_t nmemb, size_t size)=0;
static void (*free_default)(void* p)=0;
static void* (*realloc_default)(void* p, size_t size)=0;
static void* (*memalign_default)(size_t blocksize, size_t bytes)=0;
static void* (*mmap_default)(void* addr, size_t len, int prot, int flag, int fd, off_t offset)=0;
static void (*munmap_default)(void* addr, size_t len)=0;

/*
  Derived from mem.c, but intended for using using LD_PRELOAD.

  Record allocated memory and their size.  Warn if some memory is not freed.
  Notice that the pointer return by tsearch and the pointer passwd to the
  function called by twalk is the addess of the tree node, which points to the
  key. so **key is required.

  We use the backtrace to see the chain of calling when memory is
  allocated. This is better than using BASEFILE, __LINE__ which does not tell
  where the call from eventually from ans is not useful.
*/
#define MIN(a,b) (((a)<(b))?(a):(b))
#define dbg(A...) fprintf(stderr, A)
#define error(A...) ({ fprintf(stderr, "\033[01;31mFatal error\033[00;00m\t" A); exit(1);})
#define warning(A...) fprintf(stderr, "\033[01;31mWarning\033[00;00m\t" A)
#define LOCK(A) pthread_mutex_lock(&A)
#define UNLOCK(A) pthread_mutex_unlock(&A)
/*
#define LOCK(A) ({int tid=syscall(SYS_gettid); dbg("%d locking...%p, ", tid, &A); \
		pthread_mutex_lock(&A);					\
		dbg("%d locked %p\n", tid, &A);})
#define UNLOCK(A) ({int tid=syscall(SYS_gettid); dbg("%d unlocking...%p, ", tid, &A);\
		pthread_mutex_unlock(&A); \
		dbg("%d unlocked %p\n", tid, &A);})
*/
#define PNEW(A) static pthread_mutex_t A=PTHREAD_MUTEX_INITIALIZER
#define READ_ENV_INT(A,min,max)				\
    if(getenv("MAOS_"#A)){				\
	A=strtol(getenv("MAOS_"#A),NULL,10);		\
	dbg(#A"=%d\n", A);				\
	if(A>max || A<min){				\
	    error("MAOS_%s: invalid range\n", #A);	\
	}						\
    }
PNEW(mutex_mem);
void memdbg_report(void);
void memdbg_enable(int flag){
	if(flag){
		dbg("memdbg enabled\n");
		disable_memdbg=0;
	} else{
		dbg("memdbg disabled\n");
		disable_memdbg=1;
	}
}
char* get_job_progname(int pid){
	char* progname=NULL;
	char path[PATH_MAX];
	/*readlink doesn't append \0 */
	char procpath[PATH_MAX];
	snprintf(procpath, PATH_MAX, "/proc/%d/exe", pid>0?pid:(int)getpid());
	int nprog=readlink(procpath, path, PATH_MAX);
	if(nprog>0){
		path[nprog]='\0';
		char* delt=strstr(path, "(deleted)");
		if(delt){
			delt[0]='\0';
		}
		progname=strdup(path);
	} else{
		warning("Failed to readlink %s\n", procpath);
	}
	return progname;
}

char* call_addr2line(const char* buf){
	FILE* fpcmd=popen(buf, "r");
	char out[4096]; strcpy(out, "Backtrace:");
	if(!fpcmd){
		strcpy(out, "Command failed\n");
	} else{
		char ans[4096];
		while(fgets(ans, sizeof(ans), fpcmd)){
			char* tmp=strrchr(ans, '/');
			if(tmp){
				tmp++;
				char* tmp2=strchr(tmp, '\n'); tmp2[0]='\0';
				strcat(out, "->");
				strcat(out, tmp);
				if(strlen(out)>=4000) break;
			}
		}
		pclose(fpcmd);
	}
	return strdup(out);
}
/**
   Convert backtrace address to source line.
 */
void print_backtrace_symbol(void* const* buffer, int size){
	LOCK(mutex_mem);
	int disable_save=disable_memdbg;
	disable_memdbg=1;
	char add[24];
	char* progname=get_job_progname(0);/*don't free pointer. */
	if(!progname){
		warning("Unable to get progname\n");
		goto end;
	}
	char cmdstr[4096];
	strcpy(cmdstr, "addr2line -f -e ");
	strcat(cmdstr, progname);
	free(progname);
	for(int it=size-1; it>0; it--){
		snprintf(add, 24, " %p", buffer[it]);
		strcat(cmdstr, add);
	}
	char* ans=call_addr2line(cmdstr);
	dbg("%s\n", ans);
	free(ans);
	sync();
end:
	disable_memdbg=disable_save;
	UNLOCK(mutex_mem);
}
#include <execinfo.h>
void print_backtrace(){
	int size0, size1;
	size0=1024;
	void* buffer[size0];
	size1=backtrace(buffer, size0);
	print_backtrace_symbol(buffer, size1);
}

static int MEM_VERBOSE=0;
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
}T_MEMKEY;
typedef struct{
	void* p;
	void** func;
	int nfunc;
	size_t size;
}T_STATKEY;

static int stat_cmp(const void* a, const void* b){
	const T_STATKEY* pa=a;
	const T_STATKEY* pb=b;
	int nfunc=MIN(pa->nfunc, pb->nfunc);
	for(int i=0; i<nfunc; i++){
		if(pa->func[i]<pb->func[i]){
			return -1;
		} else if(pa->func[i]>pb->func[i]){
			return 1;
		}
	}
	if(pa->nfunc<pb->nfunc){
		return -1;
	} else if(pa->nfunc>pb->nfunc){
		return 1;
	} else{
		return 0;
	}
}
static void stat_usage(const void* key, VISIT which, int level){
	const T_MEMKEY* key2=*((const T_MEMKEY**)key);
	(void)level;
	if(which==leaf||which==postorder){
		void** found;
		T_STATKEY key3;
		key3.func=(void**)key2->func;
		key3.nfunc=key2->nfunc-2;
		if(!(found=tfind(&key3, &MSTATROOT, stat_cmp))){
			T_STATKEY* keynew=calloc_default(1, sizeof(T_STATKEY));
			P(keynew)=P(key2);
			keynew->func=malloc_default(key3.nfunc*sizeof(void*));
			memcpy(keynew->func, key3.func, sizeof(void*)*key3.nfunc);
			keynew->nfunc=key3.nfunc;
			keynew->size=key2->size;
			if(!tsearch(keynew, &MSTATROOT, stat_cmp)){
				error("Error inserting to tree\n");
			}
		} else{
			T_STATKEY* keynew=*found;
			keynew->size+=key2->size;
		}
	}
}
static void print_usage(const void* key, VISIT which, int level){
	const T_STATKEY* key2=*((const T_STATKEY**)key);
	(void)level;
	if(which==leaf||which==postorder){
		dbg("size %4zu B@%p\n", (key2->size), P(key2));
		print_backtrace_symbol(key2->func, key2->nfunc);
	}
}
typedef struct deinit_t{/*contains either fun or data that need to be freed. */
	void (*fun)(void);
	void* data;
	struct deinit_t* next;
}deinit_t;
deinit_t* deinit_head=NULL;

static int key_cmp(const void* a, const void* b){
	void* p1, * p2;
	if(!a||!b) return 1;
	p1=((T_MEMKEY*)a)->p;
	p2=((T_MEMKEY*)b)->p;
	if(p1<p2) return -1;
	else if(p1>p2) return 1;
	else return 0;
}
static void memkey_add(void* p, size_t size){
	if(!p||disable_memdbg) return;
	LOCK(mutex_mem);
	int disable_save=disable_memdbg;
	disable_memdbg=1;//prevent recursive calling
	T_MEMKEY* key=calloc_default(1, sizeof(T_MEMKEY));
	P(key)=p;
	key->size=size;
	key->nfunc=backtrace(key->func, DT);
	void** found=tfind(key, &MROOT, key_cmp);
	if(found){
		T_MEMKEY* key1=*found;
		warning("memkey_add: %p already exists with size %zd. new size %zd\n",
			p, key1->size, size);
		key1->size=size;
	} else{
		if(!tsearch(key, &MROOT, key_cmp)){
			warning("memkey_add: Error inserting to tree\n");
		}
		memcnt++;
		memalloc+=size;
	}
	disable_memdbg=disable_save;
	if(MEM_VERBOSE==1){
		dbg("%p malloced with %zu bytes: %s\n", p, size, found?"collision":"success");
	} else if(MEM_VERBOSE==2&&size>1024){
		dbg("Alloc:%.3f MB mem used\n", (memalloc-memfree)/1024./1024.);
	}
	UNLOCK(mutex_mem);
	if(found) print_backtrace();
}

static void memkey_del(void* p){
	if(!p||disable_memdbg) return;
	void** found=0;
	T_MEMKEY key;
	LOCK(mutex_mem);
	key.p=p;
	int disable_save=disable_memdbg;
	disable_memdbg=1;
	found=tfind(&key, &MROOT, key_cmp);
	if(found){
		T_MEMKEY* key1=*found;/*the address of allocated T_MEMKEY. */
		if(MEM_VERBOSE==1){
			dbg("%p freed with %zu bytes: success\n", p, key1->size);
		} else if(MEM_VERBOSE==2&&key1->size>1024){
			dbg("Free: %.3f MB mem used\n", (memalloc-memfree)/1024./1024.);
		}
		memfree+=key1->size;
		memcnt--;
		if(!tdelete(&key, &MROOT, key_cmp)){/*return parent. */
			warning("memkey_del: Error deleting old record\n");
		}
		free_default(key1);
	}
	disable_memdbg=disable_save;
	UNLOCK(mutex_mem);
	if(!found){
		warning("%p not found\n", p);
		print_backtrace();
	}
}

static int inited=0;
static void get_default(){
	inited=1;//Mark as false so if dlsym calls malloc/calloc won't cycle
	LOCK(mutex_mem);
	dbg("Overridding default malloc with customized one.\n");
	malloc_default=dlsym(RTLD_NEXT, "malloc");
	calloc_default=dlsym(RTLD_NEXT, "calloc");
	realloc_default=dlsym(RTLD_NEXT, "realloc");
	free_default=dlsym(RTLD_NEXT, "free");
	memalign_default=dlsym(RTLD_NEXT, "memalign");
	mmap_default=dlsym(RTLD_NEXT, "mmap");
	munmap_default=dlsym(RTLD_NEXT, "munmap");
	UNLOCK(mutex_mem);
}
static void* simple_calloc(size_t nmemb, size_t size){
	PNEW(mutex_local);
	LOCK(mutex_local);
	static char buffer[4096];
	static char* bufoff=0;
	if(!bufoff){
		bufoff=buffer;
		memset(buffer, 0, sizeof(buffer));
	}
	void* res=0;
	if(size*nmemb<sizeof(buffer)+(buffer-bufoff)){
		res=bufoff;
		bufoff+=size*nmemb;
	} else{
		error("Too much memory requested. Increase buffer size\n");
	}
	UNLOCK(mutex_local);
	return res;
}
void* malloc(size_t size){
	if(!malloc_default){//Need initializing
		if(!inited){//True first time
			get_default();
		} else{//simple allocator
			return simple_calloc(size, 1);
		}
	}
	void* p=malloc_default(size);
	if(!disable_memdbg) memkey_add(p, size);
	return p;
}

void* calloc(size_t nmemb, size_t size){
	if(!calloc_default){//Need initializing
		if(!inited){//True first time
			get_default();
		} else{//simple allocator
			return simple_calloc(nmemb, size);
		}
	}
	void* p=calloc_default(nmemb, size);
	if(!disable_memdbg) memkey_add(p, size);
	return p;
}

void* realloc(void* p0, size_t size){
	if(!disable_memdbg&&p0) memkey_del(p0);
	void* p=realloc_default(p0, size);
	if(!disable_memdbg) memkey_add(p, size);
	return p;
}

void free(void* p){
	if(p&&!disable_memdbg) memkey_del(p);
	free_default(p);
}

void* memalign(size_t blocksize, size_t bytes){
	if(!blocksize){
		if(bytes==(size_t)-1){
			memdbg_toreport=1;
			memdbg_enable(1);
		} else if(bytes==(size_t)-2){
			memdbg_enable(0);
			memdbg_report();
		} else{
			error("Invalid usage of free()\n");
		}
		return 0;
	}
	void* p=memalign_default(blocksize, bytes);
	if(!disable_memdbg) memkey_add(p, bytes);
	return p;
}

void* mmap(void* addr, size_t len, int prot, int flag, int fd, off_t offset){
	void* p=mmap_default(addr, len, prot, flag, fd, offset);
	if(!disable_memdbg) memkey_add(p, len);
	return p;
}

void munmap(void* addr, size_t len){
	if(!disable_memdbg) memkey_del(addr);
	munmap_default(addr, len);
}

/**
   Register routines to be called with mem.c is unloading (deinit).
 */
void memdbg_report(void){
	memdbg_toreport=0;
	int disable_save=disable_memdbg;
	disable_memdbg=1;
	if(MROOT){
		dbg("%ld (%.3f B) allocated memory not freed!!!\n",
			memcnt, (memalloc-memfree));
		twalk(MROOT, stat_usage);
		dbg("stat_usage is done\n");
		twalk(MSTATROOT, print_usage);
	} else{
		dbg("All allocated memory are freed.\n");
		if(memcnt>0){
			dbg("But memory count is still none zero: %ld\n", memcnt);
		}
	}
	dbg("Total allocated memory is %.3f MB\n", memalloc/1024./1024.);
	dbg("Total freed     memory is %.3f MB\n", memfree/1024./1024.);
	disable_memdbg=disable_save;
}
__attribute__((constructor)) void init(){
	unsetenv("LD_PRELOAD");//prevent ld_preload on child process.
	dbg("constructor called\n");
	/*
	  Don't assign malloc in constructor because other constructor may execute
	  earlier than this constructor and already calls malloc. This is confirmed
	  because get_default printed message earlier than init().
	*/
}
__attribute__((destructor)) void deinit(){
	if(memdbg_toreport){
		memdbg_report();
	}
}
#endif
