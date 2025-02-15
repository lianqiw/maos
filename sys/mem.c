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
#include <sys/time.h> //setitimer
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
int MEM_FUNTRACE=0;
int LOG_LEVEL=0;
static void print_mem_debug();
PNEW(mutex_deinit);

static unsigned int  memcnt=0;
static size_t memalloc=0, memfree=0;
/*max depth in backtrace */
#define DT 8 //matches funtrace_len/8

#define METHOD 1
#if METHOD == 1 //Use hash table
typedef struct{
	void *p;
	unsigned int size;
	int nfunc;
	union{
		char funtrace[funtrace_len];//new way of recording initial call location (not yet used)
		void *func[DT];
	};
}memkey_t;
memkey_t *memkey_all=NULL;
unsigned int memkey_len=0; //should be all 1 at lower bits
unsigned int memkey_thres =0;
unsigned int memkey_maxadd=0;
unsigned int memkey_maxdel=0;
unsigned int memkey_notfound=0;
//Enable or disable MEM_DEBUG
static void memkey_init(int enabled){
	if(memkey_len){
		print_mem_debug();
	}
	if(enabled && !memkey_len){
		memcnt=0;
		memalloc=0;
		memfree=0;
		memkey_maxadd=0;
		memkey_maxdel=0;
		memkey_len=0x1FFFFFF;
		dbg("initializing memkey_all for %u of %lu bytes\n", memkey_len, sizeof(memkey_t));
		memkey_all=calloc_default(memkey_len+1, sizeof(memkey_t));
		memkey_thres=memkey_len>>1;
		MEM_DEBUG=enabled;//enable after we setup.
	}else if(!enabled && memkey_len){
		MEM_DEBUG=enabled;//disable before we clean up.
		dbg("free memkey_all and disable MEM_DEBUG\n");
		memkey_len=0;
		memkey_thres=0;
		free(memkey_all);
		memkey_all=0;
	}
	
}
//convert address to index. Make sure last two bits has at least 1 set
static unsigned int addr2ind(void*p){
	unsigned long ind=(unsigned long)p;
	do{
		ind=ind>>1;
	}while(!(ind & 3));
	return (unsigned int)(ind&memkey_len);
}
//add memory pointer to list
static void memkey_add(void *p, size_t size){
	if(atomic_add_fetch(&memcnt, 1)>memkey_thres){
		memkey_thres=memcnt;
		dbg_once("memcnt=%u >= memkey_len=%u. Half slots are filled. please increase memkey_max.\n", memcnt, memkey_len);
	}
	__atomic_add_fetch(&memalloc, size, MEM_ORDER);
	unsigned int counter=0;
	for(unsigned int ind=addr2ind(p); counter<=memkey_len ; ind=(ind+1)&memkey_len, counter++){
		if(!memkey_all[ind].p){//found empty slot
			void *dummy=NULL;
			//try to occupy slot using CAS
			if(atomic_compare_exchange_n(&memkey_all[ind].p, &dummy, p)){//success
				memkey_all[ind].size=(unsigned int)size;
				if(funtrace[0]){
					memcpy(memkey_all[ind].funtrace, funtrace, funtrace_len);
				} else{
					memkey_all[ind].nfunc=backtrace(memkey_all[ind].func, DT);
				}
				if(counter>memkey_maxadd){
					memkey_maxadd=counter;
					if(memkey_maxadd>100){
						dbg("Max try is %u for add\n", memkey_maxadd);
					}
				}
				break;
			}//else: failed, occupied by another thread. Try next slot
		}
	}
}
//Memory allocated from exterior library will not be found. We limit number of tries to avoid searching too much
static size_t memkey_del(void *p){
	size_t size=0;
	unsigned int counter=0;
	for(unsigned int ind=addr2ind(p); counter<=memkey_maxadd ; ind=(ind+1)&memkey_len, counter++){
		if(memkey_all[ind].p==p){//found. no race condition can occure
			__atomic_add_fetch(&memfree, memkey_all[ind].size, MEM_ORDER);
			atomic_sub_fetch(&memcnt, 1);
			memkey_all[ind].p=0;
			size=memkey_all[ind].size;
			break;
		}
	}
	if(counter>memkey_maxdel){
		memkey_maxdel=counter;
		if(memkey_maxdel>100){
			dbg("Max try is %u for del\n", memkey_maxdel);
		}
	}
	if(counter>memkey_maxadd){
		memkey_notfound++;
		//dbg("%16p: unable to find after %u tries\n", p, memkey_maxadd);
		//print_backtrace();
	}
	return size;
}
static void print_mem_debug(){
	if(!memkey_len){
		dbg("MEM_DEBUG was not enabled, canceled.\n");
		return;
	}
	if(!memcnt){
		info3("All allocated memory are freed.\n");
	} else{
		if(memkey_notfound) info3("%u memory is not found to delete.\n", memkey_notfound);
		if(memalloc > memfree){
			info3("%u (%'lu B) allocated memory not freed. \n", memcnt, (memalloc-memfree));
		}
	}
	if(signal_caught){
		info3("Signal is caught, will not print detailed memory usage.\n");
	}else{
		unsigned int counter=0;
		int ans=0;
		for(unsigned int i=0; i<=memkey_len; i++){
			if(memkey_all[i].p){
				if(counter<50){
					counter++;
					info3("%p %9dB", memkey_all[i].p, memkey_all[i].size);
					if(memkey_all[i].nfunc){
						int offset=memkey_all[i].nfunc>3?1:0;
						if(ans||(ans=print_backtrace_symbol(memkey_all[i].func, memkey_all[i].nfunc-offset))){
							info3(" %p\n", memkey_all[i].p);
						}
					}else if(memkey_all[i].funtrace[0]){
						info3(" %s\n", memkey_all[i].funtrace);
					}
				}else{
					info3("Stop after %u prints\n", counter);
					break;
				}
			}
		}
	}
}
#elif METHOD == 2 //linked list. unsolved race condition.
#error "Implementation is not complete. race condition not resolved" 
/*
	Head counter is added by 1 when popped.

	For a node
	* A freshly added node has a counter of 10. 
	* Counter is set to 0 when marked as deleted.
	* Counter is set to 1 when locked by a thread (for deletion).

	For a stack like usage, only the head can be pushed or popped. Push is easy
	as we own the node to be pushed. Pop is more involved as we need to make
	sure both the head and the node to pop is not modified.  Since the first
	node cannot be modified when the head is not modified, we only need to make
	sure the head has not been modified (check next alone is not sufficient) .
	We put a counter on the head and any pop will increment the counter by 1.

	For a list like usage, insert node is again straightforward as we owns the
	node to be inserted. To delete a node, we need to assign the next values of
	the node to the previous node which requires atomic operation on both the
	current node and previous node. This is not possible with single CAS. We
	need to have a way to *lock* or own the node to be deleted so that its value
	is modified, then do a CAS on the previous node. We need to also make sure
	the previous node is not modified (still points to current node). Delete a
	node may fail as the previous node may be changed by another thread. So we
	first mark a node as to be deleted (counter=0), then try to recycle it. To
	actually recycle a node, we first mark it as locked (counter=1), and do CAS
	on the previous node only if the previous node is not locked (counter!=1). 

	Steps to delete a node:
	1. Mark counter to = 0
	2. CAS counter from 0 to 1. Cancel if failed
	3. CAS previous code only if its counter is not 1. Keep its counter.

	When tranversing a list, also need to check for race condition that the node
	is deleted which will have next point the wrong list 

	NOTICE: The implementation does not yet support insert into the list.
	Requires avoid inserting after a node if its counter is 0 or 1.
 */
typedef struct{
	union {
		unsigned long state;
		struct{//in small indian machine, state takes value of next when counter=0
			unsigned int next;
			unsigned int counter;
		};
	};
}memhead_t;
typedef struct{
	memhead_t list;//first element so it can be casted to memkey_head
	void *p;
	size_t size;
	union{
		char funtrace[funtrace_len];//new way of recording initial call location (not yet used)
		struct{
			long nfunc;
			void *func[DT];
		}
	};
}memkey_t;

enum{
	S_DEL =0, //to be deleted
	S_LOCK=1, //locked
	S_IDLE=2, //free node
	S_USED=3  //used ndoe
};
unsigned int memkey_count=0;
unsigned int memkey_i=1;//start with 1.
unsigned int memkey_max=0;
memkey_t *memkey_all=NULL;
memhead_t memhead_idle={0L};
memhead_t memhead_used={0L};
//Place to to the beginning of list.
static void memhead_push(memhead_t *head, unsigned int ind){
	memkey_all[ind].list.next=head->next;
	do{
		if(memkey_all[ind].list.next==ind){//why can this happen?
			warning("ind=%u, head->next=%u, p=%p, counter=%u, %s\n", ind, memkey_all[ind].list.next, memkey_all[ind].p, memkey_all[ind].list.counter, head==&memhead_idle?"idle":"used");
			return;
		}
	}
	while(!atomic_compare_exchange_n(&(head->next), &memkey_all[ind].list.next, ind));
	
}
//Remove memkey from beginning of list (only
static unsigned int memhead_pop(memhead_t *head){
	memhead_t key, key2;
	//compare both counter and pointer to make sure the node is not changed.
	key.counter=atomic_add_fetch(&(head->counter), 1);//increse the counter
	key.next=head->next;
	do{
		key2.counter=key.counter;
		key2.next=memkey_all[key.next].list.next;
	} while(!atomic_compare_exchange(&head->state, &key.state, &key2.state));
	
	return key.next;
}

/**
	Delete node i if the previous node is not locked. Returns True if success.
*/
static int memhead_recycle(memhead_t *head, unsigned int iprev, unsigned int i){
	int ans=0;
	if(i==iprev){
		error("i=%u, iprev=%u; cancel\n", i, iprev);
		return ans;
	}
	//Check and lock current node to be deleted. Important
	unsigned int dummy=S_DEL;
	if(!atomic_compare_exchange_n(&memkey_all[i].list.counter, &dummy, S_LOCK)){
		return ans;//lock node failed. cancel operation. some other thread may have done it.
	}
	//After lock success, we owns node i so can modify its values.
	memhead_t *phead;
	memhead_t key, key2;
	if(iprev==0){//head
		phead=head;
	} else{
		phead=&memkey_all[iprev].list;
	}
	key.counter=phead->counter;//save and check its value
	//Only proceed if the previous node is not locked for delete (head is never marked as such)
	if(key.counter!=S_LOCK){
		key.next=i;//check link is right
		key2.counter=key.counter;//keep its value
		key2.next=memkey_all[i].list.next;
		if(i==memkey_all[i].list.next){
			warning("invalid list: i=%u, next=%u\n", i, memkey_all[i].list.next);
		}
		if(atomic_compare_exchange(&phead->state, &key.state, &key2.state)){//success
			//we already own the node i. Can directly modify its value.
			memkey_all[i].p=NULL;
			memkey_all[i].list.counter=S_IDLE;
			//no longer owns node i after the following step
			memhead_push(&memhead_idle, i);//recycle the node 
			atomic_sub_fetch(&memkey_count, 1);
			//info("%u delete succeed memkey_count=%u\n", i, memkey_count);
			ans=1;
		}//else: phead is changed due to insertion or it self is deleted.
	}//else phead is changed due to insertion or it self is deleted.
		//info("%u is not recycled prev is locked, counter=%u\n", i, key.counter);
	//}
	if(!ans){//failed to delete. unlock the node. set is back to pending recycle
		dummy=S_LOCK;
		if(!atomic_compare_exchange_n(&memkey_all[i].list.counter, &dummy, S_DEL)){
			dbg("%u unlock node failed.\n", i);
		}
	}
	return ans;
}

//delete a node (first mark as delete then recycle nodes)
static void memkey_del(void *p){
	int done=0;
	memhead_t next={0L};
	unsigned int nretry=0;
	if(!p) info("recyling all\n");
retry:
	//mark node as delete if p is valid. Otherwise, just try recyle all deleted nodes
	for(unsigned int i=memhead_used.next, iprev=0; i&&!done; iprev=i, i=next.next){
		next=memkey_all[i].list;//copy both the conter and index
		if(next.counter==S_IDLE){
			nretry++;
			if(nretry<5){
				warning("node is deleted. Retry for %u times\n", nretry);
				goto retry;
			}
			warning("node %u is deleted. cancel for\n", i);
			break;
		}
		if(i&&next.next==i){
			error("i=%u == inext. counter=%u\n", i, next.counter);
			nretry++;
			if(nretry<5){
				goto retry;
			}
			break;
		}

		if(p && memkey_all[i].p==p){//found
			//Mark as current node as delete atomically (do we need atomic?)
			unsigned int counter=memkey_all[i].list.counter;
			while(counter>1&&!atomic_compare_exchange_n(&memkey_all[i].list.counter, &counter, S_DEL)){
				dbg("unexpected: mark delete failed: counter=%u\n", counter);
			}
			if(counter>1){//mark success.
				atomic_sub_fetch(&memcnt,1);
			}else{//Why would the following happen? counter is shown as 0
				//It seems this failure follows failure of memkey_push() which has ind==head->next which is not expected. p is the same as us.
				//It seems the same memory is used by two threads?
				warning("unexpected: mark delete %p at %u canceled: counter=%u\n", p, i, counter);
			}
			done=1;
		}
		if(!memkey_all[i].list.counter){//try recycle
			if(memhead_recycle(&memhead_used, iprev, i)){
				if(!p){
					info("Recycle %u failed\n", i);
				}
			}else{
				if(!p){
					warning("Recycle %u failed\n", i);
				}
			}
		}else if(!p){
			warning("unexpected: %u still has counter %u\n", i, memkey_all[i].list.counter);
		}

	}
	if(p && !done){//should never happen
		dbg("unexpected: %p is not found memcnt=%u, memkey_count=%u, nretry=%u\n", p, memcnt, memkey_count, nretry);
	}
}
//Add a memkey
static void memkey_add(void *p, size_t size){
	unsigned int ind;
	memalloc+=size;
	atomic_add_fetch(&memcnt, 1);
	atomic_add_fetch(&memkey_count, 1);
	int nretry=0;
retry:	
	ind=memhead_pop(&memhead_idle);
	if(!ind){
		ind=atomic_fetch_add(&memkey_i, 1);
		if(ind>=memkey_max){
			warning("All nodes are used. recycle unused\n");
			memkey_del(NULL);//recycle deleted nodes
			nretry++;
			if(nretry<2){
				goto retry;
			}else{
				error("Please increase memkey_max and recompile\n");
				return;
			}
		}
	}else{
		if(memkey_all[ind].list.counter!=S_IDLE){
			dbg("memkey_all[%u].p=%p should have counter=%u but is %u. Retry.\n", ind, memkey_all[ind].p, S_IDLE, memkey_all[ind].list.counter);
			goto retry;
		}
	}
	memkey_all[ind].p=p;
	memkey_all[ind].size=size;
	memkey_all[ind].list.counter=S_USED;
	if(funtrace[0]){
		memcpy(memkey_all[ind].funtrace, funtrace, funtrace_len);
	}else{
		memkey_all[ind].nfunc=backtrace(memkey_all[ind].func, DT);
	}
	memhead_push(&memhead_used, ind);
}

void print_mem_debug(){
	unsigned int i;
	unsigned int print_count=0;
	unsigned int print_max=10;
	memkey_del(NULL);//recycle
	if(!memcnt){
		info3("All allocated memory are freed. %u not recycled\n", memkey_count);
	}else{
		info3("%u (%'lu B) allocated memory not freed. %u not recycled!!!\n", memcnt, (memalloc-memfree), memkey_count);
	}
	for(i=memhead_used.next; i; i=memkey_all[i].list.next){
		print_count++;
		if(print_count<print_max){
			info3("%9zuB (%.40s)", memkey_all[i].size, (char*)memkey_all[i].p);
			if(memkey_all[i].funtrace[0]){
				info3(" %s i=%u p=%p counter=%u\n", memkey_all[i].funtrace, i, memkey_all[i].p, memkey_all[i].list.counter);
			} else if(memkey_all[i].nfunc){
				int offset=memkey_all[i].nfunc>3?1:0;
				if(print_backtrace_symbol(memkey_all[i].func, memkey_all[i].nfunc-offset)){
					info3("print_backtrace_symbol failed, stop.\n");
					print_count=101;
				}
			}
		} else {
			info3("Stop after %d records\n", print_count);
			break;
		}
	}
	
	info3("Memory allocation is %lu kB, still in use %lu KB.\n", memalloc>>10, (memalloc-memfree)>>10);
}
void memkey_init(){
	if(MEM_DEBUG){
		memkey_max=1000000;
		memkey_all=calloc_default(memkey_max, sizeof(memkey_t));
		memhead_used.counter=S_USED;//start with 10
		memhead_idle.counter=S_USED;
	}else{
		free(memkey_all);
		memkey_all=NULL;
		memkey_max=0;
}
#else
typedef struct{
	void *p;
	size_t size;
	union{
		char funtrace[funtrace_len];//new way of recording initial call location (not yet used)
		void *func[DT];
	};
	int nfunc;
}memkey_t;
PNEW(mutex_mem);
static void *MROOT=NULL;
static void *MSTATROOT=NULL;//reuses keys from MROOT.

typedef int(*compar)(const void *, const void *);
static int stat_cmp(const memkey_t *pa, const memkey_t *pb){
	if(!pa->nfunc&&!pb->nfunc){
		return strcmp(pa->funtrace, pb->funtrace);
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
	const memkey_t *key2=*((const memkey_t **)key);
	(void)level;
	//merge calls from the same location
	if(which==leaf||which==postorder){
		void *found;
		if(!(found=tfind(key2, &MSTATROOT, (compar)stat_cmp))){//not found
			if(!tsearch(key2, &MSTATROOT, (compar)stat_cmp)){
				warning("Error inserting to tree\n");
			}
		} else{
			(*(memkey_t **)found)->size+=key2->size;
		}
	}
}

static int print_count=0;
static void print_usage(const void *key, VISIT which, int level){
	const memkey_t *key2=*((const memkey_t **)key);
	(void)level;
	if(which==leaf||which==postorder){
		print_count++;
		if(print_count<100){
			info3("%p(%9zu)", key2->p, key2->size);
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

void print_mem_debug(){
	if(MROOT){
		info3("%u (%lu MB) allocated memory not freed!!!\n", memcnt, (memalloc-memfree)>>20);
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
			info3("Memory counter is still not zero: %u\n", memcnt);
		}
	}
	info3("Memory allocation is %lu kB, still in use %lu KB.\n", memalloc>>10, (memalloc-memfree)>>10);
}
static int key_cmp(const void *a, const void *b){
	void *p1, *p2;
	if(!a||!b) return 1;
	p1=((memkey_t *)a)->p;
	p2=((memkey_t *)b)->p;
	if(p1<p2) return -1;
	else if(p1>p2) return 1;
	else return 0;
}
static void memkey_add(void *p, size_t size){
	if(!p){
		if(size){
			warning("memory allocation for %zu failed\n", size);
		}
		return;
	}

	memkey_t *key=(memkey_t *)calloc_default(1, sizeof(memkey_t));
	key->p=p;
	key->size=size;
	
	if(funtrace[0]){//do not use backtrace
		memcpy(key->funtrace, funtrace, funtrace_len);
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
			memalloc+=size;
		}
	}
	UNLOCK(mutex_mem);
}
//Return 0 if success. 1 if failed.
static int memkey_del(void *p){
	if(!p) return 0;
	int ans=1;
	void *found=0;
	memkey_t key;
	key.p=p;
	LOCK(mutex_mem);
	found=tfind(&key, &MROOT, key_cmp);
	if(found){
		memkey_t *key1=*(memkey_t **)found;/*the address of allocated memkey_t. */
		memfree+=key1->size;
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

static void memkey_init(){

}
#endif

void *calloc_maos(size_t nbyte, size_t nelem){
	void *p=calloc_default(nbyte, nelem);
	if(!p){
		error("calloc for %ld bytes failed (%d): %s\n", nbyte*nelem, errno, strerror(errno));
	}
	if(memkey_len){
		memkey_add(p, nbyte*nelem);
	}
	if(MEM_VERBOSE){
		info("%p %s calloc with %zux%zu bytes\n", p, funtrace, nelem, nbyte);
	}
	funtrace_unset;
	return p;
}
void *malloc_maos(size_t size){
	void *p=malloc_default(size);
	if(!p){
		error("malloc for %ld bytes failed (%d): %s\n", size, errno, strerror(errno));
	}
	if(memkey_len){
		memkey_add(p, size);
	}
	if(MEM_VERBOSE){
		info("%p %s malloc with %zu bytes\n", p, funtrace, size);
	}
	funtrace_unset;
	return p;
}
void *realloc_maos(void *p0, size_t size){
	//to avoid race condition if p0 and p is different and p0 is allocated to another thread
	//We delete p0 first and then readd p.
	if(memkey_len&&p0){
		memkey_del(p0);
	}
	void *p=realloc_default(p0, size);
	if(size){
		if(!p){
			error("realloc for %ld bytes failed (%d): %s\n", size, errno, strerror(errno));
		}
		if(memkey_len){
			memkey_add(p, size);
		}
	}
	if(MEM_VERBOSE){
		info("%p %s realloc with %zu bytes\n", p, funtrace, size);
	}
	funtrace_unset;
	return p;
}
void *realloc2_maos(void *p0, size_t oldsize, size_t newsize){//zero out newly allocated segment
	void* pnew=realloc_maos(p0, newsize);
	if(newsize>oldsize){
		memset((char*)pnew+oldsize, 0, newsize-oldsize);
	}
	return pnew;
}
void free_maos(void *p){
	size_t size=0;
	if(memkey_len && p){
		size=memkey_del(p);//calls free_default when marked success.
	}
	free_default(p);//must be after memkey_del to avoid race condition (add before del)
	if(MEM_VERBOSE){
		info("%p %s freed with %zu bytes.\n", p, funtrace, size);
	}
	funtrace_unset;
}

void read_sys_env(){
	READ_ENV_INT(MEM_DEBUG, 0, 1);//enable memory accounting
	READ_ENV_INT(MEM_VERBOSE, 0, 2);//enable verbose memory operation
	READ_ENV_INT(MEM_FUNTRACE, 0, 1);//enable use funtrace instead of backtrace
	READ_ENV_INT(LOG_LEVEL, -5, 5);
}

static void init_mem(){
	if(!calloc_default){
		calloc_default=(void *(*)(size_t, size_t))dlsym(RTLD_DEFAULT, "calloc");
		malloc_default=(void *(*)(size_t))dlsym(RTLD_DEFAULT, "malloc");
		realloc_default=(void *(*)(void *, size_t))dlsym(RTLD_DEFAULT, "realloc");
		free_default=(void(*)(void *))dlsym(RTLD_DEFAULT, "free");
		read_sys_env();
		void init_process(void);
		init_process();//call before memkey_init to avoid mem counting
		memkey_init(MEM_DEBUG);
		init_hosts();
		//fpconsole=fdopen(dup(fileno(stdout)), "a");
	}
}
static __attribute__((constructor)) void init(){
	init_mem();
}
/**
 * register an external set of malloc functions
 * */
void register_malloc(void* (*ex_malloc)(size_t), void* (*ex_calloc)(size_t, size_t), void* (*ex_realloc)(void *, size_t), void(*ex_free)(void *)){
	if(ex_malloc)  malloc_default= ex_malloc;
	if(ex_calloc)  calloc_default= ex_calloc;
	if(ex_realloc) realloc_default=ex_realloc;
	if(ex_free)    free_default=   ex_free;
}
typedef struct deinit_t{/*contains either fun or data that need to be freed. */
	void (*fun)(void);
	void *data;
	struct deinit_t *next;
}deinit_t;
deinit_t *deinit_head=NULL;

/**
   Register routines to be called with mem.c is unloading (deinit).
 */
static __attribute__((destructor)) void deinit(){
	void freepath();
	void thread_pool_destroy();
	thread_pool_destroy();
	//remove files that are 365 days old.
	remove_file_older(DIRCACHE, 1, 30*24*3600);//1 month
	remove_file_older(DIRLOCK, 1, 30*24*3600);//1 month
	freepath();
	free_hosts();
	for(deinit_t *p1=deinit_head;p1;p1=deinit_head){
		deinit_head=p1->next;
		if(p1->fun) p1->fun();
		if(p1->data) myfree(p1->data);
		free_default(p1);
	}
	if(MEM_DEBUG){
		print_mem_debug();
	}
	free_process();//call after print_mem_debug which needs TEMP
	if(fplog) fclose(fplog);
}

/**
   Register a function or data to call or free upon exit
*/
void register_deinit(void (*fun)(void), void *data){
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

static int (*signal_handler)(int)=0;//If set will be called by default_signal_handler to avoid exit.
static volatile sig_atomic_t fatal_error_in_progress=0;
void default_signal_handler(int sig, siginfo_t *siginfo, void *unused){
	(void)unused;
	if(sig==0){
		dbg_time("Signal 0 caught. do nothing\n");
		return;
	}else if(sig==SIGUSR1){//Use SIGUSR1 to enable MEM_DEBUG and print memory infromation
		memkey_init(1);
		return;
	}else if(sig==SIGTERM && siginfo){
		char sender[PATH_MAX]={0};
		get_job_progname(sender, PATH_MAX, siginfo->si_pid);
		info("Code is %d, send by %d (uid=%d, %s).\n",
			siginfo->si_code, siginfo->si_pid, siginfo->si_uid, sender);
	}//SIGUSR2 is raised by error()	
	signal_caught=sig;
	info("\nSignal caught: %s (%d)\n", strsignal(sig), sig);
	//sync();
	int cancel_action=0;

	if(fatal_error_in_progress){
		info("Signal handler is already in progress. Force quit without cleanup.\n");
		raise(SIGTERM);
	}else{
		fatal_error_in_progress++;
	}
	if(sig!=SIGINT){
		if(iscrash(sig)){
			if(siginfo&&siginfo->si_addr){
				info("Memory location: %p\n", siginfo->si_addr);
			}
			struct itimerval timer={{10,0},{10,0}};
			if(setitimer(ITIMER_REAL, &timer, NULL)){
				warning("Setitimer failed\n");
			}
			//It is not safe to call backtrace in SIGSEGV, to prevent hang, we set an alarm to force quit
		}
		print_backtrace();
	}
	if(signal_handler){
		dbg_time("Call signal_handler %p\n", signal_handler);
		cancel_action=signal_handler(sig);
		dbg_time("Signal handler returns %d\n", cancel_action);
	}
	if(!cancel_action){//Propagate signal to default handler.
		dbg_time("Propagate signal\n");
		struct sigaction act={0};
		act.sa_handler=SIG_DFL;
		sigaction(sig, &act, 0);
		raise(sig);
	}else{//cancel signal, keep going
		dbg_time("Cancel action\n");
		sync();
		if(sig==SIGUSR2){//cancel error();
			fatal_error_in_progress=0;
		}
	}
}
int dummy_signal_handler(int sig){
	info2("Signal %d caught, will not quit.\n", sig);
	return 1;
}

/**
   Register signal handler
*/
void register_signal_handler(int (*func)(int)){
	struct sigaction act={0};
	act.sa_sigaction=default_signal_handler;
	act.sa_flags=SA_SIGINFO|SA_RESETHAND;//SA_RESETHAND resets the handler to default after one shot
	sigaction(SIGBUS, &act, 0);
	sigaction(SIGILL, &act, 0);
	sigaction(SIGSEGV, &act, 0);
	sigaction(SIGINT, &act, 0);
	sigaction(SIGALRM, &act, 0);
	sigaction(SIGTERM, &act, 0);
	sigaction(SIGABRT, &act, 0);
	//sigaction(SIGHUP, &act, 0);
	signal(SIGHUP, SIG_IGN);//ignore the signal.
	sigaction(SIGUSR1, &act, 0);
	sigaction(SIGUSR2, &act, 0);
	sigaction(SIGQUIT, &act, 0);
	signal_handler=func;
}
