/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "unittest.h"

static void dummy_loader(unsigned int depth, unsigned int n, int urgent);
typedef struct data_t{
	unsigned int i;
	unsigned int n;
	unsigned int tot;
	unsigned int depth;
	int urgent;
}data_t;

static void* dummy_runner(void *pdata){
	data_t *data=(data_t *)pdata;
	int do_loader=0;
	unsigned int i;
	unsigned int tot=0;
	while((i=atomic_fetch_add(&data->i, 1))<data->n){
		tot+=i;
		if(i+1==data->n&&data->depth){
			do_loader=1;
		}
	}
	atomic_fetch_add(&data->tot, tot);
	if(do_loader){
		//dbg("depth=%u\n", data->depth);
		dummy_loader(data->depth, data->n, data->urgent);
	}
	return NULL;
}
static void dummy_loader(unsigned int depth, unsigned int n, int urgent){
	//dbg("depth=%u\n", depth);
	data_t data={0,n,0,depth-1,urgent};
	tp_counter_t counter={0};
	thread_pool_queue(&counter, dummy_runner, &data, n, urgent);
	thread_pool_wait(&counter, urgent);
	assert_int_equal(data.tot, n*(n-1)/2);
	if(data.tot!=n*(n-1)/2){
		dbg("data.tot=%u, depth=%u\n", data.tot, depth);
	}
}
static void test_thread_pool(void** state){
    (void)state;
	thread_pool_init(NTHREAD);
	dummy_loader(100, 1000, 1);
	dummy_loader(20, 1000, 0);
	dummy_loader(1, 1, 0);
	thread_pool_destroy();
}
/** Memory accounting works with concurrent allocation and free from OpenMP threads. */
static void test_mem(void** state){
	(void)state;
	long count=0;
#pragma omp parallel for reduction(+:count)
	for(int i=0; i<10000; i++){
		int *a=mycalloc(100, int);
		a[0]=1;
		a[99]=1;
		count+=a[0]+a[99]+a[50];
		myfree(a);
	}
	assert_int_equal(count, 20000);
}
UNITTEST_MAIN(
	cmocka_unit_test(test_thread_pool),
	cmocka_unit_test(test_mem),
)
