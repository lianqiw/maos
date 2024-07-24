/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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



#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdint.h>
#include <cmocka.h>
#include "../lib/aos.h"

static void readstr_basic(void** state){
	(void)state;
	assert_float_equal(readstr_num("test","1+2*3", NULL), 7, 0);
	assert_float_equal(readstr_num("test","1*2*3", NULL), 6, 0);
	assert_float_equal(readstr_num("test","1*2+3", NULL), 5, 0);
	assert_float_equal(readstr_num("test","1+2+3", NULL), 6, 0);
	assert_float_equal(readstr_num("test","1+2-3", NULL), 0, 0);
	assert_float_equal(readstr_num("test","1+2/3", NULL), 1+2./3, 0.01);
	assert_float_equal(readstr_num("test","1+2+3+4*5", NULL), 1+2+3+4*5, 0.01);
	assert_float_equal(readstr_num("test", "1,", NULL), 1, 0.01);
	assert_float_equal(readstr_num("test", "1;", NULL), 1, 0.01);
	assert_float_equal(readstr_num("test", "1 1", NULL), 1, 0.01);
	assert_float_equal(readstr_num("test", "1;1", NULL), 1, 0.01);
	assert_true(isnan(readstr_num("test","/0.5+2+3+4*5", NULL)));
	assert_true(isnan(readstr_num("test","*0.5+2+3+4*5", NULL)));
	assert_true(!isnan(readstr_num("test","+0.5+2+3+4*5", NULL)));
}
static void readstr_array(void **state){
	(void)state;
	double *ret=0;
	int nrow;
	int ncol;

	//Simple read
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "[]"), 0);
	assert_int_equal(nrow, 0);
	assert_int_equal(ncol, 0);
	free(ret);ret=NULL;

	//Simple read
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "2/[1 2]"), 2);
	assert_float_equal(ret[0], 2, 0.0001);
	assert_float_equal(ret[1],1,0.0001);
	assert_int_equal(nrow, 2);
	assert_int_equal(ncol, 1);
	free(ret);ret=NULL;

	//Simple read of a matrix
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "2/[1 2; 3 4; 5 6]"), 6);
	assert_float_equal(ret[0], 2, 0.0001);
	assert_float_equal(ret[1], 1, 0.0001);
	assert_int_equal(nrow, 2);
	assert_int_equal(ncol, 3);
	free(ret);ret=NULL;

	//Simple read of a matrix; invalid entry
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "2/[1 2; 3 4; 5 ]"), -1);
	free(ret);ret=NULL;

	//Simple read without bracket. Ignore after ;
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "2 1;3"), 2);
	assert_float_equal(ret[0], 2, 0.0001);
	assert_float_equal(ret[1], 1, 0.0001);
	assert_int_equal(nrow, 2);
	assert_int_equal(ncol, 1);
	free(ret);ret=NULL;

	//Relaxed read type 1
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 10, 1, M_DBL, "test", "2/[2]"), 10);
	assert_float_equal(ret[1],1,0.0001);
	assert_float_equal(ret[9],1,0.0001);
	assert_int_equal(nrow,10);
	assert_int_equal(ncol,1);
	free(ret);ret=NULL;
	//Relaxed read type 2
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 10, 2, M_DBL, "test", "2/[1 2]"), 10);
	assert_float_equal(ret[2],1,0.0001);
	assert_float_equal(ret[9],1,0.0001);
	assert_int_equal(nrow,10);
	assert_int_equal(ncol,1);
	free(ret);ret=NULL;

	//Marix
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 0, 0, M_DBL, "test", "3.2/[1/2*3+2*4/2 1*2 1+1*2 1-2; 1 2 3 -4]*2+1*2"), 8);
	assert_int_equal(nrow, 4);
	assert_int_equal(ncol, 2);
	assert_float_equal(ret[0], 3.2/(1./2*3+2.*4/2)*2+1.*2, 0.0001);
	assert_float_equal(ret[7], 3.2/(-4)*2+1.*2, 0.0001);
	free(ret);ret=NULL;

	//If not enough data is available, a relaxed vector is read if relax is set but not 1. matrix dimension is not honored.
	assert_int_equal(readstr_numarr((void **)&ret, &nrow, &ncol, 9, 2, M_DBL, "test", "3.2/[1/2*3+2*4/2 1*2 1+1*2 1-2; 1 2 3 -4]*2+1*2"), 9);
	assert_int_equal(nrow,9);
	assert_int_equal(ncol,1);
	assert_float_equal(ret[8],ret[7],EPS);
}

int main(void){
	register_signal_handler(dummy_signal_handler);
	LOG_LEVEL=0;
	const struct CMUnitTest tests[]={
		cmocka_unit_test(readstr_basic),
		cmocka_unit_test(readstr_array),
	};

	return cmocka_run_group_tests(tests, NULL, NULL);
}
