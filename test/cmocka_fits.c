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
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdint.h>
#include <cmocka.h>
#include "../lib/aos.h"
void test_matbin(void **state){
	(void)state;
    int N=4;
    dmat *A=dnew(N,N);
	dmat *A2=NULL;
    rand_t rstat;
    seed_rand(&rstat, 1);
    drandn(A, 1, &rstat);
    A->keywords=strdup("dx=1/64\ndy=1/64\nLong=kaflasdjfl ksflksjlfjasdflkj asldfj als fals fl asfjalsdf jalsf djasldf jal fsalfkjasdflkj asldfkj asldf \nadf\nwvf=1.e-4");
    writebin(A, "A.bin");
    dfree(A2); A2=dread("A.bin");
	assert_float_equal(ddiff(A, A2), 0, 0);

	writebin(A, "A.bin.gz");
	dfree(A2); A2=dread("A.bin");
	assert_float_equal(ddiff(A, A2), 0, 0);
#if HAVE_LIBZSTD
	writebin(A, "A.bin.zst");
	dfree(A2); A2=dread("A.bin");
	assert_float_equal(ddiff(A, A2), 0, 0);
#endif
	writebin(A, "A.fits");
	dfree(A2); A2=dread("A.fits");
	assert_float_equal(ddiff(A, A2), 0, 0);

	writebin(A, "A.fits.gz");
	dfree(A2); A2=dread("A.fits");
	assert_float_equal(ddiff(A, A2), 0, 0);
#if HAVE_LIBZSTD
	writebin(A, "A.fits.zst");
	dfree(A2); A2=dread("A.fits");
	assert_float_equal(ddiff(A, A2), 0, 0);
#endif
    dfree(A);
	dfree(A2);
}
void test_cellbin(void**state){
	(void)state;
    int N=4;
    dcell *A=dcellnew(2,2);
	dcell *A2=NULL;
    P(A,0)=dnew(N,N);
    P(A,2)=dnew(N,N);
    rand_t rstat;
    seed_rand(&rstat, 1);
    drandn(P(A,0), 1, &rstat);
    drandn(P(A,2), 1, &rstat);
    A->keywords=strdup("dx=1/64\ndy=1/64\nLong=kaflasdjfl ksflksjlfjasdflkj asldfj als fals fl asfjalsdf jalsf djasldf jal fsalfkjasdflkj asldfkj asldf \nadf\nwvf=1.e-4");
    writebin(A, "A.bin");
    cellfree(A2); A2=dcellread("A.bin");
	assert_float_equal(dcelldiff(A, A2), 0, 0);
	
	writebin(A, "A.fits");
	cellfree(A2); A2=dcellread("A.fits");
	assert_float_equal(dcelldiff(A, A2), 0, 0);
return;
	writebin(A, "A.bin.gz");
	cellfree(A2); A2=dcellread("A.bin");
	assert_float_equal(dcelldiff(A, A2), 0, 0);
#if HAVE_LIBZSTD
	writebin(A, "A.bin.zst");
	cellfree(A2); A2=dcellread("A.bin");
	assert_float_equal(dcelldiff(A, A2), 0, 0);
#endif


	writebin(A, "A.fits.gz");
	cellfree(A2); A2=dcellread("A.fits");
	assert_float_equal(dcelldiff(A, A2), 0, 0);
#if HAVE_LIBZSTD
	writebin(A, "A.fits.zst");
	cellfree(A2); A2=dcellread("A.fits");
	assert_float_equal(dcelldiff(A, A2), 0, 0);
#endif
    cellfree(A);
	cellfree(A2);
}

int main(void){
    register_signal_handler(dummy_signal_handler);//captures error(). remove for debugging.
    LOG_LEVEL=-4;//set higher level to suppress print out
    const struct CMUnitTest tests[]={
        cmocka_unit_test(test_matbin),
        cmocka_unit_test(test_cellbin)
    };
    return cmocka_run_group_tests(tests, NULL, NULL);
}
