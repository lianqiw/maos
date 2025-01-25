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

void test_psd2d(){
	dmat *screen=(dmat *)genatm_simple(0.186, 0, -11./3., 1./64., 128, 1);
	dmat *psd=psd2d_aniso(screen, 1./64.);
	real inte1=sqrt(psd_inte2(psd))*1e9;
	//writebin(psd, "psd_%ld", nx);
	dfree(psd);
	dmat *extra=NULL;
	psd=psd2d(&extra, screen, 1./64.);
	real inte2=sqrt(psd_inte2(psd))*1e9;
	//writebin(psd, "psd2_%ld", nx);dfree(psd);
	real rms1=sqrt(dsumsq(screen)/PN(screen))*1e9;
	info("size is %ld rms is %g, psdint is %g, psdint 2 is %g, recon r0=%g, slope=%g\n",
		NX(screen), rms1, inte1, inte2, P(extra, 0), P(extra, 1));
	assert_float_equal(rms1, 169.648, 1e-3);
	assert_float_equal(inte1, 56.4647, 1e-3);
	assert_float_equal(inte2, 216.578, 1e-3);
	assert_float_equal(P(extra,0), 0.128266, 1e-3);
	assert_float_equal(P(extra, 1), -3.68282, 1e-3);
	dfree(extra);
	dfree(screen);
}

void test_stfun(){
	dmat *screen=(dmat *)genatm_simple(0.186, 0, -11./3., 1./64., 128, 1);
	dmat *st=stfun_batch(screen, NULL);
	info("sum=%g, max=%g\n", dsum(st), dmax(st));
	assert_float_equal(dsum(st)*1e9, 3.01473, 1e-3);
	assert_float_equal(dmax(st)*1e13, 1.15708, 1e-3);
	//dshow(st, "st");
	dfree(st);
	dfree(screen);
}
int main(void){
	register_signal_handler(dummy_signal_handler);//captures error(). remove for debugging.
	LOG_LEVEL=0;//set higher level to suppress print out
	const struct CMUnitTest tests[]={
		cmocka_unit_test(test_psd2d),
		cmocka_unit_test(test_stfun),
	};
	return cmocka_run_group_tests(tests, NULL, NULL);
}
