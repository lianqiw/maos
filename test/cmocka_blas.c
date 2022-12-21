/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>

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

static void blas_dpinv(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *A=dnew(10, 4);
	drandn(A, 1, &rstat);
	dmat *w=dnew(10, 1);
	dset(w, 2);
	dmat *Ap=dpinv(A, CELL(w));
	dmat *ApA=NULL;
	dmm(&ApA, 0, Ap, A, "nn", 1);
	dsp *spw=dspnewdiag(w->nx, P(w), 1);
	dmat *Ap2=dpinv(A, CELL(spw));
	assert_float_equal(dsum(Ap), 0.675211, 1e-6);
	assert_float_equal(ddiff(Ap, Ap2), 0, 1e-6);

	dcell *At=dcellnew_same(2,2,5,2);
	dcell *Ac=NULL;
	d2cell(&Ac, A, At);
	dcell *Acp=dcellpinv2(Ac, NULL, 1e-14, 0);
	dmat *Acp2=dcell2m(Acp);
	assert_float_equal(ddiff(Acp2, Ap2), 0, 1e-6);

	dfree(A); dfree(w); dfree(Ap); dfree(ApA); dspfree(spw); dfree(Ap2);
	dcellfree(At);
	dcellfree(Ac);
	dcellfree(Acp);
	dfree(Acp2);
}

int main(void){
	register_signal_handler(dummy_signal_handler);//captures error().
	LOG_LEVEL=0;//set higher level to suppress print out
	const struct CMUnitTest tests[]={
		cmocka_unit_test(blas_dpinv),
	};
	return cmocka_run_group_tests(tests, NULL, NULL);
}
