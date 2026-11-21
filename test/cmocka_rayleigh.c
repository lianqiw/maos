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

static real cellsumabs(const dcell *A){
	real sum=0;
	for(long i=0; i<PN(A); i++) sum+=dsumabs(P(A, i));
	return sum;
}
static real cellmaxabs(const dcell *A){
	real val=0;
	for(long i=0; i<PN(A); i++) val=MAX(val, dmaxabs(P(A, i)));
	return val;
}
static int setup_rayleigh(void **state){
	rayleigh_t *cfg=rayleigh_setup(0);
	cfg->npsf=64;//smaller than default to run fast
	*state=cfg;
	return 0;
}
static int free_rayleigh(void **state){
	rayleigh_t *cfg=(rayleigh_t *)*state;
	rayleigh_free(cfg);
	return 0;
}
/** Output has one npsf*npsf nonnegative pattern per subaperture and WFS, and is not empty. */
static void test_rayleigh(void **state){
	rayleigh_t *cfg=(rayleigh_t *)*state;
	assert_null(rayleigh(NULL));
	dcell *res=rayleigh(cfg);
	assert_non_null(res);
	assert_int_equal(NX(res), cfg->saloc->nloc);
	assert_int_equal(NY(res), PN(cfg->thetax));
	real tot=0;
	for(long i=0; i<PN(res); i++){
		dmat *pi=P(res, i);
		assert_int_equal(NX(pi), cfg->npsf);
		assert_int_equal(NY(pi), cfg->npsf);
		assert_true(dmin(pi)>=0);
		assert_true(isfinite(dsum(pi)));
		tot+=dsum(pi);
	}
	assert_true(tot>0);
	dcellfree(res);
}
/** Backscatter is linear in the scattering fraction and quadratic in the throughput (up and down). */
static void test_rayleigh_scaling(void **state){
	rayleigh_t *cfg=(rayleigh_t *)*state;
	dcell *res0=rayleigh(cfg);
	dmat *tau0=ddup(cfg->tau);
	for(long ih=0; ih<NX(cfg->tau); ih++){
		P(cfg->tau, ih, 2)*=3;
	}
	dcell *res1=rayleigh(cfg);
	dcp(&cfg->tau, tau0);
	for(long ih=0; ih<NX(cfg->tau); ih++){
		P(cfg->tau, ih, 1)*=0.5;
	}
	dcell *res2=rayleigh(cfg);
	dcp(&cfg->tau, tau0);
	real s0=dcellsum(res0);
	assert_float_equal(dcellsum(res1)/s0, 3, 1e-12);
	assert_float_equal(dcellsum(res2)/s0, 0.25, 1e-12);
	dcelladd(&res1, 1, res0, -3);
	assert_float_equal(cellmaxabs(res1), 0, 1e-12*cellmaxabs(res0));
	dcellfree(res0);
	dcellfree(res1);
	dcellfree(res2);
	dfree(tau0);
}
/** No backscatter without scattering, or when the laser beams are far outside the field of view. */
static void test_rayleigh_zero(void **state){
	rayleigh_t *cfg=(rayleigh_t *)*state;
	dmat *tau0=ddup(cfg->tau);
	for(long ih=0; ih<NX(cfg->tau); ih++){
		P(cfg->tau, ih, 2)=0;
	}
	dcell *res=rayleigh(cfg);
	assert_float_equal(cellsumabs(res), 0, 0);
	dcellfree(res);
	dcp(&cfg->tau, tau0);
	dfree(tau0);

	dset(cfg->llt_ox, 1000);//launch the beams far away
	res=rayleigh(cfg);
	assert_float_equal(cellsumabs(res), 0, 0);
	dcellfree(res);
	dzero(cfg->llt_ox);
}
int main(void){
	register_signal_handler(dummy_signal_handler);//captures error().
	LOG_LEVEL=-4;//suppress print out
	const struct CMUnitTest tests[]={
		cmocka_unit_test(test_rayleigh),
		cmocka_unit_test(test_rayleigh_scaling),
		cmocka_unit_test(test_rayleigh_zero),
	};
	return cmocka_run_group_tests(tests,setup_rayleigh, free_rayleigh);
}
