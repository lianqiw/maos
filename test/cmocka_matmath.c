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
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdint.h>
#include <cmocka.h>
#include "../lib/aos.h"
void test_embed(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	for(long nx=6; nx<8; nx++){
		for(double ang=0; ang<M_PI/2; ang+=M_PI/3){
			dmat *in=dnew(nx,nx);
			drandn(in, 1, &rstat);
			cmat *cin=cnew(nx,nx);
			cembedd(cin, in, 0);
			//dshow(in, "in");
			for(long i=0; i<10; i++){
				dmat *out=dnew(i, i);
				dembed(out, in, ang);
				//dshow(out, "out");
				if(i>0){
					assert_float_equal(P(in, nx>>1, nx>>1),P(out, i>>1, i>>1),1e-15);
				}
				dfree(out);
				cmat *cout=cnew(i, i);
				cembedd(cout, in, ang);
				if(i>0){
					assert_float_equal(P(in, nx>>1, nx>>1), creal(P(cout, i>>1, i>>1)), 1e-15);
				}
				cmat *cout2=cnew(i, i);
				cembed(cout2, cin, ang);
				assert_float_equal(csumdiffsq(cout, cout2),0, 1e-12);
				cfree(cout);
				cfree(cout2);
			}
			cfree(cin);
			dfree(in);
		}
	}
}
void test_denc(void **state){
	(void)state;
	dmat *screen=(dmat *)genatm_simple(0.186, 0, -11./3., 1./64., 129, 1);
	/*dmat *screen=dnew(129,129);
	rand_t rstat;
	seed_rand(&rstat, 1);
	drandn(screen,1,&rstat);*/
	dmat *xx, *enc;
	real ss1=0, ss2=0, ss3=0;
	xx=dlinspace(0, 1, 60); enc=denc(screen, xx, -1, 10); ss1=P(enc, 10); dfree(enc); dfree(xx);
	xx=dlinspace(0, 1, 70); enc=denc(screen, xx, -1, 10); ss2=P(enc, 10); dfree(enc); dfree(xx);
	dmat *screen2=dnew(128, 128);
	dembed(screen2, screen, 0);
	xx=dlinspace(0, 1, 70); enc=denc(screen2, xx, -1, 10); ss3=P(enc, 10); dfree(enc); dfree(xx);
	//info("ss=%g,%g,%g\n", ss1, ss2, ss3);
	assert_float_equal(ss1, ss2, 1e-9);
	assert_float_equal(ss1, ss3, 1e-9);
	dfree(enc);
	dfree(xx);
	dfree(screen);
}
void test_dcircle(void **state){
	(void)state;
	//test dcircle, dcog, dshift2center
	long nn=10;
	dmat *A=dnew(nn, nn);
	real r=3;
	real val=2;
	real cx=3;
	real cy=5;
	assert_int_equal(dcircle(NULL, cx, cy, 1, 1, r, val), -1);
	assert_int_equal(dcircle(A, cx, cy, 1, 1, -1, val), -1);
	assert_int_equal(dcircle(A, cx, cy, 1, 1, r, val), 0);
	//dshow(A, "A");
	real r2=sqrt(dsum(A)/val/M_PI);
	assert_float_equal(r, r2, 0.05);
	//test dcog 
	real grad[2];
	dcog(grad, A, 0, 0, 0, 0, 0, NULL);
	assert_float_equal(grad[0], cx-(nn-1)*0.5, 0.1);
	assert_float_equal(grad[1], cy-(nn-1)*0.5, 0.1);
	//test shift2center
	dshift2center(A, 0, 0);
	dcog(grad, A, 0, 0, 0, 0, 0, NULL);
	assert_float_equal(grad[0], 0, 0.1);
	assert_float_equal(grad[1], 0, 0.1);
	dshift2center(A, 0.1, -0.7);
	dcog(grad, A, 0, 0, 0, 0, 0, NULL);
	assert_float_equal(grad[0], 0.1, 0.1);
	assert_float_equal(grad[1], -0.7, 0.1);

	r2=sqrt(dsum(A)/M_PI);
	assert_float_equal(r+1, r2, 0.5);
	info("r=%g, r2=%g, grad=%g %g\n", r, r2, grad[0], grad[1]);
}

int main(void){
	//register_signal_handler(dummy_signal_handler);//captures error().
	LOG_LEVEL=0;//set higher level to suppress print out
	const struct CMUnitTest tests[]={
		cmocka_unit_test(test_embed),
		cmocka_unit_test(test_denc),
		cmocka_unit_test(test_dcircle),
	};
	return cmocka_run_group_tests(tests, NULL, NULL);
}
