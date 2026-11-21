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
/**
   \file unittest.h
   Common header for the cmocka based unit tests.

   cmocka.h must be included after aos.h, because its fail() macro breaks
   std::basic_ios::fail() in the C++ headers when compiled as C++. Wrap
   complex results in creal() or cabs() before passing them to assert_float_equal().

   Tests must not read or write temporary files: generate input data in memory
   and check the results with assertions.
*/
#ifndef AOS_TEST_UNITTEST_H
#define AOS_TEST_UNITTEST_H
#include "../lib/aos.h"
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C"{ //some versions of cmocka.h lack the guard
#endif
#include <cmocka.h>
#ifdef __cplusplus
}
#endif

/**
   Check |a-b|<=eps in double precision. Used by assert_float_equal().
*/
static inline void unittest_near(double a, double b, double eps, const char *sa, const char *sb, const char *file, int line){
	if(!(fabs(a-b)<=eps)){
		print_error("%s=%.17g != %s=%.17g (tolerance %g)\n", sa, a, sb, b, eps);
		_fail(file, line);
	}
}
/*Some versions of cmocka cast the arguments of assert_float_equal() to float,
  which makes small tolerances meaningless. Compare in double instead.*/
#undef assert_float_equal
#define assert_float_equal(a, b, eps) unittest_near((double)(a), (double)(b), (double)(eps), #a, #b, __FILE__, __LINE__)

/**
   Assert that two dmat (or two smat in single precision) agree to within tol in max absolute difference.
*/
#define assert_dmat_near(A, B, tol) do{\
	const dmat *A_=(A); const dmat *B_=(B);\
	assert_non_null(A_); assert_non_null(B_);\
	assert_int_equal(NX(A_), NX(B_));\
	assert_int_equal(NY(A_), NY(B_));\
	real maxd_=0;\
	for(long i_=0; i_<PN(A_); i_++){\
		real d_=fabs(P(A_,i_)-P(B_,i_));\
		if(!(d_<=maxd_)) maxd_=d_;\
	}\
	assert_float_equal(maxd_, 0, (tol));\
}while(0)

/**
   Run a group of tests with the standard settings: error() is captured by
   the signal handler and log output is suppressed.
*/
#define UNITTEST_MAIN(...) int main(void){\
	register_signal_handler(dummy_signal_handler);\
	LOG_LEVEL=-4;\
	const struct CMUnitTest tests[]={__VA_ARGS__};\
	return cmocka_run_group_tests(tests, NULL, NULL);\
}
#endif
