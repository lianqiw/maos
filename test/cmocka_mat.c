/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
static void mat_ref(void** state){
	(void)state;
	dmat *a=dnew(3,3);
	assert_int_equal(mem_isref(a->mem), 0);
  dmat *b=dref(a);
	P(b,2,2)=1;
	assert_float_equal(P(a,2,2), 1,0);
	assert_int_equal(mem_isref(b->mem), 1);
	dfree(a);
	assert_null(a);
	assert_int_equal(mem_isref(b->mem), 0);
	dfree(b);
	assert_null(b);
}

int main(void)
{
    const struct CMUnitTest tests[] = {
        cmocka_unit_test(mat_ref),
    };

    return cmocka_run_group_tests(tests, NULL, NULL);
}
