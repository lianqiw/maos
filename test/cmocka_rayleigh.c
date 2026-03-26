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

static int setup_rayleigh(void **state){
	rayleigh_t *cfg=rayleigh_setup(0);
	*state=cfg;
	return 0;
}
static void test_rayleigh(void **state){
	rayleigh_t *cfg=*state;
	dcell *res=rayleigh(cfg);
	writebin(cfg->saloc, "rayleigh_saloc");
	writebin(cfg->tau, "rayleigh_tau");
	writebin(res, "rayleigh_res.bin");
	dcellfree(res);
}
static int free_rayleigh(void **state){
	rayleigh_t *cfg=*state;
	rayleigh_free(cfg);
	return 0;
}
int main(void){
	register_signal_handler(dummy_signal_handler);//captures error().
	LOG_LEVEL=0;//set higher level to suppress print out
	const struct CMUnitTest tests[]={
		cmocka_unit_test(test_rayleigh),
	};
	return cmocka_run_group_tests(tests,setup_rayleigh, free_rayleigh);
}
