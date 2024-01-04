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

static void mat_fresnel_prop(void **state){
	(void)state;
	if(zfexist("wvf0")){
		cmat *wvf0=cread("wvf0");
		cmat *wvf1=0;
		real dxout=0;
		fresnel_prop(&wvf1, &dxout, wvf0, 1./64, 0.589e-6, 100, 1, 1);
	}
}

int main(void){
	register_signal_handler(dummy_signal_handler);//captures error().
	LOG_LEVEL=0;//set higher level to suppress print out
	const struct CMUnitTest tests[]={
		cmocka_unit_test(mat_fresnel_prop),
	};
	return cmocka_run_group_tests(tests, NULL, NULL);
}
