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


//#include <setjmp.h>
#include <signal.h>
#include <stdio.h>
#include "aos.h"
//jmp_buf exception_return;

static __attribute__((constructor)) void init(){
	fprintf(stderr, "aolib loaded\n");
	register_signal_handler(dummy_signal_handler);
	//function that calls setjmp() cannot return before you call longjmp().
	//if(setjmp(exception_return)){
//	info2("Error setting jmp_buf\n");
	// }
	//default_signal_handler(SIGUSR2, 0, 0);
}
static __attribute__((destructor)) void deinit(){
	fprintf(stderr, "aolib unloaded\n");
}
