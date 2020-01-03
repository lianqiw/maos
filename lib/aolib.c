/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include <stdio.h>
#include "aos.h"
//jmp_buf exception_return;
static void py_quitfun(const char *msg){
    info("%s", msg);
    //longjmp(exception_return, 1);
}
static void py_signal_handler(int sig){
    info("%d caught\n", sig);
    //longjmp(exception_return, sig);
}
static void(*default_handler)(int)=NULL;
static __attribute__((constructor)) void init(){
    if(!default_handler){
	default_handler=signal(SIGTERM, py_signal_handler);
    }
    quitfun=py_quitfun;
    //function that calls setjmp() cannot return before you call longjmp().
    //if(setjmp(exception_return)){
//	info("Error setting jmp_buf\n");
    // }
}
static __attribute__((destructor)) void deinit(){
    fprintf(stderr, "aolib unloaded\n");
    if(default_handler){
	signal(SIGTERM, default_handler);
    }else{
	signal(SIGTERM, SIG_DFL);
    }
}
