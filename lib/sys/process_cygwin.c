/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#if defined (__CYGWIN__)
#include "process.h"
/*Largely not implemented. */
const char *get_job_progname(void){
    return "maos";
}
int get_job_mem(void){
    return 0;
}
double get_job_launchtime(int pid){
    (void) pid;
    return 0;
}
int get_usage_running(void){
    return 0;
}
double get_usage_load(void){
    return 0;
}
double get_usage_mem(void){
    return 0;
}
double read_self_cpu(void){
    return 0;
}
int read_usage_cpu(long *user, long *tot){
    *user=0;
    *tot=0;
    return 0;
}
int get_ncpu(void){
    return 1;
}
#endif
