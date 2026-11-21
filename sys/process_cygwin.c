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
#if defined (__CYGWIN__)
#include "process.h"

/**
   Get the executable name of current process. Not implemented; always returns "maos".
   The process functions in this file are largely not implemented for cygwin.

   @param[out] res Buffer for the executable name.
   @param[in] nres Size of res.
   @param[in] pid Unused.
   @return 0.
*/
int get_job_progname(char* res, int nres, int pid){
	(void)pid;
	strscpy(res, "maos", nres); 
	return 0;
}
/**
   Get the memory usage of current process. Not implemented; returns 0.

   @return 0.
*/
size_t get_job_mem(void){
	return 0;
}
/**
   Get the launch time of a process. Not implemented; returns 0.

   @param[in] pid Unused.
   @return 0.
*/
double get_job_launchtime(int pid){
	(void)pid;
	return 0;
}
/**
   Return the number of running jobs. Not implemented; returns 0.

   @return 0.
*/
int get_usage_running(void){
	return 0;
}
/**
   Get the system load. Not implemented; returns 0.

   @return 0.
*/
double get_usage_load(void){
	return 0;
}
/**
   Get the system memory usage. Not implemented; returns 0.

   @return 0.
*/
double get_usage_mem(void){
	return 0;
}
/**
   Return CPU usage of current process. Not implemented; returns 0.

   @return 0.
*/
double read_self_cpu(void){
	return 0;
}
/**
   Read the cpu counter. Not implemented; sets user and tot to 0.

   @param[out] user Set to 0.
   @param[out] tot Set to 0.
   @return 0.
*/
int read_cpu_counter(long* user, long* tot){
	*user=0;
	*tot=0;
	return 0;
}
/**
   Return number of processors. Not implemented; returns 1.

   @return 1.
*/

int get_ncpu(void){
	return 1;
}
#endif
