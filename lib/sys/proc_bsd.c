/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#if defined(__FreeBSD__)||defined(__NetBSD__)
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/resource.h>

#include "common.h"

#include "proc.h"

const char *get_job_progname(void){
    static char *progname=NULL;
    if(!progname){ //Not tested yet.
	int mib[4];
	mib[0] = CTL_KERN;
	mib[1] = KERN_PROC;
	mib[2] = KERN_PROC_PATHNAME;
	mib[3] = -1;
	char buf[PATH_MAX];
	size_t cb = sizeof(buf);
	sysctl(mib, 4, buf, &cb, NULL, 0);
	progname=strdup(buf);
    }
    return progname;
}
int get_job_mem(void){//return in KiB
    int mem;
    struct rusage usage;
    getrusage(0, &usage);
    mem=usage.ru_idrss;
    return mem;
}
double get_job_launchtime(int pid){
    double starttime;
#error "Please implement"
}
int get_usage_running(void){
#error "Please implement"
}
double get_usage_mem(void){
#error "Please implement"
}
int read_usage_cpu(long *user, long *tot){
#error "Please implement"
}
#endif
