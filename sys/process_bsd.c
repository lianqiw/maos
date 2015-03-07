/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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




#include <limits.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <dirent.h>
#include <sys/resource.h>
#include "common.h"
#include "process.h"

int get_job_progname(char *res, int nres,int pid){
    int mib[4];
    mib[0] = CTL_KERN;
    mib[1] = KERN_PROC;
    mib[2] = KERN_PROC_PATHNAME;
    mib[3] = pid>0?pid:getpid();
    char buf[PATH_MAX];
    size_t cb = sizeof(buf);
    sysctl(mib, 4, buf, &cb, NULL, 0);
    strncpy(res, nres, path2); res[nres-1]=0;
    return 0;
}
int get_job_mem(void){/*return in KiB */
    int mem;
    struct rusage usage;
    getrusage(0, &usage);
    mem=usage.ru_idrss;
    return mem;
}
double get_job_launchtime(int pid){
#warning "Please implement"
    double starttime;
    return 0;
}
int get_usage_running(void){
#warning "Please implement"
  return 0;
}
double get_usage_mem(void){
#warning "Please implement"
  return 0;
}
int read_cpu_counter(long *user, long *tot){
#warning "Please implement"
  return 0;
}
double get_usage_load(void){
#warning "Please implement"
  return 0;
}
int get_ncpu(void){
  return sysconf(_SC_NPROCESSORS_ONLN);
}
double read_self_cpu(void){
#warning "Please implement"
  return 0;
}
#endif
