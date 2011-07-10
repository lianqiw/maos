/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#if defined(__APPLE__)
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/proc.h>
#include <sys/sysctl.h>
#include <mach-o/dyld.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#include <mach/mach_types.h>
#include <mach/task_info.h>
#include <mach/task.h>
#include <mach/vm_statistics.h>
#include <mach/vm_map.h>
#include "common.h"
#include "process.h"
#include "misc.h"
const char *get_job_progname(void){
    static char *progname=NULL;
    if(!progname){
	char path[PATH_MAX],path2[PATH_MAX];
	uint32_t size=sizeof(path);
	if(_NSGetExecutablePath(path,&size)==0){
	    if(realpath(path,path2)){
		progname=strdup0(path2);
	    }
	}
    }
    return progname;
}

int get_job_mem(void){//return in KiB
    int mem;

    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS != task_info(mach_task_self(),
				  TASK_BASIC_INFO, (task_info_t)&t_info, 
				  &t_info_count)){
	mem=0;
	// resident size is in t_info.resident_size;
	// virtual size is in t_info.virtual_size;
    }else{
	mem=t_info.resident_size/1024;
    }
    info("mem=%d KiB\n",mem);
    //struct rusage usage;
    //getrusage(0, &usage);
    //info("mem2=%d KiB\n",(int)usage.ru_maxrss);
    return mem;
}

double get_job_launchtime(int pid){
    double starttime;
    struct timeval tp;
    gettimeofday(&tp, NULL);
    //set to the current time
    starttime=tp.tv_sec+tp.tv_usec;
    //warning("Please finish it\n");
    return starttime;
}

int get_usage_running(void){
    //warning("Please implement it\n");
    return 0;
}
double get_usage_load(void){
    double load=0;
    return load;
    //warning("Please implement");
}
double get_usage_mem(void){
    double mem=0;
    //warning("Please implement");
    return mem;
}
double read_self_cpu(void){
    //warning("Please implement");
    return 0;
}
int read_usage_cpu(long *user, long *tot){
    processor_cpu_load_info_data_t *pinfo;
    mach_msg_type_number_t info_count;
    //long tot1=0;
    //long idle1=0;
    unsigned int ncpu0;
    if (host_processor_info (mach_host_self (),
			     PROCESSOR_CPU_LOAD_INFO,
			     &ncpu0,
			     (void*) &pinfo,
			     &info_count)) {
	return -1;
    }
    for (int i = 0; i < ncpu0; i++) {
	*user+=pinfo[i].cpu_ticks[CPU_STATE_USER]
	    +pinfo[i].cpu_ticks[CPU_STATE_NICE]
	    +pinfo[i].cpu_ticks[CPU_STATE_SYSTEM];
	*tot=*user+pinfo[i].cpu_ticks[CPU_STATE_IDLE];
    }
    return 0;
}

int get_ncpu(void){
    int mib[2]={CTL_HW, HW_NCPU};
    int maxproc=1;
    size_t len=sizeof(maxproc);
    if(sysctl(mib,2,&maxproc,&len,NULL,0)==-1){
	error("sysctl failed\n");
    }
    return maxproc;
}

#endif
