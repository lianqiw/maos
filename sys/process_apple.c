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

#if defined(__APPLE__)

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
#include <libproc.h>
#include "common.h"
#include "process.h"
#include "misc.h"
/**
   Get the executable name of current process.
*/
int get_job_progname(char* res, int nres, int pid){
	int ans=1;
	if(pid>0){
		char buf[PATH_MAX];
		if(proc_pidpath(pid, buf, sizeof(buf))>0){
			strncpy(res, buf, nres); res[nres-1]=0;
			ans=0;
		}
	} else{
		char path[PATH_MAX], path2[PATH_MAX];
		uint32_t size=sizeof(path);
		if(_NSGetExecutablePath(path, &size)==0){
			if(realpath(path, path2)){
				strncpy(res, path2, nres); res[nres-1]=0;
				ans=0;
			}
		}
	}
	return ans;
}
/**
   Get the memory usage of current process.
 */
int get_job_mem(void){/*return in KiB */
	int mem;

	struct task_basic_info t_info;
	mach_msg_type_number_t t_info_count=TASK_BASIC_INFO_COUNT;

	if(KERN_SUCCESS!=task_info(mach_task_self(),
		TASK_BASIC_INFO, (task_info_t)&t_info,
		&t_info_count)){
		mem=0;
		/* resident size is in t_info.resident_size; */
		/* virtual size is in t_info.virtual_size; */
	} else{
		mem=t_info.resident_size/1024;
	}
	/*struct rusage usage; */
	/*getrusage(0, &usage); */
	/*dbg("mem2=%d KiB\n",(int)usage.ru_maxrss); */
	return mem;
}
/**
   Get the launch time of current process.
*/
double get_job_launchtime(int pid){
	(void)pid;
	double starttime;
	struct timeval tp;
	gettimeofday(&tp, NULL);
	/*set to the current time */
	starttime=tp.tv_sec+tp.tv_usec;
	/*warning("Please finish it\n"); */
	return starttime;
}
/**
   Get the number of running jobs.
*/
int get_usage_running(void){
	/*warning("Please implement it\n"); */
	return 0;
}
/**
   Get the system memory usage.
*/
double get_usage_load(void){
	double load=0;
	return load;
	/*warning("Please implement"); */
}
double get_usage_mem(void){
	double mem=0;
	/*warning("Please implement"); */
	return mem;
}
double read_self_cpu(void){
	/*warning("Please implement"); */
	return 0;
}
int read_cpu_counter(long* user, long* tot){
	processor_cpu_load_info_data_t* pinfo;
	mach_msg_type_number_t info_count;
	/*long tot1=0; */
	/*long idle1=0; */
	unsigned int ncpu0;
	if(host_processor_info(mach_host_self(),
		PROCESSOR_CPU_LOAD_INFO,
		&ncpu0,
		(int**)&pinfo,
		&info_count)){
		return -1;
	}
	*user=0; *tot=0;
	for(unsigned int i=0; i<ncpu0; i++){
		long tmp=pinfo[i].cpu_ticks[CPU_STATE_USER]
			+pinfo[i].cpu_ticks[CPU_STATE_NICE]
			+pinfo[i].cpu_ticks[CPU_STATE_SYSTEM];
		*user+=tmp;
		*tot+=tmp+pinfo[i].cpu_ticks[CPU_STATE_IDLE];
	}
	vm_deallocate(mach_task_self(), (vm_address_t)pinfo, info_count);
	return 0;
}

int get_ncpu(void){
	int mib[2]={CTL_HW, HW_NCPU};
	int maxproc=1;
	size_t len=sizeof(maxproc);
	if(sysctl(mib, 2, &maxproc, &len, NULL, 0)==-1){
		error("sysctl failed\n");
	}
	return maxproc;
}

#endif
