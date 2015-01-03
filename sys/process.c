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
/**
   Obtain information about system load
   
*/
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <sys/types.h>
#include <dirent.h>
#if defined __CYGWIN__
#include <sys/cygwin.h>
#endif
#include "common.h"
#include "misc.h"
#include "process.h"
#include "daemonize.h"
#ifdef HAVE_NUMA_H
#include <numa.h>
#endif
/**
   A few routines handles process resource.
*/
int NNUMA=1; /*Number of numa nodes*/
int NCPU=0;
int NTHREAD=0;/*NTHREAD=2*NCPU when hyperthreading is enabled. */
int TCK=0;
long NMEM=0;/*Total memory in byte. */
const char *HOME=NULL;
const char *USER=NULL;
char TEMP[PATH_MAX];
char EXEP[PATH_MAX];/*absolute path of the exe.*/
/**
   Set the HOME, TEMP, USER names.
*/
void init_process(void){
#if defined(__CYGWIN__)
    cygwin_internal(CW_SYNC_WINENV);
#endif
    if(!HOME){
	USER=getenv("USER");
	if(!USER){
	    USER=getenv("USERNAME");
	}
#if defined(__CYGWIN__)
	HOME=getenv("USERPROFILE");
	const char *temp=getenv("TMP");
	if(!temp){
	    temp=getenv("TEMP");
	}
	if(!temp || !exist(temp)){
	    temp="C:/Windows/Temp";
	}
	if(!exist(temp)){
	    temp=HOME;/*set to home */
	}
	if(!exist(temp)){
	    error("Unable to determine the path to temporary files");
	}
	HOME=cygwin_create_path(CCP_WIN_A_TO_POSIX,HOME);
	temp=cygwin_create_path(CCP_WIN_A_TO_POSIX,temp);
#else
	HOME=getenv("HOME");
	const char *temp="/tmp";
#endif
	strcpy(TEMP, temp);
	strcat(TEMP, "/maos-");
	strcat(TEMP, USER);
	mymkdir("%s",TEMP);
	mymkdir("%s/.aos/",HOME);

	{/*PATH to executable*/
	    char exepath[PATH_MAX];
	    if(!get_job_progname(exepath, PATH_MAX, 0)){
		char *tmp=strrchr(exepath,'/');
		if(exepath[0]=='/' && tmp){
		    *tmp=0;
		}else{
		    strncpy(exepath, mygetcwd(), PATH_MAX); 
		}
		strncpy(EXEP, exepath, PATH_MAX-1); EXEP[PATH_MAX-1]=0;
	    }else{
		EXEP[0]=0;
	    }
	}
    }

    NCPU= get_ncpu();
    NTHREAD=sysconf( _SC_NPROCESSORS_ONLN );
#ifdef HAVE_NUMA_H
    NNUMA=numa_num_configured_nodes();
#endif
    TCK = sysconf(_SC_CLK_TCK);
#if defined(__linux__)
    FILE *fp=fopen("/proc/meminfo","r");
    if(fp && fscanf(fp, "%*s %ld %*s", &NMEM)==1){
	NMEM=NMEM*1024; 
    }else{
	NMEM=0;
    }
    fclose(fp);
#else
    NMEM=0;/*do not know. */
#endif
}

/**
   Obtain the current usage level of CPU, between 0 and 1.
*/
double get_usage_cpu(void){
    static double lasttime=0;
    double thistime=myclockd();
    static long user1, tot1;
    static double cent=1;
    long user2=0, tot2=0;
    if(thistime >=lasttime+2){/*information was too old. */
	read_cpu_counter(&user1, &tot1);
	usleep(50000);
    }
    if(thistime <=lasttime+0.1){
	return cent;
    }
    read_cpu_counter(&user2, &tot2);
    long user=user2-user1;
    long tot=tot2-tot1;
    if(tot==0) 
	cent=0;
    else
	cent=(double)user/(double)tot;
    lasttime=thistime;
    user1=user2;
    tot1=tot2;
    //cent=cent*NTHREAD/NCPU;/*discount hyperthreading. */
    return cent;
}
/**
   Return number of idle CPUs that are available to run jobs. Do not count hyperthread cores.
*/
int get_cpu_avail(void){
    int avail=0;
    double load=get_usage_load();
    double cent=get_usage_cpu();
    int nrunning=get_usage_running();
    info2("load=%g, %d%%, nrun=%d, ncpu=%d\n", load, (int)(cent*100), nrunning, NCPU);
    if(load>NCPU+1){/*don't want to put too much load on the machine. */
	return 0;
    }
    avail=(int)round((1.-cent)*NCPU);
    if(avail>NCPU-nrunning){
	avail=NCPU-nrunning;
    }
    if(avail<0) avail=0;
    /*info("CPU is %.1f%% Busy. %d running jobs. %d available.\n",cent*100, nrunning, avail); */
    return avail;
}
/**
   Wait for available CPUs in case scheduler is not available.
*/
void wait_cpu(int nthread){
    char fnlock[64];
    snprintf(fnlock,64,"%s/aos.lock", getenv("HOME"));
    int fd;
    fd=lock_file(fnlock,1,-1);
    while(get_cpu_avail()<nthread-1){
	sleep(5);
    }
    close(fd);/*remove lock */
}
#if defined(__APPLE__)
#include "process_apple.c"
#elif defined(__linux__)
#include "process_linux.c"
#elif defined(__CYGWIN__)
#include "process_cygwin.c"
#elif defined(__FreeBSD__)||defined(__NetBSD__)
#include "process_bsd.c"
#else
#error("Unknown plateform")
#endif

