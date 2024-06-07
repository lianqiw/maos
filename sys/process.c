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


/**
   Obtain information about system
*/

#include <limits.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#if defined __CYGWIN__
#include <sys/cygwin.h>
#endif
#define IN_MEM_C 1 //disable memory management in this file.
#include "common.h"
#include "misc.h"
#include "process.h"
#include "daemonize.h"
#if _OPENMP
#include <omp.h>
#endif
/**
   A few routines handles process resource.
*/
int NCPU=0;
int MAXTHREAD=0;/*MAXTHREAD=2*NCPU when hyperthreading is enabled. */
int NTHREAD=0;//Default to MAXTHREAD
int TCK=0;
int PID=0;
long NMEM=0;/*Total memory in byte. */
const char* HOME=NULL;
const char* USER=NULL;
char* HOST=NULL;
char* TEMP=NULL;//Do not put temp in /tmp as it is automatically cleaned by system.
char* CACHE=NULL;//Directory for caching files that are expensive to compute.
char* EXEP=NULL;/*absolute path of the exe.*/
char* DIRSTART=NULL;//Start up directory. HOME is replaced by ~
/**
   Set the HOME, TEMP, USER names.
*/
void init_process(void){
#if defined(__CYGWIN__)
	cygwin_internal(CW_SYNC_WINENV);
#endif
	//Local host name
	char htemp[255];
	if(gethostname(htemp, sizeof htemp)){
		warning("Unable to get hostname, set to localhost\n");
		HOST=mystrdup("localhost");
	}else{
		char *tmp;
		if((tmp=strchr(htemp, '.'))){
			tmp[0]='\0';
		}
		HOST=mystrdup(htemp);
	}

	//Get User
	USER=getenv("USER");
	if(!USER){
		USER=getenv("USERNAME");
	}
	if(!USER){
		USER="nobody";
	}
	//Get Home
#if defined(__CYGWIN__)
	HOME=getenv("USERPROFILE");
	HOME=cygwin_create_path(CCP_WIN_A_TO_POSIX, HOME);
#else
	HOME=getenv("HOME");
#endif
	if(!HOME || !exist("/tmp")){
		HOME=".";
	}
	//Get Temp directory. It must not be shared by different servers that runs MAOS.
	if(HOME){//new method
		TEMP=stradd(HOME, "/.aos/tmp/", HOST, NULL);
	}else{//deprecated
		const char *TMP=getenv("XDG_RUNTIME_DIR");
		if(TMP && exist(TMP)){
			TEMP=stradd(TMP, "/maos", NULL);
		}else if (exist("/var/tmp")){// /var/tmp is more persistent than /tmp
			TEMP=stradd("/var/tmp/maos-",USER,NULL);
		}else if (exist("/tmp")){
			TEMP=stradd("/tmp/maos-", USER, NULL);
		}else{
			TEMP=stradd(HOME, "/.aos/tmp/", HOST, NULL);
		}
	}
	//Create temporary folders
	mymkdir("%s", TEMP);
	{//preserve compatibility with old maos versions. Notice that it may be removed by automatic tmp cleanup.
		char TEMP2[100];
		snprintf(TEMP2, sizeof(TEMP2), "/tmp/maos-%s", USER);
		if(exist("/tmp") && !exist(TEMP2)){
			mysymlink(TEMP, TEMP2);
		}
	}
	if(!HOME) HOME=TEMP;
	mymkdir("%s/.aos/", HOME);
	CACHE=stradd(HOME, "/.aos/cache",NULL);
	mymkdir("%s", CACHE);
	
	DIRSTART=mygetcwd();
	if(!mystrcmp(DIRSTART, HOME)){
		DIRSTART[0]='~';
		memmove(DIRSTART+1, DIRSTART+strlen(HOME), strlen(DIRSTART)-strlen(HOME)+1);
	}
	{/*PATH to executable*/
		char exepath[PATH_MAX];
		if(!get_job_progname(exepath, PATH_MAX, 0)){
			char* tmp=strrchr(exepath, '/');
			if(exepath[0]=='/'&&tmp){
				*tmp=0;
				EXEP=mystrdup(exepath);
			} else{
				EXEP=mystrdup(DIRSTART);
			}
		}
	}

	NCPU=get_ncpu();
	MAXTHREAD=(int)sysconf(_SC_NPROCESSORS_ONLN);
#if _OPENMP
	//The openmp library may have not yet initialized, so we parse OMP_NUM_THREADS instead.
	if(getenv("OMP_NUM_THREADS")){
		int nthread=strtol(getenv("OMP_NUM_THREADS"), 0, 10);
		if(nthread<MAXTHREAD&&nthread>0){
			MAXTHREAD=nthread;
		}
	}
#endif
	if(!NTHREAD) NTHREAD=MAXTHREAD;
	TCK=(int)sysconf(_SC_CLK_TCK);
#if defined(__linux__)
	FILE* fp=fopen("/proc/meminfo", "r");
	if(fp&&fscanf(fp, "%*s %ld %*s", &NMEM)==1){
		NMEM=NMEM*1024;
	} else{
		NMEM=0;
	}
	fclose(fp);
#else
	NMEM=0;/*do not know. */
#endif
	PID=(int)getpid();
}
/**
 * free memory
 * */
void free_process(){
	/*make sure free_process is called in the end. scheduler_connect_self
	requires TEMP for backtrace printing*/
	free(TEMP); TEMP=NULL;
	free(CACHE); CACHE=NULL;
	free(EXEP); EXEP=NULL;
	free(DIRSTART); DIRSTART=NULL;
	free(HOST); HOST=NULL;
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
	if(thistime>=lasttime+2){/*information was too old. */
		read_cpu_counter(&user1, &tot1);
		mysleep(0.05);
	}
	if(thistime<=lasttime+0.1){
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
	//cent=cent*MAXTHREAD/NCPU;/*discount hyperthreading. */
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
	//info("load=%g, %d%%, nrun=%d, ncpu=%d\n", load, (int)(cent*100), nrunning, NCPU);
	if(load>NCPU+1){/*don't want to put too much load on the machine. */
		return 0;
	}
	avail=(int)round((1.-cent)*NCPU);
	if(avail>NCPU-nrunning){
		avail=NCPU-nrunning;
	}
	if(avail<0) avail=0;
	/*dbg("CPU is %.1f%% Busy. %d running jobs. %d available.\n",cent*100, nrunning, avail); */
	return avail;
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

