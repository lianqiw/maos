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
/* make a client address */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h> 
#include <errno.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/socket.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netinet/tcp.h> /*SOL_TCP */
#include <limits.h>
#include <string.h>
#include "sock.h"
#include "sockio.h"
#include "process.h"
#include "daemonize.h"
#include "common.h"
#include "misc.h"
#include "hashlittle.h"
#include "scheduler_server.h"
#include "scheduler_client.h"
#include "thread.h"
#ifdef MAOS_DISABLE_SCHEDULER
int scheduler_start(char *path, int nthread, int waiting){
  (void)path;
  (void)nthread;
  (void)waiting;
   return 0;
}
int scheduler_wait(void){
   return 0;
}
void scheduler_finish(int status){
  (void)status;
}
void scheduler_report(STATUS_T *status){
  (void)status;
}
void print_backtrace_symbol(void *const *buffer, int size){
  (void) buffer;
  (void)size;
}
void print_backtrace(int sig){
  (void) sig;
}

#else
void scheduler(void);
char *scheduler_fnlog=NULL;

uint16_t PORT=0;
uint16_t PORTMON=0;
char** hosts;
int nhost;
int hid;
static int myhostid(const char *host){
    int i;
    for(i=0; i<nhost; i++){
	if(!strncasecmp(hosts[i],host,strlen(hosts[i])))
	    break;
    }
    if(i==nhost){
	i=-1;
    }
    return i;
}


/*Initialize hosts and associate an id number */
static __attribute__((constructor))void init(){
    init_process();/*the constructor in process.c may not have been called. */
    char fn[PATH_MAX];
    snprintf(fn,PATH_MAX,"%s/.aos/jobs.log", HOME);
    scheduler_fnlog=strdup0(fn);
    snprintf(fn,PATH_MAX,"%s/.aos/port",HOME);
    PORT=0;
    {
	FILE *fp=fopen(fn,"r");
	if(fp){
	    if(fscanf(fp,"%hu", &PORT)!=1){
		/*warning3("Failed to read port from %s\n",fn); */
	    }
	    fclose(fp);
	}
    }
    if(PORT==0){
	/*user dependent PORT to avoid conflict */
	PORT= (uint16_t)((uint16_t)(hashlittle(USER,strlen(USER),0)&0x2FFF)|10000);
    }
    snprintf(fn,PATH_MAX,"%s/.aos/hosts",HOME);
    if(!exist(fn)){
	nhost=1;
	hosts=malloc(nhost*sizeof(char*));
	hosts[0]=calloc(60,sizeof(char));
	gethostname(hosts[0],60);
	register_deinit(NULL,hosts[0]);
    }else{
	nhost=64;
	hosts=malloc(nhost*sizeof(char*));
	FILE *fp=fopen(fn,"r");
	int ihost=0;
	if(fp){
	    char line[64];
	    while(fscanf(fp,"%s\n",line)==1){
		if(strlen(line)>0){
		    hosts[ihost]=strdup0(line);
		    ihost++;
		    if(ihost>=nhost){
			nhost*=2;
			hosts=realloc(hosts,nhost*sizeof(char*));
		    }
		}
	    }
	    fclose(fp);
	    hosts=realloc(hosts,ihost*sizeof(char*));
	    nhost=ihost;
	}else{
	    error("failed to open file %s\n",fn);
	}
    }
    register_deinit(NULL,hosts);

    char host[60];
    if(gethostname(host,60)) warning3("Unable to get hostname\n");
    hid=myhostid(host);
    if(hid==-1){
	hosts[nhost]=strdup0(host);/*use local machine */
	hid=nhost;
	nhost++;
	if(hid==-1){
	    warning3("Unable to determine proper hostname. Monitor may not work\n");
	}
    }
}

/**
   Launch the scheduler. We already obtained singleton lock and is in a forked process.
 */
static void scheduler_launch_do(void *junk){
    (void)junk;
#if defined(__CYGWIN__)
    char *fn_scheduler=stradd(BUILDDIR, "/bin/scheduler.exe", NULL);
#else
    char *fn_scheduler=stradd(BUILDDIR, "/bin/scheduler", NULL);
#endif
    if(exist(fn_scheduler)){
	info2("Run %s\n", fn_scheduler);
	execl(fn_scheduler, "scheduler", NULL);
    }else{/*fall back. this won't have the right argv set. */
	info2("Launch scheduler using shell\n");
	if(execlp("scheduler", "scheduler", NULL)){
	    warning("scheduler not found\n");
	    scheduler();
	}
    }
}

static void scheduler_launch(void){
    char lockpath[PATH_MAX];
    snprintf(lockpath,PATH_MAX,"%s",TEMP);
    /*launch scheduler if it is not already running. */
    single_instance_daemonize(lockpath,"scheduler", scheduler_version,
			      (void(*)(void*))scheduler_launch_do,NULL);
}

/**
   To open a port and connect to scheduler in the local host*/
static int scheduler_connect_self(int block){
    /*start the scheduler if it is not running*/
    int sock=connect_port("localhost", PORT, 0, 0);
    if(sock<0){
	scheduler_launch();
	sleep(1);
        sock=connect_port("localhost", PORT, block, 0);
    }
    if(sock<0){
	warning2("Unable to connect to port %d\n", PORT);
    }
    return sock;
}

static int psock;
static char *path_save=NULL;
static void scheduler_report_path(char *path){
    if(path){
	if(path_save){
	    path_save=realloc(path_save, strlen(path)+1);
	    strcpy(path_save, path);
	}else{
	    path_save=strdup(path);
	    register_deinit(NULL, path_save);
	}
    }else{
	if(!path_save){
	    path_save=strdup("unknown");
	    register_deinit(NULL, path_save);
	}
    }
    int cmd[2];
    cmd[0]=CMD_PATH;
    cmd[1]=getpid();
    stwriteintarr(psock, cmd, 2);
    stwritestr(psock,path_save);
}
/* called by mcao to wait for available cpu. */
int scheduler_start(char *path, int nthread, int waiting){
    psock=scheduler_connect_self(1);
    if(psock==-1){
	warning3("Failed to connect to scheduler\n");exit(0);
	return -1;
    }
    scheduler_report_path(path);
    int cmd[2];
    cmd[0]=CMD_START;
    cmd[1]=getpid();
    stwriteintarr(psock,cmd,2);
    cmd[0]=nthread;
    cmd[1]=waiting;
    stwriteintarr(psock,cmd,2);
    return 0;
}

int scheduler_wait(void){
    if(psock==-1){
	warning3("Failed to connect to scheduler\n");
	return -1;
    }
    /*read will block until clearance is received. */
    int cmd[2];
    if(read(psock, cmd, sizeof(int))==sizeof(int)){
	/*info2("Scheduler replied %d.\n",cmd[0]); */
	return 0;
    }else{
	warning("Failed to get answer from scheduler.\n");
	return -1;
    }
    /*don't close socket. */
}
/* called by mcao to notify scheduler the completion of a job */
void scheduler_finish(int status){
    if(psock==-1){
	psock=scheduler_connect_self(0);
	scheduler_report_path(NULL);
    }
    if(psock==-1) return;
    int cmd[2];
    if(status==0)
	cmd[0]=CMD_FINISH;
    else 
	cmd[0]=CMD_CRASH;
    cmd[1]=getpid();
    stwriteintarr(psock,cmd,2);
    close(psock);psock=-1;
}
/* called by sim.c to report job status */
void scheduler_report(STATUS_T *status){
    if(psock==-1){
	psock=scheduler_connect_self(0);
	scheduler_report_path(NULL);
    }
    if(psock==-1) return;
    int cmd[2];
    cmd[0]=CMD_STATUS;
    cmd[1]=getpid();
    stwriteintarr(psock,cmd,2);
    stwrite(psock,status,sizeof(STATUS_T));
    /*don't close socket. */
}

/*!defined(__INTEL_COMPILER)||1)  */
#if (_POSIX_C_SOURCE >= 2||_XOPEN_SOURCE||_POSIX_SOURCE|| _BSD_SOURCE || _SVID_SOURCE) && !defined(__CYGWIN__)
#define PRINTBACKTRACE 1
#else
#define PRINTBACKTRACE 0
#endif

void print_backtrace_symbol(void *const *buffer, int size){
    char cmdstr[BACKTRACE_CMD_LEN];
    char add[24];
    const char *progname=get_job_progname();/*don't free pointer. */
    if(!progname){
	error("Unable to get progname\n");
    }
#if PRINTBACKTRACE == 1 
    snprintf(cmdstr,BACKTRACE_CMD_LEN,"addr2line -f -e %s",progname);
#else
    snprintf(cmdstr,BACKTRACE_CMD_LEN,"%s: ",progname);
#endif
    int it;
    for(it=size-1; it>-1; it--){
	snprintf(add,24," %p",buffer[it]);
	strncat(cmdstr,add,BACKTRACE_CMD_LEN-strlen(cmdstr)-1);
    }
#if PRINTBACKTRACE == 1 
    PNEW(mutex);//Only one thread can do this.
    LOCK(mutex);
    if(psock==-1)
	psock=scheduler_connect_self(0);
    if(psock==-1) goto end;
    int cmd[2];
    cmd[0]=CMD_TRACE;
    cmd[1]=getpid();
    char *ans;
    if(stwrite(psock,cmd,sizeof(int)*2) || stwritestr(psock,cmdstr) || streadstr(psock, &ans)){
	goto end;
    }
    info2(" %s\n",ans);
    free(ans);
 end:
    UNLOCK(mutex);
#else
    info2(" %s\n",cmdstr);
#endif
}
#if !defined(__CYGWIN__) && !defined(__FreeBSD__) && !defined(__NetBSD__)
#include <execinfo.h>
void print_backtrace(int sig){
    int size0,size1;
    if(sig !=0)
	info2("Segmentation %d fault caught. print backtrace\n",sig);
    size0=1024;
    void *buffer[size0];
    size1=backtrace(buffer,size0);
    print_backtrace_symbol(buffer,size1);
    if(sig !=0)
	raise(SIGABRT);
}
#endif
#endif

/**
   Fork and launch drawdaemon.
*/
int scheduler_launch_drawdaemon(char *fifo){
    int method=0;
#if defined(__APPLE__) && 0
    char cmdopen[1024];
    /*Run the exe directly can pass the argumnents. --args is a new feature in 10.6 to do the samething with open */
    snprintf(cmdopen, 1024, "%s/scripts/drawdaemon.app/Contents/MacOS/drawdaemon %s &", SRCDIR, fifo);
    if(system(cmdopen)){
	method=0;/*failed */
	warning("%s failed\n", cmdopen);
    }else{/*succeed */
	info("%s succeeded\n", cmdopen);
	method=3;
    }
    if(method==0){
	snprintf(cmdopen, 1024, "open -n -a drawdaemon.app --args %s", fifo);
	if(system(cmdopen)){
	    warning("%s failed\n", cmdopen);
	    method=0;/*failed */
	}else{
	    info("%s succeeded\n", cmdopen);
	    method=3;
	}
    }
#endif
    char *fn=stradd(BUILDDIR, "/bin/drawdaemon",NULL);
    if(method==0){
	info2("Looking for drawdaemon in %s\n",fn);
	if(exist(fn)){
	    info2("Found drawdaemon in %s, run it.\n",fn);
	    method=1;
	}else{
	    warning3("Not found drawdaemon in %s, use bash to find and run drawdaemon.\n",fn);
	    int found=!system("which drawdaemon");
	    if(found){
		method=2;
	    }else{
		warning3("Unable to find drawdaemon\n");
	    }
	}
    }
    int ans;
    if(method==0){
	ans=1;/*failed */
    }else{
	ans=0;/*succeed */
    }
    if(method==3){
	return ans;
    }
    /*Now start to fork. */
    int pid2=fork();
    if(pid2<0){
	warning3("Error forking\n");
    }else if(pid2>0){
	/*wait the child so that it won't be a zoombie */
	waitpid(pid2,NULL,0);
	return ans;
    }
    pid2=fork();
    if(pid2<0){
	warning3("Error forking\n");
	_exit(EXIT_FAILURE);
    }else if(pid2>0){
	_exit(EXIT_SUCCESS);/*waited by parent. */
    }
    /*safe child. */
    setsid();
    fclose(stdin);
    if(method==1){
	if(execl(fn, "drawdaemon",fifo,NULL)){
	    perror("execl");
	    warning("execl failed\n");
	}
    }else if(method==2){
	if(execlp("drawdaemon","drawdaemon",fifo,NULL)){
	    warning("execlp failed\n");
	}
    }else{
	error("Invalid method.\n");
    }
    return ans;
}
/**
   Get the drawdaemon.  if direct=1, do not go through the scheduler_server
(useful for drawres, drawbin), otherwise, go through the scheduler_server(
useful for MAOS because forking a big program is expensive) */
char* scheduler_get_drawdaemon(int pid, int direct){
    int launch=0;
    static char *fifo=NULL;
    if(!fifo){
        fifo=malloc(100);
        snprintf(fifo,100,"%s/drawdaemon_%d.fifo",TEMP,pid);
	fifo=realloc(fifo,strlen(fifo)+1);
    }
    if(exist(fifo)){
        /*warning2("fifo already exist. test when drawdaemon exists\n"); */
        char fnpid[PATH_MAX];
	snprintf(fnpid, PATH_MAX, "%s/drawdaemon_%d.pid", TEMP, pid);
	FILE *fp=fopen(fnpid, "r");
	if(fp){
	    int fpid;
	    if(fscanf(fp, "%d", &fpid)!=1){
		fpid=-1;/*failed to read fpid. */
		launch=1;
	    }else{
		if(kill(fpid,0)){
		    /*warning2("Drawdaemon has exited\n"); */
		    launch=1;
		}
	    }
	    fclose(fp);
	}else{
	    /*warning2("Drawdaemon has exited\n"); */
	    launch=1;
	}
    }else{
        /*info2("make fifo\n"); */
	if(mkfifo(fifo,0700)){
	    warning3("Error making fifo\n");
	}
	launch=1;
    }
    if(launch){
#ifdef MAOS_DISABLE_SCHEDULER
	scheduler_launch_drawdaemon(fifo);
#else
	if(direct){
	    /*launch directly, used by drawres, drawbin where overhead is small. */
	    scheduler_launch_drawdaemon(fifo);
	}else{
	    /*launch indirectly through drawdaemon by maos. Slow to launch
	      directly because maos is using a lot of memory and forking it
	      takes time*/
	    int sock;
	    for(int retry=0; retry<10; retry++){
		sock=scheduler_connect_self(0);
		if(sock==-1){
		    warning2("failed to connect to scheduler\n");
		    sleep(1);
		}else{
		    break;
		}
	    }
	    if(sock==-1){
		warning2("failed to connect to scheduler\n");
		return NULL;
	    }
	    int cmd[2];
	    cmd[0]=CMD_DRAW;
	    cmd[1]=pid;
	    stwrite(sock,cmd,sizeof(int)*2);
	    /*make sure the drawdaemon appears in our DISPLAY. */
	    const char *display=getenv("DISPLAY");
	    if(strlen(display)==0){
		warning("There is no DISPLAY\n");
		return NULL;
	    }
	    const char *xauth=getenv("XAUTHORITY");
	    stwritestr(sock,display);
	    stwritestr(sock,xauth);
	    stwritestr(sock,fifo);
	    if(stread(sock,cmd,sizeof(int))) return NULL;
	    if(sock!=-1) close(sock);
	    if(cmd[0]==-1) return NULL;/*failed */
	}
#endif
	sleep(1);/*wait for drawdaemon to start. */
    }
    return fifo;
}

