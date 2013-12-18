/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <ctype.h>
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
#include "scheduler.h"
#include "scheduler_client.h"
#include "thread.h"

/**
   Contains routines that will be used to talk to the scheduler. The following usage are supported:

   1) by maos or skyc to ask permission to go or report status.
   2) by monitor to show process information or request plotting.
   3) by maos or skyc to print backtrace.

*/

int is_scheduler=0;
#ifndef MAOS_DISABLE_SCHEDULER
#define MAOS_DISABLE_SCHEDULER 0
#endif
#if MAOS_DISABLE_SCHEDULER
int scheduler_start(char *path, int nthread, int waiting){
    (void)path;
    (void)nthread;
    (void)waiting;
    return -1;
}
int scheduler_wait(void){
    return -1;
}
int scheduler_finish(int status){
    (void)status; 
    return -1;
}
int scheduler_report(STATUS_T *status){
    (void)status; return -1;
}
int scheduler_listen(void(*fun)(int)){
    return -1;
}
int scheduler_launch_exe(const char *host, int argc, const char *argv[]){
    (void)host;
    (void)argc;
    (void)argv;
    return -1;
}
int scheduler_send_socket(int sfd){
    return -1;
}
int scheduler_recv_socket(int *sfd){
    return -1;
}
#else
uint16_t PORT=0;
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
    nhost=0;
    snprintf(fn,PATH_MAX,"%s/.aos/hosts",HOME);
    if(exist(fn)){
	nhost=64;
	hosts=malloc(nhost*sizeof(char*));
	FILE *fp=fopen(fn,"r");
	int ihost=0;
	if(fp){
	    char line[64];
	    while(fscanf(fp,"%s\n",line)==1){
		if(strlen(line)>0 && line[0]!='#'){
		    hosts[ihost]=strdup0(line);
		    ihost++;
		    if(ihost>=nhost){
			nhost*=2;
			hosts=realloc(hosts,nhost*sizeof(char*));
		    }
		}
	    }
	    fclose(fp);
	    hosts=realloc(hosts,(ihost+1)*sizeof(char*));
	    nhost=ihost;
	}else{
	    error("failed to open file %s\n",fn);
	}
    }
    hid=myhostid(myhostname());
    if(hid==-1){
	hosts=realloc(hosts, sizeof(char*)*(nhost+1));
	hosts[nhost]=strdup0(myhostname());/*use local machine */
	hid=nhost;
	nhost++;
    }
    register_deinit(NULL,hosts);
}

/**
   Launch the scheduler. We already obtained singleton lock and is in a forked process.
*/
static void launch_scheduler_do(void *junk){
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
	}
    }
}

static void launch_scheduler(void){
    char lockpath[PATH_MAX];
    snprintf(lockpath,PATH_MAX,"%s",TEMP);
    /*launch scheduler if it is not already running. */
    single_instance_daemonize(lockpath,"scheduler", scheduler_version,
			      (void(*)(void*))launch_scheduler_do,NULL);
}

/**
   To open a port and connect to scheduler in the local host*/
static int scheduler_connect_self(int block){
    char fn[PATH_MAX];
    if(TEMP[0]=='/'){
	snprintf(fn, PATH_MAX, "%s/scheduler", TEMP);
    }else{
	snprintf(fn, PATH_MAX, "localhost");
    }
    int sock=connect_port(fn, PORT, 0, 0);
    if(sock<0 && block){
	/*start the scheduler if it is not running*/
	launch_scheduler();
	sock=connect_port(fn, PORT, block, 0);
    }
    return sock;
}

static int psock=-1;

static void scheduler_report_path(char *path){
    static char *path_save=NULL;
    if(psock==-1){
	return;
    }
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
#define CATCH_ERR(A) if(A){psock=-1; return -1;}

/**
   Started by maos to listen to the sock which connects to the
   scheduler for commands
*/
int scheduler_listen(void(*fun)(int)){
    if(psock!=-1 && fun){
	fun(psock);
	return 0;
    }else{
	return -1;
    }
}

/**
   Called by maos to report a job start to scheduler.
 */
int scheduler_start(char *path, int nthread, int waiting){
    psock=scheduler_connect_self(1);
    if(psock==-1){
	warning3("Failed to connect to scheduler\n");
	exit(0);
	return -1;
    }
    scheduler_report_path(path);
    int cmd[4];
    cmd[0]=CMD_START;
    cmd[1]=getpid();
    cmd[2]=nthread;
    cmd[3]=waiting;
    CATCH_ERR(stwriteintarr(psock,cmd,4));
    return 0;
}

/**
   Called by maos to wait for go signal from scheduler.
*/
int scheduler_wait(void){
    if(psock==-1){
	warning3("Failed to connect to scheduler\n");
	return -1;
    }
    /*read will block until clearance is received. */
    int cmd;
    if(streadint(psock, &cmd)){
	warning("Failed to get answer from scheduler.\n");
	return -1;
    }else{
	return 0;
    }
}
static int done=0;
/**
   Called by maos to notify scheduler the completion of a job */
int scheduler_finish(int status){
    static int last_status=0;
    if(!done){
	done=1;
    }else{
	warning("scheduler_finish is called with status %d after %d\n", 
		status, last_status);
	return 1;
    }
    last_status=status;
    if(psock==-1){
	psock=scheduler_connect_self(0);
	scheduler_report_path(NULL);
	if(psock==-1) return -1;
    }
    int cmd[2];
    if(status==0){
	cmd[0]=CMD_FINISH;
    }else{
	cmd[0]=CMD_CRASH;
    }
    cmd[1]=getpid();
    CATCH_ERR(stwriteintarr(psock,cmd,2));
    close(psock);psock=-1;
    return 0;
}

/**
   called by sim.c to report job status */
int scheduler_report(STATUS_T *status){
    if(done){
	warning("scheduler_report called after scheduler_finish\n");
	return 0;
    }	
    if(psock==-1){
	psock=scheduler_connect_self(0);
	scheduler_report_path(NULL);
	if(psock==-1) return -1;
    }
    int cmd[2];
    cmd[0]=CMD_STATUS;
    cmd[1]=getpid();
    CATCH_ERR(stwriteintarr(psock,cmd,2));
    CATCH_ERR(stwrite(psock,status,sizeof(STATUS_T)));
    return 0;
}

/**
   Ask scheduler in another machine to launch exe.
*/
int scheduler_launch_exe(const char *host, int argc, const char *argv[]){
    int ret=0;
    int sock=connect_port(host, PORT, 0, 0);
    if(sock<=-1){
	warning2("Failed to connect to %s:%d: %s\n", host, PORT, strerror(errno));
	return -1;
    }
    int cmd[2]={CMD_LAUNCH, 2};
    char *scmd=argv2str(argc, argv, " ");
    if(stwriteintarr(sock, cmd, 2)
       || stwritestr(sock, argv[0]) 
       || stwritestr(sock, scmd) 
       || streadint(sock, &ret)){
	warning2("Failed to write to scheduler at %s\n", host);
	ret=-1;
    }
    free(scmd);
    close(sock);
    return ret;
}

/**
   send a sock to the scheduler for caching
*/
int scheduler_send_socket(int sfd){
    int ans=-1;
    int ssock=scheduler_connect_self(0);
    if(ssock!=-1 && sfd!=-1){
	int cmd[2]={CMD_SOCK, 1};
	info("send socket %d to scheduler\n", sfd);
	if(stwriteintarr(ssock, cmd, 2) || stwritefd(ssock, sfd)){
	    ans=-1;
	    warning("Talk to scheduler failed\n");
	}else{
	    ans=0;
	}
	close(ssock);
    }
    return ans;
}
/**
   get a socket from the scheduler for reuse
*/
int scheduler_recv_socket(int *sfd){
    int ans=-1;
    int ssock=scheduler_connect_self(0);
    if(ssock!=-1){
	int cmd[2]={CMD_SOCK, -1};
	int ans2=-1;
	if(stwriteintarr(ssock, cmd, 2) || streadint(ssock, &ans2)){
	    warning("Talk to scheduler failed\n");
	}else if(!ans2 && !streadfd(ssock, sfd)){
	    ans=0;
	    info("received %d from scheduler\n", *sfd);
	}else{
	    info("scheduler had no valid fd\n");
	}
	close(ssock);
    }
    return ans;
}
#endif /*MAOS_DISABLE_SCHEDULER*/

/**
   Execute addr2line as specified in buf, combine the answer, and return the
   string.*/
char* call_addr2line(const char *buf){
    FILE *fpcmd=popen(buf,"r");
    char *out=NULL;
    if(!fpcmd){ 
	out=strdup("Command failed\n");
    }else{
	char ans[4096];
	while(fgets(ans, sizeof(ans), fpcmd)){
	    char *tmp=strrchr(ans,'/');
	    if(tmp){
		tmp++;
		char *tmp2=strchr(tmp,'\n'); tmp2[0]='\0';
		if(out){
		    out=stradd(out, "->", tmp, NULL);
		}else{
		    out=strdup(tmp);
		}
	    }
	}   
	pclose(fpcmd);
    }
    return out;
}

/**
   Convert backtrace address to source line.
 */
void print_backtrace_symbol(void *const *buffer, int size){
    static int connect_failed=0;
    char *cmdstr=NULL;
    char add[24];
    char *progname=get_job_progname(0);/*don't free pointer. */
    if(!progname){
	warning("Unable to get progname\n");
	return;
    }
    cmdstr=stradd("addr2line -f -e", progname, NULL);
    free(progname);
    for(int it=size-1; it>0; it--){
	snprintf(add,24," %p",buffer[it]);
	char *tmp=cmdstr;
	cmdstr=stradd(tmp, add, NULL);
	free(tmp);
    }
    if(connect_failed){
	info("%s\n", cmdstr);
	free(cmdstr);
	return;
    }
#if (_POSIX_C_SOURCE >= 2||_XOPEN_SOURCE||_POSIX_SOURCE|| _BSD_SOURCE || _SVID_SOURCE) && !defined(__CYGWIN__)
    PNEW(mutex);//Only one thread can do this.
    LOCK(mutex);
    if(MAOS_DISABLE_SCHEDULER || is_scheduler){
	info("backtrace directly\n");
	char *ans=call_addr2line(cmdstr);
	info2("%s\n", ans);
	free(ans);
    }else{//Create a new socket and ask scheduler to do popen and return answer.
#if MAOS_DISABLE_SCHEDULER == 0
	//Create a new connection.
	int sock=scheduler_connect_self(0);
	if(sock!=-1){
	    int cmd[2];
	    cmd[0]=CMD_TRACE;
	    cmd[1]=getpid();
	    char *ans=NULL;
	    if(stwrite(sock, cmd, sizeof(int)*2)){
		warning2("write cmd %d failed\n", cmd[0]);
	    }else if(stwritestr(sock,cmdstr)){
		warning2("write cmd %s failed\n", cmdstr);
	    }else if(streadstr(sock, &ans)){
		warning2("read cmd failed\n");
	    }else{
		info2(" %s\n",ans); free(ans);
	    }
	    close(sock);
	}else{
	    warning("Failed to connect to scheduler\n");
	    connect_failed=1;
	}
#else
	info2("MAOS_DISABLE_SCHEDULER is no 0\n");
#endif
    }
    UNLOCK(mutex);
#else
    info2("Please call manually: %s\n",cmdstr);
#endif
    sync();
    free(cmdstr);
}
#if !defined(__CYGWIN__) && !defined(__FreeBSD__) && !defined(__NetBSD__)
#include <execinfo.h>
void print_backtrace(){
    int size0,size1;
    size0=1024;
    void *buffer[size0];
    size1=backtrace(buffer,size0);
    print_backtrace_symbol(buffer,size1);
    sync();
}
#else
void print_backtrace(){
}
#endif
