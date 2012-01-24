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
/* make a client address */
#include <stdio.h>
#include <stdlib.h>
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h> 
#include <errno.h>
#include <arpa/inet.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/socket.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <netinet/in.h>
#include <limits.h>
#include <string.h>
#include "sockio.h"
#include "io.h"
#include "process.h"
#include "daemonize.h"
#include "common.h"
#include "misc.h"
#include "scheduler_server.h"
#include "scheduler_client.h"

static int psock;
static int scheduler_crashed;
void scheduler_shutdown(int *sock, int mode){
   if(mode>0){/*tell the server to shutdown read or write */
	int cmd[2];
	if(mode==1){
	    cmd[0]=CMD_SHUTRD;
	}else if(mode==2){
	    cmd[0]=CMD_SHUTWR;
	}
	cmd[1]=getpid();
	stwriteintarr(sock, cmd, 2);
    }
}
int init_sockaddr (struct sockaddr_in *name,
		   const char *hostname, uint16_t port){
    struct hostent *hostinfo;
    
    name->sin_family = AF_INET;
    name->sin_port = htons(port);
    hostinfo = gethostbyname (hostname);
    if (hostinfo == NULL){
	/*perror("gethostbyname"); */
	/*fprintf (stderr, "Unknown host %s.\n", hostname); */
	return -1;
    }else{
	struct in_addr *addr = (struct in_addr *) hostinfo->h_addr_list[0];
	if(addr){
	    name->sin_addr = *addr;
	    return 0;
	}else{
	    /*warning("h_addr_list is NULL for host %s\n", hostname); */
	    return -1;
	}
    }
}
/* To open a port and connect to scheduler */
int scheduler_connect_self(int block, int mode){
    /*
      mode=0: read/write
      mode=1: read only by the client. the server won't read
      mode=2: write only by the client. the server won't write
     */
    int sock;
 
    if(scheduler_crashed) {
	return -1;
    }
    int count=0;
    struct sockaddr_in servername;
    do{
	/* Create the socket. */
	sock = socket (PF_INET, SOCK_STREAM, 0);
	if (sock < 0) {
	    perror ("socket (scheduler)");
	    scheduler_crashed=1; 
	    return sock;
	}
	cloexec(sock);
	socket_tcp_keepalive(sock);
	/* Give the socket a name. */
	init_sockaddr(&servername, "localhost", PORT);
	if(connect(sock, (struct sockaddr *)&servername, sizeof (servername))<0){
	    perror("connect");
	    close(sock);
	    if(!block){
		return -1;
	    }
	    sleep(4);
	    count++;
	}else{
	    scheduler_shutdown(&sock,mode);
	    return sock;
	}
    }while(count<10);
    return -1;
}

static char *path_save=NULL;
static void scheduler_report_path(char *path){
    if(path){
	if(path_save) free(path_save);
	path_save=strdup(path);
    }else{
	if(!path_save){
	    path_save=strdup("unknown");
	}
    }
    int cmd[2];
    cmd[0]=CMD_PATH;
    cmd[1]=getpid();
    stwriteintarr(&psock, cmd, 2);
    stwritestr(&psock,path_save);
}
void scheduler_deinit(void){
    free(path_save);
}
static __attribute__((constructor)) void init(){
    register_deinit(scheduler_deinit, NULL);
}
/* called by mcao to wait for available cpu. */
int scheduler_start(char *path, int nthread, int waiting){
    psock=scheduler_connect_self(1,0);
    if(psock==-1){
	warning3("Failed to connect to scheduler\n");
	return -1;
    }
    scheduler_report_path(path);
    int cmd[2];
    cmd[0]=CMD_START;
    cmd[1]=getpid();
    stwriteintarr(&psock,cmd,2);
    cmd[0]=nthread;
    cmd[1]=waiting;
    stwriteintarr(&psock,cmd,2);
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
	psock=scheduler_connect_self(0,2);
	scheduler_report_path(NULL);
    }
    if(psock==-1) return;
    int cmd[2];
    if(status==0)
	cmd[0]=CMD_FINISH;
    else 
	cmd[0]=CMD_CRASH;
    cmd[1]=getpid();
    stwriteintarr(&psock,cmd,2);
    close(psock);psock=-1;
}
/* called by sim.c to report job status */
void scheduler_report(STATUS_T *status){
    if(psock==-1){
	psock=scheduler_connect_self(0, 2);
	scheduler_report_path(NULL);
    }
    if(psock==-1) return;
    int cmd[2];
    cmd[0]=CMD_STATUS;
    cmd[1]=getpid();
    stwriteintarr(&psock,cmd,2);
    stwrite(&psock,status,sizeof(STATUS_T));
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
    if(psock==-1)
	psock=scheduler_connect_self(0,0);
    if(psock==-1) return;
    int cmd[2];
    cmd[0]=CMD_TRACE;
    cmd[1]=getpid();
    stwrite(&psock,cmd,sizeof(int)*2);
    stwritestr(&psock,cmdstr);
    char *ans=streadstr(&psock);
    fprintf(stderr, " %s\n",ans);
    free(ans);
#else
    fprintf(stderr, " %s\n",cmdstr);
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
	if(direct){
	    /*launch directly, used by drawres, drawbin where overhead is small. */
	    scheduler_launch_drawdaemon(fifo);
	}else{
	    /*launch indirectly through drawdaemon by maos. Slow to launch
	      directly because maos is using a lot of memory and forking it
	      takes time*/
	    int sock;
	    for(int retry=0; retry<10; retry++){
		sock=scheduler_connect_self(0,0);
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
	    stwrite(&sock,cmd,sizeof(int)*2);
	    /*make sure the drawdaemon appears in our DISPLAY. */
	    const char *display=getenv("DISPLAY");
	    if(strlen(display)==0){
		warning("There is no DISPLAY\n");
		return NULL;
	    }
	    const char *xauth=getenv("XAUTHORITY");
	    stwritestr(&sock,display);
	    stwritestr(&sock,xauth);
	    stwritestr(&sock,fifo);
	    if(stread(&sock,cmd,sizeof(int))) return NULL;
	    if(sock!=-1) close(sock);
	    if(cmd[0]==-1) return NULL;/*failed */
	}
	sleep(1);/*wait for drawdaemon to start. */
    }
    return fifo;
}
