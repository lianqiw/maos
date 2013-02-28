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
/**
   Contains routines to do job scheduling and monitoring.
   
   Change log: 

   2010-01-17: Removed communicate between schedulers on different
   machines. This is complicated and not reliable. It is up to the monitor to
   connect all the servers to get status information.
   
   2010-07-26: Split scheduler to scheduler_server and scheduler_client.

   socket programming guideline: 

   1) Send heartbeat signal (TCP keepalive does this) to detect link broken
   
   2) Do not pass host id around, as it might be defined differently for the same machine.

   2012-10-25: This file only contains the routines to be used by the server.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "../sys/sys.h"

extern char *scheduler_fnlog;
/**
   Struct to hold information of running jobs.
*/
typedef struct RUN_T{
    struct RUN_T *next;
    STATUS_T status;
    double started;/*started execution. */
    double launchtime;
    int pid;
    int sock;
    int nthread;
    int time;
    char *cwd; /*Current job path*/
    char *exe; /*Path to the executable.*/
    char *path;/*Job arguments.*/
}RUN_T;

/**
  Struct to hold available monitors waiting for information.
*/
typedef struct MONITOR_T{
    int sock;
    int load;/*handle machine load information. */
    struct MONITOR_T *next;
}MONITOR_T;


static MONITOR_T *pmonitor=NULL;
static int all_done=0;
/*A linked list to store running process*/
static RUN_T *running=NULL;/*points to the begining */
static RUN_T *running_end=NULL;/*points to the last to add. */
/*A linked list to store ended process*/
static RUN_T *runned=NULL;
static RUN_T *runned_end=NULL;
static int nrun=0;
static double usage_cpu;
static RUN_T* running_add(int pid,int sock);
static RUN_T *running_get(int pid);
static RUN_T *running_get_wait(int status);

static RUN_T *runned_get(int pid);
static void runned_remove(int pid);
static int runned_add(RUN_T *irun);
static void running_remove(int pid);
static void running_update(int pid, int status);
static RUN_T *running_get_by_sock(int sock);
//static MONITOR_T *monitor_get(int hostid);
static void monitor_remove(int hostid);
static MONITOR_T *monitor_add(int hostid);
static void monitor_send(RUN_T *run,char*path);
static void monitor_send_initial(MONITOR_T *ic);
static void monitor_send_load(void);

/**
   The following runned_* routines maintains the linked list
   "runned" which contains the already finished jobs.
*/
static int runned_add(RUN_T *irun){
    /*
      Add to the end of the runned list intead of the beggining.
     */
    if(runned_get(irun->pid)){
	warning("runned_add: already exist\n");
	return 1;
    }else{
	/*move irun to runned list from running list. */
	irun->next=NULL;
	if(runned_end){
	    runned_end->next=irun;
	    runned_end=irun;
	}else{
	    runned_end=runned=irun;
	}
	irun->status.done=1;
	if(irun->status.info<10){
	    warning3("Should be a finished process\n");
	    return 1;
	}else{
	    return 0;
	}
    }
}
static void runned_remove(int pid){
    RUN_T *irun,*irun2=NULL;
    int removed=0;
    for(irun=runned; irun; irun2=irun,irun=irun->next){
	if(irun->pid==pid){
	    irun->status.info=S_REMOVE;
	    monitor_send(irun,NULL);
	    if (irun2){
		irun2->next=irun->next;
		if(irun->next==NULL)
		    runned_end=irun2;
	    }else{
		runned=irun->next;
		if(irun->next==NULL)
		    runned_end=runned;
	    }
	    free(irun->path);
	    free(irun->exe);
	    free(irun->cwd);
	    free(irun);
	    removed=1;
	    break;
	}
    }
    if(!removed){
	warning3("runned_remove: Record %s:%d not found!\n",hosts[hid],pid);
    }
}
static RUN_T *runned_get(int pid){
    RUN_T *irun;
    for(irun=runned; irun; irun=irun->next){
	if(irun->pid==pid){
	    break;
	}
    }
    return irun;
}
/*
static RUN_T *waiting_add(int pid){
    RUN_T *irun=calloc(1, sizeof(RUN_T));
    irun->pid=pid;
    if(waiting_end){
	waiting_end->next=irun;
	waiting_end=irun;
    }else{
	waiting=irun;
	waiting_end=irun;
    }
    return irun;
}
static RUN_T *waiting_get_wait(){
    RUN_T *irun=waiting;
    if(waiting){
	waiting=waiting->next;
        if(!waiting){
	    waiting_end=NULL;
	}
    }
    return irun;
    }*/
/**
   The following running_* routiens operates on running linked list
   which contains queued jobs.
*/
static RUN_T* running_add(int pid,int sock){
    RUN_T *irun;
    if(pid && (irun=running_get(pid))){
	return irun;
    }else{
	/*create the node */
	irun=calloc(1, sizeof(RUN_T));
	irun->pid=pid;
	irun->sock=sock;
	/*record the launch time */
	if(pid>0){
	    irun->launchtime=get_job_launchtime(pid);
	}else{
	    irun->launchtime=INFINITY;
	}
	irun->next=NULL;
	if(running_end){/*list is not empty*/
	    if(pid<=0 || running_end->launchtime <= irun->launchtime){
		/*append the node to the end. */
		running_end->next=irun;
		running_end=irun;
	    }else{
		/*insert the node in the middle. */
		RUN_T *jrun,*jrun2=NULL;
		for(jrun=running; jrun; jrun2=jrun,jrun=jrun->next){
		    if(jrun->launchtime >= irun->launchtime){
			irun->next=jrun;
			if(jrun2){
			    jrun2->next=irun;
			}else{
			    running=irun;
			}
			break;
		    }
		}
		if(!jrun){
		    warning("failed to insert pid %d\n",pid);
		}
	    }
	}else{
	    running_end=running=irun;
	}
	irun->status.timstart=myclocki();
	return irun;
    }
}
static void running_update(int pid, int status){
    RUN_T *irun=running_get(pid);
    if(irun){
	if(irun->status.info!=S_WAIT && status>10){
	    nrun-=irun->nthread;/*negative decrease */
	}
	irun->status.info=status;
	monitor_send(irun,NULL);
	if(status>10){/*ended */
	    irun->time=(int)myclockd();
	    if(nrun<0) nrun=0;
	    running_remove(pid);/*moved to runned. */
	}
    }
}
static void running_remove(int pid){
    RUN_T *irun,*irun2=NULL;
    int removed=0;
    for(irun=running; irun; irun2=irun,irun=irun->next){
	if(irun->pid==pid){
	    if (irun2){
		irun2->next=irun->next;
		if(irun->next==NULL)
		    running_end=irun2;
	    }else{
		running=irun->next;
		if(irun->next==NULL)
		    running_end=running;
	    }
	    if(irun->pid>0){
		irun->status.timend=myclocki();
		/*move irun to runned */
		runned_add(irun);
		{
		    const char *statusstr=NULL;
		    switch (irun->status.info){
		    case S_CRASH:
			statusstr="Crashed"; break;
		    case S_FINISH:
			statusstr="Finished";break;
		    case S_KILLED:
			statusstr="Killed";break;
		    default:
			statusstr="Unknown";break;
		    }
		    FILE *fp=fopen(scheduler_fnlog,"a");
		    if(fp){
			fprintf(fp,"[%s] %s %5d %8s '%s' nrun=%d\n",
				myasctime(),hosts[hid],pid,statusstr,irun->path,nrun);
			fclose(fp);
		    }
		}
	    }else{//The process is not running. run it and remove from monitor.
		if(irun->status.info==S_QUEUED){
		    irun->status.info=S_REMOVE;
		    warning("monitor_send S_REMOVE\n");
		    monitor_send(irun, NULL);
		    free(irun->path);
		    free(irun->exe);
		    free(irun->cwd);
		    free(irun);
		}else{//Killed. Add to runned queue
		    runned_add(irun);
		}
	    }
	    removed=1;
	    break;
	}
    }
    if(!removed){
	warning3("runned_remove %s:%d not found\n",hosts[hid],pid);
    }
}

static RUN_T *running_get(int pid){
    RUN_T *irun;
    for(irun=running; irun; irun=irun->next){
	if(irun->pid==pid){
	    break;
	}
    }
    return irun;
}

static RUN_T *running_get_wait(int status){
    RUN_T *irun;
    int jrun=0;
    for(irun=running; irun; irun=irun->next){
	jrun++;
	if(irun->status.info==status){
	    break;
	}
    }
    return irun;
}

static RUN_T *running_get_by_sock(int sock){
    RUN_T *irun;
    for(irun=running; irun; irun=irun->next){
	if(irun->sock==sock){
	    break;
	}
    }
    return irun;
}
static void check_jobs(void){
    /**
       check all the jobs. remove if any job quited.
     */
    RUN_T *irun;
 restart:
    if(running){
	for(irun=running; irun; irun=irun->next){
	    if(irun->pid>0 && kill(irun->pid,0)){
		running_update(irun->pid,S_CRASH);
		/*list is changed. need to restart the loop. */
		goto restart;
	    }
	}
    }
}
/**
 examine the process queue to start routines once CPU is available.
*/
static void process_queue(void){
    static double timestamp=0;
    if(nrun>0 && myclockd()-timestamp<10) {
	return;
    }
    timestamp=myclockd();
    if(nrun>=NCPU) return;
    if(nrun<0){
	warning("process_queue: nrun=%d. set to 0.\n", nrun);
	nrun=0;
    }
    int avail=get_cpu_avail();
    info2("process_queue: nrun=%d avail=%d\n", nrun, avail);
    if(avail<1) return;
    RUN_T *irun=running_get_wait(S_WAIT);
    if(irun){
	int nthread=irun->nthread;
	if(nrun+nthread<=NCPU && (nthread<=avail || avail >=3)){
	    nrun+=nthread;
	    /*don't close the socket. will close it in select loop. */
	    /*warning3("process %d launched. write to sock %d cmd %d\n", */
	    /*irun->pid, irun->sock, S_START); */
	    info2("process_queue: Notify %d at %d\n", irun->pid, irun->sock);
	    if(stwriteint(irun->sock,S_START)){
		perror("stwriteint");
		warning("failed to notify maos\n");
	    }
	    irun->status.timstart=myclocki();
	    irun->status.info=S_START;
	    monitor_send(irun,NULL);
	    FILE *fp=fopen(scheduler_fnlog,"a");
	    if(fp){
		fprintf(fp,"[%s] %s %5d  started '%s' nrun=%d\n",
			myasctime(),hosts[hid],irun->pid,irun->path,nrun);
	    }
	    fclose(fp);
	}
    }else{
	if(avail>1){
	    irun=running_get_wait(S_QUEUED);
	    if(!irun){
		all_done=1;
	    }else{
		info2("process_queue: no more jobs pending, process waiting list\n");
		if(launch_exe(irun->exe, irun->cwd, irun->path)){
		    warning("launch_exe %s failed\n", irun->path);
		}else{
		    warning("launch_exe %s succeed\n", irun->path);
		}
		running_remove(irun->pid);
		sleep(1);//wait for the job to actually start.
	    }
	}
    }
}
/**
   respond to client requests
*/
static int respond(int sock){
    /*
       Don't modify sock in this routine. otherwise select
       will complain Bad file descriptor
    */
    int ret=0, pid, cmd[2];
    info2("\rrespond %d start ... ", sock);
    if((ret=streadintarr(sock, cmd, 2))){
	info2("read failed");
	goto end;
    }
    pid=cmd[1];
    info2("\rrespond %d got %d %d.", sock, cmd[0], cmd[1]);
    switch(cmd[0]){
    case CMD_START://Called by maos when job starts.
	{
	    /* Data read. */
	    int nthread;
	    int waiting;
	    if(streadint(sock, &nthread) || streadint(sock, &waiting)){
		ret=-1;
		break;
	    }
	    if(nthread<1) 
		nthread=1;
	    else if(nthread>NCPU)
		nthread=NCPU;
	    RUN_T *irun=running_add(pid,sock);
	    irun->nthread=nthread;
	    if(waiting){
		irun->status.info=S_WAIT;
		all_done=0;
		info2("respond: %d queued. nrun=%d\n",pid,nrun);
	    }else{/*no waiting, no need reply. */
		info2("scheduler: nrun is increased by %d\n", nthread);
		nrun+=nthread;
		irun->status.info=S_START;
		irun->status.timstart=myclocki();
		info2("respond: %d started. nrun=%d\n",pid,nrun);
	    }
	    if(irun->path) monitor_send(irun,irun->path);
	    monitor_send(irun,NULL);
	}
	break;
    case CMD_FINISH:
	running_update(pid,S_FINISH);
	return -1;
	break;
    case CMD_STATUS:/*called by maos */
	{
	    RUN_T *irun=running_get(pid);
	    int added=0;
	    if(!irun){
		added=1;
		/*started before scheduler is relaunched. */
		irun=running_add(pid,sock);
		irun->status.info=S_START;
	    }
	    if(sizeof(STATUS_T)!=read(sock,&(irun->status), sizeof(STATUS_T))){
		warning3("Error reading\n");
	    }
	    if(added){
		irun->nthread=irun->status.nthread;
		info2("respond, STATUS: nrun is increased by %d\n", irun->nthread);
		nrun+=irun->nthread;
	    }
	    monitor_send(irun,NULL);
	}
	break;
    case CMD_CRASH:/*called by maos */
	running_update(pid,S_CRASH);
	return -1;
	break;
    case CMD_MONITOR:/*called by monitor */
	{
	    MONITOR_T *tmp=monitor_add(sock);
	    if(pid>=0x8){/*check monitor version. */
		tmp->load=1;
	    }
	    info2("respond: Monitor is connected at sock %d.\n", sock);
	}
	break;
    case CMD_PATH:{
	RUN_T *irun=running_add(pid, sock);
	if(irun->path) free(irun->path);
	if(streadstr(sock, &irun->path)){
	    ret=-1;
	    break;
	}
	info2("respond: Received path: %s\n",irun->path);
	monitor_send(irun,irun->path);
    }
	break;
    case CMD_KILL:/*called by mnitor. */
	{
	    RUN_T *irun=running_get(pid);
	    if(irun){
		if(irun->status.info!=S_QUEUED){
		    kill(pid,SIGTERM);
		    sleep(1);
		    if(!kill(pid,0)){
			kill(pid,SIGKILL);
		    }
		}
		running_update(pid,S_KILLED);
	    }
	}
	break;
    case CMD_TRACE://called by maos to print backtrace
	{
	    char *buf=NULL, *out=NULL;
	    if(streadstr(sock, &buf)||!(out=call_addr2line(buf))||stwritestr(sock,out)){
		info2("respond: CMD_TRACE failed. buf=%s, out=%s\n", buf, out);
		ret=-1; 
	    }
	    free(buf);
	    free(out);
	}
	break;
    case CMD_DRAW://not used
	break;
    case CMD_SOCK:
	break;
    case CMD_REMOVE:/*called by monitor to remove a runned object. */
	{
	    RUN_T*irun=runned_get(pid);
	    if(irun){
		runned_remove(pid);
	    }else{
		warning3("CMD_REMOVE: %s:%d not found\n",hosts[hid],pid);
	    }
	}
	break;	  
    case CMD_DISPLAY:/*called by remote display to talk to maos*/
	{
	    RUN_T *irun=running_get(pid);
	    int cmd2[2]={CMD_SOCK, 0};
	    if(stwriteintarr(irun->sock, cmd2, 2) || stwritefd(irun->sock, sock)){
		warning("Unable to pass socket %d to maos at %d\n", sock, irun->sock);
		stwriteint(sock, -1);//respond failure message.
	    }else{
		stwriteint(sock, 0);
	    }
	    ret=-1;
	}
	break;
    case CMD_VERSION://not used
	break;
    case CMD_LOAD://intended for monitor
	break;
    case CMD_UNUSED2://intended for maos
	break;
    case CMD_UNUSED3:
	break;
    case CMD_ADDHOST://only used in monitor
	break;
    case CMD_LAUNCH:{
	/*called by maos/skyc from another machine to start a job in this
	  machine*/
	char *execmd=NULL;
	char *execwd=NULL;
	char *exename=NULL;
	if(pid==3){
	    if(streadstr(sock, &exename) || streadstr(sock, &execwd)){
		warning("Unable to read exename.\n");
	    }
	}
	if(!streadstr(sock, &execmd)){
	    static long counter=0;//an negative index to retrieve this irun.
	    counter--;
	    RUN_T *irun=running_add(counter, -1);
	    irun->status.info=S_QUEUED;
	    irun->exe=exename;
	    irun->cwd=execwd;
	    irun->path=execmd;
	    char *tmp=stradd(irun->cwd, "/", irun->path, NULL);//add path.
	    monitor_send(irun, tmp); free(tmp);
	    monitor_send(irun, NULL);
	    all_done=0;
	    info2("respond: exe=%s\ncwd=%s\npath=%s\n added to waitinglist. nrun=%d\n", irun->exe, irun->cwd, irun->path, nrun);
	    ret=0;
	}else{
	    warning("Unable to read execmd\n");
	    ret=-1;
	}
	if(stwriteint(sock, ret)){
	    ret=-1;
	}
    }
	break;
    default:
	warning3("Invalid cmd: %x\n",cmd[0]);
	ret=-1;
    }
    cmd[0]=-1;
    cmd[1]=-1;
 end:
    if(ret){
	RUN_T *irun=running_get_by_sock(sock);//is maos
	if(irun && irun->status.info<10){
	    sleep(1);
	    int pid2=irun->pid;
	    if(kill(pid2,0)){
		running_update(pid2,S_CRASH);
	    }
	}
	ret=-1;
    }
    return ret;/*don't close the port yet. may be reused by the client. */
}

static void scheduler_timeout(void){
    static int lasttime3=0;
    if(!all_done){
	process_queue();
    }
    /*Report CPU usage every 3 seconds. */
    int thistime=myclocki();
    if(thistime>=(lasttime3+3)){
	if(nrun>0){
	    check_jobs();
	}
	usage_cpu=get_usage_cpu();
	monitor_send_load();
	lasttime3=thistime;
    }
}

/*The following routines maintains the MONITOR_T linked list. */
static MONITOR_T* monitor_add(int sock){
    /*info2("added monitor on sock %d\n",sock); */
    MONITOR_T *node=calloc(1, sizeof(MONITOR_T));
    node->sock=sock;
    node->next=pmonitor;
    pmonitor=node;
    monitor_send_initial(node);
    return pmonitor;
}
static void monitor_remove(int sock){
    MONITOR_T *ic,*ic2=NULL;
    /*int freed=0; */
    for(ic=pmonitor; ic; ic2=ic,ic=ic->next){
	if(ic->sock==sock){
	    if(ic2){
		ic2->next=ic->next;
	    }else{
		pmonitor=ic->next;
	    }
	    //close(ic->sock); close panics accept
	    warning3("Remove monitor at %d\n", ic->sock);
	    shutdown(ic->sock, SHUT_WR);
	    free(ic);
	    break;
	}
    }
}
/*static MONITOR_T *monitor_get(int sock){
    MONITOR_T *ic;
    for(ic=pmonitor; ic; ic=ic->next){
	if(ic->sock==sock)
	    break;
    }
    return ic;
    }*/
static int monitor_send_do(RUN_T *irun, char *path, int sock){
    int cmd[3];
    cmd[1]=hid;
    cmd[2]=irun->pid;
    if(path){/*don't do both. */
	cmd[0]=CMD_PATH;
	return (stwrite(sock, cmd, 3*sizeof(int)) || stwritestr(sock, path));
    }else{
	cmd[0]=CMD_STATUS;
	return (stwrite(sock, cmd, 3*sizeof(int)) || stwrite(sock, &irun->status, sizeof(STATUS_T)));
    }
}
/* Notify alreadyed connected monitors job update. */
static void monitor_send(RUN_T *irun,char*path){
    MONITOR_T *ic;
 redo:
    for(ic=pmonitor; ic; ic=ic->next){
	int sock=ic->sock;
	if(monitor_send_do(irun,path,sock)){
	    monitor_remove(sock);
	    goto redo;
	}
    }
}
/* Notify alreadyed connected monitors machine load. */
static void monitor_send_load(void){
    MONITOR_T *ic;
    double mem=get_usage_mem();
    int cmd[3];
    cmd[0]=CMD_LOAD;
    cmd[1]=hid;
    int memi=(int)(mem*100);
    int cpui=(int)(usage_cpu*100);
    cmd[2]=memi | (cpui << 16);
 redo:
    for(ic=pmonitor; ic; ic=ic->next){
	int sock=ic->sock;
	if(!ic->load)
	    continue;

	if(stwrite(sock,cmd,sizeof(int)*3)){
	    monitor_remove(sock);
	    goto redo;//restart
	}
    }
}
/* Notify the new added monitor all job information. */
static void monitor_send_initial(MONITOR_T *ic){
    int sock;
    RUN_T *irun;
    sock=ic->sock;
    {
	int cmd[3];
	cmd[0]=CMD_VERSION;
	cmd[1]=scheduler_version;
	cmd[2]=hid;/*Fixme: sending hid to monitor is not good is hosts does not match. */
	if(stwrite(sock,cmd,sizeof(int)*3)){
	    return;
	}
    }
    for(irun=runned; irun; irun=irun->next){
	if(monitor_send_do(irun,irun->path,sock)|| monitor_send_do(irun,NULL,sock)){
	    monitor_remove(sock);
	    return;
	}   
    }
    for(irun=running; irun; irun=irun->next){
	if(monitor_send_do(irun,irun->path,sock)|| monitor_send_do(irun,NULL,sock)){
	    monitor_remove(sock);
	    return;
	}
    }
}

int main(){
    char fn[PATH_MAX];
    snprintf(fn, PATH_MAX, "%s/scheduler", TEMP);
    remove(fn);
    extern int is_scheduler;
    is_scheduler=1;
    listen_port(PORT, fn,respond, 1, scheduler_timeout, 0);
    remove(fn);
    exit(0);
}
