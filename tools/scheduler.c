/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifdef HAVE_CONFIG_H
#include "config.h" 
#endif
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include <sys/socket.h>
#include <poll.h>
#include "../sys/sys.h"
#include "scheduler.h"

//static char* scheduler_fnlog=NULL;
static int NGPU=0;
/**
   Struct to hold information of running jobs.
*/
typedef struct RUN_T{
	struct RUN_T* next;
	char* exe; /*Path to the executable.*/
	char* path;/*Job path and Job arguments.*/
	char* path0;/*same as path, with fields separated by \n instead of space*/
	status_t status;
	time_t last_time; //time of last report
	int pid;
	int pidnew;//the new pid
	int sock;//socket for primary maos connection
	int sock_cmd;//socket for secondary maos connection for maos_command()
	int nthread;//number of threads
	int ngpu;//number of gpus requested.
}RUN_T;

/**
  Struct to hold available monitors waiting for information.
*/
typedef struct MONITOR_T{
	int sock;
	int load;/*handle machine load information. */
	int plot;/*whether plot is enabled*/
	int http;/*client wants text mode */
	int (*func)(char*buf, int nlen, int mode, void *userdata);//function to send text
	void* userdata;//parameter for func
	struct MONITOR_T* next;
}MONITOR_T;


static MONITOR_T* pmonitor=NULL;//head of linked list
static int all_done=0;
/*A linked list to store running process*/
static RUN_T* running=NULL;/*points to the begining */
static RUN_T** running_end=&running;/*points to the last to add. */
/*A linked list to store ended process*/
static RUN_T* runned=NULL;
static RUN_T** runned_end=&runned;//points to the pointer of end slot
static double usage_cpu;
static RUN_T* running_add(int pid, int sock);
static RUN_T* running_get_by_pid(int pid);
static RUN_T* running_get_by_status(int status);
static RUN_T *running_get_by_sock(int sock);

//static RUN_T* runned_get(int pid);
static int runned_add(RUN_T* irun);
static void running_remove(int pid, int status);
//static MONITOR_T *monitor_get(int hostid);
static void monitor_remove(int sock);
static void scheduler_timeout(void);
static void monitor_send(RUN_T* run, const char* path);
static void monitor_send_initial(MONITOR_T* ic);
static void monitor_send_load(void);
static long counter=-1;//an negative index for a pending run
static int nused_cpu=0;//number of CPUs being used
static int nused_gpu=0;//number of GPUs being used
static time_t lasttime3=0;//every 3 seconds
static time_t lasttime10=0;//every 10 seconds
#define MAXSOCK 1024
static const char* sockname[MAXSOCK]={0};//saving the client name of each socket
static int nrun_handle(int cmd, int pid, int nthread, int ngpu){
	switch(cmd){
	case 0:
		if(pid==1){
			return nused_gpu;
		} else{
			return nused_cpu;
		}
		break;
	case 1: //put back resource
		nused_cpu+=nthread;
		nused_gpu+=ngpu;
		dbg_time("post by %d: nused_cpu %d->%d. nused_gpu %d->%d.\n",
			pid, nused_cpu-nthread, nused_cpu, nused_gpu-ngpu, nused_gpu);
		break;
	case 2: //remove resource
		nused_cpu-=nthread;
		nused_gpu-=ngpu;
		dbg_time("remove by %d: nused_cpu %d->%d. nused_gpu %d->%d.\n",
			pid, nused_cpu+nthread, nused_cpu, nused_gpu+ngpu, nused_gpu);
		if(nused_cpu<0){
			dbg_time("nused_cpu=%d\n", nused_cpu);
			nused_cpu=0;
		}
		if(nused_gpu<0){
			dbg_time("nused_gpu=%d\n", nused_gpu);
			nused_gpu=0;
		}
		break;
	case 3: //reset
		if(nused_cpu||nused_gpu){
			nused_cpu=0;
			nused_gpu=0;
			dbg_time("reset: nused_cpu=%d, nused_gpu=%d\n", nused_cpu, nused_gpu);
		}
		break;
	}
	return nused_cpu;
}
static int nrun_get(int type){
	return nrun_handle(0, type, 0, 0);
}
static void nrun_add(int pid, int nthread, int ngpu){
	nrun_handle(1, pid, nthread, ngpu);
}
static void nrun_sub(int pid, int nthread, int ngpu){
	nrun_handle(2, pid, nthread, ngpu);
}
/**
   The following runned_* routines maintains the linked list
   "runned" which contains the already finished jobs.
*/
static int runned_add(RUN_T* irun){
	/*move irun to end of runned list from running list. */
	irun->next=NULL;
	irun->status.done=1;
	*runned_end=irun; runned_end=&irun->next;
	if(irun->status.info<10){
		warning_time("Should be a finished process\n");
		return 1;
	} else{
		return 0;
	}
}
void runned_remove(int pid){
	RUN_T *irun=NULL;
	for(RUN_T **curr=&runned; *curr;){
		irun=*curr;
		if(irun->pid==pid || pid==-1){//remove node
			if(runned_end==&irun->next){
				runned_end=curr;
			}
			*curr=irun->next;
			irun->status.info=S_REMOVE;
			monitor_send(irun, NULL);
			free(irun->path);
			free(irun->path0);
			free(irun->exe);
			free(irun);
			if(pid!=-1) break;
		}else{
			curr=&irun->next;
		}
	}
	if(!irun && pid!=-1){
		warning_time("runned_remove: Record %s:%d not found!\n", HOST, pid);
	}
}
/*static RUN_T* runned_get(int pid){
	RUN_T* irun;
	for(irun=runned; irun; irun=irun->next){
		if(irun->pid==pid){
			break;
		}
	}
	return irun;
}*/
/**
   Restart a crashed/finished job
*/
static void runned_restart(int pid){
	RUN_T *irun=NULL;
	for(RUN_T **curr=&runned; *curr;){
		irun=*curr;
		if(irun->pid==pid){
			if(irun->status.info<10){
				dbg_time("status.info=%d, cancelled.\n", irun->status.info);
				break;
			}
			if(runned_end==&irun->next){
				runned_end=curr;
			}
			*curr=irun->next;
			//Insert to the beginning of running list
			irun->next=running;
			if(running_end==&running) running_end=&irun->next;
			running=irun; 
			//update status
			memset(&irun->status, 0, sizeof(status_t));
			irun->status.info=S_QUEUED;
			irun->status.timlast=myclocki();
			//pidnew changed to be different from pid to indicate restarted job.
			irun->pidnew=--counter;
			monitor_send(irun, NULL);
			irun->pid=irun->pidnew;//sync after notify monitor.
			irun->sock=0;
			all_done=0;
			break;
		}else{
			curr=&irun->next;
		}
	}
	if(!irun){
		warning_time("Process with pid=%d is not found\n", pid);
	}
}
/**
   The following running_* routiens operates on running linked list
   which contains queued jobs.
*/
static RUN_T* running_add(int pid, int sock){
	RUN_T* irun;
	if((irun=running_get_by_pid(pid))){
		//dbg_time("PID %d is already in running.\n", pid); /*create the node */
		if(irun->sock!=sock) irun->sock=sock;
		return irun;
	} else{
		char progname[PATH_MAX];
		if(pid>0){
			if(get_job_progname(progname, PATH_MAX, pid)){//failed. Already exited.
				warning_time("PID %d already exited. will not add\n", pid);
				return NULL;
			}
		}
		dbg2_time("adding %d to running\n", pid); /*create the node */
		irun=mycalloc(1, RUN_T);
		irun->pidnew=irun->pid=pid;
		irun->sock=sock;
		irun->status.timlast=irun->status.timstart=myclocki();
		if(pid>0){
			irun->exe=strdup(progname);
		}
		/*record the launch time */
		irun->next=NULL;
		*running_end=irun; running_end=&irun->next;//add to the end
		return irun;
	}
}

static void running_remove(int pid, int status){
	//dbg_time("Removing %d from running\n", pid);
	RUN_T* irun=NULL;
	for(RUN_T** curr=&running; *curr; ){
		irun=*curr;
		if(irun->pid==pid){
			if(running_end==&irun->next){
				running_end=curr;
			}
			*curr=irun->next;
			//Only the following states are used
			//S_QUEUED, S_WAIT, S_START, S_RUNNING, S_UNKNOWN, S_FINISH
			//nrun_add is called when state changes to S_START
			//S_START may change to S_RUNNING, S_UNKNOWN before this function is called
			//Only maos can change to S_RUNNING
			if(irun->pid>0&&(irun->status.info==S_START || irun->status.info==S_RUNNING || irun->status.info==S_UNKNOWN)){
				nrun_sub(irun->pid, irun->nthread, irun->ngpu);
			}
			irun->status.timlast=myclocki();
			irun->status.info=status;
			//dbg_time("remove job %d with status %d\n", pid, irun->status.info);
			monitor_send(irun, NULL);
			runned_add(irun);
			break;
		}else{
			curr=&irun->next;
		}
	}
	if(!irun){
		dbg_time("%s:%d not found\n", HOST, pid);
	}
}

static RUN_T* running_get_by_pid(int pid){
	RUN_T* irun;
	for(irun=running; irun; irun=irun->next){
		if(irun->pid==pid||(irun->pid>0&&pid==0)){
			break;
		}
	}
	return irun;
}

static RUN_T* running_get_by_status(int status){
	RUN_T* irun;
	//int jrun=0;
	for(irun=running; irun; irun=irun->next){
		//jrun++;
		if(irun->status.info==status){
			break;
		}
	}
	return irun;
}

static RUN_T* running_get_by_sock(int sock){
	RUN_T* irun;
	for(irun=running; irun; irun=irun->next){
		if(irun->sock==sock){
			break;
		}
	}
	return irun;
}
void running_kill(int pid){
	RUN_T* irun=running_get_by_pid(pid);
	if(irun){
		if(irun->status.info!=S_QUEUED){
			kill(pid, SIGTERM);
			if(irun->status.info==S_WAIT){//wait up the process.
				stwriteint(irun->sock, S_START);
			}
			irun->status.info=S_UNKNOWN;
			monitor_send(irun, NULL);
		} else{
			running_remove(pid, S_KILLED);
		}
		//dbg_time("%5d term signal sent\n", pid);
	}
}
/**
	check all the jobs. remove if any job quited.
 */
static void check_jobs(void){
	RUN_T* irun=NULL, *irun_next=NULL;
	int nrunning=0;
	if(running){
		time_t now=myclocki();
		for(irun=running; irun; irun=irun_next){
			irun_next=irun->next;
			if(irun->pid>0){//Running job
				if(kill(irun->pid, 0)){//No longer exists
					if(irun->last_time+60<now){//allow grace period.
						dbg_time("check_jobs: Job %d no longer exists. Change status to S_CRASH\n", irun->pid);
						running_remove(irun->pid, S_CRASH);
					}else{
						nrunning++;
					}
				} else {
					nrunning++;
					if((irun->last_time+600<now)&&irun->status.info==S_RUNNING){
						dbg_time("check_jobs: Job %d does not update after %lu seconds. Change status to S_UNKNOWN\n",
								irun->pid, now-irun->last_time);
						irun->status.info=S_UNKNOWN;
						monitor_send(irun, NULL);
					}
				}
			}
		}
	}
	if(!nrunning && (nrun_get(0) || nrun_get(1))){
		warning("There are no running jobs, but ncpu=%d, ngpu=%d. reset to 0.\n", nrun_get(0), nrun_get(1));
		nused_cpu=0;
		nused_gpu=0;
	}

}
/**
 examine the process queue to start routines once CPU is available.
*/
static void process_queue(void){
	//dbg_time("process_queue: enter\n");
	if(!running) return;
	static time_t lasttime=0;
	time_t thistime=myclocki();
	if(nrun_get(0)>0&&thistime<lasttime+1){
		//dbg_time("process_queue: too son\n");
		return;
	}
	lasttime=thistime;
	if(nrun_get(0)>=NCPU||(NGPU&&nrun_get(1)>=NGPU)){
		//dbg_time("process_queue: enough jobs are running\n");
		return;
	}
	//dbg2_time("nused_cpu=%d, ngpu=%d, avail=%d\n", nrun_get(0), nrun_get(1), avail);
	RUN_T* irun=running_get_by_status(S_WAIT);
	if(irun){//There are jobs waiting.
		if(irun->sock>0){//already connected.
			int nthread=irun->nthread;
			int avail=get_cpu_avail();
			if(nrun_get(0)+nthread<=NCPU&&(nthread<=avail||avail>=3)
				&&(!NGPU||!irun->ngpu||nrun_get(1)+irun->ngpu<=NGPU)){//resource available to star the job
				irun->last_time=myclocki();
				irun->status.timlast=myclocki();
				if(stwriteint(irun->sock, S_START)){
					dbg_time("Starting Job %d at %d failed: %s\n", irun->pid, irun->sock, strerror(errno));
				} else{
					dbg_time("Starting Job %d at %d ok.\n", irun->pid, irun->sock);
				}
				//we mark the sate as running in either state and let check_job do the clean up
				irun->status.info=S_START;
				nrun_add(irun->pid, irun->nthread, irun->ngpu);
				monitor_send(irun, NULL);
			}
		} else{
			if(kill(irun->pid, 0)){
				dbg_time("Job %d already exited. irun->sock=%d, change to UNKNOWN.\n", irun->pid, irun->sock);
				irun->status.info=S_UNKNOWN;
			} else{
				dbg_time("Wait for %d to connect. irun->sock=%d\n", irun->pid, irun->sock);
			}
		}
	} else{
		if(nrun_get(0)<NTHREAD&&(!NGPU||nrun_get(1)<NGPU)){//resource available to star a new job
			//static double lasttime2=0;
			//time_t thistime2=myclocki();
			//if(thistime2>lasttime2+0.001){//wait 1ms for job to connect.
				//lasttime2=thistime2;
				irun=running_get_by_status(S_QUEUED);
				if(!irun){
					dbg_time("all jobs are done\n");
					all_done=1;
					nrun_handle(3, 0, 0, 0);
					counter=-1; //reset the counter
				} else{
					int pid;
					irun->last_time=myclocki();
					irun->status.timlast=myclocki();
					if((pid=launch_exe(irun->exe, irun->path0))<0){
						dbg_time("Launch job %d failed: %d (%s)\n", irun->pid, pid, irun->exe);
						running_remove(irun->pid, S_CRASH);
					} else{
						dbg_time("Launch job %d as pid %d (%s)\n", irun->pid, pid, irun->exe);
						//inplace update the information in monitor
						irun->status.info=S_WAIT;//will be checked again by process_queue
						irun->pidnew=pid;
						monitor_send(irun, NULL);
						irun->pid=pid;
					}
				}
			//}
		}else{
			dbg_time("no resource available. nused_cpu=%d, ngpu=%d.\n", nrun_get(0), nrun_get(1));
		}
	}
}
/** replace \n by space*/
static char* remove_endl(const char* path0){
	char* path=strdup(path0);
	for(char* tmp=path; tmp[0]; tmp++){
		if(tmp[0]=='\n'){
			tmp[0]=' ';
		}
	}
	return path;
}
static void queue_new_job(const char* exename, const char* execmd){
	RUN_T* irun=running_add(--counter, -1);
	if(!irun){
		warning_time("scheduler: running_add failed.\n");
		return;
	}
	irun->status.info=S_QUEUED;
	irun->exe=strdup(exename);
	irun->path0=strdup(execmd);
	irun->path=remove_endl(irun->path0);
	dbg_time("%d (%s)\n", irun->pid, exename);
	monitor_send(irun, irun->path);
	monitor_send(irun, NULL);
	all_done=0;
}

/**
 * @brief Sends a command to a running process (maos) and optionally passes a socket.
 *
 * This function retrieves the running process by its PID, then attempts to send a command
 * (and optionally a socket) to the process via its command socket. It also replies to the
 * provided socket (if valid) with the result of the operation.
 *
 * @param pid   The process ID of the target running process.
 * @param sock  The socket descriptor to reply to, or -1 if no reply is needed.
 * @param cmd   The command to send to the running process.
 *
 * @return 0 on success, -1 on failure.
 *
 * The function logs a warning if the operation fails, or a debug message on success.
 * If the process or its command socket is not found, or if writing the command fails,
 * the function returns -1.
 */
int maos_command(int pid, int sock, int cmd){
	RUN_T* irun=running_get_by_pid(pid);
	int cmd2[2]={cmd, 0};
	int failed=(!irun||!irun->sock_cmd||(stwriteintarr(irun->sock_cmd, cmd2, 2)));
	if(sock!=-1){
		if(!stwriteint(sock, failed?-1:0) && !failed){//scheduler reply first before maos gets the chance to start writing to the socket
			stwritefd(irun->sock_cmd, sock);
		}
	}
	if(failed){
		warning_time("Failed to send command %d and socket %d to maos (PID=%d, sock_cmd=%d)\n", cmd, sock, irun?irun->pid:-1, irun?irun->sock_cmd:-1);
	} else{
		dbg_time("Successfully send command %d and socket %d to maos (PID=%d, sock_cmd=%d)\n", cmd, sock, irun->pid, irun->sock_cmd);
	}
	return failed?-1:0;//do not keep this connection.
}
/*
 Stores drawdaemon sockets
 */
typedef struct SOCKID_M{
	int id;
	int sock;
	struct SOCKID_M* prev;
	struct SOCKID_M* next;
} SOCKID_T;
static SOCKID_T* shead=0;

/**
 * @brief Saves or updates the association between a socket and an ID.
 *
 * This function searches for an existing entry in the linked list of SOCKID_T structures
 * that matches the given socket descriptor (`sock_save`). If found, it updates the associated
 * ID with the provided `id` value and logs the update. If not found, it creates a new entry,
 * initializes it with the given socket and ID, inserts it at the head of the list, and logs the save.
 *
 * @param sock_save The socket descriptor to be saved or updated.
 * @param id        The ID to associate with the socket.
 */
static void socket_save(int sock_save, int id){
	int found=0;
	for(SOCKID_T* p=shead; p; p=p->next){
		if(p->sock==sock_save){
			found=1;
			p->id=id;
			dbg_time("Update socket %d with id %d\n", sock_save, id);
		}
	}
	if(!found){
		SOCKID_T* tmp=mycalloc(1, SOCKID_T);
		tmp->id=id;
		tmp->sock=sock_save;
		tmp->next=shead;
		tmp->prev=0;
		if(shead){
			shead->prev=tmp;
		}
		shead=tmp;
		dbg_time("Save socket %d with id %d\n", sock_save, id);
	}
}
static void socket_remove(SOCKID_T*p){
	//remove from list. make sure p->next is saved
	if(p->prev){//middle item
		p->prev->next=p->next;
	} else{//first item
		shead=p->next;
	}
	if(p->next){//not last item
		p->next->prev=p->prev;
	}
	free(p);
}
//retrieve socket from linked list based in id. 
static int socket_get(int id){
	int sock_save=-1;
	//dbg_time("Looking up socket with id %d\n", id);
	for(SOCKID_T* p_next, *p=shead; p; p=p_next){
		p_next=p->next;//save so that socket_remove can be called
		int badsock=0;
		if(p->sock==-1){//should not happen
			socket_remove(p);
		}else if((badsock=stcheck(p->sock))||p->id==id){
			if(badsock){
				dbg_time("Remove bad socket %d with id %d\n", p->sock, p->id);
				close(p->sock); //closed socket
				p->sock=-1;
				socket_remove(p);
			} else if(sock_save==-1){
				sock_save=p->sock; //valid match
				dbg_time("Get socket %d with id %d\n", p->sock, p->id);
				p->sock=-1;
				socket_remove(p);
			}else{
				dbg_time("Skip socket %d with id %d\n", p->sock, p->id);
				//check next socket for disconnection.
			}
		}else{
			dbg_time("Skip socket %d with id %d\n", p->sock, p->id);
		}
	}
	if(sock_save==-1&&id!=0){
		sock_save=socket_get(0);//if not found, also check id=0
	}
	return sock_save;	
}
static void socket_heartbeat(){
	for(SOCKID_T* p_next, *p=shead; p; p=p_next){
		p_next=p->next;
		int cmd[1]={DRAW_HEARTBEAT};
		if(stcheck(p->sock) || stwrite(p->sock, cmd, sizeof(cmd))){
			dbg_time("Remove bad socket %d with id %d\n", p->sock, p->id);
			socket_remove(p);
		}
	}
}
/*Close all connection upon scheduler exit.*/
static void socket_close(){
	for(SOCKID_T *p_next, *p=shead; p; p=p_next){
		p_next=p->next;
		close(p->sock);
		socket_remove(p);
	}
}
static void set_sockname(int sock, const char *name){
	if(sock>-1&&sock<MAXSOCK&&!sockname[sock]){
		sockname[sock]=name;
	}
}
static void unset_sockname(int sock){
	if(sock>-1&&sock<MAXSOCK){
		sockname[sock]=NULL;
	}
}
static const char* get_sockname(int sock){
	if(sock>-1 && sock<MAXSOCK && sockname[sock]){
		return sockname[sock];
	}else{
		return "";
	}
}

static int scheduler_recv_wait=-1;//>-1: there is pending scheduler_recv_socket.
static int monitor_request_draw(int sock, int pid){
	MONITOR_T* pm=0;
	int ans=-1;
	for(pm=pmonitor; pm; pm=pm->next){
		if(pm->plot){
			if(pm->http){
				char buf[4096];
				int padding=pm->func?http_padding:(int)(3*sizeof(int));
				int nlen=snprintf(buf+padding, sizeof(buf)-padding, "%d&DRAW$", pid);
				if(pm->func){//directly write to client
					ans=pm->func(buf+padding, nlen, 0, pm->userdata);
				}else{//send to proxy
					((int*)buf)[0]=DRAW_ENTRY;
					((int*)buf)[1]=nlen;
					((int*)buf)[2]=0;//0 indicate text data
					ans=stwrite(pm->sock, buf, nlen+padding);
				}
			}else{
				int moncmd[3]={MON_DRAWDAEMON,0,pid};
				ans=stwrite(pm->sock, moncmd, sizeof(moncmd));
			}
			scheduler_recv_wait=sock;
			dbg_time("(%d:%s) request monitor at %d to start drawdaemon\n", sock, get_sockname(sock), pm->sock);
			break;
		}
	}
	return ans;
}
/**
 * @brief Send socket to draw for drawing.
 * 
 * @param sock 
 * @param pid 
 * @return int 0 for success. -1 for failure
 */
int send_draw_sock(int sock, int pid){
	int ret;
	if(pid<=0){//this is for pending scheduler_recv_socket or a new idle connection
		if(scheduler_recv_wait==-1||stwriteint(scheduler_recv_wait, 0)||stwritefd(scheduler_recv_wait, sock)){
			int sock2=dup(sock);
			if(scheduler_recv_wait!=-1){
				dbg_time("(%d:%s) Failed to pass sock to draw at %d, save socket as %d for future\n", sock, get_sockname(sock), scheduler_recv_wait, sock2);
			}else{
				dbg_time("(%d:%s) No pending drawdaemon request, save socket as %d for future\n", sock, get_sockname(sock), sock2);
			}
			socket_save(sock2, abs(pid));//duplicate socket and keep it 
			set_sockname(sock2, "drawdaemon");
			ret=0;//prevent scheduler from listening to this socket.
		} else{
			dbg_time("passed sock %d to draw at %d\n", sock, scheduler_recv_wait);
			ret=0;//close socket on scheduler.
		}
		stwriteint(sock, 0);
		scheduler_recv_wait=-1;
	} else{//ask maos to start drawing with this drawdaemon
		ret=maos_command(pid, sock, MAOS_DRAW);
	}
	dbg_time("returns %d for sock=%d, pid=%d\n", ret, sock, pid);
	return ret;
}
//PNEW(mutex_sch);//respond() and scheduler_handle_ws() much lock this before preceed.
/**
   respond to client requests. The fixed header is int[2] for opcode and
   pid. Optional data follows based on opcode. The disadvanced of this is that
   new commands with additional payload cannot be easily introduced in the
   scheduler().
*/
static int respond(struct pollfd *pfd, int flag){
	if(!pfd){
		return 0;
	}
	if(flag==-1){//ask client to shutdown
		int cmd[3]={-1, 0, 0};
		stwrite(pfd->fd, cmd, 3*sizeof(int));
	}
	int sock=pfd->fd;
	/*
	   Don't modify sock in this routine. otherwise select
	   will complain Bad file descriptor
	*/
	int ret=0, pid, cmd[2]={0,0};
	//dbg_time("\rrespond %2d start ... ", sock);errno=0;
	if((ret=streadintarr(sock, cmd, 2))){
		/*if(errno!=ESRCH && errno!=ENOENT){//timeout or closed
			dbg_time("read %d failed (%d): %s,ret=%d\n", sock, errno, strerror(errno), ret);
		}*/
		//dbg_time("(%d:%s) connection lost or closed by client.\n", sock, get_sockname(sock)); 
		unset_sockname(sock);
		return ret;
	}else{
		dbg2_time("(%d:%s) cmd is %d %d. \n", sock, get_sockname(sock), cmd[0], cmd[1]);
	}
	pid=cmd[1];
	//LOCK(mutex_sch);
	switch(cmd[0]){
	case CMD_START://1: Called by maos when job starts.
	{
		set_sockname(sock, "maos");
		/* Data read. */
		int nthread;
		int ngpu;
		int waiting;
		if(streadint(sock, &nthread)||streadint(sock, &waiting)){
			ret=-1;
			break;
		}
		ngpu=waiting>>1;
		if(ngpu>NGPU){
			ngpu=NGPU;
		}
		if(ngpu)//maos restrict number of threads to ngpu+1 when gpu is in use.
			nthread=ngpu;
		else if(nthread<1)
			nthread=1;
		else if(nthread>NCPU)
			nthread=NCPU;

		RUN_T* irun=running_add(pid, sock);
		if(!irun){
			warning_time("scheduler: running_add %d failed. Exe already exited.\n", pid);
			break;
		}
		irun->nthread=nthread;
		irun->ngpu=ngpu;
		if(waiting&0x1){
			dbg_time("(%d:%s) Job %d waiting to start with nused_cpu=%d, ngpu=%d\n", sock, get_sockname(sock), pid, nthread, ngpu);
			irun->status.info=S_WAIT;
		} else{/*no waiting, no need reply. */
			dbg_time("(%d:%s) Job %d started with nused_cpu=%d, ngpu=%d\n", sock, get_sockname(sock), pid, nthread, ngpu);
			irun->status.info=S_START;
			nrun_add(pid, nthread, ngpu);
		}
		all_done=0;
		if(irun->path) monitor_send(irun, irun->path);
		monitor_send(irun, NULL);
	}
	break;
	case CMD_FINISH://2: Called by MAOS when job finishes.
		dbg_time("(%d:%s) Job %d reports finished\n", sock, get_sockname(sock), pid);
		running_remove(pid, S_FINISH);
		ret=0;//-1. Do not yet close. wait for client to close.
		break;
	case CMD_STATUS://3: Called by MAOS to report status at every time step
	{
		RUN_T* irun=running_get_by_pid(pid);
		if(!irun){/*started before scheduler is relaunched. */
			dbg_time("(%d:%s) pid=%d is running but not recorded.\n", sock, get_sockname(sock), pid);
			irun=running_add(pid, sock);
			if(!irun){
				warning_time("scheduler: running_add %d failed. Exe already exited.\n", pid);
				break;
			}
			irun->status.info=S_START;
			irun->nthread=irun->status.nthread;
			nrun_add(pid, irun->nthread, irun->ngpu);
			monitor_send(irun, NULL);
		}
		//save and restore time information
		time_t timlast=irun->status.timlast;
		time_t timstart=irun->status.timstart;
		if(sizeof(status_t)!=read(sock, &(irun->status), sizeof(status_t))){
			warning_time("Error reading status.\n");
		}
		irun->status.timlast=timlast;
		irun->status.timstart=timstart;
		irun->last_time=myclocki();
		monitor_send(irun, NULL);
	}
	break;
	case CMD_CRASH://4: called by MAOS when job crashes
		dbg_time("(%d:%s) Job %d reports crashed\n", sock, get_sockname(sock), pid);
		running_remove(pid, S_CRASH);
		ret=0;//-1; wait for client to close.
		break;
	case CMD_MONITOR://5: Called by Monitor when it connects
	{
		set_sockname(sock, "monitor");
		monitor_add(sock, pid, NULL, NULL);

		dbg_time("(%d:%s) monitor from %s is connected.\n", sock, get_sockname(sock), addr2name(socket_peer(sock)));
	}
	break;
	case CMD_PATH://6: Called by MAOS to report the PATH.
	{
		set_sockname(sock, "maos");
		RUN_T* irun=running_add(pid, sock);
		if(!irun){
			warning_time("(%d:%s) running_add %d failed. Exe already exited.\n", sock, get_sockname(sock), pid);
			break;
		}
		if(irun->path0) free(irun->path0);
		if(irun->path) free(irun->path);
		if(streadstr(sock, &irun->path0)){
			dbg_time("(%d:%s) Job %d receiving path failed. \n", sock, get_sockname(sock), pid);
			ret=-1;
			break;
		}
		irun->path=remove_endl(irun->path0);
		dbg_time("(%d:%s) Job %d received path. \n", sock, get_sockname(sock), pid);
		monitor_send(irun, irun->path);
	}
	break;
	case CMD_KILL://7: Called by Monitor to kill a task.
	{
		warning_time("(%d:%s) Received monitor command to kill %d\n", sock, get_sockname(sock), pid);
		running_kill(pid);
	}
	break;
	case CMD_TRACE://8: Called by MAOS to request a backtrace
	{
		set_sockname(sock, "trace");
		char* buf=NULL, out[200];
		if(streadstr(sock, &buf)
			||call_addr2line(out, sizeof out, buf)
			||stwritestr(sock, out)){
			warning_time("(%d:%s) CMD_TRACE failed. buf=%s, out=%s\n", sock, get_sockname(sock), buf, out);
			ret=-1;
		}else{
			//info_time("(%d:%s) CMD_TRACE success. buf=%s, out=%s\n", sock, get_sockname(sock), buf, out);
		}
		free(buf);
	}
	break;
	case CMD_PROBE://9: called by monitor to probe the connection
	{	//send back a dummy response
		int cmd2[3]={MON_VERSION, 0, 0};
		stwrite(sock, &cmd2, sizeof(cmd2));
	}
	break;
	case CMD_DRAWCLI://10:Called by draw() to cache or request a fd. Valid over UNIX socket only.
	{
		if(pid>0){//receive sock from draw() and save for later use
			set_sockname(sock, "save drawdaemon");
			int sock_save;
			if(streadfd(sock, &sock_save)){
				warning_time("(%d:%s) receive socket for saving failed\n", sock, get_sockname(sock));
				sock_save=-1;
				ret=-1;
			} else{
				dbg_time("(%d:%s) received socket %d for saving.\n", sock, get_sockname(sock), sock_save);
				socket_save(sock_save, pid);
			}
		} else if(pid<0 && pid>-10000){//send existing sock to draw()
			set_sockname(sock, "get drawdaemon");
			int sock_save=socket_get(-pid);//drawdaemon with the same session id.
			dbg_time("(%d:%s) received socket request, sock_saved=%d\n", sock, get_sockname(sock), sock_save);

			//cannot pass -1 as sock, so return a flag first. sock can be zero.
			if(stwriteint(sock, sock_save>-1?0:-1)){
				warning_time("(%d:%s) send socket status failed.\n", sock, get_sockname(sock));
				ret=-1;
			}
			if(sock_save>-1){
				if(stwritefd(sock, sock_save)){
					warning_time("(%d:%s) send socket %d failed\n", sock, get_sockname(sock), sock_save);
					ret=-1;//close connection to draw()
				} else{//socket is transferred to draw. we close it.
					dbg_time("(%d:%s) send socket %d from %s success\n", sock, get_sockname(sock), sock_save, addr2name(socket_peer(sock_save)));
				}
				close(sock_save);
			}
		} else{//pid==0 or pid<10000; request a drawdaemon using monitor. 
			set_sockname(sock, "open drawdaemon");
			if(pid<-10000) pid+=10000;
			if(monitor_request_draw(sock, pid)){
				//there is no available drawdameon. Need to create one by sending request to monitor.
				warning_time("(%d:%s) there is no monitor available to start drawdaemon\n", sock, get_sockname(sock));
				if(stwriteint(sock, -1)){
					warning_time("(%d:%s) failed to respond to draw.\n", sock, get_sockname(sock));
				}
			}
		}
	}
	break;
	case CMD_REMOVE://11: Called by Monitor to remove a finished job fron the list*/
		runned_remove(pid);
	break;
	case CMD_DRAWSER://12: called by Drawdaemon/Monitor to provide drawdaemon to draw().*/
		set_sockname(sock, "drawdaemon");
		dbg_time("(%d:%s) drawdaemon for pid=%d from %s.\n", sock, get_sockname(sock), pid, addr2name(socket_peer(sock)));
		ret=send_draw_sock(sock, pid)?-2:-1;
		break;
	case CMD_MAOSCLI://13:  for a maos client to create a client link to maos.
		set_sockname(sock, "maos client");
		dbg_time("(%d:%s) pass command to maos %d for client at %s\n", sock, get_sockname(sock), pid, addr2name(socket_peer(sock)));
		ret=maos_command(pid, sock, MAOS_VAR)?-2:-1;
		break;
	case CMD_MAOSSER://14: called by maos to save a port to run maos_command
	{
		set_sockname(sock, "maos server");
		RUN_T* irun=running_get_by_pid(pid);
		if(irun){
			irun->sock_cmd=sock;
			dbg_time("(%d:%s) maos server socket is saved\n", sock, get_sockname(sock));
		} else{
			warning_time("(%d:%s) maos server socket does not have irun. close.\n", sock, get_sockname(sock));
			ret=-1;
		}
	}
	break;
	case CMD_RESTART://15: Called by monitor to restart a job
		dbg_time("(%d:%s) restart job %d\n", sock, get_sockname(sock), pid);
		runned_restart(pid);
		break;
	case CMD_KILLED://16: called by maos to indicate that job is cancelled or killed
		running_remove(pid, S_KILLED);
		ret=0;//-1; //wait for client to close.
		break;
	case CMD_DUMMY://17; probe. no response.
		break;
	case CMD_LAUNCH://18: called from maos from another machine to start a job in this machine
	{
		set_sockname(sock, "job launcher");
		char* exename=NULL;
		char* execwd=NULL;
		char* execmd=NULL;
		if(pid>=2){
			if(streadstr(sock, &exename)){
				warning_time("(%d:%s) Unable to read exename.\n", sock, get_sockname(sock));
				ret=-1;
			}
		}
		if(ret!=-1&&pid>=3){//old method. should not be used.
			if(streadstr(sock, &execwd)){
				warning_time("(%d:%s) Unable to read execwd.\n", sock, get_sockname(sock));
				ret=-1;
			}
		}
		if(ret!=-1&&!streadstr(sock, &execmd)){
			if(execwd){//merge cwd to new
				char* tmp=execmd;
				execmd=stradd(execwd, "/", execmd, NULL);
				free(tmp);
			}
			if(strstr(execmd, "--server")){//Launch immediately in server mode
				int pidn=-1;
				if((pidn=launch_exe(exename, execmd))<0){
					warning_time("launch_exe %s failed\n", exename);
					ret=-1;
				}
			} else{//Queue in jobs.
				queue_new_job(exename, execmd);
			}
		} else{
			warning_time("(%d:%s) Unable to read execmd\n", sock, get_sockname(sock));
			ret=-1;
		}
		if(stwriteint(sock, ret)){
			ret=-1;
		}
		free(exename);
		free(execwd);
		free(execmd);
	}
	break;
	default:
		warning_time("(%d:%s) Invalid cmd: %x\n", sock, get_sockname(sock), cmd[0]);
		ret=-1;
	}
	//job is finished. process the next job. don't wait for timeout
	if((cmd[0]==CMD_KILLED||cmd[0]==CMD_FINISH||cmd[0]==CMD_CRASH)&&running){
		process_queue();
	}
	if(ret){
		RUN_T* irun=running_get_by_sock(sock);//is maos
		if(irun&&irun->status.info<10){
			//connection failed to a running maos job.
			if(kill(irun->pid, 0)){
				dbg_time("(%d:%s) Job %d no longer exists, crashed?\n", sock, get_sockname(sock), irun->pid);
				running_remove(irun->pid, S_CRASH);
			}
		}
		ret=-1;
	}
	if(ret<0){
		dbg_time("cmd=%d returned ret=%d\n", cmd[0], ret);
	}
	//if(cmd[0]!=CMD_STATUS||ret) dbg_time("(%d:%s) respond returns %d for %d.\n", sock, get_sockname(sock), ret, cmd[0]);
	return ret;//ret=-1 will close the socket.
}


static void scheduler_timeout(void){
	time_t thistime=myclocki();
	//Process job queue
	if(running){
		process_queue();
	}
	//respond to heatbeat of saved sockets

	if(thistime>=(lasttime3+3)){
		if(running){
			check_jobs();
		}
		
		lasttime3=thistime;
	}
	//Report CPU usage and drawdaemon heartbeat every 10 seconds. 
	if(thistime>=(lasttime10+10)){
		usage_cpu=get_usage_cpu();
		monitor_send_load();
		socket_heartbeat();
		lasttime10=thistime;
	}
}

/*The following routines maintains the MONITOR_T linked list. */
void monitor_add(int sock, int flag, int (*func)(char* buf, int nlen, int mode, void *userdata), void* userdata){
	//dbg_time("added monitor on sock %d\n",sock);
	MONITOR_T *node=NULL;
	for(MONITOR_T *ic=pmonitor; ic; ic=ic->next){
		if(ic->sock==sock){
			node=ic;
			break;
		}
	}
	if(!node){
		node=mycalloc(1, MONITOR_T);
		node->sock=sock;
		node->next=pmonitor;
		pmonitor=node;
		dbg_time("monitor is added: fd=%d flag=%x\n", node->sock, flag);
	}else{
		dbg_time("monitor is updated: fd=%d flag=%x\n", node->sock, flag);
	}
	
	if(flag>=0x8){//old scheme where scheduler_version at compiling is passed
		node->load=1;
	} else if(flag<0){//new scheme with bit indicating different options
		node->load=flag&1;
		node->plot=(flag>>1)&1;
		node->http=(flag>>2)&1;
	}
	node->func=func;
	node->userdata=userdata;
	monitor_send_initial(node);
}
static void monitor_remove(int sock){
	MONITOR_T* ic=NULL;
	for(MONITOR_T **curr=&pmonitor; *curr;){
		ic=*curr;
		if(ic->sock==sock){
			*curr=ic->next;
			dbg_time("Removed monitor at %d\n", ic->sock);
			shutdown(ic->sock, SHUT_WR);
			//close(ic->sock);// close panics accept
			free(ic);
			return;
		}else{
			curr=&ic->next;
		}
	}
	warning_time("Monitor is not found fd=%d\n", sock);
}
static int http_convert(RUN_T *irun, const char *path, char *buf, int size){
	int len;
	if(path){
		len=snprintf(buf, size, "%d&PATH&%s$", irun->pid, path);
	} else{
		status_t* st=&irun->status;
		struct tm* tim=localtime(&(st->timlast));
		char stime[80];
		strftime(stime, 80, "[%a %H:%M:%S]", tim);
		len=snprintf(buf, size, "%d&STATUS&%d&%d" /*pid, key, pidnew, status*/
			"&%s&%.2f&%.2f" /*start time, errhi, errlo*/
			"&%d&%d&%d&%d" /*iseed, nseed, isim, nsim*/
			"&%ld&%ld&%.3f$" /*rest, tot, step timing*/
			, irun->pid, irun->pidnew, st->info,
			stime, st->clerrhi, st->clerrlo,
			st->nseed==0?0:st->iseed+1, st->nseed, st->simend==0?0:st->isim+1, st->simend,
			st->rest, st->laps+st->rest, st->tot*st->scale);
	}
	return len;
}
static int monitor_send_do(RUN_T* irun, const char* path, MONITOR_T *ic){
	if(!irun || !ic) return -1;
	int ans=0;
	if(ic->http){
		char buf[4096];
		int padding=ic->func?http_padding:(int)(3*sizeof(int));
		int nlen=http_convert(irun, path, buf+padding, sizeof(buf)-padding);
		if(ic->func){//directly write to client
			ans=ic->func(buf+padding, nlen, 0, ic->userdata);
		}else{//send to proxy
			((int*)buf)[0]=DRAW_ENTRY;
			((int*)buf)[1]=nlen;
			((int*)buf)[2]=0;//0 indicate text data
			ans=stwrite(ic->sock, buf, nlen+padding);
		}
	}else{//write to TCP client
		int cmd[3];
		cmd[1]=irun->pidnew;
		cmd[2]=irun->pid;
		if(path){/*don't do both. */
			cmd[0]=MON_PATH;
			ans=(stwrite(ic->sock, cmd, 3*sizeof(int))||stwritestr(ic->sock, path));
		} else{
			cmd[0]=MON_STATUS;
			ans=(stwrite(ic->sock, cmd, 3*sizeof(int))||stwrite(ic->sock, &irun->status, sizeof(status_t)));
		}
	}
	return ans;
}

/* Notify alreadyed connected monitors job update. */
static void monitor_send(RUN_T* irun, const char* path){
	MONITOR_T* ic;
redo:
	for(ic=pmonitor; ic; ic=ic->next){
		if(monitor_send_do(irun, path, ic)){
			dbg_time("monitor_send to fd=%d failed, remove\n", ic->sock);
			monitor_remove(ic->sock);
			goto redo;
		}
	}
}
double get_usage_gpu(double* const gpu_mem){
	static int enabled=1;
	*gpu_mem=0;
	if(enabled){
		FILE* fp=popen("nvidia-smi --format=csv,nounits,noheader --query-gpu=memory.used,memory.total", "r");
		if(!fp){
			dbg_time("popen nvidia-smi failed, disable gpu usage\n");
			enabled=0;
		} else{
			long used0=0, total0=0;
			long used, total;
			long nuse=0, ngpu=0;
			while(fscanf(fp, "%ld, %ld\n", &used, &total)==2){
				used0+=used;
				total0+=total;
				if(used>500){
					nuse++;
				}
				ngpu++;
			}
			pclose(fp);
			if(!total0){
				enabled=0;
				dbg_time("total0=0, disable gpu usage.\n");
			}
			*gpu_mem=(double)used0/total0;
			return (double)nuse/ngpu;
		}
	}
	return 0;
}
/* Notify alreadyed connected monitors machine load. */
static void monitor_send_load(void){
	MONITOR_T* ic, * ic2;
	double mem=get_usage_mem();
	double gpu_mem=0;
	double gpu=0;
	/*if(NGPU){
		gpu=get_usage_gpu(&gpu_mem);
	}*/
	int cmd[3];
	cmd[0]=MON_LOAD;
	cmd[1]=0;
	int memi=(int)((MAX(mem, gpu_mem))*100);
	int cpui=(int)((MAX(usage_cpu, gpu))*100);
	cmd[2]=(memi&0xFFFF)|(cpui<<16);

	for(ic=pmonitor; ic; ic=ic2){
		ic2=ic->next;
		int sock=ic->sock;
		if(!ic->load)
			continue;

		if(stwrite(sock, cmd, sizeof(int)*3)){
			monitor_remove(sock);
		}
	}
}
/* Notify the new added monitor all job information. */
static void monitor_send_initial(MONITOR_T* ic){
	int sock=ic->sock;
	for(RUN_T* irun=runned; irun; irun=irun->next){
		if(monitor_send_do(irun, irun->path, ic)||monitor_send_do(irun, NULL, ic)){
			monitor_remove(sock);
			return;
		}
	}
	for(RUN_T* irun=running; irun; irun=irun->next){
		if(monitor_send_do(irun, irun->path, ic)||monitor_send_do(irun, NULL, ic)){
			monitor_remove(sock);
			return;
		}
	}
}

int main(int argc, const char* argv[]){
	int flag=0;
	//err2out=0;//send error() and dbg() to stderr
	//std2out=0;//send info() to stderr
	LOG_LEVEL=1;
	if(argc>1 && argv[1]){
		if(!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")){
			info("%s [options]\n"
			"Options:\n"
			"\t-h --help: \tprint this help\n"
			"\t-d: \t\tfork twice and launch as daemon\n"
			"\t--direct: \tskip all checking\n", argv[0]);
			return 0;
		}
		if(!strcmp(argv[1], "-d")){//daemonize
			flag=1;
		} else if(!strcmp(argv[1], "--direct")){//launched by launch_scheduler_do, already locked.
			flag=2;
		}
	}
	//Handle different ways of launching drawdaemon
	switch(flag){
	case 0:{
		char fnlock[PATH_MAX];
		snprintf(fnlock, PATH_MAX, "%s/scheduler.pid", TEMP);
		int fd=0;
		if((fd=lock_file_version(fnlock, 0, scheduler_version))<0){
			dbg_time("Failed to lock_file. exit\n");
			return 1;
		}else{
			char strpid[60];
			snprintf(strpid, 60, "%d %d\n", getpid(), scheduler_version);
			lseek(fd, 0, SEEK_SET);
			if(ftruncate(fd, 0)<0)
				warning("Unable to truncate file\n");
			if(write(fd, strpid, strlen(strpid)+1)!=(long)strlen(strpid)+1){
				warning("Write pid %d to %s failed\n", getpid(), fnlock);
			}
		}
		//do not close fd.
	}
		break;
	case 1:{//daemonize
		single_instance_daemonize(TEMP, "scheduler", scheduler_version, NULL, NULL);
	}
		break;
	case 2://do nothing
		break;
	}
	char slocal[PATH_MAX];//local port
	snprintf(slocal, PATH_MAX, "%s/scheduler", TEMP);
	//remove(slocal);
	extern int is_scheduler;
	is_scheduler=1;
	setbuf(stdout, NULL);
	{
		//Find out number of gpus.
		FILE* fpcmd=popen("nvidia-smi -L 2>/dev/null", "r");
		if(!fpcmd){
			NGPU=0;
		} else{
			char line[4096];
			while(fgets(line, sizeof(line), fpcmd)){
				if(!strncmp(line, "GPU", 3)){
					NGPU++;
				} else{
					warning_time("unknown entry: %s\n", line);
				}
			}
		}
		pclose(fpcmd);
		dbg("NGPU=%d\n", NGPU);
	}

	//Use ws_service in a separate thread does not update job status.
	//We call it in scheduler_timeout using poll. no need for LOCK in this case
#if HAS_LWS
	extern int PORT;
	start_lws(PORT+100);
#endif
	//Must acquire mutex_sch before handling run_t
	//still need time out to handle process queue.
	extern int PORT;
	listen_port((listen_opt_t){.port=PORT, .localpath=slocal, .responder=respond, .http_responder=http_handler, .timeout_fun=scheduler_timeout, .timeout_sec=10});
	runned_remove(-1);
	socket_close();
	signal_caught=0;//enable printing memory debug information
	exit(0);
}
