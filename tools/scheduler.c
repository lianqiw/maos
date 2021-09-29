/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#include <tgmath.h>
#include <errno.h>
#include <sys/socket.h>
#include "../sys/sys.h"
#if HAS_LWS
#include "scheduler_ws.h"
#endif
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
	int sock2;//socket for secondary maos connection for maos_command()
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
	struct MONITOR_T* next;
}MONITOR_T;


static MONITOR_T* pmonitor=NULL;
static int all_done=0;
/*A linked list to store running process*/
static RUN_T* running=NULL;/*points to the begining */
static RUN_T* running_end=NULL;/*points to the last to add. */
/*A linked list to store ended process*/
static RUN_T* runned=NULL;
static RUN_T* runned_end=NULL;
static double usage_cpu;
static RUN_T* running_add(int pid, int sock);
static RUN_T* running_get(int pid);
static RUN_T* running_get_by_status(int status);

static RUN_T* runned_get(int pid);
static void runned_remove(int pid);
static int runned_add(RUN_T* irun);
static void running_remove(int pid, int status);
static RUN_T* running_get_by_sock(int sock);
//static MONITOR_T *monitor_get(int hostid);
static void monitor_remove(int hostid);
static MONITOR_T* monitor_add(int hostid);
static void scheduler_timeout(void);
static void monitor_send(RUN_T* run, char* path);
static void monitor_send_initial(MONITOR_T* ic);
static void monitor_send_load(void);
static long counter=-1;//an negative index for a pending run
static int nused_cpu=0;//number of CPUs being used
static int nused_gpu=0;//number of GPUs being used
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
	/*
	  Add to the end of the runned list intead of the beggining.
	*/

	/*move irun to runned list from running list. */
	irun->next=NULL;
	if(runned_end){
		runned_end->next=irun;
		runned_end=irun;
	} else{
		runned_end=runned=irun;
	}
	irun->status.done=1;
	if(irun->status.info<10){
		warning_time("Should be a finished process\n");
		return 1;
	} else{
		return 0;
	}
}
static void runned_remove(int pid){
	RUN_T* irun, * irun2=NULL;
	int removed=0;
	for(irun=runned; irun; irun2=irun, irun=irun->next){
		if(irun->pid==pid){
			irun->status.info=S_REMOVE;
			monitor_send(irun, NULL);
			if(irun2){
				irun2->next=irun->next;
				if(irun->next==NULL)
					runned_end=irun2;
			} else{
				runned=irun->next;
				if(irun->next==NULL)
					runned_end=runned;
			}
			free(irun->path);
			free(irun->path0);
			free(irun->exe);
			free(irun);
			removed=1;
			break;
		}
	}
	if(!removed){
		warning_time("runned_remove: Record %s:%d not found!\n", HOST, pid);
	}
}
static RUN_T* runned_get(int pid){
	RUN_T* irun;
	for(irun=runned; irun; irun=irun->next){
		if(irun->pid==pid){
			break;
		}
	}
	return irun;
}
/**
   Restart a crashed/finished job
*/
static void runned_restart(int pid){
	RUN_T* irun, * irun_prev=NULL;
	for(irun=runned; irun; irun_prev=irun, irun=irun->next){
		if(irun->pid==pid){
			if(irun->status.info<10){
				dbg_time("status.info=%d\n", irun->status.info);
				break;
			}
			if(irun_prev){//not start of list
				irun_prev->next=irun->next;
			} else{//start if list
				irun_prev=runned=irun->next;
			}
			if(!irun->next){
				runned_end=irun_prev;
			}
			//Insert to the beginning of running list
			irun->next=running;
			running=irun;
			if(!running_end){
				running_end=irun;
			}
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
	if((irun=running_get(pid))){
		dbg2_time("PID %d is already in running.\n", pid); /*create the node */
		if(irun->sock!=sock) irun->sock=sock;
		return irun;
	} else{
		char progname[PATH_MAX];
		if(pid>0){
			if(get_job_progname(progname, PATH_MAX, pid)){//failed. Already exited.
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
		if(running_end){/*list is not empty*/
			running_end->next=irun;
			running_end=irun;
		} else{
			running_end=running=irun;
		}
		return irun;
	}
}

static void running_remove(int pid, int status){
	//dbg_time("Removing %d from running\n", pid);
	RUN_T* irun, * irun2=NULL;
	for(irun=running; irun; irun2=irun, irun=irun->next){
		if(irun->pid==pid){
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
			dbg_time("remove job %d with status %d\n", pid, irun->status.info);
			if(irun2){
				if(!(irun2->next=irun->next)){
					running_end=irun2;
				}
			} else{//beginning of list
				if(!(running=irun->next)){
					running_end=running;
				}
			}
			monitor_send(irun, NULL);
			runned_add(irun);
			break;
		}
	}
	if(!irun){
		dbg_time("%s:%d not found\n", HOST, pid);
	}
}

static RUN_T* running_get(int pid){
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
	int jrun=0;
	for(irun=running; irun; irun=irun->next){
		jrun++;
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
/**
	check all the jobs. remove if any job quited.
 */
static void check_jobs(void){
	RUN_T* irun, * irun2;
	int nrunning=0;
	if(running){
		time_t now=myclocki();
		for(irun=running; irun; irun=irun2){
			irun2=irun->next;
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
	int avail=get_cpu_avail();
	if(avail<1){
		dbg_time("process_queue: no CPUs are available\n");
		return;
	}
	dbg_time("nused_cpu=%d, ngpu=%d, avail=%d\n", nrun_get(0), nrun_get(1), avail);
	RUN_T* irun=running_get_by_status(S_WAIT);
	if(irun){//There are jobs waiting.
		if(irun->sock>0){//already connected.
			int nthread=irun->nthread;
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
		if(avail>1&&nrun_get(0)<NTHREAD&&(!NGPU||nrun_get(1)<NGPU)){//resource available to star a new job
			static double lasttime2=0;
			time_t thistime2=myclocki();
			if(thistime>lasttime2+0.001){
				lasttime2=thistime2;
				irun=running_get_by_status(S_QUEUED);
				if(!irun){
					dbg_time("all jobs are done\n");
					all_done=1;
					nrun_handle(3, 0, 0, 0);
					counter=-1; //reset the counter
				} else{
					int pid;
					irun->last_time=thistime2;
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
			}
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
   Pass socket and command to maos.
 */
static int maos_command(int pid, int sock, int cmd){
	RUN_T* irun=running_get(pid);
	int cmd2[2]={cmd, 0};
	if(!irun||!irun->sock2||(stwriteintarr(irun->sock2, cmd2, 2)||stwritefd(irun->sock2, sock))){
		warning_time("Unable to pass socket (%d) to maos (%d, %d)\n", sock, irun?irun->pid:-1, irun?irun->sock2:-1);
		stwriteint(sock, -1);//respond failure message.
	} else{
		dbg_time("Successfully passed socket (%d) to maos (%d, %d)\n", sock, irun->pid, irun->sock2);
		stwriteint(sock, 0);//respond succeed
	}
	return -1;//do not keep this connection.
}
typedef struct SOCKID_M{
	int id;
	int sock;
	struct SOCKID_M* prev;
	struct SOCKID_M* next;
} SOCKID_T;
static SOCKID_T* shead=0;
//save_socket to linked list. socket is uniq, but id may not be.
static void socket_save(int sock_save, int id){
	int found=0;
	for(SOCKID_T* p=shead; p; p=p->next){
		if(p->sock==sock_save){
			found=1;
			p->id=id;
		}
	}
	if(!found){
		SOCKID_T* tmp=(SOCKID_T*)malloc(sizeof(SOCKID_T));
		tmp->id=id;
		tmp->sock=sock_save;
		tmp->next=shead;
		tmp->prev=0;
		if(shead){
			shead->prev=tmp;
		}
		shead=tmp;
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
	for(SOCKID_T* p_next, *p=shead; p&&sock_save==-1; p=p_next){
		p_next=p->next;
		int badsock=0;
		if((badsock=stcheck(p->sock))||p->id==id){
			if(badsock){
				close(p->sock); //closed socket
			} else{
				sock_save=p->sock; //valid match
			}
			socket_remove(p);
		}
	}
	return sock_save;
}
static void socket_heartbeat(){
	for(SOCKID_T* p_next, *p=shead; p; p=p_next){
		p_next=p->next;
		int cmd[1]={DRAW_HEARTBEAT};
		if(stwrite(p->sock, cmd, sizeof(cmd))){
			socket_remove(p);
		}
	}
}

static int scheduler_recv_wait=-1;//>-1: there is pending scheduler_recv_socket.
//PNEW(mutex_sch);//respond() and scheduler_handle_ws() much lock this before preceed.
/**
   respond to client requests. The fixed header is int[2] for opcode and
   pid. Optional data follows based on opcode. The disadvanced of this is that
   new commands with additional payload cannot be easily introduced in the
   scheduler().
*/
static int respond(int sock){
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
		goto end;
	}
	pid=cmd[1];
	dbg3_time("respond %d got %d %d. \n", sock, cmd[0], cmd[1]);
	//LOCK(mutex_sch);
	switch(cmd[0]){
	case CMD_START://1: Called by maos when job starts.
	{
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
			dbg_time("(%d) Job %d waiting to start with nused_cpu=%d, ngpu=%d\n", sock, pid, nthread, ngpu);
			irun->status.info=S_WAIT;
		} else{/*no waiting, no need reply. */
			dbg_time("(%d) Job %d started with nused_cpu=%d, ngpu=%d\n", sock, pid, nthread, ngpu);
			irun->status.info=S_START;
			nrun_add(pid, nthread, ngpu);
		}
		all_done=0;
		if(irun->path) monitor_send(irun, irun->path);
		monitor_send(irun, NULL);
	}
	break;
	case CMD_FINISH://2: Called by MAOS when job finishes.
		dbg_time("(%d) Job %d reports finished\n", sock, pid);
		running_remove(pid, S_FINISH);
		ret=0;//-1. Do not yet close. wait for client to close.
		break;
	case CMD_STATUS://3: Called by MAOS to report status at every time step
	{
		RUN_T* irun=running_get(pid);
		if(!irun){/*started before scheduler is relaunched. */
			dbg_time("(%d) pid=%d is running but not recorded.\n", sock, pid);
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
		dbg_time("(%d) Job %d reports crashed\n", sock, pid);
		running_remove(pid, S_CRASH);
		ret=0;//-1; wait for client to close.
		break;
	case CMD_MONITOR://5: Called by Monitor when it connects
	{
		MONITOR_T* tmp=monitor_add(sock);
		if(tmp){
			if(pid>=0x8){//old scheme where scheduler_version at compiling is passed
				tmp->load=1;
			} else if(pid<0){//new scheme with bit indicating different options
				tmp->load=pid&1;
				tmp->plot=pid&(1<<1);
			}
		}
		dbg_time("(%d) Monitor is connected.\n", sock);
	}
	break;
	case CMD_PATH://6: Called by MAOS to report the PATH.
	{
		RUN_T* irun=running_add(pid, sock);
		if(!irun){
			warning_time("(%d) running_add %d failed. Exe already exited.\n", sock, pid);
			break;
		}
		free(irun->path0);
		free(irun->path);
		if(streadstr(sock, &irun->path0)){
			dbg_time("(%d) Job %d receiving path failed. \n", sock, pid);
			ret=-1;
			break;
		}
		irun->path=remove_endl(irun->path0);
		dbg_time("(%d) Job %d received path. \n", sock, pid);
		monitor_send(irun, irun->path);
	}
	break;
	case CMD_KILL://7: Called by Monitor to kill a task.
	{
		RUN_T* irun=running_get(pid);
		warning_time("(%d) Received monitor command to kill %d\n", sock, pid);
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
	break;
	case CMD_TRACE://8: Called by MAOS to request a backtrace
	{
		char* buf=NULL, out[10000];
		if(streadstr(sock, &buf)
			||call_addr2line(out, 10000, buf)
			||stwritestr(sock, out)){
			warning_time("CMD_TRACE failed. buf=%s, out=%s\n", buf, out);
			ret=-1;
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
	case CMD_SOCK://10:Called by draw() to cache or request a fd. Valid over UNIX socket only.
	{
		if(pid>0){//receive sock from draw()
			int sock_save;
			if(streadfd(sock, &sock_save)){
				warning_time("(%d) receive socket failed\n", sock);
				sock_save=-1;
				ret=-1;
			} else{
				dbg_time("(%d) received socket %d\n", sock, sock_save);
				socket_save(sock_save, pid);
			}
		} else if(pid<0){//send existing sock to draw()
			int sock_save=socket_get(-pid);//drawdaemon with the same session id.
			if(sock_save==-1){
				sock_save=socket_get(0);//available drawdaemon with no session id.
			}
			dbg_time("(%d) received socket request, sock_saved=%d\n", sock, sock_save);

			//cannot pass -1 as sock, so return a flag first. sock can be zero.
			if(stwriteint(sock, sock_save>-1?0:-1)){
				warning_time("(%d) Unable to talk to draw\n", sock);
				ret=-1;
			}
			if(sock_save>-1){
				if(stwritefd(sock, sock_save)){
					warning_time("(%d) send socket %d failed\n", sock, sock_save);
					ret=-1;//close connection to draw()
				} else{//socket is transferred to draw. we close it.
					dbg_time("(%d) send socket %d success\n", sock, sock_save);
				}
				close(sock_save);
			}
		} else{//pid==0; request a drawdaemon using monitor. 
			MONITOR_T* pm=0;
			for(pm=pmonitor; pm; pm=pm->next){
				if(pm->plot){
					int moncmd[3]={MON_DRAWDAEMON,0,0};
					stwrite(pm->sock, moncmd, sizeof(moncmd));
					scheduler_recv_wait=sock;
					dbg_time("(%d) request monitor (%d) to start drawdaemon\n", sock, pm->sock);
					break;
				}
			}
			if(!pm){
				//there is no available drawdameon. Need to create one by sending request to monitor.
				warning_time("(%d) there is no monitor available to start drawdaemon\n", sock);
				if(stwriteint(sock, -1)){
					warning_time("(%d) Failed to respond to draw.\n", sock);
				}
			}
		}
	}
	break;
	case CMD_REMOVE://11: Called by Monitor to remove a finished job fron the list*/
	{
		RUN_T* irun=runned_get(pid);
		if(irun){
			runned_remove(pid);
		} else{
			dbg_time("(%d) CMD_REMOVE: %s:%d not found\n", sock, HOST, pid);
		}
	}
	break;
	case CMD_DISPLAY://12: called by monitor enable connection of drawdaemon to draw().*/
		dbg_time("(%d) CMD_DISPLAY received with pid=%d\n", sock, pid);
		if(pid<=0){//this is for pending scheduler_recv_socket
			if(scheduler_recv_wait==-1||stwriteint(scheduler_recv_wait, 0)
				||stwritefd(scheduler_recv_wait, sock)){
				dbg_time("(%d) Failed to pass sock to draw at %d, save socket for future\n", sock, scheduler_recv_wait);
				socket_save(dup(sock), abs(pid));//duplicate socket and keep it 
				ret=-1;//prevent scheduler from listening to this socket.
			} else{
				dbg_time("(%d) passed sock to draw at %d\n", sock, scheduler_recv_wait);
				ret=-1;//close socket on scheduler.
			}
			stwriteint(sock, 0);
			scheduler_recv_wait=-1;
		} else{
			ret=maos_command(pid, sock, MAOS_DRAW);
		}
		break;
	case CMD_MAOS://13: create a live link to maos.
		dbg_time("(%d) pass command to maos %d\n", sock, pid);
		ret=maos_command(pid, sock, MAOS_VAR);
		break;
	case CMD_MAOSDAEMON://14; called by maos to save a port to run maos_command
	{
		RUN_T* irun=running_get(pid);
		if(irun){
			irun->sock2=sock;
			dbg_time("(%d) maos socket is saved\n", sock);
		} else{
			warning_time("(%d) maos irun not found\n", sock);
			ret=-1;
		}
	}
	break;
	case CMD_RESTART://15: Called by monitor to restart a job
		dbg_time("(%d) restart job %d\n", sock, pid);
		runned_restart(pid);
		break;
	case CMD_KILLED://16: called by maos to indicate that job is cancelled or killed
		running_remove(pid, S_KILLED);
		ret=0;//-1; //wait for client to close.
		break;
	case CMD_UNUSED4://17;not used
		break;
	case CMD_LAUNCH://18: called from maos from another machine to start a job in this machine
	{
		char* exename=NULL;
		char* execwd=NULL;
		char* execmd=NULL;
		if(pid>=2){
			if(streadstr(sock, &exename)){
				warning_time("(%d) Unable to read exename.\n", sock);
				ret=-1;
			}
		}
		if(ret!=-1&&pid>=3){//old method. should not be used.
			if(streadstr(sock, &execwd)){
				warning_time("(%d) Unable to read execwd.\n", sock);
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
			warning_time("(%d) Unable to read execmd\n", sock);
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
		warning_time("(%d) Invalid cmd: %x\n", sock, cmd[0]);
		ret=-1;
	}
	//job is finished. process the next job. don't wait for timeout
	if((cmd[0]==CMD_KILLED||cmd[0]==CMD_FINISH||cmd[0]==CMD_CRASH)&&running){
		process_queue();
	}
end:
	if(ret){
		RUN_T* irun=running_get_by_sock(sock);//is maos
		if(irun&&irun->status.info<10){
			//connection failed to a running maos job.
			if(kill(irun->pid, 0)){
				dbg_time("(%d) Job %d no longer exists, crashed?\n", sock, irun->pid);
				running_remove(irun->pid, S_CRASH);
			}
		}
		ret=-1;
	}
#if HAS_LWS
	ws_service(0);
#endif
	return ret;//ret=-1 will close the socket.
}
/*
  handle requests from web browser via websockets. Browser sends request over
  text messages with fields separated by '&'.
*/
void scheduler_handle_ws(char* in, size_t len){
	(void)len;
	char* sep=strchr(in, '&');
	if(!sep){
		warning_time("Invalid message: %s\n", in);
		return;
	}
	*sep='\0';
	int pid=strtol(in, 0, 10);
	sep++;
	char* tmp;
	if((tmp=strchr(sep, '&'))){
		*tmp='\0';
	}
	if((tmp=strchr(sep, ';'))){
		*tmp='\0';
	}
	//LOCK(mutex_sch);
	if(!strcmp(sep, "REMOVE")){
		RUN_T* irun=runned_get(pid);
		if(irun){
			runned_remove(pid);
		} else{
			warning_time("CMD_REMOVE: %s:%d not found\n", HOST, pid);
		}
	} else if(!strcmp(sep, "KILL")){
		RUN_T* irun=running_get(pid);
		if(irun){
			if(irun->status.info!=S_QUEUED){
				kill(pid, SIGTERM);
				if(irun->status.info==S_WAIT){//wait up the process.
					stwriteint(irun->sock, S_START);
				}
			} else{
				running_remove(pid, S_KILLED);
			}
			warning_time("HTML client send term signal to %5d term signal.\n", pid);
		}
	} else{
		warning_time("Unknown action: %s\n", sep);
	}
	//UNLOCK(mutex_sch);
}
static void scheduler_timeout(void){
	static time_t lasttime3=0;//every 3 seconds
	static time_t lasttime10=0;//every 10 seconds

	time_t thistime=myclocki();
	//Process job queue
	if(running){
		process_queue();
	}
	//respond to heatbeat of saved sockets
#if HAS_LWS
	ws_service(1000);//service http.
#endif	
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
static MONITOR_T* monitor_add(int sock){
	/*dbg_time("added monitor on sock %d\n",sock); */
	MONITOR_T* node=mycalloc(1, MONITOR_T);
	node->sock=sock;
	node->next=pmonitor;
	pmonitor=node;
	monitor_send_initial(node);
	return pmonitor;
}
static void monitor_remove(int sock){
	MONITOR_T* ic, * ic2=NULL;
	/*int freed=0; */
	for(ic=pmonitor; ic; ic2=ic, ic=ic->next){
		if(ic->sock==sock){
			if(ic2){
				ic2->next=ic->next;
			} else{
				pmonitor=ic->next;
			}
			//close(ic->sock); close panics accept
			dbg_time("Removed monitor at %d\n", ic->sock);
			shutdown(ic->sock, SHUT_WR);
			free(ic);
			break;
		}
	}
}

/**
   Convert RUN_T to string for websocket.
*/
#if HAS_LWS
static void html_push(RUN_T* irun, char* path, char** dest, size_t* plen, long prepad, long postpad){
	char temp[4096]={0};
	size_t len;
	if(path){
		len=snprintf(temp, 4096, "%d&PATH&%s$", irun->pid, path);
	} else{
		status_t* st=&irun->status;
		struct tm* tim=localtime(&(st->timlast));
		char stime[80];
		strftime(stime, 80, "[%a %H:%M:%S]", tim);
		len=snprintf(temp, 4096, "%d&STATUS&%d&%d" /*pid, key, pidnew, status*/
			"&%s&%.2f&%.2f" /*start time, errhi, errlo*/
			"&%d&%d&%d&%d" /*iseed, nseed, isim, nsim*/
			"&%ld&%ld&%.3f$" /*rest, tot, step timing*/
			, irun->pid, irun->pidnew, st->info,
			stime, st->clerrhi, st->clerrlo,
			st->nseed==0?0:st->iseed+1, st->nseed, st->simend==0?0:st->isim+1, st->simend,
			st->rest, st->laps+st->rest, st->tot*st->scale);
	}
	if(!dest){
		ws_push(temp, len);
	} else{
		*plen=len;
		*dest=(char*)malloc(len+prepad+postpad);
		memcpy(*dest+prepad, temp, len);
	}
}
static void html_push_all_do(RUN_T* irun, l_message** head, l_message** tail, long prepad, long postpad){
	l_message* node=(l_message*)calloc(1, sizeof(l_message));
	l_message* node2=(l_message*)calloc(1, sizeof(l_message));
	node->next=node2;
	node2->next=0;
	if(*tail){
		(*tail)->next=node;
	} else{
		*head=node;
	}
	*tail=node2;
	html_push(irun, irun->path, &node->payload, &node->len, prepad, postpad);
	html_push(irun, 0, &node2->payload, &node2->len, prepad, postpad);
}
//called by scheduelr_ws upon client connection. 
void html_push_all(l_message** head, l_message** tail, long prepad, long postpad){
	//LOCK(mutex_sch);
	*head=*tail=0;
	for(RUN_T* irun=runned; irun; irun=irun->next){
		html_push_all_do(irun, head, tail, prepad, postpad);
	}
	for(RUN_T* irun=running; irun; irun=irun->next){
		html_push_all_do(irun, head, tail, prepad, postpad);
	}
	//UNLOCK(mutex_sch);
}
#endif
static int monitor_send_do(RUN_T* irun, char* path, int sock){
	int cmd[3];
	cmd[1]=irun->pidnew;
	cmd[2]=irun->pid;
	if(path){/*don't do both. */
		cmd[0]=MON_PATH;
		return (stwrite(sock, cmd, 3*sizeof(int))||stwritestr(sock, path));
	} else{
		cmd[0]=MON_STATUS;
		return (stwrite(sock, cmd, 3*sizeof(int))||stwrite(sock, &irun->status, sizeof(status_t)));
	}
}
/* Notify alreadyed connected monitors job update. */
static void monitor_send(RUN_T* irun, char* path){
#if HAS_LWS
	html_push(irun, path, 0, 0, 0, 0);
#endif
	MONITOR_T* ic;
redo:
	for(ic=pmonitor; ic; ic=ic->next){
		int sock=ic->sock;
		if(monitor_send_do(irun, path, sock)){
			monitor_remove(sock);
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
	if(NGPU&&0){
		gpu=get_usage_gpu(&gpu_mem);
	}
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
	int sock;
	RUN_T* irun;
	sock=ic->sock;
	for(irun=runned; irun; irun=irun->next){
		if(monitor_send_do(irun, irun->path, sock)||monitor_send_do(irun, NULL, sock)){
			monitor_remove(sock);
			return;
		}
	}
	for(irun=running; irun; irun=irun->next){
		if(monitor_send_do(irun, irun->path, sock)||monitor_send_do(irun, NULL, sock)){
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
		if((fd=lock_file(fnlock, 0, scheduler_version))<0){
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
	ws_start(PORT+100);
#endif	
	//Must acquire mutex_sch before handling run_t
	//still need time out to handle process queue.
	extern int PORT;
	listen_port(PORT, slocal, respond, 1, scheduler_timeout, 0);
	remove(slocal);

	exit(0);
}
