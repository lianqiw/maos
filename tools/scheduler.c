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
	status_t status;
	double launch_time;
	double last_time; //time of last report
	int pid;
	int pidnew;//the new pid
	int sock;//socket for primary maos connection
	int sock2;//socket for secondary maos connection for maos_command()
	int nthread;//number of threads
	int ngpu;//number of gpus requested.
	char* exe; /*Path to the executable.*/
	char* path;/*Job path and Job arguments.*/
	char* path0;/*same as path, with fields separated by \n instead of space*/
}RUN_T;

/**
  Struct to hold available monitors waiting for information.
*/
typedef struct MONITOR_T{
	int sock;
	int load;/*handle machine load information. */
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
static RUN_T* running_get_wait(int status);

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
static int nrun_handle(int cmd, int pid, int nthread, int ngpu_used){
	static int ncpu=0;
	static int ngpu=0;
	switch(cmd){
	case 0:
		if(pid==1){
			return ngpu;
		} else{
			return ncpu;
		}
		break;
	case 1:
		ncpu+=nthread;
		ngpu+=ngpu_used;
		dbg2_time("%d: ncpu %d->%d. ngpu %d->%d.\n",
			pid, ncpu-nthread, ncpu, ngpu-ngpu_used, ngpu);
		break;
	case 2:
		ncpu-=nthread;
		ngpu-=ngpu_used;
		dbg2_time("%d: ncpu %d->%d. ngpu %d->%d.\n",
			pid, ncpu+nthread, ncpu, ngpu+ngpu_used, ngpu);
		if(ncpu<0){
			dbg2_time("ncpu=%d\n", ncpu);
			ncpu=0;
		}
		if(ngpu<0){
			dbg2_time("ngpu=%d\n", ngpu);
			ngpu=0;
		}
		break;
	}
	return ncpu;
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
		if(pid>0){
			irun->exe=strdup(progname);
		}
		/*record the launch time */
		if(pid>0){
			irun->launch_time=get_job_launchtime(pid);
		} else{
			irun->launch_time=INFINITY;
		}
		irun->next=NULL;
		if(running_end){/*list is not empty*/
			if(pid<=0||running_end->launch_time<=irun->launch_time){
			/*append the node to the end. */
				running_end->next=irun;
				running_end=irun;
			} else{
			/*insert the node in the middle. */
				RUN_T* jrun, * jrun2=NULL;
				for(jrun=running; jrun; jrun2=jrun, jrun=jrun->next){
					if(jrun->launch_time>=irun->launch_time){
						irun->next=jrun;
						if(jrun2){
							jrun2->next=irun;
						} else{
							running=irun;
						}
						break;
					}
				}
				if(!jrun){
					warning_time("failed to insert pid %d\n", pid);
				}
			}
		} else{
			running_end=running=irun;
		}
		irun->status.timstart=myclocki();
		return irun;
	}
}

static void running_remove(int pid, int status){
	//dbg_time("Removing %d from running\n", pid);
	RUN_T* irun, * irun2=NULL;
	for(irun=running; irun; irun2=irun, irun=irun->next){
		if(irun->pid==pid){
			if(irun->pid>0&&(irun->status.info==S_START
							 ||irun->status.info==S_RUNNING)){
				nrun_sub(irun->pid, irun->nthread, irun->ngpu);
			}
			if(irun->status.info<10){//mark termination time.
				irun->status.timend=myclocki();
			}
			irun->status.info=status;
			dbg_time("remove job %d with status %d\n", pid, irun->status.info);
			//remove from the running list
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

static RUN_T* running_get_wait(int status){
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
	if(running){
		double now=myclockd();
		for(irun=running; irun; irun=irun2){
			irun2=irun->next;
			if(irun->pid>0){
				if(kill(irun->pid, 0)){
					if(irun->last_time+60<now){//allow grace period.
						dbg_time("check_jobs: Job %d no longer exists\n", irun->pid);
						running_remove(irun->pid, S_CRASH);
					}
				} else if((irun->last_time+600<now)&&irun->status.info==S_RUNNING){
					dbg_time("check_jobs: Job %d does not update after %g seconds.\n",
							irun->pid, now-irun->last_time);
					irun->status.info=S_UNKNOWN;
					monitor_send(irun, NULL);
				}
			}
		}
	}
}
/**
 examine the process queue to start routines once CPU is available.
*/
static void process_queue(void){
	static double lasttime=0;
	if(nrun_get(0)>0&&myclockd()-lasttime<1){
		//dbg_time("process_queue: too son\n");
		return;
	}
	lasttime=myclockd();
	if(nrun_get(0)>=NCPU){
		//dbg_time("process_queue: enough jobs are running\n");
		return;
	}
	int avail=get_cpu_avail();
	//dbg_time("process_queue: nrun=%d avail=%d\n", nrun_get(0), avail);
	if(avail<1){
		//dbg_time("process_queue: no CPUs are available\n");
		return;
	}
	RUN_T* irun=running_get_wait(S_WAIT);
	while(irun&&irun->pid>0&&kill(irun->pid, 0)){//job exited
		running_remove(irun->pid, S_CRASH);
		irun=running_get_wait(S_WAIT);
	}
	if(irun){//There are jobs waiting.
		if(irun->sock>0){//already connected.
			int nthread=irun->nthread;
			if(nrun_get(0)+nthread<=NCPU&&(nthread<=avail||avail>=3)
				&&(!NGPU||!irun->ngpu||nrun_get(1)+irun->ngpu<=NGPU)){//resource available to star the job
				dbg_time("Start %d at %d\n", irun->pid, irun->sock);
				if(stwriteint(irun->sock, S_START)){
					warning_time("stwriteint to %d failed: %s\n", irun->sock, strerror(errno));
				}
				nrun_add(irun->pid, nthread, irun->ngpu);
				irun->status.timstart=myclocki();
				irun->status.info=S_START;
				monitor_send(irun, NULL);
			}
		} else{
			dbg_time("Wait for %d to connect. irun->sock=%d\n", irun->pid, irun->sock);
		}
	} else{
		if(avail>1&&nrun_get(0)<NTHREAD&&(!NGPU||nrun_get(1)<NGPU)){//resource available to star a new job
			static double lasttime2=0;
			double thistime=myclockd();
			if(thistime>lasttime2+0.001){
				lasttime2=thistime;
				irun=running_get_wait(S_QUEUED);
				if(!irun){
					dbg_time("all jobs are done\n");
					all_done=1;
					counter=-1; //reset the counter
				} else{
					int pid;
					if((pid=launch_exe(irun->exe, irun->path0))<0){
						warning_time("start new job failed: %s\n", irun->path);
						running_remove(irun->pid, S_CRASH);
					} else{
						dbg_time("new job is started: %d\n", pid);
						//inplace update the information in monitor
						irun->status.info=S_WAIT;
						irun->status.timstart=myclocki();
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
static void new_job(const char* exename, const char* execmd){
	RUN_T* irun=running_add(--counter, -1);
	if(!irun){
		warning_time("scheduler: running_add failed.\n");
		return;
	}
	irun->status.info=S_QUEUED;
	irun->exe=strdup(exename);
	irun->path0=strdup(execmd);
	irun->path=remove_endl(irun->path0);
	dbg_time("(%s) (%s)\n", exename, execmd);
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
static int scheduler_recv_wait=-1;//>-1: there is pending scheduler_recv_socket.
PNEW(mutex_sch);//respond() and scheduler_handle_ws() much lock this before preceed.
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
	int ret=0, pid, cmd[2];
	//dbg_time("\rrespond %2d start ... ", sock);errno=0;
	if((ret=streadintarr(sock, cmd, 2))){
		//if(errno!=EAGAIN && errno!=ENOENT){//timeout or closed
		dbg_time("read %d failed (%d): %s,ret=%d\n", sock, errno, strerror(errno), ret);
		//}
		goto end;
	}
	pid=cmd[1];
	dbg3_time("respond %d got %d %d. \n", sock, cmd[0], cmd[1]);
	LOCK(mutex_sch);
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
			irun->status.info=S_WAIT;
			all_done=0;
			//dbg_time("%5d queued.\n",pid);
		} else{/*no waiting, no need reply. */
			nrun_add(pid, nthread, ngpu);
			irun->status.info=S_START;
			irun->status.timstart=myclocki();
			//dbg_time("%5d started\n",pid);
		}
		if(irun->path) monitor_send(irun, irun->path);
		monitor_send(irun, NULL);
	}
	break;
	case CMD_FINISH://2: Called by MAOS when job finishes.
		running_remove(pid, S_FINISH);
		ret=-1;
	break;
	case CMD_STATUS://3: Called by MAOS to report status at every time step
	{
		RUN_T* irun=running_get(pid);
		if(!irun){/*started before scheduler is relaunched. */
			dbg_time("pid=%d is running but not recorded.\n", pid);
			irun=running_add(pid, sock);
			if(!irun){
				warning_time("scheduler: running_add %d failed. Exe already exited.\n", pid);
				break;
			}
			irun->status.info=S_START;
			irun->nthread=irun->status.nthread;
			nrun_add(pid, irun->nthread, irun->ngpu);
		}
		if(sizeof(status_t)!=read(sock, &(irun->status), sizeof(status_t))){
			warning_time("Error reading status.\n");
		}
		irun->last_time=myclockd();
		monitor_send(irun, NULL);
	}
	break;
	case CMD_CRASH://4: called by MAOS when job crashes
		running_remove(pid, S_CRASH);
		ret=-1;
		break;
	case CMD_MONITOR://5: Called by Monitor when it connects
	{
		MONITOR_T* tmp=monitor_add(sock);
		if(tmp&&pid>=0x8){/*check monitor version. */
			tmp->load=1;
		}
		//dbg_time("Monitor is connected at sock %d.\n", sock);
	}
	break;
	case CMD_PATH://6: Called by MAOS to report the PATH.
	{
		RUN_T* irun=running_add(pid, sock);
		if(!irun){
			warning_time("scheduler: running_add %d failed. Exe already exited.\n", pid);
			break;
		}
		free(irun->path0);
		free(irun->path);
		if(streadstr(sock, &irun->path0)){
			ret=-1;
			break;
		}
		irun->path=remove_endl(irun->path0);
		dbg_time("Received path: %s\n", irun->path);
		monitor_send(irun, irun->path);
	}
	break;
	case CMD_KILL://7: Called by Monitor to kill a task.
	{
		RUN_T* irun=running_get(pid);
		warning_time("Received monitor command to kill %d\n", pid);
		if(irun){
			if(irun->status.info!=S_QUEUED){
				kill(pid, SIGTERM);
				if(irun->status.info==S_WAIT){//wait up the process.
					stwriteint(irun->sock, S_START);
				}
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
		break;
	case CMD_SOCK://10:Called by draw() to cache or request a fd. Valid over UNIX socket only.
	{
		typedef struct SOCKID_M{
			int id;
			int sock;
			struct SOCKID_M* prev;
			struct SOCKID_M* next;
		} SOCKID_T;
		static SOCKID_T* head=0;
		if(pid>0){//receive sock from draw()
			int found=0;
			int sock_save;
			if(streadfd(sock, &sock_save)){
				warning_time("receive socket from %d failed\n", sock);
				sock_save=-1;
			} else{
				dbg_time("received socket %d from %d\n", sock_save, sock);
				for(SOCKID_T* p=head; p; p=p->next){
					if(p->id==pid){
						close(p->sock);
						p->sock=sock_save;
						found=1;
					}
				}
				if(!found){
					SOCKID_T* tmp=(SOCKID_T*)malloc(sizeof(SOCKID_T));
					tmp->id=pid;
					tmp->sock=sock_save;
					tmp->next=head;
					tmp->prev=0;
					if(head){
						head->prev=tmp;
					}
					head=tmp;
				}
			}
		} else if(pid<0){//send existing sock to draw()
			int sock_save=-1;
			for(SOCKID_T* p_next, *p=head; p&&sock_save==-1; p=p_next){
				p_next=p->next;
				int badsock=0;
				if((badsock=stcheck(p->sock))||(p->id==-pid)){
					if(badsock){
						close(p->sock); //closed socket
					} else{
						sock_save=p->sock; //valid match
					}
					//remove from list 
					if(p->prev){//middle item
						p->prev->next=p->next;
					} else{//first item
						head=p->next;
					}
					if(p->next){//not last item
						p->next->prev=p->prev;
					}
					free(p);
				}
			}
			dbg_time("received socket request from %d, sock_saved=%d\n", sock, sock_save);

			//cannot pass -1 as sock, so return a flag first. sock can be zero.
			if(stwriteint(sock, sock_save>-1?0:-1)){
				warning_time("Unable to talk to draw\n");
			}
			if(sock_save>-1){
				if(stwritefd(sock, sock_save)){
					warning_time("send socket %d to %d failed\n", sock_save, sock);
				} else{//socket is transferred to draw. we close it.
					dbg_time("send socket %d to %d success\n", sock_save, sock);
				}
				close(sock_save);
			}

		} else{//pid==0; request a drawdaemon using monitor. 
			if(pmonitor){
				//there is no available drawdameon. Need to create one by sending request to monitor.
				dbg_time("request monitor to start drawdaemon\n");
				int moncmd[3]={MON_DRAWDAEMON,0,0};
				stwrite(pmonitor->sock, moncmd, sizeof(moncmd));
				scheduler_recv_wait=sock;
				//wait until monitor opend drawdaemon.
				//continue in CMD_DISPLAY
			} else{
				warning_time("there is no monitor available to start drawdaemon\n");
				if(stwriteint(sock, -1)){
					warning_time("Failed to respond to draw at %d.\n", sock);
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
			dbg_time("CMD_REMOVE: %s:%d not found\n", HOST, pid);
		}
	}
	break;
	case CMD_DISPLAY://12: called by monitor enable connection of drawdaemon to draw().*/
		dbg_time("CMD_DISPLAY received from %d with pid=%d\n", sock, pid);
		if(pid==0){//this is for pending scheduler_recv_socket
			if(scheduler_recv_wait==-1||stwriteint(scheduler_recv_wait, 0)
				||stwritefd(scheduler_recv_wait, sock)){
				stwriteint(sock, -1);
				warning_time("Failed to pass sock %d to draw at %d\n", sock, scheduler_recv_wait);
			} else{
				stwriteint(sock, 0);//success.
				dbg_time("passed sock %d to draw at %d\n", sock, scheduler_recv_wait);
			}
			scheduler_recv_wait=-1;
			ret=-1;
		} else{
			ret=maos_command(pid, sock, MAOS_DRAW);
		}
		break;
	case CMD_MAOS://13: create a live link to maos.
		ret=maos_command(pid, sock, MAOS_VAR);
		break;
	case CMD_MAOSDAEMON://14; called by maos to save a port to run maos_command
	{
		RUN_T* irun=running_get(pid);
		if(irun){
			irun->sock2=sock;
			dbg_time("CMD_MAOSDAEMON: socket %d is saved\n", sock);
		} else{
			warning_time("CMD_MAOSDAEMON: irun not found\n");
			ret=-1;
		}
	}
	break;
	case CMD_RESTART://15: Called by maos to restart a job
		runned_restart(pid);
		break;
	case CMD_KILLED://16: called by maos to indicate that job is cancelled or killed
		running_remove(pid, S_KILLED);
		ret=-1;
		break;
	case CMD_UNUSED4://17;not used
		break;
	case CMD_LAUNCH://18: called from maos another machine to start a job in this machine
	{
		char* exename=NULL;
		char* execwd=NULL;
		char* execmd=NULL;
		if(pid>=2){
			if(streadstr(sock, &exename)){
				warning_time("Unable to read exename.\n");
			}
		}
		if(pid>=3){//old method. should not be used.
			if(streadstr(sock, &execwd)){
				warning_time("Unable to read execwd.\n");
			}
		}
		if(!streadstr(sock, &execmd)){
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
				new_job(exename, execmd);
				ret=0;
			}
		} else{
			warning_time("Unable to read execmd\n");
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
		warning_time("Invalid cmd: %x\n", cmd[0]);
		ret=-1;
	}
	cmd[0]=-1;
	cmd[1]=-1;
end:
	if(ret){
		RUN_T* irun=running_get_by_sock(sock);//is maos
		if(irun&&irun->status.info<10){
			//connection failed to a running maos job.
			sleep(1);
			int pid2=irun->pid;
			if(kill(pid2, 0)){
				running_remove(pid2, S_CRASH);
			}
		}
		ret=-1;
	}
	sync();
	UNLOCK(mutex_sch);
	return ret;/*don't close the port yet. may be reused by the client. */
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
	LOCK(mutex_sch);
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
	UNLOCK(mutex_sch);
}
static void scheduler_timeout(void){
	static double lasttime3=0;//every 3 seconds
	static double lasttime10=0;//every 10 seconds

	double thistime=myclockd();
	//Process job queue
	if(!all_done){
		process_queue();
	}

	if(thistime>=(lasttime3+3)){
		if(running){
			check_jobs();
		}
		lasttime3=thistime;
	}
	//Report CPU usage every 3 seconds. 
	if(thistime>=(lasttime10+10)){
		usage_cpu=get_usage_cpu();
		monitor_send_load();
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
void html_convert(RUN_T* irun, char* path, char** dest, size_t* plen, long prepad, long postpad){
	char temp[4096]={0};
	size_t len;
	if(path){
		len=snprintf(temp, 4096, "%d&PATH&%s;", irun->pid, path);
	} else{
		status_t* st=&irun->status;
		struct tm* tim=localtime(&(st->timstart));
		char stime[80];
		strftime(stime, 80, "[%a %H:%M:%S]", tim);
		len=snprintf(temp, 4096, "%d&STATUS&%d&%d" /*pid, key, pidnew, status*/
			"&%s&%.2f&%.2f" /*start time, errhi, errlo*/
			"&%d&%d&%d&%d" /*iseed, nseed, isim, nsim*/
			"&%ld&%ld&%.3f;" /*rest, tot, step timing*/
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
static void html_convert_all_do(RUN_T* irun, l_message** head, l_message** tail, long prepad, long postpad){
	l_message* node=(l_message*)malloc(sizeof(l_message));
	l_message* node2=(l_message*)malloc(sizeof(l_message));
	node->next=node2;
	node2->next=0;
	html_convert(irun, irun->path, &node->payload, &node->len, prepad, postpad);
	html_convert(irun, 0, &node2->payload, &node2->len, prepad, postpad);
	if(*tail){
		(*tail)->next=node;
	} else{
		*head=node;
	}
	*tail=node2;
}

void html_convert_all(l_message** head, l_message** tail, long prepad, long postpad){
	LOCK(mutex_sch);
	*head=*tail=0;
	for(RUN_T* irun=runned; irun; irun=irun->next){
		html_convert_all_do(irun, head, tail, prepad, postpad);
	}
	for(RUN_T* irun=running; irun; irun=irun->next){
		html_convert_all_do(irun, head, tail, prepad, postpad);
	}
	UNLOCK(mutex_sch);
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
	html_convert(irun, path, 0, 0, 0, 0);
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

int main(){
	char slocal[PATH_MAX];//local port
	snprintf(slocal, PATH_MAX, "%s/scheduler", TEMP);
	remove(slocal);
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
#if HAS_LWS
	//launch in a separate thread to avoid excessive idle wakeup.
	//Must acquire mutex_sch before handling run_t
	thread_new(ws_service, (void*)(long)(PORT+100));
#endif
	//still need time out to handle process queue.
	listen_port(PORT, slocal, respond, 1, scheduler_timeout, 0);
	remove(slocal);

	exit(0);
}
