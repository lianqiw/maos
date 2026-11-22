/*
  Copyright 2009-2026 Lianqi Wang
  
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
   \file scheduler_job.c Job list and CPU/GPU accounting of the scheduler. It
   does not do any socket communication except to wake up a waiting job, and
   notifies monitors through monitor_send() that is provided by scheduler.c.
*/
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#include "scheduler_job.h"

/**number of GPUs in this machine*/
int NGPU=0;
/**A linked list to store running process. points to the begining*/
static RUN_T* running=NULL;
/**points to the last to add.*/
static RUN_T** running_end=&running;
/**A linked list to store ended process*/
static RUN_T* runned=NULL;
/**points to the pointer of end slot*/
static RUN_T** runned_end=&runned;
static RUN_T* running_get_by_status(int status);
static int runned_add(RUN_T* irun);
/**an negative index for a pending run*/
static int counter=-1;
/**number of CPUs being used*/
static int nused_cpu=0;
/**number of GPUs being used*/
static int nused_gpu=0;
/**
   Query or update the number of CPUs and GPUs used by running jobs.

   @param[in] cmd 0: query, 1: add resource usage, 2: remove resource usage, 3: reset to zero.
   @param[in] pid For cmd 0, 1 queries GPUs and other values query CPUs. Otherwise the job pid used for logging.
   @param[in] nthread Number of CPUs to add or remove.
   @param[in] ngpu Number of GPUs to add or remove.
   @return For cmd 0 with pid 1, the number of GPUs in use; otherwise the number of CPUs in use.
*/
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
/**
   Get the number of CPUs or GPUs in use.

   @param[in] type 0: CPUs, 1: GPUs.
   @return Number of CPUs or GPUs in use.
*/
int nrun_get(int type){
	return nrun_handle(0, type, 0, 0);
}
/**
   Record that a job started using its nthread CPUs and ngpu GPUs. A job is counted at most once.

   @param[in,out] irun Job record.
*/
void nrun_add(RUN_T *irun){
	if(!irun->counted){
		irun->counted=1;
		nrun_handle(1, irun->pid, irun->nthread, irun->ngpu);
	}
}
/**
   Record that a job released its CPUs and GPUs if they are counted. Counts are clipped at zero.

   @param[in,out] irun Job record.
*/
void nrun_sub(RUN_T *irun){
	if(irun->counted){
		irun->counted=0;
		nrun_handle(2, irun->pid, irun->nthread, irun->ngpu);
	}
}
/**
   Append irun to the end of the "runned" linked list of finished jobs and mark
   it done. The runned_* routines maintain this list. irun must already be
   unlinked from the running list.

   @param[in,out] irun Job record. Ownership moves to the runned list.
   @return 0 if irun has a finished status, 1 (with a warning) otherwise.
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
/**
   Remove the finished job pid (all finished jobs if pid is -1) from the runned
   list, notify monitors, and free the records.

   @param[in] pid Job pid, or -1 for all finished jobs.
*/
void runned_remove(int pid){
	RUN_T *irun=NULL;
	int found=0;
	for(RUN_T **curr=&runned; *curr;){
		irun=*curr;
		if(irun->pid==pid || pid==-1){//remove node
			found=1;
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
	if(!found && pid!=-1){
		warning_time("runned_remove: Record %s:%d not found!\n", HOST, pid);
	}
}
/**
   Restart a crashed or finished job by moving it from the runned list to the
   beginning of the running list with queued status and a new negative pid.

   @param[in] pid Pid of the finished job.
*/
void runned_restart(int pid){
	RUN_T *irun=NULL;
	int found=0;
	for(RUN_T **curr=&runned; *curr;){
		irun=*curr;
		if(irun->pid==pid){
			found=1;
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
			irun->sock=-1;
			irun->sock_cmd=-1;
			break;
		}else{
			curr=&irun->next;
		}
	}
	if(!found){
		warning_time("Process with pid=%d is not found\n", pid);
	}
}
/**
   Add a job to the end of the "running" linked list, which contains queued and
   running jobs. The running_* routines operate on this list. If the job already
   exists, only its socket is updated.

   @param[in] pid Job pid. Positive for a launched process, negative for a queued job.
   @param[in] sock Socket connected to the job, or -1.
   @return The existing or newly allocated record, owned by the running list, or NULL if a positive pid has already exited.
*/
RUN_T* running_add(int pid, int sock){
	RUN_T* irun;
	if((irun=running_get_by_pid(pid))){
		//dbg_time("PID %d is already in running.\n", pid); /*create the node */
		if(irun->sock!=sock) irun->sock=sock;
		return irun;
	} else{
		char progname[PATH_MAX];
		if(pid>0){
			if(get_job_progname(progname, PATH_MAX, pid)){//failed. Already exited.
				dbg_time("PID %d already exited. will not add\n", pid);
				return NULL;
			}
		}
		dbg2_time("adding %d to running\n", pid); /*create the node */
		irun=mycalloc(1, RUN_T);
		irun->pidnew=irun->pid=pid;
		irun->sock=sock;
		irun->sock_cmd=-1;
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

/**
   Move job pid from the running list to the runned list with the given final
   status, releasing its CPU/GPU usage and notifying monitors.

   @param[in] pid Job pid.
   @param[in] status Final status, e.g. S_FINISH, S_CRASH, or S_KILLED.
*/
void running_remove(int pid, int status){
	//dbg_time("Removing %d from running\n", pid);
	RUN_T* irun=NULL;
	int found=0;
	for(RUN_T** curr=&running; *curr; ){
		irun=*curr;
		if(irun->pid==pid){
			found=1;
			if(running_end==&irun->next){
				running_end=curr;
			}
			*curr=irun->next;
			nrun_sub(irun);//nrun_add is called when state changes to S_START
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
	if(!found){
		dbg_time("%s:%d not found\n", HOST, pid);
	}
}

/**
   Find a job in the running list by pid.

   @param[in] pid Job pid. 0 matches the first job with a positive pid.
   @return Reference to the record in the running list, or NULL if not found.
*/
RUN_T* running_get_by_pid(int pid){
	RUN_T* irun;
	for(irun=running; irun; irun=irun->next){
		if(irun->pid==pid||(irun->pid>0&&pid==0)){
			break;
		}
	}
	return irun;
}

/**
   Find the first job in the running list with the given status.

   @param[in] status Status value such as S_QUEUED or S_WAIT.
   @return Reference to the record in the running list, or NULL if not found.
*/
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

/**
   Find a job in the running list by its connection socket.

   @param[in] sock Socket connected to the job.
   @return Reference to the record in the running list, or NULL if not found.
*/
RUN_T* running_get_by_sock(int sock){
	RUN_T* irun;
	for(irun=running; irun; irun=irun->next){
		if(irun->sock==sock){
			break;
		}
	}
	return irun;
}
/**
   Kill job pid: send SIGTERM if it has been started (waking it up first if it is
   waiting) and mark it as unknown until it exits, or remove it from the running
   list as killed if it is still queued.

   @param[in] pid Job pid.
*/
void running_kill(int pid){
	RUN_T* irun=running_get_by_pid(pid);
	if(irun){
		if(irun->status.info!=S_QUEUED){
			if(pid>0) kill(pid, SIGTERM);//a negative pid would signal a process group
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
   Check whether a launched job still exists. A job whose process is gone for
   more than 60 seconds is marked as crashed and moved to the runned list. A job
   that has not reported for 600 seconds is marked as unknown.

   @param[in,out] irun Job record. May be moved to the runned list.
   @return 1 if the job process is still considered alive, 0 otherwise.
*/
int check_job_exists(RUN_T *irun){
	int ans=0;
	if(irun && irun->pid>0 && irun->status.info<10){
		time_t now=myclocki();
		if(kill(irun->pid, 0)){//No longer exists
			if(irun->last_time+60<now){//allow grace period.
				info_time("check_jobs: Job %d no longer exists. Change status to S_CRASH\n", irun->pid);
				running_remove(irun->pid, S_CRASH);
			}else{//wait 1 minute before forceful cleanup
				ans=1;
			}
		} else {
			ans=1;
			if((irun->last_time+600<now)){
				//do not attempt to kill as there maybe a networking issue.
				info_time("check_jobs: Job %d does not update after %lu seconds. Change status to S_UNKNOWN\n",
						irun->pid, now-irun->last_time);
				irun->status.info=S_UNKNOWN;
				monitor_send(irun, NULL);
			}
		}
	}
	return ans;
}
/**
	check all the jobs. remove if any job quited.
 */
void check_jobs(void){
	RUN_T* irun=NULL, *irun_next=NULL;
	int nrunning=0;
	if(running){
		for(irun=running; irun; irun=irun_next){
			irun_next=irun->next;
			if(irun->pid>0){//Running job
				nrunning+=check_job_exists(irun)?1:0;
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
void process_queue(void){
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
		if(irun->sock>=0){//already connected. we start it
			int nthread=irun->nthread;
			int avail=get_cpu_avail();
			if(nrun_get(0)+nthread<=NCPU&&(nthread<=avail||avail>=3)
				&&(!NGPU||!irun->ngpu||nrun_get(1)+irun->ngpu<=NGPU)){//resource available to star the job
				irun->last_time=myclocki();
				irun->status.timlast=myclocki();
				if(stwriteint(irun->sock, S_START)){
					warning_time("Starting Job %d at %d failed: %s\n", irun->pid, irun->sock, strerror(errno));
				} else{
					info_time("Starting Job %d at %d ok.\n", irun->pid, irun->sock);
				}
				//we mark the sate as running in either state and let check_job do the clean up
				irun->status.info=S_START;
				nrun_add(irun);
				monitor_send(irun, NULL);
			}
		} else{//check job status
			check_job_exists(irun);
		}
	}
	if(nrun_get(0)<NTHREAD&&(!NGPU||nrun_get(1)<NGPU)){//resource available to star a new job
		//static double lasttime2=0;
		//time_t thistime2=myclocki();
		//if(thistime2>lasttime2+0.001){//wait 1ms for job to connect.
			//lasttime2=thistime2;
			irun=running_get_by_status(S_QUEUED);
			if(!irun){
				dbg_time("all jobs are done\n");
				//do not reset nused_cpu or counter: other jobs may still be running.
			} else{
				int pid;
				irun->last_time=myclocki();
				irun->status.timlast=myclocki();
				if((pid=launch_exe(irun->exe, irun->path0))<0){
					warning_time("Launch job %d failed: %d (%s)\n", irun->pid, pid, irun->exe);
					running_remove(irun->pid, S_CRASH);
				} else{
					info_time("Launch job %d as pid %d (%s)\n", irun->pid, pid, irun->exe);
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
/**
   Replace newline characters by spaces.

   @param[in] path0 Input string.
   @return Newly allocated copy with newlines replaced; the caller frees it.
*/
char* remove_endl(const char* path0){
	char* path=strdup(path0);
	for(char* tmp=path; tmp[0]; tmp++){
		if(tmp[0]=='\n'){
			tmp[0]=' ';
		}
	}
	return path;
}
/**
   Add a new job to the running list with queued status and notify monitors.

   @param[in] exename Path to the executable.
   @param[in] execmd Job command line with arguments separated by newline.
*/
void queue_new_job(const char* exename, const char* execmd){
	RUN_T* irun=running_add(--counter, -1);
	if(!irun){
		warning_time("scheduler: running_add failed.\n");
		return;
	}
	irun->status.info=S_QUEUED;
	irun->exe=exename?strdup(exename):NULL;
	irun->path0=strdup(execmd);
	irun->path=remove_endl(irun->path0);
	info_time("%d (%s)\n", irun->pid, exename?exename:"");
	monitor_send(irun, irun->path);
	monitor_send(irun, NULL);
}
/**
   @return The first job in the list of queued or running jobs.
*/
RUN_T* running_head(void){
	return running;
}
/**
   @return The first job in the list of finished jobs.
*/
RUN_T* runned_head(void){
	return runned;
}
