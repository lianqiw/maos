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

/*
  All routines in this thread runs in a separate threads from the main thread.
  In order to call any gdk/gtk routines, it calls to go through
  g_idle_add. The thread is waits for select() results.

  All functions here except listen_host() are declared static to avoid calling
  from other threads to avoid corruptiong the data structures. Interaction to
  this thread is via sock_main. 

*/


#include <sys/socket.h>
#include "../sys/sys.h"
#include "monitor.h"
//the following are only by the thread running listen_host() to avoid race condition
static proc_t** pproc;
static int* nproc;
static int* hsock;   //socket of each host
static time_t* htime;//last time having signal from each host.
extern int sock_main[2]; /*listens command from main*/
static fd_set active_fd_set;
//int pipe_addhost[2]={0,0};
extern double* usage_cpu;
//extern double *usage_mem, *usage_mem2;
/*
static proc_t *proc_add(int id,int pid){
	proc_t *iproc;
	if((iproc=proc_get(id,pid))) return iproc;
	iproc=mycalloc(1,proc_t);
	iproc->iseed_old=-1;
	iproc->pid=pid;
	iproc->hid=id;
	iproc->next=pproc[id];
	pproc[id]=iproc;
	nproc[id]++;
	g_idle_add((GSourceFunc)update_title, GINT_TO_POINTER(id|nproc[id]<<8));
	return iproc;
}
*/
CHECK_ARG(1)
static void sendmail(const char*format,...){
	format2fn;
	if(mailto && fn){
		dbg_time("Sending mail to %s\n", mailto);
		char cmd[1024];
		snprintf(cmd, sizeof(cmd), "sendmail -t %s", mailto);
		FILE* p=popen(cmd, "w");
		if(p){
			fprintf(p, "To:%s\n%s\n", mailto, fn);
			fclose(p);
		}else{
			warning_time("sendmail failed for message (%s)\n", fn);
		}
	}

}
static proc_t* proc_get(int id, int pid){
	proc_t* iproc;
	if(id<0||id>=nhost){
		error("id=%d is invalid\n", id);
	}
	for(iproc=pproc[id]; iproc; iproc=iproc->next){
		if(iproc->pid==pid){
			break;
		}
	}
	if(!iproc){
		//info("%s: %d is not found\n", hosts[id], pid);
		iproc=mycalloc(1, proc_t);
		//iproc->iseed_old=-1;
		iproc->pid=pid;
		iproc->hid=id;
		iproc->next=pproc[id];
		pproc[id]=iproc;
		nproc[id]++;
		g_idle_add((GSourceFunc)update_title, GINT_TO_POINTER((id|nproc[id]<<8)));
	}
	return iproc;
}


static void proc_remove_all(int id){
	proc_t* iproc, * jproc=NULL;
	for(iproc=pproc[id]; iproc; iproc=jproc){
		jproc=iproc->next;
		g_idle_add((GSourceFunc)remove_entry, iproc->row);//frees iproc
		free(iproc->path);
		free(iproc);
	}
	nproc[id]=0;
	pproc[id]=NULL;
	g_idle_add((GSourceFunc)update_title, GINT_TO_POINTER((id|nproc[id]<<8)));
}

static void proc_remove(int id, int pid){
	proc_t* iproc, * jproc=NULL;
	for(iproc=pproc[id]; iproc; jproc=iproc, iproc=iproc->next){
		if(iproc->pid==pid){
			if(jproc){
				jproc->next=iproc->next;
			} else{
				pproc[id]=iproc->next;
			}
			nproc[id]--;

			g_idle_add((GSourceFunc)remove_entry, iproc->row);
			g_idle_add((GSourceFunc)update_title, GINT_TO_POINTER((id|nproc[id]<<8)));
			free(iproc->path);
			free(iproc);
			break;
		}
	}
}

static int host_from_sock(int sock){
	if(sock<0) return -1;
	for(int ihost=0; ihost<nhost; ihost++){
		if(hsock[ihost]==sock){
			return ihost;
		}
	}
	return -1;
}

/* Record the host after connection is established*/
static void host_added(int ihost, int sock){
	htime[ihost]=myclocki();
	proc_remove_all(ihost);/*remove all entries. */
	if(sock>-1){
		hsock[ihost]=sock;
		FD_SET(sock, &active_fd_set);
		g_idle_add(host_up, GINT_TO_POINTER(ihost));
	}
	dbg_time("connected to %s\n", hosts[ihost]);
}

/*remove the host upon disconnection*/
static void host_removed(int sock, int notify){
	if(sock>-1){
		shutdown(sock, SHUT_WR);
		close(sock);
		FD_CLR(sock, &active_fd_set);
	}
	int ihost=host_from_sock(sock);
	if(ihost!=-1 && hsock[ihost]!=-1){
		hsock[ihost]=-1;
		g_idle_add(host_down, GINT_TO_POINTER(ihost));
		warning_time("Disconnected from %s\n", hosts[ihost]);
		if(notify){
			sendmail("Subject:monitor on disconnected from %s\n\nAt %s\n",
				 hosts[ihost], myasctime(0));
		}
	}
}
//connect to scheduler(host)
static int add_host(int ihost){
	if(ihost<0 || ihost>=nhost){
		dbg_time("Invalid ihost=%d\n", ihost);
	}else if(hsock[ihost]>-1){
		dbg_time("host %d is already connected\n", ihost);
	}else if(hsock[ihost]==-1){
		hsock[ihost]--;//make it -2 so no concurrent access.
		int sock=scheduler_connect(hosts[ihost]);
		if(sock>-1){
			int cmd[2];
			cmd[0]=CMD_MONITOR;
			cmd[1]=1 | plot_enabled<<1 | 0x1<<31;
			if(stwriteintarr(sock, cmd, 2)){//write failed.
				warning_time("Failed to write to scheduler at %s\n", hosts[ihost]);
				close(sock);
				sock=-1;
			} else{
				host_added(ihost, sock);
			}
		}
		if(sock<0){
			dbg_time("Cannot reach %s\n", hosts[ihost]);
			hsock[ihost]=-1;
		}
	}else{
		dbg_time("add_host is already in progress\n");
	}
	return hsock[ihost];
}

static int test_jobs(int status, int flag){
	switch(flag){
	case -1://finished
		return status==S_FINISH;
		break;
	/*case -2://skipped
		return status==S_FINISH&&frac<1;
		break;*/
	case -3://crashed
		return status==S_CRASH||status==S_KILLED||status==S_TOKILL;
		break;
	case -4://all that is not running or pending
		return status==S_FINISH||status==S_CRASH||status==S_KILLED||status==S_TOKILL;
		break;
	default:
		return 0;
	}
}
static void save_all_jobs(){
	char* fnall=NULL;
	char* tm=strtime_pid();
	for(int ihost=0; ihost<nhost; ihost++){
		if(!pproc[ihost]) continue;
		char fn[PATH_MAX];
		const char* host=hosts[ihost];
		FILE* fp[2];
		snprintf(fn, PATH_MAX, "%s/maos_%s_%s.done", HOME, host, tm);
		fp[0]=fopen(fn, "w");
		snprintf(fn, PATH_MAX, "%s/maos_%s_%s.wait", HOME, host, tm);
		fp[1]=fopen(fn, "w");

		if(fnall){
			fnall=stradd(fnall, "\n", fn, NULL);
		} else{
			fnall=strdup(fn);
		}
		char* lastpath[2]={NULL,NULL};
		for(proc_t* iproc=pproc[ihost]; iproc; iproc=iproc->next){
			char* spath=iproc->path;
			char* pos=NULL;
			int id;
			pos=strstr(spath, "/maos ");
			if(!pos){
				pos=strstr(spath, "/skyc ");
			}
			if(iproc->status.info>10){
				id=0;
			} else{
				id=1;
			}
			if(pos){
				pos[0]='\0';
				if(!lastpath[id]||strcmp(lastpath[id], spath)){//a different folder.
					free(lastpath[id]); lastpath[id]=strdup(spath);
					fprintf(fp[id], "cd %s\n", spath);
				}
				pos[0]='/';
				fprintf(fp[id], "%s\n", pos+1);
			} else{
				fprintf(fp[id], "%s\n", spath);
			}
		}
		fclose(fp[0]);
		fclose(fp[1]);
		free(lastpath[0]);
		free(lastpath[1]);
	}
	info("Jobs saved to \n%s", fnall);
	free(fnall);
	free(tm);
}
static int scheduler_cmd(int ihost, int pid, int command);
//called by listen_host to respond to scheduler
static int respond(int sock){
	int cmd[3];
	//read fixed length header info.
	if(streadintarr(sock, cmd, 3)){
		return -1;//failed
	}
	int ihost=host_from_sock(sock);
	if(ihost>=0){
		htime[ihost]=myclocki();
	}
	int pid=cmd[2];
	switch(cmd[0]){
	case -1:{//server request shutdown
		dbg_time("Received shutdown request\n");
		return -3;
	}
		break;
	case MON_CMD://called by monitor main thread
	{
		int command;
		streadint(sock, &command);
		//dbg_time("received command %d for relaying to scheduler\n", command);
		int ans=scheduler_cmd(cmd[1], cmd[2], command);
		stwriteint(sock, ans);
	}
	break;
	case MON_DRAWDAEMON://called by scheduler to open drawdaemon
	{
		dbg_time("Received drawdaemon request\n");
		scheduler_cmd(ihost, 0, CMD_DISPLAY);
	}
	break;
	case MON_STATUS://called by scheduler to pass maos status
	{
		if(ihost<0){
			dbg_time("Host not found\n");
			return -1;
		}
		proc_t* iproc=proc_get(ihost, pid);
		iproc->tlast=myclocki();
		/*if(!iproc){
		iproc=proc_add(ihost,pid);
		}*/
		int old_info=iproc->status.info;
		if(stread(sock, &iproc->status, sizeof(status_t))){
			return -1;
		}
		if(!iproc->timstart){//only set once
			iproc->timstart=iproc->status.timstart?iproc->status.timstart:myclocki();
		}
		if(iproc->status.timlast){
			iproc->timlast=iproc->status.timlast;
		}
		if(iproc->status.info==S_REMOVE){
			proc_remove(ihost, pid);
		} else{
			if(cmd[1]!=ihost&&cmd[1]!=cmd[2]){
				/*pidnew is different from pid to indicate transition from
				queued to run or vice versa. There may be a race condition if
				new status from pid is sent before the transition.*/
				proc_t *p2=proc_get(ihost, cmd[1]);
				if(p2){//race condition happens, remove the newly created entry.
					if(p2->path&&!iproc->path){
						iproc->path=p2->path;
						p2->path=NULL;
					}
					proc_remove(ihost, cmd[1]);
				}
				iproc->pid=cmd[1];
			}
			//Alert when job crashed during running.
			if(old_info && old_info!=iproc->status.info 
				&& iproc->status.info==S_CRASH 
				&& (iproc->status.isim>0 || iproc->status.iseed>0)){
				sendmail("Subject: Job %d crashed on %s\n\nOn %s\n\nPath is %s\n", 
					iproc->pid, hosts[iproc->hid], myasctime(0), iproc->path);
			}
			g_idle_add((GSourceFunc)refresh, iproc);
		}
	}
	break;
	case MON_PATH://called by scheduler to pass maos path
	{
		if(ihost<0){
			dbg_time("Host not found\n");
			return -1;
		}
		proc_t* iproc=proc_get(ihost, pid);
		/*if(!iproc){
		iproc=proc_add(ihost,pid);
		}*/
		if(streadstr(sock, &iproc->path)){
			return -1;
		}
		char* tmp=NULL;
		while((tmp=strchr(iproc->path, '\n'))){
			tmp[0]=' ';
		}
	}
	break;
	case MON_CLEARJOB://called by monitor to clear all jobs unless pid is specified
	{
		ihost=cmd[1];
		int flag=cmd[2];
		for(proc_t* iproc=pproc[ihost]; iproc; iproc=iproc->next){
			if((flag < 0 && test_jobs(iproc->status.info, flag)) || iproc->pid==flag){
				if(scheduler_cmd(ihost, iproc->pid, CMD_REMOVE)){
					warning_time("Failed to clear the job %d on host %d\n", iproc->pid, ihost);
				}
			}
		}
	}
	break;
	case MON_KILLJOB://called by monitor to kill all jobs unless pid is specified
	{
		ihost=cmd[1];
		for(proc_t* iproc=pproc[ihost]; iproc; iproc=iproc->next){
			if(iproc->status.info<11&&(!pid||pid==iproc->pid)){
				if(scheduler_cmd(ihost, iproc->pid, CMD_KILL)){
					warning_time("Failed to kill the job %d on host %d\n", iproc->pid, ihost);
				}
			}
		}
	}
	break;
	case MON_SAVEJOB: //called by monitor to save all jobs
	{
		save_all_jobs();
	}
	break;
	case MON_VERSION://send by scheduler as a dummy message
		break;
	case MON_LOAD:
	{
		if(ihost<0){
			dbg_time("Host not found\n");
			return -1;
		}
		usage_cpu[ihost]=(double)((pid>>16)&0xFFFF)/100.;
		//usage_mem[ihost]=(double)(pid & 0xFFFF)/100.;
		usage_cpu[ihost]=MAX(MIN(1, usage_cpu[ihost]), 0);
		//usage_mem[ihost]=MAX(MIN(1,usage_mem[ihost]),0);
		g_idle_add((GSourceFunc)update_progress, GINT_TO_POINTER(ihost));
	}
	break;
	case MON_ADDHOST:
		if(cmd[1]>-1&&cmd[1]<nhost){
			add_host(cmd[1]);
		} else if(cmd[1]==-2){//quit
			return -2;
		}
		break;
	default:
		warning_time("Invalid cmd %d\n", cmd[0]);
		return -1;
	}
	return 0;
}
/**
   listen_host() live in a dedicated thread, it has the following resposibilities:
   1) listening commands from the main thread to initiate connection to servers
   2) listening to connected servers for maos status event and update the display
   3) monitor connected servers for activity. Disable pages when server is disconnected.

   write to sock_main[1] will be caught by select in listen_host(). This wakes it up.*/
void* listen_host(void* pmsock){
	int msock=GPOINTER_TO_INT(pmsock);
	pproc=mycalloc(nhost+1, proc_t*);
	nproc=mycalloc(nhost+1, int);
	hsock=mycalloc(nhost+1, int);
	for(int i=0; i<=nhost; i++){
		hsock[i]=-1;
	}
	htime=mycalloc(nhost, time_t);
	FD_ZERO(&active_fd_set);
	FD_SET(msock, &active_fd_set);//listen to monitor itself
	int keep_listen=1;
	while(keep_listen){
		fd_set read_fd_set=active_fd_set;
		struct timeval timeout={5,0};
		if(select(FD_SETSIZE, &read_fd_set, NULL, NULL, &timeout)<0){
			perror("select");
			continue;
		}
		for(int i=0; i<FD_SETSIZE; i++){
			if(FD_ISSET(i, &read_fd_set)){
				int res;
				if((res=respond(i))<0){
					host_removed(i, res==-1?1:0);
				}
				if(res==-2){//quit
					keep_listen=0;
					break;
				}
			}
		}
		time_t ntime=myclocki();
		for(int ihost=0; ihost<nhost; ihost++){
			if(htime[ihost]){//only handle hosts that are ever connected
				if(hsock[ihost]<0){//disconnected, trying to reconnect
					if(ntime>htime[ihost]+60){//try every 600 seconds
						htime[ihost]=ntime;
						add_host(ihost);//do not use _add_host_wrap. It will deadlock.
					}
				} else if(htime[ihost]>0&&ntime>htime[ihost]+60){//no activity for 60 seconds. check host connection 
					dbg_time("60 seconds no respond. probing server %s.\n", hosts[ihost]);
					scheduler_cmd(ihost, 0, CMD_PROBE);
					htime[ihost]=-ntime;
				} else if(htime[ihost]<0&&ntime>-htime[ihost]+60){//probed, but not response within 60 seconds
					dbg_time("no respond. disconnect server %s.\n", hosts[ihost]);
					host_removed(hsock[ihost], 1);
				}
			}
			//check for jobs that may have hung
			for(proc_t* iproc=pproc[ihost]; iproc; iproc=iproc->next){
				if(iproc->status.info==S_RUNNING&&iproc->tlast+600<ntime){
					iproc->status.info=S_CRASH;
					iproc->status.tot=(ntime-iproc->tlast)/iproc->status.scale;
					dbg_time("proc %d in %s is not updating in %lu seconds.\n", iproc->pid, hosts[ihost], ntime-iproc->tlast);
					g_idle_add((GSourceFunc)refresh, iproc);
				}
			}
		}
	}
	for(int i=0; i<FD_SETSIZE; i++){
		if(FD_ISSET(i, &active_fd_set)){
			close(i);
			FD_CLR(i, &active_fd_set);
		}
	}
	return NULL;
}
/**
   called by monitor to let a MAOS job remotely draw on this computer
*/
static int scheduler_display(int ihost, int pid){
	int ans=1;
	if(ihost<0||ihost>=nhost) return ans;
	dbg_time("scheduler_display: %s, pid=%d\n", hosts[ihost], pid);
	/*connect to scheduler with a new port. The schedule pass the other
	  end of the port to drawdaemon so we can communicate with it.*/
	int sock=scheduler_connect(hosts[ihost]);
	if(sock==-1) return ans;
	int cmd[2]={CMD_DISPLAY, pid};
	
	if(stwriteintarr(sock, cmd, 2)||streadintarr(sock, cmd, 1)||cmd[0]){
		warning("Failed to pass sock to draw via scheduler.\n");
	} else{
		char arg1[20];
		snprintf(arg1, 20, "%d", sock);
		if(spawn_process("drawdaemon", arg1, NULL)<0){
			warning("spawn drawdaemon failed\n");
		} else{
			ans=0;
		}
	}
	close(sock);
	return ans;
}
/**
   called by monitor to talk to scheduler.
*/
static int scheduler_cmd(int ihost, int pid, int command){
	if(ihost<0) return 0;
	if(command==CMD_DISPLAY){
		return scheduler_display(ihost, pid);
	} else{
		int ans=-1;
		int sock=hsock[ihost];
		if(sock>-1){
			int cmd[2];
			cmd[0]=command;
			cmd[1]=pid;/*pid */
			ans=stwriteintarr(sock, cmd, 2);
			if(ans){/*communicated failed.*/
				host_removed(sock, 0);
			}
		}
		return ans;
	}
}
