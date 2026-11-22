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
#include <ctype.h>
#include "../sys/sys.h"
#include "scheduler.h"
#include "scheduler_job.h"

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


/**head of linked list*/
static MONITOR_T* pmonitor=NULL;
static double usage_cpu;
//static MONITOR_T *monitor_get(int hostid);
static void scheduler_timeout(void);
static void monitor_send_initial(MONITOR_T* ic);
static void monitor_send_load(void);
/**every 3 seconds*/
static time_t lasttime3=0;
/**every 10 seconds*/
static time_t lasttime10=0;
#define MAXSOCK 1024
/**saving the client name of each socket*/
static const char* sockname[MAXSOCK]={0};
/**whether a TCP socket is authenticated with CMD_AUTH*/
static char sockauth[MAXSOCK]={0};

/**
   Send a command to a running maos process via its command socket, optionally
   passing a socket to it. The reply to sock (0 on success, -1 on failure) is
   written before the socket is passed. Logs a warning on failure.

   @param[in] pid Pid of the target maos process.
   @param[in] sock Socket to reply to and pass to maos, or -1 if none.
   @param[in] cmd Command to send to maos (e.g. MAOS_DRAW, MAOS_VAR).
   @return 0 on success, -1 if the process or its command socket is not found or writing the command fails.
*/
int maos_command(int pid, int sock, int cmd){
	RUN_T* irun=running_get_by_pid(pid);
	int cmd2[2]={cmd, 0};
	int failed=(!irun||irun->sock_cmd<0||(stwriteintarr(irun->sock_cmd, cmd2, 2)));
	if(sock!=-1 && !stwriteint(sock, failed?-1:0) && !failed){//scheduler reply first before maos gets the chance to start writing to the socket
		stwritefd(irun->sock_cmd, sock);
	}
	if(failed){
		warning_time("Failed to send command %d and socket %d to maos (PID=%d, sock_cmd=%d)\n", cmd, sock, pid, irun?irun->sock_cmd:-1);
	} else{
		dbg_time("Successfully send command %d and socket %d to maos (PID=%d, sock_cmd=%d)\n", cmd, sock, pid, irun->sock_cmd);
	}
	return failed?-1:0;//do not keep this connection.
}
/**
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
   Save a drawdaemon socket with an id in the linked list, or update the id if
   the socket is already saved.

   @param[in] sock_save Socket descriptor to save. Ownership moves to the list.
   @param[in] id Id (session id) to associate with the socket.
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
/**
   Unlink p from the saved socket list and free it. Does not close the socket.

   @param[in] p Entry to remove. Freed.
*/
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
/**
   Retrieve and remove a saved drawdaemon socket by id. Bad sockets encountered
   are closed and removed. If none matches id, a socket saved with id 0 is used.

   @param[in] id Id to look up.
   @return The socket, now owned by the caller, or -1 if not found.
*/
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
/**
   Send a heartbeat to all saved drawdaemon sockets and remove those that fail.
*/
static void socket_heartbeat(){
	for(SOCKID_T* p_next, *p=shead; p; p=p_next){
		p_next=p->next;
		int cmd[1]={DRAW_HEARTBEAT};
		if(stcheck(p->sock) || stwrite(p->sock, cmd, sizeof(cmd))){
			dbg_time("Remove bad socket %d with id %d\n", p->sock, p->id);
			close(p->sock);
			socket_remove(p);
		}
	}
}
/**
   Close all saved drawdaemon sockets upon scheduler exit.
*/
static void socket_close(){
	for(SOCKID_T *p_next, *p=shead; p; p=p_next){
		p_next=p->next;
		close(p->sock);
		socket_remove(p);
	}
}
/**
   Name a socket for logging. If the socket was previously a monitor, the stale
   monitor is removed.

   @param[in] sock Socket descriptor.
   @param[in] name Name of the client. Must remain valid (string literal).
*/
static void set_sockname(int sock, const char *name){
	if(sock>-1&&sock<MAXSOCK&&(!sockname[sock]||strcmp(sockname[sock],name))){
		if(sockname[sock] && !strcmp(sockname[sock], "monitor")){
			warning_time("remove staled monitor socket %d\n", sock);
			monitor_remove(sock);//sock is reused with old monitor sock
		}
		dbg_time("sock %d is named %s\n", sock, name);
		sockname[sock]=name;
	}
}
/**
   Clear the name of a socket.

   @param[in] sock Socket descriptor.
*/
static void unset_sockname(int sock){
	if(sock>-1&&sock<MAXSOCK){
		sockname[sock]=NULL;
	}
}
/**
   Get the name of a socket for logging.

   @param[in] sock Socket descriptor.
   @return The name, or "" if not named. Do not free.
*/
static const char* get_sockname(int sock){
	if(sock>-1 && sock<MAXSOCK && sockname[sock]){
		return sockname[sock];
	}else{
		return "";
	}
}

/**>-1: there is pending scheduler_recv_socket.*/
static int scheduler_recv_wait=-1;
/**
   Ask the first monitor with plotting enabled to start a drawdaemon for pid,
   and record sock as the pending draw request.

   @param[in] sock Socket of draw() waiting for a drawdaemon.
   @param[in] pid Pid of the process requesting drawing.
   @return 0 on success, nonzero if sending failed, -1 if no monitor supports plotting.
*/
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
					((int*)buf)[2]=-1;//-1 indicate text data
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
   Hand a drawdaemon socket to draw(). With pid<=0 it is passed to the pending
   draw() request from monitor_request_draw(); otherwise it is passed to maos pid
   via maos_command().

   @param[in] sock Drawdaemon socket.
   @param[in] pid <=0: service the pending request. >0: ask maos pid to draw with sock.
   @return 0 for success, -1 for failure.
*/
int send_draw_sock(int sock, int pid){
	int ret;//return status
	set_sockname(sock, "drawdaemon");
	if(pid<=0){//this is for pending scheduler_recv_socket or a new idle connection
		if(scheduler_recv_wait==-1||stwriteint(scheduler_recv_wait, 0)||stwritefd(scheduler_recv_wait, sock)){
			if(scheduler_recv_wait!=-1){
				dbg_time("(%d:%s) failed to pass sock to draw at %d.\n", sock, get_sockname(sock), scheduler_recv_wait);
			}else{
				dbg_time("(%d:%s) no pending drawdaemon request.\n", sock, get_sockname(sock));
			}
			ret=-1;
		} else{
			dbg_time("passed sock %d to draw at %d\n", sock, scheduler_recv_wait);
			stwriteint(sock, 0);
			ret=0;
		}
		scheduler_recv_wait=-1;
	} else{//ask maos to start drawing with this drawdaemon
		RUN_T* irun=running_get_by_pid(pid);
		if(!irun){
			dbg_time("(%d:%s) maos PID=%d not longer exists.\n", sock, get_sockname(sock), pid);
			ret=-1;
		}else if(maos_command(pid, sock, MAOS_DRAW)){
			dbg_time("(%d:%s) failed to pass sock to maos PID=%d.\n", sock, get_sockname(sock), pid);
			ret=-1;
		}else{
			dbg_time("(%d:%s) passed sock to maos PID=%d.\n", sock, get_sockname(sock), pid);
			ret=0;
		}
	}
	dbg_time("returns %d for sock=%d, pid=%d\n", ret, sock, pid);
	return ret;
}
/**
   Check that a backtrace request is an addr2line command without shell special
   characters, as it is run by popen().

   @param[in] cmd The command.
   @return 1 if valid, 0 otherwise.
*/
static int trace_cmd_valid(const char *cmd){
	if(!cmd || strncmp(cmd, "addr2line ", 10)) return 0;
	for(const char *p=cmd; *p; p++){
		if(!isalnum((int)*p) && !strchr(" /._+-", *p)) return 0;
	}
	return 1;
}
/**
   Check whether a socket is a unix (local) socket.

   @param[in] sock Socket descriptor.
   @return 1 if it is a unix socket, 0 otherwise.
*/
static int sock_is_unix(int sock){
	struct sockaddr_storage addr;
	socklen_t len=sizeof(addr);
	return !getsockname(sock, (struct sockaddr*)&addr, &len) && addr.ss_family==AF_UNIX;
}
/**
   Check whether a client may send commands: unix sockets are local and are
   protected by file permissions. TCP sockets must authenticate with CMD_AUTH.

   @param[in] sock Socket descriptor.
   @return 1 if allowed, 0 otherwise.
*/
static int sock_authed(int sock){
	if(sock>=0 && sock<MAXSOCK && sockauth[sock]) return 1;
	return sock_is_unix(sock);
}
/**
   Clear references to a socket that is being closed, so that they do not
   point to another connection that later reuses the descriptor.

   @param[in] sock Socket descriptor.
*/
static void sock_closed(int sock){
	for(RUN_T* irun=running_head(); irun; irun=irun->next){
		if(irun->sock==sock) irun->sock=-1;
		if(irun->sock_cmd==sock) irun->sock_cmd=-1;
	}
	if(scheduler_recv_wait==sock) scheduler_recv_wait=-1;
	if(sock>=0 && sock<MAXSOCK) sockauth[sock]=0;
	for(MONITOR_T* ic=pmonitor; ic; ic=ic->next){
		if(ic->sock==sock){
			monitor_remove(sock);
			break;
		}
	}
	unset_sockname(sock);
}
/*
   The following respond_* functions handle one command each. They are called
   by respond() after the fixed header int[2] (command and pid) is read and
   read the optional payload of the command. The meaning of pid depends on the
   command. See the protocol description in docs/94_scheduler.md.

   Each returns 0 to keep the connection or -1 to close it.
*/
/**CMD_START (1): called by maos when job starts. payload: int nthread, int waiting (bit 0: wait for the scheduler to start the job. bit 1 and up: number of GPUs).*/
static int respond_start(int sock, int pid){
	set_sockname(sock, "maos");
	int nthread;
	int ngpu;
	int waiting;
	if(streadint(sock, &nthread)||streadint(sock, &waiting)){
		return -1;
	}
	ngpu=waiting>>1;
	if(ngpu>NGPU){
		ngpu=NGPU;
	}else if(ngpu<0){
		ngpu=0;
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
		return 0;
	}
	nrun_sub(irun);//in case the job reports start again
	irun->nthread=nthread;
	irun->ngpu=ngpu;
	if(waiting&0x1){
		dbg_time("(%d:%s) Job %d waiting to start with nused_cpu=%d, ngpu=%d\n", sock, get_sockname(sock), pid, nthread, ngpu);
		irun->status.info=S_WAIT;
	} else{/*no waiting, no need reply. */
		dbg_time("(%d:%s) Job %d started with nused_cpu=%d, ngpu=%d\n", sock, get_sockname(sock), pid, nthread, ngpu);
		irun->status.info=S_START;
		nrun_add(irun);
	}
	irun->last_time=myclocki();
	if(irun->path) monitor_send(irun, irun->path);
	monitor_send(irun, NULL);
	return 0;
}
/**CMD_FINISH (2), CMD_CRASH (4), CMD_KILLED (16): called by maos when job finishes, crashes, or is cancelled.*/
static int respond_end(int sock, int pid, int status){
	dbg_time("(%d:%s) Job %d reports %s\n", sock, get_sockname(sock), pid, status==S_FINISH?"finished":(status==S_CRASH?"crashed":"killed"));
	running_remove(pid, status);
	process_queue();//process the next job. don't wait for timeout
	return 0;//Do not yet close. wait for client to close.
}
static int respond_finish(int sock, int pid){
	return respond_end(sock, pid, S_FINISH);
}
static int respond_crash(int sock, int pid){
	return respond_end(sock, pid, S_CRASH);
}
static int respond_killed(int sock, int pid){
	return respond_end(sock, pid, S_KILLED);
}
/**CMD_STATUS (3): called by maos to report status at every time step. payload: status_t.*/
static int respond_status(int sock, int pid){
	//Read the status first so that it is available for an unrecorded job and the
	//socket stays in sync even if the job cannot be added.
	status_t status;
	if(stread(sock, &status, sizeof(status_t))){
		warning_time("Error reading status.\n");
		return -1;
	}
	RUN_T* irun=running_get_by_pid(pid);
	if(!irun){/*started before scheduler is relaunched. */
		dbg_time("(%d:%s) pid=%d is running but not recorded.\n", sock, get_sockname(sock), pid);
		irun=running_add(pid, sock);
		if(!irun){
			dbg_time("scheduler: running_add %d failed. Exe already exited.\n", pid);
			return 0;
		}
		irun->status.info=S_START;
		irun->nthread=MAX(1, MIN(status.nthread, NCPU));
		nrun_add(irun);
		monitor_send(irun, NULL);
	}
	//save and restore time information
	time_t timlast=irun->status.timlast;
	time_t timstart=irun->status.timstart;
	irun->status=status;
	irun->status.timlast=timlast;
	irun->status.timstart=timstart;
	irun->last_time=myclocki();
	monitor_send(irun, NULL);
	return 0;
}
/**CMD_MONITOR (5): called by monitor when it connects. pid is the option flag of monitor_add().*/
static int respond_monitor(int sock, int pid){
	monitor_add(sock, pid, NULL, NULL);
	dbg_time("(%d:%s) monitor from %s is connected.\n", sock, get_sockname(sock), addr2name(socket_peer(sock)));
	return 0;
}
/**CMD_PATH (6): called by maos to report the job path. payload: string.*/
static int respond_path(int sock, int pid){
	set_sockname(sock, "maos");
	RUN_T* irun=running_add(pid, sock);
	if(!irun){
		warning_time("(%d:%s) running_add %d failed. Exe already exited.\n", sock, get_sockname(sock), pid);
		return 0;
	}
	free(irun->path0); irun->path0=NULL;
	free(irun->path); irun->path=NULL;
	if(streadstr(sock, &irun->path0)||!irun->path0){
		dbg_time("(%d:%s) Job %d receiving path failed. \n", sock, get_sockname(sock), pid);
		return -1;
	}
	irun->path=remove_endl(irun->path0);
	irun->last_time=myclocki();
	dbg_time("(%d:%s) Job %d received path. \n", sock, get_sockname(sock), pid);
	monitor_send(irun, irun->path);
	return 0;
}
/**CMD_KILL (7): called by monitor to kill a job.*/
static int respond_kill(int sock, int pid){
	warning_time("(%d:%s) Received monitor command to kill %d\n", sock, get_sockname(sock), pid);
	running_kill(pid);
	return 0;
}
/**CMD_TRACE (8): called by maos to request a backtrace. payload: string addr2line command. reply: string.*/
static int respond_trace(int sock, int pid){
	(void)pid;
	set_sockname(sock, "trace");
	int ret=0;
	char* buf=NULL;
	char out[200]={0};
	if(streadstr(sock, &buf)
		||!trace_cmd_valid(buf)//buf is run by the shell
		||call_addr2line(out, sizeof out, buf)
		||stwritestr(sock, out)){
		warning_time("(%d:%s) CMD_TRACE failed. buf=%s, out=%s\n", sock, get_sockname(sock), buf?buf:"", out);
		ret=-1;
	}
	free(buf);
	return ret;
}
/**CMD_PROBE (9): called by monitor to probe the connection. reply: int[3] {MON_VERSION, 0, 0}.*/
static int respond_probe(int sock, int pid){
	(void)pid;
	int cmd2[3]={MON_VERSION, 0, 0};
	stwrite(sock, &cmd2, sizeof(cmd2));
	return 0;
}
/**
   CMD_DRAWCLI (10): called by draw() to cache or request a drawdaemon socket.
   pid>0: save the socket that follows with id pid. -10000<pid<0: reply int 0
   and the saved socket with id -pid, or -1. Otherwise, request a new
   drawdaemon from the monitor for pid+10000.
*/
static int respond_drawcli(int sock, int pid){
	int ret=0;
	if(pid>0){//receive sock from draw() and save for later use
		set_sockname(sock, "save drawdaemon");
		int sock_save;
		if(streadfd(sock, &sock_save)){
			warning_time("(%d:%s) receive socket for saving failed\n", sock, get_sockname(sock));
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
	return ret;
}
/**CMD_REMOVE (11): called by monitor to remove a finished job from the list.*/
static int respond_remove(int sock, int pid){
	(void)sock;
	runned_remove(pid);
	return 0;
}
/**CMD_DRAWSER (12): called by drawdaemon or monitor to provide a drawdaemon to draw(). pid>0: for maos pid. Otherwise for a pending request with id -pid.*/
static int respond_drawser(int sock, int pid){
	int sock2=dup(sock);
	dbg_time("(%d:%s) drawdaemon for pid=%d from %s. duplicated as %d\n", sock, get_sockname(sock), pid, addr2name(socket_peer(sock)), sock2);
	if(sock2<0){
		warning_time("(%d:%s) dup failed: %s\n", sock, get_sockname(sock), strerror(errno));
	}else if(send_draw_sock(sock2, pid)){//pass a duplicated fd to maos
		dbg_time("Save drawdaemon %d for future use\n", sock2);
		socket_save(sock2, pid==INT_MIN?0:abs(pid));
		stwriteint(sock2, 0);//respond OK to drawdaemon
	}else{
		close(sock2);//the receiver has its own copy now.
	}
	return -1;//remove this sock from scheduler receiving pollfd
}
/**CMD_MAOSCLI (13): for a maos client to create a client link to maos pid.*/
static int respond_maoscli(int sock, int pid){
	set_sockname(sock, "maos client");
	int sock2=dup(sock);//pass a duplicated fd to maos
	info_time("(%d:%s) pass command to maos %d for client at %s\n", sock, get_sockname(sock), pid, addr2name(socket_peer(sock)));
	if(sock2>=0){
		maos_command(pid, sock2, MAOS_VAR);
		close(sock2);//maos has its own copy if passing succeeded.
	}
	return -1;//remove this sock from scheduler receiving pollfd
}
/**CMD_MAOSSER (14): called by maos to save a connection to run maos_command().*/
static int respond_maosser(int sock, int pid){
	set_sockname(sock, "maos");
	RUN_T* irun=running_get_by_pid(pid);
	if(irun){
		irun->sock_cmd=sock;
		dbg_time("(%d:%s) maos server socket is saved\n", sock, get_sockname(sock));
		return 0;
	} else{
		warning_time("(%d:%s) maos server socket does not have irun. close.\n", sock, get_sockname(sock));
		return -1;
	}
}
/**CMD_RESTART (15): called by monitor to restart a job.*/
static int respond_restart(int sock, int pid){
	info_time("(%d:%s) restart job %d\n", sock, get_sockname(sock), pid);
	runned_restart(pid);
	return 0;
}
/**CMD_DUMMY (17): probe. no action and no reply.*/
static int respond_dummy(int sock, int pid){
	(void)sock; (void)pid;
	return 0;
}
/**
   CMD_LAUNCH (18): called by maos from another machine to start a job in this
   machine. pid is the number of strings that follow: 2: exename and execmd,
   3: exename, execwd (old method), and execmd. 1: execmd. reply: int 0 or -1.
*/
static int respond_launch(int sock, int pid){
	set_sockname(sock, "job launcher");
	int ret=0;
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
	if(ret!=-1&&!streadstr(sock, &execmd)&&execmd){
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
	return ret;
}
/**CMD_AUTH (19): called by a TCP client to authenticate. See scheduler_auth(). reply: string nonce. payload: string base64_sha1(nonce+password). reply: int 0 or -1.*/
static int respond_auth(int sock, int pid){
	(void)pid;
	char *pass=scheduler_password(1);
	char *nonce=generate_session_cookie();
	char *key=NULL;
	int ret=-1;
	if(pass && nonce && sock<MAXSOCK && !stwritestr(sock, nonce) && !streadstr(sock, &key) && key){
		char *src=stradd(nonce, pass, NULL);
		char key2[64];
		base64_sha1(src, key2);
		free(src);
		if(!strcmp(key, key2)){
			sockauth[sock]=1;
			ret=0;
		}
	}
	if(ret){
		warning_time("(%d:%s) authentication from %s failed.\n", sock, get_sockname(sock), addr2name(socket_peer(sock)));
	}
	if(stwriteint(sock, ret)){
		ret=-1;
	}
	free(pass);
	free(nonce);
	free(key);
	return ret;
}
#define RESPOND_PID 1  /*pid must be positive: a command from maos about its own job*/
#define RESPOND_UNIX 2 /*only valid over a unix socket (passes file descriptors)*/
#define RESPOND_AUTH 4 /*allowed before the TCP connection is authenticated*/
/**
   The handler and requirement of each command.
*/
static const struct{
	int cmd;
	int (*fun)(int sock, int pid);
	int flag;
}respond_list[]={
	{CMD_START,   respond_start,   RESPOND_PID},
	{CMD_FINISH,  respond_finish,  RESPOND_PID},
	{CMD_STATUS,  respond_status,  RESPOND_PID},
	{CMD_CRASH,   respond_crash,   RESPOND_PID},
	{CMD_MONITOR, respond_monitor, 0},
	{CMD_PATH,    respond_path,    RESPOND_PID},
	{CMD_KILL,    respond_kill,    0},
	{CMD_TRACE,   respond_trace,   RESPOND_PID},
	{CMD_PROBE,   respond_probe,   0},
	{CMD_DRAWCLI, respond_drawcli, RESPOND_UNIX},
	{CMD_REMOVE,  respond_remove,  0},
	{CMD_DRAWSER, respond_drawser, 0},
	{CMD_MAOSCLI, respond_maoscli, 0},
	{CMD_MAOSSER, respond_maosser, RESPOND_PID},
	{CMD_RESTART, respond_restart, 0},
	{CMD_KILLED,  respond_killed,  RESPOND_PID},
	{CMD_DUMMY,   respond_dummy,   0},
	{CMD_LAUNCH,  respond_launch,  0},
	{CMD_AUTH,    respond_auth,    RESPOND_AUTH},
};
//PNEW(mutex_sch);//respond() and scheduler_handle_ws() much lock this before preceed.
/**
   Respond to client requests. The fixed header is int[2] for opcode and
   pid. Optional data follows based on opcode and is read by the handler in
   respond_list. The disadvantage of this is that new commands with additional
   payload cannot be easily introduced in the scheduler().

   @param[in] pfd Poll entry of the client socket.
   @param[in] flag -1 asks the client to shut down before reading. Otherwise unused.
   @return 0 to keep the connection, -1 to close the socket.
*/
static int respond(struct pollfd *pfd, int flag){
	if(!pfd){
		return 0;
	}
	if(flag==-1){//ask client to shutdown
		int cmd[3]={-1, 0, 0};
		stwrite(pfd->fd, cmd, 3*sizeof(int));
		return 0;//do not block reading. the reply is handled when it arrives.
	}
	int sock=pfd->fd;
	/*
	   Don't modify sock in this routine. otherwise select
	   will complain Bad file descriptor
	*/
	int ret=-1, cmd[2]={0,0};
	if(!streadintarr(sock, cmd, 2)){//otherwise connection lost or closed by client.
		dbg2_time("(%d:%s) cmd is %d %d. \n", sock, get_sockname(sock), cmd[0], cmd[1]);
		int ic=0, nc=sizeof(respond_list)/sizeof(respond_list[0]);
		while(ic<nc && respond_list[ic].cmd!=cmd[0]) ic++;
		if(ic==nc){
			warning_time("(%d:%s) Invalid cmd: %d\n", sock, get_sockname(sock), cmd[0]);
		}else if(!(respond_list[ic].flag&RESPOND_AUTH) && !sock_authed(sock)){
			warning_time("(%d:%s) cmd %d from %s is rejected as the connection is not authenticated.\n",
				sock, get_sockname(sock), cmd[0], addr2name(socket_peer(sock)));
		}else if((respond_list[ic].flag&RESPOND_PID) && cmd[1]<=0){
			warning_time("(%d:%s) cmd %d with invalid pid %d is rejected.\n", sock, get_sockname(sock), cmd[0], cmd[1]);
		}else if((respond_list[ic].flag&RESPOND_UNIX) && !sock_is_unix(sock)){
			warning_time("(%d:%s) cmd %d is rejected as it requires a unix socket.\n", sock, get_sockname(sock), cmd[0]);
		}else if((ret=respond_list[ic].fun(sock, cmd[1]))){
			dbg_time("cmd=%d returned ret=%d\n", cmd[0], ret);
		}
	}
	if(ret){
		RUN_T* irun=running_get_by_sock(sock);//is maos
		if(irun&&irun->status.info<10){
			check_job_exists(irun);
		}
		sock_closed(sock);
		ret=-1;
	}
	return ret;//ret=-1 will close the socket.
}


/**
   Periodic task called by listen_port(): process the job queue, check jobs every
   3 seconds, and report load and send drawdaemon heartbeats every 10 seconds.
*/
static void scheduler_timeout(void){
	time_t thistime=myclocki();
	//Process job queue
	process_queue();
	//respond to heatbeat of saved sockets

	if(thistime>=(lasttime3+3)){
		if(running_head()){
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
/**
   Add a monitor client to the MONITOR_T linked list (or update an existing one
   for the same socket) and send it the full job list.

   @param[in] sock Socket of the monitor.
   @param[in] flag Options. >=0x8: old scheme, load reporting only. <0: bit 0 enables load, bit 1 plotting, bit 2 text (http) mode.
   @param[in] func Function to send text directly to the client, or NULL to write to sock.
   @param[in] userdata Parameter passed to func.
*/
void monitor_add(int sock, int flag, int (*func)(char* buf, int nlen, int mode, void *userdata), void* userdata){
	//dbg_time("added monitor on sock %d\n",sock);
	set_sockname(sock, "monitor");
	socket_send_timeout(sock, 5);
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
/**
   Remove the monitor connected at sock from the monitor list and shut down the
   socket for writing.

   @param[in] sock Socket of the monitor.
*/
void monitor_remove(int sock){
	MONITOR_T* ic=NULL;
	for(MONITOR_T **curr=&pmonitor; *curr;){
		ic=*curr;
		if(ic->sock==sock){
			*curr=ic->next;
			dbg_time("Removed monitor at %d\n", ic->sock);
			unset_sockname(ic->sock);
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
/**
   Format job path or status as text for http (websocket) monitors.

   @param[in] irun Job record.
   @param[in] path Job path. If NULL, the status is formatted instead.
   @param[out] buf Output buffer.
   @param[in] size Size of buf.
   @return Length of the text (as snprintf, may exceed size if truncated).
*/
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
/**
   Send job path or status to one monitor.

   @param[in] irun Job record.
   @param[in] path Job path. If NULL, the status is sent instead.
   @param[in] ic Monitor to send to.
   @return 0 on success, nonzero on failure.
*/
static int monitor_send_do(RUN_T* irun, const char* path, MONITOR_T *ic){
	if(!irun || !ic) return -1;
	int ans=0;
	if(ic->http){
		char buf[4096];
		int padding=ic->func?http_padding:(int)(3*sizeof(int));
		int nlen=http_convert(irun, path, buf+padding, sizeof(buf)-padding);
		if(nlen>=(int)sizeof(buf)-padding) nlen=sizeof(buf)-padding-1;//truncated
		if(ic->func){//directly write to client
			ans=ic->func(buf+padding, nlen, 0, ic->userdata);
		}else{//send to proxy
			((int*)buf)[0]=DRAW_ENTRY;
			((int*)buf)[1]=nlen;
			((int*)buf)[2]=-1;//-1 indicate text data. >=0: binary data for draw
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

/**
   Notify already connected monitors of a job update. Monitors that fail are
   removed.

   @param[in] irun Job record.
   @param[in] path Job path. If NULL, the status is sent instead.
*/
void monitor_send(RUN_T* irun, const char* path){
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
/**
   Query nvidia-smi for GPU usage. GPU usage query is disabled permanently after
   the first failure.

   @param[out] gpu_mem Fraction of total GPU memory used. 0 if disabled.
   @return Fraction of GPUs that use more than 500 MB of memory, or 0 if disabled.
*/
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
			if(!total0||!ngpu){
				enabled=0;
				dbg_time("total0=0, disable gpu usage.\n");
				return 0;
			}
			*gpu_mem=(double)used0/total0;
			return (double)nuse/ngpu;
		}
	}
	return 0;
}
/**
   Notify already connected monitors of the machine load.
*/
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
/**
   Send the newly added monitor information of all jobs. The monitor is removed
   on failure.

   @param[in] ic Monitor to send to.
*/
static void monitor_send_initial(MONITOR_T* ic){
	int sock=ic->sock;
	for(RUN_T* irun=runned_head(); irun; irun=irun->next){
		if(monitor_send_do(irun, irun->path, ic)||monitor_send_do(irun, NULL, ic)){
			monitor_remove(sock);
			return;
		}
	}
	for(RUN_T* irun=running_head(); irun; irun=irun->next){
		if(monitor_send_do(irun, irun->path, ic)||monitor_send_do(irun, NULL, ic)){
			monitor_remove(sock);
			return;
		}
	}
}

/**
   Scheduler: queues and launches jobs according to available CPUs and GPUs,
   keeps their status, and serves monitor, drawdaemon, and browser (http and
   websocket) clients. -d runs as daemon; --direct skips locking as it is already
   done by the launcher.

   @param[in] argc Number of arguments.
   @param[in] argv Arguments: [-h|--help] [-d] [--direct].
   @return 0 for help, 1 if the lock file cannot be acquired. Otherwise exits with 0.
*/
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
	signal(SIGPIPE, SIG_IGN);//a client that disconnects while being written to would otherwise kill the scheduler
	{
		//Find out number of gpus.
		FILE* fpcmd=popen("nvidia-smi -L 2>/dev/null", "r");
		if(!fpcmd){
			NGPU=0;
		}else{
			char line[4096];
			while(fgets(line, sizeof(line), fpcmd)){
				if(!strncmp(line, "GPU", 3)){
					NGPU++;
				} else{
					warning_time("unknown entry: %s\n", line);
				}
			}
			pclose(fpcmd);
		}
		dbg("NGPU=%d\n", NGPU);
	}

	//Use ws_service in a separate thread does not update job status.
	//We call it in scheduler_timeout using poll. no need for LOCK in this case
#if HAS_LWS
	start_lws(PORT+100);
#endif
	//Must acquire mutex_sch before handling run_t
	//still need time out to handle process queue.
	listen_port((listen_opt_t){.localpath=slocal, .port=PORT, .responder=respond, .http_responder=http_handshake, .timeout_fun=scheduler_timeout, .timeout_sec=10});
	runned_remove(-1);
	socket_close();
	signal_caught=0;//enable printing memory debug information
	exit(0);
}
