/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

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

#include <sys/types.h>
#include <fcntl.h> 
#include <errno.h>
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

#include "sock.h"
#include "sockio.h"
#include "process.h"
#include "misc.h"
#include "hashlittle.h"
#include "scheduler_client.h"
#include "daemonize.h"

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
void parse_host(char* line){
	(void)line;
}
const char* lookup_host(const char* hostname){
	return hostname;
}
void init_hosts(){}
void free_hosts(){}
void scheduler_start(char* path, int nthread, int ngpu, int waiting){
	(void)path;
	(void)nthread;
	(void)waiting;
	(void)ngpu;
}
int scheduler_wait(void){
	return 0;
}
void scheduler_finish(int status){
	(void)status;
}
void scheduler_report(STATUS_T* status){
	(void)status;
}
int scheduler_listen(thread_fun fun){
	(void)fun;
	return 0;
}
int scheduler_launch_exe(const char* host, int argc, const char* argv[]){
	(void)host;
	(void)argc;
	(void)argv;
	return -1;
}
int scheduler_send_socket(int sfd, int id){
	(void)sfd;
	(void)id;
	return -1;
}
int scheduler_recv_socket(int* sfd, int id){
	(void)sfd;
	(void)id;
	return -1;
}
int scheduler_launch_drawdaemon(void){
	return -1;
}
#else
uint16_t PORT=0;
char** hosts=0;
static char** hostsaddr=0;
int nhost=0;
PNEW(mutex_hosts);
/**
   Parse and add host info to hosts[] and hostaddr[]. The format is hostname[=hostaddr:port]
*/
void parse_host(char* line){
	static int memhost=0;
	if(strlen(line)>0&&line[0]!='#'){
		LOCK(mutex_hosts);
		if(memhost<nhost+1){
			memhost+=10;
			hosts=realloc(hosts, memhost*sizeof(char*));
			hostsaddr=realloc(hostsaddr, memhost*sizeof(char*));
		}

		char* eq=strchr(line, '=');
		if(eq){
			eq[0]='\0'; eq++;
		} else{
			eq=line;
		}
		hosts[nhost]=strdup(line);
		hostsaddr[nhost]=strdup(eq);
		nhost++;
		UNLOCK(mutex_hosts);
	}
}
/**
   free hosts[] and hostsaddr[]
 */
void free_hosts(){
	LOCK(mutex_hosts);
	for(int ihost=0; ihost<nhost; ihost++){
		free(hosts[ihost]);
		free(hostsaddr[ihost]);
	}
	free(hosts);
	free(hostsaddr);
	nhost=0;
	UNLOCK(mutex_hosts);
}
/*

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
	}*/

/*Initialize hosts and associate an id number */
//static __attribute__((constructor))
void init_hosts(){
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX, "%s/.aos/port", HOME);
	PORT=0;

	FILE* fp=fopen(fn, "r");
	if(fp){
		if(fscanf(fp, "%hu", &PORT)!=1){
		}
		fclose(fp);
	}
	if(PORT<1000){ //user dependent PORT to avoid conflict 
		PORT=(uint16_t)((uint16_t)(hashlittle(USER, strlen(USER), 0)&0x2FFF)|10000);
	}

	snprintf(fn, PATH_MAX, "%s/.aos/hosts", HOME);
	fp=fopen(fn, "r");
	if(fp){
		char line[64];
		while(fscanf(fp, "%s ", line)==1){
			parse_host(line);
		}
		fclose(fp);
	}
	/*if(myhostid(HOST)==-1){
	parse_host("localhost");//use local machine
	}*/
}
/**
   Translate hostname to based on ~/.aos/hosts
*/
const char* lookup_host(const char* hostname){
	if(!hostname) return NULL;
	for(int ihost=0; ihost<nhost; ihost++){
		if(!strcmp(hosts[ihost], hostname)){
			return hostsaddr[ihost];
		}
	}
	return hostname;
}
/**
   Launch the scheduler. We already obtained singleton lock and is in a forked process.
*/
static void launch_scheduler_do(void* junk){
	(void)junk;
#if defined(__CYGWIN__)
	char* fn_scheduler=stradd(BUILDDIR, "/bin/scheduler.exe", NULL);
#else
	char* fn_scheduler=stradd(BUILDDIR, "/bin/scheduler", NULL);
#endif
	if(exist(fn_scheduler)){
		dbg_time("Run %s as %d\n", fn_scheduler, getpid());
		execl(fn_scheduler, "scheduler", NULL);
	} else{/*fall back. this won't have the right argv set. */
		info("Launch scheduler using shell\n");
		if(execlp("scheduler", "scheduler", NULL)){
			warning("scheduler not found\n");
		}
	}
}

static void launch_scheduler(void){
	dbg_time("launch_scheduler");
	/*launch scheduler if it is not already running. */
	single_instance_daemonize(TEMP, "scheduler", scheduler_version,
		(void(*)(void*))launch_scheduler_do, NULL);
	sleep(1);
}

/**
   To open a port and connect to scheduler in the local host*/
static int scheduler_connect_self(int block){
	char fn[PATH_MAX];
	int sock=-1;
	int retry=0;
	do{
		if(TEMP[0]=='/'){//try local connection first.
			snprintf(fn, PATH_MAX, "%s/scheduler", TEMP);
			sock=connect_port(fn, PORT, 0, 0);
		}
		if(sock<0){//try socket connection
			snprintf(fn, PATH_MAX, "localhost");
			sock=connect_port(fn, PORT, 0, 0);
		}
		if(sock<0){
	        /*start the scheduler if it is not running*/
			launch_scheduler();
		}
		retry++;
	} while(sock<0&&block&&retry<3);
	return sock;
}

static int psock=-1;

static void scheduler_report_path(char* path){
	static char path_save[PATH_MAX];
	path_save[0]=0;
	if(psock==-1){
		return;
	}
	if(path){
		strcpy(path_save, path);
	} else{
		if(!path_save[0]){
			strcpy(path_save, "unknown");
		}
	}
	int cmd[2];
	cmd[0]=CMD_PATH;
	cmd[1]=getpid();
	stwriteintarr(psock, cmd, 2);
	stwritestr(psock, path_save);
}
#define CATCH_ERR(A) if(A){psock=-1;}

/**
   Started by maos to listen to the sock which connects to the
   scheduler for commands
*/
int scheduler_listen(thread_fun fun){
	if(psock!=-1&&fun){
		thread_new(fun, &psock);
		return 0;
	} else{
		return -1;
	}
}

/**
   Called by maos to report a job start to scheduler.
 */
void scheduler_start(char* path, int nthread, int ngpu, int waiting){
	psock=scheduler_connect_self(1);
	if(psock==-1){
		warning_time("Failed to connect to scheduler\n");
	}
	scheduler_report_path(path);
	int cmd[4];
	cmd[0]=CMD_START;
	cmd[1]=getpid();
	cmd[2]=nthread;
	cmd[3]=(waiting?1:0)|(ngpu<<1);;
	CATCH_ERR(stwriteintarr(psock, cmd, 4));
}

/**
   Called by maos to wait for go signal from scheduler.
*/
int scheduler_wait(void){
	if(psock==-1){
		warning_time("Failed to connect to scheduler\n");
		return -1;
	}
	/*read will block until clearance is received. */
	int cmd;
	if(streadint(psock, &cmd)){
		warning("Failed to get answer from scheduler.\n");
		return -1;
	} else{
		return 0;
	}
}
/**
   Called by maos to notify scheduler the completion of a job */
void scheduler_finish(int status){
	if(psock!=-1){
		int cmd[2];
		if(status==0){
			cmd[0]=CMD_FINISH;
		} else{
			cmd[0]=CMD_CRASH;
		}
		cmd[1]=getpid();
		CATCH_ERR(stwriteintarr(psock, cmd, 2));
		close(psock);psock=-1;
	}
}

/**
   called by sim.c to report job status */
void scheduler_report(STATUS_T* status){
	if(psock==-1){
		psock=scheduler_connect_self(0);
		scheduler_report_path(NULL);
	}
	if(psock!=-1){
		int cmd[2];
		cmd[0]=CMD_STATUS;
		cmd[1]=getpid();
		CATCH_ERR(stwriteintarr(psock, cmd, 2));
		CATCH_ERR(stwrite(psock, status, sizeof(STATUS_T)));
	}
}

/**
   Ask scheduler (maybe in another machine) to launch executable.
*/
int scheduler_launch_exe(const char* host, int argc, const char* argv[]){
	int ret=0;
	int sock=connect_port(host, PORT, 0, 0);
	if(sock<=-1){
		warning("Failed to connect to %s:%d: %s\n", host, PORT, strerror(errno));
		return -1;
	}
	int cmd[2]={CMD_LAUNCH, 2};
	char* scmd=argv2str(argc, argv, " ");
	if(stwriteintarr(sock, cmd, 2)
		||stwritestr(sock, argv[0])
		||stwritestr(sock, scmd)
		||streadint(sock, &ret)){
		warning("Failed to write to scheduler at %s\n", host);
		ret=-1;
	}
	free(scmd);
	close(sock);
	return ret;
}

/**
   send a sock to the scheduler for caching
*/
int scheduler_send_socket(int sfd, int id){
	int ans=-1;
	int ssock=scheduler_connect_self(0);
	if(ssock!=-1&&sfd!=-1){
		int cmd[2]={CMD_SOCK, id?abs(id):1};
		if(stwriteintarr(ssock, cmd, 2)||stwritefd(ssock, sfd)){
			ans=-1;
			warning("Talk to scheduler failed\n");
		} else{
			ans=0;
		}
		close(ssock);
	}
	return ans;
}
/**
   get a socket from the scheduler for reuse. return 0 when success.
*/
int scheduler_recv_socket(int* sfd, int id){
	int ans=-1;
	int ssock=scheduler_connect_self(0);
	if(ssock!=-1){
		int cmd[2]={CMD_SOCK, (id==0?0:-abs(id))};
		int ans2=-1;
		if(stwriteintarr(ssock, cmd, 2)||streadint(ssock, &ans2)){
			warning("Talk to scheduler failed\n");
		} else if(!ans2&&!streadfd(ssock, sfd)){
			ans=0;
			info("received %d from scheduler\n", *sfd);
		} else{
			//dbg("scheduler had no valid fd\n");
		}
		close(ssock);
	} else{
		warning_time("Failed to connect to scheduler\n");
	}
	return ans;
}

#endif /*MAOS_DISABLE_SCHEDULER*/

/**
   Execute addr2line as specified in buf, combine the answer, and return the
   string.*/
int call_addr2line(char* ans, int nans, const char* buf){
	ans[0]=0;
	FILE* fpcmd=popen(buf, "r");
	if(!fpcmd){
		return 1;
	} else{
		char line[4096];
		while(fgets(line, sizeof(line), fpcmd)){
			char* tmp=strrchr(line, '/');
			if(tmp){
				tmp++;
				char* tmp2=strchr(tmp, '\n'); if(tmp2) tmp2[0]='\0';
				strncat(ans, "->", nans-3);
				strncat(ans, tmp, nans-strlen(tmp)-1);
			}
		}
		pclose(fpcmd);
	}
	return 0;
}

/**
   Convert backtrace address to source line.
   Do not call routines that may call malloc as the mutex in malloc may already be locked.
 */
void print_backtrace_symbol(void* const* buffer, int size){
	//disable memory debugging as this code may be called from within malloc_dbg
#if (_POSIX_C_SOURCE >= 2||_XOPEN_SOURCE||_POSIX_SOURCE|| _BSD_SOURCE || _SVID_SOURCE) && !defined(__CYGWIN__)
	static int connect_failed=0;
	char cmdstr[PATH_MAX+40]={0};
	char add[24]={0};
	char progname[PATH_MAX+20]={0};
	if(get_job_progname(progname, sizeof progname, 0)){
		warning("Unable to get progname\n");
		return;
	}
	snprintf(cmdstr, sizeof cmdstr, "addr2line -f -e %s", progname);
	for(int it=size-1; it>0; it--){
		snprintf(add, 24, " %p", buffer[it]);
		strncat(cmdstr, add, PATH_MAX-strlen(cmdstr)-1);
	}
	if(connect_failed){
		dbg("%s\n", cmdstr);
		return;
	}
	PNEW(mutex);//Only one thread can do this.
	LOCK(mutex);
	if(MAOS_DISABLE_SCHEDULER||is_scheduler){
		dbg("backtrace directly\n");
		char ans[10000];
		if(!call_addr2line(ans, 10000, cmdstr)){
			info("%s\n", ans);
		} else{
			info("Command failed\n");
		}
	} else{//Create a new socket and ask scheduler to do popen and return answer.
#if MAOS_DISABLE_SCHEDULER == 0
	//Create a new connection.
		int sock=scheduler_connect_self(0);
		if(sock!=-1){
			int cmd[2];
			cmd[0]=CMD_TRACE;
			cmd[1]=getpid();
			char* ans=NULL;
			if(stwrite(sock, cmd, sizeof(int)*2)){
				warning("write cmd %d failed\n", cmd[0]);
			} else if(stwritestr(sock, cmdstr)){
				warning("write cmd %s failed\n", cmdstr);
			} else if(streadstr(sock, &ans)){
				warning("read cmd failed\n");
			} else{
				info(" %s\n", ans); free(ans);
			}
			close(sock);
		} else{
			warning("Failed to connect to scheduler\n");
			connect_failed=1;
		}
#else
		info("MAOS_DISABLE_SCHEDULER is set\n");
#endif
	}
	UNLOCK(mutex);
	sync();
#else
	(void)buffer; (void)size;
#endif
}
#if !defined(__CYGWIN__) && !defined(__FreeBSD__) && !defined(__NetBSD__)
#include <execinfo.h>
void print_backtrace(){
	int size0, size1;
	size0=1024;
	void* buffer[size0];
	size1=backtrace(buffer, size0);
	print_backtrace_symbol(buffer, size1);
	sync();
}
#else
void print_backtrace(){}
#endif


