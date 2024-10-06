/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#include <strings.h> //strncasecmp
#include <unistd.h> //sync

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


void wait_cpu(int nthread){
	char fnlock[64];
	snprintf(fnlock, 64, "%s/aos.lock", getenv("HOME"));
	int fd;
	fd=lock_file(fnlock, 1);
	int nvail=0;
	if(nthread>NCPU) nthread=NCPU;
	while((nvail=get_cpu_avail())+1<nthread){
		dbg("needs %d threads but only %d is available\n", nthread-1, nvail);
		sleep(5);
	}
	close(fd);
}

uint16_t PORT=0;
char** hosts=0;
static char** hostsaddr=0;
int nhost=0;
PNEW(mutex_hosts);
/**
   Parse and add host info to hosts[] and hostaddr[].
   hosts contains the shortcut name, without dot
   and hostaddr contains the FQDN and optionally the port
*/
void parse_host(const char* line /**<contains hostname[=hostaddr:port]*/
){
	static size_t memhost=0;
	if(strlen(line)>0&&line[0]!='#'){
		LOCK(mutex_hosts);
		if(memhost<(size_t)(nhost+1)){
			memhost+=10;
			hosts=realloc(hosts, memhost*sizeof(char*));
			hostsaddr=realloc(hostsaddr, memhost*sizeof(char*));
		}

		char* eq=strchr(line, '=');
		hostsaddr[nhost]=eq?strdup(eq+1):strdup(line);
		if(eq) *eq=0;//terminate hostname
		//char *dot=strchr(line, '.');
		//size_t n0=strlen(line);
		//size_t n=MIN(dot?(size_t)(dot-line):n0, eq?(size_t)(eq-line):n0);
		hosts[nhost]=strdup(line); //strndup(line, n?n:n0);
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
/**
 * Initialize hosts and associate an id number */
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
}
/**
   Translate hostname to address based on ~/.aos/hosts
*/
const char* lookup_hostaddr(const char* hostname){
	if(hostname){
		for(int ihost=0; ihost<nhost; ihost++){
			if(!strcmp(hosts[ihost], hostname)){
				return hostsaddr[ihost];
			}
		}
	}
	return hostname;
}
/**
   Translate address to hostname based on ~/.aos/hosts
*/
const char *lookup_hostname(const char *hostaddr){
	if(hostaddr){
		for(int ihost=0; ihost<nhost; ihost++){
			if(hostsaddr[ihost] && !strcmp(hostsaddr[ihost], hostaddr)){
				return hosts[ihost];
			}
		}
	}
	return hostaddr;
}
#if! MAOS_DISABLE_SCHEDULER
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
		execl(fn_scheduler, "scheduler", "--direct", NULL);
	} else{/*fall back. this won't have the right argv set. */
		info("Launch scheduler using shell\n");
		if(execlp("scheduler", "scheduler", "--direct", NULL)){
			warning("scheduler not found\n");
		}
	}
}

static void launch_scheduler(int retry){
	dbg_time("launch_scheduler\n");
	/*launch scheduler if it is not already running. */
	while(single_instance_daemonize(TEMP, "scheduler", scheduler_version,
		(void(*)(void*))launch_scheduler_do, NULL)&&retry>-1){
		retry--;
		sleep(1);
	}
}
#endif
int scheduler_connect(const char *hostname){
#if MAOS_DISABLE_SCHEDULER
	(void)hostname; return -1;
#else		
	return connect_port(hostname, PORT, 0, 0);
#endif
}
/**
   To open a port and connect to scheduler in the local host*/
int scheduler_connect_self(int block){
	int sock=-1;
#if MAOS_DISABLE_SCHEDULER
	(void) block;
#else
	char fn[PATH_MAX];
	int retry=block?100:2;
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
			launch_scheduler(block?10:0);
			mysleep(1);
			retry--;
		}
	}while(sock<0 && retry>0);
#endif	
	return sock;
}
#if! MAOS_DISABLE_SCHEDULER
static int psock=-1; //persistent socket for maos to scheduler connection
#endif
#define CATCH_ERR(A) if(A){psock=-1;}
void scheduler_report_path(const char* path){
#if MAOS_DISABLE_SCHEDULER
	(void)path;
#else
	static char path_save[PATH_MAX];
	path_save[0]=0;
	if(path){
		strncpy(path_save, path, PATH_MAX-1);
	} else if(!path_save[0]){
		return;//do nothing
	}

	if(psock<0){
		psock=scheduler_connect_self(0);//non-blocking connection
	}
	if(psock<0){
		return;
	}

	int cmd[2];
	cmd[0]=CMD_PATH;
	cmd[1]=getpid();
	CATCH_ERR(stwriteintarr(psock, cmd, 2));
	CATCH_ERR(stwritestr(psock, path_save));
#endif	
}


/**
   Called by maos to report a job start to scheduler and wait for signal before proceeding if waiting is set.
 */
void scheduler_start(int nthread, int ngpu, int waiting){
#if MAOS_DISABLE_SCHEDULER
	(void)ngpu;
	if(waiting) wait_cpu(nthread);
#else
	if(psock<0){
		psock=scheduler_connect_self(waiting?1:0);
	}
	if(psock>=0){
		int cmd[4];
		cmd[0]=CMD_START;
		cmd[1]=getpid();
		cmd[2]=nthread;
		cmd[3]=(waiting?1:0)|(ngpu<<1);;
		CATCH_ERR(stwriteintarr(psock, cmd, 4));
		if(waiting){
			int cmd2;
			if(streadint(psock, &cmd2)){
				warning_time("Failed to get answer from scheduler.\n");
			}else{
				waiting=0;
			}
		}
	}else{
		warning_time("Failed to connect to scheduler\n");
	}
	if(waiting){
		wait_cpu(nthread);
	}
#endif
}

/**
   Called by maos to notify scheduler the completion of a job */
void scheduler_finish(int status){
#if MAOS_DISABLE_SCHEDULER
	(void)status;
#else
	if(psock<0){
		psock=scheduler_connect_self(0);//non-blocking connection
	}
	if(psock<0){
		return;
	}
	int cmd[2];
	if(status==0){
		cmd[0]=CMD_FINISH;
	}else if(iscrash(status)){
		cmd[0]=CMD_CRASH;
	}else{
		cmd[0]=CMD_KILLED;
	}
	cmd[1]=getpid();
	CATCH_ERR(stwriteintarr(psock, cmd, 2));
	close(psock);psock=-1;
#endif	
}
#if !MAOS_DISABLE_SCHEDULER
pthread_t cthread=0;
static void* scheduler_connect_thread(void *data){
	(void) data;
	dbg_time("started.\n");
	psock=scheduler_connect_self(1);
	if(psock!=-1){
		scheduler_report_path(NULL);
		cthread=0;
	}
	dbg_time("finished with psock=%d.\n", psock);
	return NULL;
}
#endif
/**
   called by sim.c to report job status. Do not try to reconnect. */
void scheduler_report(status_t* status){
#if MAOS_DISABLE_SCHEDULER
	(void)status;
#else
	if(psock!=-1){
		int cmd[2];
		cmd[0]=CMD_STATUS;
		cmd[1]=getpid();
		CATCH_ERR(stwriteintarr(psock, cmd, 2));
		CATCH_ERR(stwrite(psock, status, sizeof(status_t)));
	}else if(!cthread){//launch a thread to connect
		pthread_create(&cthread, NULL, scheduler_connect_thread, NULL);
	}
#endif
}

/**
   Started by maos to listen to the sock which connects to the
   scheduler for commands.
   It requires a new connection to the scheduler to avoid data racing (?).
*/
pthread_t scheduler_listen(thread_fun fun){
	pthread_t ans=0;
#if MAOS_DISABLE_SCHEDULER
	(void) fun;
#else
	if(!fun) return 0;
	//int	sock=scheduler_connect_self(0);
	if(psock<0){
		psock=scheduler_connect_self(0);
	}
	if(psock!=-1){
		int cmd[2];
		cmd[0]=CMD_MAOSSER;
		cmd[1]=getpid();
		stwriteintarr(psock, cmd, 2);
		ans=thread_new(fun, (void *)(long)psock);
	} else{
		warning_time("Failed to connect to scheduler\n");
	}
#endif
	return ans;
}
/**
   Ask scheduler (maybe in another machine) to launch executable.
*/
int scheduler_launch_exe(const char* host, int argc, const char* argv[]){
	char *scmd=argv2str(argc, argv, " ");
	int ret=0;
#if MAOS_DISABLE_SCHEDULER
	ret=launch_exe(argv[0], scmd)<0?-1:0;
	(void)host;(void)argc;(void)argv;
#else
	int sock=connect_port(host, PORT, 0, 0);
	if(sock<=-1){
		warning("Failed to connect to %s:%d: %s\n", host, PORT, strerror(errno));
		return -1;
	}
	int cmd[2]={CMD_LAUNCH, 2};
	
	if(stwriteintarr(sock, cmd, 2)
		||stwritestr(sock, argv[0])
		||stwritestr(sock, scmd)
		||streadint(sock, &ret)){
		warning("Failed to write to scheduler at %s\n", host);
		ret=-1;
	}
	
	close(sock);
	
#endif
	free(scmd);
	return ret;
}
/**
   send a drawing sock to the scheduler for caching (dir=1) or request a sock for drawing (dir=-1).
   id should be between 1 and 9999.
*/
int scheduler_socket(int dir, int *sfd, int id){
	int ans=-1;
#if MAOS_DISABLE_SCHEDULER
	(void)dir;(void)sfd;(void)id;
#else
	int ssock=scheduler_connect_self(1);
	if(ssock==-1 || !sfd){
		return -1;
	}
	id=abs(id);
	if(id>10000 || id==0) id=1;
	if(dir==1&&*sfd!=-1){//send
		int cmd[2]={CMD_DRAWCLI, id};
		ans=(stwriteintarr(ssock, cmd, 2)||stwritefd(ssock, *sfd));
	}else if(dir==-1 || dir==0){//recv
		int ans2;
		int cmd[2]={CMD_DRAWCLI, dir==-1?(-id):(-10000-id)};
		ans=(stwriteintarr(ssock, cmd, 2)||streadint(ssock, &ans2))||ans2;
		if(!ans){
			ans=streadfd(ssock, sfd);
		}
	}else{
		dbg_time("invalid dir=%d\n", dir);
	}
	close(ssock);
	if(ans){
		dbg2_time("scheduler_socket operation for %s failed.\n", dir==1?"send":"receive");
	}
#endif
	return ans;
}

/**
   Execute addr2line as specified in buf, combine the answer, and return the
   string.*/
int call_addr2line(char* ans, int nans, const char* cmdstr){
	ans[0]=0;
	int res=1;
	FILE *fpcmd=popen(cmdstr, "r");
	if(!fpcmd){
		warning("popen failed: %s", cmdstr);
	} else{
		char line[4096];
		nans--;
		while(fgets(line, sizeof(line), fpcmd)){
			char *tmp=strrchr(line, '/');
			if(tmp){
				tmp++;
				char* tmp2=strchr(tmp, '\n'); if(tmp2) tmp2[0]='\0';
				tmp2=strchr(tmp,'('); if(tmp2) tmp2[-1]='\0';
				if(!mystrcmp(tmp, "sp")||!mystrcmp(tmp, "mat")||!mystrcmp(tmp, "cell.c")||!mystrcmp(tmp, "loc.c")||!mystrcmp(tmp, "mem.c")){
					break;
				}
				strncat(ans, "->", nans); nans-=3;
				strncat(ans, tmp, nans);nans-=strlen(tmp);
				res=0;
			}
		}
		pclose(fpcmd);
	}
	if(res){//failed. copy the progname and addresses
		char *prog=strrchr(cmdstr, '/');
		strncat(ans, prog?(prog+1):cmdstr, nans);
		res=0;
	}
	return res;
}

/**
   Convert backtrace address to source line.
   Do not call routines that may call malloc as the mutex in malloc may already be locked.
   return 0 if success otherwise -1
 */
int print_backtrace_symbol(void* const* buffer, int size){
	int ans=-1;
	//disable memory debugging as this code may be called from within malloc_dbg
//#if (_POSIX_C_SOURCE >= 2||_XOPEN_SOURCE||_POSIX_SOURCE|| _BSD_SOURCE || _SVID_SOURCE) && !defined(__CYGWIN__)
#if __linux__ //mac does not have addr2line. //TODO: use atos
	char cmdstr[PATH_MAX]={0};
	char add[24]={0};
	char progname[PATH_MAX]={0};
	if(get_job_progname(progname, sizeof progname, 0)){
		dbg("Unable to get progname\n");
		return ans;
	}
	snprintf(cmdstr, sizeof cmdstr, "addr2line -f -e %s", progname);
	for(int it=size-1; it>0; it--){
		snprintf(add, 24, " %p", buffer[it]);
		strncat(cmdstr, add, PATH_MAX-strlen(cmdstr)-1);
	}

	static int connect_failed=0;
	PNEW(mutex);//Only one thread can do this.
	LOCK(mutex);
	if(MAOS_DISABLE_SCHEDULER||is_scheduler||connect_failed){
		char line[200];
		if(!call_addr2line(line, sizeof line, cmdstr)){
			dbg("%s\n", line);
			if(strlen(line)){
				ans=0;
			}
		} else{
			char *prog=strrchr(cmdstr, '/');
			dbg("%s\n", prog?(prog+1):cmdstr);
		}
	} else{
#if MAOS_DISABLE_SCHEDULER == 0
		//Create a new connection.
		int sock=scheduler_connect_self(0);
		if(sock!=-1){
			int cmd[2];
			int len;
			cmd[0]=CMD_TRACE;
			cmd[1]=getpid();
			char line[PATH_MAX];//avoid malloc
			if(stwrite(sock, cmd, sizeof(int)*2)){
				dbg("write cmd %d failed\n", cmd[0]);
			} else if(stwritestr(sock, cmdstr)){
				dbg("write cmd %s failed\n", cmdstr);
			} else if(streadint(sock, &len) || len>=PATH_MAX || stread(sock, line, len)){
				dbg("read cmd failed\n");
			} else{
				dbg(" %s\n", line);
				if(strlen(line)){
					ans=0;
				}
			}
			close(sock);
		} else{
			dbg("Failed to connect to scheduler\n");
			connect_failed=1;
		}
#endif
	}
	UNLOCK(mutex);
	//sync();
#else
	(void)buffer; (void)size;
#endif
	return ans;
}
#if !defined(__CYGWIN__) && !defined(__FreeBSD__) && !defined(__NetBSD__)
#include <execinfo.h>
void print_backtrace(){
	int size0, size1;
	size0=1024;
	void* buffer[size0];
	size1=backtrace(buffer, size0);
	print_backtrace_symbol(buffer, size1);
	//sync();
}
#else
void print_backtrace(){}
#endif


