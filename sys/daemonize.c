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


#include <unistd.h>
#include <signal.h>
#include <netdb.h>
#include <fcntl.h> 
#include <errno.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <netinet/in.h>
#include <pthread.h>
#include <ctype.h>
#include "common.h"
#include "process.h"
#include "misc.h"
#include "daemonize.h"
int detached=0;
FILE *fplog=NULL;
/*
  \file daemonize.c

  Process handling routines.
*/
int lock_file(const char *fnlock, /**<The filename to lock on*/
	int block         /**<block on waiting. set to 0 for no waiting.*/
){
	int fd=open(fnlock, O_RDWR|O_CREAT, 0644);
	if(fd<0){
		warning("Open file %s failed: %s\n", fnlock, strerror(errno));
	}else{
		int op=block?F_LOCK:F_TLOCK;
		//lockf() works on both local and NFS
		if(lockf(fd, op, 0)){/*lock faild. another process already locked file.*/
			dbg("Lock %s failed: %s\n", fnlock, strerror(errno));
			close(fd); fd=-1;
		}else{//: success
			dbg("Lock %s success.\n", fnlock);
			char msg[128];
			snprintf(msg, sizeof(msg), "Locked by %s:%d\n", HOST, PID);
			size_t nmsg=strlen(msg)+1;
			if(write(fd, msg, nmsg)!=(ssize_t)nmsg){
				dbg("Failed to write message %s to file %s\n", msg, fnlock);
			}
		}
	}
	return fd;
}
/**
   Ensure exclusive access by opening and maintaining lock of file fn.

   Will proceed if the lock if succeed, in which case non-negative number of the
   fd is returned. If other process already locked, will return negative
   number. If version is specified and is bigger than the value contained in
   existing fnlock, will kill the old process that created fnlock.

   returns:
   -PID (<-1): if a running process already locked file
   -1: failed to lock
   fd (>-1): successfully locked
*/
int lock_file_version(const char* fnlock, /**<The filename to lock on*/
	int block,         /**<block on waiting. set to 0 for no waiting.*/
	int version        /**<The version of the software that locks the file, primarily for managing scheduler. May be zero.*/
){
	int fd;
	int count=0;
	if((count++)>10){
		return -1;
	}
	fd=open(fnlock, O_RDWR|O_CREAT, 0644);
retry:
	if(fd<0){
		warning_time("Open file %s failed: %s\n", fnlock, strerror(errno));
	}else{//check file content first
		int op=LOCK_EX;
		if(!block) op|=LOCK_NB;
		//lockf does not inherit to fork(), so we have to use flock here
		if(flock(fd, op)){/*lock faild. another process already locked file.*/
			if(block){/*In block mode, we should never fail. */
				dbg_time("Blocked Lock failed: %s\n", strerror(errno));
			}
			//lock failed. check version of the running process. kill it if old version
			int pid=0, version_old=0;
			int ans1=0, ans2=-1;
			FILE* fp=fdopen(dup(fd), "r");
			if(fp){//check existing instance
				ans1=fscanf(fp, "%d %d", &pid, &version_old);
				fclose(fp); fp=NULL;
				if(ans1>0 && pid>0){
					if((ans2=kill(pid, 0))){//already exited or in a different host
						pid=0;
					}
				}
				if(pid){//already running
					if(ans1>1 && version && version_old && version>version_old){
						ans2=kill(pid, SIGTERM);
						if(ans2){//signal send failed
							dbg_time("%d failed to send TERM signal to old executive %d, return -1.\n", getpid(), pid);
							return -pid;
						}else{
							dbg_time("%d successfully sends TERM signal to old executive %d\n", getpid(), pid);
							sleep(1);
							goto retry;
						}
					}else{
						dbg_time("Already running as PID %d.\n", pid);
						return -pid;//keep existing drawdaemon
					}
				}
				
			}
		} else{/*lock succeed. write pid. the pid will be updated with correct pid later. */
			char strpid[60];
			snprintf(strpid, 60, "%d %d\n", getpid(), version);
			lseek(fd, 0, SEEK_SET);
			if(ftruncate(fd, 0)<0){
				warning_time("Unable to truncate the file\n");
			}
			if(write(fd, strpid, strlen(strpid)+1)!=(long)strlen(strpid)+1){
				warning_time("Write pid %d to %s failed\n", getpid(), fnlock);
			}
			fsync(fd);/*don't close file. maintain lock. */
		}
	}
	return fd;
}
/**
   Launch a singleton daemon. First try to create a lock file progname.pid
   in folder lockfolder_in. If lock failes, it means other daemon already
   exists, will not proceed. If lock succeed, will write pid and version in the
   file.

   the stdin and stderr stream will then be directed to progname.log in folder
   lockfolder_in. buffering is disabled to ensure real time display by tail.

   if daemon_func is true, will fork and run it with argument
   daemon_arg. otherwise the routine returns.
*/
int single_instance_daemonize(const char* lockfolder,
	const char* progname, int version, void(*daemon_func)(void*), void* daemon_arg){
	int fd=-1;
	char fnlock[PATH_MAX];
	snprintf(fnlock, PATH_MAX, "%s/%s.pid", lockfolder, progname);
	
	fd=lock_file_version(fnlock, 0, version);//non-blocking lock
	if(fd<0){
		/*lock failed. daemon already running. no need to start the daemon. */
		dbg_time("failed to lock_file %s. return\n", fnlock);
		if(daemon_func){
			if(fd==-1){
				return -1;//failed to lock for some reason.
			}else{
				return 0;//already running
			}
		} else{
			_exit(EXIT_SUCCESS);
		}
	}
	/*lock success, forking and put to background. */
	/*fork to detach from terminal. Do the fist fork*/
	pid_t pid=fork();
	if(pid<0){
		dbg_time("failed to fork. return\n");
		if(daemon_func){
			return -1;
		}else{
			exit(EXIT_FAILURE);
		}
	} else if(pid>0){
		close(fd);/*release lock in parent process. */
		mysleep(0.1);
		waitpid(pid, NULL, 0);/*prevent child (parent of another fork) from defunct*/
		if(daemon_func){/*The process that launched this routine will return. */
			return 0;
		} else{
			_exit(EXIT_SUCCESS);
		}
	}
	//child process keeps on

	/* setsid so that this process is the session leader and
	   therefore detached from the parent session and will
	   not hang when the parent session terminates. As a
	   session leader it can open /dev/tty and obtain a
	   controller terminal. We don't want this to happen. So
	   we fork again after setsid.*/
	if(setsid()==-1) warning("Error setsid\n");
	if(chdir(TEMP)) warning("Error chdir\n");
	umask(0077);
	/*We fork again. after this fork, the process is not the
	  session leader (the leader has exited) and no
	  controlling tty can every happen.*/
	pid=fork();
	if(pid<0){
		exit(EXIT_FAILURE);
	} else if(pid>0){/*exit first parent. */
		close(fd);/*close and release lock */
		_exit(EXIT_SUCCESS);
	}
	PID=getpid();
	/*redirect stdin/stdout. */
	if(!freopen("/dev/null", "r", stdin)) warning("Error closing stdin\n");
	//not good to redirect both stdout and stderr to the same file. out of order.
	char fnerr[PATH_MAX];
	char fnlog[PATH_MAX];
	
	snprintf(fnlog, PATH_MAX, "%s/%s.log", lockfolder, progname);
	if(exist(fnlog)){
		snprintf(fnerr, PATH_MAX, "%s/%s.log.%s", lockfolder, progname, myasctime(0));
		(void)rename(fnlog, fnerr);
	}
	snprintf(fnerr, PATH_MAX, "%s/%s.err", lockfolder, progname);
	if(!freopen(fnlog, "w", stdout)) warning("Error redirect stdout\n");
	if(!freopen(fnerr, "w", stderr)) warning("Error redirect stderr\n");
	setbuf(stdout, NULL);/*disable buffering. */

	char strpid[60];
	snprintf(strpid, 60, "%d %d\n", getpid(), version);
	lseek(fd, 0, SEEK_SET);
	if(ftruncate(fd, 0)<0)
		warning("Unable to truncate file\n");
	if(write(fd, strpid, strlen(strpid)+1)!=(long)strlen(strpid)+1){
		warning("Write pid %d to %s failed\n", getpid(), fnlock);
	}
	if(daemon_func){
		daemon_func(daemon_arg);
		exit(EXIT_SUCCESS);/*make sure we don't return. */
	} else{
		return 0;
	}
}
/**
  Redirect output.
  If we are in detached mode, will output to file, otherwise will output to both file and screen.
*/
void redirect(void){
	extern int disable_save;
	if(disable_save) return;
	if(fplog){
		fclose(fplog); fplog=NULL;
	}
	char fnlog[PATH_MAX];
	snprintf(fnlog, PATH_MAX, "run_%s_%ld.log", HOST, (long)getpid());
	if(detached){
		if(!freopen("/dev/null", "r", stdin)||!freopen(fnlog, "w", stdout)||dup2(fileno(stdout), fileno(stderr))==-1){
			dbg("Unable to redirect stdin, stdout, or stderr\n");
		}
	}else{
		fplog=fopen(fnlog, "w");
		setbuf(fplog, NULL);
	}
	setbuf(stdout, NULL);/*disable buffering. */
}
/**
   Daemonize a process by fork it and exit the parent. no need to fork twice since the parent exits.
*/
void daemonize(void){ /* Fork off the parent process */
	pid_t pid=fork();
	if(pid<0){
		perror("fork");
		exit(EXIT_FAILURE);
	}
	if(pid>0){/*exit first parent. */
		_exit(EXIT_SUCCESS);
	}
	/* Create a new SID for the child process */
	if(setsid()==-1) error("Error setsid\n");
	PID=getpid();
	umask(0077);
	detached=1;/*we are in detached mode, disable certain print outs.*/
	redirect();
}
/**
   fork and launch exe as specified in cmd. cmd should composed of the path to
   start the exe, exe name, and parameters. used by maos/skyc.

   if cwd is NULL, cmd must include path information.
   if cwd is not NULL cmd must only contain arguments, without path information.
 */
pid_t launch_exe(const char* exepath, const char* cmd){
	const char* exes[]={"maos", "skyc"}; int nexe=2;
	const char* args[3]={NULL};
	const char* exename=NULL;
	char* cwd=NULL;
	if(exepath){
		exename=strrchr(exepath, '/');
		if(exename){
			exename++;
		} else{
			exename=exepath;
		}
		nexe=1;
		exes[0]=exename;
	}
	/*
	  extract cwd from cmd;
	*/
	char* cmd2;
	if(cmd[0]=='~'){
		cmd2=stradd(HOME, cmd+1, NULL);
	} else{
		cmd2=strdup(cmd);
	}
	//Break cmd into path (cwd), exename (stexe), and arguments (args).
	const char* cmd2end=cmd2+strlen(cmd2);
	for(int iexe=0; iexe<nexe; iexe++){
		const char* exe=exes[iexe];
		char* cmd3=cmd2;
		char* stexe=NULL;
		do{
			stexe=strstr(cmd3, exe);
			cmd3=cmd3+strlen(exe);
		} while(stexe&&cmd3<cmd2end&&(*(stexe-1)!='/'||!isspace((int)stexe[strlen(exe)])));
		if(stexe){
			exename=exe;
			stexe[-1]='\0';
			cwd=cmd2;
			if(!exepath){
				exepath=exe;//has to use shell.
			}
			args[0]=exename;
			args[1]=stexe+strlen(exe)+1;
			break;
		}
	}
	pid_t ans=-1;
	if(args[0]){
		ans=spawn_process(exepath, args, cwd);
	} else{
		warning("Unable to interpret %s\n", cmd);
	}
	free(cmd2);
	return ans;
}
/**
   Find an executable from predetermined locations and return the absolute path.
 */
char* find_exe(const char* name){
	if(exist(name)){
		return strdup(name);
	}
	const char* dirs[]={BUILDDIR, "/usr", "/usr/local", "/opt", "/opt/homebrew"};
	for(unsigned int i=0; i<sizeof(dirs)/sizeof(dirs[0]); i++){
		char* fn=stradd(dirs[i], "/bin/", name, NULL);
		if(exist(fn)){
			return fn;
		}else{
			free(fn);
		}
	}
	return NULL;
}
/**
   fork twice and launch exename, with arguments args. Returns PID of the grandchild process.
*/
int spawn_process(const char *exename, const char *const *args, const char *path){
	int pipfd[2];//pipe to obtain pid of grand child.
	if(pipe(pipfd)){
		warning("unable to create pipe\n");
		return -1;
	}
	pid_t pid=fork();
	pid_t ans=0;
	if(pid<0){//fork failed.
		return -1;
	} else if(pid>0){//parant
		close(pipfd[1]);
		if(read(pipfd[0], &ans, sizeof(pid_t))!=sizeof(pid_t)){
			ans=-1;
		}
		waitpid(pid, NULL, 0);/*wait child*/
		close(pipfd[0]);
	} else{//child
		//redirect stdout, stderr to /dev/null to avoid polluting parent output log.
		//child will redirect it to their own file.
		/*if(!freopen("/dev/null", "w", stdout) || !freopen("/dev/null", "w", stderr)){
			warning("redirect stdout or stderr to /dev/null failed\n");
		}*/
		close(pipfd[0]);
		detached=1;
		setenv("MAOS_DIRECT_LAUNCH", "1", 1);
		pid_t pid2=fork();//fork twice to avoid zoombie process.
		if(pid2){//parent, or fail.
			if(write(pipfd[1], &pid2, sizeof(pid_t))!=sizeof(pid_t)){
				warning("Report pid(%d) failed\n", (int)pid2);
			}
			waitpid(pid2, NULL, WNOHANG);
			_exit(EXIT_SUCCESS);
		} else{//child
			close(pipfd[1]);
			if(setsid()==-1) warning("Error setsid\n");
			if(path&&chdir(path)) warning("Error chdir to %s\n", path);
			//Do not use putenv as it does not copy the string
			char* fn=find_exe(exename);
			if(fn){
				execv(fn, (char *const*)args);
			} else{
				execvp(exename, (char* const*)args);
			}
			warning("exec failed to run %s\n", exename);
			_exit(EXIT_FAILURE);
		}
	}
	return ans;
}
