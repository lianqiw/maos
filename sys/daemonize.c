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

/*
  \file daemonize.c

  Process handling routines.
*/
/**
   Ensure exclusive access by opening and maintaining lock of file fn.

   Will proceed if the lock if succeed, in which case non-negative number of the
   fd is returned. If other process already locked, will return negative
   number. If version is specified and is bigger than the value contained in
   existing fnlock, will kill the old process that created fnlock.
*/
int lock_file(const char* fnlock, /**<The filename to lock on*/
	long block,         /**<block on weighting. set to 0 for no waiting.*/
	long version        /**<The version of the software that locks the file, primarily for managing scheduler. May be zero.*/
){
	int fd;
	int count=0;
retry:
	if((count++)>10){
		fd=-1;
		return fd;
	}
	fd=open(fnlock, O_RDWR|O_CREAT, 0644);
	if(fd>=0){
		//fcntl(fd, F_SETFD, FD_CLOEXEC);//do not set CLOEXEC to hold the lock after fork/exec.
		int op=LOCK_EX;
		if(!block) op|=LOCK_NB;
		if(flock(fd, op)){/*lock faild. */
			perror("flock");
			if(block){/*In block mode, we should never fail. */
				error("Lock failed\n");
			}
			long pid=0;
			FILE* fp=fdopen(fd, "r");
			if(!fp||fscanf(fp, "%ld", &pid)==EOF){
				warning("Error reading pid\n");
				fclose(fp);
				fd=-1;
				sleep(1);
				goto retry;
			} else{
				if(kill(pid, 0)){
					warning_time("Unknown process %ld already locks file %s. (NFS mounted system?)\n", pid, fnlock);
					return -1;
				} else{/*already running. check version */
					/*warning("Process %ld already locks file %s\n",pid,fnlock); */
					long version_old=0;
					if(version>0&&(fscanf(fp, "%ld", &version_old)==EOF||version_old<version)){
						info("%d is sending TERM signal to old executive\n", getpid());
						if(!kill(pid, SIGTERM)){//signal sent
							sleep(5);
							if(!kill(pid, 0)){//still running
								warning_time("%d is sending KILL signal to the old executive.\n", getpid());
								if(!kill(pid, SIGKILL)){//signal sent
									sleep(5);
									if(!kill(pid, 0)){
										warning_time("Failed to kill the old executive\n");
										return -pid;
									}
								}
							}
						}
						fclose(fp);/*closed fd also. */
						goto retry;
					} else{
						fclose(fp);
						return -pid;
					}
				}
			}
		} else{/*lock succeed. write pid. */
			char strpid[60];
			snprintf(strpid, 60, "%d %ld\n", getpid(), version);
			lseek(fd, 0, SEEK_SET);
			if(ftruncate(fd, 0)<0)
				warning_time("Unable to truncate the file\n");
			if(write(fd, strpid, strlen(strpid)+1)!=(long)strlen(strpid)+1){
				warning_time("Write pid %d to %s failed\n", getpid(), fnlock);
			}
			fsync(fd);/*don't close file. maintain lock. */
		}
	} else{
		warning_time("Open file failed. This should never happen\n");
		sleep(1);
		goto retry;
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
void single_instance_daemonize(const char* lockfolder_in,
	const char* progname, long version,
	void(*daemon_func)(void*),
	void* daemon_arg){
	int fd=-1;
	char* lockfolder=expand_filename(lockfolder_in);
	char fnlock[PATH_MAX];
	char fnlog[PATH_MAX];
	snprintf(fnlock, PATH_MAX, "%s/%s.pid", lockfolder, progname);
	snprintf(fnlog, PATH_MAX, "%s/%s.log", lockfolder, progname);
	free(lockfolder);
	
	fd=lock_file(fnlock, 0, version);//non-blocking lock
	if(fd<0){
	/*lock failed. daemon already running. no need to start the daemon. */
		if(daemon_func){
			dbg_time("failed to lock_file. return\n");
			return;
		} else{
			_exit(EXIT_SUCCESS);
		}
	}
	/*lock success, forking and put to background. */
	/*fork to detach from terminal. Do the fist fork*/
	long pid=fork();
	if(pid<0){
		exit(EXIT_FAILURE);
	} else if(pid>0){
		close(fd);/*release lock in this process that is not daemon. */
		usleep(1e5);
		waitpid(pid, NULL, 0);/*prevent child from defunct*/
		if(daemon_func){/*The process that launched this routine will return. */
			sleep(3);
			return;
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
	if(!freopen(fnlog, "w", stdout)||dup2(1, 2)==-1) warning("Error redirect stdout or stderr\n");
	setbuf(stdout, NULL);/*disable buffering. */

	char strpid[60];
	snprintf(strpid, 60, "%d %ld\n", getpid(), version);
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
		return;
	}
}
int detached=0;
typedef struct{
	int pfd;//in
	FILE** fps;
	int nfp;
}dup_stdout_t;
/**
  Replicate stream written to stdout to both stdout and file.
* */
static void* dup_stdout(dup_stdout_t* data){
	int fd=data->pfd;//read.
	FILE** fpout=data->fps;
	int nfp=data->nfp;
	char buf[400];
	int nactive=nfp;
	while(nactive){
		int len=read(fd, buf, sizeof buf);
		if(len<=0&&errno==EINTR){
			continue;
		}
		if(len<=0){/*pipe closed. parent exited. we exit also.*/
			break;
		}
		nactive=0;
		for(int i=0; i<2; i++){
			if(fpout[i]){
				if(fwrite(buf, len, 1, fpout[i])!=1){
					fpout[i]=NULL;
				} else{
					nactive++;
				}
			}
		}
	}
	return 0;
}

/**
  Redirect output.
  If we are in detached mode, will output to file, otherwise will output to both file and screen.
*/
FILE* fpconsole=0;
void redirect(void){
	extern int disable_save;
	if(disable_save) return;
	char fn[PATH_MAX];
	snprintf(fn, PATH_MAX, "run_%s_%ld.log", HOST, (long)getpid());
	if(detached){//only output to file
		if(!freopen(fn, "w", stdout)||dup2(1, 2)==-1) warning("Error redirecting stdout or stderr.\n");
		//don't close stdin to prevent fd=0 from being used by file.
		if(!freopen("/dev/null", "r", stdin)) warning("Error redirecting stdin\n");
	} else{
	/* output to both file and screen. we first keep a reference to our
	   console output fd. The stdout and stderr is redirected to one of of
	   the pipe and the other end of the pipe is output to the screen and file.
	 */
		int pfd[2];
		if(pipe(pfd)){//fail to create pipe.
			warning("pipe failed, failed to redirect stdout.\n");
		} else{
			fpconsole=fdopen(dup(fileno(stdout)), "w");
			dup_stdout_t* data=mycalloc(1, dup_stdout_t);
			data->nfp=2;
			data->fps=mycalloc(data->nfp, FILE*);
			data->fps[0]=fpconsole;
			data->fps[1]=fopen(fn, "w");
			data->pfd=pfd[0];//read
			//spawn a thread to duplicate output to both console and file.
			pthread_t thread;
			//child thread read from pfd[0] and write to stdout.
			pthread_create(&thread, NULL, (void* (*)(void*))dup_stdout, data);
			//master threads redirects stderr and stdout to pfd[1]
			if(dup2(pfd[1], 1)==-1||dup2(pfd[1], 2)==-1){
				warning("Error redirecting stdout or stderr");
			}
		}
	}
	setbuf(stdout, NULL);
}
/**
   Daemonize a process by fork it and exit the parent. no need to fork twice since the parent exits.
*/
void daemonize(void){ /* Fork off the parent process */
	pid_t pid=fork();
	if(pid<0){
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
	const char* args=NULL;
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
			cwd=strdup(cmd2);
			if(!exepath){
				exepath=exe;//has to use shell.
			}
			args=stexe;
			break;
		}
	}
	pid_t ans=-1;
	if(args){
		args=strstr(args, exename);
		if(!args){
			warning("args=(%s) has wrong format\n", args);
			goto end;
		}
		args+=strlen(exename)+1;
		ans=spawn_process(exepath, args, cwd);
	} else{
		warning("Unabel to interpret %s\n", cmd);
	}
end:
	free(cwd);
	free(cmd2);
	return ans;
}
/**
   Find an executable from predetermined locations and return the absolute path.
 */
char* find_exe(const char* name){
	char* fn=stradd(BUILDDIR, "/bin/", name, NULL);
	if(!exist(fn)){
		free(fn);
		char* cwd=mygetcwd();
		fn=stradd(cwd, "/", name, NULL);
	}
	if(!exist(fn)){
		free(fn); fn=NULL;
	}
	return fn;
}
/**
   fork twice and launch exename, with arguments args. Returns PID of the grandchild process.
*/
int spawn_process(const char* exename, const char* args, const char* path){
	int pipfd[2];
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
				execl(fn, exename, args, NULL);
			} else{
				execlp(exename, exename, args, NULL);
			}
			_exit(EXIT_FAILURE);
		}
	}
	return ans;
}
