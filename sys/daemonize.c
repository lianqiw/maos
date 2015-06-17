/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <sys/types.h>
#include <fcntl.h> 
#include <errno.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <netinet/in.h>
#include <sys/wait.h>
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
int lock_file(const char *fnlock, /**<The filename to lock on*/
	      long block,         /**<block on weighting. set to 0 for no waiting.*/
	      long version        /**<The version of the software that locks the file, primarily for managing scheduler. May be zero.*/
	      ){
    int fd;
    int count=0;
 retry:
    if(count>10){
	fd=-1;
	return fd;
    }
    count++;
    fd=open(fnlock,O_RDWR|O_CREAT,0644);
    if(fd>=0){
	int op=LOCK_EX;
	if(!block) op |= LOCK_NB;
	if(flock(fd, op)){/*lock faild. */
	    perror("flock");
	    if(block){/*In block mode, we should never fail. */
		error("Lock failed\n");
	    }
	    long pid=0;
	    FILE *fp=fdopen(fd,"r");
	    if(!fp || fscanf(fp,"%ld",&pid)==EOF){
		warning("Error reading pid\n");
		fclose(fp);
		fd=-1;
		sleep(1);
		goto retry;
	    }else{
		if(kill(pid,0)){
		    warning("Unknown process %ld already locks file %s. (NFS mounted system?)\n" ,
			    pid, fnlock);
		    return -1;
		}else{/*already running. check version */
		    /*warning("Process %ld already locks file %s\n",pid,fnlock); */
		    long version_old=0;
		    if(version>0 && (fscanf(fp,"%ld",&version_old)==EOF || version_old < version)){
			info2("%d is sending TERM signal to old executive\n", getpid());
			if(!kill(pid,SIGTERM)){//signal sent
			    sleep(5);
			    if(!kill(pid,0)){//still running
				warning_time("%d is sending KILL signal to the old executive.\n", getpid());
				if(!kill(pid,SIGKILL)){//signal sent
				    sleep(5);
				    if(!kill(pid,0)){
					warning_time("Failed to kill the old executive\n");
					return -pid;
				    }
				}
			    }
			}
			fclose(fp);/*closed fd also. */
			goto retry;
		    }else{
			fclose(fp);
			return -pid;
		    }
		}
	    }
	}else{/*lock succeed. write pid. */
	    char strpid[60];
	    snprintf(strpid,60,"%d %ld\n",getpid(),version);
	    lseek(fd,0,SEEK_SET);
	    if(ftruncate(fd,0)<0)
		warning("Unable to truncate the file\n");
	    if(write(fd,strpid,strlen(strpid)+1)!=strlen(strpid)+1){
		warning("Write pid %d to %s failed\n",getpid(),fnlock);
	    }
	    fsync(fd);/*don't close file. maintain lock. */
	}
    }else{
	warning("Open file failed. This should never happen\n");
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
void single_instance_daemonize(const char *lockfolder_in, 
			       const char *progname, long version,
			       void(*daemon_func)(void*), 
			       void* daemon_arg){
    int fd=-1;
    char *lockfolder=expand_filename(lockfolder_in);
    char *fnlock0=stradd(lockfolder,"/",progname,".pid",NULL);
    char *fnlog0=stradd(lockfolder,"/",progname,".log",NULL);
    /*Make those stack variable so that we don't have to free them. */
    char *fnlock=alloca(strlen(fnlock0)+1);
    strcpy(fnlock,fnlock0);
    char *fnlog=alloca(strlen(fnlog0)+1);
    strcpy(fnlog,fnlog0);
    free(lockfolder);
    free(fnlog0);
    free(fnlock0);

    fd=lock_file(fnlock,0,version);
    if(fd<0){
	/*lock failed. daemon already running. no need to start the daemon. */
	if(daemon_func){
	    return;
	}else{
	    _exit(EXIT_SUCCESS);
	}
    }
    /*lock success, forking and put to background. */
    /*fork to detach from terminal. Do the fist fork*/
    long pid = fork();
    if(pid<0){
	exit(EXIT_FAILURE);
    }else if (pid > 0) {
	close(fd);/*release lock in this process that is not daemon. */
	sleep(1);
	waitpid(pid,NULL,0);/*prevent child from defunct*/
	if(daemon_func){/*The process that launched this routine will return. */
	    return;
	}else{
	    exit(EXIT_SUCCESS);
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
    if(chdir("/")) warning("Error chdir\n");
    umask(0077);
    /*We fork again. after this fork, the process is not the
      session leader (the leader has exited) and no
      controlling tty can every happen.*/
    pid = fork();
    if(pid<0){
	exit(EXIT_FAILURE);
    }else if (pid > 0) {/*exit first parent. */
	close(fd);/*close and release lock */
	exit(EXIT_SUCCESS);
    }
  
    /*redirect stdin/stdout. */
    if(!freopen("/dev/null","r",stdin)) warning("Error closing stdin\n");
    if(!freopen(fnlog, "w", stdout)) warning("Error redirect stdout\n");
    if(!freopen(fnlog, "w", stderr)) warning("Error redirect stderr\n");
    setbuf(stdout,NULL);/*disable buffering. */
    setbuf(stderr,NULL);
 
    char strpid[60];
    snprintf(strpid,60,"%d %ld\n",getpid(),version);
    lseek(fd,0,SEEK_SET);
    if(ftruncate(fd,0)<0)
	warning("Unable to truncate file\n");
    if(write(fd,strpid,strlen(strpid)+1)!=strlen(strpid)+1){
	warning("Write pid %d to %s failed\n",getpid(),fnlock);
    }
    if(daemon_func){
	daemon_func(daemon_arg);
	exit(EXIT_SUCCESS);/*make sure we don't return. */
    }else{
	return;
    }
}
int detached=0;
typedef struct{
    int pfd;
    int stdoutfd;
    const char *fn;
}redirect_t;

static void* fputs_stderr(redirect_t *data){
    int fd=data->pfd;
    int stdoutfd=data->stdoutfd;
    const char *fn=data->fn;
    free(data);
    FILE *fpout[2];
    fpout[0]=fdopen(stdoutfd, "w");
    fpout[1]=fopen(fn, "w");
    setbuf(fpout[0], NULL);
    setbuf(fpout[1], NULL);
    if(!fpout[0] || !fpout[1]) {
	printf("open stderr failed\n");
	return 0;
    }
    char buf[400];
    while(fpout[0] || fpout[1]){
	int len=read(fd, buf, sizeof buf);
	if(len<=0 && errno == EINTR){
	    continue;
	}
	if(len<=0){/*pipe closed. parent exited. we exit also.*/
	    break;
	}
	for(int i=0; i<2; i++){
	    if(fpout[i] && fwrite(buf, len, 1, fpout[i])!=1){
		fpout[i]=NULL;
	    }
	}
    }
    return 0;
}
/**
   Redirect stdout (1) and stderr (2) to fd
 */
static void redirect2fd(int fd){
    if(dup2(fd, 1)<0 || dup2(fd, 2)<0){
	warning("Error redirecting stdout/stderr\n");
    }
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
}
/*
   Redirect stdout and stderr to fn
 */
static void redirect2fn(const char *fn){
    if(!freopen(fn, "w", stdout)){
	warning("Error redirecting stdout/stderr\n");
    }
    //Redirect stderr to stdout 
    dup2(fileno(stdout), fileno(stderr));
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
}
/**
  Redirect output. 
  If we are in detached mode, will output to file, otherwise will output to both file and screen.
*/
void redirect(void){
    extern int disable_save;
    if(disable_save) return;
    char fn[PATH_MAX];
    snprintf(fn, PATH_MAX, "run_%s_%ld.log", myhostname(), (long)getpid());
    (void)remove(fn);
    if(detached){//only output to file
	redirect2fn(fn);
	if(!freopen("/dev/null","r",stdin)) warning("Error redirectiont stdin\n");
    }else{
	/* output to both file and screen. we first keep a reference to our
	   console output fd. The stdout and stderr is redirected to one of of
	   the pipe and the other end of the pipe is output to the screen and file.
	 */
	int stdoutfd=dup(fileno(stdout));
	int pfd[2];
	if(pipe(pfd)){//fail to create pipe.
	    warning("pipe failed\n");
	}else{
	    redirect_t *data=calloc(1, sizeof(redirect_t));
	    data->pfd=pfd[0];//read
	    data->stdoutfd=stdoutfd;//write to console
	    data->fn=fn;
	    //spawn a thread to handle output.
	    pthread_t thread;
	    //child thread read from pfd[0] and write to stdout.
	    pthread_create(&thread, NULL, (void *(*)(void *))fputs_stderr, data);
	    //master threads redirects stderr and stdout to pfd[1]
	    redirect2fd(pfd[1]);
	}
    }
    remove("run_recent.log");
    mysymlink(fn, "run_recent.log");

}
/**
   Daemonize a process by fork it and exit the parent. no need to fork twice since the parent exits.
*/
void daemonize(void){ /* Fork off the parent process */       
    pid_t pid = fork();
    if (pid < 0) {
	exit(EXIT_FAILURE);
    }
    if (pid > 0) {/*exit first parent. */
	_exit(EXIT_SUCCESS);
    }
    /* Create a new SID for the child process */
    if(setsid()==-1) error("Error setsid\n");
    pid=getpid();
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
pid_t launch_exe(const char *exepath, const char *cmd){
    const char *exes[]={"maos", "skyc"}; int nexe=2;
    const char *args=NULL;
    const char *exename=NULL;
    char *cwd=NULL;
    if(exepath){
	exename=strrchr(exepath, '/');
	if(exename){
	    exename++;
	}else{
	    exename=exepath;
	}
	nexe=1;
	exes[0]=exename;
    }
    /*
      extract cwd from cmd;
    */
    char *stexe=NULL;
    char *cmd2;
    if(cmd[0]=='~'){
	cmd2=stradd(HOME, cmd+1, NULL);
    }else{
	cmd2=strdup(cmd);
    }
    const char *cmd2end=cmd2+strlen(cmd2);
    for(int iexe=0; iexe<nexe; iexe++){
	const char *exe=exes[iexe];
	const char *cmd3=cmd2;
	do{
	    stexe=strstr(cmd3, exe);
	    cmd3=cmd3+strlen(exe);
	}while(stexe && cmd3<cmd2end && (*(stexe-1)!='/'||!isspace((int)stexe[strlen(exe)])));
	if(stexe){
	    exename=exe;
	    stexe[-1]='\0';
	    cwd=strdup(cmd2);
	    if(!exepath) {
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
	int pipfd[2];
	if(pipe(pipfd)){
	    warning("unable to create pipe\n");
	    goto end;
	}
	pid_t pid=fork();
	if(pid<0){
	    goto end;//unable to fork
	}else if(pid>0){//parant
	    close(pipfd[1]);
	    waitpid(pid, NULL, 0);/*wait child*/
	    if(read(pipfd[0], &ans, sizeof(pid_t))!=sizeof(pid_t)){
		ans=-1;
	    }
	    close(pipfd[0]);
	}else{//child
	    close(pipfd[0]);
	    detached=1;
	    pid_t pid2=fork();
	    if(pid2<0){//error forking
		if(write(pipfd[1], &pid2, sizeof(pid_t))!=sizeof(pid_t)){
		    warning("Report pid(%d) failed\n", (int)pid2);
		}
		_exit(EXIT_FAILURE);//unable to fork
	    }else if(pid2>0){//parent
		usleep(1000000);
		waitpid(pid2, NULL, WNOHANG);
		//Must waitpid before running kill because otherwise child maybe zoombie.
		if(kill(pid2, 0)){//exec failed.
		    warning("child pid %d not found. exec failed?\n", pid2);
		}
		if(write(pipfd[1], &pid2, sizeof(pid_t))!=sizeof(pid_t)){
		    warning("Report pid(%d) failed\n", (int)pid2);
		}
		_exit(EXIT_SUCCESS);
	    }else{//child
		close(pipfd[1]);
		if(setsid()==-1) warning("Error setsid\n");
		if(chdir(cwd)) error("Error chdir to %s\n", cwd);
		//Do not use putenv as it does not copy the string
		setenv("MAOS_DIRECT_LAUNCH", "1", 1);
		if(execlp(exepath, exepath, args, NULL)){
		    error("Unable to exec: %s\n", strerror(errno));
		}
		_exit(EXIT_FAILURE);
	    }
	}
    }else{
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
char *find_exe(const char *name){
    char *fn=stradd(BUILDDIR, "/bin/", name,NULL);
    if(!exist(fn)){
	free(fn);
	char *cwd=mygetcwd();
	fn=stradd(cwd, "/", name, NULL);
    }
    if(!exist(fn)){
	free(fn); fn=NULL;
    }
    return fn;
}
/**
   fork and launch drawdaemon
*/
int spawn_drawdaemon(int sock){
    pid_t pid;
    if((pid=fork())<0){
	warning("forked failed\n");
	close(sock);
	return -1;
    }else if(pid>0){
	close(sock);
	waitpid(pid, NULL, 0);
	return 0;
    }else{
	detached=1;
	setsid();
	pid=fork();
	if(pid<0){
	    exit(1);
	}else if(pid>0){
	    close(sock);
	    exit(0);
	}
	char arg1[20];
	snprintf(arg1, 20, "%d", sock);
	char *fn=find_exe("drawdaemon");
	if(fn){
	    execl(fn, "drawdaemon", arg1,  NULL);
	}else{
	    execlp("drawdaemon", "drawdaemon", arg1, NULL);
	}
	exit(0);//in case child comes here. quit.
    }
}
