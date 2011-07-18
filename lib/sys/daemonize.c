/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
// Check whether scheduler is already running.
#include <stdio.h>
#include <stdlib.h>
#include <netdb.h>
#include <unistd.h>
#include <sys/types.h>
#include <fcntl.h> 
#include <errno.h>
#include <arpa/inet.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/file.h>
#include <sys/stat.h>
#include <netinet/in.h>
#include <sys/wait.h>
#include "common.h"
#include "misc.h"
#include "daemonize.h"

/**
   Open and maintain lock of file fn.
   return negative number if failed.
*/
int lock_file(const char *fnlock, /**<The filename to lock on*/
	      long block,         /**<block on weighting. set to 0 for no waiting.*/
	      long version        /**<The version of the software that locks the file, primarily for managing scheduler*/
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
    //cloexec(fd);//do not do this. need to carry over to forked routines
    if(fd>=0){
	int op=LOCK_EX;
	if(!block) op |= LOCK_NB;
	if(flock(fd, op)){//lock faild.
	    if(block){//In block mode, we should never fail.
		perror("flock");
		error("Lock failed\n");
	    }
	    long pid=0;
	    FILE *fp=fdopen(fd,"r");
	    if(fscanf(fp,"%ld",&pid)==EOF){
		warning("Error reading pid\n");
		fclose(fp);
		fd=-1;
		sleep(1);
		goto retry;
	    }else{
		if(kill(pid,0)){
		    warning("Process %ld already locks file, but it already exited\n",pid);
		    fclose(fp);
		    remove(fnlock);
		    goto retry;
		}else{//already running. check version
		    //warning("Process %ld already locks file %s\n",pid,fnlock);
		    long version_old=0;
		    if(version>0 && (fscanf(fp,"%ld",&version_old)==EOF || version_old < version)){
			info("Killing old executive\n");
			if(kill(pid,SIGTERM) && (sleep(2) && kill(pid,0))){
			    warning3("Failled to kill the old executive, send KILL signal\n");
			    if(kill(pid,SIGKILL) && (sleep(2) && kill(pid,0))){
				warning3("Stilled failled to kill the old executive\n");
				return -pid;
			    }
			}
			fclose(fp);//closed fd also.
			goto retry;
		    }else{
			fclose(fp);
			return -pid;
		    }
		}
	    }
	}else{//lock succeed. write pid.
	    char strpid[60];
	    snprintf(strpid,60,"%ld %ld\n",(long)getpid(),version);
	    lseek(fd,0,SEEK_SET);
	    if(ftruncate(fd,0)<0)
		warning("Unable to truncate the file\n");
	    if(write(fd,strpid,strlen(strpid)+1)!=strlen(strpid)+1){
		warning("Write pid %d to %s failed\n",(int)getpid(),fnlock);
	    }else{
		//info("Write pid %d to %s success\n",(int)getpid(),fnlock);
	    }
	    fsync(fd);//don't close file. maintain lock.
	}
    }else{
	warning("Open file failed. This should never happen\n");
	sleep(1);
	goto retry;
    }
    return fd;
}

void single_instance_daemonize(const char *lockfolder_in, 
			       const char *progname, long version,
			       void(*daemon_func)(void*), 
			       void* daemon_arg){
    /*
      if daemon_func is true, will run the daemon_func in forked daemon 
      and never return. otherwise, the function will return.
      daemon_arg is the arguments to daemon_func
    */
    int fd=-1;
    char *lockfolder=NULL;
    expand_filename(&lockfolder,lockfolder_in);
    char *fnlock0=stradd(lockfolder,"/",progname,".pid",NULL);
    char *fnlog0=stradd(lockfolder,"/",progname,".log",NULL);
    //Make those stack variable so that we don't have to free them.
    char *fnlock=alloca(strlen(fnlock0)+1);
    strcpy(fnlock,fnlock0);
    char *fnlog=alloca(strlen(fnlog0)+1);
    strcpy(fnlog,fnlog0);
    free(lockfolder);
    free(fnlog0);
    free(fnlock0);

    fd=lock_file(fnlock,0,version);
    if(fd<0){//lock failed. daemon already running. no need to start the daemon.
	if(daemon_func){
	    return;
	}else{
	    _exit(EXIT_SUCCESS);
	    return;//just in case
	}
    }
    //lock success, forking and put to background.
    /*fork to detach from terminal. Do the fist fork*/
    long pid = fork();
    if(pid<0){
	exit(EXIT_FAILURE);
    }else if (pid > 0) {//exit first parent.
	close(fd);//release lock in this process that is not daemon.
	waitpid(pid,NULL,0);//prevent child from defunct.
	if(daemon_func){//The process that launched this routine will return.
	    return;
	}else{
	    exit(EXIT_SUCCESS);
	}
    }
    /* setsid so that this process is the session leader and
       therefore detached from the parent session and will
       not hang when the parent session terminates. As a
       session leader it can open /dev/tty and obtain a
       controller terminal. We don't want this to happen. So
       we fork again after setsid.*/
    if(setsid()==-1) warning("Error setsid\n");
    if(chdir("/")) warning("Error chdir\n");
    umask(0077);
    //redirect stdin/stdout.
    if(!freopen("/dev/null","r",stdin)) warning("Error closing stdin\n");
    if(!freopen(fnlog, "w", stdout)) warning("Error redirect stdout\n");
    if(!freopen(fnlog, "w", stderr)) warning("Error redirect stderr\n");
    setbuf(stdout,NULL);//disable buffering.
    setbuf(stderr,NULL);

    /*We fork again. after this fork, the process is not the
      session leader (the leader has exited) and no
      controlling tty can every happen.*/
    pid = fork();
    if(pid<0){
	exit(EXIT_FAILURE);
    }else if (pid > 0) {//exit first parent.
	close(fd);//close and release lock
	exit(EXIT_SUCCESS);
    }
  
    char strpid[60];
    snprintf(strpid,60,"%ld %ld\n",(long)getpid(),version);
    lseek(fd,0,SEEK_SET);
    if(ftruncate(fd,0)<0)
	warning("Unable to truncate file\n");
    if(write(fd,strpid,strlen(strpid)+1)!=strlen(strpid)+1){
	warning("Write pid %d to %s failed\n",(int)getpid(),fnlock);
    }else{
	info("Write pid %d to %s success\n",(int)getpid(),fnlock);
    }
    
    //fsync(fd);
    //maintain lock. don't close fd
    //return to main program and continue in the newly process.
    if(daemon_func){
	daemon_func(daemon_arg);
	exit(EXIT_SUCCESS);//make sure we don't return.
    }else{
	return;
    }
}
int detached=0;
/**
   Daemonize a process by fork it and exit the parent. no need to fork twice since the parent exits.
*/
void daemonize(void){
    /* Fork off the parent process */       
    pid_t pid = fork();
    if (pid < 0) {
	exit(EXIT_FAILURE);
    }
    if (pid > 0) {//exit first parent.
	/*give enough time for the job to communicate with the scheduler so that
	  our jobs are in order.*/
	usleep(50000);
	_exit(EXIT_SUCCESS);
    }
    /* Create a new SID for the child process */
    if(setsid()==-1) error("Error setsid\n");
    
    umask(0077);
    detached=1;/*we are in detached mode, disable certain print outs.*/
    /*
      Closing stdin makes a file opened after freopen for stdout and stderr
      to have an fd of 2, which is the stderr fd. then after another
      freopen, the file ended put being the stderr stream out put. messes up
      everything.
    */
    /*It may be the following producing garbage characters.*/
    if(!freopen("/dev/null","r",stdin)) warning("Error redirectiont stdin\n");
    char fn[256];
    pid=getpid();
    snprintf(fn,256,"run_%d.log",pid);
    if(!freopen(fn, "w", stdout)) warning("Error redirecting stdout\n");
    snprintf(fn,256,"run_%d.log",pid);
    if(!freopen(fn, "w", stderr)) warning("Error redirecting stderr\n");
	
    mysymlink(fn, "run_recent.log");
    //turns off file buffering
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
    snprintf(fn,256,"kill_%d",pid);
    FILE *fp=fopen(fn,"w");
    fprintf(fp,"#!/bin/sh\nkill %d && rm $0 -rf \n", pid);
    fclose(fp);
    chmod(fn,00700);
}
