/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file load.c

   This routine is a matlab script helper. It launches matlab scripts and
   supervise it. The job is controlled by scheduler.
 */
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <limits.h>
#include "../sys/sys.h"
long matlab_pid;
static int kill_matlab(int sig){
    if(sig!=0){
	kill(matlab_pid,SIGTERM);
	sleep(10);
	if(!kill(matlab_pid,0)){
	    kill(matlab_pid,SIGKILL);
	}
    }
    return 0;
}
int main(int argc, char **argv){
    if(argc!=2){
	error("Usage:load script.m or load script\n");
    }
    char *cwd=mygetcwd();
    int slen=strlen(cwd)+3+strlen(argv[1]);
    char scmd[slen];
    if(!mystrcmp(cwd,HOME)){
	strcpy(scmd,"~");
	strcat(scmd,cwd+strlen(HOME));
    }else{
	strcpy(scmd,cwd);
    }
    free(cwd);
    strcat(scmd,"/");
    strcat(scmd,argv[1]);
    if(!strcmp(argv[1]+strlen(argv[1])-2,".m")){
	argv[1][strlen(argv[1])-2]='\0';
    }
    info("Will launch %s\n",scmd);
    long pid=fork();
    if(pid<0){
	exit(EXIT_FAILURE);
    }else if(pid>0){
	waitpid(pid,NULL,0);
	exit(EXIT_SUCCESS);
    }
    if(setsid()==-1) error("Error setsid\n");
    umask(0077);
    pid=fork();
    if(pid<0){
	exit(EXIT_FAILURE);
    }else if(pid>0){
	exit(EXIT_SUCCESS);
    }

    /*Now we are in the background. */
    /*dbg("Waiting for the scheduelr\n"); */
	scheduler_report_path(scmd);
    scheduler_start(1,0,1);
    /*dbg("Ready to go\n"); */
    
    pid=fork();
    if(pid<0){
	exit(EXIT_FAILURE);
    }else if(pid>0){
	matlab_pid=pid;
	register_signal_handler(kill_matlab);
	int status;
	double starttime=myclocki();
	waitpid(pid,&status,WUNTRACED | WCONTINUED);
	dbg("status=%d\n",status);
	status_t *st=mycalloc(1,status_t);
	st->laps=0;
	st->rest=myclocki()-starttime;
	st->iseed=0;
	st->nseed=1;
	st->isim=-1;
	st->simend=1;
	if(WIFEXITED(status)){
	    st->info=S_FINISH;
	    st->isim=0;
	    st->laps=myclocki()-starttime;
	    st->rest=0;
	}else if(WIFSIGNALED(status) || WIFSTOPPED(status)){
	    st->info=S_KILLED;	
	}else{
	    st->info=S_CRASH;
	}
	scheduler_report(st);
	free(st);
	/*scheduler_finish(status); */
	exit(EXIT_SUCCESS);
    }else{
	char fnlog[PATH_MAX];
	snprintf(fnlog,PATH_MAX,"load_%s_%d.log",argv[1],(int)getpid());
	if(!freopen("/dev/null", "r", stdin)) warning("Error redirect stdin\n");
	if(!freopen(fnlog, "w", stdout)) warning("Error redirect stdout\n");
	if(!dup2(fileno(stdout), fileno(stderr))) warning("Error redirect stderr\n");
	setbuf(stdout,NULL);/*disable buffering. */
	unsetenv("DISPLAY");
	execlp("matlab","matlab","-nojvm","-nodesktop","-nosplash","-r",argv[1],NULL);
    }
}
