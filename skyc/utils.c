/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <getopt.h>
#include "skyc.h"
#include "parms.h"
#include "utils.h"
char *dirsetup;
/**
   Print out usage information.
 */
static void print_usage(void){
    info2(
	  "Usage: skyc [OPTION...] [FILE]...\n"
	  "skyc is a simulation tool developed to do sky coveragepostprocessing\n\n"
	  "Options: \n"
	  "-h, --help        to print out this message\n"
	  "-d, --detach      to detach from terminal and run in background\n"
	  "-f, --force       force starting simulation without scheduler\n"
	  "-n N, --nthread=N Use N threads, default is 1\n"
	  "-o DIR, --output=DIR\n"
	  "                  output the results to DIR.\n"
	  "-c FILE.conf, --conf=FILE.conf\n"
	  "                  Use FILE.conf as the baseline config instead of maos.conf\n"
	  );
    exit(0);
}
/**
   Parse command line arguments argc, argv
 */
ARG_S *parse_args(int argc, const char *argv[]){
    ARG_S *arg=calloc(1, sizeof(ARG_S));
    char *host=NULL; int local=0;
    ARGOPT_T options[]={
	{"help", 'h', T_INT, 2, print_usage, NULL},
	{"detach", 'd',T_INT, 0, &arg->detach, NULL},
	{"override",'O',T_INT,0, &arg->override, NULL},
	{"force",  'f',T_INT, 0, &arg->force, NULL},
	{"output", 'o',T_STR, 1, &arg->dirout, NULL},
	{"nthread",'n',T_INT, 1, &arg->nthread,NULL},
	{"conf",   'c',T_STR, 1, &arg->conf, NULL},
	{"path",   'P',T_STR, 3, addpath, NULL},
	{"run",    'r',T_STR, 1, &host, NULL},
	{"local",  'l',T_INT, 0, &local, NULL},
	{NULL, 0,0,0, NULL, NULL}
    };
    char *cmds=strnadd(argc-1, argv+1, " ");
    parse_argopt(cmds, options);
    if(!host && !arg->detach){//foreground running
	arg->force=1;
    }else if(local || getenv("MAOS_DIRECT_LAUNCH")){
	/*lanched through scheduler to run locally. We are already detached, so
	  don't daemonize again.*/
	arg->detach=0;
	arg->force=0;
	detached=1;
    }else{
#ifndef MAOS_DISABLE_SCHEDULER
	/*Detached version. Always launch through scheduler if available.*/
	if(!host){
	    host=strdup(myhostname());
	}
	if(scheduler_launch_exe(host, argc, argv)){
	    error2("Unable to launch skyc at server %s\n", host);
	}else{
	    exit(EXIT_SUCCESS);
	}
#endif
    }
    free(host);
    /*do not use NTHREAD here. causes too much cpu load*/
    if(arg->nthread>NCPU || arg->nthread<=0){
	arg->nthread=NCPU;
    }else{
        NTHREAD=arg->nthread;
    }
    char fntmp[PATH_MAX];
    snprintf(fntmp,PATH_MAX,"%s/skyc_%ld.conf",TEMP,(long)getpid());
    FILE *fptmp=fopen(fntmp,"w");
    if(cmds){
	fputs(cmds, fptmp);
	fclose(fptmp);
	free(cmds); cmds=NULL;
    }
    arg->confcmd=strdup(fntmp);
    if(!arg->dirout){
	arg->dirout=strtime();
    }
    if(!arg->conf){ /*If -c is not specifid in path, will use maos.conf*/
	arg->conf=strdup("maos.conf");
    }
    addpath(".");
    mymkdir("%s",arg->dirout);
    if(chdir(arg->dirout)){
	error("Unable to chdir to %s\n", arg->dirout);
    }
    return arg;
}
/**
   Rename the log files when simulation exits.
 */
void rename_file(int sig){
    if(sig==0){
	rename("run_recent.log", "run_done.log");
	mysymlink("run_done.log", "run_recent.log");
	rename("skyc_recent.conf", "skyc_done.conf");
	mysymlink("skyc_done.conf", "skyc_recent.conf");
    }
}
/**
   Handles signals. 
 */
int skyc_signal_handler(int sig){
    psignal(sig, "skyc");
    rename_file(sig);/*handles signal */
    scheduler_finish(sig);
    return 0;
}
