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
#include <unistd.h>
#include <getopt.h>
#include "skyc.h"
#include "parms.h"
#include "utils.h"
char* dirsetup;
/**
   Print out usage information.
 */
static void print_usage(void){
	info(
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
void free_arg(arg_s** parg){
	if(!parg) return;
	arg_s *arg=*parg;
	free(arg->conf); 
	free(arg->confcmd);
	free(arg->dirout); 
	free(arg);
}
/**
   Parse command line arguments argc, argv
 */
arg_s* parse_args(int argc, const char* argv[]){
	arg_s* arg=mycalloc(1, arg_s);
	char* host=NULL; int local=0;
	argopt_t options[]={
	{"help",  'h', M_INT, 0, 1, (void*)print_usage, NULL},
	{"detach", 'd',M_INT, 0, 0, &arg->detach, NULL},
	{"override",'O',M_INT,0, 0, &arg->override, NULL},
	{"force",  'f',M_INT, 0, 0, &arg->force, NULL},
	{"output", 'o',M_STR, 1, 0, &arg->dirout, NULL},
	{"nthread",'n',M_INT, 1, 0, &arg->nthread,NULL},
	{"conf",   'c',M_STR, 1, 0, &arg->conf, NULL},
	{"path",   'P',M_STR, 1, 1, (void*)addpath, NULL},
	{"run",    'r',M_STR, 1, 0, &host, NULL},
	{"local",  'l',M_INT, 0, 0, &local, NULL},
	{NULL, 0,0,0,0, NULL, NULL}
	};
	arg->confcmd=strnadd(argc-1, argv+1, " ");
	parse_argopt(arg->confcmd, options);
	if(!host&&!arg->detach){//foreground running
		arg->force=1;
	} else if(local||getenv("MAOS_DIRECT_LAUNCH")){
	/*lanched through scheduler to run locally. We are already detached, so
	  don't daemonize again.*/
		arg->detach=0;
		arg->force=0;
		detached=1;
	} else{
#ifndef MAOS_DISABLE_SCHEDULER
	/*Detached version. Always launch through scheduler if available.*/
		if(scheduler_launch_exe(host, argc, argv)<0){
			warning("Launch with scheduler failed. Restart without scheduler.\n");
		} else{
			free_arg(&arg);
			dbg3("Launched using scheduler.\n");
			exit(EXIT_SUCCESS);
		}
#endif
	}
	free(host);
	/*do not use NTHREAD here. causes too much cpu load*/
	if(arg->nthread>NCPU||arg->nthread<=0){
		arg->nthread=NCPU;
	} else{
		NTHREAD=arg->nthread;
	}
	if(!arg->dirout){
		arg->dirout=strtime_pid();
	}
	addpath2(2, ".");
	mymkdir("%s", arg->dirout);
	if(chdir(arg->dirout)){
		error("Unable to chdir to %s\n", arg->dirout);
	}
	return arg;
}
/**
   Rename the log files when simulation exits.
 */
void skyc_final(int sig){
	rename_log(sig, "skyc");
}
/**
   Handles signals.
 */
int skyc_signal_handler(int sig){
	info("skyc: %s", strsignal(sig));
	skyc_final(sig);/*handles signal */
	scheduler_finish(sig);
	return 0;
}
