/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "skyc.h"
#include "parms.h"
#include "utils.h"
#include "skysim.h"
#include "version.h"
/**
   \file skyc.c This file contains main routine about skycoverage.*/
char *dirstart;

void skyc_version(void){
    info("SRC: %s v%s %s\n", SRCDIR, PACKAGE_VERSION, GIT_VERSION);
    info("BUILD: %s by %s on %s %s", BUILDDIR, COMPILER, __DATE__, __TIME__);
#if USE_CUDA
#if CUDA_DOUBLE
    info(" +CUDA(double)");
#else
    info(" +CUDA(single)");
#endif
#else
    info(" -CUDA");
#endif
#ifdef __OPTIMIZE__
    info(" +optimization.\n");
#else
    info(" -optimization\n");
#endif
    info("Launched at %s in %s with PID %ld.\n",myasctime(),HOST, (long)getpid());
#if HAS_LWS
    extern uint16_t PORT;
    info("The web based job monitor can be accessed at http://localhost:%d\n", 1+PORT);
#endif
}
/**
   The main(). It parses the command line, setup the parms, ask the scheduler
   for signal to proceed, and then starts skysim to do sky coverage.
 */
int main(int argc, const char *argv[]){
    dirstart=mygetcwd();
    char *scmd=argv2str(argc, argv, " ");
    ARG_S* arg=parse_args(argc,argv);
    /*In detach mode send to background and disable drawing*/
    if(arg->detach){
	daemonize();
    }else{
	redirect();
    }
    info("%s\n", scmd);
    info("Output folder is '%s'. %d threads\n",arg->dirout, arg->nthread);
    skyc_version();
    /*register signal handler */
    register_signal_handler(skyc_signal_handler);
    /*
      Ask job scheduler for permission to proceed. If no CPUs are available,
      will block until ones are available.  if arg->force==1, will run
      immediately.
    */
    scheduler_start(scmd,arg->nthread,0,!arg->force);
    /*setting up parameters before asking scheduler to check for any errors. */
    dirsetup=stradd("setup",NULL);
    PARMS_S * parms=setup_parms(arg);
    if(parms->skyc.dbg){
	mymkdir("%s",dirsetup);
    }
    if(!arg->force){
	info("Waiting start signal from the scheduler ...\n");
	/*Failed to wait. fall back to own checking.*/
	int count=0;
	while(scheduler_wait()&& count<60){
	    warning_time("failed to get reply from scheduler. retry\n");
	    sleep(10);
	    count++;
	    scheduler_start(scmd,arg->nthread,0,!arg->force);
	}
	if(count>=60){
	    warning_time("fall back to own checker\n");
	    wait_cpu(arg->nthread);
	}
    }
    info("Simulation started at %s in %s.\n",myasctime(),HOST);
    free(scmd);
    free(arg->dirout);
    free(arg);
    THREAD_POOL_INIT(parms->skyc.nthread);
    /*Loads the main software*/
    OMPTASK_SINGLE skysim(parms);
    free_parms(parms);
    free(dirsetup);
    free(dirstart);
    rename_file(0);
    scheduler_finish(0);
    info("End:\t%.2f MiB\n",get_job_mem()/1024.);
    info("Simulation finished at %s in %s.\n",myasctime(),HOST);
    return 0;
}
