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
#include "skyc.h"
#include "parms.h"
#include "utils.h"
#include "skysim.h"
/**
   \file skyc.c This file contains main routine about skycoverage.*/
char *dirstart;
/**
   The main(). It parses the command line, setup the parms, ask the scheduler
   for signal to proceed, and then starts skysim to do sky coverage.
 */
int main(int argc, char **argv){
    dirstart=mygetcwd();
    char *scmd=argv2str(argc, argv);
    //Strip out PATH information from the command line.
    char*fn=mybasename(argv[0]);
    strcpy(argv[0],fn);
    free(fn);
    
    ARG_S* arg=parse_args(argc,argv);
    /*In detach mode send to background and disable drawing*/
    if(arg->detach){
	disable_draw=1;//disable drawing.
	daemonize();
    }
    info2("Output folder is '%s' %d threads\n",arg->dirout, arg->nthread);

#ifdef SVN_REV
    if(strlen(SVN_REV)>1 && strcmp(SVN_REV,"exported")){
	info2("MAOS Version %s, Revision %s,",PACKAGE_VERSION,SVN_REV);
    }
#else
    info2("MAOS Version %s,",PACKAGE_VERSION);
#endif
    info2("Launched at %s in %s.\n",myasctime(),myhostname());

    fprintf(stderr, "Compiled on %s %s by %s ",
	    __DATE__, __TIME__, __VERSION__);
#ifdef __OPTIMIZE__
    fprintf(stderr, "with optimization.\n");
#else
    fprintf(stderr, "without optimization!!!\n");
#endif
    //register signal handler
    register_signal_handler(skyc_signal_handler);
    /*
      Ask job scheduler for permission to proceed. If no CPUs are available,
      will block until ones are available.  if arg->force==1, will run
      immediately.
    */
    scheduler_start(scmd,arg->nthread,!arg->force);
    //setting up parameters before asking scheduler to check for any errors.
    dirsetup=stradd("setup",NULL);
    PARMS_S * parms=setup_parms(arg);
    if(parms->skyc.dbg){
	mymkdir("%s",dirsetup);
    }
    if(!arg->force){
	info2("Waiting start signal from the scheduler ...\n");
	/*Failed to wait. fall back to own checking.*/
	int count=0;
	while(scheduler_wait()&& count<60){
	    warning3("failed to get reply from scheduler. retry\n");
	    sleep(10);
	    count++;
	    scheduler_start(scmd,arg->nthread,!arg->force);
	}
	if(count>=60){
	    warning3("fall back to own checker\n");
	    wait_cpu(arg->nthread);
	}
    }
    info2("Simulation started at %s in %s.\n",myasctime(),myhostname());
 
    free(arg);
    free(scmd);
    THREAD_POOL_INIT(parms->skyc.nthread);
    /*Loads the main software*/
    skysim(parms);
    free_parms(parms);
    free(dirsetup);
    free(dirstart);
    rename_file(0);
    scheduler_finish(0);
    exit_success=1;
    info2("End:\t%.2f MiB\n",get_job_mem()/1024.);
    info2("Simulation finished at %s in %s.\n",myasctime(),myhostname());
    return 0;
}
