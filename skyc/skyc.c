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
#include <unistd.h>

#include "skyc.h"
#include "parms.h"
#include "utils.h"
#include "skysim.h"

/**
   \file skyc.c This file contains main routine about skycoverage.*/
char* dirstart;

/**
   The main(). It parses the command line, setup the parms, ask the scheduler
   for signal to proceed, and then starts skysim to do sky coverage.
 */
int main(int argc, const char* argv[]){
	dirstart=mygetcwd();
	char* scmd=argv2str(argc, argv, " ");
	arg_s* arg=parse_args(argc, argv);
	/*In detach mode send to background and disable drawing*/
	if(arg->detach){
		daemonize();
	} else{
		redirect();
	}
	info("%s\n", scmd);
	info("Output folder is '%s'. %d threads\n", arg->dirout, arg->nthread);
	print_version();
	/*register signal handler */
	register_signal_handler(skyc_signal_handler);
	/*
	  Ask job scheduler for permission to proceed. If no CPUs are available,
	  will block until ones are available.  if arg->force==1, will run
	  immediately.
	*/
	scheduler_report_path(scmd);
	/*setting up parameters before asking scheduler to check for any errors. */
	dirsetup=stradd("setup", NULL);
	parms_s* parms=setup_parms(arg);
	if(parms->maos.nseed){
		if(parms->skyc.dbg){
			mymkdir("%s", dirsetup);
		}
		scheduler_start(arg->nthread, 0, !arg->force);
		info("Simulation started at %s in %s.\n", myasctime(0), HOST);
		THREAD_POOL_INIT(parms->skyc.nthread);
		/*Loads the main software*/
OMPTASK_SINGLE 
		skysim(parms);
	}
	free(scmd);
	free_arg(&arg);
	free_parms(parms);parms=NULL;
	free(dirsetup);
	free(dirstart);
	skyc_final(0);
	scheduler_finish(signal_caught);
	print_mem("End");
	info("Simulation finished at %s in %s.\n", myasctime(0), HOST);
	return 0;
}
