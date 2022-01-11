/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "version.h"
/**
   \file skyc.c This file contains main routine about skycoverage.*/
char* dirstart;

void skyc_version(void){
	info2("SRC: %s v%s %s\n", SRCDIR, PACKAGE_VERSION, GIT_VERSION);
	char exe[PATH_MAX];
	if(!get_job_progname(exe, PATH_MAX, 0)){
		info2("BUILT: %s by %s on %s", BUILDDIR, COMPILER, myasctime(fmtime(exe)));
	} else{
		info2("BUILT: %s by %s on %s %s", BUILDDIR, COMPILER, __DATE__, __TIME__);//__DATE__ and __TIME__ is only applicable to this specific file
	}
#if CPU_SINGLE 
	info2(" CPU(single)");
#else
	info2(" CPU(double)");
#endif
#if USE_CUDA
#if CUDA_DOUBLE
	info2(" +CUDA(real)");
#else
	info2(" +CUDA(single)");
#endif
#else
	info2(" -CUDA");
#endif
#ifdef __OPTIMIZE__
	info2(" +optimization.\n");
#else
	info2(" -optimization\n");
#endif
	info2("Launched at %s in %s with PID %ld.\n", myasctime(0), HOST, (long)getpid());
#if HAS_LWS
	extern uint16_t PORT;
	info2("The web based job monitor can be accessed at http://localhost:%d\n", 100+PORT);
#endif
}
/**
   The main(). It parses the command line, setup the parms, ask the scheduler
   for signal to proceed, and then starts skysim to do sky coverage.
 */
int main(int argc, const char* argv[]){
	dirstart=mygetcwd();
	char* scmd=argv2str(argc, argv, " ");
	ARG_S* arg=parse_args(argc, argv);
	/*In detach mode send to background and disable drawing*/
	if(arg->detach){
		daemonize();
	} else{
		redirect();
	}
	info("%s\n", scmd);
	info("Output folder is '%s'. %d threads\n", arg->dirout, arg->nthread);
	skyc_version();
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
	PARMS_S* parms=setup_parms(arg);
	if(parms->maos.nseed){
		if(parms->skyc.dbg){
			mymkdir("%s", dirsetup);
		}
		scheduler_start(arg->nthread, 0, !arg->force);
		info("Simulation started at %s in %s.\n", myasctime(0), HOST);
		THREAD_POOL_INIT(parms->skyc.nthread);
		/*Loads the main software*/
		OMPTASK_SINGLE skysim(parms);
	}
	free(scmd);
	free(arg->conf); 
	free(arg->confcmd);
	free(arg->dirout); 
	free(arg);
	free_parms(parms);parms=NULL;
	free(dirsetup);
	free(dirstart);
	rename_file(0);
	scheduler_finish(signal_caught);
	print_mem("End");
	info("Simulation finished at %s in %s.\n", myasctime(0), HOST);
	return 0;
}
