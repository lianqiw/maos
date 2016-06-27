/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "common.h"
#include "sim_utils.h"
#include "maos.h"
#include "version.h"
static void maos_daemon(int sock){
    thread_block_signal();
    int cmd[2];
    while(!streadintarr(sock, cmd, 2)){
	switch(cmd[0]){
	case MAOS_DRAW://Starting draw to received fd.
	    {
		int fd;
		if(streadfd(sock, &fd)){
		    warning("unable to read fd from %d\n", sock);
		    continue;
		}else{
		    warning("got fd=%d\n", fd);
		}
		draw_add(fd);
		PARMS_T *parms=(PARMS_T*)global->parms;//cast away constness
		parms->plot.setup=1;
		parms->plot.run=1;
		if(global->setupdone){//setup is already finished. request plot setup.
		    plot_setup(global->parms, global->powfs, global->aper, global->recon);
		}
	    }break;
	default:
	    warning("unknown cmd %d\n", cmd[0]);
	    break;
	}
    }
}
void maos_version(void){
    info2("MAOS Version %s. Compiled on %s %s by %s, %d bit", PACKAGE_VERSION, __DATE__, __TIME__, __VERSION__, (int)sizeof(long)*8);
#if USE_CUDA
#if CUDA_DOUBLE
    info2(", w/t CUDA(double)");
#else
    info2(", w/t CUDA(single)");
#endif
#else
    info2(", w/o CUDA");
#endif
#ifdef __OPTIMIZE__
    info2(", w/t optimization.\n");
#else
    info2(", w/o optimization\n");
#endif
    info2("Source: %s %s\n", SRCDIR, GIT_VERSION);
    info2("BUILD: %s\n", BUILDDIR);
    info2("Launched at %s in %s with PID %ld.\n",myasctime(),HOST, (long)getpid());
#if HAS_LWS
    extern uint16_t PORT;
    info2("The web based job monitor can be accessed at http://localhost:%d\n", 1+PORT);
#endif
}

/**
   This is the standard entrance routine to the program.  It first calls
   setup_parms() to setup the simulation parameters and check for possible
   errors. It then waits for starting signal from the scheduler if in batch
   mode. Finally it hands the control to maos() to start the actual simulation.

   Call maos with overriding *.conf files or embed the overriding parameters in
   the command line to override the default parameters, e.g.

   <p><code>maos base.conf save.setup=1 'powfs.phystep=[0 100 100]'</code><p>

   Any duplicate parameters will override the pervious specified value. The
   configure file nfiraos.conf will be loaded as the master .conf unless a -c
   switch is used with another .conf file. For scao simulations, call maos with
   -c switch and the right base .conf file.
   
   <p><code>maos -c scao_ngs.conf override.conf</code><p>

   for scao NGS simulations 

   <p><code>maos -c scao_lgs.conf override.conf</code><p>

   for scao LGS simulations.  With -c switch, nfiraos.conf will not be read,
   instead scao_ngs.conf or scao_lgs.conf are read as the master config file.
   Do not specify any parameter that are not understood by the code, otherwise
   maos will complain and exit to prevent accidental mistakes.
       
   Generally you link the maos executable into a folder that is in your PATH
   evironment or into the folder where you run simulations.

   Other optional parameters:
   \verbatim
   -d          do detach from console and not exit when logged out
   -s 2 -s 4   set seeds to [2 4]
   -n 4        launch 4 threads.
   -f          To disable job scheduler and force proceed
   \endverbatim
   In detached mode, drawing is automatically disabled.
   \callgraph
*/
int main(int argc, const char *argv[]){
    char *scmd=argv2str(argc,argv," ");
    ARG_T* arg=parse_args(argc,argv);/*does chdir */
    if(arg->detach){
	daemonize();
    }else{
	redirect();
    }
    /*Launch the scheduler if it is not running and report about our process */
    int ngpu;
#if USE_CUDA
    ngpu=arg->ngpu;
    if(!ngpu) ngpu=0xFFFFFF;
#else
    ngpu=0;
#endif
    scheduler_start(scmd,NTHREAD,ngpu,!arg->force);
    info2("%s\n", scmd);
    info2("Output folder is '%s'. %d threads\n",arg->dirout, NTHREAD);
    maos_version();
    /*setting up parameters before asking scheduler to check for any errors. */
    PARMS_T *parms=setup_parms(arg->conf, arg->confcmd, arg->override);
    free(arg->conf); arg->conf=0;
    if(arg->confcmd){
	remove(arg->confcmd); free(arg->confcmd); arg->confcmd=0;
    }
    info2("After setup_parms:\t %.2f MiB\n",get_job_mem()/1024.);
    /*register signal handler */
    register_signal_handler(maos_signal_handler);
 
    if(!arg->force){
	/*
	  Ask job scheduler for permission to proceed. If no CPUs are
	  available, will block until ones are available.
	  if arg->force==1, will run immediately.
	*/
	info2("Waiting start signal from the scheduler ...\n");
	int count=0;
	while(scheduler_wait()&& count<60){
	    /*Failed to wait. fall back to own checking.*/
	    warning_time("failed to get reply from scheduler. retry\n");
	    sleep(10);
	    count++;
	    scheduler_start(scmd,NTHREAD,ngpu,!arg->force);
	}
	if(count>=60){
	    warning_time("fall back to own checker\n");
	    wait_cpu(NTHREAD);
	}
    }
    thread_new((thread_fun)scheduler_listen, (void*)maos_daemon);
    setup_parms_gpu(parms, arg->gpus, arg->ngpu);
    free(scmd);
    free(arg->dirout);
    free(arg->gpus);
    free(arg);
    /* do not use prallel single in maos(). It causes blas to run single threaded
     * during preparation. Selective enable parallel for certain setup functions
     * that doesn't use blas*/
    maos(parms);
    rename_file(0);
    scheduler_finish(0);
    return 0;
}
