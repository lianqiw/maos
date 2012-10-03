/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "maos.h"
#include "setup_powfs.h"
#include "setup_recon.h"
#include "setup_aper.h"
#include "sim.h"
#include "sim_utils.h"
#include "setup_surf.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
/**
   \file maos.c
   Contains main() and the entry into simulation maos()
*/
const PARMS_T *curparms=NULL;
int curiseed=0;
int use_cuda=0;
/**
   This is the routine that calls various functions to do the simulation. maos()
   calls setup_aper(), setup_powfs(), and setup_recon() to set up the aperture
   (of type APER_T), wfs (of type POWFS_T), and reconstructor (of type RECON_T)
   structs and then hands control to sim(), which then stars the simulation.
   \callgraph */
void maos(const PARMS_T *parms){    
    TIC;tic;
    APER_T  * aper=NULL;
    POWFS_T * powfs=NULL;
    RECON_T * recon=NULL;
    aper  = setup_aper(parms);
    info2("After setup_aper:\t%.2f MiB\n",get_job_mem()/1024.);
    if(!parms->sim.evlol){
	powfs = setup_powfs(parms, aper);
	info2("After setup_powfs:\t%.2f MiB\n",get_job_mem()/1024.);
	if(parms->dbg.wfslinearity!=-1){
	    int iwfs=parms->dbg.wfslinearity;
	    assert(iwfs>-1 || iwfs<parms->nwfs);
	    info2("Studying wfslineariy for WFS %d\n", iwfs);
	    wfslinearity(parms, powfs, iwfs);
	    rename_file(0);
	    scheduler_finish(0);
	    exit_success=1;/*tell mem.c to print non-freed memory in debug mode. */
	    exit(0);
	}
	recon = setup_recon(parms, powfs, aper);
	info2("After setup_recon:\t%.2f MiB\n",get_job_mem()/1024.);
    }
    setup_surf(parms, aper, powfs, recon);/*setting up M1/M2/M3 surface OPD. */
    
#if USE_CUDA
    if(!parms->sim.evlol && (parms->gpu.wfs || parms->gpu.tomo)){
	gpu_wfsgrad_init(parms, powfs);
    }
    if(parms->gpu.evl){
	gpu_perfevl_init(parms, aper);
    }
    if(parms->gpu.wfs && powfs){
	gpu_wfssurf2gpu(parms, powfs);
    }
    if(parms->gpu.evl){
	gpu_evlsurf2gpu(aper);
    }
    if(!parms->sim.evlol && (parms->gpu.tomo || parms->gpu.fit)){
	gpu_setup_recon(parms, powfs, recon);
    }
#endif

    if(!parms->sim.evlol && parms->recon.mvm){
	setup_recon_mvm(parms, recon, powfs);
    }
    /*
      Before entering real simulation, make sure to delete all variables that
      won't be used later on to save memory.
    */
#if USE_MKL==1
    int one=1;
    omp_set_num_threads(&one);/*only allow 1 thread after svd/chol is done. */
#endif
    free_recon_unused(parms, recon);
    toc2("Presimulation");
    sim(parms, powfs, aper, recon);
    /*Free all allocated memory in setup_* functions. So that we
      keep track of all the memory allocation.*/
    free_recon(parms, recon); recon=NULL;
    free_powfs(parms, powfs); powfs=NULL;
    free_aper(aper, parms); aper=NULL;
#if USE_CUDA
    if(use_cuda){
	gpu_cleanup();
    }
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
int main(int argc, char **argv){
    char *fn=mybasename(argv[0]);
    strcpy(argv[0], fn);
    free(fn);
    char *scmd=argv2str(argc,argv);
    ARG_T* arg=parse_args(argc,argv);/*does chdir */

    if(arg->detach){
	daemonize();
    }else{
	redirect();
	/*foreground task will start immediately. */
	arg->force=1;
    }
    /*Launch the scheduler and report about our process */
    scheduler_start(scmd,arg->nthread,!arg->force);
    info2("%s\n", scmd);
    info2("MAOS Version %s. Compiled on %s %s by %s, %d bit", PACKAGE_VERSION, __DATE__, __TIME__, __VERSION__, (int)sizeof(long)*8);
#if USE_CUDA
    info2(", w/t CUDA");
#else
    info2(", w/o CUDA");
#endif
#ifdef __OPTIMIZE__
    info2(", w/t optimization.\n");
#else
    info2(", w/o optimization\n");
#endif
    info2("Source: %s\n", SRCDIR);
    info2("Launched at %s in %s with %d threads\n",myasctime(),myhostname(),arg->nthread);
    info2("Output folder is '%s'.\n",arg->dirout);
    /*setting up parameters before asking scheduler to check for any errors. */
    PARMS_T *parms=setup_parms(arg);
    
#if USE_MKL==1
    if(arg->nthread==1){
	int one=1;
	/*only allow 1 thread to make debug easy. */
	omp_set_num_threads(&one);
    }
#endif

    curparms = parms;
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

    info2("\n*** Simulation started at %s in %s. ***\n\n",myasctime(),myhostname());
    setup_parms_running(parms, arg);
    free(scmd);
    free(arg->dirout);
    free(arg->gpus);
    free(arg);
    THREAD_POOL_INIT(parms->sim.nthread);
    dirsetup=stradd("setup",NULL);
    if(parms->save.setup || parms->save.recon){
	mymkdir("%s",dirsetup);
    }
    if(parms->sim.skysim){
	dirskysim=stradd("skysim",NULL);
	mymkdir("%s",dirskysim);
    }else{
	dirskysim=strdup(".");
    }

    /*Loads the main software*/
    maos(parms);
    free_parms(parms);
    free(dirsetup);
    free(dirskysim);
    info2("Job finished at %s\n",myasctime());
    rename_file(0);
    scheduler_finish(0);
    exit_success=1;/*tell mem.c to print non-freed memory in debug mode. */
    return 0;
}
