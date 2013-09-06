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
GLOBAL_T *global=NULL;//record for convenient access.
int use_cuda=0;
const char *dirsetup=NULL;
const char *dirskysim=NULL;
volatile int maos_server_fd=-1;
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
	powfs = setup_powfs_init(parms, aper);
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
	recon = setup_recon_init(parms);
	/*Setup DM fitting parameters so we can flatten the DM*/
	setup_recon_dm(recon, parms);
	/*setting up M1/M2/M3, Instrument, Lenslet surface OPD. DM Calibration, WFS bias.*/
	setup_surf(parms, aper, powfs, recon);
	/*set up physical optics powfs data*/
	setup_powfs_phy(parms, powfs);
    }

    if(!parms->sim.evlol){
	setup_recon(recon, parms, powfs, aper);
	info2("After setup_recon:\t%.2f MiB\n",get_job_mem()/1024.);
    }
    global->powfs=powfs;
    global->aper=aper;
    global->recon=recon;
    global->setupdone=1;
    if(parms->plot.setup){
	plot_setup(parms, powfs, aper, recon);
    }
#if USE_CUDA
    extern int cuda_dedup;
    cuda_dedup=1;
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
    cuda_dedup=0;
#endif

    if(!parms->sim.evlol && parms->recon.mvm){
	setup_recon_mvm(parms, recon, powfs);
    }
    /*
      Before entering real simulation, make sure to delete all variables that
      won't be used later on to save memory.
    */

    free_recon_unused(parms, recon);
    toc2("Presimulation");sync();
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
static void maos_server(PARMS_T *parms){
    if(maos_server_fd<0){
	error("Invalid maos_server_fd\n");
	EXIT;
    }
    warning("maos running in server mode\n");
    int msglen;
    int sock=maos_server_fd;
    while(!streadint(sock, &msglen)){
	int ret=0;/*acknowledgement of the command. 0: accepted. otherwise: not understood*/
	int *msg=alloca(sizeof(int)*msglen);
	
	if(streadintarr(sock, msg, msglen)){
	    break;
	}
	switch(msg[0]){
	case MAOS_ASSIGN_WFS:{/*Specifies which WFS to be handled*/
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		parms->wfs[iwfs].sock=msg[iwfs+1]?-1:0;
	    }
	}break;
	case MAOS_ASSIGN_EVL:{/*Specifies which EVL to be handled*/
	    if(!parms->evl.sock){
		parms->evl.sock=calloc(parms->evl.nevl, sizeof(int));
	    }
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		parms->evl.sock[ievl]=msg[ievl+1]?-1:0;
	    }
	}break;
	case MAOS_ASSIGN_RECON:{/*Specifies whether recon should be handled*/
	    parms->recon.sock=msg[1]?-1:0;
	}break;
	case MAOS_ASSIGN_DONE:{/*Now start configuration*/
	    error("Not completed\n");
	}break;
	default:
	    ret=1;
	}/*switch*/

	if(stwriteint(sock, ret)){
	    break;
	}
    }
    warning("maos_server exited\n");
    //todo:listening on socket for commands.
    //maos client mode: start maos server mode via scheduler.
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
    /*Launch the scheduler and report about our process */
    scheduler_start(scmd,arg->nthread,!arg->force);
    info2("%s\n", scmd);
    info2("Output folder is '%s'. %d threads\n",arg->dirout, arg->nthread);
    maos_version();
    /*setting up parameters before asking scheduler to check for any errors. */
    PARMS_T *parms=setup_parms(arg);
    global=calloc(1, sizeof(GLOBAL_T));
    global->parms=parms;
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
    if(!disable_save){
	//Make the symlinks for running job only.
	char fnpid[PATH_MAX];
	snprintf(fnpid, PATH_MAX, "maos_%d.conf", (int)getpid());
	mysymlink(fnpid, "maos_recent.conf");
	snprintf(fnpid, PATH_MAX, "run_%d.log", (int)getpid());
	mysymlink(fnpid, "run_recent.log");
    }
    info2("\n*** Simulation started at %s in %s. ***\n\n",myasctime(),myhostname());
    thread_new((thread_fun)scheduler_listen, maos_daemon);
    THREAD_POOL_INIT(parms->sim.nthread);
    setup_parms_running(parms, arg);
    if(arg->server){
	while(maos_server_fd<0){
	    warning("Waiting for fd\n");
	    sleep(1);
	}
	maos_server(parms);
	EXIT;
    }
    free(scmd);
    free(arg->dirout);
    free(arg->gpus);
    free(arg);
    if(parms->save.setup){
	dirsetup="setup";
	mymkdir("%s",dirsetup);
    }else{
	dirsetup=".";
    }
    if(parms->sim.skysim){
	dirskysim="skysim";
	mymkdir("%s",dirskysim);
    }else{
	dirskysim=".";
    }

    /*Loads the main software*/
#if _OPENMP>=200805
#pragma omp parallel
#pragma omp single nowait
#endif
    maos(parms);
    free_parms(parms);
    info2("Job finished at %s\n",myasctime());
    rename_file(0);
    scheduler_finish(0);
    draw_final(0);
    exit_success=1;/*tell mem.c to print non-freed memory in debug mode. */
    return 0;
}
