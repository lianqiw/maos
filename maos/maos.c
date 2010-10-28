/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include "setup_powfs.h"
#include "setup_recon.h"
#include "setup_aper.h"
#include "sim.h"
#include "sim_utils.h"
/**
   \file maos.c
   Contains main() and the entry into simulation maos()
*/
/**
    \page page1 Introduction

    \section sect-intro Introduction

    This code is a complete new implementation of the MCAO simulator LAOS, which
    was written in MATLAB language. The motivation to develop this software is
    to create a multi-conjugate adaptive optics simulator that runs fast,
    consumes less memory, and does the job without MATLAB, which is proprietary
    and has a large memory footprint. The code is written in C language with a
    function oriented design. The code is completely configurable through
    configuration files. The code tries its best to check the configuration for
    any apparent errors or conflicts. The configuration files are very readable
    and easy to maintain. The final development goal is to make an efficient,
    easy to use, general purpose adaptive optics simulator to help the
    development of adaptive optics systems, particularly advanced, very high
    order systems for TMT and other ELTs.

    Atmospheric turbulence is represented as one or multiple discrete screen(s)
    at different heights that evolves according to frozen flow with given wind
    velocity. The resulting aberrations are sensed by one or multiple natural
    guide star (NGS) or laser guide star (LGS) Shack-Hartmann wavefront
    sensor(s) (WFS). The sensors can be simulated as idealized wavefront
    gradient sensors, best Zernike fit tilt sensors, or a physical optic WFS
    using user specified pixel characteristics and a matched filter estimation
    algorithm. The laser guide star wavefront sensing model includes the impact
    of guide star elongation for a specified sodium layer profile and (an
    optional) polar coordinate CCD. The tomographic wavefront reconstruction
    then estimates the turbulence in one or several different heights from the
    pseudo open loop gradients measured by the WFS, using a computationally
    efficient implementation of a minimum variance reconstruction algorithm as
    described in [Ellerbroek, 2002]. These reconstructed turbulence screens are
    then fit to the actuators on one or several deformable mirrors (DMs) to
    achieve best possible correction over a specified field of view
    (FoV). Performance evaluation is done in terms of RMS wavefront error,
    Strehl ratio, and point spread function (PSFs) at a few science objects in
    the target FoV, which might be different from the FoV used for the fitting
    step. A range of additional specified features are implemented, such as
    telescope wavefront aberrations and the “split tomography” control
    algorithm.

    This software is written in the C language (revision 99), with external
    dependent libraries of FFTW version 3 (www.fftw.org) and blas/lapack
    (http://www.netlib.org/lapack/). The code contains a local copy of the
    package Arpack (http://www.caam.rice.edu/software/ARPACK/, to find eigen-
    value/vectors of sparse matrices, written in f77, so a FORTRAN compiler is
    needed) and Cholmod (http://www.cise.ufl.edu/research/sparse/SuiteSparse/,
    to do Cholesky factorization, copied from suite-sparse). An optimized blas
    library such ATLAS, GOTO-BLAS, or Intel MKL is necessary to get good
    performance for Cholesky decompositions of the tomography and/or fitting
    matrix..

    A C99 compliant compiler is required to compile the code, such as GCC 4 or
    Intel C++ compiler. The code has been successfully compiled on 64 bit
    GNU/Linux using GNU GCC and the Intel C++ compiler ICC (in C mode), and in
    MAC OS X 10.5 using the GCC compiler. For machines with Intel CPUs, the code
    generally runs faster when compiled with ICC.

    This software also contains two optional executable (drawdaemon, monitor)
    for plotting and job monitoring (jobs can be monitored on several different
    machines). These two executable requires GTK+ and the Cairo development and
    run-time library.

    The multi-threading is achieved using Posix Threads (no MPI is
    implemented). This parallelism works fairly well in multi-core, shared
    memory machines. For example, in a recently purchased desktop with Intel
    Core i5 quad-core processor at 2.66 GHz, the software (configured in a
    NFIRAOS split tomography mode with physical optics LGS and geometric best
    fit Zernike tilt NGS), compiled with Intel C++ compiler version 10, can run
    at 10 seconds per simulation time step with 1 thread and 3 seconds per time
    step with 4 threads. The same setup in a dual Intel W5590 server runs 7
    seconds per simulation time step with 1 thread and 4 seconds per step with 2
    threads. The scaling with number of cores is fairly good. The memory
    consumption in the NFIRAOS baseline case is about 2 G.


    
    The code contains three main directories. The lib, maos, and skyc
    directory. The directory lib contains routines that are not specific to the
    simulations and can be adopted by other routines. The directory maos
    contains routines that are specific the adaptive optics simulations. The
    directory skyc contains routines for sky coverage postprocessing.


 */

/**
   \page page2 Requirements 

   \section sect-requirement Requirements

    - C99 compliant compiler: GCC4 or ICC.

    - FORTRAN compiler: gfortran or Intel IFC

    - FFTW version 3. Can usually be installed through the linux distribution
    official repository or using source from www.fftw.org
    
    - Optimized blas and lapack. Blas and lapack can usually be installed
    through the linux distribution official repository. For 64 bit linux
    machines, it is also possible to use a redistributed version of the Intel
    MKL by appending --enable-mkl when you run configure
    
    - (Optional) GTK+ and Cairo. For drawdaemon and monitor.
*/

/**
   \page page3 Run the code

   \section sect-run Run the code
    \subsection sect-prep Preparing the folders and compiling

    We recommend using three different folders to 1) store the source tree, 2)
    compile the code, and 3) run the simulation.

    The source tree is where the source is extracted. Obtain the software in the
    form of the compressed tar ball: maos_version.tar.gz. Extract it into some
    directory and name the directory src_dir. Example:

    \verbatim
    cd ~/work/programming
    tar zxvf maos_version.tar.gz 
    cd maos_version
    export src_dir=$(pwd)
    \endverbatim
    
    Next we create another folder, where we are going to compile the code. example:

    \verbatim
    cd ~/work/maos
    mkdir comp_optim && cd comp_optim
    $src_dir/configure
    make
    \endverbatim

    which will config the compiling folder with GCC compiler and default
    optimization flags and compile the code.
    
    To use intel C++ compiler, do the following

    \verbatim
    cd ~/work/maos
    mkdir comp_icc && cd comp_icc
    $src_dir/configure --enable-icc
    make
    \endverbatim
    
    To use the redistributed Intel MKL library (for fast cholesky decomposition)
    \verbatim
    $src_dir/configure --enable-icc --enable-mkl
    \endverbatim

    To enable debug mode
    \verbatim
    $src_dir/configure --enable-debug
    \endverbatim
    
    To specify your own compiler using gcc4.1 which is not the system default
    \verbatim
    CC=gcc4.1 $src_dir/configure
    \endverbatim
    
    The compiled executable is maos in the sub-folder “bin” of the compiling
    folder. You do not have to do "make install" to run the simulations.

    \subsection sect-monitor Using the job monitor
    
    When reasonably recent GTK+ and Cairo libraries are present in the system,
    two additional executives will be compiled, “drawdaemon” and “monitor”. The
    plotting utility “drawdaemon” is located in sub-folder “lib”. It is launched
    automatically when plotting commands are issued in the software. The
    monitoring utility “monitor” is located in sub-folder “bin”.  Link monitor to
    a folder that is in your PATH and type “monitor” in the command line to use
    it.

    In order to use monitor to monitor jobs in different machines, follow these steps:

    Create a folder .aos in your home folder, create two files in there. Put
    these files in all the machines you want to monitor.
    \verbatim
~/.aos/hosts  Each line contains a machine's hostname you want to monitor. 
              They should be the names returned by command "hostname", 
              and you should be able to log into those machines 
              with these names. Set up /etc/hosts properly to ensure this.
~/.aos/port   contains the port you want the scheduler to bind to 
              so that the monitor can connect to it. 
              Any non-previledged port numbers are fine (like 10000). 
              Different users should have different port number. 
	      The port has to be the same on all the machines. 
    \endverbatim

    \subsection sect-simu-setup Defining simulation setup files

    All user modifiable options are configured through these config files with
    suffix “.conf”. The supplied default setup files are stored in sub-folder
    “config/maos”. The config files are organized in a modular way. One .conf file
    can include another .conf file with include=”filename.conf”. Most of the
    .conf files in “config/maos” file are self-documented.

    All the options are in the form of key=value. The usually convention is: if
    the entry in parms called parms.wfs.thetax, the entry in .conf file will be
    in the form wfs.thetax=value. Array values are embraced by[], string are
    embraced by double quote “”. a line ended with \ will be concatenated with
    the next line (no space is allowed after \). Anything after comment string #
    will be ignored during reading. The config key=value can be specified
    multiple times for any key, in which case the latter value will override the
    previous value. The software will complain and exit if any key=value pair is
    not understood or required but not provided.

    The usual way to organize config files for a specific problem is to keep a
    baseline .conf file in sub-folder “config/maos” in the source folder that specify
    all the necessary key=value entries through optionally including other .conf
    files, and one or multiple overriding file (in the simulation folder, not
    inside “config/maos”) override the values you want to change from the baseline.

    In the “config/maos” sub-folder in the source folder, there are three baseline
    .conf files: 1) nfiraos.conf for NFIRAOS mcao simulations, 2) scao_lgs.conf
    for single conjugate LGS simulation, and 2) scao_ngs.conf for single
    conjugate NGS simulation. Each of this file include a few other .conf files
    by specifying include=”filename.conf”, for example, all the three baseline
    .conf files include 1) “sim.conf” to specify closed loop information such as
    loop gain, number of time steps, etc, 2) “recon.conf” to specify
    reconstruction parameters in closed loop, 3) “dbg.conf” to specify debugging
    parameters. The MCAO baseline “nfiraos.conf” includes “atm_mk13n50p.conf” to
    specify MK13N median seeing turbulence profile, while the two SCAO baselines
    include “atm_scao.conf” that specifies the turbulence information with only
    a single turbulence layer.

    In the MCAO case, you can put include=”atm_srd.conf” in your overriding
    .conf file to use the TMT SRD profile instead of the default MK13N median
    seeing turbulence profile.

    The .conf files in the “config/maos” folder are only modified during code
    development and should not be modified by the end user. Use overriding .conf
    files in simulation folders to modify the values.

    \subsection sect-run Run the simulation

    The simulation is usually done outside of the source and compilation folder
    by linking the maos executable into the simulation folder or link maos into
    a direction that is in the PATH. Example

    \verbatim
cd ~/work/maos
mkdir NFIRAOS && cd NFIRAOS
ln ~/work/aos/comp_optim/bin/maos . -s
./maos
    \endverbatim
    
    You can also link “maos” into your PATH and run “maos” instead of “./maos”

    \verbatim
    cd ~/local/bin #make sure this is in your PATH
    ln ~/work/aos/comp_optim/bin/maos . -s 
    #Make sure you use -s to make symbolic link
    \endverbatim

    Use the following command to get the help message:
    \verbatim
    maos -h  #will print the help message and exit.
    \endverbatim

    If maos is launched without any parameters, the baseline NFIRAOS
    configuration nfiraos.conf will be used. The valid command line arguments
    are as follows

    \verbatim
-c scao_lgs.conf   Change the baseline config from nfiraos.conf to scao_lgs.conf
-d     deteach from the console and run in background. You can then close
       the terminal and log out without killing the job
-n 4   use 4 threads. (default is 1)
-s 2   to override the seed (normally specified in sim.seed) to 2. 
       repeat to have multiple seeds (-s1 -s2 -s3)
-f     Force starting immediately. Otherwise only proceed when cpu is available
       as controlled by the built-in "scheduler"
a.conf Append a.conf to the config file list
-o ABC Put the results in folder ABC (create if not already exists). If not 
       specified the default output folder name will be composed of the data 
       and time string in the form of 2010-03-05-161539.
    \endverbatim

    Example:
    \verbatim
    maos -c scao_lgs.conf override.conf -n 4 -d -o test 
    \endverbatim

    Will launch the software using “scao_lgs.conf” (Single conjugate, LGS
    system) as the baseline and override some parameters by override.conf. It
    will launch 4 threads in simulation, detach from the console, and output the
    results into folder “test”.

*/
/**
   \page page4 Guidelines

   \section sect-guide Guidelines

   When modifying the code, please adhere to the following guidelines as
   strictly as possible.

   Memory allocation:

   - When declaring a struct, always initialize it to zero by calling calloc or
   memset, unless you have a good reason not to.

   - Use calloc instead of malloc to initialize the memory unless for large
   arrays that you will initialize immediately.

   Functions:
      
      
   - The input parameters that is also the output should be grouped in the
   beginning.

   - In any function, there should usually be at most 1 return statement unless
   NULL or error code are returned before matural.

   - One utility function should be handle only one major mask to maximize
   reusability.

   - Input arguments to utility functions (in lib folder) should be basic types
   to maximize resuability. Input arguments to simulation functions can be
   wrapped (like simu) to maximum readability.
   
   Others:
   
   - Avoid non constant static variables/pointers unless you have a good reason
   to do so. static variables are bad for multi-threading and hard to free
   after usage.

   - whenever modify something temporarily for debugging, make a warning. let
   it easily identifiable. not make a hidden bug.
   
   - Do not hard code adjustable parameters.
   
   - Always declare unchanged variables to constants. declare pointers as
   restrict if possible.

   - Do not include system headers in header files unless necessary. include
   thoese headers in .c file. 

 */
/**
   \page page5 Profiling MAOS 

   MAOS is a pthread enabled software, and in the current stage, the gprof based
   profiling does not generate useful information. Instead, we have to reply on
   system wide profilers, such as the oprofile.

   Prepare:

   <code>opcontrol --callgraph=16</code><p>
   <code>opcontrol --start </code><p>
   
   Then Run MAOS
   
   View profile:

   <code>opcontrol --stop</code><p>
   <code>opcontrol --dump</code><p>
   
   <code>opreport -cgf -l maos | gprof2dot.py -f oprofile | dot -Tpng -o output.png</code><p>

   or 

   <code>opgprof maos | gprof2dot.py | dot -Tpng -o output.png</code><p>

   gprof2dot.py can be found in the scripts folder.
*/

/**
   \page page 6 C Fundementals

   The user is required to have a good understanding of the C fundementals, like
   pointer usage, to study and improve the code. Here we list a few to remind
   the user.
   \section sect-array Arrays
   
   If you declare an array:

   <code>double a[1024];</code> the symbol actually contains the address of a
   1024x8byte memory region that is automatically allocated in the stack (malloc,
   calloc allocations memory in the heap). Suppost the address of the memory region of 0xFF00. Assignment 
   <code>a[0]=1;</code>puts 1 in memory started at 0xFF00
   <code>a[1]=1;</code>puts 1 in memory started at 0xFF08.

   \section sect-pointer Pointers A pointer is a variable that stores the
   address of a memory region.  <code>void *p;</code> is a minimal pointer that
   just stores the address and can not do any thing about the memory it points
   to. Its value can be assigned to any other non const vectr or from any vector
   (1d or 2d).

   However, a normal pointer like <code>double *p;</code> also contains
   information about the memory it points to. More specifically, the compiler
   know the length of each element so that it can compute the address of the elements, just like arrays. For example.
   <code>double *p;</code>if p contains value of 0xFF00.<p>
   <code>p[0]=1;</code>will store a double 1 into memory region at 0xFF00<p>
   <code>p[1]=1;</code>will store a double 1 into memory region at 0xFF08 because a double requires 8 bytes.<p>
   <code>p=p+1;</code>p will now contain 0xFF08, i.e., points to the next element, instead of 0xFF01<p>
   <code>p[0]=1;</code> is strictly equivalent to <code>*p=1;<code> *p means visit the memory pointed by p instead of p itself.
   <code>p[1]=1;</code> is strictly equivalent to <code>*(p+1)=1;<code><p>

   
   \subsection sect-1d-pointer 1D Pointers
   The simplist usage of pointer is to create and index an 1-D array:
   <code>double *p;</code>A pointer stores the address of some memory region<p>
   <code>p=calloc(1024, sizeof(double));</code>Allocate a memory block of 1024x8 byte and store its address in p.<p>
   <code>p[0]=1;</code>We assign 1 to the first number in the array<p>

   This is equivalent as using the arrays
   <code>double a[1024];</code>Declare an 1-D array.<p>
   <code>a[0]=1;</code>We assign 1 to the first number in the array<p>
   The symbol a is like a constant pointer. We can assign it's value to a constant pointer:
   <code>double *const p=a;</code><p>
   <code>p[0]=1;</code><p>

   \subsection sect-const Const pointers vs. pointers point to const memory
   Beware about a const pointer and pointer that points to constant memory regions:
   <code>const double *p;</code>p is a pointer that points to a memory region that can not be changed.<p>
   Thus <code>p[0]=1;</code> is illegal.<p>
   <code>double * const p;</code>p is a constant pointer that can not have its value changed.<p>
   Thus <code>p=a</code> is illegal, unless you assign at declaration :<code>double *const p=a</code><p>
   However, it is legal to modify the values in the memory p points to:
   <code>p[0]=1;</code> is perfectly ok.

   \subsection sect-2d-pointer 2D Pointers A 2D pointer is like a 1D pointer. It
   stores the address of a memory region. However, unlike the 1D pointer, the
   compiler knows some additional infos about the 2D pointer. The compiler knows
   the length of the first dimension so that is can compute the memory location
   given 2D index like a[2][3];
   <code>double (*p)[2]=calloc(3*2,sizeof(double))</code> can be used in the same way as <code>double a[3][2]</code>
   
   The pointer p in our dmat (or cmat) struct points to the allocated memory. It
   can be cast to a two-d vector to facilitize visiting of the elements. The
   cast is conveniently done by PDMAT(a,pa) for dmat a. This simply declares a
   2-d pointer and assign a->p to it: <code>double
   (pa*)[2]=(void*)a->p</code>. a->p is first cast to void * because it is a
   double * and can not be assigned to double (*)[2] directly. The array
   elements (ix,jy) can be accessed by pa[iy][ix]. Notice that ix changes fastest.
 */
/** 
    \mainpage MAOS Documentation

    Multi-Thread Adaptive Optics Simulator (MAOS) is a end-to-end adaptive
    optics system simulator. It has the capability to simulation many different
    type of adaptive optics systems, including conventional single conjugate AO
    (SCAO), multi-conjugate AO (MCAO), laser tomography AO (LTAO), multi-object
    AO (MOAO), and ground layer AO (GLAO).

    <p>

    - \ref page1 <p>

    - \ref page2 <p>

    - \ref page3 <p>

    - \ref page4 <p>
  
    - \ref page5 <p>

    \section sect-structure Code Structure    

    - main() is the entrance into the program. It calls
        - setup_parms(): to setup the parameters of types PARMS_T and check for errors.
        - maos(): calls the following to do simulation.
            - setup_aper():  to setup the aperture (of type APER_T)
            - setup_powfs(): to setup the wfs type information (of type POWFS_T)
            - setup_recon(): to setup the wavefront reconstructor and dm fit structs (of type RECON_T)
            - sim(): to start the simulation. It then calls
                - save_skyc(): (optional) to save information for sky coverage postproc
                - FOR EACH SEED
                    - init_simu(): to initialize the run time structs
                    - genscreen(): (in closed loop) to generate atmosphere turbulence
                    - FOR EACH TIME STEP
                         - sim_update_etf(): (optional) to update the sodium profile
                         - genscreen(): (in open loop) to generate atmosphere turbulence
                         - perfevl(): (in closed loop) to evaluate the performance ins cience field
                         - wfsgrad(): to compute WFS gradients
                         - tomofit(): to do tomography and DM fit
                         - filter(): to do temporal filtering of DM commands.
                         - moao_recon(): (optional) to compute MOAO DM commands.
                         - perfevl(): (in open loop) to performance evaluation
                         - save_simu(): to save simulation telemetry data.
                         - print_progress(): to display progress data.
                    - free_simu(): to delete run time structs and close files

    \author Please write to Lianqi Wang <lianqiw at gmail.com> to obtain the source code.

 */


/**
   This is the routine that calls various functions to do the simulation. maos()
   calls setup_aper(), setup_powfs(), and setup_recon() to set up the aperture
   (of type APER_T), wfs (of type POWFS_T), and reconstructor (of type RECON_T)
   structs and then hands control to sim(), which then stars the simulation.
   \callgraph */
void maos(const PARMS_T *parms){    
    APER_T  * aper;
    POWFS_T * powfs;
    RECON_T * recon;

    aper  = setup_aper(parms);
    info2("After setup_aper:\t%.2f MiB\n",get_job_mem()/1024.);
    powfs = setup_powfs(parms, aper);
    info2("After setup_powfs:\t%.2f MiB\n",get_job_mem()/1024.);
    
    recon = setup_recon(parms, powfs, aper);
    info2("After setup_recon:\t%.2f MiB\n",get_job_mem()/1024.);
    /*
      Before entering real simulation, make sure to delete all variables that
      won't be used later on to save memory.
    */
    if(parms->dbg.evlol){
	sim_evlol(parms, powfs, aper, recon);
    }else{
	sim(parms, powfs, aper, recon);
    }
    /*Free all allocated memory in setup_* functions. So that we
      keep track of all the memory allocation.*/
    free_recon(parms, recon); recon=NULL;
    free_powfs(parms, powfs); powfs=NULL;
    free_aper(aper); powfs=NULL;
#if USE_PTHREAD == 2
    if(parms->sim.nthread>1)
	thr_pool_destroy(default_pool);
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

    char*fn=mybasename(argv[0]);
    if(!strncmp(fn, "scheduler",9)){//launch the scheduler.
	scheduler();
	exit(0);
    }
    char *scmd=argv2str(argc,argv);
    strcpy(argv[0],fn);
    free(fn);

    ARG_T* arg=parse_args(argc,argv);
    /*In detach mode send to background and disable drawing*/
    if(arg->detach){
	disable_draw=1;//disable drawing.
	//info2("Sending to background\n");
	daemonize();
	fprintf(stderr, "%s\n", scmd);
    }
    info2("MAOS Version %s. ", PACKAGE_VERSION);
#ifdef SVN_REV
    if(strlen(SVN_REV)>1 && strcmp(SVN_REV,"exported")){
	info2("Revision %s. ", SVN_REV);
    }
#endif
    info2("Launched at %s in %s.\n",myasctime(),myhostname());
    info2("Compiled on %s %s by %s ", __DATE__, __TIME__, __VERSION__);
#ifdef __OPTIMIZE__
    info2("with optimization.\n");
#else
    info2("without optimization!!!\n");
#endif
    //register signal handler
    register_signal_handler(maos_signal_handler);
  

    scheduler_start(scmd,arg->nthread,!arg->force);
    //setting up parameters before asking scheduler to check for any errors.
    PARMS_T * parms=setup_parms(arg);
    info2("After setup_parms:\t %.2f MiB\n",get_job_mem()/1024.);
    if(!lock_seeds(parms)){
	warning("There are no seed to run. Exit\n");
	maos_done(0);
	return 1;
    }
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
    free(scmd);
    free(arg->seeds);
    free(arg->dirout);
    free(arg->conf);
    free(arg);
    
#if USE_PTHREAD == 2
    //Create thread pool.
    if(parms->sim.nthread>1)
	default_pool=thr_pool_create(1,parms->sim.nthread,3600,NULL);
#endif
    dirsetup=stradd("setup",NULL);
    if(parms->save.setup){
	mymkdir("%s",dirsetup);
	addpath(dirsetup);
    }
    if(parms->sim.skysim){
	dirskysim=stradd("skysim",NULL);
	mymkdir("%s",dirskysim);
    }else{
	dirskysim=mystrdup(".");
    }

    /*Loads the main software*/
    maos(parms);
    maos_done(0);
    free_parms(parms);
    free(dirsetup);
    free(dirskysim);
    info2("Job finished at %s\t%.2f MiB\n",myasctime(),get_job_mem()/1024.);
    exit_success=1;//tell mem.c to print non-freed memory in debug mode.
    return 0;
}
