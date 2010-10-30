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

    - \ref page6 <p>

    - \ref page7 <p>

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


    \author Lianqi Wang <lianqiw@gmail.com> at TMT Corporation http://www.tmt.org
    
    The source code can be obtained in http://github.com/lianqiw/maos

 */
/**
    \page page1 Introduction to MAOS

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
    Cholmod (http://www.cise.ufl.edu/research/sparse/SuiteSparse/, to do
    Cholesky factorization, copied from suite-sparse). It will be compile if
    cholmod is not found in the system. An optimized blas library such ATLAS,
    GOTO-BLAS, or Intel MKL is necessary to get good performance for Cholesky
    decompositions of the tomography and/or fitting matrix. Use --enable-mkl
    when running configure to download a redistributed copy of MKL.

    A C99 compliant compiler is required to compile the code, such as GCC 4 or
    Intel C++ compiler. The code has been successfully compiled on GNU/Linux
    using GNU GCC and the Intel C++ compiler ICC (in C mode), and in MAC OS X
    10.5 using the GCC compiler. For machines with Intel CPUs, the code
    generally runs faster when compiled with ICC and linked with MKL.

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
   \page page2 Compile and Run the Code

   \section sect-requirement Requirements

    - C99 compliant compiler: GCC4 or ICC.

    - FFTW version 3. Can usually be installed through the linux distribution
    official repository or using source from www.fftw.org
    
    - Optimized blas and lapack. Blas and lapack can usually be installed
    through the linux distribution official repository. For 64 bit linux
    machines, it is also possible to use a redistributed version of the Intel
    MKL by appending --enable-mkl when you run configure
    
    - (Optional) GTK+ and Cairo. For drawdaemon and monitor.

    - (Optional) cholmod. If not found in the system, will compile and use the
      shipped cholmod.

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
   \page page3 Defining Simulation Parameters
   
    \section sect-simu-setup Defining Simulation Parameters

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

    The recommanded way to organize config files for a specific problem is to
    keep a baseline .conf file in sub-folder “config/maos” in the source folder
    that specify all the necessary key=value entries through optionally
    including other .conf files. Embed the options in the command line (embrace
    it with quote if the option contains spaces, like arrays) or put in one or
    multiple overriding file (in the simulation folder, not inside
    “config/maos”) to override the values you want to change from the baseline.

    The .conf files in the “config/maos” folder are only modified during code
    development and should not be modified by the end user. Use overriding .conf
    files in simulation folders or embed configurations in the command line to
    modify the values.

    In the “config/maos” sub-folder in the source folder, there are three baseline
    .conf files: 
    
    - nfiraos.conf for NFIRAOS mcao simulations,
    - scao_lgs.conf for single conjugate LGS simulation
    - scao_ngs.conf for single conjugate NGS simulation. 

    Each of this file include a few other .conf files by specifying
    include=”filename.conf”, for example, all the three baseline .conf files
    include 1) “sim.conf” to specify closed loop information such as loop gain,
    number of time steps, etc, 2) “recon.conf” to specify reconstruction
    parameters in closed loop, 3) “dbg.conf” to specify debugging
    parameters. The MCAO baseline “nfiraos.conf” includes “atm_mk13n50p.conf” to
    specify MK13N median seeing turbulence profile, while the two SCAO baselines
    include “atm_scao.conf” that specifies the turbulence information with only
    a single reconstruction layer.

    \ref conf0 shows the full content of nfiraos.conf
*/


/**
   \page page4 Interpreting Simulation Results
   
   MAOS will generate output in binary format \c .bin or zipped \c .bin.gz files
   as described below. There are two MATLAB mex routines \c read and \c write
   that can read and write \c .bin or \c .bin.gz files in MATLAB. The source of
   these mex routines are located in sub-folderx \c mex. Cd to this folder and
   issue \c make to compile these mex routines.  Copy \c write.mexa64 and \c
   read.mexa64 into a folder that is in your matlab path, such as
   $HOME/matlab. The usage in MATLAB is as follows:

\verbatim 
   cle=read('Res_1'); 
or cle=read('Res_1.bin'); 
or cle=read('Res_1.bin.gz');
   each wil try to read Res_1.bin or Res_1.bin.gz 
   
   write(cle,'Res_1');     will write the data to Res_1.bin.gz;
   write(cle,'Res_1.bin'); will write the data to Res_1.bin without compression.
\endverbatim

There will be several files created during simulation in the result folder. The
number after underscore _ is the seed. For example, with seed 1 the following
files are produced. Read in these .bin or .bin.gz files using provided mex
function \c read in MATLAB.  

 - \c Res_1.bin: A binary file containing a cell array that include the main
results. i.e. res=read('Res_1'); 
     - \c res{1} contains the open loop wavefront variance (WFV in units of \f$m^2\f$) in row vectors. The rows are 
          - Piston removed WFV
          - WFV in tip/tilt modes
          - Piston/tip/tilt removed WFV
     - \c res{2} contains the residual wavefront variance after the tomography phase screen is directly applied as the correction in open loop mode. Empty is evl.tomo is zero.
     - \c res{3} contains the closed loop wavefront variance in row vectors. The rows are
          - Piston removed WFV
          - WVF in tip/tilt modes
          - Piston/tip/tilt removed WFV
     - \c res{4} (Only in split tomography) contains the closed loop wavefront variance. The rows are
          - WFV in LGS contains modes
          - WFV in NGS Tip/Tilt modes
          - WFV in NGS modes (includeing Tip/Tilt and additional modes controlled by NGS (On-Instrument WFS)

 - \c Resolep_1.bin: Open loop wavefront variance for each science field point. Each cell represent a science field point. The format is similar to res{1} above.

 - \c Resolmp_1.bin: Open loop wavefront Zernike (not normalized wrt radius)
   modes defined on not-normalized coordinate on the aperture. The format is
   similar as \c Resolep_1.bin
 
 - \c Resclep_1.bin: Closed loop wavefront variance for each science field
   point, in the same format as Resolep_1.bin

 - \c Resclmp_1.bin: Close loop wavefront Zernike modes, in the same format as
   Resolmp_1.bin

 - \c Resclemp_1.bin: (Only in split tomography) LGS/TT/NGS mode wavefront error
   for each direction.

 - \c RescleNGSm_1.bin: (Only in split tomography) contains a row vector array
   of either the residual NGS modes (in radian like unit) or the applied NGS
   mode correction if ideal NGS mode correction is made. Mainly used to save the
   ideal NGS modes applied in skycoverage pre-simulation.

 - \c maos_975.conf: The final effective arrays of the configurations for MAOS
   run with PID 975. Can be used to reproduce the simulation or check the
   configurations.

 - \c sanea_sim_1.bin: When wavefront sensors are running in physical optics
   mode, the average gradinet measurement error for that wavefront sensor is
   saved (in order) in a cell in this file. Each cell is a column vector with
   elements twice the number of subaperture. The first half is for x (or radial)
   gradient, and the second half is for y (or azimuthal) gradient. They are in
   units of \f$rad^2\f$.

 - \c evlpsfmean_1.bin: When evl.psfmean is 1, contains the time averaged
   psf. if is a cell array of \f$n_{wvl}\times n_{evl}\f$. Each cell contains
   a center cut of the science PSF.

 - \c evlpsfhist_1_ievlx.bin: When evl.psfhist is 1, each of these files
   contains the time history of the complex PSF of evaluation direction x.

 - \c Resuptcmd_1.bin: Each cell contains the uplink tip/tilt command (only none
   empty for LGS WFS) history in unit of radian.

 - \c Resupterr_1.bin: Each cell contains the uplink tip/tilt error history in
   unit of radian.
 */
/**
   \page page5 Programming Guidelines

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
   \page page6 Profiling MAOS 

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
   \page page7 C Fundementals

   The user is required to have a good understanding of the C fundementals, like
   the pointer, to study and improve the code. Here we list a few to remind the
   user. 

   \section sect-array Arrays
   
   If you declare an array: <code>double a[1024];</code> the symbol actually
   contains the address of a 1024x8 byte memory region that is automatically
   allocated in the stack (malloc, calloc allocations memory in the
   heap). Suppost the address of the memory region is 0xFF00. Assignment
   <code>a[0]=1;</code> puts 1 in memory started at 0xFF00 <code>a[1]=1;</code>
   puts 1 in memory started at 0xFF08.

   Stack is the memory region allocated dynamically when calling a function. All
   variables reside in stack are automatically deleted when the function
   returns. 
   
   Heap is the memory reservoir. The functions malloc, calloc, realloc reserves
   memory in the heap while free removes the reservation. 

   \section sect-pointer Pointers 

   A pointer is a variable that stores the address of a memory region.
   <code>void *p;</code><p> is a basic pointer that just stores the address and can
   not do any thing about the memory it points to. Its value can be assigned to
   any other non const vectr or from any vector (1d or 2d).

   However, a normal pointer like <code>double *p;</code> also contains
   information about the memory it points to. More specifically, the compiler
   know the length of each element so that it can compute the address offset of
   the elements, just like arrays. 

   For example: <code>double *p;</code> if p contains value of 0xFF00.<p>
   <code>p[0]=1;</code>will store a double 1 into memory region at 0xFF00<p>
   <code>p[1]=1;</code>will store a double 1 into memory region at 0xFF08
   because a double requires 8 bytes.<p> <code>p=p+1;</code>p will now contain
   0xFF08, i.e., points to the next element, instead of 0xFF01<p>
   <code>p[0]=1;</code> is strictly equivalent to <code>*p=1;</code> *p means
   visit the memory pointed by p instead of p itself.  <code>p[1]=1;</code> is
   strictly equivalent to <code>*(p+1)=1;</code><p>

   
   \subsection sect-1d-pointer 1D Pointers

   The simplist usage of pointer is to create and index an 1-D array:
   <code>double *p;</code>.  

   <code>p=calloc(1024, sizeof(double));</code> Allocate a memory block of
   1024x8 byte and store its address in p.<p> <code>p[0]=1;</code> will assign 1 to
   the first number in the array<p>

   This is equivalent to using the arrays.

   <code>double a[1024];</code> Declare an 1-D array.

   <code>a[0]=1;</code > will assign 1 to the first number in the array.  The
   symbol a is like a constant pointer. We can assign it's value to a constant
   pointer: 

   <code>double *const p=a;</code><p> <code>p[0]=1;</code><p>

   \subsection sect-const Const pointers vs. pointers point to const memory

   Beware about a const pointer and pointer that points to constant memory
   regions: 

   <code>const double *p;</code>p is a pointer that points to a memory
   region that can not be changed. Thus <code>p[0]=1;</code> is illegal.

   <code>double * const p;</code>p is a constant pointer that can not have its
   value changed. Thus <code>p=a</code> is illegal, unless you assign at
   declaration: <code>double *const p=a</code>. However, it is legal to modify
   the values in the memory p points to: <code>p[0]=1;</code>.

   \subsection sect-2d-pointer 2D Pointers 

   A 2D pointer is like a 1D pointer. It stores the address of a memory
   region. However, unlike the 1D pointer, the compiler knows some additional
   infos about the 2D pointer. The compiler knows the length of the first
   dimension so that is can compute the memory location given 2D index like
   a[2][3]; <code>double (*p)[2]=calloc(3*2,sizeof(double))</code> can be used
   in the same way as <code>double a[3][2]</code>. The last index changes
   fastest.
   
   The pointer p in our dmat (or cmat) struct points to the allocated memory. It
   can be cast to a two-d vector to facilitize visiting of the elements. The
   cast is conveniently done by PDMAT(a,pa) for dmat a, which simply declares a
   2-d pointer and assign a->p to it: 
   
   <code>double (pa*)[2]=(void*)a->p</code>. a->p is first cast to void *
   because it is a double * and can not be assigned to double (*)[2]
   directly. The array elements (ix,jy) can be accessed by pa[iy][ix]. Notice
   that ix changes fastest.
 */


/**
    \page conf0 nfiraos.conf

    The following shows the content of nfiraos.conf, which is the baseline
    configuration of TMT NFIRAOS system

    \include config/maos/nfiraos.conf

    
    
*/
