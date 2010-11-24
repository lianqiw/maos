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
