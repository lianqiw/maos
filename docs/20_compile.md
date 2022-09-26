\page page20_compile  Compile the Code

\tableofcontents
<!--
This step can be skipped if the user choose to use the pre-compiled binary files,
which are available for all the popular platforms: Windows, Mac, and
Linux. Please notice that the job scheduling and monitoring utilities are not
included in the released binaries (the binaries may not be up to date).
-->
# Requirements

- Posix compatatible platform: Linux, macOS, or WSL for windows. BSDs may work
  but are not actively verified.

- C11 compliant compiler: GCC, Intel ICC, or clang, and make.

- FFTW version 3. Will download from MAOS site if not available in the system.

- Optimized blas and lapack. Blas and lapack can usually be installed through
  the linux distribution official repository. Will download from MAOS site if
  not available in the system. use --enable-mkl to force using Inte MKL.

- Cholmod. For sparse matrix cholesky factorization. Will download from MAOS
  site if not available in the system.

- If compile from git repo, autotools (autoconf, automake, libtool) are needed.

- (Optional) CUDA. For GPU acceleration.

- (Optional) GTK+. For drawdaemon and monitor.

- (Optional) libwebsockets. For web based job monitoring. It requires cmake to compile.

- (Optional) CMocka. For unit testing.

- (Optional) Doxygen. For generating documentation (this website)



# Preparing the folders and compiling

We recommend using different folders to 1) store the source code, 2) compile the
code (different directory can be used for different compiling options), and 3)
run simulations.

## Prerequisites

To build the source code that does not come with `configure`, a few additional utilities are needed: ``` automake, autoconf, libtool, make, cmake ```. 
In linux, use the native package manage to install. In macOS, we recommend using homebrew (https://brew.sh/) to install. 

## Download the code
### Option 1:
The preferred method is to use the always up-to-data Git version:
```
cd ~/work/programming
git clone https://github.com/lianqiw/maos.git maos
cd maos
./autogen.sh
export src_dir=$(pwd)
```

### Option 2:
Download the tag version that is not formally released:

```
cd ~/work/programming
curl -L https://github.com/lianqiw/maos/archive/refs/tags/v2.8.tar.gz -o maos.tar.gz
tar xvf maos.tar.gz
cd maos-2.8
./autogen.sh
export src_dir=$(pwd)
```


### Option 3:
There are released version that comes with `configure`. Use this if the utilities mentioned above are not available. In this case, do not try to update the `configure.ac` file which will trigger regeneration of `configure` and require these utilities.
```
cd ~/work/maos
#curl -L https://github.com/lianqiw/maos/releases/download/v2.8/maos_linux_amd64-2.8.0.tar.bz2 #for linux
#curl -L https://github.com/lianqiw/maos/releases/download/v2.8/maos_macos_amd64-2.8.0.tar.bz2 #for macos 
#curl -L https://github.com/lianqiw/maos/releases/download/v2.8/maos_macos_aarch64-2.8.0.tar.bz2 #for macos apple silicon
tar xvf maos_*.bz2
```

### Option 4:
Instead of compiling the code, there are also pre-compiled binary that can be downloaded. This option works the best when there is no desire to modify the code.
```
cd ~/work/programming
curl -L https://github.com/lianqiw/maos/archive/refs/tags/maos-2.8.0.tar.gz
cd maos_2.8.0
export src_dir=$(pwd)
```
## Compile the Code
Next we create another folder, where we are going to compile the code. example:

```
cd ~/work/maos
mkdir optim && cd optim
$src_dir/configure
make -j4 #this will compile using 4 threads.
```

The compiled executable is maos in the sub-folder `bin` of the compiling
folder. You do not have to do `make install` to run the simulations.

\subsection sect-compiler Compiler options

The default compiler for the system is used with optimization flags. There
are a few optional switches or additional components that can be enabled when
calling configure.

To use a different compiler or enable debugging:
```
$src_dir/configure CC=icc    #ICC compiler
$src_dir/configure CC=clang  #clang compiler. Also used for compiling cuda code if cuda is detected.
$src_dir/configure CC=gcc4.5 #use gcc4.5
$src_dir/configure --enable-debug #enable debugging (debug is not enabled by default)
$src_dir/configure --disable-openmp #Use pthreads instead of openmp (openmp is the default)
$src_dir/configure --disable-double #Use single precision (double is the default)
$src_dir/configure --disable-scheduler #disable the built-in scheduler (for CPU cluster with job management)
```

\subsection sect-cuda GPU acceleration

To enable GPU computing, first install nvidia graphics driver and cuda
software development toolkit (SDK). Issue the following command to ensure
graphics configuration is correct: \c nvidia-smi. You should see all your
nvidia graphics gpus are listed.

The configure can automatically detect cuda if nvcc is in the path, or if
cuda is installed in /usr/local/cuda. For other locations, you
need to pass the cuda directory to configure:

`$src_dir/configure --with-cuda=$cuda_dir`

When GPU compute is compiled in maos and Nvidia GPUs are available, maos
will make use of all suitable GPUs to run simulation. If that is not desired, the
GPU compute can be disabled at run time by passing -g-1 to maos. Different
set of GPU can also be selected by passing one or multiple -gi or -Gn:

```
$maos -g-1    #disable GPU compute
$maos -g0 -g2 #Only use gpu 0 and 2
$maos -G2     #Automatically choose 2 gpus.
```

\subsection sect-mex Matlab Mex Routines

To compile mex routines for matlab, pass the matlab directory to configure if it is not in \c $PATH:

```
$src_dir/configure --with-matlab=$matlab_dir
```

The compiled mex routines are in \c mex folder in the compilation
directory.

- The \c read() and \c write() are for handling .bin and .fits files.
- The \c aolib() exports many useful routines for use in \c matlab. Run \c aolib() to see the list of available functions.
- The \c maos() can be used to run maos within matlab. Run \c maos() to see instructions.

\subsection sect-mac-gtk Installing GTK+ in MAC OS and Compile Monitor, Drawdaemon

For macOS, it is simplest to install gtk and gtk-mac-integration. using macport
or homebrew. No additional configuration is required. For manual install, follow
the instructions on [page](https://www.gtk.org/docs/installations/) to install
gtk+ for osx and configure the envirionmental parameters. Version 2 and 3 are
supported. Make sure pkg-config is in the path.

Now rerun MAOS configure and make. Monitor and drawdaemon should appear in bin folder now. To select gtk version when multiple ones are available:
```
$src_dir/configure --disable-gtk-3 #disables gtk-3
$src_dir/configure --disable-gtk-2 #disable gtk-2 (gtk-3 is used by default)
```

# Graphical User Interface

When GTK+ libraries are found in the system, additional executables will
be compiled, \c drawdaemon, \c drawres, \c drawbin and \c monitor in the \c bin folder.

- The plotting utility \c drawdaemon is launched automatically when plotting commands are issued in simulation (`maos plot.all=1`) or in post processing with \c drawres or \c drawbin

- The monitoring utility \c monitor can be used to monitor jobs in this and linked machines. 

- The \c drawres can be used to automatically plot the results folder.

- The \c drawbin can be used to automatically plot the OPDs or maps.

## Monitor

In order to use monitor to monitor jobs in different machines, follow these steps:

Create a folder .aos in your home folder, create two files in there. Put
these files in all the machines you want to monitor.

- \c ~/.aos/hosts Each line contains a machine's hostname you want to monitor. It also supports the
full format \c prettyname=full_hostname:port

- \c ~/.aos/port contains the port you want the \c scheduler to bind to so that
  the monitor can connect to it.  Any non-privileged port numbers are fine (like
  10000).  Different users should have different port number.  The port has to
  be the same on all the machines unless explicitly specified in \c hosts. If
  not set, the number will be automatically computed based on username. Make
  sure the firewall allows them. 

- Please be aware that no authentication is performed. The user should make sure such ports are protected.

## Drawdaemon

Once the monitor is running, it will display running jobs on each connected host.
Right click on a running job and select `Plot selected jobs` will launch drawdaemon
and start plotting the telemetry data in realtime. This works for local and remote host.
Multiple drawdaemon windows can be started for the same or different jobs.

An alternative way to launch drawdaemon is to include \c plot.setup=1 \c
plot.run=1 or \c plot.all=1 in maos command line.

Remote plotting: when a job is running remotely or in a headless machine, and
`drawres`, `drawbin`, or `maos` plotting is issued, the `monitor` can open a
`drawdaemon` in the local computer and plot remotely.

## Plotting results

`$drawres result_folder` #can be used to plot the results. Wildcards are supported to match many folders. 
It will launch drawdaemon using an active monitor in a remote machine if the current shell does not have DISPLAY set. 
This can be used to view results remotely.


# Python Scripts 

MAOS comes with a few useful python routines location in \c scripts/ folder:

- `readbin.py` contains `readbin` that can read `bin` or `fits` files.

- `aolib.py` contains routines to read and process \c MAOS results. It also imports \c libaos.py and \c draw.py

- `libaos.py` contains routines that wraps \c MAOS library routines via \c aolib.so and \c ctypes 

- `lib2py.py` generates \c libaos.py

- `draw.py` used to draw opd defined on coordinate loc. \c draw(loc, opd) or \c draw(map) where \c map is a 2d array.

- `maos_client` provides an interface to a running maos session. See next subsection.

For `libaos.py` to work, an environment variable \c MAOS_AOLIB need to set to
the path of `aolib.so` or `aolib.dylib`

## Interface to MAOS

It is possible to retrieve or override MAOS internal data when it is running. See the following steps:

First, set up shell environment so that the libaos.py can find it.

- `export MAOS_AOLIB=\path\to\bin\aolib.so`
  
Then open \c ipython (or \c python) and run the following

- \c import \c maos_client  #import the module
- \c maos_client.connect(host,port) #the host and port is the same as used by monitor
- \c maos_client.get_list() #get list of available variables
- \c maos_client.get_var(name)  #get variable with the name
- \c maos_client.get_all()  #get all variables
- \c maos_client.pause(0)   #unpause maos
- \c maos_client.pause(1)   #pause or step maos.
