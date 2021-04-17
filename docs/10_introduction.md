\page page10_intro Introduction

The Multithreaded Adaptive Optics Simulator (MAOS) started as a reimplementation
in C of the algorithms in the original MATLAB based linear adaptive optics
simulator (LAOS), which suffers performance and multi-threading limitations. The
goal is to create an efficient time domain general purpose adaptive
optics simulator to help the development of adaptive optics systems,
particularly advanced, very high order systems for TMT and other ELTs. MAOS now
is fairly feature complete in modeling astronomical adpative optics systems such
as, multiconjugate AO (MCAO), laster tomography AO (LTAO), ground layer AO
(GLAO), multi-object AO (MOAO) and Pyramid based classic AO.

Simulations are configured through configuration files and command line options,
which are checked for any apparent errors or conflicts. The configuration files
are very human readable and easy to maintain. Default parameters are maintained
in the default configuration files and only parameters that need to be overriden
need to be specified by the end user.

Atmospheric turbulence is represented as one or multiple discrete optical path
difference (OPD) screen(s) at different altitudes that evolve according to
frozen flow with given wind velocities. The turbulence strength $r_0$ can also
vary with time. The turbulence is sensed by one or multiple natural guide star
(NGS) or laser guide star (LGS) Shack-Hartmann wavefront sensor(s) (WFS).  

The telescope aperture is modeled with a gray pixel amplitude map that can
simulate segmented primary mirrors and secondary mirror support vignetting.
Aberration of the telescope or instrument optics are modeled as static OPD maps
with various conjugation. Non-common path aberrations that originate from the WFS
or science optics are handled by automaticaly biasing the WFS measurements.

The wavefront sensors (WFS) can be simulated as idealized wavefront gradient
sensors, best Zernike fit tilt sensors, or physical optic WFS using user
specified detector pixel characteristics and center of gravity or matched filter
pixel processing algorithm. Pyramid WFS is also implemented for the NGS with
modulation and optical gain optimization. The laser guide star wavefront sensing
model includes the impact of guide star elongation for a specified sodium layer
profile and an optional polar coordinate CCD. Automatic subaperture or optical
gain optimization can be done by dithering the spot in a circular motion and
compare the input and measured motion.

One or multiple deformable mirrors are simulated with grid of actuators that
couple together with a predefined coupling coefficient. The actuator hysteresis,
stroke limit, inter-actuator stroke limit, stuck or floating dead actuators can
also be modeled. Both zonal and modal control are supported.

The tomographic wavefront reconstruction estimates the turbulence in several
different heights from the pseudo open loop gradients (formed by adding the
deformation mirror (DM) correction back to the gradients) measured by the WFS,
using a computationally efficient implementation of a minimum variance
reconstruction algorithm as described in [Ellerbroek, 2002]. These reconstructed
turbulence screens are then fit to the actuators on one or several deformable
mirrors (DMs) to achieve best possible correction over a specified field of view
(FoV). For LGS based MCAO with separate low order WFS controlling tip/tilt and
plate scale, so called `split tomography` control algorithm is used that
reconstructs the low order modes separately from the high order reconstruction
with additional sanitization to decouple the two. For classic AO, least squares
control is also implemented. A simple integrator with gain update is used to
control the DM. For low order control, an LQG based kalman filter is also
implemented.

Performance evaluation is done in terms of RMS wavefront error, Strehl ratio,
and/or point spread functions (PSFs) at a few science objects in the target FoV,
which might be different from the FoV used for the DM fitting step.

Various knobs and switches are implemented to probe the contributions of various
effects, such as DM fitting error, tomography error, telescope aberration, DM
hysteresis error, failed actuators, etc. These simulations are used to optimize
system parameters and/or build the error budget.

This software is written in the C language (revision 99), with external
dependent libraries of FFTW version 3 (www.fftw.org), blas/lapack
(www.netlib.org/lapack/, Openblas, or Intel MKL), Cholmod
(www.cise.ufl.edu/research/sparse/SuiteSparse/ for Cholesky factorization), and
optionally CUDA SDK (for GPU acceleration). Automake, autoconf, libtool and make
are needed to configure the code and compile. It also requires POSIX compliant C
library for a few system calls. 

A C99 compliant compiler is required to compile the code. The cost has been
tested with GNU GCC, Intel ICC, and CLANG. For machines with Intel CPUs, the
code generally runs faster when compiled with ICC and linked with MKL. Nvidia
CUDA SDK is needed to compile the optional cuda code which runs in Nvidia GPUs
(>2.0 compute capability).

This software also contains two optionally compiled graphical executables
(drawdaemon, monitor) for plotting and job monitoring (jobs can be monitored
on several different machines). These two executable require GTK+.

The multi-threading is achieved using Posix Threads (no MPI is implemented) or
OpenMP. This parallelism works fairly well in multi-core, shared memory
machines. For example, in a recently purchased desktop as of this writing with
Intel Core i5 quad-core processor, the software (configured in a NFIRAOS split
tomography mode with physical optics LGS and geometric best fit Zernike tilt
NGS), can run at 7 seconds per simulation time step with 1 thread and ~1 seconds
per time step with 8 threads. The scaling with number of cores is fairly good.
Each time step takes only ~0.06 seconds on a 8x GTX 580 or dual RTX 2080 Ti GPU
setup (in single precision mode). The memory consumption in the NFIRAOS baseline
case is about 2 GB.

The code contains a few directories. The `sys` directory contains the lower
levels libraries that often deals with the system for input/output, and process
management. The `math` directory contains mathmatical data structures and
associated functions. The `lib` directory contains adaptive optics related data
structures and functions. The `maos` directory contains the end to end adaptive
optics simulation implementation. The `skyc` is contains routines for sky
coverage postprocessing that takes `maos` results as input. The `scripts` folder
contains useful Python, Matlab, and IDL routines that can be used to read/write
the simulation results or call maos internal routines. 

Please see \ref page20_compile on how to compile and use MAOS.

