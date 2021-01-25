/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/**
    \page page10_intro Introduction

    MAOS is a complete reimplementation of the algorithms in the original MATLAB
    based simulator LAOS, which suffers performance and multi-threading
    limitations. The motivation to develop this software is to create an
    adaptive optics simulator that runs fast, consumes less memory, and does the
    job without MATLAB, which is proprietary and has a large memory
    footprint. The code is mostly written in C language with a function oriented
    design. Lately it acquired the ability to employ Nvidia GPUs for ultra fast
    execution using the CUDA programming model.

    MAOS is configured through configuration files and command line
    options. The code checks the configuration for any apparent errors or
    conflicts. The configuration files are very human readable and easy to
    maintain. The development goal is to make an efficient, easy to use, general
    purpose adaptive optics simulator to help the development of adaptive optics
    systems, particularly advanced, very high order systems for TMT and other
    ELTs.

    Atmospheric turbulence is represented as one or multiple discrete screen(s)
    at different heights that evolve according to frozen flow with given wind
    velocity. The resulting aberrations are sensed by one or multiple natural
    guide star (NGS) or laser guide star (LGS) Shack-Hartmann wavefront
    sensor(s) (WFS).  The sensors can be simulated as idealized wavefront
    gradient sensors, best Zernike fit tilt sensors, or a physical optic WFS
    using user specified pixel characteristics and center of gravity or matched
    filter estimation algorithm. Pyramid WFS is also implemented for the NGS
    with modulation and optical gain optimization which improves the
    magnitude limit.

    The laser guide star wavefront sensing model includes the impact
    of guide star elongation for a specified sodium layer profile and an
    optional polar coordinate CCD. The tomographic wavefront reconstruction
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
    (www.netlib.org/lapack/). The code contains a local copy of the Cholmod
    (www.cise.ufl.edu/research/sparse/SuiteSparse/, to do Cholesky
    factorization, copied from suite-sparse).  An optimized blas library such
    ATLAS, GOTO-BLAS, or Intel MKL is necessary to get good performance for
    Cholesky decompositions of the tomography and/or fitting matrix. Use
    --enable-mkl when running configure to download a redistributed copy of
    MKL. Automake, autoconf, libtool and make are needed to configure the code for
    compiling.

    A C99 compliant compiler is required to compile the code. It usually
    compiles smoothly using GNU GCC, Intel ICC, and CLANG. For machines with
    Intel CPUs, the code generally runs faster when compiled with ICC and linked
    with MKL. Nvidia CUDA SDK is needed to compile the optional cuda code which
    runs in Nvidia GPUs (>2.0 compute capability).

    This software also contains two optionally compiled graphical executables
    (drawdaemon, monitor) for plotting and job monitoring (jobs can be monitored
    on several different machines). These two executable requires GTK+ and the
    Cairo development and run-time library.

    The multi-threading is achieved using Posix Threads (no MPI is implemented)
    or OpenMP. This parallelism works fairly well in multi-core, shared memory
    machines. For example, in a recently purchased desktop as of this writing
    with Intel Core i5 quad-core processor at 2.66 GHz, the software (configured
    in a NFIRAOS split tomography mode with physical optics LGS and geometric
    best fit Zernike tilt NGS), compiled with Intel C++ compiler version 10, can
    run at 10 seconds per simulation time step with 1 thread and 3 seconds per
    time step with 4 threads. The same setup in a dual Intel W5590 server runs 7
    seconds per simulation time step with 1 thread and 4 seconds per step with 2
    threads. The scaling with number of cores is fairly good. The memory
    consumption in the NFIRAOS baseline case is about 2 G.

    The code contains a few directories. The \c sys directory contains the lower
    levels libraries that often deals with the system for input/output, and
    process management. The \c math directory contains mathmatical data
    structures and associated fucntions. The \c lib directory contains adaptive
    optics related structures and functions. The \c maos directory contains the
    MAOS simulation implementation. The \c skyc is contains routines for sky
    coverage postprocessing routine that takes \c maos results as input.

    \ref page20_compile

 */
