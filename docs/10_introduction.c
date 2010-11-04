/**
    \page page1 Introduction

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