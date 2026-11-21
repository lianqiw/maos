\page page10_intro Introduction

\tableofcontents

The Multithreaded Adaptive Optics Simulator (MAOS) started as a C
re-implementation of the algorithms in the MATLAB based Linear Adaptive Optics
Simulator (LAOS), which was limited in speed and multi-threading. The goal is an
efficient, general purpose, time domain adaptive optics simulator to help the
development of adaptive optics systems, in particular the advanced, very high
order systems for TMT and other extremely large telescopes (ELTs). MAOS is now
fairly feature complete for modeling astronomical adaptive optics systems such
as multi-conjugate AO (MCAO), laser tomography AO (LTAO), ground layer AO
(GLAO), multi-object AO (MOAO), and classic single conjugate AO with
Shack-Hartmann or pyramid wavefront sensors.

# What is modeled {#intro-model}

## Configuration {#intro-config}

Simulations are configured through configuration files and command line
options. The configuration files are plain text `key=value` lists that are easy
to read and maintain. Default values of all parameters are kept in the default
configuration files, so a user only needs to specify the parameters that differ.
MAOS checks the configuration for apparent errors and conflicts, and stops on
any unknown key. See \ref page30_run and \ref page34_conf.

## Atmosphere {#intro-atm}

Atmospheric turbulence is represented by one or more discrete optical path
difference (OPD) screens at different altitudes that evolve by frozen flow with
the given wind velocities. The turbulence strength \f$r_0\f$ can also vary with
time (`atm.r0evolve`). The turbulence is sensed by one or more natural guide
star (NGS) or laser guide star (LGS) wavefront sensors (WFS).

## Telescope {#intro-telescope}

The telescope aperture is modeled by a gray pixel amplitude map that can
represent segmented primary mirrors and the vignetting by the secondary mirror
supports. Aberrations of the telescope or instrument optics are modeled as
static OPD maps at their conjugation altitudes. Non-common path aberrations
(NCPA) that originate in the WFS or science optics are handled by automatically
biasing the WFS measurements.

## Wavefront sensors {#intro-wfs}

A WFS can be simulated as an idealized wavefront gradient sensor, a best
Zernike fit tilt sensor, or a physical optics WFS with user specified detector
pixel characteristics and center of gravity or matched filter pixel processing.
The pyramid WFS is implemented with modulation and optical gain optimization.
The LGS model includes spot elongation from a specified sodium layer profile
and an optional polar coordinate CCD. Subaperture or optical gain can be
optimized during the simulation by dithering the spot in a circular motion and
comparing the input and measured motion.

## Deformable mirrors {#intro-dm}

One or more deformable mirrors (DMs) are simulated as grids of actuators with a
given inter-actuator coupling. Actuator hysteresis, stroke limit,
inter-actuator stroke limit, and stuck or floating actuators can also be
modeled. Both zonal and modal control are supported.

## Reconstruction and control {#intro-recon}

The tomographic wavefront reconstruction estimates the turbulence at several
altitudes from the pseudo open loop gradients (the measured gradients plus the
contribution of the current DM commands). It uses a computationally efficient
implementation of the minimum variance reconstructor described in [Ellerbroek,
2002]. The reconstructed turbulence screens are then fit to the actuators of
one or more DMs to give the best correction over a specified field of view
(FoV). For LGS MCAO with separate low order WFS for tip/tilt and plate scale,
the so called split tomography algorithm reconstructs the low order modes
separately from the high order modes, with additional steps to decouple the
two. For classic AO, least squares reconstruction is also implemented. The DMs
are controlled by integrators whose gains can be optimized during the
simulation. For low order control, an LQG based Kalman filter is also
available.

## Performance evaluation {#intro-perf}

Performance is evaluated in terms of RMS wavefront error, Strehl ratio, and
point spread functions (PSFs) at science directions in the target FoV, which
can differ from the FoV used for DM fitting.

Many switches are available to isolate the contributions of individual
effects, such as DM fitting error, tomography error, telescope aberration, DM
hysteresis, or failed actuators. These simulations are used to optimize system
parameters and to build error budgets.

# Implementation {#intro-impl}

MAOS is written in C (C11 with GNU extensions; the code also compiles as C++)
and requires a POSIX compliant platform. It depends on the following external
libraries:

- FFTW version 3 (https://www.fftw.org) for Fourier transforms.
- BLAS and LAPACK (reference LAPACK, OpenBLAS, or Intel MKL).
- CHOLMOD from SuiteSparse
  (https://github.com/DrTimothyAldenDavis/SuiteSparse) for sparse Cholesky
  factorization.
- zlib for compressed files.
- Optional: the CUDA toolkit for GPU acceleration on Nvidia GPUs, GTK and cairo
  for the graphical tools, MATLAB for the mex interface, CMocka for unit tests,
  and zstd, OpenSSL, and libwebsockets for extra file and network features.

The configure script downloads FFTW, CHOLMOD, and some of the optional
libraries automatically when they are not found in the system. Autoconf,
automake, libtool, and make are needed to build from the git repository. MAOS
is built with GCC and Clang in continuous integration; the Intel compilers (icx
and the classic icc) are also recognized by the build system. See \ref
page20_compile for details.

Parallelism uses POSIX threads and OpenMP (no MPI), so a simulation runs
within one shared memory machine. It scales well with the number of cores. For
example, on a quad-core Intel Core i5 desktop of the time, the NFIRAOS
configuration (split tomography with physical optics LGS and geometric best
fit Zernike tilt NGS) ran at 7 seconds per time step with 1 thread and about 1
second per step with 8 threads. The same case took about 0.06 seconds per step
on 8 GTX 580 GPUs or 2 RTX 2080 Ti GPUs (in single precision). The memory
usage in the NFIRAOS baseline case is about 2 GB.

Two optional graphical programs are also included: `drawdaemon` for plotting
and `monitor` for monitoring jobs on one or more machines. They require GTK
(version 3 or 4; GTK 2 is still accepted). Jobs can also be monitored and plotted in a web browser
without GTK.

# Source tree {#intro-source}

The source code is organized in layers, each using only the ones above it:

- `sys`: the operating system layer: memory accounting, threads and thread
  pool, sockets, file input/output (`.bin` and `.fits`, optionally compressed),
  and the job scheduler client.
- `math`: numerical data structures (dense and sparse matrices, cell arrays,
  coordinates `loc_t`, maps `map_t`) and their functions, with wrappers for
  BLAS, LAPACK, FFTW, and CHOLMOD.
- `lib`: adaptive optics building blocks such as turbulence generation, ray
  tracing, WFS and DM models, and the configuration file parser. It is also
  built as `aolib.so` for use from Python.
- `cuda`: the GPU versions of the simulation and reconstruction code.
- `maos`: the end-to-end simulator.
- `skyc`: the sky coverage post-processor that takes `maos` results as input.
  See \ref skycoverage.
- `tools`: the job scheduler, the `monitor` (GTK and browser versions), the
  `drawdaemon` plotting program, and utilities such as `drawres`, `drawbin`,
  and `bin2fits`.
- `mex`: the MATLAB mex interface.
- `scripts`: Python, MATLAB, IDL, and Julia routines to read and write
  simulation results or to call `maos` library routines, plus test scripts.
- `config`: the configuration files.
- `test`: unit tests.

See \ref page93_maos for the structure of the simulator, and \ref
page20_compile for how to compile and use MAOS.
