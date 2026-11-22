\page page93_maos Architecture

This page describes how the code is organized and how a \c maos simulation
runs. See \ref page91_basics for the data types and \ref page95_profiling for
debugging.

\section arch-layers Library Layering

The source is organized as layered libraries. Each layer depends only on the
layers above it (this is also the \c SUBDIRS build order in \c Makefile.am):

1. \c sys/: operating system layer. Memory accounting (\c mem.c), threads and
   the thread pool (\c thread.c, \c thread_pool.c), sockets (\c sock.c,
   \c sockio.c), file IO with the \c bin and \c fits formats and gzip/zstd
   compression (\c bin.c, \c zstdio.c), the scheduler client
   (\c scheduler_client.c), and per OS process information (\c process_*.c).
2. \c math/: numerical containers (#dmat, #cmat, #smat, #zmat, #lmat, sparse
   #dsp/#csp, #cell arrays) with BLAS, LAPACK, FFTW, and
   CHOLMOD wrappers. The functions are type generic (see
   \ref sect-type-generic).
3. \c lib/: physics and AO building blocks reusable outside of \c maos, built
   into \c libaos together with \c sys and \c math: the point and grid
   types #loc_t, #pts_t, and #map_t derived from #dmat (\c loc.c, \c map.c),
   turbulence generation,
   ray tracing (\c accphi.c), the \c mkg, \c mkh, \c mkz, \c mkw operators,
   PCG, servo, PSD, pyramid WFS, Shack-Hartmann WFS subaperture geometry and
   CoG (\c shwfs.c), the laser launch telescope (\c llt.c), NEA conversions
   (\c nea.c), dithering (\c dither.c), the Fourier domain preconditioner of
   the tomography (\c fdpcg.c), the configuration file parser
   (\c readcfg.c), and the plotting client (\c draw.c). Functions here take
   basic types, not \c maos structs. \c lib/aolib.c builds \c bin/aolib.so,
   which is used by the Python (\c scripts/libaos.py) and Matlab interfaces.
4. \c cuda/: GPU versions of the simulation and reconstruction, organized as
   \c cuda/math, \c cuda/recon, and \c cuda/sim. They mirror the CPU code
   paths and are compiled only when CUDA is found (\c USE_CUDA).
5. \c maos/: the simulator.
6. \c skyc/: sky coverage post processing, a separate executable with its own
   parameters (see \ref skycoverage).
7. \c tools/: the job \c scheduler, the \c monitor (GTK and browser
   interfaces), \c drawdaemon (GTK plotting), \c drawres, \c bin2fits, and
   other utilities.
8. \c mex/ and \c scripts/: Matlab and Python bindings and scripts.
   \c scripts/lib2mex.py and \c scripts/lib2py.py parse the headers to
   generate \c mex/aolib.c and \c scripts/libaos.py, so the formatting of
   structs and doc comments in those headers matters.

\section sect-structure Simulation Flow

The following is the call sequence of \c maos. Function names link to their
documentation.

- main() (\c maos/maos.c) parses the command line, then
    - setup_parms() (\c maos/parms.c) reads the configuration files into
      #parms_t and checks for errors:
        - The \c readcfg_*() functions read every key, for example
          \c readcfg_sim(), \c readcfg_powfs(), \c readcfg_wfs(),
          \c readcfg_dm(), \c readcfg_tomo(), \c readcfg_gpu().
        - All keys that are read are saved to \c maos_HOST_PID.conf (linked
          as \c maos_recent.conf) in the output folder so that the simulation
          can be reproduced.
        - The \c setup_parms_postproc_*() functions derive dependent values
          and validate the parameters.
    - scheduler_start() waits for permission from the scheduler to start (see
      \ref arch-scheduler).
    - setup_parms_gpu() selects the GPUs and sets the \c parms->gpu flags.
    - \c THREAD_POOL_INIT() initializes the thread pool or OpenMP with
      \c NTHREAD threads (option \c -n).
    - maos_setup() (\c maos/setup.c) sets up the AO system geometry and data
      structs. Its main steps are:
        - setup_aper(): the aperture grid and amplitude map (#aper_t).
        - setup_powfs_init(): the WFS subaperture geometry (#powfs_t).
        - setup_recon_prep(): the reconstruction grids (ploc, aloc, floc,
          xloc) and the ray tracing operators (#recon_t).
        - setup_surf(): surface errors and NCPA.
        - setup_shwfs_phy(): physical optics WFS data (dtf, etf, i0, matched
          filter, cog).
        - setup_powfs_neasim(): gradient noise model for geometric WFS.
        - setup_powfs_calib(): NCPA calibration.
        - setup_recon_GA(), setup_recon_GF(), twfs_setup_GR(): actuator,
          focus, and radial mode to gradient operators.
        - ngsmod_prep(): low order NGS mode geometry for split tomography.
        - setup_recon_saneai(): noise covariance.
        - setup_recon_dither_dm(): DM dithering.
        - setup_recon_control(): tomography, DM fitting, or least squares
          reconstructor.
        - setup_recon_moao(): MOAO, if enabled.
        - setup_powfs_fit(): fitting of turbulence onto the WFS subapertures
          (\c sim.wfsalias=2 or \c sim.idealwfs=2).
        - setup_recon_mvm() or gpu_setup_recon_mvm(): the control matrix when
          \c recon.mvm is enabled.
        - gpu_wfsgrad_init(), gpu_setup_recon(), gpu_perfevl_init(),
          gpu_wfssurf2gpu(): copy data to the GPU when the corresponding
          \c parms->gpu flag is set.
        - setup_recon_misc(): PSD and other miscellaneous data.
        - plot_setup(): plot the geometry when \c plot.setup is set.
    - maos_sim() (\c maos/sim.c) runs the simulation for each seed in
      \c sim.seeds:
        - maos_iseed() calls sim_init() to create the runtime data
          (#sim_t) and generates the frozen flow turbulence with genatm().
        - In the parallel (\c PARALLEL=1) and serial (\c PARALLEL=0) modes,
          maos_isim() is called for each time step from \c sim.start to
          \c sim.end. It calls \c prepare_isim() (turbulence for non frozen
          flow, sodium profile update, DM projection), then
            - with \c PARALLEL=1, queues reconstruct(), perfevl(), and
              wfsgrad() to run concurrently, waits for them, then calls
              shift_grad() and filter_dm(). The GPU versions of perfevl and
              wfsgrad are queued first so that they start early.
            - with \c PARALLEL=0 in closed loop, calls perfevl(), wfsgrad(),
              reconstruct(), shift_grad(), and filter_dm() in sequence.
            - with \c PARALLEL=0 in open loop, calls wfsgrad(), shift_grad(),
              reconstruct(), filter_dm(), and perfevl() in sequence.
            - \c postproc_isim() updates the reconstructor when requested
              (for example by Cn2 estimation or dithering).
            - print_progress() reports the progress.
        - In the event driven mode (\c PARALLEL=2), three loops run
          concurrently for all time steps and synchronize with each other
          through counters:
            - perfevl() and print_progress(),
            - wfsgrad() and shift_grad(),
            - reconstruct(), filter_dm(), and \c postproc_isim(), which
              updates the reconstructor when requested (for example by Cn2
              estimation).
        - sim_free() frees the runtime data of the seed.
    - maos_reset() and free_parms() free the remaining data.

The main functions of each time step are:

- perfevl() (\c maos/perfevl.c): ray traces the turbulence and DM shapes to
  the science directions and computes the wavefront error and PSFs.
- wfsgrad() (\c maos/wfsgrad.c): computes the WFS OPD, the subaperture images
  in physical optics mode, and the gradients.
- shift_grad(): copies the new gradients to \c gradlastcl/\c gradlastol for
  the reconstructor.
- reconstruct() (\c maos/recon.c): computes the DM error signal from the
  gradients.
- filter_dm() (\c maos/filter.c): applies the servo filter to compute the DM
  commands (\c dmcmd) and the actual DM shape (\c dmreal), including
  hysteresis (see \ref hysteresis).

In closed loop, reconstruct() uses the gradients of the previous time step
(\c simu->reconisim is \c isim-1), which is what allows it to run concurrently
with wfsgrad() and perfevl().

\c PARALLEL is read from the environment variable \c MAOS_PARALLEL (default 1;
see read_sim_env()). maos_setup() changes it to 0 in open loop
(\c sim.closeloop=0) or when \c evl.tomo is set. \c PARALLEL=2 requires OpenMP
and does not support \c sim.dmproj; otherwise it falls back to 1. The
environment variables \c MAOS_NO_WFS, \c MAOS_NO_EVL, and \c MAOS_NO_RECON
skip the corresponding calls, which is useful for timing.

\section arch-structs Key Data Structures

- #parms_t (\c maos/parms.h): all the input parameters, read once by
  setup_parms() and treated as constant afterwards. The members mirror the
  configuration keys, for example \c parms->powfs[ipowfs] (#powfs_cfg_t),
  \c parms->wfs[iwfs] (#wfs_cfg_t), \c parms->dm[idm] (#dm_cfg_t),
  \c parms->sim (#sim_cfg_t), and \c parms->gpu (#gpu_cfg_t).
- #aper_t (\c maos/types.h): the telescope aperture: the evaluation grid
  (\c locs), the amplitude, and the PSF computation data (\c locfft).
- #powfs_t: data shared by all WFS of the same type (geometry, amplitude, and
  physical optics data such as the subaperture OTFs and matched filters).
- #recon_t: the reconstruction grids and operators (\c ploc, \c aloc,
  \c xloc, \c GA, \c GX, \c HXF, ...), the tomography and fitting parameters,
  and the control matrix.
- #sim_t: the runtime data of one seed: atmosphere, DM commands and shapes,
  gradients, performance results, saved telemetry, and the synchronization
  state for the parallel modes.
- #global_t: holds pointers to the above (\c global->parms, \c global->powfs,
  \c global->aper, \c global->recon, \c global->simu) for use by the
  scheduler listener, the signal handler, and the Matlab/Python interfaces.

\section arch-gpu GPU Code Paths

When \c maos is compiled with CUDA, the \c cuda/ library provides GPU versions
of the WFS (\c gpu_wfsgrad*), performance evaluation (\c gpu_perfevl*), and
reconstruction (\c gpu_setup_recon, \c gpu_tomo, \c gpu_fit, ...). Each task
uses the GPU when its flag in \c parms->gpu is set: \c gpu.wfs, \c gpu.evl,
\c gpu.tomo, \c gpu.fit, \c gpu.lsr, \c gpu.psf, and \c gpu.moao (defaults in
\c config/maos/sim.conf). setup_parms_gpu() turns off the flags that do not
apply to the configuration, sets \c gpu.recon, initializes the GPUs with
gpu_init(), and clears all the flags if no GPU is used. The command line option \c -g selects
the GPUs (\c -g-1 disables GPU) and \c -G sets the number of GPUs. The CPU
code calls the \c gpu_* functions inside \c \#if \c USE_CUDA blocks
depending on these flags.

\section arch-param Adding a Parameter

Every runtime option is a \c key=value entry in the \c config/maos/*.conf
files and maps 1:1 to a field of #parms_t: \c wfs.thetax becomes
\c parms->wfs.thetax (or \c parms->wfs[iwfs].thetax for per WFS keys). To add
a parameter:

1. Add the field with a doc comment to the struct in \c maos/parms.h.
2. Read it in the matching \c readcfg_*() function in \c maos/parms.c with one
   of the \c READ_* macros, for example \c READ_INT(sim.foo),
   \c READ_DBL(sim.foo), \c READ_STR(sim.foo), or \c READ_DMAT(sim.foo), or
   \c READ_POWFS(dbl,foo), \c READ_WFS(dbl,foo), \c READ_DM(dbl,foo) for array
   parameters that have one value per powfs, WFS, or DM.
3. Derive dependent values and validate the input in the appropriate
   \c setup_parms_postproc_*() function.
4. Add the key with its default value and a comment to the appropriate
   \c .conf file, for example \c sim.conf, \c recon.conf, \c dbg.conf, or
   \c powfs_common.conf.
5. Free any allocated memory in free_parms().

The configuration parser (\c lib/readcfg.c) checks the keys when it is closed
in setup_parms(): a key that is supplied by the user (on the command line or
in a user \c .conf file) but never read is an error, and so is a key that is
read more than once. Unused keys in the default configuration files only
produce a warning. Set \c MAOS_PERMISSIVE=1 to turn these errors into
warnings. A key that is read but not found produces a warning and is assumed
to be 0, so always provide a default in the \c .conf files. Conf files
\c include= each other; later values override earlier ones, and \c key+=value
appends to an array.

\section arch-thread Threading

\c maos parallelizes both across the main tasks of a time step (see above)
and inside each task. The macros in \c sys/thread.h provide a common
interface:

- \c QUEUE(group, fun, arg, njob, urgent) queues \c njob calls of
  \c fun(arg); \c WAIT(group, urgent) waits for the queued jobs in the group;
  \c CALL() queues and waits.
- \c QUEUE_THREAD(), \c CALL_THREAD() run an array of \c thread_t, which
  divides an index range among the threads (see \c thread_prep()).
- \c OMP_FOR(), \c OMP_TASK_FOR() parallelize loops with OpenMP.

When compiled with OpenMP 4.5 or newer, these macros use OpenMP tasks.
Otherwise they use the thread pool in \c sys/thread_pool.c, which keeps a
fixed number of threads that take jobs from a queue; a thread that waits for
its jobs runs other pending jobs in the meantime. main() runs maos_sim()
inside an OpenMP \c parallel/\c single region (\c OMPTASK_SINGLE) so that the
tasks can be scheduled.

\section arch-mem Memory Accounting

\c sys/mem.h replaces \c malloc, \c calloc, \c realloc, and \c free (in C code)
with wrappers that record every allocation when memory debugging is enabled.
It is enabled by default in debug builds (\c --enable-debug) and can be set
with the environment variable \c MAOS_MEM_DEBUG=1 or 0. At exit, the
unreleased memory is printed together with the backtrace of its allocation.
\c MAOS_MEM_VERBOSE=1 prints every allocation and free. Use \c mycalloc(n,
type) and \c mymalloc(n, type) to allocate typed memory. See
\ref page95_profiling.

\section arch-log Logging

The logging macros are defined in \c sys/common.h and filtered by the
environment variable \c MAOS_LOG_LEVEL (default 0, range -5 to 5). The output
goes to the terminal (unless detached) and to the log file.

- \c info3(), \c info2(), \c info(): information of decreasing importance,
  shown at the default level.
- \c dbg(), \c dbg2(), \c dbg3(): debugging information, shown only when
  \c MAOS_LOG_LEVEL is at least 1, 2, or 3, respectively.
- \c warning(): a warning with the file and line number.
- \c error(): prints the message with the file and line number and raises
  \c SIGUSR2. The signal handler (\c sys/mem.c) prints a backtrace and
  calls the handler registered by the program, which normally terminates it.
- Variants such as \c info_once(), \c warning_once(), \c dbg_time(), and
  \c info_time() print only once or add a time stamp.

\section arch-scheduler Scheduler, Monitor, and Drawdaemon

- The \c scheduler (\c tools/scheduler.c) is a per user daemon that queues
  jobs and keeps their status. It listens on a TCP port that is derived from
  the user name or read from \c ~/.aos/port. When compiled with
  libwebsockets, it also serves the browser monitor (\c monitor.html) on that
  port plus 100. \c maos launches the scheduler if it is not running.
- \c maos talks to the scheduler through \c sys/scheduler_client.c:
  scheduler_report_path() reports the job, scheduler_start() waits until
  enough CPUs are available (the \c -f option skips the wait),
  scheduler_report() sends the status every step, and scheduler_finish()
  reports the end. scheduler_listen() starts a thread in \c maos that handles
  requests forwarded by the scheduler, such as a request to start drawing.
- The \c monitor (\c tools/monitor*.c, or \c monitor.html in a browser)
  connects to the schedulers on the hosts listed in \c ~/.aos/hosts, shows the
  jobs, and can kill jobs or request plots.
- Plotting uses \c lib/draw.c. When plotting is requested (\c plot.* options or
  from the monitor), \c maos sends the data over a socket to \c drawdaemon
  (\c tools/drawdaemon*.c), which is either launched locally or connected
  through the monitor and the scheduler.

Configure with \c --disable-scheduler (implied by \c --enable-release) to build
without the scheduler.
