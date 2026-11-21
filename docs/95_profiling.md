\page page95_profiling Debugging

\section debug-build Debug Build

Configure a separate build directory with \c --enable-debug to compile without
optimization (\c -g3 \c -O0) and with the \c DEBUG macro set to 1. In a debug
build:

- \c P() checks that array indices are within bounds and calls error()
  otherwise (see \ref sect-maos-math).
- Memory accounting is enabled by default (see below).

\section debug-tools Debugging Tools

The GNU \c gdb can be used for basic troubleshooting. error() raises
\c SIGUSR2, so \c gdb stops at the location of the error. When the program
crashes or calls error(), the signal handler prints a backtrace. If \c gdb is
not able to find the problem, \c valgrind is good at detecting out of bounds
memory access and use after free.

For CUDA related debugging, use \c cuda-gdb.

Increase the verbosity with the environment variable \c MAOS_LOG_LEVEL (from
-5 to 5, default 0). Levels 1 to 3 enable the \c dbg(), \c dbg2(), and
\c dbg3() output.

\section debug-mem Memory Accounting

Set the environment variable \c MAOS_MEM_DEBUG=1 before launching \c maos to
track memory allocation and free. Unreleased memory is printed at exit
together with the backtrace of its allocation, which helps to find memory
leaks. Memory accounting is enabled by default in debug builds; set
\c MAOS_MEM_DEBUG=0 to disable it. Related options:

- \c MAOS_MEM_VERBOSE=1 prints every allocation and free.
- \c MAOS_MEM_FUNTRACE=1 records the file, line, and function of the
  allocation instead of the backtrace.
- Sending \c SIGUSR1 to a running process (<code>kill -USR1 pid</code>)
  enables memory accounting from that point on.

The accounting covers memory allocated through \c malloc, \c calloc,
\c realloc, and \c free in C code, which \c sys/mem.h redirects to the
accounting wrappers.

\section debug-coverage Code Coverage

Configure with \c --enable-gcov (which also enables \c --enable-debug) to
compile with gcov instrumentation. Then:

- \c make \c check runs the unit tests. If \c lcov and \c genhtml are
  installed, the \c check-local target in \c test/Makefile.am collects the
  coverage of the whole build and writes an HTML report to
  \c test/coverage/ in the build directory.
- \c make \c doc runs the unit tests from a clean coverage state, and
  \c test/unittest_coverage.sh summarizes the line coverage of \c sys,
  \c math, and \c lib using \c gcov. The summary is included in the unit test
  page of the manual (\ref page96_unittest).
- To measure the coverage of a \c maos run instead, run \c maos and then
  \c scripts/gcov.sh at the root of the build directory, which uses \c lcov
  and \c genhtml to write the report to \c out/index.html.

\section debug-profiling Profiling

MAOS is multi-threaded, and the \c gprof based profiling does not generate
useful information. Instead, use a system wide sampling profiler, such as
\c perf, OProfile, or Intel VTune, to find the hotspots. Setting
\c MAOS_PARALLEL=0 disables running the main tasks of each time step in
parallel, which helps to obtain more accurate timing of each part, and
\c MAOS_NO_WFS=1, \c MAOS_NO_EVL=1, or \c MAOS_NO_RECON=1 skip the WFS,
performance evaluation, or reconstruction. The option \c -n sets the number
of threads.

For example, with OProfile:

    operf ./maos [options]
    opreport -cgf -l maos | gprof2dot.py -f oprofile | dot -Tpng -o output.png

\c gprof2dot.py can be found in the \c scripts folder. With \c perf, use

    perf record -g ./maos [options]
    perf report
