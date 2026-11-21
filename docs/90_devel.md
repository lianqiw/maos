\page page90_devel Development

This section is for developers who want to understand, modify, or extend the
code.

- \subpage page91_basics -- C fundamentals, MAOS data types, the \c bin file
  format, and coding guidelines.
- \subpage page93_maos -- library layering, the simulation flow, key data
  structures, threading, and how to add a parameter.
- \subpage page95_profiling -- debugging, memory checking, code coverage, and
  profiling.
- \subpage page96_unittest -- results of the unit tests and code coverage.

The physics and numerical algorithms are described in \ref algorithm.

\section devel-build Build Directories

The code uses autotools and is always built out of the source tree. Several
build directories configured with different options can coexist next to the
source directory, for example:

```
cd maos            #the source directory
./autogen.sh       #needed for a fresh clone or after editing configure.ac or Makefile.am
mkdir ../optim && cd ../optim && ../maos/configure                  #optimized build
mkdir ../debug && cd ../debug && ../maos/configure --enable-debug   #debug build
make -j8           #executables are placed in bin/ of the build directory
```

Useful \c configure options for development are \c --enable-debug,
\c --enable-gcov (code coverage, implies \c --enable-debug),
\c --disable-double (single precision), \c --enable-release (portable build
without the scheduler), and \c CC=clang++ (compile the C code as C++). See
\ref page20_compile for all the options.

\section devel-unittest Unit Tests

The unit tests are in \c test/cmocka_*.c and use the
<a href="https://cmocka.org">CMocka</a> framework. Each test includes
\c test/unittest.h, which provides helpers such as a double precision
\c assert_float_equal(). Tests must not read or write temporary files:
generate the input data in memory and check the results with assertions. To
add a test program, add it to \c bin_PROGRAMS in \c test/Makefile.am.

Run the tests from the build directory:

```
make check                                  #all tests (requires cmocka)
make -C test check TESTS=cmocka_mat         #a single test
test/cmocka_mat                             #run a test program directly
```

Each run updates the status report \c docs/unittest_status.md in the build
directory. \c make \c doc runs the unit tests, generates the report (with the
line coverage of \c sys, \c math, and \c lib when configured with
\c --enable-gcov), and builds the manual with doxygen, where the report
replaces the placeholder \ref page96_unittest.

\section devel-regression End-to-end Regression

\c scripts/maos_check.sh runs \c maos for 100 steps in many configurations
(NGS and LGS SCAO, GLAO, LTAO, MOAO, MCAO with various tomography algorithms,
pyramid WFS, NCPA calibration, and others) and compares the mean
wavefront error with reference values. Run it from the \c bin folder of a
build directory:

```
cd ../optim/bin
../../maos/scripts/maos_check.sh 5        #5 m telescope
../../maos/scripts/maos_check.sh 5 -g-1   #force CPU
```

The first argument is the telescope diameter: 30 uses \c tmt.conf, 10 uses
\c keck.conf, and other values set \c aper.d. The remaining arguments are
passed to \c maos. A configuration fails if the result differs from the
reference by more than about 10%; the script exits with the number of
failures. The reference values are the \c REF_GPU and \c REF_CPU arrays in the
script, one pair per diameter. When a change intentionally alters the results,
update these arrays: each run appends a line with the new
\c REF_GPU=(...) or \c REF_CPU=(...) array to \c maos_check.ref in the current
directory, which can be copied into the script. The results and logs are saved
to \c maos_check_*.* files.

\section devel-ci Continuous Integration

The GitHub workflows in \c .github/workflows run on every push and pull
request:

- \c linux.yml: builds with gcc (\c --enable-release with CUDA, without a
  GPU so that the simulation runs on the CPU), runs \c make \c check and
  \c make \c doc, publishes the documentation to the \c gh-pages branch,
  builds with clang and runs \c make \c check, and runs
  \c scripts/maos_check.sh \c 5.
- \c macos.yml: builds on macOS with \c --enable-release and runs
  \c make \c check.
- \c macos-cpp.yml: builds on macOS with \c CC=clang++ and runs
  \c make \c check. The code must therefore also compile as C++.

\section devel-workflow Contribution Checklist

- Follow the guidelines in \ref sect-guide.
- Add or update unit tests for changes in \c sys, \c math, and \c lib, and
  run \c make \c check.
- Run \c scripts/maos_check.sh for changes that affect the simulation, and
  update its reference values if the change is intended to alter the results.
- Document new parameters in the \c .conf files (see \ref arch-param) and
  new functions with doxygen comments. \c make \c doc builds the manual.
