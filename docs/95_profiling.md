\page page95_profiling Debugging

\section Debugging

The GNU gdb can be used for basic bug troubleshooting. In case it is not able
to detect the problem, \c valgrind can be used which is good at
detecting memory out of bound access or use after free.

Set macro MAOS_MEM_DEBUG=1 before launching maos will also enable it to track
the memory allcation/free and print unrealeased memory in the end that can be
used to detector any potential memory leakage.

For CUDA related debugging, use cuda-gdb.

\section Code Coverage

Compile maos with --enable-gcov to enable code coverage reporting.

\section Profiling

MAOS is a pthread enabled software, and in the current stage, the gprof based
profiling does not generate useful information. Instead, we can use system
wide profilers, such as the Oprofile or Intel vtune Amplifier to detect
hotspots and optimize the code. Set env MAOS_PARALLEL=0 disables internal
parallelization that may help to obtain more accurate timing.

Prepare oprofile:

    opcontrol --callgraph=16
    opcontrol --start

Then Run MAOS. To view profile:

    opcontrol --stop
    opcontrol --dump
    opreport -cgf -l maos | gprof2dot.py -f oprofile | dot -Tpng -o output.png

or 

    opgprof maos | gprof2dot.py | dot -Tpng -o output.png

    \c gprof2dot.py can be found in the scripts folder.

