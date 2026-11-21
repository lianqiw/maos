Overview {#mainpage}
===

The Multithreaded Adaptive Optics Simulator (MAOS) is an end-to-end adaptive
optics (AO) simulator. It can simulate many kinds of astronomical AO systems,
including single conjugate AO (SCAO), multi-conjugate AO (MCAO), laser
tomography AO (LTAO), multi-object AO (MOAO), and ground layer AO (GLAO). It is
written in C, runs on multi-core CPUs, and can optionally use Nvidia GPUs.

# User guide

- \ref page10_intro &ndash; what MAOS models and how it is organized.
- \ref page20_compile &ndash; download, dependencies, configure options, GPU and
  MATLAB support, the job monitor, and the Python scripts.
- \ref page30_run &ndash; command line options, configuration files, and environment
  variables.
- \ref page33_example &ndash; walk-through of the TMT NFIRAOS MCAO and SCAO example
  configurations.
- \ref page34_conf &ndash; the complete list of configuration files.
- \ref page40_results &ndash; the output files and how to read and plot them.
- \ref page43_nfiraos &ndash; performance of the TMT NFIRAOS baseline.
- \ref algorithm &ndash; notes on DM influence functions, DM hysteresis, and physical
  optics beam propagation.
- \ref skycoverage &ndash; the sky coverage post-processor `skyc`.

# Developer guide

- \ref page90_devel
  - \ref page91_basics &ndash; coding conventions and the basic data types.
  - \ref page93_maos &ndash; the structure of the `maos` simulator.
  - \ref page95_profiling &ndash; debugging and profiling.
  - \ref page96_unittest &ndash; the latest unit test results.

# Publications

- <a href="https://github.com/lianqiw/files/blob/master/Wang_AO4ELT2_MAOS.pdf">AO4ELT2 paper on MAOS</a>
- <a href="https://github.com/lianqiw/files/blob/master/Wang_SPIE2012_GPU.pdf">SPIE 2012 paper on MAOS and a GPU based RTC</a>
- More references are on <a href="https://www.researchgate.net/profile/Lianqi-Wang">ResearchGate</a>.

# Historical benchmarks

Time per simulation step for the TMT NFIRAOS baseline (30 m aperture, 6 LGS,
dual order 60x60 DMs), measured on hardware of the time:

- Dual Intel Xeon W5590 at 3.33 GHz: 1.3 s with 8 threads.
- Intel Core i7-2600 at 3.4 GHz: 1.9 s with 4 threads.
- Intel Core i7-2600 at 3.4 GHz: 1.67 s with 8 hyper-threads (4 physical cores).
- One Nvidia GTX 580 GPU: 0.2 s.
- Eight Nvidia GTX 580 GPUs: 0.03 s.

# Source code and author

The source code is on GitHub: <a href="https://github.com/lianqiw/maos">lianqiw/maos</a>.
This manual is at <a href="https://lianqiw.github.io/maos">lianqiw.github.io/maos</a>.

\author Lianqi Wang <a href="https://www.linkedin.com/in/lianqiw/">(LinkedIn)</a> at <a href="https://www.tmt.org">TMT International Observatory</a>.
