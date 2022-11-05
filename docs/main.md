Overview {#mainpage}
===

Multi-Thread Adaptive Optics Simulator (MAOS) is an end-to-end adaptive
optics system simulator that is capable to simulate many kinds of
astronominal adaptive optics systems, including conventional single
conjugate AO (SCAO), multi-conjugate AO (MCAO), laser tomography AO (LTAO),
multi-object AO (MOAO), and ground layer AO (GLAO). Please follow the
following links for relevant information.

<p>

- \ref page10_intro 
- \ref page20_compile 
- \ref page30_run
- \ref page33_example
- \ref page40_results
- \ref page43_nfiraos 
- \ref algorithm
- \ref page90_devel

<a href="https://github.com/downloads/lianqiw/files/maos_gpu.pdf">AO4ELT2 Paper on MAOS</a>

<p>

Some benchmarking results using the TMT NFIRAOS (30 m aperture, 6 LGS, dual order 60x60 DM):
- Dual Intel Xeon W5590 at 3.33 Ghz: 1.3 s per time step with 8 threads.
- Intel Core i7-2600 at 3.4 Ghz: 1.9s per time step with 4 threads.
- Intel Core i7-2600 at 3.4 Ghz: 1.67s per time step with 8 hyper threads (4 physical cores).
- Nvidia GTX 580 GPU: 0.2s per time step
- 8x Nvidia GTX 580 GPU: 0.03s per time step


\author Lianqi Wang at TMT International Observatory www.tmt.org 

The source code can be obtained in <a href="http://github.com/lianqiw/maos">http://github.com/lianqiw/maos </a>
