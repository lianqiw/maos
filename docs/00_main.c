/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
    \mainpage MAOS Documentation

    Multi-Thread Adaptive Optics Simulator (MAOS) is a end-to-end adaptive
    optics system simulator. It has the capability to simulate all kinds of
    astronominal adaptive optics systems, including conventional single
    conjugate AO (SCAO), multi-conjugate AO (MCAO), laser tomography AO (LTAO),
    multi-object AO (MOAO), and ground layer AO (GLAO).

    <p>

    - \ref page10 

    - \ref page20 

    - \ref page25

    - \ref simulation 

    - \ref page40 
    
    - \ref page41 
    
    - \ref page42 

    - \ref page50 

    - \ref page60 

    - \ref page70 

    <a href="https://github.com/downloads/lianqiw/files/maos_gpu.pdf">AO4ELT2 Paper on MAOS</a>

    <p>
    
    Some benchmarking results using the TMT NFIRAOS (30 m aperture, 6 LGS, dual order 60x60 DM):
    - Dual Intel Xeon W5590 at 3.33 Ghz: 1.3 s per time step with 8 threads.
    - Intel Core i7-2600 at 3.4 Ghz: 1.9s per time step with 4 threads.
    - Intel Core i7-2600 at 3.4 Ghz: 1.67s per time step with 8 hyper threads (4 physical cores).
    - Nvidia GTX 580 GPU: 0.2s per time step
    - 8x Nvidia GTX 580 GPU: 0.03s per time step
  

    \author Lianqi Wang <lianqiw@gmail.com> at TMT Corporation www.tmt.org 
    
    The source code can be obtained in <a href="http://github.com/lianqiw/maos">http://github.com/lianqiw/maos </a>
 */
