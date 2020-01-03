/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \page simulation Configuring Simulation
   \section  Misregistration
   

   It is tricky to implement the misregistration of the deformable mirror and the WFS.
   
   \subsection sect-mis-dm Misregistration of the deformable mirror.  
   
   Denote the grid of actuators as saloc, and with misregistration, it becomes
   salocm. 
   
   - We should use saloc for DM fitting and pseudo open loop gradient
   computations because we are not aware of the misregistration. 

   - We should use salocm for raytracing to WFS and performance evaluation
   because this is what happens in the system, no matter whether we know the
   misregistration or not.

   \subsection sect-mis-wfs Misregistration of the WFS.
   
   Denote the grid of the WFS entry pupil as loc, and with misregistration, it
   becomes locm. Project the telescope amplitude map onto loc, we get amp, which
   is what we know. Project the telescope amplitude map onto locm, we get ampm,
   which is what happens but we don't know.

   - We should use loc and amp to compute the reconstructor.

   - We should use locm and ampm for tracing to WFS and performance evaluation.
   
   - We should use loc and ampm for matched filter, averaging gradient, or
   zernike tip/tilt operator computation, since this is the process that happens
   in the system, and the matched filter is updated through dithering to the
   true amplitude map. We use loc because it is our model within the WFS.
   
   - We should use ptsm->area (real area) to select subapertures for use in
     reconstruction, because in real system, the subapertures are selected by
     their illumination level, and is affected by the real amplitude.

  */
