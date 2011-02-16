/**
   \page page11 Misregistration

   It is tricky to implement the misregistration of the deformable mirror and the WFS.
   
   \section sect-mis-dm Misregistration of the deformable mirror.  
   
   Denote the grid of actuators as saloc, and with misregistration, it becomes
   salocm. 
   
   - We should use saloc for DM fitting and pseudo open loop gradient
   computations because we are not aware of the misregistration. 

   - We should use salocm for raytracing to WFS and performance evaluation
   because this is what happens in the system, no matter whether we know the
   misregistration or not.

   \section sect-mis-wfs Misregistration of the WFS.
   
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
