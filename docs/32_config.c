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
   \page page30_run Run simulations
   \tableofcontents
   \section advanced Advanced configuration
 
   \subsection sect-surface Specifying Surface OPDs

   We can optionally setup one or more static surfaces that cover science fields
   and/or wavefront sensors. Each surface file must contain a 2-d array of the
   OPD with a header specifying the following keys:

   \verbatim
   dx  the sampling in x direction (first dimension).
   dy  the sampling in y direction (second dimension).
   ox  the origin in x direction.
   oy  the origin in y direction.
   h   the height conjugation of this surface
   vx  the frozen flow wind speed along x of this surface (0 for static)
   vy  the frozen flow wind speed along y of this surface (0 for static)
   \endverbatim

   Optionally, the header can also include the following keys to indicate its coverage

   \verbatim
   SURFNAME #M1 or M2 for the primary or secondary mirror
   SURFEVL=[0 1 2] #covers science evaluation directions 0 1 2
   SURFWFS=[0 1 2] #covers WFS 0 1 2
   \endverbatim
    
   Use the \c write mex routine, \c writebin.m, or \c writebin.py to write the bin file:

       write(OPD, header, 'opd1.bin')
    
   Put the list of surface `.bin` files in key \c surf.

       maos surf=['opd1.bin','opd2.bin']

   The amplitude map of the telescope can be specified with
   `aper.fnamp=aper.bin` with a similar data format, with OPD replaced by
   amplitude map between 0 and 1.

   We can also setup one or more tilted M3 surfaces that are common to science
   fields and wavefront sensors. Put the list of surface files in key \c
   tsurf. Each surface file must contain a 2-d array with a header specifying
   the following keys:

   \verbatim
   dx    the sampling in x direction (first dimension).
   dy    the sampling in y direction (second dimension).
   ox    the origin in x direction.
   oy    the origin in y direction.
   txdeg the x tilt angle in degrees wrt beam (90 is perp)
   tydeg the y tilt angle in degrees wrt beam (90 is perp)
   ftel  the telescope effective focal length
   fexit the distance between the exit pupil and the focus
   fsurf the distance between the center of the M3 surface and the focus.
   \endverbatim

   \subsection sect-wfs WFS Configuration

   The WFS configuration are mostly supplied in \c powfs, which applies to all
   the wfs belonging to this type. For wfs specific parameters, use \c wfs. If
   there is a single number in \c powfs.[key], it applies to all powfs. The
   follow lists a few common options

   \verbatim
   powfs.noisy=[0]         #set noise free for all powfs
   powfs.noisy=[1 0 0]     #set noise free only for powfs 1 and 2.
   powfs.phystep=[0 -1 -1] #use physical optics wfs for the first type (high order), and geomtric for the remaining (t/t/f)
   powfs.phytype=[1 2 2]   #pixel processing 1: use matched filter, 2: cog
   powfs.siglev=[900 100 100] #specify signal level for physical optics wfs mode at 800 Hz per pixel.
   powfs.bkgrnd=[10 0 0]   #specify background level at 800 Hz per pixel
   powfs.rne=[3 1 1]       #specify read out noise for physical optics wfs mode
   powfs.nearecon=[20 2 1] #specify noise equivalent angle in milli-arcsecond for geometric wfs mode
   powfs.neasim=[20 2 1]   #specify nea for simulation. -1 to match nearecon.
   \endverbatim

   WFS specific parameters usually include WFS coordinate
   \verbatim
   wfs.thteax=[]           #specify the x coordinate in arcsec
   wfs.thetay=[]           #specify the y coordinate in arcsec
   \endverbatim

   \subsection sect-perfevl Point Spread Function

   Maos only computes RMS averaged wavefront errors by default, which are saved
   to \c Res_[seed].bin. When desired, PSFs computing can be enabled for some or all
   of the science evaluation directions. Modify \c evl.thetax , \c evl.thetay , and \c
   evl.wt to modify evaluation directions and relative weighting for field
   averaging. Use the following parameters to enable PSF computation.
   \verbatim
   maos evl.psfisim=20 evl.psfmean=1 # start averaging from step 20, save averaged PSF once in the end.
   maos evl.psfmean=1 evl.psfol=1   # also include open loop PSF
   maos evl.psfmean=1 evl.psfsize=[1024,1024]     # only save center part 1024x1024 of the PSF
   maos evl.psfmean=1 evl.psfgridsize=[4096,4096] # specify PSF FFT grid size to control the sampling.
   \endverbatim

   \subsection sect-act Actuator Slaving

   Due to pupil obscuration, there may be actuators that are not sensed by the
   WFS and therefore cannot be accurately reconstructed. The actuator coupling
   coefficient, computed using DM fitting ray tracing (HA in Minimum variance
   reconstruction) or gradient interaction matrix (GA in least squares
   reconstruction), is used to determine how accurately an actuator is
   sensed/controlled. Activated by
   
       lsr.actslave=1 #or 2, for lsr
       fit.actslave=1 #or 2, for mvr
   
   - \c actslave=1 enables implementation with threshold \c .actthres.  When an
   actuator is outside of the pupil, its coupling coefficient is usually very
   small. When it is below lsr.acthres, its value is slaved to its neighboring
   actuators.
   - \c actslave=2 in addition enables implementation with threshold \c
   .actthres2.  When an actautor is hidden by secondary mirror support struts,
   it may still have relatively strong coupling, but regions separated by the
   structs may have different piston values, so called island effect. The
   mitigation method is to limit the value gap around such actautors whoes
   coupling is below .actthres2.

*/
