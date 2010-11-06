/**
   \page page4 Simulation Results
   
   \section simu-results Simulation Results
   MAOS will generate output in binary format \c .bin or zipped \c .bin.gz files
   as described below. There are two MATLAB mex routines \c read and \c write
   that can read and write \c .bin or \c .bin.gz files in MATLAB. The source of
   these mex routines are located in sub-folderx \c mex. Cd to this folder and
   issue \c make to compile these mex routines.  Copy \c write.mexa64 and \c
   read.mexa64 into a folder that is in your matlab path, such as
   $HOME/matlab. The usage in MATLAB is as follows:

\verbatim 
   cle=read('Res_1'); 
or cle=read('Res_1.bin'); 
or cle=read('Res_1.bin.gz');
   each wil try to read Res_1.bin or Res_1.bin.gz 
   
   write(cle,'Res_1');     will write the data to Res_1.bin.gz;
   write(cle,'Res_1.bin'); will write the data to Res_1.bin without compression.
\endverbatim

There will be several files created during simulation in the result folder. The
number after underscore _ is the seed. For example, with seed 1 the following
files are produced. Read in these files using provided mex function \c read in
MATLAB. Notice that the suffix .bin or .bin.gz has been removed. You do not need
to add the suffix when use \c read.

 - \c Res_1: A binary file containing a cell array that include the main
results. i.e. res=read('Res_1'); 
     - \c res{1} contains the open loop wavefront variance (WFV in units of \f$m^2\f$) in row vectors. The rows are 
          - Piston removed WFV
          - WFV in tip/tilt modes
          - Piston/tip/tilt removed WFV
     - \c res{2} contains the residual wavefront variance after the tomography phase screen is directly applied as the correction in open loop mode. Empty is evl.tomo is zero.
     - \c res{3} contains the closed loop wavefront variance in row vectors. The rows are
          - Piston removed WFV
          - WVF in tip/tilt modes
          - Piston/tip/tilt removed WFV
     - \c res{4} (Only in split tomography) contains the closed loop wavefront variance. The rows are
          - WFV in LGS contains modes
          - WFV in NGS Tip/Tilt modes
          - WFV in NGS modes (includeing Tip/Tilt and additional modes controlled by NGS (On-Instrument WFS)

 - \c Resolep_1: Open loop wavefront variance for each science field point. Each cell represent a science field point. The format is similar to res{1} above.

 - \c Resolmp_1: Open loop wavefront Zernike (not normalized wrt radius)
   modes defined on not-normalized coordinate on the aperture. The format is
   similar as \c Resolep_1
 
 - \c Resclep_1: Closed loop wavefront variance for each science field
   point, in the same format as Resolep_1

 - \c Resclmp_1: Close loop wavefront Zernike modes, in the same format as
   Resolmp_1

 - \c Resclemp_1: (Only in split tomography) LGS/TT/NGS mode wavefront error
   for each direction.

 - \c RescleNGSm_1: (Only in split tomography) contains a row vector array
   of either the residual NGS modes (in radian like unit) or the applied NGS
   mode correction if ideal NGS mode correction is made. Mainly used to save the
   ideal NGS modes applied in skycoverage pre-simulation.

 - \c maos_975.conf: The final effective arrays of the configurations for MAOS
   run with PID 975. Can be used to reproduce the simulation or check the
   configurations.

 - \c sanea_sim_1: When wavefront sensors are running in physical optics
   mode, the average gradinet measurement error for that wavefront sensor is
   saved (in order) in a cell in this file. Each cell is a column vector with
   elements twice the number of subaperture. The first half is for x (or radial)
   gradient, and the second half is for y (or azimuthal) gradient. They are in
   units of \f$rad^2\f$.

 - \c evlpsfcl_1: When evl.psfmean is 1, contains the time averaged science
   closed loop psf. if is a cell array of \f$n_{wvl}\times n_{evl}\f$. Each cell
   contains a center cut of the science PSF.

 - \c evlpsfdl_1: diffraction limited, PSF.

 - \c evlpsfol_1: frame and field averaged open loop science PSF.

 - \c evlpsfhist_1_ievlx: When evl.psfhist is 1, each of these files
   contains the time history of the complex PSF of evaluation direction x.


 - \c Resuptcmd_1: Each cell contains the uplink tip/tilt command (only none
   empty for LGS WFS) history in unit of radian.

 - \c Resupterr_1: Each cell contains the uplink tip/tilt error history in
   unit of radian.


  

 */
