\page page40_results Simulation Results

\tableofcontents

\section drawres Plotting Results

The wavefront error time history can be processed using the built in tool 
`drawres` and plotted using the Gtk based `drawdaemon`. 

  drawres folder1 folder2 ...

If the specified folder contains subfolders with results, they will be plotted 
together. When using on a text terminal with `DISPLAY` not set, if there is
`monitor` running in a connected computer, it will launch `drawdaemon`.

\section bin Reading Results

MAOS generates output in binary format \c .bin files (which maybe gzipped)
as described below. PSFs are written into fits files with extensions.

The data are saved to the \c .bin files in a modular way. For a simple matrix
(like double matrix or cell matrix), we first write a magic number (4 byte
unsigned int) representing the data type, and then x and y dimension (two 8
byte unsigned int), and finally the data itself (if it is cell, recursively
write the data contained in the cell with the same method).

Each block of data may optionally be proceeded with a metadata header that
describes the data, for example, the sampling of a loc_t grid. The surf OPD
contains header information about the sampling, origin and height of the OPD
array. 

We use \c .bin format instead of \c .fits because it can be conveniently
memory mapped for best efficiency in reading or saving.

\subsection read-matlab MATLAB

There are two MATLAB \c mex routines \c read and \c write
that can read and write \c .bin or \c .fits files in MATLAB. The source of
these mex routines are located in source sub-folder \c mex. If \c mex is found in
the system, this mex routines will be compiled automatically in the \c mex
subfolder of the compiling folder. If it doesn't get compiled, please goto
the \c mex subfolder in the source folder and type ./compile.sh to compile
them. Copy \c write.mexa64 and \c read.mexa64 into a folder that is in your
matlab path, such as \c $HOME/matlab. The usage in \c MATLAB is as follows:

``` 
cle=read('Res_1');     #will try .bin or .fits
cle=read('Res_1.bin'); #only use .bin

write(cle,'Res_1');     will write the data to \c Res_1.bin
write(cle,'Res_1.bin'); will write the data to \c Res_1.bin without compression.
```

If there is any problem in compiling \c read.c and \c write.c, there are also
two matlab scripts in the \c scripts folder, \c readbin.m and \c writebin.m,
that have most of functionality of \c read.c, \c write.c albeit a bit slower
and cannot handle compressed files conveniently.

\subsection sect-python Python

There is \c readbin.py in \c scripts folder for reading \c .bin or \c .fits
files in python. Alternatively, it is also possible to call \c C routines in
python, including \c readbin and \c writebin as \c read and \c write. To do
so, first set an environment variable \c MAOS_AOLIB poing to the \c aolib.so
file in the \c .bin folder in the compiling directory. Then import \c
scripts/aoslib.py to your python.

\subsection sect-idl IDL

There are two idl scripts in the sub-folder \c script. They are \c readbin.pro
for reading \c .bin files, and \c writebin.pro for writing \c .bin files.

\subsection sect-fits FITS

Certain results (e.g., PSF) are saved in \c .fits format. Any standard fits
reader (e.g., \c fv, \c ds9) should be able to read them. There is also an
executable \c bin2fits that can convert \c .bin to \c. fits files.

\section sect-interpret Result Files

There will be several files created during simulation in the result
folder. The number after underscore _ is the seed. For example, with seed 1
the following files are produced. Read in these files using provided mex
function \c read in MATLAB. Notice that the suffix \c .bin or \c .bin.gz has
been removed. You do not need to add the suffix when use \c read.

- \c Res_1: A binary file containing a cell array that include the main
results. i.e. res=read('Res_1'); 
  - \c res{1} contains the open loop wavefront variance (WFV in units of \f$m^2\f$) in row vectors. The rows are 
      - Piston removed WFV
      - WFV in tip/tilt modes
      - Piston/tip/tilt removed WFV
  - \c res{2} contains the residual wavefront variance after the tomography phase screen is directly applied as the correction in open loop mode. Empty if evl.tomo is zero.
  - \c res{3} contains the closed loop wavefront variance in row vectors. The rows are
      - Piston removed WFV
      - WVF in tip/tilt modes
      - Piston/tip/tilt removed WFV
  - \c res{4} (Only in split tomography) contains the closed loop wavefront variance. The rows are
      - WFV in LGS contains modes
      - WFV in NGS Tip/Tilt modes
      - WFV in NGS modes (including Tip/Tilt and additional modes controlled by NGS (On-Instrument WFS))

- \c Resolep_1: Open loop wavefront variance for each science field point. Each cell represent a science field point. The format is similar to res{1} above.

- \c Resolmp_1: Open loop wavefront Zernike (not normalized wrt radius)
modes defined on not-normalized coordinate on the aperture. The format is
similar as \c Resolep_1

- \c Resclep_1: Closed loop wavefront variance for each science field
point, in the same format as \c Resolep_1

- \c Resclmp_1: Close loop wavefront Zernike modes, in the same format as \c Resolmp_1

- \c Resclemp_1: (Only in split tomography) LGS/TT/NGS mode wavefront error
for each direction.

- \c RescleNGSm_1: (Only in split tomography) contains a row vector array
of either the residual NGS modes (in radian like unit) or the applied NGS
mode correction if ideal NGS mode correction is made. Mainly used to save the
ideal NGS modes applied in skycoverage pre-simulation.

- \c maos_[PID].conf: The final effective arrays of the configurations for MAOS
run with pid PID. Can be used to reproduce the simulation or check the
configurations. Lines begin with # are defaults configurations while the rest
have been overridden.

- \c sanea_sim_1: When wavefront sensors are running in physical optics
mode, the average gradient measurement error for that wavefront sensor is
saved (in order) in a cell in this file. Each cell is a column vector with
elements twice the number of subaperture. The first half is for \f$x\f$ (or radial)
gradient, and the second half is for \f$y\f$ (or azimuthal) gradient. They are in
units of \f$rad^2\f$.

- \c evlpsfcl_1: When \c evl.psfmean is 1, contains the time averaged
science closed loop psf. if is a cell array of \f$n_{wvl}\times
n_{evl}\f$. Each cell contains a center cut of the science PSF. 

The PSF total FOV is \c evl.wvl/evl.dx but may be trimmed by setting \c
evl.psfsize to a smaller number (of pixels). The PSF sampling is usually \c
evl.wvl/(2*aper.d) but may be \c evl.wvl/(evl.psfgridsize*evl.dx) if \c
evl.psfgridsize is set. The exact number of those can be checked in the fits
header.


- \c evlpsfdl_1: diffraction limited, PSF.

- \c evlpsfol_1: frame and field averaged open loop science PSF.

- \c evlpsfhist_1_ievlx: When \c evl.psfhist is 1, each of these files contains
the time history of the complex PSF of evaluation direction \f$x\f$.


- \c Resuptcmd_1: Each cell contains the uplink tip/tilt command (only none
empty for LGS WFS) history in unit of radian.

- \c Resupterr_1: Each cell contains the uplink tip/tilt error history in
unit of radian.

