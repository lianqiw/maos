\page page40_results Simulation Results

\tableofcontents

\section maosres RMS WFE

MAOS provides a few convenient routines to quickly process the simulation
results and obtain time averaged wavefront error. 

## Python

Use the following (see \c aolib.py):

    res,fds=maos_res(dir[,seeds] [,iframe] [,iframe2]) #gather field averaged WFE. [] indicate optional
    res,fds=maos_res(dir*) #gather field averaged WFE from multiple folders recursively
    res,fds=maos_res_each(dir*) #gather results for each field

## Matlab

Use the following:

    res=maos_res(dir[,seeds] [,range]) %gather field averaged WFE. [] indicate optional
    res,fds=maos_reseach(dir[, seeds] [,range]) #gather results for each field

# Plotting Results

The wavefront error time history can be processed using the built in tool 
`drawres` and plotted using the Gtk based `drawdaemon`. 

    drawres folder1 folder2 ...

If the specified folder contains subfolders with results, they will be plotted 
together. When using on a text terminal with `DISPLAY` not set, if there is
`monitor` running in a connected computer, it will launch `drawdaemon`.

# Reading Results

MAOS generates output in binary format `.bin` files (which maybe gzipped)
as described below. PSFs are written into `.fits` files with extensions.

## .bin file format
The data are saved to the \c .bin files in a simple modular way. For a simple
matrix (like double matrix or cell matrix), we first write a magic number (4
byte unsigned int) representing the data type, and then x and y dimension (two 8
byte unsigned int), and finally the data itself (if it is cell, recursively
write the data contained in the cell with the same method).

Each block of data may optionally be proceeded with a metadata header that
describes the data, for example, the sampling of a `loc_t` grid. The surf OPD
contains header information about the sampling, origin and height of the OPD
array. 

We use `.bin` format instead of `.fits` because it can be conveniently
memory mapped for best efficiency in reading or saving.

For implementation details, please refer to \ref sect-bin.

## MATLAB

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
executable \c bin2fits that can convert \c .bin to \c .fits files.

# Result Files

There will be several files created during simulation in the result folder and
subfolder `extra` (those results are not saved when \c save.extra=0). The number
after underscore _ is the seed. For example, with seed 1 the following files are
produced. Read in these files using provided mex function \c read in MATLAB.
Notice that the suffix \c .bin or \c .fits has been omitted. You do not need to
add the suffix when use \c read.
## Wavefront error
- \c Res_1: A binary file containing a cell array that include the main results. 
  - Use res=read('Res_1'); 
  - \c res[0] contains the open loop wavefront variance (WFV in units of \f$m^2\f$) in row vectors. When sim.rmax=1, the rows are 
      - Piston removed WFV
      - WFV in tip/tilt modes
      - Piston/tip/tilt removed WFV
       
    When sim.rmax>1, the rows are
      - Piston removed WFV
      - Piston/Tip removed WFE
      - Piston/Tip/Tilt removed WFE
      - ...
  - \c res[1] no longer used.
  - \c res[2] contains the closed loop wavefront variance in row vectors. The format is the same as res[0]
  - \c res[3] (Only in split tomography) contains the closed loop wavefront variance. The rows are (not depend on sim.rmax)
      - WFV in LGS contains modes
      - WFV in NGS Tip/Tilt modes
      - WFV in NGS modes (including Tip/Tilt and additional modes controlled by NGS (On-Instrument WFS))

- \c Resp_1: \c resp: Results for each direction. Combines the following four files previous saved individually: Resolmp_1, Resclmp_1, Resolep_1, Resclep_1. 
    - Use \c resp=read('Resp_1)
    - \c resp[0,idir]: Open loop wavefront Zernike modes, for evaluation direction idir,
      defined on un-normalized coordinate on the aperture, in the Noll's Zernike order of piston/tip/tilt ...

    - \c resp[1,idir]: Close loop wavefront Zernike modes, in the same format as \c resp[0,idir].

    - \c resp[2,idir]: Open loop wavefront variance for each science field point in the same format as res[0].

    - \c resp[3,idir]: Closed loop wavefront variance for each science field point in the same format as res[0].


## Split tomography

- \c Resclemp_1: LGS/TT/NGS mode wavefront error for each direction.

- \c RescleNGSm_1: contains a row vector array of
  either the residual NGS modes (in radian like unit) or the applied NGS mode
  correction if ideal NGS mode correction is made. Mainly used to save the ideal
  NGS modes applied in skycoverage pre-simulation.

- \c RescleNGSmp_1: This is NGS mode before averaging over directions.

- \c ResoleNGSm_1 and \c ResoleNGSmp_1: Like above, but for open loop.

## Log files
- \c maos_[HOST]_[PID].conf: The final effective configuration for MAOS run with
  pid PID. Can be used to reproduce the simulation or check the configurations.
  Lines begin with # are defaults configurations while the rest have been
  overridden. File for most recent run is linked to \c maos_recent.conf

- \c run_[HOST]_[PID].log: Log file. File for most recent run is linked to \c
  run_recent.log

- \c sanea_sim_1: When wavefront sensors are running in physical optics mode,
  the average gradient measurement error for that wavefront sensor is saved (in
  order) in a cell in this file. Each cell is a column vector with elements
  twice the number of subaperture. The first half is for \f$x\f$ (or radial)
  gradient, and the second half is for \f$y\f$ (or azimuthal) gradient. They are
  in units of \f$rad^2\f$.

## PSF

The PSF total FOV is \c evl.wvl/evl.dx but may be trimmed by setting \c
evl.psfsize to a smaller number (of pixels). The PSF sampling is usually 
`evl.wvl/(2*aper.d)` but may be `evl.wvl/(evl.psfgridsize*evl.dx)` if \c
evl.psfgridsize is set. The exact number of those can be checked in the fits
header.

- \c evlpsfcl_1: When \c evl.psfmean is 1, contains the time averaged science
  closed loop psf. if is a cell array of \f$n_{wvl}\times n_{evl}\f$. Each cell
  contains a center cut of the science PSF. 

- \c evlpsfdl_1: diffraction limited, PSF.

- \c evlpsfol_1: frame and field averaged open loop science PSF.

- \c evlpsfhist_1_ievlx: When \c evl.psfhist is 1, each of these files contains
the time history of the complex PSF of evaluation direction \f$x\f$.

## Other

- \c Resfsmcmd_1: Each cell contains the FSM tip/tilt command (only none empty
  for LGS WFS) history in unit of radian.

- \c Resfsmerr_1: Each cell contains the FSM tip/tilt error history in unit of
  radian.

- \c psdcl_1 and \c psdol_1: Closed loop and pseudo open loop PSD for high order
  control obtained from telemetry. 

- `ResCG_1`: Tomography and DM Fitting conjugate gradient residual time history.

- `Reszoompos_1`: For LGS only, the trombine position. One cell per WFS.

\section geometry Geometry Data

When \c save.setup=1, MAOS will save the geometry data created before
entering simulation to files in folder \c setup. When \c save.recon=1, the
tomography and DM fitting matrices are also saved, which takes space and are
not enabled by \c save.setup=1.

The following explains each file. Not all of them may exist, depending on set
of options. The suffix of the file names are removed.


File name | Description
---------|--------------
\c actslave  |The slaving operator
\c ahst_GM   |ad hoc split tomography, model to gradient operator
\c ahst_MCC  |ahst, model cross coupling matrix.
\c ahst_Mdm  |ahst, the NGS mode defined on DM grid
\c ahst_Pngs |ahst, ngs mode removal operator from DM commands
\c ahst_Ptt  |ahst, tip/tilt removal operator from DM commands
\c ahst_Rngs |ahst, ngs mode reconstructor
\c ahst_Wa   |ahst, the weighting using science field.
\c aloc      |DM actuator grid.
\c ploc      |The coarse sampled grid on aperture for reconstruction
\c xloc      |The tomography grid.
\c TT        |The global tip/tilt modes for wfs
\c PTT       |The global tip/tilt removal operator for wfs
\c saneai    |The inverse of nea^2 used in tomography (rad^2)
\c W0        |The W0 weighting defined on ploc. 
\c W1        |The W1 weighting defined on ploc.
\c aper_locs |The fine sampled grid on telescope aperture.
\c aper_amp  |Telescope aperture amplitude defined on aper_locs
\c aper_mcc  |modal cross-coupling matrix for modes defined on aper_locs
\c FLM       |Fit left hand side operator, sparse matrix
\c FLU       |Fit left hand side operator, Low rank U matrix
\c FLV       |Fit left hand side operator, Low rank V matrix
\c FRM       |Fit right hand side operator, sparse matrix
\c FRU       |Fit right hand side operator, Low rank U matrix
\c FRV       |Fit right hand side operator, Low rank V matrix
\c GA        |Gradient operator from aloc.
\c GP        |Gradient operator from ploc.
\c GX        |Gradient operator from xloc.
\c HA        |Ray tracing from aloc to ploc along fit directions
\c HXF       |Ray tracing from xloc to ploc along fit directions
\c HXW       |Ray tracing from xloc to ploc along wfs directions
\c L2        |Laplacian regularization on xloc
\c NW        |Low rank terms in fitting
\c NW2       |Adjusted low rank terms by slaving.
\c powfs0_area         |The subaperture area
\c powfs0_dtf0_nominal |The nominal of DTF
\c powfs0_dtf0_si      |The si of DTF
\c powfs0_etfprep0_2d  |The elongation transfer function
\c powfs0_GP           |The gradient operator from ploc.
\c powfs0_GS0          |The gradient operator from aper_locs.
\c powfs0_i0           |The subaperture time averaged image for matched filter
\c powfs0_gx           |The pixel by pixel x gradient of i0
\c powfs0_gy           |The pixel by pixel y gradient of i0
\c powfs0_imcc         |The inverse of model cross-coupling matrix of piston/tip/tilt modes
\c powfs0_loc          |The grid for all subapertures (grouped by subapertures)
\c powfs0_amp          |The amplitude defined on powfs0_loc
\c powfs0_llt_loc      |The aperture grid of the uplink laser launch telescope (LLT)
\c powfs0_llt_amp      |The aperture amplitude of LLT defined on powfs0_llt_loc
\c powfs0_llt_imcc     |The inverse of model cross-coupling matrix of p/t/t modes for LLT
\c powfs0_srot         |The orientation of each subaperture wrt LLT
\c powfs0_srsa         |The distance of each subaperture from the LLT
\c powfs0_mtche        |The matched filter gradient estimator
\c powfs0_pts          |The lower left grid point of each subaperture.
\c powfs0_saloc        |The lower left corner of each subaperture
\c powfs0_sanea        |The subaperture noise equivalent angle(nea) in rad^2
\c powfs0_saneaxy      |The nea along x/y directions in rad^2
\c powfs0_saneaixy     |The inverse of nea along x/y directions in rad^-2
\c powfs0_saneara      |The subaperture nea along r/a directions in rad^-2
\c powfs0_sepsf        |The subaperture short exposure psf (tip/tilt removed)
\c powfs0_sodium       |The sodium layer profile
\c powfs0_sprint       |which subapertures to print
\c powfs0_GS0          |The averaging gradient operator from pts.
\c powfs1_ZS0          |The ztilt best fit gradient operator from pts.

\section telemetry Telemetry Data

Depending on parameters enabled, the simulation telemetry data will be saved
to files in the simulation folder. The following describes them in detail.
Notice that save.grad, save.gradgeom, save.ints and save.wfsopd takes values
of 1, 2 or 3, where 1 means saving for all wfs, 2 means saving for only high
order wfs and 3 means saving for lower order wfs.

Name convention: the last number after underscore is the seed. The following
shows results for seed1. When there are multiple wfs or science fields, we
only show the results for the first. 

Many of the results contains a cell array of length \c nsim (number of
simulation steps), whenever that makes sense.

The data regarding the DM commands are all defined on the actuator grid \c aloc.

The second column of the tables shows which option in \c dbg.conf enables the
save of this data.

The suffix of the file names are removed.


File name          |Option to enable | Description
-------------------|-----------------|-------------
\c atm_1           |save.atm   |The atmosphere 
\c dmerr_hi_1      |save.dm   |The DM error signal for high order wfs
\c dmfit_hi_1      |save.dm   |The DM fit result for high order wfs
\c Merr_lo_1       |save.dm   |The low order mode error signal (split tomography)
\c Mint_lo_1       |save.dm   |The low order mode integrator output (split tomography)
\c dmint_hi_1      |save.dm   |The DM integrator output of high order wfs output (split integrator only)
\c dmint_1         |save.dm   |The DM integrator output (command integrator for both high/low wfs)
\c dmreal_1        |save.dm   |The real, effective DM commands applied
\c dmpttr_1        |save.dmpttr   |The piston/tip/tilt removed, effective, DM commands.
\c evl0_opdcl_1    |save.evlopd   |The closed loop opd for science fields
\c evl0_opdol_1    |save.evlopd   |The open loop opd for science fields
\c gcov_wfs0_5_10_1      |save.gcov   |The covariance between gradients of wfs 0 and 5 saved at time step 10
\c opdr_1          |save.opdr   |The tomography output, defined on xloc
\c opdx_1          |save.opdx   |The atmosphere projected onto xloc (direct fitting)
\c wfs0_gradcl_1   |save.grad  |The wfs gradient measurement.
\c wfs0_gradnf_1   |save.grad  |The wfs noise free gradient.
\c wfs0_gradol_1 |save.grad  |The wfs pseudo open loop gradient
\c wfs0_gradgeom_1 |save.gradgeom   |The wfs geometric gradient (in physical optics wfs simulations)
\c wfs0_intsnf_1   |save.ints   |The wfs subaperture images (noise free)
\c wfs0_intsny_1   |save.ints   |The wfs subaperture images (noisy)
\c wfs0_lltopd_1   |save.wfsopd   |The wfs laser launch telescope OPD.
\c wfs0_opd_1      |save.wfsopd   |The wfs OPD.
------
