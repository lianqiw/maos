\page page40_results Simulation Results

\tableofcontents

MAOS writes all results to the output folder given by `-o`. Numerical results
use the `.bin` format described in \ref sect-bin; PSFs and science OPDs use
FITS. The number after the last underscore in a file name is the seed (for
example `Res_1.bin` is for seed 1). When reading a file with the tools below,
the suffix (`.bin`, `.fits`, with optional `.gz`) may be omitted.

# RMS WFE

## Python

`scripts/maos_result.py` contains routines that read `Res_*.bin` or
`Resp_*.bin` from one or more result folders and average the wavefront
variance over a range of time steps and over the seeds found in each folder.
Add the `scripts` folder to `PYTHONPATH` so that `readbin.py` can be imported.

    from maos_result import *
    res,fds=maos_res(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1)      #field averaged WFE
    res,fds=maos_res_tot(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1)  #high order + low order total
    res,fds=maos_res_hi(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1)   #high order only
    res,fds=maos_res_each(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1) #each evaluation direction

- `fds` is a folder name or a list of folder names. Wildcards are expanded
  with `glob`. Folders without valid results are skipped; the returned `fds`
  lists the folders that are used, one per row of `res`.
- `seeds` is a list of seeds. By default all `Res_*.bin` files in each folder
  are used.
- `iframe1` and `iframe2` select the time steps to average. A value less than
  1 (up to 1 for `iframe2`) is a fraction of the number of completed steps;
  a larger value is a step index. The default discards the first 20% of the steps. With
  `iframe1=0` and `iframe2=0` the time history is returned instead of the
  average.
- The results are wavefront variance in \f$nm^2\f$. Use `np.sqrt(res)` to get
  the RMS WFE in nm.
- `maos_res` returns one row per folder with columns: high order WFE,
  tip/tilt WFE, and low order total WFE (plus focus when split tomography
  includes a focus mode). For integrated tomography, the high order column is
  the piston/tip/tilt removed WFE and the last two columns are both tip/tilt.
- `maos_res_each` returns an array of dimension `nfolder*ndirection*nmod` with
  the closed loop piston removed, tip/tilt, and piston/tip/tilt removed WFE
  for each evaluation direction (the last part of `Resp_*.bin`).

## MATLAB

There is no MATLAB counterpart of `maos_result.py`. Read `Res_1.bin` with the
`read` mex routine (see \ref sect-mex) and process the cells as described in
\ref sect-res.

# Plotting Results

The wavefront error time history can be plotted with the built in tool
`drawres`, which sends the plots to the GTK based `drawdaemon`.

    drawres [-s seed] folder1 folder2 ...

- Without a folder, the current directory is used. Subfolders of each folder
  (except `skysim`) are searched recursively, and every folder that contains
  `Res_*.bin` is plotted together, so that results can be compared.
- `-s seed` limits the plots to the given seed. It may be repeated.
- The time axis is in seconds when `sim.dt` is found in `maos_done.conf` or
  `maos_recent.conf` of the folder, and in steps otherwise.
- The closed loop results (tab `CL`) are shown as total, high order, low
  order, tip/tilt, plate scale, and focus WFE (the last two only for split
  tomography) for each seed and averaged over the seeds. Open loop results
  are in tab `OL`. When there are several evaluation directions, the WFE as a
  function of field angle is also plotted from `Resp_*.bin`.
- Set environment variable `MAOS_DRAWRES_HI=1` (or call the program as
  `drawreshi`) to plot only the high order closed loop results.
- `drawres` also plots sky coverage results (`Res*_*.bin` written by `skyc`),
  but maos and skyc results cannot be mixed in one call.

`drawres` connects to an idle `drawdaemon` or asks the scheduler to launch
one. A `monitor` that is connected to the scheduler can display the plots.
Otherwise `drawdaemon` is launched locally when `DISPLAY` or
`WAYLAND_DISPLAY` is set.

# Reading Results

## .bin file format

Each array is saved as an optional header, a 4 byte magic number that
identifies the data type, two 8 byte unsigned integers for the `nx` and `ny`
dimensions, and the data in column major order (`x` changes fastest). A cell
array has its own magic number and dimensions, followed by each of its
elements written in the same way. The header is a text string with
`key=value` entries, for example the sampling of a `loc_t` grid or the
`sim.dt` and `evl.dx` of `Res_1.bin`.

Files are written uncompressed unless the file name ends with `.gz` (gzip) or
`.zst` (zstd). Readers detect the compression from the file content. FITS
files are also supported for reading and writing. We use `.bin` by default
because it is simple, supports cell arrays, and can be memory mapped.

For implementation details, please refer to \ref sect-bin.

\subsection sect-mex MATLAB

The source folder `mex` contains the mex routines `read` and `write` that
read and write `.bin` or `.fits` files in MATLAB. They are built in the
`mex` subfolder of the build folder if MATLAB is found by `configure` (see
`--with-matlab`). Alternatively, run `compile.sh` in the `mex` subfolder of
the source folder. Copy `read.mexa64` and `write.mexa64` to a folder that is
in your MATLAB path, such as `$HOME/matlab`. Usage:

```
res=read('Res_1');          %tries Res_1.bin, Res_1.bin.gz, and Res_1.fits
[res, header]=read('Res_1');%also returns the header
res=read('Res_1', n, start);%for a cell array, reads n cells starting at cell start (1 based)
write(res, 'Res_1');        %writes to Res_1.bin
write(res, header, 'Res_1');%writes with a header
```

The mex `write` compresses `.bin` files with gzip (keeping the `.bin`
suffix) and writes `.fits` files uncompressed. The mex `read` handles gzip
but not zstd compressed files.

If the mex routines cannot be compiled, the MATLAB scripts `readbin.m` and
`writebin.m` in the `scripts` folder provide most of the functionality,
albeit slower.

\subsection sect-python Python

`scripts/readbin.py` reads `.bin` and `.fits` files (uncompressed or
gzip) into numpy arrays:

    from readbin import readbin
    res=readbin('Res_1')   #tries Res_1.bin, Res_1.bin.gz, Res_1.fits, Res_1.fits.gz

Because numpy uses row major order, an array of size `nx*ny` in MAOS is
returned with shape `(ny, nx)`. A cell array is returned as a numpy object
array with the same convention (a cell with `ny=1` becomes a 1-d array). The
header of each array is available in its `header` attribute. A FITS file with
several images is returned as an array of the images.

Alternatively, the C routines in `aolib.so` (built in the `bin` folder of
the build folder) can be called from Python. Set the environment variable
`MAOS_AOLIB` to the path of `aolib.so` and import `scripts/libaos.py`. It
provides `read` and `write` that use the C implementation (with zstd
support) as well as many routines from `lib`.

\subsection sect-idl IDL

The `scripts` folder contains `readbin.pro` (with `readdata.pro`) for
reading `.bin` files in IDL. There is also a Julia reader `readbin.jl` that
uses `aolib.so`.

\subsection sect-fits FITS

Certain results (PSFs, science OPDs) are saved in `.fits` format. Any
standard FITS reader (e.g., `fv`, `ds9`, `astropy`) can read them. A cell
array is saved as a sequence of FITS images (extensions). The executable
`bin2fits` converts `.bin` files to `.fits` files.

# Result Files

The main results are saved in the output folder. Additional results are saved
in the subfolder `extra` when `save.extra=1` (the default in `dbg.conf`).
The file names below are for seed 1. Most of the results are saved
incrementally during the simulation, so they can be read while the simulation
is running; time steps that have not been simulated yet contain NaN or 0.

\subsection sect-res Wavefront error

- `Res_1.bin`: A cell array of 4 elements with the main results. The header
  contains `sim.dt` and `evl.dx`. Each element is an array with one column
  per time step. Values are wavefront variance in \f$m^2\f$, averaged over the
  science evaluation directions using the weights `evl.wt`.
  - `Res[0]` contains the open loop wavefront variance. When `evl.rmax=1`
    (the default), the rows are
      - Piston removed WFV
      - WFV in tip/tilt modes
      - Piston/tip/tilt removed WFV

    When `evl.rmax>1`, there are `(evl.rmax+1)*(evl.rmax+2)/2` rows. Row `i`
    is the WFV with the first `i+1` Zernike modes removed:
      - Piston removed WFV
      - Piston/tip removed WFV
      - Piston/tip/tilt removed WFV
      - ...
  - `Res[1]` is no longer used and is empty.
  - `Res[2]` contains the closed loop wavefront variance in the same format as
    `Res[0]`.
  - `Res[3]` (only when `evl.split` is set, e.g., in split tomography)
    contains the closed loop wavefront variance split into (independent of
    `evl.rmax`)
      - WFV in LGS controlled modes (total minus NGS modes)
      - WFV in NGS tip/tilt modes
      - WFV in all NGS modes (tip/tilt and the additional low order modes
        controlled by the NGS WFS)
      - WFV in NGS focus mode (only when focus is one of the NGS modes)

- `Resp_1.bin`: Results for each science evaluation direction. A cell array
  of `nevl` by 4. The header contains `sim.dt`, `evl.dx`, `evl.thetax`, and
  `evl.thetay` (in arcsec). In MATLAB, use `resp{idir,k}`; in Python use
  `resp[k,idir]`.
  - `k=0`: Open loop wavefront modes. When `evl.rmax=1`, the rows are
    piston (m), tip (rad), and tilt (rad). When `evl.rmax>1`, the rows are
    the coefficients (m) of the orthonormalized Zernike modes in Noll's order
    of piston/tip/tilt/...
  - `k=1`: Closed loop wavefront modes in the same format.
  - `k=2`: Open loop wavefront variance in the same format as `Res[0]`.
  - `k=3`: Closed loop wavefront variance in the same format as `Res[0]`.

## Split tomography results

These files are saved in `extra` when `evl.split` is set. The NGS modes are
the low order modes defined in `ahst_Modes` (tip/tilt first, then plate
scale, focus, etc. depending on the configuration).

- `Resclemp_1.bin`: For each direction, a 3 by `nsim` array of the closed
  loop WFV in LGS controlled modes, in NGS tip/tilt modes (field averaged),
  and in the NGS modes.

- `RescleNGSm_1.bin`: The closed loop NGS mode coefficients (field averaged),
  one column per time step. With `tomo.ahst_idealngs` set, it contains the
  NGS modes that are corrected ideally. It is the input of the sky coverage
  post processing (see \ref skycoverage).

- `RescleNGSmp_1.bin`: For each direction, the inner product of the closed
  loop OPD with the NGS modes (before averaging over directions).

- `ResoleNGSm_1.bin` and `ResoleNGSmp_1.bin`: Like above, but for open loop.

- `RescorrNGSm_1.bin`: The NGS mode correction (integrator output) applied
  (only when `recon.split=1` and `sim.fuseint=0`).

## Log files

- `maos_[HOST]_[PID].conf`: The final effective configuration of the MAOS
  process with pid PID. The first line includes the default configuration
  file. Lines beginning with # are default values while the rest have been
  overridden. It can be used to reproduce the simulation or check the
  configuration. When the simulation finishes successfully, it is renamed to
  `maos_done.conf` (an older one is renamed to `maos_done.1.conf`, etc.).
  `maos_recent.conf` links to the file of the most recent run.

- `run_[HOST]_[PID].log`: The log file. It is renamed to `run_done.log` when
  the simulation finishes successfully, or to `run_[HOST]_[PID].err` if it
  fails. `run_recent.log` links to the file of the most recent run.

- `Res_1.done`: Created when seed 1 finishes. MAOS skips seeds that are
  already done in the output folder unless `-O` (override) is given.

## PSF

PSFs are computed when `evl.psfmean` or `evl.psfhist` is set, starting from
time step `evl.psfisim`, for directions with `evl.psf` set. Each PSF has
FITS keywords that describe the field location, zenith angle, \f$r_0\f$,
\f$L_0\f$, wavelength, OPD sampling, FFT grid size, PSF sampling (`dp` in
arcsec), total intensity, and exposure time. The PSF sampling is
`evl.wvl/(nfft*evl.dx)`, where the FFT size `nfft` is set by
`evl.psfgridsize` or, by default, is about twice the aperture diameter, giving
a sampling of about `evl.wvl/(2*aper.d)`. The PSF size is `nfft` pixels unless
trimmed to `evl.psfsize` pixels. The suffix `x0_y0` is the field location
in arcsec (followed by the source height for sources not at infinity).

- `evlpsfcl_1_x0_y0.fits`: The closed loop science PSF, cumulatively averaged
  since `evl.psfisim`. Every `evl.psfmean` steps, one PSF per wavelength is
  appended, so the last `evl.nwvl` images are the average over the whole
  simulation.

- `evlpsfcl_ngsr_1_x0_y0.fits`: Like above, but with the NGS modes removed
  (`evl.psf=2` or 3 in split tomography).

- `evlpsfol_1.fits`: When `evl.psfol` is set, the open loop science PSF
  on axis (`evl.psfol=1`) or averaged over all directions (`evl.psfol=2`),
  in the same format.

- `evlpsfdl.fits`: The diffraction limited PSF, one per wavelength.

- `evlpsfhist_1_x0_y0.fits`: When `evl.psfhist` is set, the time history of
  the complex PSF.

- `evlopdcov_1_x0_y0.fits` and `evlopdmean_1_x0_y0.fits`: When `evl.cov` or
  `evl.opdmean` is set, the cumulative average of the covariance and mean of
  the closed loop science OPD. Open loop versions are `evlopdcovol_1.fits`
  and `evlopdmeanol_1.fits`.

## Subaperture Pixels

The subaperture pixels can be saved by specifying `save.setup=1` for models
computed during the setup process, and `save.ints=1` for runtime data. The
pixels are grouped by subapertures, with the subaperture ordering defined in
`saloc` (see `setup/powfs0_saloc.bin`).

- `setup/powfs0_i0.bin`: The WFS model image used to build the matched filter
  or calibrate CoG. `powfs0_gx.bin` and `powfs0_gy.bin` are its x and y
  derivatives.
- `wfs0_intsnf_1.bin`: The noise free WFS pixels for every frame.
- `wfs0_intsny_1.bin`: The noisy WFS pixels for every frame.

The config input `powfs.bkgrndfn` uses the same format.

## Other

These files are saved in `extra` when `save.extra=1`. Cell arrays with one
element per WFS have empty elements for WFS that do not use the feature.

- `Resfsmerr_1.bin` and `Resfsmcmd_1.bin`: For each LGS WFS (or WFS with
  `powfs.dither=1`), the uplink FSM tip/tilt error and command time history in
  radian.

- `Resltpm_1.bin`: For each LGS powfs, the common path pointing mirror
  command of the laser launch telescope.

- `Reszoompos_1.bin`: For LGS WFS when `sim.mffocus` is set, the trombone
  (zoom) position.

- `Resfocuserrs_1.bin`: For LGS WFS, the LGS focus error time history (raw
  and high pass filtered).

- `Respsdcl_hi_1.bin` and `Respsdol_hi_1.bin`: Closed loop and pseudo open
  loop PSD of the high order control computed from telemetry when
  `recon.psd` is set. `Respsdcl_lo_1.bin` and `Respsdol_lo_1.bin` are for the
  low order control in split tomography. `Resgain_1.bin` contains the
  updated control gains.

- `ResCG_1.bin`: The conjugate gradient residual time history of the
  tomography and DM fitting.

- `Resditheramp_1.bin`: For WFS with `powfs.dither` set, the dithering
  statistics (estimated signal amplitude and gain adjustment) at each update.
  Saved with `save.extra` or `save.dither`, as is `Resdithererr_wfs*_1.bin`,
  the common path dithering error history.

- `Rescn2_1.bin`: The estimated \f$r_0\f$ and \f$C_n^2\f$ profile when the
  \f$C_n^2\f$ estimation (`cn2.pair`) is enabled.

- `Restwfs_1.bin`: The radial mode error measured by the truth WFS at each
  of its outputs, when there is a truth WFS.

- `telws_1.bin`: The telescope wind shake tip/tilt time series generated from
  `sim.wspsd`.

- `dmhist_1.bin`: The histogram of the DM commands of each actuator when
  `dm.hist` is set.

\section geometry Geometry Data

When `save.setup=1`, MAOS saves the geometry data created before entering
simulation to files in folder `setup`. When `save.recon=1`, the tomography and
DM fitting matrices are also saved, which takes space and are not enabled by
`save.setup=1`. Some files are only saved with `save.setup=2` or higher, and
`save.setup>1` also enables `save.recon`, `save.mvst`, `save.ncpa`, and
`save.fdpcg`.

The following explains the common files. Not all of them exist, depending on
the options. The suffix of the file names are removed.

File name | Description
---------|--------------
`aper_locs` |The fine sampled grid on the telescope aperture for performance evaluation.
`aper_amp`  |Telescope aperture amplitude defined on aper_locs
`aper_mcc`  |Modal cross-coupling matrix for piston/tip/tilt modes defined on aper_locs
`ploc`      |The coarse sampled grid on the aperture for reconstruction
`floc`      |The grid on the aperture for DM fitting
`xloc`      |The tomography grid, one per reconstruction layer.
`aloc`      |DM actuator grid, one per DM.
`amap`      |DM actuator grid as a map, one per DM.
`W0`        |The W0 weighting defined on floc.
`W1`        |The W1 weighting defined on floc.
`saneai`    |The inverse of the NEA squared used in tomography (\f$rad^{-2}\f$)
`TT`        |The global tip/tilt modes for WFS (with save.recon)
`PTT`       |The global tip/tilt removal operator for WFS
`TTF`       |The global tip/tilt and focus modes for WFS
`PTTF`      |The global tip/tilt and focus removal operator for WFS
`GA`        |Gradient operator from aloc (with save.recon)
`GP`        |Gradient operator from ploc (with save.recon)
`HXW`       |Ray tracing from xloc to ploc along WFS directions (with save.recon)
`HA`        |Ray tracing from aloc to floc along fitting directions (with save.recon)
`L2`        |Laplacian regularization on xloc (with save.setup>1)
`FLM`, `FLU`, `FLV` |Fit left hand side operator: sparse matrix and low rank U and V matrices (with save.recon)
`FRM`, `FRU`, `FRV` |Fit right hand side operator: sparse matrix and low rank U and V matrices (with save.recon)
`fit_actcpl` |Actuator coupling used to treat inactive actuators in DM fitting
`fit_actextrap` |Actuator extrapolation operator
`tomo_RRM`, `tomo_RRU`, `tomo_RRV` |Tomography right hand side operator (with save.recon)
`tomo_RLM`, `tomo_RLU`, `tomo_RLV` |Tomography left hand side operator (with save.recon)
`MVM`       |The control matrix when `recon.mvm` is set
`ahst_MCC`  |ahst, NGS mode cross coupling matrix.
`ahst_Modes` |ahst, the NGS modes defined on DM grid
`ahst_GM`   |ahst, NGS modes to gradient operator
`ahst_Pngs_wt1` |ahst, NGS mode removal operator from DM commands (suffix is tomo.ahst_wt)
`ahst_Rngs` |ahst, NGS mode reconstructor
`focus_RFlgsg` |LGS focus reconstructor from LGS gradients
`psd_Herr`, `psd_eloc` |Ray tracing operator and grid used to compute the PSD from telemetry
`powfs0_loc`          |The grid for all subapertures (grouped by subapertures)
`powfs0_amp`          |The amplitude defined on powfs0_loc
`powfs0_pts`          |The lower left grid point of each subaperture.
`powfs0_saloc`        |The lower left corner of each subaperture
`powfs0_saa`          |The subaperture area (relative to a full subaperture)
`powfs0_GS0`          |The averaging gradient operator from powfs0_loc (with save.setup>1).
`powfs0_neasim`       |The Cholesky factor of the NEA covariance used to add noise to gradients
`powfs0_gradoff`      |The gradient reference vector
`powfs0_llt_loc`      |The aperture grid of the uplink laser launch telescope (LLT)
`powfs0_llt_amp`      |The aperture amplitude of LLT defined on powfs0_llt_loc
`powfs0_llt_imcc`     |The inverse of model cross-coupling matrix of p/t/t modes for LLT
`powfs0_llt_psf`      |The LLT uplink PSF
`powfs0_srot`         |The orientation of each subaperture wrt LLT
`powfs0_srsa`         |The distance of each subaperture from the LLT
`powfs0_sprint`       |Which subapertures to print
`powfs0_sodium`       |The sodium layer profile
`powfs0_focus`        |The additional LGS focus error time series due to sodium range variation (powfs0_llt.fnrange)
`powfs0_dtf0_nominal` |The nominal of DTF for wavelength 0 (with save.setup>1)
`powfs0_dtf0_si`      |The si of DTF for wavelength 0 (with save.setup>1)
`powfs0_etfprep0`     |The elongation transfer function for wavelength 0 (with save.setup>2)
`powfs0_sepsf`        |The subaperture short exposure PSF (tip/tilt removed, with save.setup>1)
`powfs0_i0`           |The subaperture time averaged image for matched filter
`powfs0_gx`           |The pixel by pixel x gradient of i0
`powfs0_gy`           |The pixel by pixel y gradient of i0
`powfs0_mtche`        |The matched filter gradient estimator
`powfs0_sanea`        |The subaperture noise equivalent angle (NEA) covariance in \f$rad^2\f$

\section telemetry Telemetry Data

Depending on the parameters enabled, the simulation telemetry data are saved to
files in the simulation folder. The following describes them in detail.

The options `save.ints`, `save.grad`, `save.gradnf`, `save.gradpsol`,
`save.gradgeom`, and `save.wfsopd` take either a scalar value (0 or 1, applied
to all WFS) or a vector with one value per WFS, where 1 enables saving for that
WFS.

Name convention: the last number after underscore is the seed. The following
shows results for seed 1. When there are multiple WFS or science fields, we
only show the results for the first.

Most of these files contain a cell array of length `nsim` (number of
simulation steps). Elements for steps without output (e.g., WFS that run at a
lower rate) are empty.

The data regarding the DM commands are all defined on the actuator grid `aloc`.
The gradients are in radian, with the x gradients of all subapertures followed
by the y gradients.

The second column of the tables shows which option in `dbg.conf` enables the
save of this data.

The suffix of the file names are removed.

File name          |Option to enable | Description
-------------------|-----------------|-------------
`atm_1`           |save.atm   |The atmosphere screens
`atm_winddir_1`   |save.atm   |The wind direction of each layer
`dmerr_1`         |save.dm   |The DM error signal
`dmrecon_1`       |save.dm   |The DM reconstruction (tomography and fitting) output
`Merr_lo_1`       |save.dm   |The low order mode error signal (split tomography)
`Mint_lo_1`       |save.dm   |The low order mode integrator output (split tomography, sim.fuseint=0)
`dmint_1`         |save.dm   |The DM integrator output
`dmcmd_1`         |save.dm   |The DM commands
`dmreal_1`        |save.dm   |The real, effective DM commands applied
`dmproj_1`        |save.dm   |The atmosphere projected onto the DM (sim.dmproj)
`ttmreal_1`       |save.dm   |The tip/tilt mirror commands (when sim.lpttm is set)
`wfs0_moaofit_1`  |save.dm   |The MOAO DM commands for WFS 0
`evl0_moaofit_1`  |save.dm   |The MOAO DM commands for science direction 0
`evl0_opdcl_1.fits` |save.evlopd   |The closed loop OPD for science fields (every save.evlopd steps)
`evl0_opdol_1.fits` |save.evlopd   |The open loop OPD for science fields (every save.evlopd steps)
`gcov_1_wfs0_5_10` |save.gcov   |The covariance between the PSOL gradients of WFS 0 and 5 averaged over the first 10 steps (saved every save.gcovp steps)
`ecov_1_9`        |save.ecov   |The covariance of the DM error vector accumulated until time step index 9 (saved every save.ecov steps)
`opdr_1`          |save.opdr   |The tomography output, defined on xloc
`opdx_1`          |save.opdx   |The atmosphere projected onto xloc (direct fitting)
`wfs0_gradcl_1`   |save.grad  |The WFS gradient measurement.
`wfs0_gradnf_1`   |save.gradnf  |The WFS noise free gradient.
`wfs0_gradpsol_1` |save.gradpsol  |The WFS pseudo open loop gradient
`wfs0_gradgeom_1` |save.gradgeom   |The WFS geometric gradient (in physical optics WFS simulations)
`wfs0_intsnf_1`   |save.ints   |The WFS subaperture images (noise free)
`wfs0_intsny_1`   |save.ints   |The WFS subaperture images (noisy)
`wfs0_opd_1`      |save.wfsopd   |The WFS OPD.
`wfs0_opdol_1`    |save.wfsopd   |The WFS OPD before DM correction.
`wfs0_lltopd_1`   |save.wfsopd   |The WFS laser launch telescope OPD.

## Summary of save options

All the `save.*` options are documented in `config/maos/dbg.conf`. The
following table summarizes their effects.

Option | Effect
-------|-------
`save.extra`  | Saves the additional results in subfolder `extra` (enabled by default).
`save.all`    | Sets `save.setup` and `save.run` to the same value.
`save.setup`  | Saves the setup data in subfolder `setup` (see \ref geometry). Values greater than 1 save more and also enable `save.recon`, `save.mvst`, `save.ncpa`, and `save.fdpcg`.
`save.recon`  | Saves the reconstructor matrices (large).
`save.mvst`   | Saves the MVST intermediate matrices.
`save.ncpa`   | Saves the NCPA surface OPD on the aperture and WFS.
`save.fdpcg`  | Saves the FDPCG preconditioner matrices.
`save.mvm`    | Saves the MVM control matrix (`MVM.bin`).
`save.mvmi`, `save.mvmf` | Save the intermediate results of the MVM control matrix assembly on GPU (`mvmi.bin`, `mvmf.bin`).
`save.atm`    | Saves the atmosphere. Not enabled by `save.run` or `save.all`.
`save.run`    | Saves the RTC telemetry: sets `save.dm`, `save.extra`, `save.dither`, `save.gradoff`, `save.ints`, `save.grad`, `save.gradnf`, `save.gradpsol`, and `save.gradgeom`. When greater than 1, also sets `save.opdr`, `save.opdx`, `save.evlopd`, and `save.wfsopd`, which use a lot of disk space and slow down the simulation.
`save.dm`     | Saves the DM commands and error signals.
`save.ints`, `save.grad`, `save.gradnf`, `save.gradpsol`, `save.gradgeom`, `save.wfsopd` | Save the WFS pixels, gradients, and OPDs (see the table above).
`save.evlopd` | Saves the science OPD every this many time steps.
`save.opdr`, `save.opdx` | Save the tomography output and the atmosphere projected onto `xloc`.
`save.gcov`, `save.gcovp` | Save the PSOL gradient covariance for the pairs of WFS listed in `save.gcov` every `save.gcovp` steps.
`save.ecov`   | Saves the covariance of the DM error vector every this many steps.
`save.dither` | Saves the matched filter and gradient scaling estimated from dithering (in `extra`).
`save.gradoff` | Saves the gradient reference vector updates (in `extra`).
