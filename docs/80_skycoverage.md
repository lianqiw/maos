\page skycoverage Sky Coverage

\tableofcontents

Sky coverage of an LGS MCAO system is evaluated in two steps. First, `maos`
runs a presimulation with ideal correction of the NGS controlled modes and
saves the information needed for NGS WFS modeling. Then, the separate program
`skyc` uses these data to simulate the NGS WFS and the low order control loop
for many random star fields and asterisms. The skyc configuration files are in
`config/skyc`.

\section presim MAOS presimulation

Setting `sim.skysim=1` enables the presimulation. This changes the maos
configuration as follows:
- Split tomography (`recon.split=1`) with ideal NGS mode correction
  (`tomo.ahst_idealngs=1` unless already set) is used, so that the residual
  NGS modes are open loop turbulence that the NGS loop needs to correct.
- `powfs.psfout` is set for low order powfs, so that the complex PSFs of the
  NGS WFS subapertures are saved.
- `save.extra` is set.
- An output folder (`-o`) is required.

`config/maos/skyc_10.conf` (or `skyc_20.conf`) adds `sim.skysim=1`,
`sim.end=5000`, and a grid of NGS WFS spaced by 10" (or 20") that covers the
patrol field, for a configuration with one LGS powfs (6 WFS) followed by two
NGS powfs (e.g., the TTF and TT OIWFS of NFIRAOS). The following example runs
the presimulation for NFIRAOS at three turbulence percentiles and several
zenith angles:
```
for prof in 50 25 75 ;do
	for za in 0 30 45 50 55 60;do
		maos -c nfiraos_lgs.conf skyc_10.conf atm_mk13n${prof}p.conf sim.zadeg=$za sim.seeds=[1] -o ${prof}p_za${za}
	done
done
```

The files needed by `skyc` are saved in subfolder `skysim` of the output
folder:
- `maos.conf`: The `maos.*` keys that describe the presimulation (`maos.r0z`,
  `maos.dt`, `maos.zadeg`, `maos.D`, `maos.hc`, `maos.hs`, `maos.nmod`,
  `maos.wvl`, `maos.ngsgrid`, `maos.npowfs`, `maos.msa`, `maos.nsa`,
  `maos.seeds`, `maos.nstep`, the file names below, etc.). It ends with
  `include="skyc.conf"` so that it can be used directly as the skyc
  configuration.
- `RescleNGSm_1.bin` and `RescleNGSmp_1.bin`: Hard links to the files in
  `extra`. The residual NGS modes (field averaged and for each evaluation
  direction) that the NGS loop has to correct (see \ref page40_results).
- `MCC_za0.bin` and `MCC_OA_za0.bin`: The NGS mode cross coupling matrices,
  field averaged and on axis, used to convert NGS mode to wavefront error. The
  number is the zenith angle.
- `powfs1_loc.bin`, `powfs1_amp.bin`, `powfs1_saloc.bin`: The pupil grid,
  amplitude, and subaperture locations of each NGS powfs.
- `wvfout/wvfout_seed1_sa2_x0_y0.bin`: For each NGS WFS on the grid, the time
  history of the complex PSF of each subaperture and wavelength. `sa2` is the
  subaperture order (2 for 2x2) and `x0_y0` is the location in arcsec.
- `ztiltout/ztiltout_seed1_sa2_x0_y0.bin`: The corresponding time history of
  the geometric (z-tilt) gradients.
- `gradoff/gradoff_sa2_x0_y0.bin`: The gradient offset of each NGS WFS, when
  it exists.

\section postsim Sky coverage post processing simulation

Change directory to the `skysim` folder of the presimulation and run `skyc`:
```
cd ${prof}p_za${za}/skysim
skyc -d -o APD_JHKs_typeImr skyc.multirate=1 skyc.seed=1
```
Options:
- `-o DIR`: Output folder (relative to the current folder). A time stamped
  folder is used if omitted.
- `-c FILE.conf`: Use `FILE.conf` as the baseline configuration. By default
  `maos.conf` in the current folder is used.
- `-d`: Run in the background (detached) through the scheduler.
- `-n N`: Number of threads.
- `-O`: Rerun seeds that are already done.
- `key=value` or `FILE.conf`: Override parameters.

`skyc` reads `maos.conf`, which includes `config/skyc/skyc.conf` for the
default `skyc.*` parameters. Zenith angle dependent parameters (detector
pixels per subaperture and implementation error) are in
`config/skyc/skyc_za*.conf`, which can be added to the command line, e.g.,
`skyc_za30.conf`. The main parameters are:

Parameter | Description
----------|------------
`skyc.npowfs`, `skyc.nwfsmax` | Number of NGS powfs (must equal `maos.npowfs`) and maximum number of WFS of each powfs in an asterism.
`skyc.pixtheta`, `skyc.pixpsa`, `skyc.pixblur`, `skyc.pixguard`| Detector pixel size (arcsec), number of pixels per subaperture, pixel blur, and guard window of each powfs.
`skyc.rne`, `skyc.excess` | Read out noise (a negative value selects a frame rate dependent model) and excess photon noise factor.
`skyc.zb.*` | Zero magnitude flux, sky background, detector QE, throughput, and excess background for each wavelength.
`skyc.patfov`, `skyc.keepout`, `skyc.minrad` | Diameter of the patrol field, probe keep out distance, and minimum star radius (arcsec).
`skyc.nsky`, `skyc.start` | Number of star fields and index of the first one.
`skyc.lat`, `skyc.lon`, `skyc.catscl` | Galactic latitude and longitude of the Besancon star catalog (`config/skyc/besancon`) used to generate random star fields, and the star count scaling.
`skyc.stars` | Load star fields from this file instead of generating them. It is a cell array with one element per star field, each containing a (2+nwvl) by nstar array of star location (radian) and magnitudes.
`skyc.maxstar`, `skyc.maxaster`, `skyc.maxdtrat` | Maximum number of (brightest) stars per field, asterisms per field, and sampling rates per asterism to simulate.
`skyc.dtrats`, `skyc.multirate`, `skyc.snrmin_fast`, `skyc.snrmin_slow` | Candidate sampling periods (in units of `maos.dt`, in descending order). With `skyc.multirate=1`, each OIWFS can run at a different rate, selected using the minimum SNR thresholds.
`skyc.servo`, `skyc.gsplit` | Servo type (1: integrator, 2: type II, -1: LQG) and separate gain for each mode.
`skyc.psd_ws`, `skyc.addws` | Telescope wind shake PSD (skip it if the maos presimulation already includes `sim.wspsd`) and how it is added.
`skyc.evlstart`, `skyc.phystart` | Time step to start performance evaluation and physical optics simulation.
`skyc.phytype` | Pixel processing: 1: matched filter, 2: CoG, 3: correlation.
`skyc.imperrnm`, `skyc.imperrnmb` | Additional implementation error on axis and at the edge of the patrol field (nm).
`skyc.estimate` | Only estimate the performance with servo analysis, without time domain simulation.
`skyc.seed` | Seed for the star fields and random numbers.
`skyc.save` | Save the selected asterism of each star field for a maos simulation.

For each star field, `skyc` forms the asterisms from the available stars,
estimates their performance with servo analysis, and runs time domain physical
optics simulations of the NGS WFS (using the saved PSFs) for the most
promising asterisms and sampling rates. The best asterism determines the
result of the star field.

Before the first run, `skyc` processes the saved PSFs and stores the results in
the `skysim` folder so that later runs can reuse them: subfolders `pistat`
(time averaged subaperture PSFs), `neaspec`, and `phygrad`, and the WFS
nonlinearity `powfs0_nonlin.bin`, `powfs1_nonlin.bin`.

\section skyres Sky coverage results

The following file names are given for maos seed 1 (`maos.seeds`) and skyc
seed 1 (`skyc.seed`), in the output folder of `skyc`. The results are mean
square wavefront error in \f$m^2\f$. Use `drawres` to plot the sky coverage,
i.e., the cumulative distribution of the WFE of the star fields (see \ref
page40_results).

- `Res1_1.bin`: An array with one column per star field and the following
  rows:
	- Total NGS mode WFE (including residual wind shake when it is not part of
	  the simulated time series)
	- NGS mode WFE due to turbulence
	- Tip/tilt WFE due to turbulence
	- Residual wind shake tip/tilt WFE
	- Estimated WFE from the servo analysis
	- On axis NGS mode WFE due to turbulence
	- On axis tip/tilt WFE due to turbulence
	- Sampling frequency (Hz) of the selected asterism
	- Column index of the selected asterism in `Res1_1_aster`
	- Index of the selected asterism within the star field

  Star fields without a usable asterism have the open loop WFE.
- `Res1_1_aster.bin`: All asterisms that have physical optics simulation
  results. Each column is for an asterism. The rows are:
	- WFE in physical optics simulations
	- WFE from servo analysis
	- Star field index
	- Number of WFS
	- For each WFS:
		- sampling frequency (Hz) (0 if not used)
		- signal to noise ratio at this sampling frequency
		- x coordinate (radian)
		- y coordinate (radian)
		- magnitude in each wavelength (`maos.wvl`)
- `Res1_1_stars.bin`: The star fields used, in the same format as `skyc.stars`.
- `Res1_1_mres.bin`: The residual NGS mode time history of the selected
  asterism for each star field.
- `Res1_1_gain.bin`: The servo gains of the selected asterism for each star
  field (when `skyc.servo>0` and `skyc.multirate=0`).
- `Res1_1_minput.bin`: The input NGS modes (with wind shake added when
  `skyc.addws` is set).
- `psds_m2.bin`: The PSD of the input NGS modes.
- `Res1_1.done`: Created when the seed is done. It is skipped in later runs
  unless `-O` is given.
- `skyc_[HOST]_[PID].conf` and `run_[HOST]_[PID].log`: The effective
  configuration and the log. Renamed to `skyc_done.conf` and `run_done.log`
  when finished successfully.
- `Res1_1_maos/sky0/`: When `skyc.save=1`, the PSFs, NEA, residual modes, and
  a `base.conf` of the selected asterism of each star field, for running maos
  with that asterism.
