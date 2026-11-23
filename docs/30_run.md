\page page30_run Run simulations
\tableofcontents
\section sect-run Usage

We assume that `maos` is already in the PATH so we can type `maos` to launch it.

Use the following command to get the help message:
```
maos -h  #will print the help message and exit.
```

The general form of the command line is

```
maos [OPTION...] [-c base.conf] [override.conf ...] [key=value ...] [-o DIR]
```

Options, `.conf` files, and `key=value` entries can appear in any order. A
value that contains spaces, such as an array `[1 2 3]`, does not need to be
quoted: maos joins all the arguments and keeps text inside `[]` or quotes
together. Quote it anyway if your shell would otherwise expand it.

Examples:
```
maos -o result    # Run the default configuration and save results to folder result.

maos -c scao_ngs.conf sim.seeds=[1 10] -d -o scao_ngs override.conf
                  # Run a single conjugate natural guide star case, with seed 1 and 10,
                  # detach from the terminal, output results to folder scao_ngs,
                  # and read overriding parameters from override.conf.
```

\subsection sect-run-options Command line options

| Option | Description |
|--------|-------------|
| `-h`, `--help` | Print the help message and exit. |
| `-c FILE.conf`, `--conf=FILE.conf` | Use `FILE.conf` as the baseline (main) configuration instead of `default.conf`. |
| `-o DIR`, `--output=DIR` | Save results to `DIR` (created if needed). maos changes into `DIR` before running. Without `-o`, nothing is saved, and options that require saving (`save.setup`, `save.all`, `sim.skysim`, `evl.psfmean`, `evl.psfhist`) cause an error. |
| `-n N`, `--nthread=N` | Use `N` threads. The default is the number of physical CPU cores (hyperthreads are not counted), or `OMP_NUM_THREADS` if it is set and smaller than the number of online processors. The best value depends on the configuration: more threads than needed add synchronization overhead. |
| `-g i`, `--gpu=i` | Use GPU `i` (in `nvidia-smi` order). Repeat to use several GPUs, e.g. `-g0 -g2`. `-g-1` disables the GPU. Default: automatic. |
| `-G N`, `--ngpu=N` | Use `N` GPUs, chosen automatically. When both `-g` and `-G` are given, `-g` takes priority. |
| `-d`, `--detach` | Detach from the terminal and run in the background. Requires `-o`. |
| `-f`, `--force` | Start immediately without waiting for the scheduler to allocate CPUs. Jobs running in the foreground always start immediately. |
| `-O`, `--override` | Rerun seeds that are already finished (have a `Res_<seed>.done` file). See \ref sect-run-seeds. |
| `-r host`, `--run=host` | Run the job on another host. Requires the scheduler (`scheduler`) on that host. |
| `-p DIR`, `--path=DIR` | Add `DIR` to the search path for `.conf` and data files. |

The help message printed by `maos -h` also lists a few environment variables.
See \ref sect-run-env for the complete list.

\subsection sect-run-sched Scheduler and detached runs

When maos runs in the foreground (no `-d` and no `-r`), it starts right away and
reports its progress to the scheduler, which is launched automatically if it is
not already running. The `monitor` program (GTK or the browser interface)
connects to the scheduler to show the progress of all jobs and to plot results.

With `-d` or `-r host`, maos hands the job to the scheduler, which queues it
and launches it when enough CPUs (and GPUs) are free. If the scheduler cannot be
reached, maos warns and runs the job by itself in the background. Without
`-f`, such a job still waits until enough CPUs are idle. When maos is built with
`--disable-scheduler` (or `--enable-release`), detached jobs start immediately.

\subsection sect-run-gpu GPU selection

When maos is compiled with CUDA, it uses the GPUs by default. Without `-g` or
`-G`, maos uses all usable GPUs that have enough free memory, the least busy
ones first. Use `-g i` (repeatable) to select GPUs by their `nvidia-smi` index,
`-G N` to use `N` automatically selected GPUs, or `-g-1` to run on CPU only.
Setting the environment variable `MAOS_GPU_DISABLE_<i>=1` excludes GPU `i` from
use.

Which parts of the simulation run on the GPU is controlled by the `gpu.*` keys
in \ref sim, all enabled (`1`) by default: `gpu.wfs` (wavefront sensors),
`gpu.evl` (performance evaluation), `gpu.tomo` (tomography), `gpu.fit` (DM
fitting), `gpu.lsr` (least squares reconstruction), `gpu.psf` (PSF
accumulation), and `gpu.moao`. maos turns off some of them automatically when
the configuration is not supported on GPU, for example `gpu.evl` when
`evl.rmax>1`, or `gpu.tomo` when `tomo.cxxalg` is not 0.

\subsection sect-run-seeds Seeds, restarting, and multiple runs

`sim.seeds` is a list of random seeds (default `[1]`). Each seed is a separate
realization of the turbulence and noise, and the seeds are simulated one after
another in the same process:

```
maos sim.seeds=[1 2 3 4] -o mydir
```

When a seed finishes, maos creates `Res_<seed>.done` in the output folder.
When maos starts, it skips every seed that already has a `.done` file and every
seed that another maos process is running in the same folder (detected with a
lock file). Therefore:

- To resume after an interruption, rerun the same command. Finished seeds are
  skipped; an unfinished seed is restarted from the beginning.
- Several maos processes can share an output folder and a seed list; each seed
  is run only once.
- Use `-O` to rerun finished seeds anyway. If there is no seed left to run,
  maos prints `There are no seed to run. Use -O to override.` and exits.

`sim.start` and `sim.end` (default `0` and `500` in \ref sim) set the first and
last (exclusive) time step. `sim.dt` (default `1/600` s) is the loop time step.

\subsection sect-run-output Output folder

With `-o DIR`, maos writes the following in `DIR`:

- `Res_<seed>.bin` and `Resp_<seed>.bin`: the wavefront error time history and
  other results. See \ref page40_results for the content and how to read them.
- `Res_<seed>.done`: marks a finished seed.
- `extra/`: additional results when `save.extra=1` (the default).
- `maos_<host>_<pid>.conf`: every configuration key with the value used. It
  starts with `include=` the main `.conf` file followed by all the keys; keys
  that keep the value from the main configuration are commented out. Rerun the
  same configuration with
  `maos -c maos_done.conf -o newdir`.
- `run_<host>_<pid>.log`: the log of the run.
- When the run finishes normally, the two files are renamed to `maos_done.conf`
  and `run_done.log` (an existing pair from a previous run is renamed to
  `maos_done.N.conf` and `run_done.N.log`). When the run is killed or fails, the
  log is renamed to `run_<host>_<pid>.err`. The symbolic links `maos_recent.conf`
  and `run_recent.log` point to the most recent files.
- `setup/`: data structures prepared during setup when `save.setup=1`.
- `skysim/`: inputs for the sky coverage post processor when `sim.skysim=1`.
- PSF `.fits` files (`evlpsfcl_<seed>_x<thetax>_y<thetay>.fits`, etc.) when
  `evl.psfmean` or `evl.psfhist` is set.

Other `save.*` keys in \ref dbg save telemetry such as gradients, DM commands,
and OPDs.

\subsection sect-run-env Environment variables

The following environment variables are read by maos. Variables starting with
`MAOS_` can also be given as `MAOS_X=value` in a `.conf` file or on the command
line, which sets the environment variable before it is used. This does not
work for `MAOS_CONFIG_PATH`, `MAOS_PERMISSIVE`, and the low level flags at the
end of this section, which are read before the configuration is loaded; set
them in the shell.

| Variable | Default | Description |
|----------|---------|-------------|
| `MAOS_CONFIG_PATH` | (unset) | Folder that contains the configuration files: `$MAOS_CONFIG_PATH/config/maos` or `$MAOS_CONFIG_PATH/maos`. See \ref sect-config-path. |
| `MAOS_LOG_LEVEL` | 0 | Amount of output, from -5 to 5. -3: errors and warnings only; -2: also the most important information; -1: also useful information; 0: all information; 1, 2, 3: increasing amount of debugging output. |
| `MAOS_PERMISSIVE` | 0 | Set to 1 to turn the errors for unrecognized keys (or keys used more than once) into warnings. |
| `MAOS_PARALLEL` | 1 | 1: run WFS, reconstruction, and performance evaluation of each time step in parallel. 0: run them one after another. 2: event driven parallel mode (requires OpenMP; falls back to 1 otherwise). |
| `MAOS_NO_WFS` | 0 | Set to 1 to skip all WFS calls (for debugging and timing). |
| `MAOS_NO_EVL` | 0 | Set to 1 to skip performance evaluation calls. |
| `MAOS_NO_RECON` | 0 | Set to 1 to skip reconstruction calls. |
| `MAOS_KEEP_MEM` | 0 | Set to 1 to keep some temporary memory between time steps. |
| `MAOS_TOMOSCALE` | 1e-12 | Scaling applied to the turbulence regularization terms in tomography to keep numbers in range, needed for single precision. |
| `MAOS_MEM_DEBUG` | 0 (1 in debug builds) | Set to 1 to enable malloc/free accounting, which reports memory leaks at exit. |
| `MAOS_MEM_VERBOSE` | 0 | Set to 1 or 2 to print every malloc/free call. |
| `MAOS_MEM_FUNTRACE` | 0 | Set to 1 to record the allocating function name instead of a backtrace in memory accounting. |
| `MAOS_GPU_DISABLE_<i>` | (unset) | Set to a nonzero value to exclude GPU `<i>` (nvidia-smi index). |
| `MAOS_MISREG_SEQ` | 1 | Random seed used to randomize WFS misregistration set by `powfs.misregx`, `powfs.misregy`, and `powfs.misregc`. |
| `MAOS_DRAW_NOREUSE` | 0 | Set to 1 to start a new `drawdaemon` window instead of reusing an existing one. |
| `OMP_NUM_THREADS` | (unset) | When smaller than the number of online processors, sets the default number of threads (see `-n`). |
| `CUDA_LAUNCH_BLOCKING` | (unset) | When nonzero, maos uses only the default CUDA stream. |
| `MAOS_DIRECT_LAUNCH` | (internal) | Set by the scheduler when it launches a job. Do not set it. |

The following variables are for testing and developing maos:
`MAOS_PYWFS_DEBUG`, `MAOS_PYWFS_TT_DUAL`, `MAOS_PYWFS_PSIZE` (default 3),
`MAOS_PYWFS_FULL` (pyramid WFS), `MAOS_LGS_SPH_PSD`, `MAOS_SODIUM_FIT`,
`MAOS_DITHER_ND`, `MAOS_DITHER_MD2`, `MAOS_DITHER_MD3`, and the low level flags
`MAOS_PROP_GRID_MAP_OPTIM`, `MAOS_LOC_MAP_EXTEND`, `MAOS_SVD_SDD_THRES`, and
`MAOS_CN2EST_NO_NEGATIVE`.

\section sect-config Configuration Files

All user modifiable options are configured through config files with suffix
\c .conf or `key=value` entries on the command line. The supplied default setup files
are stored in sub-folder \c config/maos. The config files are organized in a
modular way. One \c .conf file can include another \c .conf file with \c
include=filename.conf. The \c .conf files in \c config/maos
are self-documented. \ref page34_conf lists all of them.

All the options are in the form of \c key=value. The usual convention is:
if the entry in `parms` is called `parms->wfs.thetax`, the entry in the .conf file will
be in the form `wfs.thetax=value`.

\subsection sect-config-syntax Syntax

- Anything after `#` is a comment.
- A line ending with `\` continues on the next line.
- `include=file.conf` (or a line with only `file.conf`) reads another file at
  that place. `path=DIR` adds `DIR` (relative to the folder of the current file)
  to the search path.
- `key=value` sets a value. The same key can be set many times. Within the same
  level (see below), the value read last wins.
- `key+=[...]` appends to an existing array and `key-=[...]` removes the given
  values from the end of the array. Both values must be enclosed in `[]`.
- `key=ref:other.key` uses the value of another key, e.g.,
  `sim.aplo=ref:sim.aphi`.
- Numbers can be written as simple expressions without spaces, such as `1/600`,
  `0.5*2`, or `30+5`. `inf` and `nan` are accepted.
- Numerical arrays are enclosed in `[]` with entries separated by spaces or
  commas: `[1 2 3]` or `[1, 2, 3]`.
  - `n@m` repeats the number `n` `m` times: `[0.5@3]` is `[0.5 0.5 0.5]`.
  - A factor can be applied to the whole array: `[1 2 3]*2`, `2*[1 2 3]`,
    `[1 2 3]/2`, `2/[1 2 3]`, or `[1 2 3]+1`.
  - `;` starts a new column: `[1 2 3; 4 5 6]` has 3 rows and 2 columns.
    A trailing `'` transposes: `[1 2 3; 4 5 6]'` has 2 rows and 3 columns.
  - A key that takes a numerical array can also take a file name (`.bin` or
    `.fits`, optionally gzipped) that contains the array.
- Strings may be enclosed in single or double quotes. Quotes are needed only
  when the string contains spaces or commas. String arrays use `[]` with
  entries separated by commas: `["a.bin", "b.bin"]` or `[a.bin,b.bin]`.
  `[]` is an empty array; `[""]` is an array with one empty string.
- Many per-type arrays, such as most `powfs.*` keys, accept either one value per
  element or a single value that is used for all elements.

Old key names that have been renamed are listed in \ref changes and are
converted automatically, with a warning if the key was renamed with `=>`.

\subsection sect-config-override Overriding and checking keys

Configuration is read in three levels, from low to high priority:

1. The main `.conf` file (`-c FILE.conf`, or `default.conf`) and all files it includes.
2. `.conf` files given on the command line, and the files they include.
3. `key=value` entries given on the command line.

A value from a higher level always overrides a lower one, no matter in which
order they appear on the command line. Within the same level, the value read
last wins.

maos checks that every supplied key is used:

- A key in level 2 or 3 that maos does not recognize is an error, so a
  misspelled key cannot silently be ignored. In level 1, and for `dbg.*`
  keys, it is only a warning. Set `MAOS_PERMISSIVE=1` to make all of these
  warnings.
- A key that maos reads but is not defined anywhere gives a warning and is
  assumed to be 0 or empty. The supplied `.conf` files define every key, so this
  only happens with incomplete custom main `.conf` files.
- Invalid values or array lengths are errors.

The full list of keys and values actually used is saved to `maos_done.conf` in
the output folder (see \ref sect-run-output).

\subsection sect-config-path Search path

The configuration folder is the first of the following that exists:

1. `config/maos` next to the maos executable,
2. `$MAOS_CONFIG_PATH/config/maos` or `$MAOS_CONFIG_PATH/maos`,
3. `config/maos` in the source tree that maos was compiled from,
4. `~/.aos/config-<version>/maos`.

`.conf` files and data files named in them (e.g. `aper.fnamp`) are searched in
the folder where maos is launched (before changing into the `-o` folder),
then in the configuration folder and its `bin`, `atm`, and `examples`
sub-folders, then in folders added by `-p` or `path=`. An absolute path is
used as is.

Without `-c`, the main file is `~/.aos/default.conf` if it exists, otherwise
`default.conf` found through the search path (so a `default.conf` in the
launch folder takes precedence).

\subsection sect-config-org Organization

The recommended way to organize config files for a specific problem is to
keep a baseline \c .conf file that specifies all the necessary \c key=value
entries, usually by including other \c .conf files. The default entry is
`default.conf` unless `-c file.conf` is supplied on the command line.

Put the values you want to change on the command line (quote them if your shell
requires it) or in one or more overriding `.conf` files in the simulation
folder.

The \c .conf files in the \c config/maos folder are only modified during
code development and should not be modified by the end user. Use overriding
\c .conf files in simulation folders or embed configurations on the command
line to modify the values.

In the \c config/maos folder, there are a few
baseline \c .conf files that are complete for each configuration:

- \ref mcao_lgs for LGS MCAO simulations (the default through \ref default),
- \ref mcao_ngs for NGS MCAO simulations,
- \ref scao_lgs for single conjugate LGS simulation,
- \ref scao_ngs for single conjugate NGS simulation,
- \ref scao_pywfs for single conjugate NGS simulation with a pyramid WFS,
- \ref ltao for LGS LTAO simulations,
- \ref glao for GLAO simulations.

Each of these files includes a few other \c .conf files by specifying \c
include=`filename.conf`. For example, all of them
include the following common `.conf` files that are independent of the
simulation mode or geometry.

- \ref sim to specify general simulation parameters such as the telescope
  aperture (it includes `keck.conf` from `examples/`), loop gain, number of time steps, etc.,
- \ref recon to specify reconstruction parameters in closed loop,
- \ref dbg to specify debugging, plotting, saving, and loading parameters.

The following list the exclusive `.conf` files used to specify a particular
set of configurations. Each group should only appear once.

- `atm_*.conf` (in `atm/`): a preconfigured turbulence profile.
- `dm_*.conf`: the deformable mirror configuration.
- `fov_*.conf`: the science FoV (`evl_*.conf`) and DM fitting FoV (`fit_*.conf`).
- `wfs_*.conf`: the set of WFS types. Each is composed of `powfs_*.conf` files.
  \ref powfs_none (included first) clears the WFS list, and each `powfs_*.conf`
  appends one WFS type with `+=`.
- `examples/*.conf`: telescope apertures (`keck.conf`, `tmt.conf`, `gemini.conf`, `elt.conf`) and
  complete configurations of existing or planned AO systems.

For example, the default configuration \ref mcao_lgs includes \ref mcao_base and
\ref wfs_lgs_ttf_tt. \ref mcao_base includes \ref sim, \ref dm_dual for a dual DM
system, \ref atm_mk13n50p for the Maunakea 13N median seeing turbulence profile,
\ref fov_sq34 for the science FoV (\ref evl_sq34, Simpson weighting on 3x3
points in the square 34x34" FoV) and the DM fitting FoV (\ref fit_sq34, the same
points), \ref recon, and \ref dbg. \ref wfs_lgs_ttf_tt defines six LGS WFS, one
tip/tilt/focus NGS WFS, and two tip/tilt NGS WFS.

For an MCAO configuration, the DM fitting FoV can be larger than the science FoV
in order to sharpen diffraction limited TTF NGS. For an LTAO configuration, use
\ref fit_oa to set only on axis DM fitting FoV. A GLAO system can be simulated
by having a single DM but with a wide DM fitting FoV, as an alternative to
setting `recon.glao=1` which averages WFS measurements before reconstruction.

See \ref page33_example for more detailed explanations.

\subsection sect-config-params Guide to the main parameters

The table below lists the most frequently changed keys in each group. The
defaults are in the listed file; the comment next to each key in that file
describes it.

| Group | File | Main keys |
|-------|------|-----------|
| `aper` | \ref sim, `examples/keck.conf` | `aper.d` (outer and inner diameter, `[11.0 2.4]` from `keck.conf`), `aper.fnamp` (pupil amplitude map), `aper.rotdeg` |
| `sim` | \ref sim | `sim.seeds`, `sim.end` (500), `sim.dt` (1/600), `sim.zadeg` (30), `sim.closeloop` (1), `sim.ephi` (0.5, high order integrator gain), `sim.alhi` (extra latency), `sim.eplo`, `sim.fcfocus`, `sim.idealtomo`, `sim.evlol`, `sim.skysim` |
| `atm` | \ref sim, `atm/atm_*.conf` | `atm.r0z` (0.186 m in \ref atm_mk13n50p), `atm.L0` (30 m), `atm.ht`, `atm.wt`, `atm.ws`, `atm.wddeg`, `atm.frozenflow`, `atm.size`, `atm.dx` (1/64 m) |
| `atmr` | \ref sim | `atmr.ht`, `atmr.wt` (reconstructed layers, default follow `atm`), `atmr.os` |
| `dm` | \ref sim, `dm_*.conf` | `dm.ht`, `dm.dx` (0.5 m actuator pitch), `dm.offset`, `dm.iac`, `dm.stroke`, `dm.actstuck`, `dm.actfloat` |
| `powfs` | `powfs_*.conf`, \ref powfs_common | per WFS type: `powfs.type` (0: SHWFS, 1: pyramid), `powfs.dsa`, `powfs.nwfs`, `powfs.hs`, `powfs.wvl`, `powfs.dtrat`, `powfs.siglev`, `powfs.bkgrnd`, `powfs.rne`, `powfs.noisy`, `powfs.phystep`, `powfs.phytype_recon`, `powfs.pixtheta`, `powfs.pixpsa`, `powfs.nearecon`, `powfs.neasim`, `powfs.fnllt` |
| `wfs` | `powfs_*.conf` | per WFS: `wfs.thetax`, `wfs.thetay` (arcsec, concatenated for all WFS types), `wfs.siglev` |
| `llt` | \ref llt_CL, \ref llt_SL | laser launch telescope and sodium layer: `llt.d`, `llt.ox`, `llt.oy`, `llt.fnprof`, `llt.fnrange` |
| `evl` | \ref sim, `evl_*.conf` | `evl.thetax`, `evl.thetay`, `evl.wt`, `evl.fov`, `evl.wvl`, `evl.psfmean`, `evl.psfisim`, `evl.psfsize`, `evl.psfgridsize` |
| `fit` | \ref recon, `fit_*.conf` | `fit.thetax`, `fit.thetay`, `fit.wt`, `fit.fov`, `fit.alg`, `fit.actslave` |
| `recon`, `tomo`, `lsr` | \ref recon | `recon.alg` (0: MVR, 1: LSR), `recon.split` (1), `recon.glao`, `recon.mvm`, `recon.modal`, `recon.psd`, `tomo.alg`, `tomo.precond`, `tomo.maxit`, `lsr.alg`, `lsr.actslave` |
| `ncpa` | \ref sim | `ncpa.surf`, `ncpa.calib` (1), `ncpa.thetax`, `ncpa.thetay` |
| `moao` | \ref sim | `moao.dx`, `moao.gdm`, `evl.moao`, `powfs.moao` |
| `gpu` | \ref sim | `gpu.wfs`, `gpu.evl`, `gpu.tomo`, `gpu.fit`, `gpu.lsr`, `gpu.psf`, `gpu.moao` (see \ref sect-run-gpu) |
| `save`, `plot`, `load`, `dbg` | \ref dbg | `save.setup`, `save.run`, `save.all`, `save.extra`, `plot.all`, `plot.setup`, `plot.run`, `load.*` |

Values that have physical units use SI units unless noted otherwise. Angles
such as `wfs.thetax` and `evl.thetax` are in arcsec. When `evl.fov` (or
`fit.fov`) is not 0 or 1, `evl.thetax` and `evl.thetay` are given in units of
`evl.fov` (the full width in arcsec) and scaled by it. Similarly, when
`powfs.astscale` is not 0 or 1, `wfs.thetax` and `wfs.thetay` of that powfs
are multiplied by it, so the asterism is given with a radius of 0.5 and
`powfs.astscale` is its diameter in arcsec (maos checks the radius unless all
WFS of the powfs are on axis).

\section sect-exe Sample Runs

To run predefined AO modes:

```
maos -o mcao                # default is LGS MCAO with 2 DMs. Save results to folder mcao.
maos -c mcao_ngs.conf       # NGS MCAO
maos -c scao_ngs.conf       # NGS SCAO
maos -c scao_pywfs.conf     # NGS SCAO using a pyramid WFS
maos -c scao_lgs.conf       # LGS SCAO
maos -c ltao.conf           # LGS LTAO
maos -c glao.conf           # GLAO
maos -c mcao_lgs.conf dm_single.conf evl.moao=0 moao.dx=[1/2] # LGS LTAO with MOAO for each science field.
maos -c mcao_lgs.conf dm_triple.conf # LGS MCAO with 3 DMs
```

To run simulation for sample AO systems (in `config/maos/examples`):

```
maos -c keck_ngs.conf          # Keck NGS AO
maos -c keck_kapa.conf         # Keck KAPA LGS LTAO
maos -c nfiraos_lgs.conf       # TMT NFIRAOS LGS MCAO
maos -c nfiraos_lgs_pyttf.conf # TMT NFIRAOS LGS MCAO with a pyramid WFS as the TTF sensor
maos -c nfiraos_ngs.conf       # TMT NFIRAOS NGS AO with a pyramid WFS
```

The default case starts with `default.conf`, which includes \ref mcao_lgs. Put
your own `~/.aos/default.conf` to change the default (see \ref sect-config-path).

To customize the number of components

```
maos dm_single.conf wfs_lgs_only.conf recon.glao=1 # LGS GLAO
maos dm_triple.conf dm.dx=0.3  # triple DM with 0.3 m actuator pitch
maos wfs_lgs_ttf_tt_twfs.conf powfs.dtrat=[1 1 1 100] # with LGS, TTF, TT, and truth WFS. The truth WFS must have dtrat>1.
maos evl_x.conf                # evaluate performance along the x axis.
```

Change aperture

```
maos aper.d=[10]    # 10 meter circular
maos aper.d=[10, 1] # 10 meter annular diameter with 1 meter inner obscuration
maos aper.d=[10, 1] aper.fnamp=aper.bin # with supplied amplitude map.
maos tmt.conf       # TMT aperture. Also keck.conf, gemini.conf, elt.conf.
```

When `aper.fnamp` is given, its diameter must agree with `aper.d` (within 1 m),
otherwise the amplitude map is not used. Set `aper.d=[0]` to take the
diameter from the amplitude map.

Change turbulence

```
maos atm_mk13n25p.conf # use another predefined profile. The default is atm_mk13n50p.conf
maos atm.r0z=0.1     # change r0 at zenith to 0.1m
maos sim.zadeg=45    # change zenith angle to 45 degrees (default is 30)
maos atm_single.conf # use single ground layer turbulence
```

Configuring WFS. Adjust the number of elements depending on how many powfs are
in use. A powfs is a type of WFS. Each type can have multiple WFS. The default
configuration has three: LGS, TTF NGS, and TT NGS.

```
maos powfs.noisy=[0]     # use noise free WFS for all wfs types
maos powfs.noisy=[1 0 0] # use noise free WFS for wfs type 1 and 2.
maos powfs.phystep=[-1]  # use geometric wfs instead of physical optics wfs for all
```

Adjust the controller

```
maos sim.ephi=0.3 # Change gain of the main integrator to 0.3
```

\section sect-ebb Error Budget Breakdown

There are options to allow probing for error budget items, such as DM fitting error, tomography error, etc.

```
maos -c scao_ngs.conf sim.closeloop=0 atm.frozenflow=1 sim.idealtomo=1 # step1: single DM fitting on axis
maos -c mcao_ngs.conf sim.closeloop=0 atm.frozenflow=1 sim.idealtomo=1 # step1b: Multiconjugate DM fitting over a field
maos -c scao_ngs.conf sim.closeloop=0 atm.frozenflow=1 powfs.phystep=[-1] powfs.noisy=[0] # step2: noise free classic AO simulation with geometric SHWFS.
maos -c mcao_ngs.conf sim.closeloop=0 atm.frozenflow=1 powfs.phystep=[-1] powfs.noisy=[0] # step2b: noise free MCAO with geometric SHWFS.
maos -c mcao_ngs.conf powfs.phystep=[-1] powfs.noisy=[0] # step3: noise free MCAO with geometric SHWFS in closed loop simulation.
maos -c mcao_ngs.conf powfs.noisy=[0] # step4: noise free MCAO with physical optics WFS.
maos -c mcao_lgs.conf # step5: noisy MCAO with physical optics LGS WFS.
```

Note that the option `sim.closeloop=0 atm.frozenflow=1` uses open loop
correction (no servo lag) with frozen flow turbulence. The option
`sim.idealtomo=1` enables DM fitting directly from turbulence rather than from
tomography output. For open loop simulations, truth WFS should not be included
or should be disabled with `sim.eptwfs=0`, and non-iterative solvers
(`tomo.alg=0`, `fit.alg=0`) usually work better than CG.

The following error budget items are computed:

- `DM fitting error`: step1 gives DM fitting error which depends on the actuator spacing.
- `DM projection error`: the quadrature difference between step1 and step1b gives DM projection error which depends on the field of view and number of DMs.
- `WFS aliasing error`: the quadrature difference between step1 and step2 gives WFS aliasing error.
- `Tomography error`: the quadrature difference between the RSS of the three terms above and step2b is the tomography error which depends on the asterism geometry.
- `Servo lag error`:  the quadrature difference between step2b and step3 gives servo lag error.
- `WFS nonlinearity error`: the quadrature difference between step3 and step4 gives WFS non-linearity error.
- `WFS noise effect`: the quadrature difference between step4 and step5 gives WFS noise error.

Other keys in \ref sim that help isolate error terms:

- `sim.evlol=1` evaluates the uncorrected (open loop) turbulence only and skips
  WFS and reconstruction.
- `sim.idealwfs=1` uses only the turbulence in the DM correctable space for the
  WFS (removes aliasing). It conflicts with `sim.wfsalias`.
- `sim.wfsalias=1` studies the WFS aliasing by projecting the turbulence onto
  the space orthogonal to the DM vector space.
- `sim.idealevl=1` evaluates the performance only in the DM correctable space,
  ignoring DM fitting and projection error.
- `sim.noatm=1` disables the atmosphere, e.g., to evaluate telescope surface
  errors only.

\section advanced Advanced configuration

\subsection sect-surface Specifying Surface OPDs

We can optionally set up one or more static surfaces that cover science fields
and/or wavefront sensors. Each surface file must contain a 2-d array of the
OPD with a header specifying the following keys. If `dx` is missing, it is
computed from `D` (the diameter) or assumed to be 1/64 m. If `ox` or `oy` is
missing, the OPD is centered on the pupil.

```
dx  the sampling in x direction (first dimension).
dy  the sampling in y direction (second dimension). Default is dx.
ox  the origin in x direction.
oy  the origin in y direction.
h   the height conjugation of this surface (default 0)
vx  the frozen flow wind speed along x of this surface (0 for static)
vy  the frozen flow wind speed along y of this surface (0 for static)
```

Optionally, the header can also include the following keys to indicate its coverage

```
SURFNAME=name       #name of the surface. M1 or M2 for the primary or secondary mirror.
                    #M1 and M2 are rotated by aper.rotdeg. M1 is also shifted by aper.misreg.
SURFEVL=[1 1 1 ...] #length: nevl. 1: enabled for this science evaluation direction (assume all 1 if omitted)
SURFWFS=[1 1 1 ...] #length: 1, nwfs, or npowfs. 1: enabled for this WFS (assume all 1 if omitted)
```

Use the \c write mex routine, \c writebin.m, or \c aolib.writebin to write the bin file:
```
   write(OPD, header, 'opd1.bin')
```
Or simply use fits format. Put the list of surface file names in key \c ncpa.surf (\c surf is an accepted alias).
```
   maos ncpa.surf=['opd1.bin','opd2.bin', 'opd3.fits']
```

It is also possible to generate a surface on the fly with an inline
description. The keys are `r0` (in m) or `rms` (in nm), `slope` (PSD power law,
default -11/3), `L0` (default `nx*dx`), `nx` (default 2048), `ny`, `dx`
(default 1/64), `seed`, `ht`, `vx`, `vy`, plus the coverage keys above. `mode=n`
generates the Zernike mode `n` with the given `rms` instead.
```
maos ncpa.surf=["'r0=0.1;slope=-4;SURFEVL=1;', 'r0=0.2; slope=-4;SURFWFS=1 1 1 1 1 1 0 0 0'"]
maos ncpa.surf=["'r0=0.1;slope=-4;L0=30;nx=2048;dx=1./64;SURFEVL=1;SURFWFS=0'"]
```
The double quote is necessary here to group the single quoted entries together.

The amplitude map of the telescope can be specified with
`aper.fnamp=aper.bin` with a similar data format, with OPD replaced by
amplitude map between 0 and 1. By default, NCPA are calibrated (`ncpa.calib=1`
uses all DMs, 2 uses only the ground DM). Use `ncpa.calib=0` to disable calibration.

<!---
We can also setup one or more tilted M3 surfaces that are common to science
fields and wavefront sensors. Put the list of surface files in key \c
tsurf. Each surface file must contain a 2-d array with a header specifying
the following keys:

```
dx    the sampling in x direction (first dimension).
dy    the sampling in y direction (second dimension).
ox    the origin in x direction.
oy    the origin in y direction.
txdeg the x tilt angle in degrees wrt beam (90 is perp)
tydeg the y tilt angle in degrees wrt beam (90 is perp)
ftel  the telescope effective focal length
fexit the distance between the exit pupil and the focus
fsurf the distance between the center of the M3 surface and the focus.
```
-->
\subsection sect-wfs WFS Configuration

The WFS configuration is mostly supplied in \c powfs, which applies to all
the wfs belonging to this type. For wfs specific parameters, use \c wfs. If
there is a single number in \c powfs.[key], it applies to all powfs. The
following lists a few common options

```
powfs.noisy=[0]         #set noise free for all powfs
powfs.noisy=[1 0 0]     #set noise free only for powfs 1 and 2.
powfs.phystep=[0 -1 -1] #use physical optics wfs for the first type (high order), and geometric for the remaining (t/t/f)
powfs.phytype_recon=[1 2 2] #pixel processing 1: use matched filter, 2: cog
powfs.siglev=[900 100 100]  #signal level per subaperture per frame at sim.dtref for physical optics wfs mode.
powfs.bkgrnd=[10 0 0]   #background level per pixel per frame at sim.dtref
powfs.rne=[3 1 1]       #read out noise for physical optics wfs mode
powfs.nearecon=[20 2 1] #noise equivalent angle in milli-arcsecond for geometric wfs mode
powfs.neasim=[20 2 1]   #nea for simulation. -1 to match nearecon.
```

WFS specific parameters usually include WFS coordinates
```
wfs.thetax=[]           #specify the x coordinate in arcsec
wfs.thetay=[]           #specify the y coordinate in arcsec
```

\subsection sect-perfevl Point Spread Function

MAOS only computes RMS averaged wavefront errors by default, which are saved to
\c Res_[seed].bin and \c Resp_[seed].bin. When desired, PSF computation can be
enabled for some or all of the science evaluation directions (`evl.psf`). Modify \c
evl.thetax , \c evl.thetay , and \c evl.wt to modify evaluation directions and
relative weighting for field averaging. Notice that this will significantly slow
down the simulation. Use the following parameters to enable PSF computation.
```
maos evl.psfisim=20 evl.psfmean=1 # start averaging from step 20, save averaged PSF once in the end.
maos evl.psfmean=100            # also save the averaged PSF every 100 steps.
maos evl.psfmean=1 evl.psfol=1  # also include open loop PSF
maos evl.psfmean=1 evl.psfsize=[1024]     # only save center part 1024x1024 of the PSF
maos evl.psfmean=1 evl.psfgridsize=[4096] # specify PSF FFT grid size to control the sampling.
```
`evl.psfsize` and `evl.psfgridsize` take one value for all wavelengths or one
value per wavelength in `evl.wvl`.

\subsection sect-act Actuator Slaving

Due to pupil obscuration, there may be actuators that are not sensed by the
WFS and therefore cannot be accurately reconstructed. The actuator coupling
coefficient, computed using DM fitting ray tracing (HA in minimum variance
reconstruction) or gradient interaction matrix (GA in least squares
reconstruction), is used to determine how accurately an actuator is
sensed/controlled. Activated by

      lsr.actslave=1 #or 2, for lsr (default 1)
      fit.actslave=1 #or 2, for mvr (default 1)

- \c actslave=1 enables implementation with threshold \c .actthres
(`lsr.actthres=0.7`, `fit.actthres=0.1` by default). When an
actuator is outside of the pupil, its coupling coefficient is usually very
small. When it is below the threshold, its value is slaved to its neighboring
actuators.
- \c actslave=2 in addition enables implementation with threshold \c
.actthres2 (`lsr.actthres2=0.8`, `fit.actthres2=0.7`). When an actuator is hidden by secondary mirror support struts,
it may still have relatively strong coupling, but regions separated by the
struts may have different piston values, the so called island effect. The
mitigation method is to limit the value gap around such actuators whose
coupling is below .actthres2.

\subsection sect-sodium Sodium range variation

The sodium range variation is specified through the LLT configuration key \c
llt.fnrange in \ref llt_CL or \ref llt_SL, which is selected by
`powfs.fnllt`. To override it, use \c powfs0_llt.fnrange, where `0` is the index
of the LGS powfs. The file should contain additional change (wrt fnprof) of
sodium height in meters with dimension nx1 or nxnLGS where \c n can be less
than, equal to or greater than the length of simulation and nLGS is the
number of LGS in this powfs (6 for NFIRAOS).

The range variation is simulated by adding corresponding focus mode to LGS WFS wavefront.

\subsection pixel-background Rayleigh backscatter

The Rayleigh backscatter effect can be introduced using `powfs.bkgrndfn=[filename,]`
which contains the background level for each pixel per frame, defined at
`sim.dtref`. The file contains a cell array of `nsa*1` (shared by all WFS) or
`nsa*nwfs`, where `nwfs` is the number of `wfs` for this `powfs`, and each cell
contains a `npixpsax*npixpsay` numerical array. The subapertures are ordered as
in `powfs0_saloc.bin`. `powfs.bkgrndfnc` specifies the part that is subtracted
(calibrated).

\section maos_skycoverage Sky coverage

The sky coverage simulation is done in two parts. \c maos is run first to prepare NGS mode and WFS PSF time series:

\c maos -c nfiraos_lgs.conf skyc_10.conf -o 50p_za0_ncpa

The second step is to run \c skyc in the sky sim folder
\c cd 50p_za0_ncpa/skysim
\c skyc -d -o APD_JHKs_typeImr

The results are stored in \c Res1_1.bin for \c maos seed 1, \c skyc seed 1. There are \c skyc.nsky (500 by default) columns with the following row definitions:

1. Total residual error
2. Residual atmospheric (included in \c maos simulations) error in NGS modes
3. Residual atmosphere error only in tip/tilt modes.
4. Residual windshake if `skyc.psd_ws` is set and `skyc.addws=0`.
5. Estimated error.
6. Residual atmospheric error in NGS modes on axis only.
7. Residual atmospheric error in tip/tilt modes on axis only.
8. Selected NGS WFS sampling frequency.
9. Index into \c Res1_1_aster.bin.
10. Index of the selected asterism.

See \ref skycoverage for details and \ref page40_results for further explanation and results interpretation.
