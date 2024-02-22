\page page30_run Run simulations
\tableofcontents
\section sect-run Usage

We assume that `maos` is already in the PATH so we can type `maos` to launch it.

Use the following command to get the help message:
```
maos -h  #will print the help message and exit.
```

The valid command line arguments are listed as follows

```
Usage: maos [OPTION...] [FILE]...
maos is a simulation tool developed to adaptive optics systems

Examples:
maos -o result # Run the default configuration from default.conf and save results to \c result folder

maos -c scao_ngs.conf sim.seeds=[1 10] -d -o scao_ngs override.conf
         # Run a single conjugate natural guide star case, with seed 1 and 10
         # detach from the terminal and output results to folder scao_ngs
         # and read in overriding parameters stored in override.conf and chol.conf

Options:
-h, --help        to print out this message
-d, --detach      to detach from terminal and run in background
-f, --force       force starting simulation without scheduler
-n, --nthread=N   Use N threads, default is 1
-o, --output=DIR  output the results to DIR.
-c, --conf=FILE.conf   Use FILE.conf as the baseline config instead of nfiraos.conf
-p, --path=dir    Add dir to the internal PATH
-g, --gpu=i       Use the i'th gpu. 0 for the first. -1 to disable. default: automatic
-G, --ngpu=N      Use a total of N gpus.
-r, --run=host    Run the job in another host. Requires scheduler to be running.

The following environment variables are supported
MAOS_TOMOSCALE=1e12 Rebalance tomography terms for single precision calculation
MAOS_PARALLEL=1      Set to 0 to disable parallel launch
MAOS_NO_WFS=0        Set to 1 to disable all WFS calls
MAOS_NO_EVL=0        Set to 1 to disable evaluation calls
MAOS_NO_RECON=0      Set to 1 to disable reconstruction calls
MAOS_KEEP_MEM=0      Set to 1 to keep some temporary memory between steps.
MAOS_MEM_DEBUG=0     Set to 1 to enable malloc/free accounting. Useful to detect memory leak.
MAOS_MEM_VERBOSE=0   Set to 1 to print detailed malloc/free info
MAOS_LOG_LEVEL=0     Set logging level. -3: error and warning only,
                           -2:  essential info, -1  useful info, 0:  all info,
                           1:  debugging info, 2: more debugging info, 3: everything.
MAOS_GPU_DISABLE_0=1         Do not utilize GPU 0 (sudo nvidia-smi -i 0 -c 2 has the same effect)
```

\section sect-config Configuration Files

All user modifiable options are configured through config files with suffix
\c .conf or key=value in the command line. The supplied default setup files
are stored in sub-folder \c config/maos. The config files are organized in a
modular way. One \c .conf file can include another \c .conf file with \c
include=filename.conf.  The \c .conf files in \c config/maos file
are self-documented.

All the options are in the form of \c key=value. The usually convention is:
if the entry in `parms` called `parms->wfs.thetax`, the entry in .conf file will
be in the form `wfs.thetax=value`.

Some of the rules are:

- array values must be embraced by \c [] and separated by \c ,
- string are embraced by double or single quotes ('' or "").
- A line ended with \ will be concatenated with the next line.
- Anything after comment string # will be ignored during reading.

The config `key=value` can be specified multiple times for any key, in which
case the latter value will override the previous value, except when \c
key+=value is used which causes it to append instead. The software will
complain and exit if any `key=value` pair is not understood or required but
not provided.

The recommended way to organize config files for a specific problem is to
keep a baseline \c .conf file in sub-folder \c config/maos in the source
folder that specify all the necessary \c key=value entries through
optionally including other \c .conf files. The default entry is
`default.conf` unless `-c file.conf` switch is supplied in the command line.

Embed the overriding options in the command line (embrace it with quote if
the option contains spaces, like arrays) or put in one or multiple
overriding file (in the simulation folder, not inside -c config/maos) to
override the values you want to change from the baseline.

The \c .conf files in the \c config/maos folder are only modified during
code development and should not be modified by the end user. Use overriding
\c .conf files in simulation folders or embed configurations in the command
line to modify the values.

In the \c config/maos folder, there are a few
sample \c .conf files that is complete for each configuration:

- \ref mcao_lgs for NFIRAOS mcao simulations (the default),
- \ref mcao_ngs for NGS mcao simulations
- \ref scao_lgs for single conjugate LGS simulation
- \ref scao_ngs for single conjugate NGS simulation.
- \ref scao_pywfs for Pyramid WFS NGS simulation.

Each of this file include a few other \c .conf files by specifying \c
include=`filename.conf`, for example, all these baseline \c .conf files
include the following common `.conf` files that is independent on the
simulation mode or geometry.

- \ref sim to specify general simulation parameters such as loop gain, number of time steps, etc,
- \ref recon to specify reconstruction parameters in closed loop,
- \ref dbg to specify debugging parameters.

The following list the exclusive `.conf` files used to specified particular
set of configurations. Each group should only appear once.

- `atm_*.conf`: set to a different preconfigured turbulence profile.
- `dm_*.conf`: set deformable mirror configuration.
- `fov_*.conf`: set the science FoV.
- `wfs_*.conf`: set the wfs group type

For example, the baseline configuration for TMT MCAO system `mcao_lgs.conf`
includes `atm_mk13n50p.conf` to specify MK13N median seeing turbulence profile,
`dm_dual.conf` to specify dual DM system, `wfs_lgs_ttf_tt.conf` to specify WFS
configuration as LGS WFS plus TTF NGS WFS plus TT NGS WFS, `evl_sq34.conf` to
specify the science FoV as simpson weighting on 3x3 points in the square 34x34"
FoV, and `fit_sq34.conf` to sepcify the DM fitting FoV as simpson weighting on
3x3 points in the square 34x34" FoV. 

For an MCAO configuration, the DM fitting FoV can be larger than the science FoV
in order to sharpen diffraction limited TTF NGS. For an LTAO configuration, use
`fit_oa.conf` to set only on axis DM fitting FoV. A GLAO system can be simulated
by having a single DM but with a wide DM fitting FoV, as an alternative to
setting `recon.glao=0` which averages WFS measurements before reconstruction.

See \ref page33_example for more detailed explanations.

\section sect-exe Sample Runs

To run predefined AO modes:

```
maos -o mcao # default is MCAO with 2 dms. Save results to folder mcao. 
maos -c mcao_ngs.conf  # NGS MCAO
maos -c scao_ngs.conf  # NGS SCAO
maos -c scao_pwfs.conf # NGS SCAO using PWFS
maos -c scao_lgs.conf  # LGS SCAO
maos -c examples/keck_lgs.conf #KECK LGS AO
maos -c mcao_lgs.conf dm_single.conf evl.moao=0 moao.dx=[1/2] #LGS LTAO with MOAO for each science field.
maos -c mcao_lgs.conf dm_triple.conf #LGS MCAO with 3 DMs
```
The default case starts with default.conf. It can be overriden by the environment variable MAOS_DEFAULT.

To customize number of components

```
maos dm_single.conf wfs_lgs_only.conf recon.glao=1 # LGS GLAO
maos dm_triple.conf dm.dx=0.3  # triple DM with 0.3 actuator pitch
maos wfs_lgs_ttf_tt_twfs.conf  # with LGS, TTF, TT, and TWFS
maos evl_x.conf #evaluate performance along x axis.
```

Change aperture

```
maos aper.d=[10]    # 10 meter circular
maos aper.d=[10, 1] # 10 meter annular diameter with 1 meter inner obscuration
maos aper.d=[10, 1] aper.fnamp=aper.bin # with supplied amplitude map. Diameters need to match.
```

Change turbulence

```
maos atm_mk13n25p.conf #use another predefined profile. The default is atm_mk13n50p.conf
maos atm.r0z=0.1     # change r0 at zenith to 0.1m
maos sim.zadeg=30    # change zenith angle to 30 degrees
maos atm_single.conf # use single ground layer turbulence
```

Configuring WFS. Adjust the number of elements depending on how many powfs is in use. A powfs
is a type of WFS. Each type can have multiple WFS.

```
maos powfs.noisy=[0]     # use noise free WFS for all wfs types
maos powfs.noisy=[1 0 0] # use noise free WFS for wfs type 1 and 2.
maos powfs.phystep=[-1]  # use geometric wfs instead of physical optics wfs for all
```



Adjust the controller

```
maos sim.ephi=0.3 #Change gain of the main integrator to 0.3
```
\section sect-ebb Error Budget Breakdown

There are options to allow probing for error budget items, such as DM fitting error, tomography error, etc. 

```
maos scao_ngs.conf sim.closeloop=0 atm.frozenflow=1 sim.idealfit=1 # step1: single DM fitting on axis
maos mcao_ngs.conf sim.closeloop=0 atm.frozenflow=1 sim.idealfit=1 # step1b: Multiconjugate DM fitting over a field
maos scao_ngs.conf sim.closeloop=0 atm.frozenflow=1 powfs.phystep=[-1] powfs.noisy=[0] # step2: noise free classic AO simulation with geometric SHWFS.
maos mcao_ngs.conf sim.closeloop=0 atm.frozenflow=1 powfs.phystep=[-1] powfs.noisy=[0] # step2b: noise free MCAO with geometric SHWFS.
maos mcao_ngs.conf powfs.phystep=[-1] powfs.noisy=[0] # step3: noise free MCAO with geometric SHWFS in closed loop simulation.
maos mcao_ngs.conf powfs.noisy=[0] # step4: noise free MCAO with physical optics WFS.
maos mcao_lgs.conf powfs.noisy=[1] # step5: noisy MCAO with physical optics LGS WFS.

```

Note that the option `sim.closeloop=0 atm.frozenflow=1` uses open loop correction (no servo lag) with frozen flow turbulence. The option `sim.idealfit=1` enables DM fitting directly from turbulence rather than from tomography output. For open loop simulations, truth wfs should not be included or disabled with sim.eptwfs=0, and non-iterative solvers (tomo.alg=0, fit.alg=0) usually work better than CG.

The following error budget items are computed:

- `DM fitting error`: step1 gives DM fitting error which depends on the actuator spacing.
- `DM projection error`: the quadrature difference between step1 and step1b gives DM projecting error which depends on the field of view and number of DMs.
- `WFS aliasing error`: the quadrature difference between step1 and step2 gives WFS aliasing error.
- `Tomography error`: the quadrature difference between the RSS of the three terms above and step2b is the tomography error which depends on the asterism geometry. 
- `Servo lag error`:  the quadrature difference between step2b and step3 gives servo lag error.
- `WFS nonlinearity error`: the quadrature difference between step3 and step4 gives WFS non-linearity error.
- `WFS noise effect`: the quadrature difference between step4 and step5 gives WFS noise error.


\section advanced Advanced configuration

\subsection sect-surface Specifying Surface OPDs

We can optionally setup one or more static surfaces that cover science fields
and/or wavefront sensors. Each surface file must contain a 2-d array of the
OPD with a header specifying the following keys. If header is not available,
the OPD is assumed to have 1/64 sampling and centered on the pupil.

```
dx  the sampling in x direction (first dimension).
dy  the sampling in y direction (second dimension).
ox  the origin in x direction.
oy  the origin in y direction.
h   the height conjugation of this surface
vx  the frozen flow wind speed along x of this surface (0 for static)
vy  the frozen flow wind speed along y of this surface (0 for static)
```

Optionally, the header can also include the following keys to indicate its coverage

```
SURFNAME=name       #name of the surface. M1 or M2 for the primary or secondary mirror
SURFEVL=[1 1 1 ...] #length: nevl. 1: enabled for this science evaluation direction (assume all 1 if omitted)
SURFWFS=[1 1 1 ...] #length: nwfs. 1: enabled for this WFS (assume all 1 if omitted)
```

Use the \c write mex routine, \c writebin.m, or \c aolib.writebin to write the bin file:
```
   write(OPD, header, 'opd1.bin')
```
Or simply use fits format. Put the list of surface file names in key \c surf.
```
   maos surf=['opd1.bin','opd2.bin', 'opd3.fits']
```

It is also possible to specify NCPA with inline configuration:
```
maos surf=["'r0=0.1;slope=-4;SURFEVL=1;', 'r0=0.2; slope=-4;SURFWFS=1 1 1 1 1 1 0 0 0'"]
maos surf=["'r0=0.1;slope=-4;L0=30;nx=2048;dx=1./64;SURFEVL=1;SURFWFS=0'"]
```
The double quote is necessary here to group the single quoted entries together.

The amplitude map of the telescope can be specified with
`aper.fnamp=aper.bin` with a similar data format, with OPD replaced by
amplitude map between 0 and 1. By default, NCPA are calibrated, use `sim.ncpa_method=0` to disable calibration

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

\subsection sect-wfs WFS Configuration

The WFS configuration are mostly supplied in \c powfs, which applies to all
the wfs belonging to this type. For wfs specific parameters, use \c wfs. If
there is a single number in \c powfs.[key], it applies to all powfs. The
follow lists a few common options

```
powfs.noisy=[0]         #set noise free for all powfs
powfs.noisy=[1 0 0]     #set noise free only for powfs 1 and 2.
powfs.phystep=[0 -1 -1] #use physical optics wfs for the first type (high order), and geometric for the remaining (t/t/f)
powfs.phytype=[1 2 2]   #pixel processing 1: use matched filter, 2: cog
powfs.siglev=[900 100 100] #specify signal level for physical optics wfs mode at sim.dtref Hz per pixel.
powfs.bkgrnd=[10 0 0]   #specify background level at sim.dtref Hz per pixel
powfs.rne=[3 1 1]       #specify read out noise for physical optics wfs mode
powfs.nearecon=[20 2 1] #specify noise equivalent angle in milli-arcsecond for geometric wfs mode
powfs.neasim=[20 2 1]   #specify nea for simulation. -1 to match nearecon.
```

WFS specific parameters usually include WFS coordinate
```
wfs.thteax=[]           #specify the x coordinate in arcsec
wfs.thetay=[]           #specify the y coordinate in arcsec
```

\subsection sect-perfevl Point Spread Function

Maos only computes RMS averaged wavefront errors by default, which are saved
to \c Res_[seed].bin. When desired, PSFs computing can be enabled for some or all
of the science evaluation directions. Modify \c evl.thetax , \c evl.thetay , and \c
evl.wt to modify evaluation directions and relative weighting for field
averaging. Use the following parameters to enable PSF computation.
```
maos evl.psfisim=20 evl.psfmean=1 # start averaging from step 20, save averaged PSF once in the end.
maos evl.psfmean=1 evl.psfol=1   # also include open loop PSF
maos evl.psfmean=1 evl.psfsize=[1024,1024]     # only save center part 1024x1024 of the PSF
maos evl.psfmean=1 evl.psfgridsize=[4096,4096] # specify PSF FFT grid size to control the sampling.
```

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

\subsection sect-sodium Sodium range variation

The sodium range variation is specified through the LLT configuration key \c
llt.fnrange in \ref llt_CL or \ref llt_SL. When supplied by user, use \c
powfs0_llt.fnrange. The file should contain additional change (wrt fnprof) of
sodium height in meters with dimension nx1 or nxnLGS where \c n can be less
than, equal to or greater than the length of simulation and nLGS is the
number of LGS in this powfs (6 for NFIRAOS).

The range variation is simulated by adding corresponding focus mode to LGS WFS wavefront.

\subsection pixel-background Rayleigh backscatter

The Rayleigh backscatter effect can be introduced using `powfs.bkgrndfn=[filename,]` which contains simulation background level for each pixel defined at `sim.dtref`. The file contains a cell array of `nsa*nlgswfs` while each cell contains a `npixpisax*npixpsay` numerical array. `nlgswfs` is the number of `wfs` for this `powfs`. 

\section skycoverage Sky coverage

The sky coverage simulation is done in two parts. \c maos is run first used to prepare NGS mode and WFS PSF time series:

\c maos -c nfiraos_lgs.conf skyc_10.conf -o 50p_za0_ncpa

The second step is to run \c skyc in the sky sim folder
\c cd 50p_za0_ncpa/skysim
\c skyc -d -o APD_JHKs_typeImr

The results are stored in \c Res1_1.bin for \c maos seed 1, \c skyc seed 1. There are \c nsky=500 columns with the following row definitions:

1. Total residual error
2. Residual atmpheric (included in \c maos simulations) error
3. Residual atmosphere error only in tip/tilt modes.
4. Residual windshake if skyc.psd_is is set and skyc.addws=0. (deprecated use)
5. 0 (not used)

===

See \ref page40_results for further explanation and results interpretation.
