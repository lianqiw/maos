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
   \section sect-run Usage
   
   We assume that maos is already in the PATH so we can type "maos" to launch it.
   
   Use the following command to get the help message:
   \verbatim
   maos -h  #will print the help message and exit.
   \endverbatim

   The valid command line arguments are listed as follows

   \verbatim
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
   -r, --run=host    Run the job in another host.

   The following environment variables are supported
   MAOS_TOMO_SCALE=1e12 Rebalance tomography terms for single precision calculation
   MAOS_PARALLEL=1      Set to 0 to disable parallel launch
   MAOS_NO_WFS=0        Set to 1 to disable all WFS calls
   MAOS_NO_EVL=0        Set to 1 to disable evaluation calls
   MAOS_NO_RECON=0      Set to 1 to disable reconstruction calls
   MAOS_KEEP_MEM=0      Set to 1 to keep temporary memory between steps
   MAOS_MEM_DEBUG=0     Set to 1 to enable malloc/free accounting
   MAOS_MEM_VERBOSE=0   Set to 1 to print detailed malloc/free info
   MAOS_LOG_LEVEL=0     Set logging level. -3: error only, -2:  warning, -1  essential info, 0:  useful info,
                        1:  debugging info, 2: debugging info, 3: everything.

   \endverbatim

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

   The recommanded way to organize config files for a specific problem is to
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
   - \ref scao_pwfs for Pyramid WFS NGS simulation.

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
   includes `atm_mk13n50p.conf` to specify MK13N median seeing turbulence
   profile, `dm_dual.conf` to specify dual DM system, `wfs_lgs_ttf_tt.conf` to
   specify WFS configuration as LGS WFS plus TTF NGS WFS plus TT NGS WFS,
   `fov_sq34.conf` to specify the science FoV as simpson weighting on 3x3
   points in the square 34x34" FoV.

   See page31_example for more detailed explanations.

   \section sect-exe Sample Runs

   To run predefined AO modes
   
   \verbatim
   maos -o mcao # default is dual conjugate AO, save results to folder mcao
   maos -c mcao_ngs.conf  # NGS MCAO
   maos -c scao_ngs.conf  # NGS SCAO
   maos -c scao_pwfs.conf # NGS SCAO using PWFS
   maos -c scao_lgs.conf  # LGS SCAO
   \endverbatim
   
   To customize 

   \verbatim
   maos dm_single.conf wfs_lgs_only.conf recon.glao=1 # LGS GLAO
   maos dm_triple.conf dm.dx=0.3  # triple DM with 0.3 actuator pitch
   maos wfs_lgs_ttf_tt_twfs.conf  # with LGS, TTF, TT, and TWFS
   maos evl_x.conf #evaluate performance along x axis.
   \endverbatim
   
   Change aperture

   \verbatim
   maos aper.d=[10]    # 10 meter circular
   maos aper.d=[10, 1] # 10 meter annular diameter with 1 meter inner obscuration
   maos aper.d=[10, 1] aper.fnamp=aper.bin # with supplied amplitude map. Diameters need to match.
   \endverbatim

   Change turbulence

   Override keys listed in \c atm_mk13n50p.conf

   \verbatim
   maos atm.r0z=0.1     # change r0 at zenith to 0.1m 
   maos sim.zadeg=30    # change zenith angle to 30 degrees
   maos atm_single.conf # use single ground layer turbulence
   \endverbatim

   Configuring WFS

   \verbatim
   maos powfs.noisy=[0]     # use noise free WFS for all wfs types
   maos powfs.noisy=[1 0 0] # use noise free WFS for wfs type 1 and 2.
   maos powfs.phystep=[-1]  # use geometric wfs instead of physical optics wfs for all
   \endverbatim

   Adjust the number of elements depending on how many powfs is in use. A powfs
   is a type of WFS. Each type can have multiple WFS.

   Adjust the controller
   
   \verbatim
   maos sim.ephi=0.3 #Change gain of the main integrator to 0.3
   \endverbatim
 
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

   \subsection Sodium range variation

   The sodium range variation is specified through the LLT configuration key \c
   llt.fnrange in \ref llt_CL or \ref llt_SL. When supplied by user, use \c
   powfs0_llt.fnrange. The file should contain additional change (wrt fnprof) of
   sodium height in meters with dimension nx1 or nxnLGS where \c n can be less
   than, equal to or greater than the length of simulation and nLGS is the
   number of LGS in this powfs (6 for NFIRAOS).

===

   See \ref page40_results for further explanation and results interpretation.
*/
