/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \page page30 Simulation Parameters
   
    All user modifiable options are configured through these config files with
    suffix \c .conf. The supplied default setup files are stored in sub-folder
    \c config/maos. The config files are organized in a modular way. One \c
    .conf file can include another \c .conf file with \c
    include=”filename.conf”.  Most of the \c .conf files in \c config/maos file
    are self-documented.

    All the options are in the form of \c key=value. The usually convention is:
    if the entry in `parms` called `parms->wfs.thetax`, the entry in .conf file will
    be in the form `wfs.thetax=value`. 

    Some of the rules are:

    -# Array values are preferred to be embraced by []
    -# string are embraced by double or single quotes ('' or "").  
    -# A line ended with \ will be concatenated with the next line. 
    -# Anything after comment string # will be ignored during reading. 
    
    The config `key=value` can be specified multiple times for any key, in which
    case the latter value will override the previous value. The software will
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
    
    In the \c config/maos sub-folder in the source folder, there are a few
    envelop \c .conf files that is complete for each configuration:
    
    -# \c `nfiraos.conf` for NFIRAOS mcao simulations,
    -# \c `mcao_ngs.conf` for NGS mcao simulations
    -# \c `scao_lgs.conf` for single conjugate LGS simulation
    -# \c `scao_ngs.conf` for single conjugate NGS simulation. 

    Each of this file include a few other \c .conf files by specifying \c
    include=`filename.conf`, for example, all these baseline \c .conf files
    include the following common `.conf` files that is independent on the
    simulation mode or geometry.

    -# \c `sim.conf` to specify general simulation parameters such as loop gain, number of time steps, etc,
    -# \c `recon.conf` to specify reconstruction parameters in closed loop, 
    -# \c `dbg.conf` to specify debugging parameters. 

    The following list the exclusive `.conf` files used to specified particular
    set of configurations. Each group should only appear once.
    
    - `atm_*.conf`: set to a different preconfigured turbulence profile.
    
    - `dm_*.conf`: set deformable mirror configuration.

    - `fov_*.conf`: set the science FoV.

    - `wfs_*.conf`: set the wfs type (excludes `wfs_common.conf`)

    For example, the baseline configuration for TMT MCAO system `nfiraos.conf`
    includes `atm_mk13n50p.conf` to specify MK13N median seeing turbulence
    profile, `dm_dual.conf` to specify dual DM system, `wfs_lgs_ttf_tt.conf` to
    specify WFS configuration as LGS WFS plus TTF NGS WFS plus TT NGS WFS,
    `fov_iris.conf` to specify the science FoV as simpson weighting on 3x3
    points in the square FoV.

    \subsection sect-surface Specifying Surface OPDs

    We can optionally setup one or more surfaces that can cover science fields
    or wavefront sensors. Put the list of surface `.bin` files in key \c
    surf. Each surface file must contain a 2-d array of the OPD with a header
    specifying the following keys:

    \verbatim
    dx  the sampling in x direction (first dimension).
    dy  the sampling in y direction (second dimension).
    ox  the origin in x direction.
    oy  the origin in y direction.
    h   the height conjugation of this surface
    vx  the frozen flow wind speed along x of this surface (0 for static)
    vy  the frozen flow wind speed along y of this surface (0 for static)
    \endverbatim
    
    Use the write mex routine or writebin.m to write the bin file:

    `write(OPD, header, 'file.bin')`
    
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

    \ref nfiraos Shows the example configuration file for NFIRAOS.
    
    \ref SCAO_NGS shows the example configuration file of a single conjugate NGS system.

    \ref SCAO_LGS shows the example configuration file of a single conjugate LGS system.
*/
