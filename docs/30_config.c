/**
   \page page3 Simulation Parameters
   
    All user modifiable options are configured through these config files with
    suffix \c .conf. The supplied default setup files are stored in sub-folder
    \c config/maos. The config files are organized in a modular way. One \c
    .conf file can include another \c .conf file with \c
    include=”filename.conf”.  Most of the \c .conf files in \c config/maos file
    are self-documented.

    All the options are in the form of \c key=value. The usually convention is:
    if the entry in parms called parms.wfs.thetax, the entry in .conf file will
    be in the form wfs.thetax=value. 

    Some of the rules are:

    -# Array values are preferred to be embraced by []
    -# string are embraced by double or single quotes ('' or "").  
    -# A line ended with \ will be concatenated with the next line. 
    -# Anything after comment string # will be ignored during reading. 
    
    The config key=value can be specified multiple times for any key, in which
    case the latter value will override the previous value. The software will
    complain and exit if any key=value pair is not understood or required but
    not provided.

    The recommanded way to organize config files for a specific problem is to
    keep a baseline \c .conf file in sub-folder \c config/maos in the source
    folder that specify all the necessary \c key=value entries through
    optionally including other \c .conf files. 

    Embed the overriding options in the command line (embrace it with quote if
    the option contains spaces, like arrays) or put in one or multiple
    overriding file (in the simulation folder, not inside -c config/maos) to
    override the values you want to change from the baseline.

    The \c .conf files in the \c config/maos folder are only modified during
    code development and should not be modified by the end user. Use overriding
    \c .conf files in simulation folders or embed configurations in the command
    line to modify the values.
    
    In the \c config/maos sub-folder in the source folder, there are three
    baseline \c .conf files:
    
    -# \c nfiraos.conf for NFIRAOS mcao simulations,
    -# \c scao_lgs.conf for single conjugate LGS simulation
    -# \c scao_ngs.conf for single conjugate NGS simulation. 

    Each of this file include a few other \c .conf files by specifying
    \c include=”filename.conf”, for example, all the three baseline \c .conf files
    include

    -# \c “sim.conf” to specify general simulation parameters such as loop gain, number of time steps, etc,
    -# \c “recon.conf” to specify reconstruction parameters in closed loop, 
    -# \c “dbg.conf” to specify debugging parameters. 

    The MCAO baseline “nfiraos.conf” includes “atm_mk13n50p.conf” to specify
    MK13N median seeing turbulence profile, while the two SCAO baselines include
    “atm_scao.conf” that specifies the turbulence information with only a single
    reconstruction layer.

    We can optionally setup one or more surfaces that can cover science fields
    or wavefront sensors. Put the list of surface files in key \c surf. Each
    surface file must contain a 2-d array with a header specifying the following
    keys:

    dx  the sampling
    ox  the origin in x direction.
    oy  the origin in y direction.
    h   the height conjugation of this surface
    vx  the frozen flow wind speed along x of this surface (0 for static)
    vy  the frozen flow wind speed along y of this surface (0 for static)

    
    We can also setup one or more tilted M3 surfaces that are common to science
    fields and wavefront sensors. Put the list of surface files in key \c
    tsurf. Each surface file must contain a 2-d array with a header specifying
    the following keys:
    
    dx    the sampling in x direction (first dimension).
    dy    the sampling in y direction (second dimension)
    ox    the origin in x direction.
    oy    the origin in y direction.
    txdeg the x tilt angle in degrees wrt beam (90 is perp)
    tydeg the y tilt angle in degrees wrt beam (90 is perp)
    ftel  the telescope effective focal length
    fexit the distance between the exit pupil and the focus
    fsurf the distance between the center of the M3 surface and the focus.

    \ref nfiraos Shows the example configuration file for NFIRAOS.
    
    \ref SCAO_NGS shows the example configuration file of a single conjugate NGS system.

    \ref SCAO_LGS shows the example configuration file of a single conjugate LGS system.
*/
