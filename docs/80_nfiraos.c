/**
   \page nfiraos Example: NFIRAOS
   The following shows the content of nfiraos.conf, which is the baseline
    configuration of TMT NFIRAOS system
    
    \include config/maos/nfiraos.conf
    
    It includes the following files to define all the parameters.
    \c "sim.conf" to define some general parameters:
    \include config/maos/sim.conf
    
    \c "aper_tmt.conf" to define the TMT aperture
    \include config/maos/aper_tmt.conf

    \c "dm_dual.conf" to define the dual deformable mirrors
    \include config/maos/dm_dual.conf

    \c "atm_mk13n50p.conf" to define the turbulence profiles using the Mauna Kea 13N median seeing
    \include config/maos/atm_mk13n50p.conf

    \c "wfs_multi_lgs.conf" to define the wavefront sensors.
    \include config/maos/wfs_multi_lgs.conf
    \include config/maos/llt_CL.conf

    \c "fov_iris.conf" to define the science fov for performance evaluation with IRIS.
    \include config/maos/fov_iris.conf

    \c "recon.conf" to define the parameters about turbulence reconstruction and DM fitting
    \include config/maos/recon.conf

    \c "dbg.conf" to define debugging and saving/loading geometry and telemetry.
    \include config/maos/dbg.conf

*/
