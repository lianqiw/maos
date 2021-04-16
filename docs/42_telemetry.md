\page page40_results
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

