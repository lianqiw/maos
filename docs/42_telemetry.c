/**
   \page page42 Telemetry Data

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

<table border="0" cellspacing="0" cellpadding="2">
           <tr><td>File name                </td><td>Option to enable</td><td>Description</td>
</td></tr><tr><td>\c atm_1              </td><td>save.atm   </td><td>The atmosphere 
</td></tr><tr><td>\c dmerr_hi_1      </td><td>save.dm   </td><td>The DM error signal for high order wfs
</td></tr><tr><td>\c dmfit_hi_1      </td><td>save.dm   </td><td>The DM fit result for high order wfs
</td></tr><tr><td>\c Merr_lo_1       </td><td>save.dm   </td><td>The low order mode error signal (split tomography)
</td></tr><tr><td>\c Mint_lo_1       </td><td>save.dm   </td><td>The low order mode integrator output (split tomography)
</td></tr><tr><td>\c dmint_hi_1      </td><td>save.dm   </td><td>The DM integrator output of high order wfs output (split integrator only)
</td></tr><tr><td>\c dmint_1         </td><td>save.dm   </td><td>The DM integrator output (command integrator for both high/low wfs)
</td></tr><tr><td>\c dmreal_1        </td><td>save.dm   </td><td>The real, effective DM commands applied
</td></tr><tr><td>\c dmpttr_1        </td><td>save.dmpttr   </td><td>The piston/tip/tilt removed, effective, DM commands.
</td></tr><tr><td>\c evl0_opdcl_1    </td><td>save.evlopd   </td><td>The closed loop opd for science fields
</td></tr><tr><td>\c evl0_opdol_1    </td><td>save.evlopd   </td><td>The open loop opd for science fields
</td></tr><tr><td>\c gcov_wfs0_5_10_1      </td><td>save.gcov   </td><td>The covariance between gradients of wfs 0 and 5 saved at time step 10
</td></tr><tr><td>\c opdr_1          </td><td>save.opdr   </td><td>The tomography output, defined on xloc
</td></tr><tr><td>\c opdx_1          </td><td>save.opdx   </td><td>The atmosphere projected onto xloc (direct fitting)
</td></tr><tr><td>\c wfs0_gradcl_1   </td><td>save.grad  </td><td>The wfs gradient measurement.
</td></tr><tr><td>\c wfs0_gradnf_1   </td><td>save.grad  </td><td>The wfs noise free gradient.
</td></tr><tr><td>\c wfs0_gradpsol_1 </td><td>save.grad  </td><td>The wfs pseudo open loop gradient
</td></tr><tr><td>\c wfs0_gradgeom_1 </td><td>save.gradgeom   </td><td>The wfs geometric gradient (in physical optics wfs simulations)
</td></tr><tr><td>\c wfs0_intsnf_1   </td><td>save.ints   </td><td>The wfs subaperture images (noise free)
</td></tr><tr><td>\c wfs0_intsny_1   </td><td>save.ints   </td><td>The wfs subaperture images (noisy)
</td></tr><tr><td>\c wfs0_lltopd_1   </td><td>save.wfsopd   </td><td>The wfs laser launch telescope OPD.
</td></tr><tr><td>\c wfs0_opd_1      </td><td>save.wfsopd   </td><td>The wfs OPD.
</td></tr>
*/
