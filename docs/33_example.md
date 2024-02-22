\page page33_example Examples
\tableofcontents
- \subpage nfiraos_example
- \subpage scao_base_example
- \subpage scao_lgs_example
- \subpage scao_ngs_example
- \subpage page34_conf


\page nfiraos_example TMT NFIRAOS LGS MCAO

The default configuration \ref default uses \ref mcao_lgs which describes the
baseline configuration of TMT NFIRAOS. It includes the following files to
define all the parameters.

- \ref sim to define some general parameters. It deletes all previous entries,

- \ref dm_dual to define the dual deformable mirrors

- \ref atm_mk13n50p to define the turbulence profiles using the Mauna Kea 13N median seeing

- \ref wfs_lgs_ttf_tt to define the wavefront sensors which is composed of the following

    - \ref powfs_common
    - \ref powfs_none
    - \ref powfs_shlgs
    - \ref llt_CL
    - \ref powfs_shttf
    - \ref powfs_shtt

- \ref fov_sq34 which is composed of the following
    - \ref fit_sq34 to define the DM fitting fov
    - \ref evl_sq34 to define the science fov for performance evaluation.

- \ref recon to define the parameters about turbulence reconstruction and DM fitting

- \ref dbg to define debugging and saving/loading geometry and telemetry.


\page scao_base_example SCAO

The clasic AO configuration sans WFS is defined in \ref scao_base. It includes 

- \ref sim to define some general parameters. It deletes all previous entries,

- \ref atm_mk13n50p to define the turbulence profile

- \ref atmr_single to force single layer reconstruction

- \ref dm_single to define a single DM on the ground

- \ref fov_oa to define only on axis FoV

- \ref recon to define the parameters about turbulence reconstruction and DM fitting

- \ref dbg to define debugging and saving/loading geometry and telemetry.

\page scao_ngs_example SCAO NGS

The classic NGS AO is defined in \ref scao_ngs which supplements the basic \ref scao_base_example definition with NGS WFS \ref powfs_shngs. 

\page scao_lgs_example SCAO LGS

The classic LGS AO is defined in \ref scao_lgs which supplements the basic \ref scao_base_example definition with an LGS WFS \ref powfs_shlgs and a Tip/Tilt/Focus NGS WFS \ref powfs_shttf. 



