\page page33_example Examples
\tableofcontents
- \subpage nfiraos_example
- \subpage scao_base_example
- \subpage scao_lgs_example
- \subpage scao_ngs_example
- \subpage page34_conf

\section sect-example-list Available configurations

The baseline configurations in `config/maos` are complete and can be used with
`-c`:

- \ref default (includes \ref mcao_lgs): LGS MCAO, the default.
- \ref mcao_ngs -- NGS MCAO.
- \ref scao_lgs, \ref scao_ngs, \ref scao_pywfs -- single conjugate AO with LGS,
  Shack-Hartmann NGS, or pyramid NGS WFS.
- \ref ltao -- LGS LTAO with a single DM and on axis DM fitting.
- \ref glao -- GLAO with a single DM.

The folder `config/maos/examples` contains telescope apertures, which are
included by other files (\ref keck, \ref tmt, \ref gemini, \ref elt), and
configurations of existing or planned AO systems, for example:

- \ref keck_ngs -- Keck NGS AO.
- \ref keck_kapa -- Keck KAPA LGS LTAO.
- \ref strata and the `strata_*.conf` variants: Keck STRATA with the adaptive
  secondary mirror (KASM) for different instruments and fields of view.
- \ref nfiraos_lgs, \ref nfiraos_lgs_noimpl, \ref nfiraos_lgs_pyttf,
  \ref nfiraos_ngs -- TMT NFIRAOS.
- \ref gems_lgs -- a test setup of the Gemini South MCAO system GeMS.
- \ref raven_lgs, \ref raven_ngs -- MOAO configurations on the Gemini aperture
  modeled after Raven.
- \ref morfeo_lgs -- LGS MCAO on the ELT aperture.
- \ref mcao_2dm30, \ref mcao_2dm60, \ref mcao_3dm60 and the Keck variants
  \ref keck_2dm30, \ref keck_2dm60, \ref keck_3dm60 -- wide field MCAO studies.

For example, `maos -c nfiraos_lgs.conf -o nfiraos` runs TMT NFIRAOS in LGS mode.
\ref page34_conf shows the content of every file.


\page nfiraos_example LGS MCAO (default)

The default configuration \ref default uses \ref mcao_lgs, which is a generic
LGS MCAO system. It includes \ref mcao_base and \ref wfs_lgs_ttf_tt, which in
turn include the following files to define all the parameters.

- \ref sim to define some general parameters. It includes `keck.conf` (\ref keck)
  for the telescope aperture.

- \ref dm_dual to define the dual deformable mirrors (conjugated to 0 and 11.8 km)

- \ref atm_mk13n50p to define the turbulence profiles using the Maunakea 13N median seeing

- \ref fov_sq34 which is composed of the following
    - \ref evl_sq34 to define the science fov for performance evaluation.
    - \ref fit_sq34 to define the DM fitting fov

- \ref recon to define the parameters about turbulence reconstruction and DM fitting

- \ref dbg to define debugging and saving/loading geometry and telemetry.

- \ref wfs_lgs_ttf_tt to define the wavefront sensors which is composed of the following

    - \ref powfs_none to clear the WFS list. It includes \ref powfs_common for the common WFS parameters.
    - \ref powfs_shlgs for six Shack-Hartmann LGS WFS. It uses \ref llt_CL for the laser launch telescope (`powfs.fnllt`).
    - \ref powfs_shttf for one tip/tilt/focus NGS WFS (it includes \ref powfs_shtt).
    - \ref powfs_shtt for two tip/tilt NGS WFS (it includes \ref powfs_shngs).

The TMT NFIRAOS configuration \ref nfiraos_lgs is \ref mcao_lgs with the TMT
aperture (\ref tmt) and NFIRAOS specific parameters.


\page scao_base_example SCAO

The classic AO configuration without WFS is defined in \ref scao_base. It includes

- \ref sim to define some general parameters,

- \ref atm_mk13n50p to define the turbulence profile

- \ref atmr_single to force single layer reconstruction

- \ref dm_single to define a single DM on the ground

- \ref fov_oa to define only on axis FoV

- \ref recon to define the parameters about turbulence reconstruction and DM fitting

- \ref dbg to define debugging and saving/loading geometry and telemetry.

It also sets `fit.alg=0` and `sim.cachedm=0`.

\page scao_ngs_example SCAO NGS

The classic NGS AO is defined in \ref scao_ngs which supplements the basic \ref scao_base_example definition with an NGS WFS \ref powfs_shngs.

The pyramid WFS version \ref scao_pywfs uses \ref powfs_pyngs (with \ref pywfs)
instead, together with least squares reconstruction (`recon.alg=1`).

\page scao_lgs_example SCAO LGS

The classic LGS AO is defined in \ref scao_lgs which supplements the basic \ref scao_base_example definition with an LGS WFS \ref powfs_shlgs and a Tip/Tilt/Focus NGS WFS \ref powfs_shttf, both on axis.
