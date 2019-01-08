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
   \page SCAO_NGS Example: SCAO NGS

   The following shows the content of scao_ngs.conf, which is the baseline
   configuration of a single conjugate NGS AO system on TMT.
   
   \include config/maos/scao_ngs.conf

   It includes "scao_base.conf" that defines a basic SCAO system on TMT
   \include config/maos/scao_base.conf

   and "wfs_single_ngs.conf" that defines parameters for a single high order NGS WFS.
*/

/**
   \page SCAO_LGS Example: SCAO LGS

   The following shows the content of scao_ngs.conf, which is the baseline
   configuration of a single conjugate LGS AO system on TMT.
   
   \include config/maos/scao_lgs.conf

   It includes "scao_base.conf" that defines a basic SCAO system on TMT
   \include config/maos/scao_base.conf
   
   and "wfs_single_lgs.conf" that defines parameters for a single high order NGS WFS.
*/

/**
   \page SCAO_BASE Example: SCAO BASE 

   The following shows the content of scao_base.conf, which defines a basic SCAO
   system on TMT without information about wfs.

   \include config/maos/scao_base.conf
   
   It includes the following files to define all the parameters

   \c "sim.conf" for general parameters

   \include config/maos/sim.conf

   \c "aper_tmt.conf" for the aperture of TMT

   \include config/maos/aper_tmt.conf

   \c "atm_scao.conf" for the atmosphere.

   \include config/maos/atm_scao.conf

   \c "dm_single.conf" for a single deformable mirror

   \include config/maos/dm_single.conf

   \c "fov_oa.conf" for a single on axis science and fit field.

   \include config/maos/fov_oa.conf

   \c "recon.conf" for tomography and DM fit.

   \include config/maos/recon.conf

   \c "dbg.conf" for debugging, gemetry and telemetry data

   \include config/maos/dbg.conf

   Finally it overrides the tomography and DM fit algorithm to cholesky
backsolve and disables cache of DM commands.

 */
