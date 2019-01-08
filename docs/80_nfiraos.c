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
   \page nfiraos Example: NFIRAOS
   The following shows the content of nfiraos.conf, which is the baseline
    configuration of TMT NFIRAOS system
    
    \include config/maos/nfiraos.conf
    
    It includes the following files to define all the parameters.
    \c "sim.conf" to define some general parameters:
    \include config/maos/sim.conf
    
    \c "dm_dual.conf" to define the dual deformable mirrors
    \include config/maos/dm_dual.conf

    \c "atm_mk13n50p.conf" to define the turbulence profiles using the Mauna Kea 13N median seeing
    \include config/maos/atm_mk13n50p.conf

    \c "wfs_lgs_ttf_tt.conf" to define the wavefront sensors.
    \include config/maos/wfs_common.conf
    \include config/maos/wfs_lgs_ttf_tt.conf
    \include config/maos/llt_CL.conf

    \c "fov_iris.conf" to define the science fov for performance evaluation with IRIS.
    \include config/maos/fov_iris.conf

    \c "recon.conf" to define the parameters about turbulence reconstruction and DM fitting
    \include config/maos/recon.conf

    \c "dbg.conf" to define debugging and saving/loading geometry and telemetry.
    \include config/maos/dbg.conf

*/
