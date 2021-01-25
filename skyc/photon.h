/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef SKYC_PHOTON_H
#define SKYC_PHOTON_H
typedef struct ZB_S{
    real ZJ, ZH, ZK;   /**<zero magnitude flux*/
    real BJ, BH, BK;   /**<sky background*/
}ZB_S;
void photon_flux(const ZB_S *zb, real *Np, real *Nptot, real *Nbtot, real *QCSNR, real *QCNEA,
		 int nwvl, real* wvls, real *mags, 
		 real dxsa, int iscircle, real pixtheta, 
		 real dt, real za, real *strehl, real imperrnm,
		 real *thruput, real *qe, real rne);
#endif
