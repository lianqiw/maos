/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


#ifndef AOS_POWFS_H
#define AOS_POWFS_H
#include "common.h"
powfs_t * setup_powfs_init(const parms_t *parms, aper_t *aper);
void setup_powfs_misreg_tel(powfs_t *powfs, const parms_t *parms, aper_t *aper, int ipowfs);
void setup_powfs_misreg_dm(powfs_t *powfs, const parms_t *parms, aper_t *aper, int ipowfs);
void setup_shwfs_phy(const parms_t *parms,  powfs_t *powfs);
void setup_powfs_neasim(const parms_t *parms,  powfs_t *powfs);
void setup_powfs_calib(const parms_t *parms, powfs_t *powfs);
void free_powfs_unused(const parms_t *parms, powfs_t *powfs);
void free_powfs(const parms_t *parms, powfs_t *powfs);
//void test_powfs(const parms_t *parms, powfs_t *powfs);
void setup_shwfs_etf(powfs_t *powfs, const parms_t *parms, int ipowfs, int mode, int icol, 
  real deltah, real thresh);
void wfspupmask(const parms_t *parms, loc_t *loc, dmat *amp, int iwfs);

#endif
