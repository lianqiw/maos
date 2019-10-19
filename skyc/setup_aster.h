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

#ifndef SKYC_SETUP_ASTER_H
#define SKYC_SETUP_ASTER_H
#include "parms.h"
#include "types.h"
ASTER_S *setup_aster_comb(int *naster, const STAR_S *star, int nstar, const PARMS_S *parms);
ASTER_S *setup_aster_sample(void);
void setup_aster(ASTER_S *aster, POWFS_S *powfs, const PARMS_S *parms, SIM_S *simu);
void free_aster(ASTER_S *aster, int naster, const PARMS_S *parms);
void setup_aster_copystar(ASTER_S *aster, STAR_S *star, const PARMS_S *parms);
void setup_aster_gm(ASTER_S *aster, STAR_S *star, const PARMS_S *parms);
void setup_aster_controller(SIM_S *simu, ASTER_S *aster, const PARMS_S *parms);
void setup_aster_read(ASTER_S *aster, const PARMS_S *parms, int seed);
int  setup_aster_select(real *result, ASTER_S *aster, int naster, STAR_S *star, real maxerror,
			const PARMS_S *parms);
void setup_aster_wvf(ASTER_S *aster, STAR_S *star, const PARMS_S *parms);
void setup_aster_ztilt(ASTER_S *aster, STAR_S *star, const PARMS_S *parms);
#endif
