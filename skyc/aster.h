/*
  Copyright 2009-2026 Lianqi Wang
  
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
aster_s *setup_aster_comb(int *naster, const star_s *star, int nstar, const parms_s *parms);
aster_s *setup_aster_sample(void);
void setup_aster(aster_s *aster, powfs_s *powfs, const parms_s *parms, sim_s *simu);
void free_aster(aster_s *aster, int naster, const parms_s *parms);
void setup_aster_controller(sim_s *simu, aster_s *aster, star_s *star, const parms_s *parms);
void setup_aster_read(aster_s *aster, const parms_s *parms, int seed);
int  setup_aster_select(real *result, aster_s *aster, int naster, star_s *star, real maxerror, const parms_s *parms);
void setup_aster_wvf(aster_s *aster, star_s *star);
void setup_aster_ztilt(aster_s *aster, star_s *star);
#endif
