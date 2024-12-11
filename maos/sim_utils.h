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



#include "common.h"
#ifndef AOS_SIM_UTILS_H
#define AOS_SIM_UTILS_H
void atm2xloc(dcell **opdx, const sim_t *simu);
void sim_update_etf(sim_t *simu);
void update_wfsflags(sim_t *simu);
void seeding(sim_t *simu);
sim_t* init_simu(const parms_t *parms,powfs_t *powfs, aper_t *aper,recon_t *recon,int iseed);
void free_simu(sim_t *simu);
void print_progress(sim_t *simu);
void save_skyc(powfs_t *powfs, recon_t *recon, const parms_t *parms);
void blend_screen_side(map_t *atm1, map_t *atm2, long overx, long overy);
void genatm(sim_t *simu);
#endif
