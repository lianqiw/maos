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
/**
   \file twfs.h

   Truth WFS related
*/
#ifndef AOS_TWFS_H
#define AOS_TWFS_H
#include "common.h"
void twfs_setup_GR(recon_t* recon, const parms_t* parms);
void twfs_setup_RR(recon_t* recon, const parms_t* parms);
void twfs_recon(sim_t* simu);
#endif