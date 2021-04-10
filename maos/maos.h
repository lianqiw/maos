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
/**
   \file maos.h
   
   Contains main() and the entry into simulation maos(). The main() is separated
   into main.c so that maos.c can be a part of libaos.la which is callable by
   MATLAB.
*/
#ifndef __MAOS_H
#define __MAOS_H
#include "common.h"
#include "setup_powfs.h"
#include "setup_recon.h"
#include "setup_aper.h"
#include "sim_utils.h"
#include "setup_surf.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
void maos_setup(const parms_t *parms);
void maos_reset();
void maos_sim();
void maos_isim(int isim);
sim_t *maos_iseed(int iseed);
#endif

