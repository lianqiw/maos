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



#ifndef __MAOS_H
#define __MAOS_H
#include "common.h"
#include "powfs.h"
#include "recon.h"
#include "aper.h"
#include "sim_utils.h"
#include "surf.h"
#include "pywfs.h"
#if USE_CUDA
#include "../cuda/gpu.h"
#endif
//also used by mex
void maos_setup(const parms_t *parms);
void maos_reset();
void maos_sim();
void maos_isim(int isim);
sim_t *maos_iseed(int iseed);
#endif

