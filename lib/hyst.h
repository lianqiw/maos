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
   \file hyst.h
   Routines for hysteresis modeling
*/

#ifndef SKYC_HYST_H
#define SKYC_HYST_H

#include "../math/mathdef.h"

typedef struct hyst_t hyst_t;
hyst_t *hyst_new(real hysteresis, real alpha, real stroke, long nact);
void hyst_reset(hyst_t *hyst);
void hyst_free(hyst_t *in);
void hyst_dmat(hyst_t *hyst, dmat *dmreal, const dmat *dmcmd);
void hyst_dcell(hyst_t **hyst, dcell *dmreal, const dcell *dmcmd);
dmat* hyst_test(real hysteresis, real hyst_alpha, real hyst_stroke, const dmat *dmcmd);
#endif
