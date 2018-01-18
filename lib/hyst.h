/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

typedef struct HYST_T HYST_T;
HYST_T *hyst_new(dmat *coeff, int naloc);
void hyst_reset(HYST_T *hyst);
void hyst_free(HYST_T *in);
void hyst_dmat(HYST_T *hyst, dmat *dmreal, const dmat *dmcmd);
void hyst_dcell(HYST_T **hyst, dcell *dmreal, const dcell *dmcmd);
void hyst_calib(HYST_T *hyst, int i);
#endif
