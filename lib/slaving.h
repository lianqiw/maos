/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_LIB_SLAVING_H
#define AOS_LIB_SLAVING_H
#include "type.h"
#include "imat.h"
spcell *slaving(loc_t **aloc, spcell *HA, dmat *W1, dcell *NW, icell *actstuck, icell *actfloat, double thres, double scl);
void act_stuck(loc_t **aloc, spcell *HA, dcell *HB, icell *stuck);
void act_zero(loc_t **aloc, dcell *HB, icell *dead);
void act_float(loc_t **aloc, spcell **HA, dcell *HB, icell *actfloat);
void act_stuck_cmd(loc_t **aloc, dcell *adm, icell *stuck);
spcell* act_float_interp(loc_t **aloc, icell *actfloat);
#endif
