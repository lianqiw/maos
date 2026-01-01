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
#ifndef AOS_LIB_SLAVING_H
#define AOS_LIB_SLAVING_H
#include "../math/mathdef.h"
/**
   \file slaving.h
   Routines to compute interpolation matrices for stuck and floating actuators.
*/
dcell *genactcpl(const_anyarray HA, const dmat *W1);
dspcell *slaving(loccell *aloc, const dcell *actcpl, const lcell *actstuck, const lcell *actfloat, real thres, real scale, int mode);
void act_stuck(loccell *aloc, anyarray HA, const lcell *stuck);
void act_zero(loccell *aloc, const dcell *HB, const lcell *dead);
void act_float(loccell *aloc, dspcell **HA, const dcell *HB, const lcell *actfloat);
void act_stuck_cmd(loccell *aloc, const dcell *adm, const lcell *stuck);
dspcell* act_float_interp(loccell *aloc, const lcell *actfloat);
dsp* act_extrap_each(loc_t *aloc,const dmat *actcpl, real scl);
dspcell* act_extrap(loccell *aloc,const dcell *actcpl, real scl, int lor);
#endif
