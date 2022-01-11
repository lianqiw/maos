/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file cure.h
*/
#include "../math/mathdef.h"
void cure1d(dmat **phix, const dmat*gx, const dmat *gy, real dx);
void cure(dmat **phix, const dmat*gx, const dmat *gy, real dx);
void cure_loc(dmat **phix, const dmat*grad, const loc_t *gloc);
void cure_loc2(dmat **phix, const dmat*grad, const loc_t *saloc);
