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


#ifndef AOS_LIB_MKW_H
#define AOS_LIB_MKW_H
#include "../math/mathdef.h"
#include "../math/mathdef.h"
/**
   \file mkw.h
   Contains functions that computes the bilinear weighting function W0, W1
*/
void mkw_amp(loc_t *loc,real *amp,dsp **W0,dmat **W1);
void mkw_circular(loc_t *loc,real cx, real cy,
		   real r,dsp **W0,dmat **W1);
void mkw_annular(loc_t *loc, real cx, real cy, real cri, real cro, 
		 dsp **W0, dmat **W1);
#endif
