/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_LIB_MKG_H
#define AOS_LIB_MKG_H
#include "../math/mathdef.h"
/**
   \file mkg.h
   Contains function that creates average gradient operator
*/
dsp *mkg(loc_t* xloc, loc_t *ploc, dmat *amp, loc_t *saloc, dmat *saa, real saat,
	real scale, real dispx, real dispy, int do_partial);
dsp *mkgt(loc_t *xloc, loc_t *ploc, dmat *amp, loc_t *saloc, dmat *saa, real saat,
	real scale, real dispx, real dispy, int do_partial);
dmat* mkg_ngs(dmat *sacent, real thetax, real thetay, real hc, real hs, 
	int indps, int indastig, int indfocus, int ahstfocus);
#endif
