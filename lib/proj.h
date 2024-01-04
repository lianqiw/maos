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



#ifndef AOS_LIB_PROJ_H
#define AOS_LIB_PROJ_H
/**
   \file proj.h
   Project from a tilt surface to flat surface (M3).
 */
void proj_rect_grid(real* phiout, const rmap_t *mapin, real thetax, real thetay,
		    const loc_t *locout,const real ratiox, const real ratioy,
		    const real *ampout,
		    real sc, real hs, real ht,
		    real betax, real betay);
void m3proj(const rmap_t *tsurf, dmat *opd, const loc_t *locin, real thetax, real thetay, real hs);
dmat *m3proj2(dmat *mapin_0, const char *keywords, const loc_t *locout, real thetax, real thetay, real hs);
#endif
