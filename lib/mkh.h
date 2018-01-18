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

#ifndef AOS_LIB_MKH_H
#define AOS_LIB_MKH_H
#include "../math/mathdef.h"
#include "../math/mathdef.h"
/**
   \file mkh.h
   Contains functions that create ray tracing operator
*/
dsp *mkhb_cubic(loc_t *locin, loc_t *locout, double displacex, double displacey, double scale, double iac);
dsp *mkh_cubic(loc_t *locin, loc_t *locout, double displacex, double displacey, double scale, double iac);
dsp * mkhb(loc_t *locin, loc_t *locout, double displacex, double displacey,double scale);
dsp * mkh(loc_t *locin, loc_t *locout, double displacex, double displacey,double scale);
dsp *mkhbin1d(dmat *xin, dmat *xout);		   
#endif
