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
#ifndef AOS_LIB_LLT_H
#define AOS_LIB_LLT_H
#include "loc.h"
/**
   \file llt.h
   Pupil of the laser launch telescope (LLT) of a laser guide star WFS: the
   points, the Gaussian beam amplitude, and the uplink focus and tip/tilt/focus
   removal of its OPD.
*/
pts_t* llt_pts(real lltd, real dsa, real dx);
dmat* llt_amp(const pts_t* pts, real lltd, real widthp, real misregx, real misregy);
void llt_add_focus(dcell* ncpa, const loc_t* loc, const dmat* amp, real lltd, real focus);
void llt_remove_ttf(dcell* ncpa, const loc_t* loc, const dmat* amp, real lltd, int ttfr, dmat** ppttf, dmat** pproj);
#endif
