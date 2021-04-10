/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file genseotf.h

   Contains routines to generate mean short exposure (tip/tilt
   removed) pixel intensities. Mostly used for LGS pixel intensity for its
   matched filter. Structure functions from kolmogorov spectrum is used. Not
   able to take into account outerscale yet.

   \todo find ways to factor in outerscale effect (use von karman spectrum
   instead of kolmogorov) 
*/
#ifndef AOS_GENSEOTF_H
#define AOS_GENSEOTF_H
#include "common.h"
void gensei(const parms_t *parms, powfs_t *powfs, int ipowfs);
void genseotf(const parms_t *parms, powfs_t *powfs, int ipowfs);
void genselotf(const parms_t *parms, powfs_t *powfs, int ipowfs);
void gensepsf(const parms_t *parms, powfs_t *powfs, int ipowfs);
void genmtch(const parms_t *parms, powfs_t *powfs, const int ipowfs);
#endif
