/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_GENSEOTF_H
#define AOS_GENSEOTF_H
#include "common.h"
void gensei(const PARMS_T *parms, POWFS_T *powfs, int ipowfs);
void genseotf(const PARMS_T *parms, POWFS_T *powfs, int ipowfs);
void genselotf(const PARMS_T *parms, POWFS_T *powfs, int ipowfs);
void gensepsf(const PARMS_T *parms, POWFS_T *powfs, int ipowfs);
void genmtch(const PARMS_T *parms, POWFS_T *powfs, const int ipowfs);
#endif
