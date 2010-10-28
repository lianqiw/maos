/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_LOCBIN_H
#define AOS_LOCBIN_H
#include "bin.h"
#include "loc.h"
void locwrite(const LOC_T *loc, const char *format,...) CHECK_ARG(2);
LOC_T *locread(const char *format,...) CHECK_ARG(1);
void locarrwrite(LOC_T ** loc, int nloc, const char *format,...) CHECK_ARG(3);
LOC_T ** locarrread(int *nloc, const char *format,...) CHECK_ARG(2);

MAP_T *sqmapread(const char *format,...) CHECK_ARG(1);
RECTMAP_T *rectmapread(const char *format,...) CHECK_ARG(1);
void sqmapwrite(const MAP_T *map, const char *format,...) CHECK_ARG(2);

void sqmaparrwrite(MAP_T ** map, int nmap, const char *format,...) CHECK_ARG(3);
MAP_T **sqmaparrread(int*nlayer, const char *format,...) CHECK_ARG(2);
#endif
