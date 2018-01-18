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

#ifndef AOS_LIB_LOCBIN_H
#define AOS_LIB_LOCBIN_H
#include "loc.h"
/**
   \file locbin.h
   i/o functions for loc_t, map_t.
*/
loc_t *locreaddata(file_t *fp, header_t* header);
void locwritedata(file_t *fp, const loc_t *loc);
void locwrite(const loc_t *loc, const char *format,...) CHECK_ARG(2);
loc_t *locread(const char *format,...) CHECK_ARG(1);
void locarrwrite(loc_t ** loc, int nloc, const char *format,...) CHECK_ARG(3);
loc_t ** locarrread(int *nloc, const char *format,...) CHECK_ARG(2);

map_t *mapread(const char *format,...) CHECK_ARG(1);
map_t *mapreaddata(file_t *fp, header_t *header);
rmap_t *rmapread(const char *format,...) CHECK_ARG(1);
rmap_t **rmaparrread(int *nlayer, const char *format, ...);
void mapwritedata(file_t *fp, map_t *map);
void maparrwrite(map_t ** map, int nmap, const char *format,...) CHECK_ARG(3);
map_t **maparrread(int*nlayer, const char *format,...) CHECK_ARG(2);

map_t* d2map(dmat *in);
mapcell *dcell2map(dcell *in);
rmap_t* d2rmap(dmat *in);
rmap_t **dcell2rmap(int *nlayer, dcell *in);

#endif
