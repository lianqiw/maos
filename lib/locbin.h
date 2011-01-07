/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
loc_t *locreaddata(file_t *fp, uint32_t magic, char *header0);
void locwritedata(file_t *fp, const loc_t *loc);
void locwrite(const loc_t *loc, const char *format,...) CHECK_ARG(2);
loc_t *locread(const char *format,...) CHECK_ARG(1);
void locarrwrite(loc_t ** loc, int nloc, const char *format,...) CHECK_ARG(3);
loc_t ** locarrread(int *nloc, const char *format,...) CHECK_ARG(2);

map_t *sqmapread(const char *format,...) CHECK_ARG(1);
map_t *sqmapreaddata(file_t *fp, uint32_t magic, char *header);
rectmap_t *rectmapread(const char *format,...) CHECK_ARG(1);
rectmap_t *rectmapreaddata(file_t *fp);
void sqmapwrite(const map_t *map, const char *format,...) CHECK_ARG(2);

void sqmaparrwrite(map_t ** map, int nmap, const char *format,...) CHECK_ARG(3);
map_t **sqmaparrread(int*nlayer, const char *format,...) CHECK_ARG(2);
#endif
