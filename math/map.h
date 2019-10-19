/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_LIB_MAP_H
#define AOS_LIB_MAP_H
#include "mathdef.h"
#include "type.h"
/**
   \file map.h
   This file defines functions relates to map_t, etc.
*/

#define rmapfree(A) ({dfree_do((dmat*)A);A=NULL;})
#define mapfree(A) ({dfree_do((dmat*)A);A=NULL;})
map_t *mapnew(long nx, long ny, real dx, real dy);
map_t *mapnew2(map_t *A);
map_t *mapref(map_t *A);
void mapcircle(map_t *map, real r, real val);
void mapcircle_symbolic(map_t *map, real r);
void map_d_din(map_t *map, real *d, real *din);
void create_metapupil(map_t **map, long* nx, long* ny, dmat *dirs, real D, 
		      real ht, real dx, real dy, real offset,real guard, 
		      long ninx, long niny, int pad,int square);
dmat *mkcirmap(long nx, long ny, real cx, real cy, real r);

map_t* d2map(const dmat *in);
mapcell *dcell2map(const dcell *in);
rmap_t* d2rmap(const dmat *in);
rmap_t **dcell2rmap(int *nlayer, const dcell *in);

void map_header(map_t *map);
void rmap_header(rmap_t *map);

//map_t *mapreaddata(file_t *fp, header_t *header);
//rmap_t *rmapreaddata(file_t *fp, header_t *header);
//void mapwritedata(file_t *fp, map_t *map);


#endif
