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
#ifndef AOS_LIB_MAP_H
#define AOS_LIB_MAP_H
#include "../math/mathdef.h"
/**
   \file map.h
   This file defines functions relates to map_t, etc.
*/
/**
   OPD or Amplitude map defined on square/rectangular grids. with equal spacing
   on x/y. Can be casted to dmat.
*/
struct map_t{
    /*The OPD, takes the same form of dmat so can be casted. */
    union{
        MATARR(real,);
		dmat dmat[1]; //base[1] is conveniently used to return a pointer as cell without casting which is error prone
    };
    real ox;      /**<Origin in x*/
    real oy;      /**<Origin in y*/
    real dx;      /**<Sampling along x*/
    real dy;      /**<Sampling along y*/
    real ht;       /**<Heigh conjugation of this surface*/
    real vx;      /**Wind velocity. Useful for atmospheric grid*/
    real vy;      /**Wind velocity. Useful for atmospheric grid*/
    real iac;     /**<Inter-actuator coupling. >0: use cubic influence function*/
	real dratio;  /**<M1 to DM demagnification ratio (>1 for shrinking beam)*/
	real aoi;	  /**<Angle of incidence in radian for ray tracing input */
};

/**
   Map with different x/y sampling. Can be cased to dmat
*/
struct rmap_t{
    union{
        MATARR(real,);
        dmat dmat[1];
    };
    real ox;      /**<Origin in x*/
    real oy;      /**<Origin in y*/
    real dx;      /**<Sampling along x (first dimension)*/
    real dy;      /**<Sampling along y (second dimension)*/
    real txdeg;   /**<the x tilt angle in degree wrt beam (90 is prep), */
    real tydeg;   /**<the y tilt angle in degree wrt beam (90 is prep), */
    real ftel;    /**<Effective focal length of the telescope*/
    real fexit;   /**<The distance between the exit pupil and the focus*/
    real fsurf;   /**<The distance between the tilted surface (M3) and the focus*/
};

//map_t and rmap_t can be used as dmat
TYPE_ASSERT(offsetof(map_t, dmat)==0 && offsetof(map_t, ox)>=sizeof(dmat));
TYPE_ASSERT(offsetof(rmap_t, dmat)==0 && offsetof(rmap_t, ox)>=sizeof(dmat));

map_t* map_convert(dmat* A);
cell* map_convert_cell(cell* A);
rmap_t* rmap_convert(dmat* A);
#define mapwrite(out, A...) writebin(dmat_cast(out), A)
#define mapread(A...)    map_convert(dmat_cast(readbin_id(M_MAP, 0, A)))

#define mapcellread(A...) (mapcell*)cellconvert(readbin_id(M_MAP, 1, A), map_convert_cell)
#define mapcellread_mmap(A...) (mapcell*)cellconvert((cell*)dcellread_mmap(A), map_convert_cell)
#define mapcellnew (mapcell*)cellnew
#define mapccellnew (mapccell*)cellnew

#define rmapread(A...)    rmap_convert(dmat_cast(readbin_id(M_RMAP, 0, A)))
#define rmapwrite(out, A...) writebin(dmat_cast(out), A)
#define rmapcellnew  (rmapcell*)cellnew
#define rmapccellnew (rmapccell*)cellnew

#define rmapfree(A) cellfree(A) //if(A){dfree_do((dmat*)A);A=NULL;}
#define mapfree(A) cellfree(A) // if(A){dfree_do((dmat*)A);A=NULL;}
map_t *mapnew(long nx, long ny, real dx, real dy);
map_t *mapnew2(map_t *A);
map_t *mapref(const map_t *A);
map_t *mapdup(const map_t *A);
void mapcircle(map_t *map, real r, real val);
void map_d_din(real *d, real *din, map_t *map);
void create_metapupil(map_t **map, long* nx, long* ny, dmat *dirs, real D, 
		      real ht, real dx, real dy, real offset,real guard, 
		      long ninx, long niny, int pad,int square);

void map_make_keywords(cell* map);
void map_parse_keywords(map_t* in);
void rmap_make_keywords(cell* map);
void rmap_parse_keywords(rmap_t* in);
//map_t *mapreaddata(file_t *fp, header_t *header);
//rmap_t *rmapreaddata(file_t *fp, header_t *header);
//void mapwritedata(file_t *fp, map_t *map);
void map_blend(map_t *atm1, map_t *atm2, long overx, long overy);

#endif
