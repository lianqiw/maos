/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


#ifndef AOS_LIB_TURBULENCE_H
#define AOS_LIB_TURBULENCE_H
/**
   \file turbulence.h
   Contains routines to generate atmospheric turbulence screens
*/
#include "../math/mathdef.h"
typedef struct genatm_t{
    rand_t *rstat;   /**<The random stream*/
    real *wt;      /**<The layer weights*/
    real r0;       /**<The Fried Parameter*/
    real *L0;       /**<The outer scale*/
    real dx;       /**<The sampling*/
    real fmin;     /**<Minimum spatial frequency*/
    real fmax;     /**<Maximum spatial frequency*/
    real slope;    /**<Power slope of PSD. -11/3 for Von Karman, -4 for Biharmonic. -1 for fractal. */
    long nx;         /**<Number of pixels along x*/
    long ny;         /**<Number of pixels along y*/
    long nlayer;     /**<The number of layers*/
    long ninit;      /**<In Fractal method, the size of initial screen*/
    long share;      /**<Use file backend for sharing of atmosphere*/
    dmat *r0logpsds; /**<Spatial PSD of log(r0) (m)=beta*f^alpha. [alpha, beta, minfreq, maxfreq]*/
    /*The following are private data. do not set when call. */
    mapcell *screen;  /**<The destination screen pointer*/
}genatm_t;
map_t *genatm_simple(real r0, real L0, real dx, long nx);
dmat *genatm_loc(loc_t *loc, real r0, real L0);
//mapcell* genatm_from_spect(genatm_t *data);
mapcell* genscreen(genatm_t *data);
mapcell* genscreen_str(const char *header);
dmat* turbcov(dmat *r, real rmax, real r0, real L0);
void spatial_psd(dmat **out, long nx, long ny, real dx, real strength, real L0,
		 real fmin, real fmax, real slope, real power);
dmat* turbpsd(long nx, long ny, real dx, real r0, real L0, real slope, real power);

real calc_aniso(real r0, int nht, real *ht, real *wt);
real calc_greenwood(real r0, int nps, real *ws, real *wt);
real calc_aniso2(real r0, int nht, real *ht, real *wt, real hc1, real hc2);
#endif
