/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
typedef struct GENATM_T{
    rand_t *rstat;   /**<The random stream*/
    double *wt;      /**<The layer weights*/
    double r0;       /**<The Fried Parameter*/
    double L0;       /**<The outer scale*/
    double dx;       /**<The sampling*/
    long nx;         /**<Number of pixels along x*/
    long ny;         /**<Number of pixels along y*/
    long nlayer;     /**<The number of layers*/
    long ninit;      /**<In Fractal method, the size of initial screen*/
    long share;      /**<Use file backend for sharing of atmosphere*/
    dmat *r0logpsds; /**<Spatial PSD of log(r0) (m)*/
    /*The following are private data. do not set when call. */
    mapcell *screen;  /**<The destination screen pointer*/
    dmat *spect;     /**<The turbulence spectrum, sqrt of PSD*/
    long method;     /**<The method*/
}GENATM_T;
map_t *genatm_simple(double r0, double L0, double dx, double nx);
dmat *genatm_loc(loc_t *loc, double r0, double L0);
mapcell* genatm_from_spect(GENATM_T *data);
mapcell* vonkarman_screen(GENATM_T *data);
mapcell* biharmonic_screen(GENATM_T *data);
mapcell *fractal_screen(GENATM_T *data);
dmat* turbcov(dmat *r, double rmax, double r0, double L0);
dmat *spatial_psd(long nx, long ny, double dx, double strength, 
		  double fmin, double fmax, double slope, double power);
INLINE double r02strength(double r0){
    return 0.0229*pow(r0,-5./3.)*pow((0.5e-6)/(2.*M_PI),2);
}
INLINE dmat* turbpsd(long nx, long ny, double dx, double r0, double L0, double slope, double power){
    return spatial_psd(nx, ny, dx, r02strength(r0), 1./L0, INFINITY, slope, power);
}

double calc_aniso(double r0, int nht, double *ht, double *wt);
double calc_greenwood(double r0, int nps, double *ws, double *wt);
double calc_aniso2(double r0, int nht, double *ht, double *wt, double hc1, double hc2);
#endif
