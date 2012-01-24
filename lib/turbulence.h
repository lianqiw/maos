/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_TURBULENCE_H
#define AOS_TURBULENCE_H

#include "loc.h"
#include "random.h"
#include "type.h"
typedef struct GENSCREEN_T{
    rand_t *rstat;   /**<The random stream*/
    double *wt;      /**<The layer weights*/
    double r0;       /**<The Fried Parameter*/
    double l0;       /**<The outerscale*/
    double dx;       /**<The sampling*/
    long nx;         /**<Number of pixels along x*/
    long ny;         /**<Number of pixels along y*/
    long nlayer;     /**<The number of layers*/
    long ninit;      /**<In Fractal method, the size of initial screen*/
    long share;      /**<Use file backend for sharing of atmosphere*/
    long nthread;    /**<Number of threads to use*/
    /*The following are private data. do not set when call. */
    map_t **screen;  /**<The destination screen pointer*/
    dmat *spect;     /**<The turbulence spectrum, sqrt of PSD*/
    long ilayer;     /**<Current layer*/
    long method;     /**<The method*/
    pthread_mutex_t mutex_ilayer;/**<Mutex lock*/
}GENSCREEN_T;

extern int genscreen_keep_unused;
map_t** genscreen_from_spect(GENSCREEN_T *data);
map_t** vonkarman_screen(GENSCREEN_T *data);
map_t** biharmonic_screen(GENSCREEN_T *data);
map_t **fractal_screen(GENSCREEN_T *data);
dmat* turbcov(dmat *r, double rmax, double r0, double L0);
dmat *turbpsd_full(long nx, long ny, double dx, double r0, double L0, double slope, double power);
#define turbpsd(nx, ny, dx, r0, L0, power) turbpsd_full(nx, ny, dx, r0, L0, -11./6., power);

#endif
