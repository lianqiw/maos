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
#ifndef AOS_CUDA_TYPES_H
#define AOS_CUDA_TYPES_H
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include <cuComplex.h>
#define fcomplex cuFloatComplex
#define dcomplex cuDoubleComplex

typedef struct{
    float *p;
    int nx;
    int ny;
    int *nref;
}curmat;

typedef struct{
    curmat **p;
    int nx;
    int ny;
}curcell;
typedef struct{
    fcomplex *p;
    int nx;
    int ny;
    int *nref;
}cucmat;

typedef struct{
    cucmat **p;
    int nx;
    int ny;
}cuccell;

typedef struct{
    int *p;
    int *i;
    float *x;
    int nx;
    int ny;
    int nzmax;
}cusp;
typedef struct{
    cusp **p;
    int nx;
    int ny;
}cuspcell;
typedef struct{
    cuspcell *Mt;
    curcell *U;
    curcell *V;
}cumuv_t;

typedef struct{
    float (*loc)[2];/*in device. */
    float dx;
    int nloc;
}culoc_t;
/*
  We use a single map_t to contain all layers instead of using an array of map_t
  because we want to use layered texture. This preference can be retired since
  the speed is largely the same with layered texture or flat memory.
 */
struct cumap_t{
    cudaArray *ca;/*3D array. for layered texture */
    float **p;/*float array. */
    float *ht;
    float *vx;
    float *vy;
    float *ox;
    float *oy;
    float *dx;
    float *iac;
    int* cubic;
    int* nx;
    int* ny;
    int nlayer;
    float **cc;/*coefficients for cubic dm. */
};


#endif
