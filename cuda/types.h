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
    int igpu;//which gpu we are on.
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
    int ref;
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
    float (*loc)[2];//in device.
    float dx;
    int nloc;
}culoc_t;
/*
  We use a single map_t to contain all layers instead of using an array of map_t
  because we want to use layered texture. This preference can be retired since
  the speed is largely the same with layered texture or flat memory.
 */
struct cumap_t{
    cudaArray *ca;//3D array. for layered texture
    float **p;//float array.
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
};


#endif
