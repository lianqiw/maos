#ifndef AOS_CUDA_TYPES_H
#define AOS_CUDA_TYPES_H
typedef struct{
    float *p;
    int nx;
    int ny;
    int ref;
}curmat;
/**
   curcell: call curcell is created by curcellnew2, p0 contains the allocated
   memory. p contains points inside the memory. curcell can therefore be casted to
   a curmat and do dot/copy operations.
*/
typedef struct{
    curmat **p;
    int nx;
    int ny;
}curcell;

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
typedef struct{
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
}cumap_t;


#endif
