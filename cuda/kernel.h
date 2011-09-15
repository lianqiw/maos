#ifndef AOS_CUDA_KERNEL_H
#define AOS_CUDA_KERNEL_H
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
__global__ void set_do(float *a, float alpha, int n);
__global__ void show_do(float *a, int nx, int ny);
__global__ void add_ptt_do(float *restrict opd, float (*restrict loc)[2], int n, float pis, float tx, float ty);
__global__ void add_ngsmod_do(float *restrict opd, float (*restrict loc)[2], int n, 
			      float m0, float m1, float m2, float m3, float m4,
			      float thetax, float thetay, float scale, float ht, float MCC_fcp, float alpha );
#endif
