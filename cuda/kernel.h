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
#ifndef AOS_CUDA_KERNEL_H
#define AOS_CUDA_KERNEL_H
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "common.h"
#include <cuComplex.h>
#define fcomplex cuFloatComplex
#define dcomplex cuDoubleComplex
const int NG1D=64;
const int NG2D=8;
const int WRAP_SIZE=32; /*The wrap size is currently always 32 */
const int REDUCE_WRAP=8;
const int REDUCE_WRAP_LOG2=3;
const int DIM_REDUCE=WRAP_SIZE*REDUCE_WRAP; /*dimension to use in reduction. */
const int REDUCE_STRIDE=WRAP_SIZE+WRAP_SIZE/2+1;
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ <200
__device__ inline float atomicAdd(float* address, float val)
{
    float old = *address;
    float assumed;
    do {
	assumed = old;
	old = __int_as_float( atomicCAS((unsigned int*)address,
					__float_as_int(assumed),
					__float_as_int(val + assumed)));
    } while (assumed != old);
    return old;
}
#endif
__device__ inline float atomicMax(float* address, float val)
{
    float old = *address;
    float assumed;
    do {
	assumed = old;
	old = __int_as_float( atomicCAS((unsigned int*)address,
					__float_as_int(assumed),
					__float_as_int(val>assumed?val:assumed)));
    } while (assumed != old);
    return old;
}

__device__ inline float CABS2(fcomplex r){
    const float a=cuCrealf(r);
    const float b=cuCimagf(r);
    return a*a+b*b;
}

__global__ void set_do(float *a, float alpha, int n);
__global__ void scale_do(float *restrict in, int n, float alpha);
__global__ void add_ptt_do(float *restrict opd, float (*restrict loc)[2], int n, float pis, float tx, float ty);
__global__ void add_focus_do(float *restrict opd, float (*restrict loc)[2], int n, float focus);
__global__ void add_ngsmod_do(float *restrict opd, float (*restrict loc)[2], int n, 
			      float m0, float m1, float m2, float m3, float m4, float focus,
			      float thetax, float thetay, float scale, float ht, float alpha);

__global__ void add_do(float *vec, float beta, int n);
__global__ void addcabs2_do(float *restrict a, float alpha, const fcomplex *restrict b, float beta, int n);
__global__ void add_do(float *restrict a, float *alpha1, float alpha2, const float *restrict b,  int n);
__global__ void add_do(float *restrict a, const float *restrict b, float *beta1, float beta2, int n);
__global__ void add_do(float *restrict a, float *alpha1, float alpha2, 
		       const float *restrict b, float *beta1, float beta2, int n);

__global__ void max_do(float *restrict res, const float *a, const int n);
__global__ void sum_do(float *restrict res, const float *a, const int n);
__global__ void sum2_do(float *restrict res, const float *a, const int n);
__global__ void inn_do(float *res_add, const float *a, const float *b, const int n);

inline void inn_wrap(float *res_add, const float *a, const float *b, const int n, cudaStream_t stream){
    inn_do<<<DIM(n, DIM_REDUCE), DIM_REDUCE*sizeof(float), stream>>>(res_add, a, b, n);
}
inline static void sum_wrap(float *res, const float * a, const int n, cudaStream_t stream){
    sum_do<<<DIM(n, DIM_REDUCE), DIM_REDUCE*sizeof(float), stream>>>(res,a,n);
}
inline static void sum2_wrap(float *res, const float * a, const int n, cudaStream_t stream){
    sum2_do<<<DIM(n, DIM_REDUCE), 0, stream>>>(res,a,n);
}
inline static void max_wrap(float *res, const float * a, const int n, cudaStream_t stream){
    max_do<<<DIM(n, DIM_REDUCE), DIM_REDUCE*sizeof(float), stream>>>(res,a,n);
}
__global__ void embed_do(fcomplex *out, float *in, int nx);
__global__ void extract_do(float *out, fcomplex *in, int nx);
__global__ void perm_f_do(fcomplex *restrict out, const fcomplex *restrict in, int *restrict perm, int nx);
__global__ void perm_i_do(fcomplex *restrict out, const fcomplex *restrict in, int *restrict perm, int nx);
__global__ void perm_f_do(float *restrict out, const float *restrict in, int *restrict perm, int nx);
__global__ void perm_i_do(float *restrict out, const float *restrict in, int *restrict perm, int nx);
__global__ void embed_wvf_do(fcomplex *restrict wvf, 
			     const float *restrict opd, const float *restrict amp, 
			     const int *embed, const int nloc, const float wvl);
__global__ void corner2center_do(fcomplex *restrict out, int noutx, int nouty,
				 const fcomplex *restrict in, int ninx, int niny);
__global__ void corner2center_abs2_do(float *restrict out, int noutx, int nouty,
				      const fcomplex *restrict in, int ninx, int niny);
__global__ void corner2center_abs2_atomic_do(float *restrict out, int noutx, int nouty,
					     const fcomplex *restrict in, int ninx, int niny);
__global__ void fftshift_do(fcomplex *wvf, const int nx, const int ny);
__global__ void add_tilt_do(float *opd, int nx, int ny, float ox, float oy, float dx, float ttx, float tty);
__global__ void cwm_do(fcomplex *dest, float *from, int n);
__global__ void unwrap_phase_do(fcomplex *wvf, float *opd, int *embed, int n, float wvl);
#endif
