/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_CUDA_UTILS_H
#define AOS_CUDA_UTILS_H
#include "types.h"

/*void gpu_set(int igpu);
  int  gpu_next(void);*/
void dbl2flt(float * restrict *dest, double *src, int n);
void spint2int(int * restrict *dest, spint *src, int n);

void cp2gpu(cumap_t **dest, map_t **source, int nps);
void cp2gpu(cusp **dest, const dsp *src, int tocsr);
void cp2gpu(cusp **dest, const spcell *src, int tocsr);
void cp2gpu(cuspcell **dest, const spcell *src, int tocsr);
void cp2gpu(float (* restrict *dest)[2], const loc_t *src);
void cp2gpu(float * restrict *dest, const double *src, int n);
void cp2gpu(float * restrict *dest, const dmat *src);
void cp2gpu(fcomplex * restrict *dest, const dcomplex *src, int n);
void cp2gpu(fcomplex * restrict *dest, const cmat *src);

void cp2gpu(curmat *restrict *dest, const dmat *src);
void cp2gpu(curcell *restrict *dest, const dcell *src);
void cp2gpu(cucmat *restrict *dest, const cmat *src);
void cp2gpu(cuccell *restrict *dest, const ccell *src);

void cp2gpu(int * restrict *dest, const long *src, int n);
void cp2gpu(int * restrict *dest, const spint *src, int n);
void cp2gpu(int * restrict *dest, const int *src, int n);
void cp2gpu(curmat *restrict *dest, const float *src, int nx, int ny, cudaStream_t stream);
inline void cp2gpu(curmat *restrict *dest, const float *src, int nx, cudaStream_t stream){
    cp2gpu(dest, src, nx, 1, stream);
}
inline void cp2gpu(curmat *restrict *dest, const smat *src, cudaStream_t stream=0){
    if(src){
	cp2gpu(dest, src->p, src->nx, src->ny, stream);
    }
}

void cuspmul (float *y, cusp *A, const float *x, int ncol, char trans,
	      float alpha, cusparseHandle_t handle);

void gpu_write(const float *p, int nx, int ny, const char *format, ...);
void gpu_write(const fcomplex *p, int nx, int ny, const char *format, ...);
void gpu_write(const int *p, int nx, int ny, const char *format, ...);
W01_T *gpu_get_W01(dsp *R_W0, dmat *R_W1);
void cp2cpu(double * restrict *dest, double alpha, float *src, double beta, int n, cudaStream_t stream, pthread_mutex_t* mutex=0);
void cp2cpu(dmat **out, double alpha, const curmat *in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void cp2cpu(smat **out, const curmat *in, cudaStream_t stream);
void cp2cpu(zmat **out, const cucmat *in, cudaStream_t stream);
void cp2cpu(dcell **out, double alpha, const curcell *in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void cp2cpu(scell **out, const curcell *in, cudaStream_t stream);
void cp2cpu(zcell **out, const cuccell *in, cudaStream_t stream);
inline void cp2cpu(dmat **out, const curmat *in, cudaStream_t stream){
    cp2cpu(out, 0, in, 1, stream);
}
inline void cp2cpu(dcell **out, const curcell *in, cudaStream_t stream){
    cp2cpu(out, 0, in, 1, stream);
}
inline void add2cpu(dmat **out, const curmat *in, cudaStream_t stream, pthread_mutex_t* mutex=0){
    cp2cpu(out, 1, in, 1, stream, mutex);
}
inline void add2cpu(dcell **out, const curcell *in, cudaStream_t stream, pthread_mutex_t* mutex=0){
    cp2cpu(out, 1, in, 1, stream, mutex);
}
void cellarr_cur(struct cellarr *ca, int i, const curmat *A, cudaStream_t stream);
void cellarr_cuc(struct cellarr *ca, int i, const cucmat *A, cudaStream_t stream);
void cellarr_curcell(struct cellarr *ca, int i, const curcell *A, cudaStream_t stream);
void cellarr_cuccell(struct cellarr *ca, int i, const cuccell *A, cudaStream_t stream);
#endif
