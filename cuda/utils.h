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
#ifndef AOS_CUDA_UTILS_H
#define AOS_CUDA_UTILS_H
#include "types.h"
#include <cublas_v2.h>
#include <cusparse.h>
#include <cufft.h>
#include "recon.h"
#include "kernel.h"
#include "wfs.h"
extern int gpu_recon;
extern int NGPU;
extern int *GPUS;
typedef struct cudata_t{ 
    /**<for accphi */
    cumap_t **atm;   /**<atmosphere: array of cumap_t */
    cumap_t **dmreal;/**<DM: array of cumap_t */
    cumap_t **dmproj;/**<DM: array of cumap_t */
    int nps; /**<number of phase screens*/
    int ndm; /**<number of DM.*/
    /*for perfevl */
    float  (*plocs)[2];
    float   *pamp;
    int    **embed;
    cuccell *evlwvf;
    curcell *surfevl;
    curcell *evlopd;
    curcell *evlpsfol;
    curcell *evlpsfcl;
    curcell *evlpsfcl_ngsr;
    curcell *evlopdcov;
    curmat *evlopdcovol;
    curcell *evlopdcov_ngsr;
    curcell *evlopdmean;
    curmat *evlopdmeanol;
    curcell *evlopdmean_ngsr;
    /*for wfsgrad */
    cuwloc_t *powfs;
    cuwfs_t *wfs;
    /*for recon */
    curecon_t *recon;
    /*for moao*/
    curcell *dm_wfs;
    curcell *dm_evl;
    /*for mvm*/
    curmat *mvm_m;
    curmat *mvm_a;
    curmat *mvm_g;
    stream_t *mvm_stream;
    pthread_mutex_t mvm_mutex;
    cudata_t(){
	memset(this, 0, sizeof(cudata_t));
    }
}cudata_t;
#ifdef __APPLE__
extern pthread_key_t cudata_key;
inline cudata_t* _cudata(){
    return (cudata_t*)pthread_getspecific(cudata_key);
}
#define cudata _cudata()
#else
extern __thread cudata_t *cudata;
#endif
extern cudata_t *cudata_all;/*use pointer array to avoid misuse. */

void gpu_print_mem(const char *msg);
long gpu_get_mem(void);
/**
   switch to the next GPU and update the pointer.
*/
inline void gpu_set(int igpu){
    igpu=igpu%NGPU;
    cudaSetDevice(GPUS[igpu]);
#ifdef __APPLE__
    pthread_setspecific(cudata_key, &cudata_all[igpu]);
#else
    cudata=&cudata_all[igpu];
#endif
}
/**
   returns next available GPU. Useful for assigning GPUs to particular wfs, evl, etc.
*/
inline int gpu_next(){
    static int cur=-1;
    return cur=(cur+1)%NGPU;
}
/*void gpu_set(int igpu);
  int  gpu_next(void);*/
void dbl2flt(float * restrict *dest, double *src, int n);
void spint2int(int * restrict *dest, spint *src, int n);

void cp2gpu(cumap_t ***dest, map_t **source, int nps);
void cp2gpu(cusp **dest, dsp *src);
void cp2gpu(cuspcell **dest, spcell *src);
void cp2gpu(float (* restrict *dest)[2], loc_t *src);
void cp2gpu(float * restrict *dest, double *src, int n);
void cp2gpu(float * restrict *dest, dmat *src);
void cp2gpu(fcomplex * restrict *dest, dcomplex *src, int n);
void cp2gpu(fcomplex * restrict *dest, cmat *src);

void cp2gpu(curmat *restrict *dest, dmat *src);
void cp2gpu(curmat *restrict *dest, smat *src, cudaStream_t stream=0);
void cp2gpu(curcell *restrict *dest, dcell *src);
void cp2gpu(cucmat *restrict *dest, cmat *src);
void cp2gpu(cuccell *restrict *dest, ccell *src);

void cp2gpu(int * restrict *dest, long *src, int n);
void cp2gpu(int * restrict *dest, spint *src, int n);
void cp2gpu(int * restrict *dest, int *src, int n);
void cp2gpu(cumuv_t *out, MUV_T *in);


#if MYSPARSE
void cuspmul (float *y, cusp *A, float *x, float alpha, cudaStream_t stream);
void cusptmul(float *y, cusp *A, float *x, float alpha, cudaStream_t stream);
#else
void cuspmul (float *y, cusp *A, float *x, float alpha, cusparseHandle_t handle);
void cusptmul(float *y, cusp *A, float *x, float alpha, cusparseHandle_t handle);
#endif
void gpu_write(float *p, int nx, int ny, const char *format, ...);
void gpu_write(fcomplex *p, int nx, int ny, const char *format, ...);
void gpu_write(int *p, int nx, int ny, const char *format, ...);

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
void cellarr_cur(struct cellarr *ca, const curmat *A, cudaStream_t stream);
void cellarr_cuc(struct cellarr *ca, const cucmat *A, cudaStream_t stream);
void cellarr_curcell(struct cellarr *ca, const curcell *A, cudaStream_t stream);
void cellarr_cuccell(struct cellarr *ca, const cuccell *A, cudaStream_t stream);
#endif
