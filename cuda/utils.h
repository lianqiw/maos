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
extern int NGPU;
extern int *GPUS;
extern int NG1D;
extern int NG2D;
typedef struct{ 
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
    curcell *surfevl;
    curcell *evlopd;
    curcell *evlpsfol;
    curcell *evlpsfcl;
    curcell *evlpsfcl_ngsr;
    curcell *evlopdcov;
    curcell *evlopdcov_ngsr;
    /*for wfsgrad */
    cusparseMatDescr_t wfsspdesc;
    cuwloc_t *powfs;
    cuwfs_t *wfs;
    /*for recon */
    curecon_t *recon;
    /*for moao*/
    curcell *moao_wfs;
    curcell *moao_evl;
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
extern cudata_t **cudata_all;/*use pointer array to avoid misuse. */

void gpu_print_mem(const char *msg);
size_t gpu_get_mem(void);
/**
   switch to the next GPU and update the pointer.
*/
inline void gpu_set(int igpu){
    igpu=igpu%NGPU;
    cudaSetDevice(GPUS[igpu]);
#ifdef __APPLE__
    pthread_setspecific(cudata_key, cudata_all[igpu]);
#else
    cudata=cudata_all[igpu];
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
void gpu_map2dev(cumap_t ***dest, map_t **source, int nps);
void gpu_sp2dev(cusp **dest, dsp *src);
void gpu_spcell2dev(cuspcell **dest, spcell *src);
void gpu_calc_ptt(double *rmsout, double *coeffout, 
		  const double ipcc, const dmat *imcc,
		  const float (*restrict loc)[2], 
		  const int nloc,
		  const float *restrict phi,
		  const float *restrict amp,
		  cudaStream_t stream
		  );
void gpu_calc_ngsmod(double *pttr_out, double *pttrcoeff_out,
		     double *ngsmod_out, int nmod,
		     double MCC_fcp, double ht, double scale,
		     double thetax, double thetay,
		     const double ipcc, const dmat *imcc,
		     const float (*restrict loc)[2], 
		     const int nloc,
		     const float *restrict phi,
		     const float *restrict amp,
		     cudaStream_t stream);
void gpu_loc2dev(float (* restrict *dest)[2], loc_t *src);
void gpu_dbl2dev(float * restrict *dest, double *src, int n);
void gpu_cmp2dev(fcomplex * restrict *dest, dcomplex *src, int n);
void gpu_dmat2dev(float * restrict *dest, dmat *src);
void gpu_dmat2cu(curmat *restrict *dest, dmat *src);
void gpu_dcell2cu(curcell *restrict *dest, dcell *src);
void gpu_cmat2dev(fcomplex * restrict *dest, cmat *src);
void gpu_dbl2flt(float * restrict *dest, double *src, int n);
void gpu_long2dev(int * restrict *dest, long *src, int n);
void gpu_spint2dev(int * restrict *dest, spint *src, int n);
void gpu_int2dev(int * restrict *dest, int *src, int n);
void gpu_spint2int(int * restrict *dest, spint *src, int n);
void gpu_dev2dbl(double * restrict *dest, double alpha, float *src, double beta, int n, cudaStream_t stream);
#if MYSPARSE
void cuspmul (float *y, cusp *A, float *x, float alpha, cudaStream_t stream);
void cusptmul(float *y, cusp *A, float *x, float alpha, cudaStream_t stream);
#else
void cuspmul (float *y, cusp *A, float *x, float alpha, cusparseHandle_t handle);
void cusptmul(float *y, cusp *A, float *x, float alpha, cusparseHandle_t handle);
#endif
__global__ void fscale_do(float *v, int n, float alpha);
void gpu_writeflt(float *p, int nx, int ny, const char *format, ...);
void gpu_writefcmp(fcomplex *p, int nx, int ny, const char *format, ...);
void gpu_writeint(int *p, int nx, int ny, const char *format, ...);
void gpu_muv2dev(cumuv_t *out, MUV_T *in);
void gpu_cur2d(dmat **out, double alpha, const curmat *in, double beta, cudaStream_t stream);
void gpu_cur2s(smat **out, const curmat *in, cudaStream_t stream);
void gpu_cuc2z(zmat **out, const cucmat *in, cudaStream_t stream);
void gpu_curcell2d(dcell **out, double alpha, const curcell *in, double beta, cudaStream_t stream);
void gpu_curcell2s(scell **out, const curcell *in, cudaStream_t stream);
void gpu_cuccell2z(zcell **out, const cuccell *in, cudaStream_t stream);
void cellarr_cur(struct cellarr *ca, const curmat *A, cudaStream_t stream);
void cellarr_cuc(struct cellarr *ca, const cucmat *A, cudaStream_t stream);
void cellarr_curcell(struct cellarr *ca, const curcell *A, cudaStream_t stream);
void cellarr_cuccell(struct cellarr *ca, const cuccell *A, cudaStream_t stream);

__device__ inline float CABS2(fcomplex r){
    const float a=cuCrealf(r);
    const float b=cuCimagf(r);
    return a*a+b*b;
}
inline void gpu_inn(float *restrict res, const float *a, const float *b, const int n, cudaStream_t stream){
    inn_do<<<DIM(n, DIM_REDUCE), DIM_REDUCE*sizeof(float), stream>>>(res, a, b, n);
}
/*somehow I must test both CUDA_ARCH existance and version.*/
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ <200
static __inline__ __device__ float atomicAdd(float* address, float val)
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
#endif
