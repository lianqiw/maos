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
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "curmat.h"
#include "pcg.h"

#include "recon.h" /*for  debugging */

#define PRINT_RES 0
#undef TIMING
#define TIMING 0
#if !TIMING
#undef TIC
#undef tic
#undef toc
#define TIC
#define tic
#define toc(A)
#endif
/* dest = a/b */
__global__ static void div_do(float *restrict dest, const float *restrict a, const float *restrict b){
    dest[0]=a[0]/b[0];
}
/* dest = a/b;then b=a*/
__global__ static void div_assign_do(float *restrict dest, const float *restrict a, float *restrict b){
    dest[0]=a[0]/b[0];
    b[0]=a[0];
}
#if PRINT_RES
__global__ static void div_sqrt_do(float *restrict dest, const float *restrict a,  const float *restrict b){
    dest[0]=sqrt(a[0]/b[0]);
}
#endif

/**
   res points to a scalar in device memory. **The caller ha to zero it**
*/
void curcellinn2(float *restrict res, const curcell *A, const curcell *B, cudaStream_t stream){
    //cudaMemsetAsync(res, 0,sizeof(float), stream);
    if(A->m && B->m){
	const int n=A->m->nx*A->m->ny;
	gpu_inn(res, A->m->p, B->m->p, n, stream);
    }else{
	for(int i=0; i<A->nx*A->ny; i++){
	    const curmat *a=A->p[i];
	    const curmat *b=B->p[i];
	    const int n=a->nx*a->ny;
	    gpu_inn(res,a->p,b->p,n,stream);
	}
    }
}
/**
   The PCG algorithm. Copy of lib/pcg.c, but replacing dcell with curcell.
   Timing: 
   curcellinn implemented with blas takes 0.457 ms each call.
   curcellinn implemented with kernel takes 0.193 ms each call with 256x32 split. 0.158 with 128x16 split. 0.137 with 64x64 split.
   curcelladd takes 0.048 ms per call.
   return non-zero during unconvergence.
 */

int gpu_pcg(curcell **px, 
	    G_CGFUN Amul, const void *A, 
	    G_PREFUN Mmul, const void *M, 
	    const curcell *b, int warm, int maxiter,
	    cudaStream_t stream){
    TIC;tic;
    int ans=0;
    curcell *r0=NULL;
    curcell *x0=NULL;/*The initial vector. equals to *px*/
    curcell *z0=NULL;/*Is reference or preconditioned value. */
    float *store;
    float *current;
    int ntot=maxiter*2+2;
#if PRINT_RES == 1
    ntot+=maxiter+1;
#endif
    DO(cudaMalloc(&store, ntot*sizeof(float)));current=store;
    DO(cudaMemsetAsync(store, 0, ntot*sizeof(float),stream));
#if PRINT_RES == 1
    float *r0z0=store;   current++;
    float *r0r0=current; current+=maxiter;
#endif
    float *r0z1=current; current++;
    float *bk=current;   current++;
    float *r0z2=current; current+=maxiter;
    float *ak=current;   current+=maxiter;
 
#if PRINT_RES == 1
    curcellinn2(r0z0, b, b, stream);
    float diff[maxiter+1];
#endif
    /*computes r0=b-A*x0 */
    curcellcp(&r0, b, stream);
    if(!*px){
	*px=curcellnew(b);
    }
    x0=*px;
    if(warm){
	Amul(&r0, A, x0, -1);/*r0=r0+(-1)*A*x0 */
    }else{
	curcellzero(x0, stream);
    }
    curcell *p0=NULL;
    if(Mmul){
	Mmul(&z0,M,r0,stream);
    }else{
	z0=r0;
    }
    curcellcp(&p0, z0, stream);
    curcellinn2(r0z1, r0, z0, stream);
    curcell *Ap=NULL;
    for(int k=0; k<maxiter; k++){
#if PRINT_RES == 1
	if(Mmul){
	    /*res->r0r0=r0'*r0; */
	    curcellinn2(r0r0+k, r0, r0, stream);
	    /*r0r0=sqrt(r0r0/r0z0); */
	    div_sqrt_do<<<1,1,0,stream>>>(r0r0+k, r0r0+k, r0z0);
	}else{
	    /*r0r0=sqrt(r0z1/r0z0); */
	    div_sqrt_do<<<1,1,0,stream>>>(r0r0+k, r0z1, r0z0);
	}
	cudaMemcpyAsync(&diff[k], r0r0+k, sizeof(float), MEMCPY_D2D, stream);
#endif
	if(Ap){
	    curcellzero(Ap, stream);
	}
	Amul(&Ap, A, p0, 1);
	/*ak=r0z1/(p0'*Ap); */
	curcellinn2(ak+k, p0, Ap, stream);
	div_do<<<1,1,0,stream>>>(ak+k, r0z1, ak+k);
	/*put here helps to remove the spikes in performance/wfs. why necessary? */
	CUDA_SYNC_STREAM;
	/*x0=x0+ak*p0 */
	curcelladd2(&x0, p0, ak+k, 1, stream);
	if(k+1==maxiter) break;
	/*r0=r0-ak*Ap */
	curcelladd2(&r0, Ap, ak+k, -1, stream);
	/*preconditioner */
	if(Mmul) Mmul(&z0,M,r0, stream);
	/*r0z2=r0'*z0 */
	curcellinn2(r0z2+k, r0, z0, stream);
	/*bk=r0z2/r0z1; r0z1=r0z2*/
	div_assign_do<<<1,1,0,stream>>>(bk, r0z2+k, r0z1);
	/*p0=bk*p0+z0 */
	curcelladd3(&p0, bk, z0, stream);
	toc("cg");
    }
    /* Instead of check in the middle, we only copy the last result. Improves performance by 20 nm !!!*/
    CUDA_SYNC_STREAM;
#if PRINT_RES == 1
    info2("CG %2d: %.5f ==> %.5f\n", maxiter, diff[0], diff[maxiter-1]);
#endif
    curcellfree(r0); 
    if(Mmul){
	curcellfree(z0);
    }
    curcellfree(Ap);
    curcellfree(p0);
    cudaFree(store);
    return ans;
}
