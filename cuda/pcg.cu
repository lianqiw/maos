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
    if( !isfinite(dest[0])) dest[0]=0;
}
/* dest = a/b;then b=a*/
__global__ static void div_assign_do(float *restrict dest, const float *restrict a, float *restrict b){
    dest[0]=a[0]/b[0];
    b[0]=a[0];
    if (!isfinite(dest[0])) dest[0]=0;
}


/**
   res points to a scalar in device memory. **The caller has to zero it**
*/
void curcellinn_add(float *restrict res, const curcell *A, const curcell *B, cudaStream_t stream){
    //cudaMemsetAsync(res, 0,sizeof(float), stream);
    if(A->m && B->m){
	const int n=A->m->nx*A->m->ny;
	inn_wrap(res, A->m->p, B->m->p, n, stream);
    }else{
	for(int i=0; i<A->nx*A->ny; i++){
	    const curmat *a=A->p[i];
	    const curmat *b=B->p[i];
	    const int n=a->nx*a->ny;
	    inn_wrap(res,a->p,b->p,n,stream);
	}
    }
}
#if PRINT_RES
/* dest = sqrt(a/b); */
__global__ static void div_sqrt_do(float *restrict dest, const float *restrict a, 
				   const float *restrict b){
    dest[0]=sqrt(a[0]/b[0]);
}
static void pcg_residual(float *r0r0, float *r0z1, float *r0z0, curcell *r0, 
			 int precond, cudaStream_t stream){
    if(precond){
	/*r0r0=r0'*r0; */
	curcellinn_add(r0r0, r0, r0, stream);
	/*r0r0=sqrt(r0r0/r0z0); */
	div_sqrt_do<<<1,1,0,stream>>>(r0r0, r0r0, r0z0);
    }else{
	/*r0r0=sqrt(r0z1/r0z0); */
	div_sqrt_do<<<1,1,0,stream>>>(r0r0, r0z1, r0z0);
    }
}
#endif
/**
   The PCG algorithm. Copy of lib/pcg.c, but replacing dcell with curcell.
   Timing: 
   curcellinn implemented with blas takes 0.457 ms each call.
   curcellinn implemented with kernel takes 0.193 ms each call with 256x32 split. 0.158 with 128x16 split. 0.137 with 64x64 split.
   curcelladd takes 0.048 ms per call.
   return non-zero during unconvergence.
 */

float gpu_pcg(curcell **px, 
	    G_CGFUN Amul, const void *A, 
	    G_PREFUN Mmul, const void *M, 
	    const curcell *b, int warm, int maxiter,
	    cudaStream_t stream){
    TIC;tic;
    curcell *r0=NULL;
    curcell *x0=NULL;/*The initial vector. equals to *px*/
    curcell *z0=NULL;/*Is reference or preconditioned value. */
    float *store;
    float *current;
    float residual=-1;/*Only useful if PRINT_RES is set*/
    int ntot=maxiter*2+2;
#if PRINT_RES 
    ntot+=maxiter+1;
    float diff[maxiter];
#endif
    DO(cudaMalloc(&store, ntot*sizeof(float)));current=store;
    DO(cudaMemsetAsync(store, 0, ntot*sizeof(float),stream));
#if PRINT_RES 
    float *r0z0=store;   current++;
    float *r0r0=current; current+=maxiter;
#endif
    float *r0z1=current; current++;
    float *bk=current;   current++;
    float *r0z2=current; current+=maxiter;
    float *ak=current;   current+=maxiter;
 
#if PRINT_RES
    curcellinn_add(r0z0, b, b, stream);
#if PRINT_RES == 2
    fprintf(stderr, "CG %d:", maxiter);
#endif
#endif
    /*computes r0=b-A*x0 */
    if(!*px){
	*px=curcellnew(b);
    }else if(!warm){
	curcellzero(*px, stream);
    }
    x0=*px;
    curcell *p0=NULL;
    curcell *Ap=NULL;
    for(int k=0; k<maxiter; k++){
	if(k%50==0){/*initial or re-start every 50 steps*/
	    /*computes r0=b-A*x0 */
	    curcellcp(&r0, b, stream);/*r0=b; */
	    CUDA_SYNC_STREAM;
	    Amul(&r0, 1, A, x0, -1);/*r0=r0+(-1)*A*x0 */
	    if(Mmul){/*z0=M*r0*/
		Mmul(&z0,M,r0,stream);
	    }else{
		z0=r0;
	    }
	    curcellcp(&p0, z0, stream);
	    if(k!=0){
		cudaMemsetAsync(r0z1, 0, sizeof(float), stream);
	    }
	    curcellinn_add(r0z1, r0, z0, stream);
	}
#if PRINT_RES 
	pcg_residual(r0r0+k, r0z1, r0z0, r0, Mmul?1:0, stream);
	cudaMemcpyAsync(&diff[k], r0r0+k, sizeof(float), MEMCPY_D2H, stream);
#if PRINT_RES == 2
	info2("%.5f ", diff[k]);
#endif	
#endif
	CUDA_SYNC_STREAM;//needed before entering Amul
	Amul(&Ap, 0, A, p0, 1);
	/*ak=r0z1/(p0'*Ap); */
	curcellinn_add(ak+k, p0, Ap, stream);
	div_do<<<1,1,0,stream>>>(ak+k, r0z1, ak+k);
	/*x0=x0+ak*p0 */
	curcelladd(&x0, p0, ak+k, 1, stream);
	if(k+1==maxiter) {
	    /* Compute residual. Use bk as a temporary*/
	    /*pcg_residual(bk, r0z1, r0z0, r0, Mmul?1:0, stream);
	    CUDA_SYNC_STREAM;
	    cudaMemcpy(&residual, bk, sizeof(float), MEMCPY_D2H);*/
#if PRINT_RES 
	    residual=diff[k];
#endif
	    break;
	}
	/*r0=r0-ak*Ap */
	curcelladd(&r0, Ap, ak+k, -1, stream);
	/*preconditioner */
	if(Mmul) {/*z0=M*r0*/
	    Mmul(&z0,M,r0, stream);
	}
	/*r0z2=r0'*z0 */
	curcellinn_add(r0z2+k, r0, z0, stream);
	/*bk=r0z2/r0z1; r0z1=r0z2*/
	div_assign_do<<<1,1,0,stream>>>(bk, r0z2+k, r0z1);
	/*p0=bk*p0+z0 */
	curcelladd(&p0, bk, z0, stream);

	toc("cg");
	/*{
	    curcellwrite(Ap, "Ap");
	    curcellwrite(p0, "p0");
	    curcellwrite(z0, "z0");
	    curcellwrite(r0, "r0");
	    exit(0);
	    }*/
    }
    /* Instead of check in the middle, we only copy the last result. Improves performance by 20 nm !!!*/
    CUDA_SYNC_STREAM;
#if PRINT_RES == 1
    fprintf(stderr, "CG %2d: %.5f ==> %.5f\n", maxiter, diff[0], diff[maxiter-1]);
#elif PRINT_RES==2
    fprintf(stderr, "\n");
#endif
    curcellfree(r0); 
    if(Mmul){
	curcellfree(z0);
    }
    curcellfree(Ap);
    curcellfree(p0);
    cudaFree(store);
    return residual;
}
