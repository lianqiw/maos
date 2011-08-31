extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "curmat.h"
#include "pcg.h"
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
__global__ static void assign_do(float *dest, const float *restrict res){
    dest[0]=res[0];
}
__global__ static void div_do(float *dest, const float *restrict a, const float *restrict b){
    dest[0]=a[0]/b[0];
}
#if PRINT_RES
__global__ static void div2_do(float *dest,  const float *restrict b){
    dest[0]/=b[0];
}
__global__ static void div_sqrt_do(float *dest,  const float *restrict b){
    dest[0]=sqrt(dest[0]/b[0]);
}
#endif
__global__ static void div_left_do(float *dest,  const float *restrict b){
    dest[0]=b[0]/dest[0];
}
__global__ static void scale_do(float *dest,  float b){
    dest[0]*=b;
}
/**
   The PCG algorithm. Copy of lib/pcg.c, but replacing dcell with curcell.
   Timing: 
   curcellinn implemented with blas takes 0.457 ms each call.
   curcellinn implemented with kernel takes 0.193 ms each call with 256x32 split. 0.158 with 128x16 split. 0.137 with 64x64 split.
   curcelladd takes 0.048 ms per call.
 */

void gpu_pcg(curcell **px, 
	     G_CGFUN Amul, const void *A, 
	     G_PREFUN Mmul, const void *M, 
	     const curcell *b, int warm, int maxiter,
	     cudaStream_t stream){
    TIC;tic;
    curcell *r0=NULL;
    curcell *x0=NULL;//The initial vector.
    curcell *z0=NULL;//Is reference or preconditioned value.
    typedef struct{
	float r0z0;
	float r0z1;
	float r0z2;
	float ak;
	float bk;
	float tmp;
    }CGRES_T;
    CGRES_T *res;
    DO(cudaMalloc((float**)&res, sizeof(CGRES_T)));
    //computes r0=b-A*x0
    curcellcp(&r0, b, stream);
    if(!*px || !warm){//start from zero guess.
	x0=curcellnew2(b);
	if(!*px) *px=curcellnew(x0->nx, x0->ny);
    }else{
	curcellcp(&x0, *px, stream);
	CUDA_SYNC_STREAM;
	Amul(&r0, A, x0, -1);//r0=r0+(-1)*A*x0
    }
    curcell *p0=NULL;
    if(Mmul){
	Mmul(&z0,M,r0);
    }else{
	z0=r0;
    }
    curcellcp(&p0, z0, stream);
    curcellinn2(&res->r0z1, r0, z0, stream);
    curcell *Ap=NULL;
#if PRINT_RES == 1
    curcellinn2(&res->r0z0, b, b, stream);
    float diff[maxiter+1];
    if(Mmul){
	//double tmp=curcellinn(r0, r0, stream);
	curcellinn2(&res->tmp, r0, r0, stream);
	div2_do<<<1,1,0,stream>>>(&res->tmp, &res->r0z0);
	//diff[0]=sqrt(tmp/r0z0);
	cudaMemcpyAsync(&diff[0], &res->tmp, 1, sizeof(float), stream);
    }else{
	diff[0]=sqrt(r0z1/r0z0);
    }
    int kres=0;
    warning2("Step %d, res=%g\n", kres, diff[kres]);
#endif
    for(int k=0; k<maxiter; k++){
	if(Ap) curcellzero(Ap, stream);
	Amul(&Ap, A, p0, 1);
	//ak=r0z1/curcellinn(p0,Ap,stream);
	curcellinn2(&res->ak, p0, Ap, stream);
	div_left_do<<<1,1,0,stream>>>(&res->ak, &res->r0z1);
	//x0=x0+ak*p0
	//r0=r0-ak*Ap
	curcelladd2(&x0, p0, &res->ak, stream);//x0=x0+ak*p0
	scale_do<<<1,1,0,stream>>>(&res->ak, -1);
	curcelladd2(&r0, Ap, &res->ak, stream);//r0=r0-ak*Ap
	
	if(Mmul){
	    Mmul(&z0,M,r0);
	}else{
	    z0=r0;
	}
	//r0z2=r0'*z0
	curcellinn2(&res->r0z2, r0, z0, stream);
	//bk=r0z2/r0z1;
	div_do<<<1,1,0,stream>>>(&res->bk, &res->r0z2, &res->r0z1);
	//p0=bk*p0+z0
	curcelladd3(&p0, &res->bk, z0, stream);
	//r0z1=r0z2;
	assign_do<<<1,1,0,stream>>>(&res->r0z1, &res->r0z2);
#if PRINT_RES == 1
	if(Mmul){
	    //diff[k+1]=sqrt(r0'*r0/r0z0)
	    curcellinn2(&res->tmp, r0, r0, stream);
	    div_sqrt_do(&res->tmp, &res->tmp, &res->r0z0);
	    cudaMemcpyAsync(&diff[k+1], &res->tmp, sizeof(float), stream);
	}else{
	    //diff[k+1]=sqrt(r0z2/r0z0);
	    div_sqrt_do(&res->tmp, &res->r0z2, &res->r0z0);
	    cudaMemcpyAsync(&diff[k+1], &res->tmp, sizeof(float), stream);
	}
	warning2("Step %d, res=%g\n", k+1, diff[k+1]);
#endif
	toc("cg");
    }
    /* Instead of check in the middle, we only copy the last result. Improves performance by 20 nm !!!*/
    curcellcp(px, x0, stream);
#if PRINT_RES == 1
    info("Solution found at step %d with residual %g\n", kres+1, diff[kres+1]);
#endif
    curcellfree(r0); 
    if(Mmul){
	curcellfree(z0);
    }
    curcellfree(x0);
    curcellfree(Ap);
    curcellfree(p0);
    cudaFree(res);
}
