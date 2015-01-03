/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "gpu.h"
}
#include "curmat.h"
#include "cucmat.h"
#include "utils.h"
#include "kernel.h"

void curset(curmat *A, Real alpha, cudaStream_t stream){
    if(A && A->p){
	set_do<<<DIM(A->nx*A->ny,256),0,stream>>>(A->p, alpha, A->nx*A->ny);
    }
}

void curcp(curmat **out, const curmat *in, cudaStream_t stream){
    if(!in){
	curzero(*out, stream);
    }else{
	if(!*out){
	    *out=curnew(in->nx, in->ny);
	}else{
	    assert((*out)->nx * (*out)->ny == in->nx * in->ny);
	}
	cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(Real), MEMCPY_D2D, stream);
    }
}
void curcp(curmat **out, const curmat *in){
    if(!in){
	curzero(*out);
    }else{
	if(!*out){
	    *out=curnew(in->nx, in->ny);
	}else{
	    assert((*out)->nx * (*out)->ny == in->nx * in->ny);
	}
	cudaMemcpy((*out)->p, in->p, in->nx*in->ny*sizeof(Real), MEMCPY_D2D);
    }
}

/**
   out=out*beta+in*alpha;
*/
void curadd(curmat **out, Real alpha, curmat *in, Real beta, cudaStream_t stream){
    if(!in) return;
    if(!*out || Z(fabs)(alpha)<EPS){
	curcp(out, in, stream);
	if(Z(fabs)(beta-1)>EPS){
	    scale_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
		((*out)->p, in->nx*in->ny, beta);
	}
    }else{
	assert(in->nx*in->ny==(*out)->nx*(*out)->ny);
	add_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
	    ((*out)->p, NULL, alpha, in->p, NULL, beta, in->nx*in->ny);
    }
}
/**
   out=out*beta+abs2(in)*alpha;
*/
void curaddcabs2(curmat **out, Real alpha, cucmat *in, Real beta, cudaStream_t stream){
    if(!*out){
	*out=curnew(in->nx,in->ny);
    }else if(Z(fabs)(alpha)<EPS){
	curzero(*out, stream);
    }
    addcabs2_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
	((*out)->p, alpha, in->p, beta, in->nx*in->ny);
}
void curscale(curmat *in, Real alpha, cudaStream_t stream){
    if(!in) return;
    if(Z(fabs)(alpha)<EPS) {
	curzero(in, stream);
    }else if(Z(fabs)(alpha-1.f)>EPS){
	int n=in->nx*in->ny;
	scale_do<<<DIM(n,256), 0, stream>>>(in->p, n, alpha); 
    }
}

/**
   Computes C = alpha * C + beta * op(A) * B ;
*/
void curmm(curmat **C, Real alpha, const curmat *A, const curmat *B, const char trans[2], Real beta, cublasHandle_t handle){
    int m,n,k,k2;
    cublasOperation_t transa, transb;
    if(trans[0]=='t'){
	m=A->ny;
	k=A->nx;
	transa=CUBLAS_OP_T;
    }else{
	m=A->nx;
	k=A->ny;
	transa=CUBLAS_OP_N;
    }
    if(trans[1]=='t'){
	n=B->nx;
	k2=B->ny;
	transb=CUBLAS_OP_T;
    }else{
	n=B->ny;
	k2=B->nx;
	transb=CUBLAS_OP_N;
    }
    if(!*C){
	*C=curnew(m,n);
    }else{
	assert((*C)->nx==m && (*C)->ny==n);
    }
    if(k!=k2) error("Matrix mismatch\n");
    DO(CUBL(gemm)(handle, transa, transb, m,n,k,
		 &beta, A->p, A->nx, B->p, B->nx, &alpha, (*C)->p, (*C)->nx));
}
/**
   Computes C = alpha * C + beta * op(A) * B ;
*/
void curmv(Real *c, Real alpha, const curmat *A, const Real *b, char trans, Real beta, cublasHandle_t handle){
    CUBL(gemv)(handle, (trans=='t'||trans==1)?CUBLAS_OP_T:CUBLAS_OP_N, A->nx, A->ny, &beta, A->p, A->nx, b, 1, &alpha, c, 1);
}
void curcellmm(curcell **C0, double alpha, const curcell *A, const curcell *B, 
	       const char trans[2], const double beta, cublasHandle_t handle){
    if(!A || !B) return;
    int ax, az;
    int nx,ny,nz;
    int bz, by;
    if(trans[0]=='n'||trans[0]=='N'){
	nx=A->nx; 
	ax=1; az=A->nx;
	nz=A->ny;
    }else{ 
	nx=A->ny;
	az=1; ax=A->nx;
	nz=A->nx;
    }
    if(trans[1]=='n'||trans[0]=='N'){
	ny=B->ny; 
	bz=1; by=B->nx;
	if(nz!=B->nx) error("miX(mat)ch\n");
    }else{
	ny=B->nx;
	by=1; bz=B->nx;
	if(nz!=B->ny) error("miX(mat)ch\n");
    }
    if(!*C0){
	*C0=curcellnew(nx,ny);
    }else{
	assert((*C0)->nx==nx && (*C0)->ny==ny);
	cudaStream_t stream;
	cublasGetStream(handle, &stream);
	if(Z(fabs)(alpha)<EPS){
	    curcellzero(*C0, stream);
	}else if(Z(fabs)(alpha-1)>EPS){
	    curcellscale(*C0, alpha, stream);
	}
    }
    curcell *C=*C0;
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    for(int iz=0; iz<nz; iz++){
		if(A->p[ix*ax+iz*az]&&B->p[iz*bz+iy*by]){
		    curmm(&C->p[ix+iy*nx],1.,A->p[ix*ax+iz*az], 
			  B->p[iz*bz+iy*by],trans,beta,handle);
		}
	    }
	}
    }
}
/*Transpose a matrix in naive way. Faster way is to use shared memory and handle
  a block each time.*/
__global__ void transpose(Real *restrict out, const Real *restrict in, int nx, int ny){
    const int stepx=blockDim.x * gridDim.x;
    const int stepy=blockDim.y * gridDim.y;
    const int ix0=threadIdx.x+blockDim.x*blockIdx.x;
    const int iy0=threadIdx.y+blockDim.y*blockIdx.y;
    for(int iy=iy0; iy<ny; iy+=stepy){
	for(int ix=ix0; ix<nx; ix+=stepx){
	    out[iy+ix*ny]=in[ix+iy*nx];
	}
    }
}
/*Transpose a matrix*/
template <>
curmat* curmat::trans(stream_t &stream){
    curmat *B=curnew(ny, nx);
    transpose<<<dim3(16,16),dim3(16,16),0,stream>>>
	(B->p, p, nx, ny);
    return B;
}

/*
  A=A*beta+B*alpha;
*/
void curcelladd(curcell **A, Real beta, const curcell *B, Real alpha, cudaStream_t stream){
    if(!B) return;
    if(!*A){
	*A=curcellnew(B);
    }else{
	assert((*A)->nx==B->nx && (*A)->ny==B->ny);
    }
    if((*A)->m && B->m){
	curadd(&(*A)->m, beta, B->m, alpha, stream);
    }else{
	for(int i=0; i<B->nx*B->ny; i++){
	    curadd(&((*A)->p[i]), beta, B->p[i], alpha,stream);
	}
    }
}

void curadd(curmat *A, Real beta, cudaStream_t stream){
    const int n=A->nx*A->ny;
    add_do<<<DIM(n, 256), 0, stream>>>(A->p, beta, n);
}
/**
   add a vector to another, scaled by alpha and beta. all in device memory.
   a=a+b*alpha*beta;
*/

/**
   out=out+in*alpha; beta, alpha lives in device memory.
*/
void curadd(curmat **out, const curmat *in, Real *alpha, Real alpha2, cudaStream_t stream){
    if(!*out){
	*out=curnew(in->nx, in->ny);
    }
    add_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
	((*out)->p, in->p, alpha, alpha2, in->nx*in->ny);
}


/**
   A=A*beta+B*alpha; beta, alpha lives in device memory.
*/
void curcelladd(curcell **A, const curcell *B, Real* alpha, Real alpha2, cudaStream_t stream){
    if(!B) return;
    if(!*A){
	*A=curcellnew(B);
    }else{
	assert((*A)->nx==B->nx && (*A)->ny==B->ny);
    }
    if((*A)->m && B->m){
	curadd(&(*A)->m, B->m, alpha, alpha2, stream);
    }else{
	for(int i=0; i<B->nx*B->ny; i++){
	    curadd((*A)->p+i, B->p[i], alpha, alpha2,  stream);
	}
    }
}

/**
   out=out*beta+in; beta, alpha lives in device memory.
*/
void curadd(curmat **out, Real *alpha1, const curmat *in, cudaStream_t stream){
    if(!*out){
	*out=curnew(in->nx, in->ny);
    }
    add_do<<<DIM(in->nx*in->ny, 256),0,stream>>>
	((*out)->p, alpha1, 1.f, in->p, in->nx*in->ny);
}

/**
   A=A*alpha1+B*alpha; alpha1, alpha lives in device memory.
*/
void curcelladd(curcell **A, Real* alpha1, const curcell *B, cudaStream_t stream){
    if(!B) return;
    if(!*A){
	*A=curcellnew(B);
    }else{
	assert((*A)->nx==B->nx && (*A)->ny==B->ny);
    }
    if((*A)->m && B->m){
	curadd(&(*A)->m, alpha1, B->m, stream);
    }else{
	for(int i=0; i<B->nx*B->ny; i++){
	    curadd(&((*A)->p[i]), alpha1, B->p[i],  stream);
	}
    }
}



Real curinn(const curmat *a, const curmat *b, cudaStream_t stream){
    Real *res;
    Real out;
    DO(cudaMalloc(&res, sizeof(Real)));
    cudaMemsetAsync(res, 0, sizeof(Real), stream);
    inn_wrap(res, a->p, b->p, a->nx*a->ny, stream);
    cudaMemcpyAsync(&out, res, sizeof(Real), cudaMemcpyDeviceToHost, stream);
    cudaFree(res);//this command is not synchrnous.
    return out;
}

/**
   Sum all the elements in an array.
 */
void cursum2(Real *restrict res, const curmat *a, cudaStream_t stream){
    cudaMemsetAsync(res, 0,sizeof(Real), stream);
    sum_wrap(res, a->p, a->nx*a->ny, stream);
}


/**
   Find the maximum value
*/
Real curmax(const curmat *a, cudaStream_t stream){
    Real out;
    Real *res;
    DO(cudaMalloc(&res, sizeof(Real)));
    cudaMemsetAsync(res, 0, sizeof(Real), stream);
    int n=a->nx*a->ny;
    max_wrap(res, a->p, n, stream);
    cudaMemcpyAsync(&out, res, sizeof(Real), cudaMemcpyDeviceToHost, stream);
    cudaFree(res);
    return out;
}

/**
   Find the maximum value
*/
Real curmaxabs(const curmat *a, cudaStream_t stream){
    Real out;
    Real *res;
    DO(cudaMalloc(&res, sizeof(Real)));
    cudaMemsetAsync(res, 0, sizeof(Real), stream);
    int n=a->nx*a->ny;
    maxabs_wrap(res, a->p, n, stream);
    cudaMemcpyAsync(&out, res, sizeof(Real), cudaMemcpyDeviceToHost, stream);
    cudaFree(res);
    return out;
}
/**
   Find the maximum value
*/
Real curcellmax(const curcell *a, cudaStream_t stream){
    int n=a->nx*a->ny;
    Real out;
    Real *res;
    DO(cudaMalloc(&res, (n+1)*sizeof(Real)));
    cudaMemsetAsync(res, 0,(n+1)*sizeof(Real), stream);
    for(int i=0; i<n; i++){
	int m=a->p[i]->nx*a->p[i]->ny;
	max_wrap(&res[i], a->p[i]->p, m, stream);
    }
    if(n>1) {
	max_wrap(&res[n], res, n, stream);
    }
    cudaMemcpyAsync(&out, &res[n>1?n:0], sizeof(Real), cudaMemcpyDeviceToHost, stream);
    cudaFree(res);
    return out;
}
/**
   Find the maximum value
*/
Real curcellmaxabs(const curcell *a, cudaStream_t stream){
    int n=a->nx*a->ny;
    Real out;
    Real *res;
    DO(cudaMalloc(&res, (n+1)*sizeof(Real)));
    cudaMemsetAsync(res, 0,(n+1)*sizeof(Real), stream);
    for(int i=0; i<n; i++){
	int m=a->p[i]->nx*a->p[i]->ny;
	maxabs_wrap(&res[i], a->p[i]->p, m, stream);
    }
    if(n>1) {
	maxabs_wrap(&res[n], res, n, stream);
    }
    cudaMemcpyAsync(&out, &res[n>1?n:0], sizeof(Real), cudaMemcpyDeviceToHost, stream);
    cudaFree(res);
    return out;
}
/**
   Scale elements
*/
void curcellscale(curcell *A, Real alpha, cudaStream_t stream){
    if(!A) return;
    if(A->m){
	curscale(A->m, alpha, stream);
    }else{
	for(int i=0; i<A->nx*A->ny; i++){
	    curscale(A->p[i], alpha, stream);
	}
    }
}
