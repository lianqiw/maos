/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "curmat.h"
#include "cucmat.h"
#include "utils.h"
#include "kernel.h"

void curset(curmat &A, Real alpha, cudaStream_t stream){
    if(A && A.P()){
	set_do<<<DIM(A.Nx()*A.Ny(),256),0,stream>>>(A.P(), alpha, A.Nx()*A.Ny());
    }
}

void curcp(curmat &out, const curmat &in, cudaStream_t stream){
    if(!in){
	cuzero(out, stream);
    }else{
	if(!out){
	    out=curmat(in.Nx(), in.Ny());
	}else{
	    assert(out.N() == in.N());
	}
	cudaMemcpyAsync(out.P(), in.P(), in.N()*sizeof(Real), MEMCPY_D2D, stream);
    }
}
void curcp(curmat &out, const curmat &in){
    if(!in){
	cuzero(out);
    }else{
	if(!out){
	    out=curmat(in.Nx(), in.Ny());
	}else{
	    assert(out.N() == in.N());
	}
	cudaMemcpy(out.P(), in.P(), in.N()*sizeof(Real), MEMCPY_D2D);
    }
}

/**
   out=out*beta+in*alpha;
*/
void curadd(curmat &out, Real alpha, const curmat &in, Real beta, cudaStream_t stream){
    if(!in) return;
    if(!out || alpha==0){
	curcp(out, in, stream);
	if(Z(fabs)(beta-(Real)1)>EPS){
	    scale_do<<<DIM(in.Nx()*in.Ny(), 256),0,stream>>>
		(out.P(), in.Nx()*in.Ny(), beta);
	}
    }else{
	assert(in.N()==out.N());
	add_do<<<DIM(in.Nx()*in.Ny(), 256),0,stream>>>
	    (out.P(), NULL, alpha, in.P(), NULL, beta, in.Nx()*in.Ny());
    }
}
/**
   out=out*beta+abs2(in)*alpha;
*/
void curaddcabs2(curmat &out, Real alpha, const cucmat &in, Real beta, cudaStream_t stream){
    if(!out){
	out=curmat(in.Nx(),in.Ny());
    }else if(alpha==0){
	cuzero(out, stream);
    }
    addcabs2_do<<<DIM(in.Nx()*in.Ny(), 256),0,stream>>>
	(out.P(), alpha, in.P(), beta, in.Nx()*in.Ny());
}
void curscale(curmat &in, Real alpha, cudaStream_t stream){
    if(!in) return;
    if(alpha==0) {
	cuzero(in, stream);
    }else if(Z(fabs)(alpha-(Real)1)>EPS){
	int n=in.Nx()*in.Ny();
	scale_do<<<DIM(n,256), 0, stream>>>(in.P(), n, alpha); 
    }
}

/**
   Computes C = alpha * C + beta * op(A) * B ;
*/
void curmm(curmat &C, Real alpha, const curmat &A, const curmat &B, const char trans[2], Real beta, cublasHandle_t handle){
    int m,n,k,k2;
    cublasOperation_t transa, transb;
    if(trans[0]=='t'){
	m=A.Ny();
	k=A.Nx();
	transa=CUBLAS_OP_T;
    }else{
	m=A.Nx();
	k=A.Ny();
	transa=CUBLAS_OP_N;
    }
    if(trans[1]=='t'){
	n=B.Nx();
	k2=B.Ny();
	transb=CUBLAS_OP_T;
    }else{
	n=B.Ny();
	k2=B.Nx();
	transb=CUBLAS_OP_N;
    }
    if(!C){
	C=curmat(m,n);
    }else{
	assert((C).Nx()==m && (C).Ny()==n);
    }
    if(k!=k2) error("Matrix mismatch\n");
    DO(CUBL(gemm)(handle, transa, transb, m,n,k,
		  &beta, A.P(), A.Nx(), B.P(), B.Nx(), &alpha, C.P(), C.Nx()));
}
/**
   Computes C = alpha * C + beta * op(A) * B ;
*/
void curmv(Real *c, Real alpha, const curmat &A, const Real *b, char trans, Real beta, cublasHandle_t handle){
    CUBL(gemv)(handle, (trans=='t'||trans==1)?CUBLAS_OP_T:CUBLAS_OP_N, A.Nx(), A.Ny(), &beta, A.P(), A.Nx(), b, 1, &alpha, c, 1);
}
void curcellmm(curcell &C, Real alpha, const curcell &A, const curcell &B, 
	       const char trans[2], const double beta, cublasHandle_t handle){
    if(!A || !B) return;
    int ax, az;
    int nx,ny,nz;
    int bz, by;
    if(trans[0]=='n'||trans[0]=='N'){
	nx=A.Nx(); 
	ax=1; az=A.Nx();
	nz=A.Ny();
    }else{ 
	nx=A.Ny();
	az=1; ax=A.Nx();
	nz=A.Nx();
    }
    if(trans[1]=='n'||trans[0]=='N'){
	ny=B.Ny(); 
	bz=1; by=B.Nx();
	if(nz!=B.Nx()) error("miX(mat)ch\n");
    }else{
	ny=B.Nx();
	by=1; bz=B.Nx();
	if(nz!=B.Ny()) error("miX(mat)ch\n");
    }
    if(!C){
	C=curcell(nx,ny);
    }else{
	assert(C.Nx()==nx && C.Ny()==ny);
	cudaStream_t stream;
	cublasGetStream(handle, &stream);
	if(alpha==0){
	    cuzero(C, stream);
	}else if(Z(fabs)(alpha-(Real)1)>EPS){
	    curcellscale(C, alpha, stream);
	}
    }
    for(int iy=0; iy<ny; iy++){
	for(int ix=0; ix<nx; ix++){
	    for(int iz=0; iz<nz; iz++){
		if(A[ix*ax+iz*az]&&B[iz*bz+iy*by]){
		    curmm(C[ix+iy*nx],1.,A[ix*ax+iz*az], 
			  B[iz*bz+iy*by],trans,beta,handle);
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
curmat curmat::trans(stream_t &stream){
    curmat B=curmat(ny, nx);
    transpose<<<dim3(16,16),dim3(16,16),0,stream>>>
	(B.P(), p, nx, ny);
    return B;
}

/*
  A=A*beta+B*alpha;
*/
void curcelladd(curcell &A, Real beta, const curcell &B, Real alpha, cudaStream_t stream){
    if(!B) return;
    if(!A){
	A=B.New();
    }else{
	assert(A.Nx()==B.Nx() && A.Ny()==B.Ny());
    }
    if(A.M() && B.M()){
	curadd(A.M(), beta, B.M(), alpha, stream);
    }else{
	for(int i=0; i<B.Nx()*B.Ny(); i++){
	    curadd(A[i], beta, B[i], alpha,stream);
	}
    }
}

void curadd(curmat &A, Real beta, cudaStream_t stream){
    const int n=A.Nx()*A.Ny();
    add_do<<<DIM(n, 256), 0, stream>>>(A.P(), beta, n);
}
/**
   add a vector to another, scaled by alpha and beta. all in device memory.
   a=a+b*alpha*beta;
*/

/**
   out=out+in*alpha; beta, alpha lives in device memory.
*/
void curadd(curmat &out, const curmat &in, Real *alpha, Real alpha2, cudaStream_t stream){
    if(!out){
	out=curmat(in.Nx(), in.Ny());
    }
    add_do<<<DIM(in.Nx()*in.Ny(), 256),0,stream>>>
	(out.P(), in.P(), alpha, alpha2, in.Nx()*in.Ny());
}


/**
   A=A*beta+B*alpha; beta, alpha lives in device memory.
*/
void curcelladd(curcell &A, const curcell &B, Real* alpha, Real alpha2, cudaStream_t stream){
    if(!B) return;
    if(!A){
	A=B.New();
    }else{
	assert(A.Nx()==B.Nx() && A.Ny()==B.Ny());
    }
    if(A.M() && B.M()){
	curadd(A.M(), B.M(), alpha, alpha2, stream);
    }else{
	for(int i=0; i<B.N(); i++){
	    curadd(A[i], B[i], alpha, alpha2,  stream);
	}
    }
}

/**
   out=out*beta+in; beta, alpha lives in device memory.
*/
void curadd(curmat &out, Real *alpha1, const curmat &in, cudaStream_t stream){
    if(!out){
	out=curmat(in.Nx(), in.Ny());
    }
    add_do<<<DIM(in.Nx()*in.Ny(), 256),0,stream>>>
	(out.P(), alpha1, 1.f, in.P(), in.Nx()*in.Ny());
}

/**
   A=A*alpha1+B*alpha; alpha1, alpha lives in device memory.
*/
void curcelladd(curcell &A, Real* alpha1, const curcell &B, cudaStream_t stream){
    if(!B) return;
    if(!A){
	A=B.New();
    }else{
	assert(A.Nx()==B.Nx() && A.Ny()==B.Ny());
    }
    if(A.M() && B.M()){
	curadd(A.M(), alpha1, B.M(), stream);
    }else{
	for(int i=0; i<B.Nx()*B.Ny(); i++){
	    curadd(A[i], alpha1, B[i],  stream);
	}
    }
}



Real curinn(const curmat &a, const curmat &b, cudaStream_t stream){
    Real *res;
    Real out;
    DO(cudaMalloc(&res, sizeof(Real)));
    cudaMemsetAsync(res, 0, sizeof(Real), stream);
    inn_wrap(res, a.P(), b.P(), a.Nx()*a.Ny(), stream);
    CUDA_SYNC_STREAM;
    cudaMemcpy(&out, res, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaFree(res);//this command is not synchrnous.
    return out;
}

/**
   Sum all the elements in an array.
 */
void cursum2(Real *restrict res,/**<Result in GPU*/
	     const curmat &a,   /**<Source in GPU*/
	     cudaStream_t stream){
    cudaMemsetAsync(res, 0, sizeof(Real), stream);
    sum_wrap(res, a.P(), a.Nx()*a.Ny(), stream);
}
/**
   Sum all the elements in an array, and return a value.
*/
Real cursum(const curmat &a, cudaStream_t stream){
    Real out;//result in CPU.
    Real *res;//result in GPU
    DO(cudaMalloc(&res, sizeof(Real)));
    cudaMemsetAsync(res, 0, sizeof(Real), stream);
    sum_wrap(res, a.P(), a.Nx()*a.Ny(), stream);
    CUDA_SYNC_STREAM;
    cudaMemcpy(&out, res, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaFree(res);
    return out;
}

/**
   Find the maximum value
*/
Real curmax(const curmat &a, cudaStream_t stream){
    Real out;
    Real *res;
    DO(cudaMalloc(&res, sizeof(Real)));
    cudaMemsetAsync(res, 0, sizeof(Real), stream);
    int n=a.Nx()*a.Ny();
    max_wrap(res, a.P(), n, stream);
    CUDA_SYNC_STREAM;
    cudaMemcpy(&out, res, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaFree(res);
    return out;
}

/**
   Find the maximum value
*/
Real curmaxabs(const curmat &a, cudaStream_t stream){
    Real out;
    Real *res;
    DO(cudaMalloc(&res, sizeof(Real)));
    cudaMemsetAsync(res, 0, sizeof(Real), stream);
    int n=a.Nx()*a.Ny();
    maxabs_wrap(res, a.P(), n, stream);
    CUDA_SYNC_STREAM;
    cudaMemcpy(&out, res, sizeof(Real), cudaMemcpyDeviceToHost);
    cudaFree(res);
    return out;
}
/**
   Find the maximum value
*/
Real curcellmax(const curcell &a, cudaStream_t stream){
    int n=a.Nx()*a.Ny();
    Real out;
    Real *res;
    DO(cudaMalloc(&res, (n+1)*sizeof(Real)));
    cudaMemsetAsync(res, 0,(n+1)*sizeof(Real), stream);
    for(int i=0; i<n; i++){
	int m=a[i].N();
	max_wrap(&res[i], a[i].P(), m, stream);
    }
    if(n>1) {
	max_wrap(&res[n], res, n, stream);
    }
    CUDA_SYNC_STREAM;
    cudaMemcpy(&out, &res[n>1?n:0], sizeof(Real), cudaMemcpyDeviceToHost);
    cudaFree(res);
    return out;
}
/**
   Find the maximum value
*/
Real curcellmaxabs(const curcell &a, cudaStream_t stream){
    int n=a.N();
    Real out;
    Real *res;
    DO(cudaMalloc(&res, (n+1)*sizeof(Real)));
    cudaMemsetAsync(res, 0,(n+1)*sizeof(Real), stream);
    for(int i=0; i<n; i++){
	int m=a[i].N();
	maxabs_wrap(&res[i], a[i].P(), m, stream);
    }
    if(n>1) {
	maxabs_wrap(&res[n], res, n, stream);
    }
    CUDA_SYNC_STREAM;
    cudaMemcpy(&out, &res[n>1?n:0], sizeof(Real), cudaMemcpyDeviceToHost);
    cudaFree(res);
    return out;
}
/**
   Scale elements
*/
void curcellscale(curcell &A, Real alpha, cudaStream_t stream){
    if(!A) return;
    if(A.M()){
	curscale(A.M(), alpha, stream);
    }else{
	for(int i=0; i<A.Nx()*A.Ny(); i++){
	    curscale(A[i], alpha, stream);
	}
    }
}
