extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"

/**
   Createa curmat object.
*/
curmat *curnew(int nx, int ny){
    curmat *out;
    out=(curmat*)calloc(1, sizeof(curmat));
    DO(cudaMalloc(&(out->p), nx*ny*sizeof(float)));
    DO(cudaMemset(out->p, 0, nx*ny*sizeof(float)));
    out->nx=nx;
    out->ny=ny;
    return out;
}
void curfree(curmat *A){
    if(A){
	if(A->p){
	    cudaFree(A->p);
	}
	free(A);
    }
}
void curzero(curmat *A, cudaStream_t stream){
    if(A && A->p){
	DO(cudaMemsetAsync(A->p, 0, A->nx*A->ny*sizeof(float), stream));
    }
}
void curcp(curmat **out, const curmat *in, cudaStream_t stream){
    if(!in){
	curzero(*out, stream);
    }else{
	if(!*out){
	    *out=curnew(in->nx, in->ny);
	}else{
	    assert((*out)->nx==in->nx && (*out)->ny==in->ny);
	}
	cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(float), cudaMemcpyDefault, stream);
    }
}
void curwritedata(const curmat *A, file_t *fp){
    if(A && A->nx >0 && A->ny>0){
	cudaDeviceSynchronize();
	float *tmp=(float*)malloc(A->nx*A->ny*sizeof(float));
	cudaMemcpy(tmp, A->p, A->nx*A->ny*sizeof(float), cudaMemcpyDefault);
	cudaDeviceSynchronize();
	do_write(fp, 0, sizeof(float), M_FLT, tmp, A->nx, A->ny);
	free(tmp);
    }else{
	do_write(fp, 0, sizeof(float), M_FLT, NULL, 0, 0);
    }
}
void curwrite(const curmat *A, const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn, "wb");
    curwritedata(A, fp);
    zfclose(fp);
}
/**
   out=out*beta+in*alpha;
*/
void curadd(curmat **out, float beta, curmat *in, float alpha, cublasHandle_t handle){
    if(fabsf(beta-1)>1e-6 && *out){
	cublasSscal(handle, (*out)->nx*(*out)->ny, &beta, (*out)->p, 1);
    }
    if(!*out) *out=curnew(in->nx, in->ny);
    cublasSaxpy(handle, in->nx*in->ny, &alpha, in->p, 1, (*out)->p, 1);
}
__global__ static void scale_do(float *restrict in, int n, float alpha){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	in[i]*=alpha;
    }
}
void curscale(curmat *in, float alpha, cudaStream_t stream){
    int n=in->nx*in->ny;
    scale_do<<<MAX(MIN(n/256, 32), 1), MIN(n, 256), 0, stream>>>(in->p, n, alpha); 
}
/**
   Computes y = alpha * y + beta * op(A) * x; x, y are vectors
*/
void curmv(curmat **C, float alpha, const curmat *A, const curmat *B, char trans, float beta, cublasHandle_t handle){
    if(!*C){
	*C=curnew(A->nx, 1);
    }else{
	assert((*C)->nx==A->nx);
    }
    cublasSgemv(handle, trans=='t'?CUBLAS_OP_T:CUBLAS_OP_N, A->nx, A->ny, &beta, A->p, A->nx, B->p, 1, &alpha, (*C)->p, 1);
}
/**
   Computes C = alpha * C + beta * op(A) * B ;
*/
void curmm(curmat **C, float alpha, const curmat *A, const curmat *B, char trans[2], float beta, cublasHandle_t handle){
    if(!*C){
	*C=curnew(A->nx, B->ny);
    }else{
	assert((*C)->nx==A->nx && (*C)->ny==B->ny && A->ny==B->nx);
    }
    if(B->ny==1){
	cublasSgemv(handle, trans[0]=='t'?CUBLAS_OP_T:CUBLAS_OP_N, A->nx, A->ny, &beta, A->p, A->nx, B->p, 1, &alpha, (*C)->p, 1);
    }else{
	cublasSgemm(handle, trans[0]=='t'?CUBLAS_OP_T:CUBLAS_OP_N, trans[1]=='t'?CUBLAS_OP_T:CUBLAS_OP_N, A->nx, B->ny, A->ny, &beta, A->p, A->nx, B->p, B->nx, &alpha, (*C)->p, (*C)->nx);
    }
}

curcell* curcellnew(int nx, int ny){
    curcell *out=(curcell*)calloc(1, sizeof(curcell));
    out->p=(curmat**)calloc(nx*ny, sizeof(void*));
    out->nx=nx;
    out->ny=ny;
    return out;
}

curcell *curcellnew2(const curcell *in){
    curcell *out=curcellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	out->p[i]=curnew(in->p[i]->nx, in->p[i]->ny);
    }
    return out;
}

void curcellfree(curcell *A){
    if(A){
	if(A->p) free(A->p);
	free(A);
    }
}

void curcellwrite(const curcell *A, const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn, "wb");
    write_magic(MCC_ANY, fp);
    if(A){
	uint64_t nx=A->nx;
	uint64_t ny=A->ny;
	zfwritelarr(fp, 2, &nx, &ny);
	for(int i=0; i<nx*ny; i++){
	    curwritedata(A->p[i], fp);
	}
    }else{
	uint64_t zero=0;
	zfwritelarr(fp, 2, &zero, &zero);
    }
    zfclose(fp);	
}
void curcellzero(curcell *A, cudaStream_t stream){
    if(!A) return;
    for(int i=0; i<A->nx*A->ny; i++){
	curzero(A->p[i], stream);
    }
}
void curcellcp(curcell **A, const curcell *B, cudaStream_t stream){
    if(!B)
	curcellzero(*A, stream);
    else{
	if(!*A){
	    *A=curcellnew(B->nx, B->ny);
	}else{
	    assert((*A)->nx==B->nx && (*A)->ny==B->ny);
	}
	for(int i=0; i<B->nx*B->ny; i++){
	    curcp(&(*A)->p[i], B->p[i], stream);
	}
    }
}
/*
  A=A*beta+B*alpha;
*/
void curcelladd(curcell **A, float beta, const curcell *B, float alpha, cublasHandle_t handle){
    if(!B) return;
    if(!*A){
	*A=curcellnew(B->nx, B->ny);
    }else{
	assert((*A)->nx==B->nx && (*A)->ny==B->ny);
    }
    for(int i=0; i<B->nx*B->ny; i++){
	curadd(&((*A)->p[i]), beta, B->p[i], alpha, handle);
    }
}
