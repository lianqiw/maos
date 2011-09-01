#ifndef AOS_CUDA_CURMAT_H
#define AOS_CUDA_CURMAT_H

#include "utils.h"

typedef struct{
    float *p;
    int nx;
    int ny;
    int ref;
}curmat;
/**
   curcell: call curcell is created by curcellnew2, p0 contains the allocated
   memory. p contains points inside the memory. curcell can therefore be casted to
   a curmat and do dot/copy operations.
*/
typedef struct{
    curmat **p;
    int nx;
    int ny;
}curcell;

typedef struct{
    int *p;
    int *i;
    float *x;
    int nx;
    int ny;
    int nzmax;
}cusp;
typedef struct{
    cusp **p;
    int nx;
    int ny;
}cuspcell;
typedef struct{
    cuspcell *Mt;
    curcell *U;
    curcell *V;
}cumuv_t;
cuspcell *cuspcellnew(int nx, int ny);
curcell *curcellnew(int nx, int ny);
curcell *curcellnew2(const curcell *in);
void curcellfree(curcell *A);
curmat *curnew(int nx, int ny);
curmat *curnew(int nx, int ny, float *p);
void curfree(curmat *A);
void curwrite(const curmat *A, const char *format, ...);
void curcellwrite(const curcell *A, const char *format, ...);
void curset(curmat *A, float alpha, cudaStream_t stream);
void curshow(curmat *A, cudaStream_t stream);
void curzero(curmat *A, cudaStream_t stream);
void curcp(curmat **out, const curmat *in, cudaStream_t stream);
void curadd(curmat **out, float beta, curmat *in, float alpha, cublasHandle_t handle);
void curscale(curmat *in, float alpha, cudaStream_t stream);
void curmv(curmat **C, float alpha, const curmat *A, const curmat *B, char trans, float beta, cublasHandle_t handle);
void curmm(curmat **C, float alpha, const curmat *A, const curmat *B, char trans[2], float beta, cublasHandle_t handle);
/**
   Compute the inner product as vectors. A.*B
*/
/*
inline void curinn(float *result, const curmat *A, const curmat *B, cublasHandle_t handle){
    cublasSdot(handle, A->nx*A->ny, A->p, 1, B->p, 1, result);
}
inline float curcellinn(const curcell *A, const curcell *B, cublasHandle_t handle){
    float sum[A->nx*A->ny];
    for(int i=0; i<A->nx*A->ny; i++){
	curinn(&sum[i], A->p[i], B->p[i], handle);//notice that result is replaced, not accumulated.
    }
    cudaStream_t stream;
    cublasGetStream(handle, &stream);
    cudaStreamSynchronize(stream);
    float result=0;
    for(int i=0; i<A->nx*A->ny; i++){
	result+=sum[i];
    }
    return result;
}
*/

void curcellzero(curcell *A, cudaStream_t stream);
void curcellcp(curcell **A, const curcell *B, cudaStream_t stream);
void curcelladd(curcell **A, float beta, const curcell *B, float alpha, cublasHandle_t handle);
__global__ void adds_do(float *vec, float *palpha, float beta, int n);
__global__ void add_do(float *restrict a, const float *restrict b, const float *restrict b_sc1, float b_sc2, int n);
void curadd2(curmat **out, const curmat *in, float *alpha, cudaStream_t stream);
void curadd3(curmat **out, float *beta, const curmat *in, cudaStream_t stream);
void curcelladd2(curcell **A, const curcell *B, float* alpha, cudaStream_t stream);
void curcelladd3(curcell **A, float* beta, const curcell *B, cudaStream_t stream);

/**
   Routine that does reduction.
*/
__global__ void inn_do(float *restrict res, const float *a, const float *b, const int n);
__global__ void reduce_do(float *res, const float *a, const int n);
float curinn(const curmat *a, const curmat *b, cudaStream_t stream);
float curcellinn(const curcell *A, const curcell *B, cudaStream_t stream);
void curinn2(float *restrict res, const curmat *a, const curmat *b, cudaStream_t stream);
void curcellinn2(float *restrict res, const curcell *A, const curcell *B, cudaStream_t stream);
void cursum2(float *restrict, const curmat *a, cudaStream_t stream);

#endif
