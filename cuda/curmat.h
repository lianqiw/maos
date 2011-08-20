#ifndef AOS_CUDA_CURMAT_H
#define AOS_CUDA_CURMAT_H
#include "utils.h"
typedef struct{
    float *p;
    int nx;
    int ny;
}curmat;
typedef struct{
    curmat **p;
    int nx;
    int ny;
}curcell;

curcell *curcellnew(int nx, int ny);
curcell *curcellnew2(const curcell *in);
void curcellfree(curcell *A);
curmat *curnew(int nx, int ny);
void curfree(curmat *A);
void curwrite(const curmat *A, const char *format, ...);
void curcellwrite(const curcell *A, const char *format, ...);

void curzero(curmat *A, cudaStream_t stream);
void curcp(curmat **out, const curmat *in, cudaStream_t stream);
void curadd(curmat **out, float beta, curmat *in, float alpha, cublasHandle_t handle);
void curscale(curmat *in, float alpha, cudaStream_t stream);
void curmv(curmat **C, float alpha, const curmat *A, const curmat *B, char trans, float beta, cublasHandle_t handle);
void curmm(curmat **C, float alpha, const curmat *A, const curmat *B, char trans[2], float beta, cublasHandle_t handle);
/**
   Compute the inner product as vectors. A.*B
*/
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
void curcellzero(curcell *A, cudaStream_t stream);
void curcellcp(curcell **A, const curcell *B, cudaStream_t stream);
void curcelladd(curcell **A, float beta, const curcell *B, float alpha, cublasHandle_t handle);
#endif
