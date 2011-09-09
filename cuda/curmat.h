#ifndef AOS_CUDA_CURMAT_H
#define AOS_CUDA_CURMAT_H

#include "utils.h"
#include "types.h"
cuspcell *cuspcellnew(int nx, int ny);
curcell *curcellnew(int nx, int ny);
curcell *curcellnew2(const curcell *in);
void curcellfree(curcell *A);
curmat *curnew(int nx, int ny);
curmat *curnew(int nx, int ny, cudaStream_t stream);
void curfree(curmat *A);
void curwrite(const curmat *A, const char *format, ...);
void curcellwrite(const curcell *A, const char *format, ...);
void curset(curmat *A, float alpha, cudaStream_t stream);
void curshow(curmat *A, cudaStream_t stream);
void curzero(curmat *A, cudaStream_t stream);
void curcp(curmat **out, const curmat *in, cudaStream_t stream);
void curadd(curmat **out,float alpha,curmat *in,float beta,cudaStream_t stream);
void curaddcabs2(curmat **out, float alpha, cucmat *in, float beta, cudaStream_t stream);
void curscale(curmat *in, float alpha, cudaStream_t stream);
void curmv(curmat **C, float alpha, const curmat *A, const curmat *B, char trans, float beta, cublasHandle_t handle);
void curmm(curmat **C, float alpha, const curmat *A, const curmat *B, char trans[2], float beta, cublasHandle_t handle);

void curcellzero(curcell *A, cudaStream_t stream);
void curcellcp(curcell **A, const curcell *B, cudaStream_t stream);
void curcelladd(curcell **A, float beta, const curcell *B, float alpha, cudaStream_t stream);
__global__ void adds_do(float *vec, float *palpha, float beta, int n);
__global__ void add2_do(float *restrict a, const float *restrict b, const float *restrict b_sc1, float b_sc2, int n);
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

__global__ void addptt_do(float *restrict opd, float (*restrict loc)[2], int n, float tx, float ty);
/**
   Add tip/tilt to OPD
*/
inline void curaddptt(curmat *opd, float (*loc)[2], float tx, float ty, cudaStream_t stream){
    addptt_do<<<DIM(opd->nx*opd->ny, 256), 0, stream>>>(opd->p, loc, opd->nx*opd->ny, tx, ty);
}
#endif
