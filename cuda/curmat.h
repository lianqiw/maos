/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_CUDA_CURMAT_H
#define AOS_CUDA_CURMAT_H
#include "types.h"
#include "kernel.h"
#include "cumat.h"
#define curcellcp    cucellcp<Real>

void curset(curmat &A, Real alpha, cudaStream_t stream);
void curshow(curmat &A, cudaStream_t stream);
void curcp(curmat &out, const curmat *in);
void curcp(curmat &out, const curmat &in, cudaStream_t stream);
void curadd(curmat &out,Real alpha,const curmat &in,Real beta,cudaStream_t stream);
void curaddcabs2(curmat &out, Real alpha, const cucmat &in, Real beta, cudaStream_t stream);
void curscale(curmat &in, Real alpha, cudaStream_t stream);
void curmv(Real *c, Real alpha, const curmat &A, const Real *b, 
	   char trans, Real beta, cublasHandle_t handle);
void curmm(curmat &C, Real alpha, const curmat &A, const curmat &B, 
	   const char trans[2], Real beta, cublasHandle_t handle);
void curcellmm(curcell &C0, Real alpha, const curcell &A, const curcell &B, 
	       const char trans[2], const double beta, cublasHandle_t handle);
void curcelladd(curcell &A, Real beta, const curcell &B, Real alpha, cudaStream_t stream);
__global__ void add_do(Real *vec, Real *palpha, Real beta, int n);
__global__ void add_do(Real *restrict a, const Real * b, const Real *restrict b_sc1, Real b_sc2, int n);
void curadd(curmat &out, const curmat &in, Real *alpha, Real alpha2, cudaStream_t stream);
void curadd(curmat &out, Real *beta, const curmat &in, cudaStream_t stream);
void curcelladd(curcell &A, const curcell &B, Real* alpha, Real alpha2, cudaStream_t stream);
void curcelladd(curcell &A, Real* beta, const curcell &B, cudaStream_t stream);
void curadd(curmat &A, Real beta, cudaStream_t stream);

/**
   Routine that does reduction.
*/
Real curinn(const curmat &a, const curmat &b, cudaStream_t stream);
void cursum2(Real *restrict, const curmat &a, cudaStream_t stream);
Real cursum(const curmat &a, cudaStream_t stream);
void curcellscale(curcell &A, Real alpha, cudaStream_t stream);
Real curmax(const curmat &a, cudaStream_t stream);
Real curmaxabs(const curmat &a, cudaStream_t stream);
Real curcellmax(const curcell &a, cudaStream_t stream);
Real curcellmaxabs(const curcell &a, cudaStream_t stream);
/**
   Add tip/tilt to OPD
*/
inline void curaddptt(curmat &opd, Real (*loc)[2], Real pis, Real tx, Real ty, cudaStream_t stream){
    add_ptt_do<<<DIM(opd.N(), 256), 0, stream>>>(opd.P(), loc, opd.N(), pis, tx, ty);
}
inline void curaddptt(curmat &opd, Real (*loc)[2], Real *ptt, Real pis, Real tx, Real ty,  cudaStream_t stream){
    add_ptt_do<<<DIM(opd.N(), 256), 0, stream>>>(opd.P(), loc, opd.N(), ptt, pis, tx, ty);
}

#endif
