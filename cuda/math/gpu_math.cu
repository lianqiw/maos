/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>

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

/**
 * \file gpu_math.cu
 * 
 * Wraps cuda routines for CPU data type.
 * */
#include "cublas.h"
#include "utils.h"
#include "gpu_math.h"
#if CUDA_VERSION>10000
/**
 * Compute svd
 * */
void gpu_dsvd(dmat **U_, dmat **S_, dmat **Vt_, const dmat *A_){
  int cuda_dedup_save=cuda_dedup; cuda_dedup=0;
  //NumArray<real, Gpu> U, S, Vt, A;//match cpu precision (double is slow)
  NumArray<Real, Gpu> U, S, Vt, A;//match GPU precision
  stream_t stream;
  cp2gpu(A, A_, stream);
  cusvd(U, S, Vt, A, stream); 
  cp2cpu(U_, U, stream);
  cp2cpu(S_, S, stream);
  //dmat *V_=NULL, *Vt2=NULL;
	//cp2cpu(&V_, V);
  //Vt2=dtrans(V_); dfree(V_);
  cp2cpu(Vt_, Vt, stream);
  //dfree(Vt2);
  cuda_dedup=cuda_dedup_save;
}
/**
 * Invert matrix (pow=-1) or raise power of a matrix with svd.
 * */
void gpu_dsvd_pow(dmat *A_, real pow, real thres){
  if(thres<0){
    error("negative thres is not supported\n");
  }
  int cuda_dedup_save=cuda_dedup; cuda_dedup=0;
  NumArray<Real, Gpu> A;//match GPU precision
  stream_t stream;
  cp2gpu(A, A_, stream);
  cusvd_pow(A, (Real)pow, (Real)thres, stream);
  cp2cpu(&A_, A, stream);
  cuda_dedup=cuda_dedup_save;
}
/**
 * matrix multplication in gpu
 */
void gpu_dgemm(dmat **C_, const real beta, const dmat *A_, const dmat *B_, const char trans[2], const real alpha){
  int cuda_dedup_save=cuda_dedup; cuda_dedup=0;
  NumArray<Real, Gpu>A,B,C;
  stream_t stream;
  cp2gpu(A, A_, stream);
  cp2gpu(B, B_, stream);
  if(*C_) cp2gpu(C, *C_, stream); 
  cugemm(C, (Real)beta, A, B, trans, (Real)alpha, stream);
  cp2cpu(C_,C,stream);
  cuda_dedup=cuda_dedup_save;
}
void gpu_dgemm_test(){
	dmat *A=dnew(2000, 2000);
	dmat *B=dref(A);
	rand_t rstat; seed_rand(&rstat, 1);
	drandn(A, 10, &rstat);
	daddI(A, 1);
	dmat *C=NULL;
	dmm(&C, 0, A, B, "nn", 1);
	dmat *C2=NULL;
	gpu_dgemm(&C2, 0, A, B, "nn", 1);
	real diffu=ddiff(C, C2);
	real CS=dsumabs(C);
	real C2S=dsumabs(C2);
	dmm(&C, 1, A, B, "nn", 1);
	gpu_dgemm(&C2, 1, A, B, "nn", 1);
	real diffu2=ddiff(C, C2);
	dfree(A); dfree(B); dfree(C); dfree(C2);
	info("dmm and gpu_dgemm diff are %g %g. sum are %g, %g\n", diffu, diffu2, CS, C2S);
}
void gpu_dsvd_test(){//test
	//gpu_dsvd is faster than dsvd for matrix larger than 500x500 (test is on cassiopeia)
	dmat *A=dnew(2000, 2000);
	rand_t rstat; seed_rand(&rstat, 1);
	drandn(A, 10, &rstat);
	daddI(A, 1);
	dmat *U=0, *Vt=0, *S=0;
	dmat *U2=0, *Vt2=0, *S2=0;
	TIC;
	tic;
	dsvd(&U, &S, &Vt, A);
	toc("dsvd");tic;
	gpu_dsvd(&U2, &S2, &Vt2, A);
	toc("gpu_dsvd");
	real diffu=ddiff(U, U2);
	real diffs=ddiff(S, S2);
	real diffv=ddiff(Vt, Vt2);
	info("dsvd and gpu_svd diff are %g %g %g\n", diffu, diffs, diffv);
	if(diffu>0.01 || diffs>0.01 || diffv>0.01){
		writebin(U, "U");
		writebin(S, "S");
		writebin(Vt, "Vt");
		writebin(U2, "U2");
		writebin(S2, "S2");
		writebin(Vt2, "Vt2");
		writebin(A, "A");
	}
}
#endif
