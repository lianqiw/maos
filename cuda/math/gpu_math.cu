/*
  Copyright 2009-2026 Lianqi Wang
  
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
#define GPU_SVD(X,R,TG)\
void gpu_##X##svd(X##mat **U_, R##mat **S_, X##mat **Vt_, const X##mat *A_){\
	int cuda_dedup_save=cuda_dedup; cuda_dedup=0;\
	NumArray<TG, Gpu> U, Vt, A;\
	NumArray<Real, Gpu> S;\
	stream_t stream;\
	cp2gpu(A, A_, stream);\
	cusvd(U, S, Vt, A, stream);\
	cp2cpu(U_, U, stream);\
	cp2cpu(S_, S, stream);\
	cp2cpu(Vt_, Vt, stream);\
	cuda_dedup=cuda_dedup_save;\
}
GPU_SVD(d, d, Real);
GPU_SVD(c, d, Comp);
GPU_SVD(s, s, Real);
GPU_SVD(z, s, Comp);

/**
 * Invert matrix (pow=-1) or raise power of a matrix with svd.
 * */
#define GPU_SVD_POW(X,TG, RG)\
void gpu_##X##svd_pow(X##mat *A_, real pow, real thres1, real thres2){\
	int cuda_dedup_save=cuda_dedup; cuda_dedup=0;\
	NumArray<TG, Gpu> A;\
	stream_t stream;\
	cp2gpu(A, A_, stream);\
	cusvd_pow(A, (RG)pow, (RG)thres1, (RG)thres2, stream);\
	cp2cpu(&A_, A, stream);\
	cuda_dedup=cuda_dedup_save;\
}
GPU_SVD_POW(d, Real, Real)
GPU_SVD_POW(c, Comp, Real)
GPU_SVD_POW(s, Real, Real)
GPU_SVD_POW(z, Comp, Real)

/**
 * matrix multplication in gpu
 */
#define GPU_GEMM(X,TC,TG)\
void gpu_##X##gemm(X##mat **C_, TC beta, const X##mat *A_, const X##mat *B_, const char trans[2], TC alpha){\
	int cuda_dedup_save=cuda_dedup; cuda_dedup=0;\
	NumArray<TG, Gpu>A, B, C;\
	stream_t stream;\
	cp2gpu(A, A_, stream);\
	cp2gpu(B, B_, stream);\
	if(*C_) cp2gpu(C, *C_, stream);\
	TG beta2, alpha2;\
	type_convert(&beta2, &beta, 1);\
	type_convert(&alpha2, &alpha, 1);\
	cugemm(C, beta2, A, B, trans, alpha2, stream);\
	cp2cpu(C_, C, stream);\
	cuda_dedup=cuda_dedup_save;\
}
GPU_GEMM(d, double, Real);
GPU_GEMM(c, dcomplex, Comp);
GPU_GEMM(s, float, Real);
GPU_GEMM(z, fcomplex, Comp);

#define GPU_FFT(X, TG, FFT_T)\
void gpu_##X##fft2(X##mat *A_, int direction){\
	int cuda_dedup_save=cuda_dedup; cuda_dedup=0;\
	NumArray<TG, Gpu>A;\
	stream_t stream;\
	cp2gpu(A, A_, stream);\
	cufftHandle plan;\
	DO(cufftPlan2d(&plan, A.Nx(), A.Ny(), FFT_T));\
	CUFFT(plan, A(), direction==1?CUFFT_FORWARD:CUFFT_INVERSE);\
	cp2cpu(&A_, A, stream);\
	cufftDestroy(plan);\
	cuda_dedup=cuda_dedup_save;\
}
GPU_FFT(c, Comp, CUFFT_C2C);
GPU_FFT(z, Comp, CUFFT_C2C);

void gpu_ext_assign(){
#define CPU_ASSIGN_BLAS(X,R,T)\
	extern void (*X##svd_ext)(X##mat **U, R##mat **S, X##mat **Vt, const X##mat *A);\
	extern void (*X##svd_pow_ext)(X##mat *A_, real pow, real thres1, real thres2);\
	extern void (*X##gemm_ext)(X##mat **out, T beta, const X##mat *A, const X##mat *B, const char trans[2], T alpha);\
	X##svd_ext=gpu_##X##svd;\
	X##svd_pow_ext=gpu_##X##svd_pow;\
	X##gemm_ext=gpu_##X##gemm;

	dbg("Using GPU for fft2, svd, svd_pow and gemm of large arrays.\n");
	CPU_ASSIGN_BLAS(d, d, double);
	CPU_ASSIGN_BLAS(c, d, dcomplex);
	CPU_ASSIGN_BLAS(s, s, float);
	CPU_ASSIGN_BLAS(z, s, fcomplex);

#define CPU_ASSIGN_FFT(X)\
	extern void (*X##fft2_ext)(X##mat *A, int direction);\
	X##fft2_ext=gpu_##X##fft2
	CPU_ASSIGN_FFT(c);
	CPU_ASSIGN_FFT(z);
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
