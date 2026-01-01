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
#ifndef AOS_CUDA_CUBLAS_H
#define AOS_CUDA_CUBLAS_H
#include <typeinfo>
#include "utils.h"
/**
 * \file cublas.h
*/
void cuspmul(curmat &y, const cusp &A, const curmat &x, long ncol, char trans, Real alpha, stream_t &stream);
//int cusvd(curmat &U, curmat &S, curmat &VT, curmat &A);
/**<
 * C=A*diag(x)
 * */
template<typename T>
cublasStatus_t cublasGdgmm(cublasHandle_t handle, cublasSideMode_t mode,
                          int m, int n,
                          const T *A, int lda,
                          const T *x, int incx,
                          T *C, int ldc);

/**
 * C=A*S
 * Sd=diagonal(S)
*/
template<typename T>
void cudgmm(NumArray<T, Gpu>C, NumArray<T, Gpu>A, NumArray<T, Gpu>Sd, stream_t &stream){
    DO(cublasGdgmm(stream.blas(), CUBLAS_SIDE_RIGHT, A.Nx(), A.Ny(), A(), A.Nx(), Sd(), 1, C(), C.Nx()));
}
template<typename T, typename TS>
void cudgmm(NumArray<T, Gpu>C, NumArray<T, Gpu>A, NumArray<TS, Gpu>Sd, stream_t &stream){
	NumArray<T, Gpu>Sd2;
	Copy(Sd2, Sd, stream);
	DO(cublasGdgmm(stream.blas(), CUBLAS_SIDE_RIGHT, A.Nx(), A.Ny(), A(), A.Nx(), Sd2(), 1, C(), C.Nx()));
}
/**
 * C=alpha*op(A)*op(B)+beta*C
 * */
template<typename T>
cublasStatus_t cublasGgemm(cublasHandle_t handle,
                           cublasOperation_t transa, cublasOperation_t transb,
                           int m, int n, int k,
                           const T *alpha,
                           const T *A, int lda,
                           const T *B, int ldb,
                           const T *beta,
                           T *C, int ldc);

/**
   Computes C = beta * C + alpha * op(A) * B ;
*/
template<typename T>
void cugemm(NumArray<T, Gpu> &C, T beta, const NumArray<T, Gpu> &A, const NumArray<T, Gpu> &B, const char trans[2], T alpha, stream_t &stream){
    int m, n, k, k2;
    cublasOperation_t transa, transb;
    if(trans[0]=='t'){
        m=A.Ny();
        k=A.Nx();
        transa=CUBLAS_OP_T;
    } else{
        m=A.Nx();
        k=A.Ny();
        transa=CUBLAS_OP_N;
    }
    if(trans[1]=='t'){
        n=B.Nx();
        k2=B.Ny();
        transb=CUBLAS_OP_T;
    } else{
        n=B.Ny();
        k2=B.Nx();
        transb=CUBLAS_OP_N;
    }
    if(!C){
        C.init(m, n);
    } else{
        assert(C.Nx()==m&&C.Ny()==n);
    }
    if(k!=k2) error("Matrix mismatch\n");
    DO(cublasGgemm(stream.blas(), transa, transb, m, n, k,
                   &alpha, A(), A.Nx(), B(), B.Nx(), &beta, C(), C.Nx()));
}
/**
 * SVD decomposition of a dense matrix using QR algorithm. Deprecated
 * A = U * S * Vt
 * */
#if CUDA_VERSION >= 10000
#define CHECK(status)\
	switch(status){\
		case CUSOLVER_STATUS_SUCCESS: break;\
		case CUSOLVER_STATUS_NOT_INITIALIZED:\
			error("Library was not initialized\n"); break;\
		case CUSOLVER_STATUS_INVALID_VALUE:\
			error("Invalid parameters were passed\n");	break;\
		case CUSOLVER_STATUS_INTERNAL_ERROR:\
			error("An internal operation failed.\n"); break;\
		default: \
			error("Unknown error happened: %d\n", (int) status); break;\
	}\

template<typename T, typename R>
int cusvd(NumArray<T, Gpu> &U, NumArray<R, Gpu> &S, NumArray<T, Gpu> &Vt, NumArray<T, Gpu> &A, stream_t &stream){
    cusolverDnParams_t params; cusolverDnCreateParams(&params);
    cusolverStatus_t status;
    S.init(MIN(A.Nx(), A.Ny()), 1);
    U.init(A.Nx(), A.Nx());
    Vt.init(A.Ny(), A.Ny());
    size_t wDevice, wHost;
    status=cusolverDnXgesvd_bufferSize(stream.dn(), params, 'A', 'A', A.Nx(), A.Ny(),
      A.dtype(), A(), A.Nx(), S.dtype(), S(), U.dtype(), U(), U.Nx(), Vt.dtype(), Vt(), Vt.Nx(), A.dtype(), &wDevice, &wHost);
    CHECK(status);
    //Allocate work space
    RefP<char, Gpu> tmpDevice(wDevice);
    RefP<char> tmpHost(wHost);
    RefP<int, Gpu>tmpans(1);
    int ans;
    status=cusolverDnXgesvd(stream.dn(), params, 'A', 'A', A.Nx(), A.Ny(),
      A.dtype(), A(), A.Nx(), S.dtype(), S(), U.dtype(), U(), U.Nx(), Vt.dtype(), Vt(), Vt.Nx(), A.dtype(), tmpDevice(), wDevice, tmpHost(), wHost, tmpans());
    CHECK(status);
    DO(cudaMemcpy(&ans, tmpans(), sizeof(int), D2H));

    if(ans){
        if(ans<0){
            warning("the %d-th parameter is wrong", -ans);
        } else{
            warning("%d superdiagonals of an intermediate bidiagonal form did not converge to zero", ans);
        }
    }
    return ans;
}
/**
 * SVD decomposition of a dense matrix using polar decomposition. It is much
 * faster than QR algorithm.
 *
 * A = U * S * Vt
 *
 * cusolverDnXgesvdp combines polar decomposition in [14] and cusolverDnXsyevd
 * to compute SVD. It is much faster than cusolverDnXgesvd which is based on QR
 * algorithm. However polar decomposition in [14] may not deliver a full unitary
 * matrix when the matrix A has a singular value close to zero. To workaround
 * the issue when the singular value is close to zero, we add a small
 * perturbation so polar decomposition can deliver the correct result. The
 * consequence is inaccurate singular values shifted by this perturbation. The
 * output parameter h_err_sigma is the magnitude of this perturbation. In other
 * words, h_err_sigma shows the accuracy of SVD.
 * */
template<typename T, typename TS>
int cusvdp(NumArray<T, Gpu> &U, NumArray<TS, Gpu> &S, NumArray<T, Gpu> &V, NumArray<T, Gpu> &A, stream_t &stream){
    cusolverDnParams_t params; cusolverDnCreateParams(&params);
    cusolverEigMode_t jobz=CUSOLVER_EIG_MODE_VECTOR;
    cusolverStatus_t status;
    S.init(MIN(A.Nx(), A.Ny()), 1);
    U.init(A.Nx(), A.Nx());
    V.init(A.Ny(), A.Ny());
    size_t wDevice, wHost;

    status=cusolverDnXgesvdp_bufferSize(stream.dn(), params, jobz, 0, A.Nx(), A.Ny(),
      A.dtype(), A(), A.Nx(), S.dtype(), S(), U.dtype(), U(), U.Nx(), V.dtype(), V(), V.Nx(), A.dtype(), &wDevice, &wHost);
    CHECK(status);
    //Allocate work space
    RefP<char, Gpu> tmpDevice(wDevice);
    RefP<char> tmpHost(wHost);
    int ans;
    RefP<int, Gpu>tmpans(1);
    double h_err_sigma;
    status=cusolverDnXgesvdp(stream.dn(), params, jobz, 0, A.Nx(), A.Ny(),
      A.dtype(), A(), A.Nx(), S.dtype(), S(), U.dtype(), U(), U.Nx(), V.dtype(), V(), V.Nx(), A.dtype(), tmpDevice(), wDevice, tmpHost(), wHost, tmpans(), &h_err_sigma);
    CHECK(status);
    DO(cudaMemcpy(&ans, tmpans(), sizeof(int), D2H));
    if(ans){
        if(ans<0){
            warning("the %d-th parameter is wrong", -ans);
        } else{
            warning("%d superdiagonals of an intermediate bidiagonal form did not converge to zero", ans);
        }
    }
    if(h_err_sigma){
        warning("small singular values are perturbed by %g amount.\n", h_err_sigma);
    }
    return ans;
}
/**
 * Use SVD decomposition to compute power (usually inverse) of A. 
 * A=U*S*Vt --> A^-1 = Vt'*S^-1*U'
 * */
template<typename T, typename ST>
int cusvd_pow(NumArray<T, Gpu> &A, ST power, ST thres1, ST thres2, stream_t &stream){
    NumArray<T, Gpu> U, Vt;
	NumArray<ST, Gpu> S;
    cusvd(U, S, Vt, A, stream);
	int nstop=S.N();
	singular_gap<<<REDUCE(S.N()), DIM_REDUCE*sizeof(int), stream>>>
		(S(), S.N(), nstop, thres1, thres2);
	//if(nstop<S.N()){
	//CUDA_SYNC_STREAM;
	//dbg("gpu: crop at %d out of %ld singular values, thres=%g, %g\n", nstop, S.N(), thres1, thres2);
	//}
	singular_pow<<<DIM(S.N(), 256), 0, stream>>>(S(), S.N(), power, nstop);
    //U=U*S
   	cudgmm(U, U, S, stream);
	T alpha, beta;
	ST alpha2=0, beta2=1;
	type_convert(alpha, alpha2);
	type_convert(beta, beta2);
    cugemm(A, alpha, Vt, U, "tt", beta, stream);
    return 0;
}
#endif
#endif
