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
#ifndef AOS_CUDA_CUBLAS_H
#define AOS_CUDA_CUBLAS_H
#include <typeinfo>
#include "types.h"
#include <cusolverDn.h>
/**
 * \file cublas.h
*/
void cuspmul(Real *y, const cusp &A, const Real *x, int ncol, char trans, Real alpha, stream_t &stream);
//int cusvd(curmat &U, curmat &S, curmat &VT, curmat &A);
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

template<typename T>
int cusvd(NumArray<T, Gpu> &U, NumArray<T, Gpu> &S, NumArray<T, Gpu> &Vt, NumArray<T, Gpu> &A){
  cusolverDnHandle_t handle; cusolverDnCreate(&handle);
  cusolverDnParams_t params; cusolverDnCreateParams(&params);
  cusolverStatus_t status;
  S.init(MIN(A.Nx(), A.Ny()), 1);
  U.init(A.Nx(), A.Nx());
  Vt.init(A.Ny(), A.Ny());
  size_t wDevice, wHost;
  status=cusolverDnXgesvd_bufferSize(handle, params, 'A', 'A', A.Nx(), A.Ny(),
    A.datatype(), A(), A.Nx(), S.datatype(), S(), U.datatype(), U(), U.Nx(), Vt.datatype(), Vt(), Vt.Nx(), A.datatype(), &wDevice, &wHost);
  CHECK(status);
  //Allocate work space
  RefP<char, Gpu> tmpDevice(wDevice);
  RefP<char, Cpu> tmpHost(wHost);
  RefP<int, Gpu>tmpans(1);
  int ans;
  status=cusolverDnXgesvd(handle, params, 'A', 'A', A.Nx(), A.Ny(),
    A.datatype(), A(), A.Nx(), S.datatype(), S(), U.datatype(), U(), U.Nx(), Vt.datatype(), Vt(), Vt.Nx(), A.datatype(), tmpDevice(), wDevice, tmpHost(), wHost, tmpans());
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
template<typename T>
int cusvdp(NumArray<T, Gpu> &U, NumArray<T, Gpu> &S, NumArray<T, Gpu> &V, NumArray<T, Gpu> &A){
  cusolverDnHandle_t handle; cusolverDnCreate(&handle);
  cusolverDnParams_t params; cusolverDnCreateParams(&params);
  cusolverEigMode_t jobz=CUSOLVER_EIG_MODE_VECTOR;
  cusolverStatus_t status;
  S.init(MIN(A.Nx(), A.Ny()), 1);
  U.init(A.Nx(), A.Nx());
  V.init(A.Ny(), A.Ny());
  size_t wDevice, wHost;

  status=cusolverDnXgesvdp_bufferSize(handle, params, jobz, 0, A.Nx(), A.Ny(),
    A.datatype(), A(), A.Nx(), S.datatype(), S(), U.datatype(), U(), U.Nx(), V.datatype(), V(), V.Nx(), A.datatype(), &wDevice, &wHost);
  CHECK(status);
  //Allocate work space
  RefP<char, Gpu> tmpDevice(wDevice);
  RefP<char, Cpu> tmpHost(wHost);
  int ans;
  RefP<int, Gpu>tmpans(1);
  double h_err_sigma;
  status=cusolverDnXgesvdp(handle, params, jobz, 0, A.Nx(), A.Ny(),
    A.datatype(), A(), A.Nx(), S.datatype(), S(), U.datatype(), U(), U.Nx(), V.datatype(), V(), V.Nx(), A.datatype(), tmpDevice(), wDevice, tmpHost(), wHost, tmpans(), &h_err_sigma);
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
#endif
#endif
