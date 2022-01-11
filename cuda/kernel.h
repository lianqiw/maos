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

#ifndef AOS_CUDA_KERNEL_H
#define AOS_CUDA_KERNEL_H
#include "common.h"
/*
__device__ static inline float atomicAdd(float* address, float val)
{
	float old = *address;
	float assumed;
	do {
	assumed = old;
	old = __int_as_float( atomicCAS((unsigned int*)address,
					__float_as_int(assumed),
					__float_as_int(val + assumed)));
	} while (assumed != old);
	return old;
}
#endif
*/
#if __CUDA_ARCH__ < 600  && CUDA_DOUBLE
__device__ static inline double myAtomicAdd(double* address, double val){
	double old=*address;
	double assumed;
	do{
		assumed=old;
		old=__longlong_as_double(atomicCAS((unsigned long long*)address,
			__double_as_longlong(assumed),
			__double_as_longlong(val+assumed)));
	} while(assumed!=old);
	return old;
}
#define atomicAdd myAtomicAdd
#endif
__device__ static inline float atomicMax(float* address, float val){
	float old=*address;
	float assumed;
	do{
		assumed=old;
		old=__int_as_float(atomicCAS((unsigned int*)address,
			__float_as_int(assumed),
			__float_as_int(val>assumed?val:assumed)));
	} while(assumed!=old);
	return old;
}
__device__ static inline double atomicMax(double* address, double val){
	double old=*address;
	double assumed;
	do{
		assumed=old;
		old=__longlong_as_double(atomicCAS((unsigned long long*)address,
			__double_as_longlong(assumed),
			__double_as_longlong(val>assumed?val:assumed)));
	} while(assumed!=old);
	return old;
}
__device__ static inline Real CABS2(Comp r){
	const Real a=Z(cuCreal)(r);
	const Real b=Z(cuCimag)(r);
	return a*a+b*b;
}

__global__ void set_do(Real* a, Real alpha, int n);
__global__ void scale_do(Real* restrict in, int n, Real alpha);
__global__ void scale_do(Comp* restrict in, int n, Real alpha);
__global__ void add_ptt_do(Real* restrict opd, Real(*restrict loc)[2], int n, Real pis, Real tx, Real ty);
__global__ void add_ptt_do(Real* restrict opd, Real(*restrict loc)[2], int n, Real* ptt, Real pis, Real tx, Real ty);
__global__ void add_focus_do(Real* restrict opd, Real(*restrict loc)[2], int n, Real focus);
__global__ void add_ngsmod_do(Real* restrict opd, Real(*restrict loc)[2], int n,
	Real ttx, Real tty,
	Real ps1, Real ps2, Real ps3,
	Real astig1, Real astig2, Real focus,
	Real thetax, Real thetay, Real scale, Real ht, Real alpha);

__global__ void add_do(Real* vec, Real beta, int n);
__global__ void addcabs2_do(Real* restrict a, Real alpha, const Comp* restrict b, Real beta, int n);
__global__ void add_do(Real* restrict a, Real* alpha1, Real alpha2, const Real* restrict b, int n);
__global__ void add_do(Real* restrict a, const Real* restrict b, Real* beta1, Real beta2, int n);
__global__ void add_do(Real* restrict a, Real* alpha1, Real alpha2,
	const Real* restrict b, Real* beta1, Real beta2, int n);

__global__ void max_do(Real* restrict res, const Real* a, const int n);
__global__ void maxabs_do(Real* restrict res, const Real* a, const int n);
__global__ void sum_do(Real* restrict res, const Real* a, const int n);
__global__ void sum2_do(Real* restrict res, const Real* a, const int n);
__global__ void inn_do(Real* res_add, const Real* a, const Real* b, const int n);

static inline void inn_wrap(Real* res_add, const Real* a, const Real* b, const int n, cudaStream_t stream){
	inn_do<<<REDUCE(n), DIM_REDUCE*sizeof(Real), stream>>>(res_add, a, b, n);
}
static inline void sum_wrap(Real* res, const Real* a, const int n, cudaStream_t stream){
	sum_do<<<REDUCE(n), DIM_REDUCE*sizeof(Real), stream>>>(res, a, n);
}
static inline void sum2_wrap(Real* res, const Real* a, const int n, cudaStream_t stream){
	sum2_do<<<REDUCE(n), 0, stream>>>(res, a, n);
}
static inline void max_wrap(Real* res, const Real* a, const int n, cudaStream_t stream){
	max_do<<<REDUCE(n), DIM_REDUCE*sizeof(Real), stream>>>(res, a, n);
}
static inline void maxabs_wrap(Real* res, const Real* a, const int n, cudaStream_t stream){
	maxabs_do<<<REDUCE(n), DIM_REDUCE*sizeof(Real), stream>>>(res, a, n);
}
__global__ void embed_do(Comp* out, Real* in, int nx);
__global__ void extract_do(Real* out, Comp* in, int nx);
__global__ void perm_f_do(Comp* restrict out, const Comp* restrict in, int* restrict perm, int nx);
__global__ void perm_i_do(Comp* restrict out, const Comp* restrict in, int* restrict perm, int nx);
__global__ void perm_f_do(Real* restrict out, const Real* restrict in, int* restrict perm, int nx);
__global__ void perm_i_do(Real* restrict out, const Real* restrict in, int* restrict perm, int nx);
__global__ void embed_wvf_do(Comp* restrict wvf,
	const Real* restrict opd, const Real* restrict amp,
	const int* embed, const int nloc, const Real wvl);
__global__ void corner2center_do(Comp* restrict out, int noutx, int nouty,
	const Comp* restrict in, int ninx, int niny);
__global__ void corner2center_abs2_do(Real* restrict out, int noutx, int nouty,
	const Comp* restrict in, int ninx, int niny);
__global__ void corner2center_abs2_atomic_do(Real* restrict out, int noutx, int nouty,
	const Comp* restrict in, int ninx, int niny);
__global__ void fftshift_do(Comp* wvf, const int nx, const int ny);
__global__ void add_tilt_do(Real* opd, int nx, int ny, Real ox, Real oy, Real dx, Real ttx, Real tty);
template<typename D, typename F>
__global__ void cwm_do(D* dest, F* from, long n){
	for(int i=threadIdx.x+blockIdx.x*blockDim.x; i<n; i+=blockDim.x*gridDim.x){
		dest[i]*=from[i];
	}
}

__global__ void cwm_do(Comp* dest, Comp* from, int lda, int ldb, int nx, int ny);
__global__ void cwm_do(Comp* dest, Comp* from1, Comp* from2, int lda, int ldb, int nx, int ny);
__global__ void unwrap_phase_do(Comp* wvf, Real* opd, int* embed, int n, Real wvl);
__global__ void mvm_do(const Real* restrict mvm, Real* restrict a, const Real* restrict g, int nact, int ng);
__global__ void multimv_do(const Real* restrict mvm, Real* restrict a, const Real* restrict g, int nact, int ng);
#endif
