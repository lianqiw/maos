/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
template <typename T>
__global__ void set_do(T *a, T alpha, int n){
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		a[i]=alpha;
	}
}
template <typename T, typename S>
__global__ void scale_do(T *restrict in, int n, S alpha){
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		in[i]*=alpha;
	}
}
/**
   add a vector to another, scaled by alpha and beta. all in device memory.
   a=a*(*alpha1)*alpha2+b*(*beta1)*beta2;
*/
template <typename T>
__global__ void add_do(T *restrict a, T *alpha1, T alpha2,
	const T *restrict b, T *beta1, T beta2, int n){
	T alpha=alpha1?(*alpha1*alpha2):alpha2;
	T beta=beta1?(*beta1*beta2):beta2;
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		a[i]=a[i]*alpha+b[i]*beta;
	}
}
template <typename T>
__global__ void add_do(T *restrict a, const T *restrict b, T *beta1, T beta2, int n){
	T beta=beta1?(*beta1*beta2):beta2;
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		a[i]+=b[i]*beta;
	}
}
template <typename T>
__global__ void add_do(T *restrict a, T *alpha1, T alpha2, const T *restrict b, int n){
	T alpha=alpha1?(*alpha1*alpha2):alpha2;
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		a[i]=a[i]*alpha+b[i];
	}
}

/**
Add a scalar to a vector.
*/
template <typename T>
__global__ void add_do(T *vec, T beta, int n){
	const int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=step){
		vec[i]+=beta;
	}
}
//__global__ void set_do(Real* a, Real alpha, int n);
//__global__ void scale_do(Real* restrict in, int n, Real alpha);
//__global__ void scale_do(Comp* restrict in, int n, Real alpha);
__global__ void add_ptt_do(Real* restrict opd, Real(*restrict loc)[2], int n, Real pis, Real tx, Real ty);
__global__ void add_ptt_do(Real* restrict opd, Real(*restrict loc)[2], int n, Real* ptt, Real pis, Real tx, Real ty);
__global__ void add_focus_do(Real* restrict opd, Real(*restrict loc)[2], int n, Real focus);
__global__ void add_ngsmod_do(Real* restrict opd, Real(*restrict loc)[2], int n,
	Real ttx, Real tty,
	Real ps1, Real ps2, Real ps3,
	Real astig1, Real astig2, Real focus,
	Real thetax, Real thetay, Real scale, Real ht, Real alpha);


__global__ void addcabs2_do(Real *restrict a, Real alpha, const Comp *restrict b, Real beta, int n);
__global__ void addcabs2_do(Real *restrict a, const Comp *restrict b, Real beta, int n);
__global__ void max_do(Real* restrict res, const Real* a, const int n);
__global__ void maxabs_do(Real* restrict res, const Real* a, const int n);
__global__ void sum_do(Real* restrict res, const Real* a, const int n);
__global__ void sum2_do(Real* restrict res, const Real* a, const int n);
__global__ void inn_do(Real* res_add, const Real* a, const Real* b, const int n);

static inline void inn_wrap(Real* res_add, const Real* a, const Real* b, const int n, cudaStream_t stream){
	inn_do<<<REDUCE(n), DIM_REDUCE*sizeof(double), stream>>>(res_add, a, b, n);
}
static inline void inn_multi_wrap(Real* res_add, const Real* a, const int acol, const Real* b, const int n, cudaStream_t stream){
	for(int icol=0; icol<acol; icol++){
		inn_do<<<REDUCE(n), DIM_REDUCE*sizeof(double), stream>>>(res_add+icol, a+icol*n, b, n);
	}
}
static inline void sum_wrap(Real* res, const Real* a, const int n, cudaStream_t stream){
	sum_do<<<REDUCE(n), DIM_REDUCE*sizeof(double), stream>>>(res, a, n);
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
__global__ void cwm_do(D *dest, const F *from, long n){
	for(long i=threadIdx.x+blockIdx.x*blockDim.x; i<n; i+=blockDim.x*gridDim.x){
		dest[i]*=from[i];
	}
}
template<typename D, typename F>
__global__ void copy_do(D *dest, const F *from, long n){
	for(long i=threadIdx.x+blockIdx.x*blockDim.x; i<n; i+=blockDim.x*gridDim.x){
		//dest[i]=(D)from[i];
		type_convert(dest[i], from[i]);
	}
}

__global__ void cwm_do(Comp* dest, Comp* from, int lda, int ldb, int nx, int ny);
__global__ void cwm_do(Comp* dest, Comp* from1, Comp* from2, int lda, int ldb, int nx, int ny);
__global__ void unwrap_phase_do(Comp* wvf, Real* opd, int* embed, int n, Real wvl);
__global__ void mvm_do(const Real* restrict mvm, Real* restrict a, const Real* restrict g, int nact, int ng);
__global__ void multimv_do(const Real* restrict mvm, Real* restrict a, const Real* restrict g, int nact, int ng);

template<typename T>
__global__ void singular_pow(T *restrict S, int n, T power, int& nstop){
	/*if(threadIdx.x==0 && blockIdx.x==0){
		printf("gpu: crop at %d out of %d\n", nstop, n);
	}*/
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=blockDim.x*gridDim.x){
		if(i<nstop){
			S[i]=pow(S[i], power);
		} else{
			S[i]=0;
		}
	}
}
/**
	Determine the index to stop compute*/
template<typename T>
__global__ void singular_gap(T *restrict S, int n, int &nstop, T thres1, T thres2){
	extern __shared__ int nb[];//index 
	nb[threadIdx.x]=n;
	thres1*=S[0];
	//First, find the index where S[i]>S[i-1]*thres2 is no longer true
	//First, use shared memory to process for each element
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<n; i+=blockDim.x*gridDim.x){
		if(i>0 && (S[i]<thres1 || S[i]<S[i-1]*thres2)){
			if(nb[threadIdx.x]>i){
				nb[threadIdx.x]=i;
			}
			//printf("nb[%d]=%d, i=%ld\n", threadIdx.x, nb[threadIdx.x], i);
			break;
		}
	}
	//Then reduce all threads within each block
	for(int step=(blockDim.x>>1);step>0;step>>=1){
		__syncthreads();
		if(threadIdx.x<step){
			if(nb[threadIdx.x]>nb[threadIdx.x+step]){
				nb[threadIdx.x]=nb[threadIdx.x+step];
				//printf("nb[%d]=%d\n", threadIdx.x, nb[threadIdx.x]);
			}
		}
	}
	if(threadIdx.x==0){
		atomicMin(&nstop, nb[threadIdx.x]);
	}
}

/*Transpose a matrix in naive way. Faster way is to use shared memory and handle
  a block each time.*/
template <typename T>
__global__ void transpose(T *restrict out, const T *restrict in, int nx, int ny){
	const int stepx=blockDim.x*gridDim.x;
	const int stepy=blockDim.y*gridDim.y;
	const int ix0=threadIdx.x+blockDim.x*blockIdx.x;
	const int iy0=threadIdx.y+blockDim.y*blockIdx.y;
	for(int iy=iy0; iy<ny; iy+=stepy){
		for(int ix=ix0; ix<nx; ix+=stepx){
			out[iy+ix*ny]=in[ix+iy*nx];
		}
	}
}
#endif
