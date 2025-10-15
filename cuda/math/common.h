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
#ifndef AOS_CUDA_COMMON_H
#define AOS_CUDA_COMMON_H
#define AOS_CUDA_H //prevents numtype.h from defining complex number operations
/**
 * Notice that multiple routines uses stream=0 as default parameter to be used during preparation. When stream==0, all other stream operation should be blocked to ensure data integrity. So
 * 		1. CUDA_API_PER_THREAD_DEFAULT_STREAM should not be defined
 * 		2. stream should be created with default behavior, not unblocking.
 * */

//#define CUDA_API_PER_THREAD_DEFAULT_STREAM 1 #DO NOT Enable
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#define DISABLE_CUSPARSE_DEPRECATED
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <cufft.h>
#include <cuComplex.h>
#include <cuda_profiler_api.h>
#include <cusolverDn.h>

typedef double2 dcomplex;
typedef float2 fcomplex;
#include "gpu_math.h"
#define NG1D 128
#define NG2D 8
#define WRAP_SIZE 32
#define REDUCE_WRAP 8
#define REDUCE_WRAP_LOG2 3
#define DIM_REDUCE (WRAP_SIZE*REDUCE_WRAP)
#define REDUCE_STRIDE (WRAP_SIZE+WRAP_SIZE/2+1)
#if CUDA_DOUBLE
typedef double2 Comp;
typedef double  Real;
//typedef dmat rmat;
#define CUSP(A) cusparseD##A
#define CUBL(A) cublasD##A
#define CUDA_R CUDA_R_64F
#define X(A) d##A
#define XC(A) c##A
#define Z(A) A
#define FFT_T_C2C CUFFT_Z2Z
#define FFT_T_C2R CUFFT_Z2D
#define FFT_T_R2C CUFFT_D2Z
#define FFT_C2C cufftExecZ2Z
#define FFT_C2R cufftExecZ2D
#define FFT_R2C cufftExecD2Z
#define MCU_REAL M_DBL
#define MCU_COMP M_CMP
#else
typedef float2 Comp;
typedef float Real;
//typedef smat rmat;
#define CUSP(A) cusparseS##A
#define CUBL(A) cublasS##A
#define CUDA_R CUDA_R_32F
#define X(A) s##A
#define XC(A) z##A
#define Z(A) A##f
#define FFT_T_C2C CUFFT_C2C
#define FFT_T_C2R CUFFT_C2R
#define FFT_T_R2C CUFFT_R2C
#define FFT_C2C cufftExecC2C
#define FFT_C2R cufftExecC2R
#define FFT_R2C cufftExecR2C
#define MCU_REAL M_FLT
#define MCU_COMP M_ZMP
#endif
#if defined (DLONG) && 0
#define CUSPARSE_INDEX CUSPARSE_INDEX_64I
typedef spint Spint;
#else
#define CUSPARSE_INDEX CUSPARSE_INDEX_32I //CUDA < 10.0 can only use integer indices
typedef int Spint;
#endif
#undef EPS
#if CUDA_DOUBLE && !CPU_SINGLE
#define EPS 1.e-16
#else
#define EPS 1.e-6
#endif
typedef Real Real2[2];
int mycudaFree(void* p);//Unreference deduplicated memory
int mycudaMalloc(void** p, size_t size);
/*static int tot_mem=0; */
#undef cudaMalloc
#undef cudaFree
#define DEBUG_MEM 0
#if DEBUG_MEM
#define cudaMalloc(p,size) ({dbg("%ld cudaMalloc for %s: %9lu Byte\n",pthread_self(),#p, size);mycudaMalloc((void**)(void*)p,size);})
#define cudaFree(p)        ({dbg("%ld cudaFree   for %s\n", pthread_self(),#p);mycudaFree((void*)p);})
#else
#define cudaMalloc(p,size) mycudaMalloc((void**)(void*)p,size)
#define cudaFree(p) mycudaFree((void*)p)
#endif

int current_gpu();
#define TO_IMPLEMENT error("Please implement")
static inline __host__ __device__ float2 operator*(const float2& a, const float2& b){
	return cuCmulf(a, b);
}
static inline __host__ __device__ float2 operator+(const float2& a, const float2& b){
	return cuCaddf(a, b);
}
static inline __host__ __device__ float2& operator*=(float2& a, const float2& b){
	a=cuCmulf(a, b);
	return a;
}
static inline __host__ __device__ float2 operator*(const float2& a, const float b){
	float2 tmp;
	tmp.x=a.x*b;
	tmp.y=a.y*b;
	return tmp;
}
static inline __host__ __device__ float2& operator*=(float2& a, const float b){
	a.x*=b;
	a.y*=b;
	return a;
}
static inline __host__ __device__ double2 operator*(const double2& a, const double2& b){
	return cuCmul(a, b);
}
static inline __host__ __device__ double2 operator+(const double2& a, const double2& b){
	return cuCadd(a, b);
}
static inline __host__ __device__ double2& operator*=(double2& a, const double2& b){
	a=cuCmul(a, b);
	return a;
}
static inline __host__ __device__ double2 operator*(const double2& a, const double b){
	double2 tmp;
	tmp.x=a.x*b;
	tmp.y=a.y*b;
	return tmp;
}
static inline __host__ __device__ double2& operator*=(double2& a, const double b){
	a.x*=b;
	a.y*=b;
	return a;
}
template<typename M, typename N>
__device__ __host__ void type_convert(M& out, const N in){
	out=static_cast<M>(in);
}
template<>
__device__ __host__ inline void type_convert<float2, double2>(float2& out, const double2 in){
	out.x=static_cast<float>(in.x);
	out.y=static_cast<float>(in.y);
}
template<>
__device__ __host__ inline void type_convert<double2, float2>(double2 &out, const float2 in){
	out.x=static_cast<double>(in.x);
	out.y=static_cast<double>(in.y);
}
template<>
__device__ __host__ inline void type_convert<float2, float>(float2& out, const float in){
	out.x=in;
	out.y=0;
}
template<>
__device__ __host__ inline void type_convert<float2, double>(float2 &out, const double in){
	out.x=in;
	out.y=0;
}
template<>
__device__ __host__ inline void type_convert<double2, float>(double2 &out, const float in){
	out.x=in;
	out.y=0;
}
template<>
__device__ __host__ inline void type_convert<double2, double>(double2 &out, const double in){
	out.x=in;
	out.y=0;
}
extern int NULL_STREAM;
#if DEBUG
#define CUDA_CHECK_ERROR DO(cudaGetLastError())
#else
#define CUDA_CHECK_ERROR
#endif
#define DO(A...) ({int _ans=(int)(A); if(_ans!=0&& _ans!=cudaErrorNotReady){print_backtrace(); error("GPU %d error %d, %s\n", current_gpu(), _ans, cudaGetErrorString((cudaError_t)_ans));}})
#define DORELAX(A...) ({int _ans=(int)(A); static int _counter=0; if(_ans!=0&& _ans!=cudaErrorNotReady){_counter++; if(_counter>5) error("GPU %d error %d, %s\n", current_gpu(), _ans, cudaGetErrorString((cudaError_t)_ans));else warning("GPU %d error %d, %s\n", current_gpu(), _ans, cudaGetErrorString((cudaError_t)_ans));}})
#define CUDA_SYNC_STREAM if(stream!=0){CUDA_CHECK_ERROR;DORELAX(cudaStreamSynchronize(stream));}
#define CUDA_SYNC_DEVICE ({CUDA_CHECK_ERROR;DORELAX(cudaDeviceSynchronize());})
//STREAM_NEW creates stream that does not synchronize with the legacy stream 0.
#define STREAM_NEW(stream) if(NULL_STREAM) {stream=0; info("Warning NULL stream\n");} else DO(cudaStreamCreateWithFlags(&stream,cudaStreamDefault))
#define STREAM_DONE(stream) if(!NULL_STREAM) DO(cudaStreamSynchronize(stream),cudaStreamDestroy(stream))
#define HANDLE_NEW(handle,stream) ({DO(cublasCreate(&handle)); DO(cublasSetStream(handle, stream));})
#define SPHANDLE_NEW(handle,stream) ({DO(cusparseCreate(&handle)); DO(cusparseSetStream(handle, stream));})
#define HANDLE_DONE(handle) DO(cublasDestroy(handle))
#define SPHANDLE_DONE(sphandle) DO(cusparseDestroy(sphandle))
#define adpind(A,i) ((A)->nx>1?(A)->p[i]:(A)->p[0])

#define DIM(nsa,nb) MIN((nsa+nb-1)/nb,NG1D),MIN((nsa),nb)
#define REDUCE(nsa) MIN((nsa+DIM_REDUCE-1)/DIM_REDUCE,NG1D), DIM_REDUCE
//Launch kernel configured to handled nz nx*ny array indexing, 
#define DIM3(nx,ny,nb,nz) dim3(MIN((nx+nb-1)/(nb),NG2D),MIN((ny+nb-1)/(nb),NG2D),nz),dim3(MIN(nx,nb),MIN(ny,nb))
#define DIM2(nx,ny,nb) DIM3(nx,ny,nb,1)
#define D2D cudaMemcpyDefault
#define D2H cudaMemcpyDeviceToHost
#define H2D cudaMemcpyHostToDevice
/*
  Notice that the CUDA FFT 4.0 is not thread safe!. Our FFT is a work around of
  the problem by using mutex locking to makesure only 1 thread is calling FFT. */
#if CUDA_VERSION < 4010
extern pthread_mutex_t cufft_mutex;
#define LOCK_CUFFT LOCK(cufft_mutex)
#define UNLOCK_CUFFT UNLOCK(cufft_mutex)
#else /*cufft>=4.1 is thread safe. no need lock.*/
#define LOCK_CUFFT
#define UNLOCK_CUFFT
#endif

#if TIMING == 1
#define ctoc(A)						\
    if(itoc==ntoc) error("please enlarge ntoc");	\
    cudaEventCreate(&events[itoc]);			\
    nametoc[itoc]=A;					\
    cudaEventRecord(events[itoc], stream);		\
    itoc++;					

#define ctoc_init(nin)				\
    const int ntoc=nin;				\
    cudaEvent_t events[ntoc]={0};		\
    const char *nametoc[ntoc]={0};		\
    int itoc=0;					\
    ctoc("init");				\

#define ctoc_final(A...)						\
    CUDA_SYNC_STREAM;							\
    float ms;								\
    char msg[1024];							\
    int ct=snprintf(msg, sizeof(msg), A);				\
    ct+=snprintf(msg+ct, sizeof(msg)-ct, ":");				\
    for(int ii=1; ii<itoc; ii++){					\
	DO(cudaEventElapsedTime(&ms, events[ii-1], events[ii]));	\
	ct+=snprintf(msg+ct,sizeof(msg)-ct," %s %.1f,",nametoc[ii],ms); \
    }									\
    DO(cudaEventElapsedTime(&ms, events[0], events[itoc-1]));		\
    info("%s Tot %.1f ms\n", msg, ms);					\
    for(int ii=0; ii<itoc; ii++){					\
	DO(cudaEventDestroy(events[ii]));				\
	events[ii]=0;							\
    }

#else
#define ctoc(A)
#define ctoc_init(A)
#define ctoc_final(A...)
#endif

extern const char* cufft_str[];
static inline void CUFFT2(cufftHandle plan, Comp* in, Comp* out, int dir){
	LOCK_CUFFT;
	int _ans=FFT_C2C(plan, in, out, dir);
	UNLOCK_CUFFT;
	if(_ans){
		error("cufft failed: %s\n", cufft_str[_ans]);
	}
}
static inline void CUFFTR2C(cufftHandle plan, const Real* in, Comp* out){
	LOCK_CUFFT;
	int _ans=FFT_R2C(plan, (Real*)in, out);
	UNLOCK_CUFFT;
	if(_ans){
		error("cufft failed: %s\n", cufft_str[_ans]);
	}
}
static inline void CUFFTC2R(cufftHandle plan, const Comp* in, Real* out){
	LOCK_CUFFT;
	int _ans=FFT_C2R(plan, (Comp*)in, out);
	UNLOCK_CUFFT;
	if(_ans){
		error("cufft failed: %s\n", cufft_str[_ans]);
	}
}
#define CUFFT(plan,in,dir) CUFFT2(plan,in,in,dir)
class stream_t{
	cudaStream_t stream;
	cublasHandle_t handle;
	cusparseHandle_t sphandle;
	cusolverDnHandle_t dnhandle; 

public:
	stream_t(){
		init();
	}
	void init(){
		STREAM_NEW(stream);//this takes a few seconds for each gpu for the first time.
		HANDLE_NEW(handle, stream);
		SPHANDLE_NEW(sphandle, stream);
		DO(cusolverDnCreate(&dnhandle)); DO(cusolverDnSetStream(dnhandle, stream));
	}
	~stream_t(){
		deinit();
	}
	void deinit(){
		cusolverDnDestroy(dnhandle);
		SPHANDLE_DONE(sphandle);
		HANDLE_DONE(handle);
		STREAM_DONE(stream);
	}
	void reset(){//to place on correct gpu.
		deinit();
		init();
	}
	void sync(){
		DO(cudaStreamSynchronize(stream));
	}
	operator cudaStream_t(){
		return stream;
	}
	/*
		the handle_t are pointers. 
	  	cuda 10.0 complains multiple conversion function to void * when multiple operators are defined.
	}*/
	cublasHandle_t blas(){
		return handle;
	}
	cusparseHandle_t sparse(){
		return sphandle;
	}
	cusolverDnHandle_t dn(){
		return dnhandle;
	}
private://do not allow copy.
	stream_t& operator=(const stream_t& in);
/*
	deinit();
	stream=in.stream;
	handle=in.handle;
	sphandle=in.sphandle;
	return *this;
*/
	stream_t(const stream_t&);

};
typedef struct event_t{
	cudaEvent_t event;
	event_t(unsigned int flag=cudaEventDefault){
		DO(cudaEventCreateWithFlags(&event, flag));
	}
	~event_t(){
		cudaEventDestroy(event);
	}
	void record(cudaStream_t stream){
		DO(cudaEventRecord(event, stream));
	}
	operator cudaEvent_t(){
	//assert(this);
		return event;
	}
}event_t;
static inline void X(pagelock)(X(mat)* A, ...){
	va_list ap;
	va_start(ap, A);
	do{
		cudaHostRegister(A->p, A->nx*A->ny*sizeof(Real), cudaHostRegisterPortable);
		A=va_arg(ap, X(mat)*);
	} while(A);
	va_end(ap);
}
static inline void X(pageunlock)(X(mat)* A, ...){
	va_list ap;
	va_start(ap, A);
	do{
		cudaHostUnregister(A->p);
		A=va_arg(ap, X(mat)*);
	} while(A);
	va_end(ap);
}
#endif
