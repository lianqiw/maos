/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#define CUDA_API_PER_THREAD_DEFAULT_STREAM 1
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <cufft.h>
#include <cuComplex.h>
typedef double2 dcomplex;
typedef float2 fcomplex;
#include "gpu.h"
#define NG1D 128
#define NG2D 8
#define WRAP_SIZE 32
#define REDUCE_WRAP 8
#define REDUCE_WRAP_LOG2 3
#define DIM_REDUCE (WRAP_SIZE*REDUCE_WRAP)
#define REDUCE_STRIDE (WRAP_SIZE+WRAP_SIZE/2+1)
#if CUDA_DOUBLE == 1
typedef double2 Comp;
typedef double  Real;
typedef dmat rmat;
#define CUSP(A) cusparseD##A
#define CUBL(A) cublasD##A
#define X(A) d##A
#define C(A) c##A
#define Z(A) A
#define zfarr_mat zfarr_dmat
#define FFT_T_C2C CUFFT_Z2Z
#define FFT_T_C2R CUFFT_Z2D
#define FFT_T_R2C CUFFT_D2Z
#define FFT_C2C cufftExecZ2Z
#define FFT_C2R cufftExecZ2D
#define FFT_R2C cufftExecD2Z
#define M_REAL M_DBL
#define M_COMP M_CMP
#else
typedef float2 Comp;
typedef float Real;
typedef smat rmat;
#define CUSP(A) cusparseS##A
#define CUBL(A) cublasS##A
#define X(A) s##A
#define C(A) z##A
#define Z(A) A##f
#define zfarr_mat zfarr_smat
#define FFT_T_C2C CUFFT_C2C
#define FFT_T_C2R CUFFT_C2R
#define FFT_T_R2C CUFFT_R2C
#define FFT_C2C cufftExecC2C
#define FFT_C2R cufftExecC2R
#define FFT_R2C cufftExecR2C
#define M_REAL M_FLT
#define M_COMP M_ZMP
#endif
extern "C"{
    void cudaProfilerStart(void);
    void cudaProfilerStop(void);
}
#undef EPS
#define EPS 1.e-5 //Float has limited, 6 digit, resolution.
typedef Real Real2[2];
int cuda_free(void *p);
/*static int tot_mem=0; */
#undef cudaMalloc
#undef cudaFree
inline int CUDAMALLOC(void **p, size_t size){
    return cudaMalloc(p,size);
}
inline int CUDAFREE(void *p){
    return cuda_free(p);
}
#define DEBUG_MEM 0
#if DEBUG_MEM
#define cudaMalloc(p,size) ({info("%ld cudaMalloc for %s: %9lu Byte\n",pthread_self(),#p, size);CUDAMALLOC((void**)(void*)p,size);})
#define cudaFree(p)        ({info("%ld cudaFree   for %s\n", pthread_self(),#p);CUDAFREE((void*)p);})
#else
#define cudaMalloc(p,size) CUDAMALLOC((void**)(void*)p,size)
#define cudaFree(p) CUDAFREE((void*)p)
#endif

int current_gpu();
#define DO(A...) ({int _ans=(int)(A); if(_ans!=0&& _ans!=cudaErrorNotReady){error("GPU %d error %d, %s\n", current_gpu(), _ans, cudaGetErrorString((cudaError_t)_ans));}})
#define DORELAX(A...) ({int _ans=(int)(A); static int counter=0; if(_ans!=0&& _ans!=cudaErrorNotReady){counter++; if(counter>5) error("GPU %d error %d, %s\n", current_gpu(), _ans, cudaGetErrorString((cudaError_t)_ans));else warning("GPU %d error %d, %s\n", current_gpu(), _ans, cudaGetErrorString((cudaError_t)_ans));}})
#define cudaCallocHostBlock(P,N) ({DO(cudaMallocHost(&(P),N)); memset(P,0,N);})
#define cudaCallocBlock(P,N)     ({DO(cudaMalloc(&(P),N));     DO(cudaMemset(P,0,N)); CUDA_SYNC_DEVICE;})
#define cudaCallocHost(P,N,stream) ({DO(cudaMallocHost(&(P),N)); DO(cudaMemsetAsync(P,0,N,stream));})
#define cudaCalloc(P,N,stream) ({DO(cudaMalloc(&(P),N));DO(cudaMemsetAsync(P,0,N,stream));})
#define TO_IMPLEMENT error("Please implement")
__inline__ static __host__ __device__ float2 operator*(const float2 &a, const float2 &b){
    return cuCmulf(a,b);
}
__inline__ static __host__ __device__ float2 operator+(const float2 &a, const float2 &b){
    return cuCaddf(a,b);
}
__inline__ static __host__ __device__ float2&operator*=(float2 &a, const float2 &b){
    a=cuCmulf(a,b);
    return a;
}
__inline__ static __host__ __device__ float2 operator*(const float2 &a, const float b){
    float2 tmp;
    tmp.x=a.x*b;
    tmp.y=a.y*b;
    return tmp;
}
__inline__ static __host__ __device__ float2&operator*=(float2 &a, const float b){
    a.x*=b;
    a.y*=b;
    return a;
}
__inline__ static __host__ __device__ double2 operator*(const double2 &a, const double2 &b){
    return cuCmul(a,b);
}
__inline__ static __host__ __device__ double2 operator+(const double2 &a, const double2 &b){
    return cuCadd(a,b);
}
__inline__ static __host__ __device__ double2&operator*=(double2 &a, const double2 &b){
    a=cuCmul(a,b);
    return a;
}
__inline__ static __host__ __device__ double2 operator*(const double2 &a, const double b){
    double2 tmp;
    tmp.x=a.x*b;
    tmp.y=a.y*b;
    return tmp;
}
__inline__ static __host__ __device__ double2&operator*=(double2 &a, const double b){
    a.x*=b;
    a.y*=b;
    return a;
}
inline void* malloc4async(size_t N){
    void *tmp;
    cudaMallocHost(&tmp, N);
    return tmp;
}
inline void free4async(void *P){
    cudaFreeHost(P);
}
extern int NULL_STREAM;
#if DEBUG
#define CUDA_CHECK_ERROR DO(cudaGetLastError())
#else
#define CUDA_CHECK_ERROR
#endif
#define CUDA_SYNC_STREAM ({CUDA_CHECK_ERROR;DORELAX(cudaStreamSynchronize(stream));})
#define CUDA_SYNC_DEVICE ({CUDA_CHECK_ERROR;DORELAX(cudaDeviceSynchronize());})
#define STREAM_NEW(stream) if(NULL_STREAM) {stream=0; info2("Warning NULL stream\n");} else DO(cudaStreamCreate(&stream))
#define STREAM_DONE(stream) if(!NULL_STREAM) DO(cudaStreamSynchronize(stream),cudaStreamDestroy(stream))
#define HANDLE_NEW(handle,stream) ({DO(cublasCreate(&handle)); DO(cublasSetStream(handle, stream));})
#define SPHANDLE_NEW(handle,stream) ({DO(cusparseCreate(&handle)); DO(cusparseSetStream(handle, stream));})
#define HANDLE_DONE(handle) DO(cublasDestroy(handle))
#define SPHANDLE_DONE(sphandle) DO(cusparseDestroy(sphandle))
#define adpind(A,i) ((A)->nx>1?(A)->p[i]:(A)->p[0])

#define DIM(nsa,nb) MIN((nsa+nb-1)/nb,NG1D),MIN((nsa),nb)
#define REDUCE(nsa) MIN((nsa+DIM_REDUCE-1)/DIM_REDUCE,NG1D), DIM_REDUCE
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ <200
#define NTH2 16
#else
#define NTH2 32
#endif
#define DIM2(nx,ny,nb) dim3(MIN((nx+nb-1)/(nb),NG2D),MIN((ny+nb-1)/(nb),NG2D)),dim3(MIN(nx,nb),MIN(ny,nb))
#define DIM3(nx,ny,nb,nbz) dim3(MIN((nx+nb-1)/(nb),NG2D),MIN((ny+nb-1)/(nb),NG2D),nbz),dim3(MIN(nx,nb),MIN(ny,nb))

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 200
#define MEMCPY_D2D cudaMemcpyDeviceToDevice
#else
#define MEMCPY_D2D cudaMemcpyDefault
#endif
#define MEMCPY_D2H cudaMemcpyDeviceToHost
#define MEMCPY_H2D cudaMemcpyHostToDevice
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

/*For timing asynchronous kernels*/
#define EVENT_INIT(n)				\
    const int NEVENT=n;				\
    Real times[NEVENT];			\
    cudaEvent_t event[NEVENT]={0};		\
    for(int i=0; i<NEVENT; i++){		\
	DO(cudaEventCreate(&event[i]));		\
    }
#define EVENT_TIC(i) DO(cudaEventRecord(event[i], stream))
#define EVENT_TOC			       \
    stream.sync();times[0]=0;		       \
    for(int i=1; i<NEVENT; i++){	       \
	DO(cudaEventElapsedTime		       \
	   (&times[i], event[i-1], event[i])); \
	times[i]*=1e3;			       \
	times[0]+=times[i];		       \
    }						
    
#define EVENT_DEINIT				\
    for(int i=0; i<NEVENT; i++){		\
	 DO(cudaEventDestroy(event[i]));	\
    }

extern const char *cufft_str[];
INLINE void CUFFT2(cufftHandle plan, Comp *in, Comp *out, int dir){
    LOCK_CUFFT;						
    int _ans=FFT_C2C(plan, in, out, dir);		
    UNLOCK_CUFFT;						
    if(_ans){						
	error("cufft failed: %s\n", cufft_str[_ans]);	
    }							
}
INLINE void CUFFTR2C(cufftHandle plan, const Real *in, Comp *out){	
    LOCK_CUFFT;						
    int _ans=FFT_R2C(plan, (Real*)in, out);		
    UNLOCK_CUFFT;						
    if(_ans){						
	error("cufft failed: %s\n", cufft_str[_ans]);	
    }							
}
INLINE void CUFFTC2R(cufftHandle plan, const Comp *in, Real *out){	
    LOCK_CUFFT;						
    int _ans=FFT_C2R(plan, (Comp*)in, out);		
    UNLOCK_CUFFT;						
    if(_ans){						
	error("cufft failed: %s\n", cufft_str[_ans]);	
    }							
}
#define CUFFT(plan,in,dir) CUFFT2(plan,in,in,dir)
typedef struct stream_t{
    cudaStream_t stream;
    cublasHandle_t handle;
    cusparseHandle_t sphandle;
    stream_t(){
	init();
    }
    void init(){
	STREAM_NEW(stream);//this takes a few seconds for each gpu for the first time.
	HANDLE_NEW(handle, stream);
	SPHANDLE_NEW(sphandle, stream);
    }
    ~stream_t(){
	deinit();
    }
    void deinit(){
	SPHANDLE_DONE(sphandle);
	HANDLE_DONE(handle);
	STREAM_DONE(stream);
    }
    void reset(){//to place on correct gpu.
	deinit();
	init();
    }
    void sync(){
	//assert(this);
	DO(cudaStreamSynchronize(stream));
    }
    operator cudaStream_t&(){
	//assert(this);
	return stream;
    }
    operator cublasHandle_t&(){
	//assert(this);
	return handle;
    }
    operator cusparseHandle_t&(){
	//assert(this);
	return sphandle;
    }
private:
    stream_t(const stream_t &);
    stream_t & operator=(const stream_t &);
}stream_t;
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
inline void X(pagelock)(X(mat) *A, ...){
    va_list ap;
    va_start(ap, A);
    do{
	cudaHostRegister(A->p, A->nx*A->ny*sizeof(Real), cudaHostRegisterPortable);
	A=va_arg(ap, X(mat) *);
    }while(A);
    va_end(ap);
}
inline void X(pageunlock)(X(mat) *A, ...){
    va_list ap;
    va_start(ap, A);
    do{
	cudaHostUnregister(A->p);
	A=va_arg(ap, X(mat) *);
    }while(A);
    va_end(ap);
}
#endif
