/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

extern "C"
{
#include "gpu.h"
}
#include <cuda.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cufft.h>
extern "C"{
    void cudaProfilerStart(void);
    void cudaProfilerStop(void);
}
#undef EPS
#define EPS 1.e-5 //Float has limited, 6 digit, resolution.
#define DEBUG_MEM 0
#if DEBUG_MEM
/*static int tot_mem=0; */
#undef cudaMalloc
inline int CUDAMALLOC(float **p, size_t size){
    return cudaMalloc((float**)p,size);
}
inline int CUDAFREE(float *p){
    return cudaFree(p);
}
#define cudaMalloc(p,size) ({info("%ld cudaMalloc for %s: %9lu Byte\n",pthread_self(),#p, size);CUDAMALLOC((float**)p,size);})
#define cudaFree(p)        ({info("%ld cudaFree   for %s\n", pthread_self(),#p);CUDAFREE((float*)p);})
#endif
#define DO(A...) ({int ans=(int)(A); if(ans!=0) error("(cuda) %d: %s\n", ans, cudaGetErrorString((cudaError_t)ans));})
#define cudaCallocHostBlock(P,N) ({DO(cudaMallocHost(&(P),N)); memset(P,0,N);})
#define cudaCallocBlock(P,N)     ({DO(cudaMalloc(&(P),N));     DO(cudaMemset(P,0,N)); CUDA_SYNC_DEVICE;})
#define cudaCallocHost(P,N,stream) ({DO(cudaMallocHost(&(P),N)); DO(cudaMemsetAsync(P,0,N,stream));})
#define cudaCalloc(P,N,stream) ({DO(cudaMalloc(&(P),N));DO(cudaMemsetAsync(P,0,N,stream));})
#define TO_IMPLEMENT error("Please implement")
__host__ __device__ static __inline__ void cuCscalef(cuFloatComplex x, float a){
    x.x*=a;
    x.y*=a;
}
inline void* malloc4async(int N){
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ <200
    void *tmp;
    cudaMallocHost(&tmp, N);
    return tmp;
#else
    return malloc(N);
#endif
}
inline void free4async(void *P){
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ <200
    cudaFreeHost(P);
#else
    free(P);
#endif
}

#define CONCURRENT 0
#if CONCURRENT
#define CUDA_SYNC_STREAM				\
    while(cudaStreamQuery(stream)!=cudaSuccess){	\
	if(THREAD_RUN_ONCE){/* no jobs to do*/		\
	    cudaStreamSynchronize(stream);		\
	    break;					\
	}						\
    }
#else
#define CUDA_SYNC_STREAM DO(cudaStreamSynchronize(stream))
#endif

#define CUDA_SYNC_DEVICE DO(cudaDeviceSynchronize())
#define TIMING 0
#if TIMING == 1
extern int nstream;
#define STREAM_NEW(stream) ({DO(cudaStreamCreate(&stream));info2("nstream=%d\n",lockadd(&nstream,1)+1);})
#define STREAM_DONE(stream) ({DO(cudaStreamSynchronize(stream),cudaStreamDestroy(stream));info2("nstream=%d\n",lockadd(&nstream,-1)-1);})
#else
#define STREAM_NEW(stream) DO(cudaStreamCreate(&stream))
#define STREAM_DONE(stream) DO(cudaStreamSynchronize(stream),cudaStreamDestroy(stream))
#endif
#undef TIMING
#define HANDLE_NEW(handle,stream) ({DO(cublasCreate(&handle)); DO(cublasSetStream(handle, stream));})
#define SPHANDLE_NEW(handle,stream) ({DO(cusparseCreate(&handle)); DO(cusparseSetKernelStream(handle, stream));})
#define HANDLE_DONE(handle) DO(cublasDestroy(handle))
#define SPHANDLE_DONE(sphandle) DO(cusparseDestroy(sphandle))
#define adpind(A,i) ((A)->nx>1?(A)->p[i]:(A)->p[0])
#define MYSPARSE 0

#define DIM(nsa,nb) MIN((nsa+nb-1)/nb,NG1D),MIN((nsa),nb)
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
  Notice that the CUDA FFT 4.0 is not thread safe!. Our FFT is a walk around of
  the problem by using mutex locking to makesure only 1 thread is calling FFT. */
#if CUDA_VERSION < 4010
extern pthread_mutex_t cufft_mutex;
#define LOCK_CUFFT LOCK(cufft_mutex)
#define UNLOCK_CUFFT UNLOCK(cufft_mutex)
#else
/*cufft 4.1 is thread safe. no need lock.*/
#define LOCK_CUFFT
#define UNLOCK_CUFFT
#endif

/*For timing asynchronous kernels*/
#define EVENT_INIT(n)				\
    const int NEVENT=n;				\
    float times[NEVENT];			\
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
	 DO(cudaEventDestroy(&event[i]));	\
    }

extern const char *cufft_str[];
#define CUFFT2(plan,in,out,dir) do{				\
	LOCK_CUFFT;						\
	int ans=cufftExecC2C(plan, in, out, dir);		\
	UNLOCK_CUFFT;						\
	if(ans){						\
	    error("cufft failed: %s\n", cufft_str[ans]);	\
	}							\
    }while(0)
#define CUFFTR2C(plan,in,out) do{				\
	LOCK_CUFFT;						\
	int ans=cufftExecR2C(plan, in, out);			\
	UNLOCK_CUFFT;						\
	if(ans){						\
	    error("cufft failed: %s\n", cufft_str[ans]);	\
	}							\
    }while(0)
#define CUFFTC2R(plan,in,out) do{				\
	LOCK_CUFFT;						\
	int ans=cufftExecC2R(plan, in, out);			\
	UNLOCK_CUFFT;						\
	if(ans){						\
	    error("cufft failed: %s\n", cufft_str[ans]);	\
	}							\
    }while(0)
#define CUFFT(plan,in,dir) CUFFT2(plan,in,in,dir)
typedef struct stream_t{
    cudaStream_t stream;
    cublasHandle_t handle;
    cusparseHandle_t sphandle;
    stream_t(){
	STREAM_NEW(stream);//this takes a few seconds for each gpu for the first time.
	HANDLE_NEW(handle, stream);
	SPHANDLE_NEW(sphandle, stream);
    }
    ~stream_t(){
	SPHANDLE_DONE(sphandle);
	HANDLE_DONE(handle);
	STREAM_DONE(stream);
    }
    void sync(){
	assert(this);
	DO(cudaStreamSynchronize(stream));
    }
    operator cudaStream_t(){
	assert(this);
	return stream;
    }
    operator cublasHandle_t(){
	assert(this);
	return handle;
    }
    operator cusparseHandle_t(){
	assert(this);
	return sphandle;
    }
private:
    stream_t(const stream_t &);
    stream_t & operator=(const stream_t &);
}stream_t;
inline void spagelock(smat *A, ...){
    va_list ap;
    va_start(ap, A);
    do{
	cudaHostRegister(A->p, A->nx*A->ny*sizeof(float), cudaHostRegisterPortable);
	A=va_arg(ap, smat *);
    }while(A);
    va_end(ap);
}
inline void spageunlock(smat *A, ...){
    va_list ap;
    va_start(ap, A);
    do{
	cudaHostUnregister(A->p);
	A=va_arg(ap, smat *);
    }while(A);
    va_end(ap);
}
#endif
