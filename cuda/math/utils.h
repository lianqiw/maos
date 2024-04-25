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

#ifndef AOS_CUDA_UTILS_H
#define AOS_CUDA_UTILS_H
#include <typeinfo>
#include <map>
#include "common.h"
#include "types.h"
#include "cublas.h"
#include "kernel.h"
#include "cucmat.h"
#include "curmat.h"
class lock_t{
	pthread_mutex_t& mutex;
	int enable;
public:
	lock_t(pthread_mutex_t& _mutex, int _enable=1):mutex(_mutex), enable(_enable){
		if(enable) LOCK(mutex);
	}
	~lock_t(){
		if(enable) UNLOCK(mutex);
	}
};
class cumemcache_t{
	public:
	std::map<uint64_t, void *> memhash;/*For reuse constant GPU memory.*/
	std::map<void *, int> memcount; /*Store count of reused memory*/
	void *memcache;/*For reuse temp array for type conversion.*/
	long nmemcache;
	pthread_mutex_t memmutex;
	cumemcache_t() :memcache(NULL), nmemcache(0){
		pthread_mutex_init(&memmutex, 0);
	}
	~cumemcache_t(){
		free(memcache); memcache=NULL; nmemcache=0;
	}
};
extern cumemcache_t cumemcache;
extern int cuda_dedup; //Set to 1 during setup and 0 during simulation
/**
   Without type conversion. Enable asynchrous transfer. It is asynchrous only if
   called allocated pinned memory.
*/

template<typename M, typename N>
void type_convert(M* out, const N* in, long nx){
	for(int i=0; i<nx; i++){
		out[i]=static_cast<M>(in[i]);
	}
}
template<>
inline void type_convert<float2, double2>(float2* out, const double2* in, long nx){
	for(int i=0; i<nx; i++){
		out[i].x=static_cast<float>(in[i].x);
		out[i].y=static_cast<float>(in[i].y);
	}
}
template<>
inline void type_convert<double2, float2>(double2* out, const float2* in, long nx){
	for(int i=0; i<nx; i++){
		out[i].x=static_cast<double>(in[i].x);
		out[i].y=static_cast<double>(in[i].y);
	}
}

template<typename M, typename N>
void cp2gpu(M* dest, const N* src, long nx, long ny, cudaStream_t stream=0){
	/*{
		static int same_size=0, size_last=0;
		if(size_last==nx*ny*sizeof(N)){
			same_size++;
			if(same_size==10){
				error("cp2gpu: %d copies to GPU of size %ld KiB\n", same_size, size_last>>10);
			}
		}else{
			same_size=0;
		}
		dbg("cp2gpu: %d copies to GPU of size %ld KiB\n", same_size, size_last>>10);
		size_last=nx*ny*sizeof(N);
	}*/
	M* from=0;
	int free_from=0;
	if(cumemcache.memcount.count(dest)&&cumemcache.memcount[dest]>1){
		error("Should not copy to deduped pointer %p. Count=%d\n", dest,
			cumemcache.memcount[dest]);
	}
	if(sizeof(M)!=sizeof(N)){
		long memsize=nx*ny*sizeof(M);
		if(memsize>20000000){//Too large. Don't cache.
			from=(M*)malloc(memsize);
			free_from=1;
		} else{
			LOCK(cumemcache.memmutex);
			if(cumemcache.nmemcache<memsize){
				cumemcache.nmemcache=memsize;
				/*dbg("GPU%d: Enlarge mem cache to %ld: %p->",
				  current_gpu(), memsize, cumemcache.memcache);*/
				cumemcache.memcache=realloc(cumemcache.memcache, cumemcache.nmemcache);
				//dbg("%p\n", cumemcache.memcache);
			}
			from=(M*)cumemcache.memcache;
		}
		//TIC;tic;
		type_convert(from, src, nx*ny);
		//toc("converting %ld elements", nx*ny);
	} else{
		from=(M*)(src);
	}
	DO(cudaMemcpyAsync(dest, from, sizeof(M)*nx*ny, cudaMemcpyHostToDevice, stream));
	if(free_from){
		if(stream!=0){
			CUDA_SYNC_STREAM;
		}
		free(from);
	} else if(sizeof(M)!=sizeof(N)){
		UNLOCK(cumemcache.memmutex);
	}
}
/*Async copy does not make sense here because malloc pinned memory is too expensive.*/
template<typename M, typename N>
void cp2gpu_dedup(M** dest, const N* src, long nx, long ny, cudaStream_t stream=0){
	if(!src){
		error("src=null\n");
	}
	//dbg("cp2gpu_dedup: copy %p to GPU.\n", src);
	uint64_t key=0;
	int skip_copy=0;
	int record_mem=0;
	if(cuda_dedup&&!*dest){
		key=hashlittle(src, nx*ny*sizeof(N), nx*ny);
		int igpu;
		cudaGetDevice(&igpu);
		key=hashlittle(&igpu, sizeof(int), key);//put GPU index as part of fingerprint.
		lock_t tmp(cumemcache.memmutex);
		if(cumemcache.memhash.count(key)){
			*dest=(M*)cumemcache.memhash[key];
			if(cumemcache.memcount[*dest]){//valid memory
				cumemcache.memcount[*dest]++;
				skip_copy=1;//no need to copy again
				//dbg("cp2gpu_dedup: increase reference to data: %p\n", *dest);
			} else{
				cumemcache.memhash.erase(key);
				cumemcache.memcount.erase(*dest);
				warning("cp2gpu_dedup: remove invalid reference data: %p\n", *dest);
				*dest=0;
			}
		}
	} else if(!cuda_dedup&&*dest){
		//Avoid overriding previously referenced memory
		lock_t tmp(cumemcache.memmutex);
		if(cumemcache.memcount.count(*dest)&&cumemcache.memcount[*dest]>1){
			warning("cp2gpu_dedup: deferencing data: %p\n", *dest);
			cumemcache.memcount[*dest]--;
			*dest=0;
		}
	}

	if(!*dest){
		if(cuda_dedup) record_mem=1;
		*dest=(M*)new Gpu<M>[nx*ny];
	}
	if(!skip_copy){
		cp2gpu(*dest, src, nx, ny, stream);
	}
	if(record_mem){
		//dbg("cp2gpu_dedup: record reference to data: %p\n", *dest);
		lock_t tmp(cumemcache.memmutex);
		cumemcache.memhash[key]=*dest;
		cumemcache.memcount[*dest]=1;
	}
}

template<typename M, typename N> static inline void
cp2gpu(NumArray<M, Gpu>& dest, const N* src, long nx, long ny, cudaStream_t stream=0){
	if(!src||!nx||!ny) return;
	if(dest){
		if(dest.N()!=nx*ny){
			error("Array is %ldx%ld, input is %ldx%ld\n", dest.Nx(), dest.Ny(), nx, ny);
		}
		cp2gpu(dest(), src, nx, ny, stream);
	} else{
		M* tmp=0;
		cp2gpu_dedup(&tmp, src, nx, ny, stream);
		dest=NumArray<M, Gpu>(nx, ny, tmp, 1);
	}
}
static inline void cp2gpu_dedup(Real** dest, const dmat* src, cudaStream_t stream=0){
	if(!src) return;
	cp2gpu_dedup(dest, src->p, src->nx, src->ny, stream);
}
//#if CPU_SINGLE==0
/*template <typename T, typename S>
static inline void cp2gpu(NumArray<T, Gpu>& dest, NumArray<S, Cpu>& src, cudaStream_t stream=0){
	cp2gpu(dest, src(), src.Nx(), src.ny(), stream);
}*/
static inline void cp2gpu(curmat& dest, const dmat* src, cudaStream_t stream=0){
	if(!src) return;
	cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
static inline void cp2gpu(cucmat &dest, const cmat *src, cudaStream_t stream=0){
	if(!src) return;
	cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
static inline void cp2gpu(curmat &dest, const smat *src, cudaStream_t stream=0){
	if(!src) return;
	cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
static inline void cp2gpu(cucmat &dest, const zmat *src, cudaStream_t stream=0){
	if(!src) return;
	cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
#if CPU_SINGLE==0
static inline void cp2gpu(cudmat &dest, const dmat *src, cudaStream_t stream=0){
	if(!src) return;
	cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
static inline void cp2gpu(cuzmat &dest, const cmat *src, cudaStream_t stream=0){
	if(!src) return;
	cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
static inline void cp2gpu(cudmat &dest, const smat *src, cudaStream_t stream=0){
	if(!src) return;
	cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
static inline void cp2gpu(cuzmat &dest, const zmat *src, cudaStream_t stream=0){
	if(!src) return;
	cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
#endif
void cp2gpu(cumapcell& dest, const mapcell* source);
void cp2gpu(cusp& dest, const dsp* src, int tocsr);
void cp2gpu(cusp& dest, const dspcell* src, int tocsr);
void cp2gpu(cuspcell& dest, const dspcell* src, int tocsr);
void cp2gpu(curmat& dest, const loc_t* src);
void cp2gpu(curcell& dest, const dcell* src);
void cp2gpu(cuccell& dest, const ccell* src);

void gpu_write(const Real* p, long nx, long ny, const char* format, ...) CHECK_ARG(4);
void gpu_write(const Comp* p, long nx, long ny, const char* format, ...) CHECK_ARG(4);
void gpu_write(const int* p, long nx, long ny, const char* format, ...) CHECK_ARG(4);
void add2cpu(float* restrict* dest,Real alpha, Real* src, Real beta, long n, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(smat** out, float alpha, const curmat& in, float beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(zmat** out, float alpha, const cucmat& in, float beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(scell** out, float alpha, const curcell& in, float beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(zcell **out, float alpha, const cuccell &in, float beta, cudaStream_t stream, pthread_mutex_t *mutex=0);
void add2cpu(dcell **out, real alpha, const curcell &in, real beta, cudaStream_t stream, pthread_mutex_t *mutex=0);
void add2cpu(ccell **out, real alpha, const cuccell &in, real beta, cudaStream_t stream, pthread_mutex_t *mutex=0);
#if CPU_SINGLE==0
void add2cpu(real* restrict* dest, real alpha, Real* src, real beta, long n, cudaStream_t stream, pthread_mutex_t* mutex=0);
#endif
void add2cpu(dmat** out, real alpha, const curmat& in, real beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(cmat** out, real alpha, const cucmat& in, real beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void cp2cpu(dmat** out, const curmat& in, cudaStream_t stream=0);
void cp2cpu(cmat **out, const cucmat &in, cudaStream_t stream=0);
void cp2cpu(smat **out, const curmat &in, cudaStream_t stream=0);
void cp2cpu(zmat **out, const cucmat &in, cudaStream_t stream=0);
#if CPU_SINGLE==0
void cp2cpu(dmat **out, const cudmat &in, cudaStream_t stream=0);
void cp2cpu(cmat **out, const cuzmat &in, cudaStream_t stream=0);
void cp2cpu(smat **out, const cudmat &in, cudaStream_t stream=0);
void cp2cpu(zmat **out, const cuzmat &in, cudaStream_t stream=0);
#endif
void cp2cpu(scell** out, const curcell& in, cudaStream_t stream=0);
void cp2cpu(zcell** out, const cuccell& in, cudaStream_t stream=0);
void cp2cpu(dcell** out, const curcell& in, cudaStream_t stream=0);
void cp2cpu(ccell** out, const cuccell& in, cudaStream_t stream=0);

void zfarr_push_scale(struct zfarr *ca, int i, const curmat &A, Real scale, cudaStream_t stream=0);
void zfarr_push_scale(struct zfarr *ca, int i, const cucmat &A, Real scale, cudaStream_t stream=0);
void zfarr_push_scale(struct zfarr *ca, int i, const curcell &A, Real scale, cudaStream_t stream=0);
void zfarr_push_scale(struct zfarr *ca, int i, const cuccell &A, Real scale, cudaStream_t stream=0);
void drawopdamp_gpu(const char* fig, loc_t* loc, const curmat& opd,  cudaStream_t stream,
	const dmat* amp, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(10);
void drawpsf_gpu(const char* fig, curmat& psf, int count, cudaStream_t stream, int log, real psfmin,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...) CHECK_ARG(10);
void curdraw_gpu(const char *fig, curmat &psf, int count, cudaStream_t stream, int log,
	const char *title, const char *xlabel, const char *ylabel,
	const char *format, ...) CHECK_ARG(9);
void cucdraw_gpu(const char *fig, cucmat &psf, int count, cudaStream_t stream, int log,
	const char *title, const char *xlabel, const char *ylabel,
	const char *format, ...) CHECK_ARG(9);
#endif
