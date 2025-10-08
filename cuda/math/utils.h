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
#ifndef AOS_CUDA_UTILS_H
#define AOS_CUDA_UTILS_H
#include <typeinfo>
#include <map>
#include "types.h"
#include "cublas.h"
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
	std::map<uint64_t, void *> memhash;/*hash is mapped to memory address*/
	std::map<void *, uint32_t[2]> memcount; /*Store usage count of memory*/
	pthread_mutex_t mutex_hash;
	long nsave=0;
	void *memcache=NULL;/*For reuse temp array for type conversion.*/
	long nmemcache=0;
	pthread_mutex_t mutex_cache;
	
	cumemcache_t(){
		pthread_mutex_init(&mutex_cache, 0);
		pthread_mutex_init(&mutex_hash, 0);
	}
	~cumemcache_t(){
		free(memcache); memcache=NULL; nmemcache=0;
	}
};
extern cumemcache_t cumemcache;
extern int cuda_dedup; //Set to 1 during setup and 0 during simulation

template<typename M, typename N>
void type_convert(M *out, const N *in, long nx){
	for(int i=0; i<nx; i++){
		type_convert(out[i], in[i]);
	}
}

template<typename M, typename N>
void cp2gpu(M* dest, const N* src, long nx, long ny, cudaStream_t stream=0){
	M* from=0;
	int free_from=0;
	if(cumemcache.memcount.count(dest)&&cumemcache.memcount[dest][0]>1){
		error("Should not copy to deduped pointer %p. Count=%d\n", dest,
			cumemcache.memcount[dest][0]);
	}
	if(sizeof(M)!=sizeof(N)){//use cache memory for type conversion
		long memsize=nx*ny*sizeof(M);
		if(memsize>20000000){//Too large. Don't cache.
			from=(M*)malloc(memsize);
			free_from=1;
		} else{
			LOCK(cumemcache.mutex_cache);
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
	//since we are not using pinned memory, memcpyAsync is sync to host code.
	DO(cudaMemcpyAsync(dest, from, sizeof(M)*nx*ny, cudaMemcpyHostToDevice, stream));
	if(sizeof(M)!=sizeof(N)){
		if(free_from){
			if(stream!=0){
				CUDA_SYNC_STREAM;
			}
			free(from);
		} else {
			UNLOCK(cumemcache.mutex_cache);
		}
	}
}
/*Convert and copy CPU data to GPU with deduplication if cuda_dedup is set.*/
template<typename M, typename N>
void cp2gpu_dedup(M** dest, const N* src, long nx, long ny, cudaStream_t stream=0){
	if(!src){
		error("src=null\n");
	}
	//dbg("cp2gpu_dedup: copy %p to GPU.\n", src);
	uint32_t key=0;
	int skip_copy=0;
	lock_t tmp(cumemcache.mutex_hash, ((cuda_dedup&&!*dest) || (!cuda_dedup&&*dest)));
	if(cuda_dedup&&!*dest){
		int igpu;
		cudaGetDevice(&igpu);
		key=hashlittle(src, nx*ny*sizeof(N), nx*ny);
		key=hashlittle(&igpu, sizeof(int), key);//put GPU index as part of fingerprint.	
		if(cumemcache.memhash.count(key)){//memory is already in the correct gpu
			*dest=(M*)cumemcache.memhash[key];//memory address
			if(cumemcache.memcount.count(*dest) && cumemcache.memcount[*dest][0]>0){//sanity check
				cumemcache.memcount[*dest][0]++;
				skip_copy=1;//no need to copy again
				cumemcache.nsave+=nx*ny*sizeof(N);
				//dbg("cp2gpu_dedup: increase reference to data: %p, memory saved is %ld MB.\n", *dest, cumemcache.nsave>>20);
			} else{
				warning("cp2gpu_dedup: erase invalid deduplicated data: %p (key=%x, exist=%d, count=%d, size is %ld)\n", 
					*dest, key, (int)cumemcache.memcount.count(*dest), cumemcache.memcount[*dest][0], nx*ny*sizeof(N));
				cumemcache.memhash.erase(key);
				cumemcache.memcount.erase(*dest);
				*dest=0;
			}
		}//else: first time
	} else if(!cuda_dedup&&*dest){
		//Avoid overriding previously deduplicate memory
		if(cumemcache.memcount.count(*dest)){
			if(cumemcache.memcount[*dest][0]>1){
				warning("cp2gpu_dedup: allocate new memory for deduplicated data: %p (count=%d)\n", *dest, cumemcache.memcount[*dest][0]);
				cumemcache.memcount[*dest][0]--;
				*dest=0;
			}else{
				warning("cp2gpu_dedup: reuse memory for singular deduplicated data: %p\n", *dest);
				cumemcache.memcount.erase(*dest);
			}
		}
	}//else: regular copy to existing target

	if(!*dest){
		*dest=(M*)new Gpu<M>[nx*ny];
		if(cuda_dedup){
			cumemcache.memhash[key]=*dest;
			/*if(nx*ny*sizeof(N)==149448){
				info("key=%x, malloc %p\n", key, *dest);
				print_backtrace();
			}*/
			cumemcache.memcount[*dest][0]=1;
			cumemcache.memcount[*dest][1]=key;
		}
	}
	if(!skip_copy){
		cp2gpu(*dest, src, nx, ny, stream);
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
void cp2gpu(cumapcell& dest, const mapcell* source);
void cp2gpu(cusp& dest, const dsp* src, int tocsr);
void cp2gpu(cusp& dest, const dspcell* src, int tocsr);
void cp2gpu(cuspcell& dest, const dspcell* src, int tocsr);
//void cp2gpu(curmat& dest, const loc_t* src);
void cp2gpu(curmat &dest, const_anyarray src_, cudaStream_t stream=0);
void cp2gpu(cucmat &dest, const_anyarray src_, cudaStream_t stream=0);
template <typename T>
void cp2gpu(NumCell<T, Gpu>& dest, const_anycell src_, cudaStream_t stream=0){
	const cell* src=src_.c;
	if(!src){
		dest.Zero();
		return;
	}
	if(!dest){
		long nc=src->nx*src->ny;
		NumArray<long> nx(nc);
		NumArray<long> ny(nc);
		for(long i=0; i<nc; i++){
			if(src->p[i]){
				nx[i]=src->p[i]->nx;
				ny[i]=src->p[i]->ny;
			} else{
				nx[i]=0;
				ny[i]=0;
			}
		}
		dest=NumCell<T, Gpu>(src->nx, src->ny, nx(), ny());
	} else if(dest.Nx()!=src->nx||dest.Ny()!=src->ny){
		error("Mismatch: %ldx%ld vs %ldx%ld\n",
			dest.Nx(), dest.Ny(), src->nx, src->ny);
	}
	if(src->m){
		//warning("cp2gpu: use m to M\n");
		cp2gpu(dest.M(), src->m, stream);
	} else{
		for(int i=0; i<src->nx*src->ny; i++){
			cp2gpu(dest[i], src->p[i], stream);
		}
	}
}
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
