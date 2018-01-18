/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "common.h"
#include "types.h"
#include "cudata.h"
class lock_t{
    pthread_mutex_t &mutex;
    int enable;
public:
    lock_t(pthread_mutex_t &_mutex, int _enable=1):mutex(_mutex),enable(_enable){
	if(enable) LOCK(mutex);
    }
    ~lock_t(){
	if(enable) UNLOCK(mutex);
    }
};

extern int cuda_dedup; //Set to 1 during setup and 0 during simulation
/**
   Without type conversion. Enable asynchrous transfer. It is asynchrous only if
   called allocated pinned memory.
*/

template<typename M, typename N> 
inline void type_convert(M *out, const N* in, int nx){
    for(int i=0; i<nx; i++){
	out[i]=static_cast<M>(in[i]);
    }
}

template<>
inline void type_convert<float2, double2>(float2*out, const double2* in, int nx){
    for(int i=0; i<nx; i++){
	out[i].x=static_cast<float>(in[i].x);
	out[i].y=static_cast<float>(in[i].y);
    }
}
template<>
inline void type_convert<double2, float2>(double2*out, const float2* in, int nx){
    for(int i=0; i<nx; i++){
	out[i].x=static_cast<double>(in[i].x);
	out[i].y=static_cast<double>(in[i].y);
    }
}
/*template<>
inline void type_convert<float2, dcomplex>(float2* out, const dcomplex* in, int nx){
    type_convert(out, (const double2*)in, nx);
}
*/
template<typename M, typename N>
void cp2gpu(M*dest, const N*src, int nx, int ny, cudaStream_t stream=0){
    M* from=0;
    int free_from=0;
    if(sizeof(M)!=sizeof(N)){
	long memsize=nx*ny*sizeof(M);
	if(memsize>20000000 || !cudata){//Too large. Don't cache.
	    from=(M*)malloc(memsize);
	    free_from=1;
	}else{
	    LOCK(cudata->memmutex);
	    if(cudata->nmemcache<memsize){
		cudata->nmemcache=memsize;
		/*info("GPU%d: Enlarge mem cache to %ld: %p->", 
		  current_gpu(), memsize, cudata->memcache);*/
		cudata->memcache=realloc(cudata->memcache, cudata->nmemcache);
		//info("%p\n", cudata->memcache);
	    }
	    from=(M*)cudata->memcache;
	}
	type_convert(from, src, nx*ny);
    }else{
	from=(M*)(src);
    }
    if(stream==(cudaStream_t)0){
	DO(cudaMemcpy(dest, from, sizeof(M)*nx*ny, cudaMemcpyHostToDevice));
    }else{
	DO(cudaMemcpyAsync(dest, from, sizeof(M)*nx*ny, cudaMemcpyHostToDevice, stream));
    }
    if(free_from){
	if(stream!=(cudaStream_t)0){
	    CUDA_SYNC_STREAM;
	}
	free(from);
    }else if(sizeof(M)!=sizeof(N)){
	UNLOCK(cudata->memmutex);
    }
}
/*Async copy does not make sense here because malloc pinned memory is too expensive.*/
template<typename M, typename N>
int cp2gpu(M**dest, const N*src, int nx, int ny, cudaStream_t stream=0){
    if(!src) return -1;
    uint64_t key=0;
    int own=0;//whether the returned process owns the pointer (yes if created).
    if(cuda_dedup && !*dest){
	key=hashlittle(src, nx*ny*sizeof(N), 0);
	key=(key<<32) | (nx*ny);
	lock_t tmp(cudata->memmutex);
	if(cudata->memhash.count(key)){
	    *dest=(M*)cudata->memhash[key];
	    cudata->memcount[*dest]++;
	    return 0;
	}
    }else if(!cuda_dedup && *dest && cudata){
	//Avoid overriding previously referenced memory
	lock_t tmp(cudata->memmutex);
	if(cudata->memcount.count(*dest) && cudata->memcount[*dest]>1){
	    //info("Deferencing data: %p\n", *dest);
	    cudata->memcount[*dest]--;
	    *dest=0;
	}
    }
    if(!*dest){
	own=1;
	DO(cudaMalloc(dest, nx*ny*sizeof(M)));
	if(cuda_dedup){
	    lock_t tmp(cudata->memmutex);
	    cudata->memhash[key]=*dest;
	    cudata->memcount[*dest]=1;
	}
    }
    cp2gpu(*dest, src, nx, ny, stream);
    return own;
}

template<typename M, typename N> inline void
cp2gpu(cumat<M>& dest, const N*src, int nx, int ny, cudaStream_t stream=0){
    if(!src || !nx || !ny) return;
    if(dest){
	if(dest.N()!=nx*ny){
	    error("cumat is %ldx%ld, input is %dx%d\n", dest.Nx(), dest.Ny(), nx, ny);
	}
	cp2gpu(dest.P(), src, nx, ny, stream);
    }else{
	M *tmp=0;
	int own=cp2gpu(&tmp, src, nx, ny, stream);
	dest=cumat<M>(nx, ny, tmp, own);
    }
}
inline void cp2gpu(Real**dest, const dmat*src, cudaStream_t stream=0){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
inline void cp2gpu(curmat &dest, const dmat*src, cudaStream_t stream=0){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
inline void cp2gpu(curmat &dest, const smat*src, cudaStream_t stream=0){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny, stream);
}

inline void cp2gpu(cucmat &dest, const cmat*src, cudaStream_t stream=0){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
inline void cp2gpu(cucmat &dest, const zmat*src, cudaStream_t stream=0){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny, stream);
}

/* A few special cases to avoid N match to cell*/
template<typename M>
void cp2gpu(M**dest, const dmat*src){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny);
}

template<typename M>
void cp2gpu(M**dest, const cmat*src){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny);
}

template<typename M>
void cp2gpu(M**dest, const zmat*src){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny);
}
void cp2gpu(cumapcell &dest, const mapcell *source);
void cp2gpu(cusp &dest, const dsp *src, int tocsr);
void cp2gpu(cusp &dest, const dspcell *src, int tocsr);
void cp2gpu(cuspcell &dest, const dspcell *src, int tocsr);
void cp2gpu(curmat &dest, const loc_t *src);
void cp2gpu(curcell &dest, const dcell *src);
void cp2gpu(cuccell &dest, const ccell *src);

void cuspmul (Real *y, cusp &A, const Real *x, int ncol, char trans,
	      Real alpha, cusparseHandle_t handle);

void gpu_write(const Real *p, int nx, int ny, const char *format, ...);
void gpu_write(const Comp *p, int nx, int ny, const char *format, ...);
void gpu_write(const int *p, int nx, int ny, const char *format, ...);
void add2cpu(double * restrict *dest, double alpha, Real *src, double beta, int n, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(float * restrict *dest, double alpha, Real *src, double beta, int n, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(smat **out, float alpha, const curmat &in, float beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(dmat **out, double alpha, const curmat &in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(zmat **out, float alpha, const cucmat &in, float beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(cmat **out, double alpha, const cucmat &in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(dcell **out, double alpha, const curcell &in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(scell **out, float alpha, const curcell &in, float beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void cp2cpu(dmat **out, const curmat &in, cudaStream_t stream);
void cp2cpu(smat **out, const curmat &in, cudaStream_t stream);
void cp2cpu(cmat **out, const cucmat &in, cudaStream_t stream);
void cp2cpu(zmat **out, const cucmat &in, cudaStream_t stream);
void cp2cpu(scell **out, const curcell &in, cudaStream_t stream);
void cp2cpu(dcell **out, const curcell &in, cudaStream_t stream);
void cp2cpu(zcell **out, const cuccell &in, cudaStream_t stream);
void cp2cpu(ccell **out, const cuccell &in, cudaStream_t stream);
void zfarr_cur(struct zfarr *ca, int i, const curmat &A, cudaStream_t stream);
void zfarr_cuc(struct zfarr *ca, int i, const cucmat &A, cudaStream_t stream);
void zfarr_curcell(struct zfarr *ca, int i, const curcell &A, cudaStream_t stream);
void zfarr_cuccell(struct zfarr *ca, int i, const cuccell &A, cudaStream_t stream);
void drawopdamp_gpu(const char *fig, loc_t *loc, const curmat &opd, cudaStream_t stream, 
		    const double *amp, double *zlim,
		    const char *title, const char *xlabel, const char *ylabel,
		    const char* format,...) CHECK_ARG(10);
#endif
