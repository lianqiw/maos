/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "common.h"
#include <typeinfo>
#include "types.h"
#include "cudata.h"

/**
   Without type conversion. Enable asynchrous transfer. It is asynchrous only if
   called allocated pinned memory.
*/
/*
template<typename M> inline
void cp2gpu(M**dest, const M*src, int nx, int ny, cudaStream_t stream=0){
    if(!src) return;
    uint64_t key=0;
    if(cuda_dedup){
	TIC;tic;
	key=hashlittle(src, nx*ny*sizeof(M))<<32 | (nx*ny);
	toc2("hashlittle");
	iterator it=cudata->memhash.find(key);
	if(it!=cudata->memhash.end()){//find
	    info2("memory at %p is already at %p in gpu.\n", src, (void*)it);
	    if(*dest){
		warning("Free %p\n", *dest);
		cudaFree(*dest);
	    }
	    *dest=(M*)it;
	    return;
	}
    }
    if(!*dest){
	DO(cudaMalloc(dest, nx*ny*sizeof(M)));
    }
    if(cuda_dedup){
	cudata->memhash.insert(key, *dest);
	info2("memory at %p is copied to %p in gpu.\n", src, *dest);
    }

    if(stream==(cudaStream_t)0){
	DO(cudaMemcpy(*dest, src, sizeof(M)*nx*ny, cudaMemcpyHostToDevice));
    }else{
	DO(cudaMemcpyAsync(*dest, src, sizeof(M)*nx*ny, cudaMemcpyHostToDevice, stream));
    }
}
*/
/**
   With type conversion
*/
template<typename M, typename N> 
inline void type_convert(M *out, const N* in, int nx){
    for(int i=0; i<nx; i++){
	out[i]=(M)in[i];
    }
}

template<>
inline void type_convert<float2, dcomplex>(float2* out, const dcomplex* in, int nx){
    double (*tmp)[2]=(double(*)[2])in;
    for(int i=0; i<nx; i++){
	out[i].x=(float)tmp[i][0];
	out[i].y=(float)tmp[i][1];
    }
}
template<> /*This one should never be called. float complex is defined by C99*/
inline void type_convert<float2, float complex>(float2* out, const float complex* in, int nx){
    memcpy(out, in, sizeof(float2)*nx);
}
template<>
inline void type_convert<float2, double2>(float2* out, const double2* in, int nx){
    for(int i=0; i<nx; i++){
	out[i].x=in[i].x;
	out[i].y=in[i].y;
    }
}
/*Async copy does not make sense here because malloc pinned memory is too expensive.*/
template<typename M, typename N>
void cp2gpu(M**dest, const N*src, int nx, int ny, cudaStream_t stream=0){
    if(!src) return;
    uint64_t key=0;
    extern int cuda_dedup;
    if(cuda_dedup && !*dest){
	key=hashlittle(src, nx*ny*sizeof(N), 0);
	key=(key<<32) | (nx*ny);
	if(cudata->memhash->count(key)){
	    info2("memory at %p is already at %p in gpu.\n", src, 
		  (*cudata->memhash)[key]);
	    *dest=(M*)(*cudata->memhash)[key];
	    return;
	}
    }
    M* from;
    if(sizeof(M)!=sizeof(N)){
	from=(M*)malloc(sizeof(M)*nx*ny);
	type_convert(from, src, nx*ny);
    }else{
	if(typeid(M)!=typeid(N)){
	    warning("No convert from %s to %s\n", typeid(M).name(), typeid(N).name());
	}
	from=(M*)(src);
    }
    //don't call previous cp2gpu() as it may call itself.
    if(!*dest){
	DO(cudaMalloc(dest, nx*ny*sizeof(M)));
    }
    if(cuda_dedup){
	(*cudata->memhash)[key]=*dest;
	info2("memory at %p is copied to %p in gpu.\n", src, *dest);
    }
    if(stream==(cudaStream_t)0){
	DO(cudaMemcpy(*dest, from, sizeof(M)*nx*ny, cudaMemcpyHostToDevice));
    }else{
	DO(cudaMemcpyAsync(*dest, from, sizeof(M)*nx*ny, cudaMemcpyHostToDevice, stream));
    }
    if((void*)from !=(void*)src) free(from);
}
/*template<typename M> inline
void cp2gpu(cumat<M>**dest, const M*src, int nx, int ny, cudaStream_t stream=0){
    if(!src) return;
    if(!*dest){
	*dest=new cumat<M>(nx, ny);
    }
    cp2gpu(&((*dest)->p), src, nx, ny, stream);
}*/
template<typename M, typename N> inline void
cp2gpu(cumat<M>**dest, const N*src, int nx, int ny, cudaStream_t stream=0){
    if(!src) return;
    if(*dest){
	cp2gpu(&((*dest)->p), src, nx, ny, stream);
    }else{
	M *p=NULL;
	cp2gpu(&p, src, nx, ny, stream);
	*dest=new cumat<M>(nx, ny, p);
    }

}
inline void cp2gpu(float**dest, const smat*src, cudaStream_t stream=0){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
inline void cp2gpu(curmat**dest, const smat*src, cudaStream_t stream=0){
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
void cp2gpu(cumap_t **dest, map_t **source, int nps);
void cp2gpu(cusp **dest, const dsp *src, int tocsr);
void cp2gpu(cusp **dest, const spcell *src, int tocsr);
void cp2gpu(cuspcell **dest, const spcell *src, int tocsr);
void cp2gpu(float (* restrict *dest)[2], const loc_t *src);
void cp2gpu(curcell *restrict *dest, const dcell *src);
void cp2gpu(cuccell *restrict *dest, const ccell *src);

void cuspmul (float *y, cusp *A, const float *x, int ncol, char trans,
	      float alpha, cusparseHandle_t handle);

void gpu_write(const float *p, int nx, int ny, const char *format, ...);
void gpu_write(const fcomplex *p, int nx, int ny, const char *format, ...);
void gpu_write(const int *p, int nx, int ny, const char *format, ...);
void cp2cpu(double * restrict *dest, double alpha, float *src, double beta, int n, cudaStream_t stream, pthread_mutex_t* mutex=0);
void cp2cpu(dmat **out, double alpha, const curmat *in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void cp2cpu(smat **out, const curmat *in, cudaStream_t stream);
void cp2cpu(zmat **out, const cucmat *in, cudaStream_t stream);
void cp2cpu(dcell **out, double alpha, const curcell *in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void cp2cpu(scell **out, const curcell *in, cudaStream_t stream);
void cp2cpu(zcell **out, const cuccell *in, cudaStream_t stream);
inline void cp2cpu(dmat **out, const curmat *in, cudaStream_t stream){
    cp2cpu(out, 0, in, 1, stream);
}
inline void cp2cpu(dcell **out, const curcell *in, cudaStream_t stream){
    cp2cpu(out, 0, in, 1, stream);
}
inline void add2cpu(dmat **out, const curmat *in, cudaStream_t stream, pthread_mutex_t* mutex=0){
    cp2cpu(out, 1, in, 1, stream, mutex);
}
inline void add2cpu(dcell **out, const curcell *in, cudaStream_t stream, pthread_mutex_t* mutex=0){
    cp2cpu(out, 1, in, 1, stream, mutex);
}
void cellarr_cur(struct cellarr *ca, int i, const curmat *A, cudaStream_t stream);
void cellarr_cuc(struct cellarr *ca, int i, const cucmat *A, cudaStream_t stream);
void cellarr_curcell(struct cellarr *ca, int i, const curcell *A, cudaStream_t stream);
void cellarr_cuccell(struct cellarr *ca, int i, const cuccell *A, cudaStream_t stream);
#endif
