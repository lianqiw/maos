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
#include <typeinfo>
#include "common.h"
#include "types.h"
#include "cudata.h"
extern int cuda_dedup; //Set to 1 during setup and 0 during simulation
/**
   Without type conversion. Enable asynchrous transfer. It is asynchrous only if
   called allocated pinned memory.
*/

template<typename M, typename N> 
inline void type_convert(M *out, const N* in, int nx){
    for(int i=0; i<nx; i++){
	out[i]=(M)in[i];
    }
}

template<>
inline void type_convert<float2, double2>(float2*out, const double2* in, int nx){
    for(int i=0; i<nx; i++){
	out[i].x=in[i].x;
	out[i].y=in[i].y;
    }
}
template<>
inline void type_convert<float2, dcomplex>(float2* out, const dcomplex* in, int nx){
    type_convert(out, (const double2*)in, nx);
}

/*Async copy does not make sense here because malloc pinned memory is too expensive.*/
template<typename M, typename N>
void cp2gpu(M**dest, const N*src, int nx, int ny, cudaStream_t stream=0){
    if(!src) return;
    uint64_t key=0;
    if(cuda_dedup && !*dest){
	key=hashlittle(src, nx*ny*sizeof(N), 0);
	key=(key<<32) | (nx*ny);
	if(cudata->memhash->count(key)){
	    *dest=(M*)(*cudata->memhash)[key];
	    return;
	}
    }
    if(!*dest){
	DO(cudaMalloc(dest, nx*ny*sizeof(M)));
    }
    if(cuda_dedup){
	(*cudata->memhash)[key]=*dest;
    }
    M* from=0;
    if(sizeof(M)!=sizeof(N)){
	if(!cuda_dedup && cudata->memcache->count((void*)(*dest))){
	    /*We cache the array used for the conversion. It is important that
	     * no gpu memory is allocated/freed at every cycle*/
	    from=(M*)(*cudata->memcache)[(void*)(*dest)];
	}else{
	    from=(M*)malloc(sizeof(M)*nx*ny);
	    if(!cuda_dedup){
		(*cudata->memcache)[(void*)*dest]=(void*)from;
	    }
	}
	type_convert(from, src, nx*ny);
    }else{
	from=(M*)(src);
    }

    if(stream==(cudaStream_t)0){
	DO(cudaMemcpy(*dest, from, sizeof(M)*nx*ny, cudaMemcpyHostToDevice));
    }else{
	DO(cudaMemcpyAsync(*dest, from, sizeof(M)*nx*ny, cudaMemcpyHostToDevice, stream));
    }
    if((void*)from !=(void*)src && cuda_dedup) {
	free(from);
    }
}

template<typename M, typename N> inline void
cp2gpu(cumat<M>**dest, const N*src, int nx, int ny, cudaStream_t stream=0){
    if(!src) return;
    if(*dest){
	if((*dest)->nx*(*dest)->ny!=nx*ny){
	    error("cumat is %ldx%ld, input is %dx%d\n", (*dest)->nx, (*dest)->ny, nx, ny);
	}
	cp2gpu(&((*dest)->p), src, nx, ny, stream);
    }else{
	M *p=NULL;
	cp2gpu(&p, src, nx, ny, stream);
	*dest=new cumat<M>(nx, ny, p);
    }

}
inline void cp2gpu(Real**dest, const X(mat)*src, cudaStream_t stream=0){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx, src->ny, stream);
}
inline void cp2gpu(curmat**dest, const X(mat)*src, cudaStream_t stream=0){
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
void cp2gpu(Real (* restrict *dest)[2], const loc_t *src);
void cp2gpu(curcell *restrict *dest, const dcell *src);
void cp2gpu(cuccell *restrict *dest, const ccell *src);


void gpu2gpu(cumap_t **dest0, cumap_t *source, int nps);

void cuspmul (Real *y, cusp *A, const Real *x, int ncol, char trans,
	      Real alpha, cusparseHandle_t handle);

void gpu_write(const Real *p, int nx, int ny, const char *format, ...);
void gpu_write(const Comp *p, int nx, int ny, const char *format, ...);
void gpu_write(const int *p, int nx, int ny, const char *format, ...);
void add2cpu(double * restrict *dest, double alpha, Real *src, double beta, int n, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(float * restrict *dest, double alpha, Real *src, double beta, int n, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(smat **out, float alpha, const curmat *in, float beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(dmat **out, double alpha, const curmat *in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(zmat **out, float alpha, const cucmat *in, float beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(cmat **out, double alpha, const cucmat *in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(dcell **out, double alpha, const curcell *in, double beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void add2cpu(scell **out, float alpha, const curcell *in, float beta, cudaStream_t stream, pthread_mutex_t* mutex=0);
void cp2cpu(dmat **out, const curmat *in, cudaStream_t stream);
void cp2cpu(smat **out, const curmat *in, cudaStream_t stream);
void cp2cpu(cmat **out, const cucmat *in, cudaStream_t stream);
void cp2cpu(zmat **out, const cucmat *in, cudaStream_t stream);
void cp2cpu(scell **out, const curcell *in, cudaStream_t stream);
void cp2cpu(dcell **out, const curcell *in, cudaStream_t stream);
void cp2cpu(zcell **out, const cuccell *in, cudaStream_t stream);
void cp2cpu(ccell **out, const cuccell *in, cudaStream_t stream);
void cellarr_cur(struct cellarr *ca, int i, const curmat *A, cudaStream_t stream);
void cellarr_cuc(struct cellarr *ca, int i, const cucmat *A, cudaStream_t stream);
void cellarr_curcell(struct cellarr *ca, int i, const curcell *A, cudaStream_t stream);
void cellarr_cuccell(struct cellarr *ca, int i, const cuccell *A, cudaStream_t stream);
#endif
