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
#include "utils.h"
#include "curmat.h"
#include <pthread.h>
const char* cufft_str[]={
	"success",
	"invalid plan",
	"allocation failed",
	"",
	"invalid value",
	"internal errlr",
	"exec failed (error elsewhere caused cufft to fail)",
	"setup failed",
	"invalid size"
};
#ifndef I
#define I (__extension__ 1.0iF)
#endif
#if CUDA_VERSION < 4010
pthread_mutex_t cufft_mutex=PTHREAD_MUTEX_INITIALIZER;
#endif
int cuda_dedup=0;//1: allow memory deduplication. Useful during setup for const memory.
cumemcache_t cumemcache;
/*
The generic one:
template<typename T, template<typename> class Dev>
template<typename S, template<typename> class Dev2>
void NumArray<T, Dev>::Copy(const NumArray<S, Dev2> &in, cudaStream_t stream){}
//Note that you cannot partially spatially a member function.
*/

/**
   Copy map_t to cumap_t. if type==1, use cudaArray, otherwise use Real
   array. Allow multiple calling to override the data.  */
void cp2gpu(cumapcell& dest, const mapcell* source){
	if(source->nx==0) return;
	if(!dest){
		dest=cumapcell(source->nx, 1);
	}
	//initalize or override parameters.
	for(int ips=0; ips<source->nx; ips++){
		dest[ips]=(source->p[ips]);
	}
	//The time is spend on loading data from disk/network. parallel does not help
	for(int ips=0; ips<source->nx; ips++){
		cp2gpu(dest[ips].p, (const dmat*)source->p[ips]);
	}
}

/*
  Convert a host dsp array to GPU sprase array.g
*/
cusp::cusp(const dsp* src_csc, /**<Source dsp in CSC*/
	int tocsr,        /**<0: Keep in CSC. 1: Convert to CSR */
	int transp        /**<1: transpose*/
){
	if(!src_csc) return;
	dsp* src_trans=0;
	const dsp* src=0;
	if(tocsr!=transp){
		src_trans=dsptrans(src_csc);
		src=src_trans;
	} else{
		src=src_csc;
	}
	ref=new cusp_ref;
	cp2gpu_dedup(&ref->p, src->pp, src->ny+1, 1);
	cp2gpu_dedup(&ref->i, src->pi, src->nzmax, 1);
	cp2gpu_dedup(&ref->x, src->px, src->nzmax, 1);
	ref->nx=src_csc->nx;
	ref->ny=src_csc->ny;
	ref->nzmax=src->nzmax;
	ref->type=tocsr?SP_CSR:SP_CSC;
	ref->count=1;
	if(src_trans){
		dspfree(src_trans);
	}

#if CUDA_VERSION >= 10000
	int nrow=(ref->type==SP_CSR?ref->nx:ref->ny);
	int ncol=(ref->type==SP_CSR?ref->ny:ref->nx);
	cusparseCreateCsr(&ref->desc, nrow, ncol, ref->nzmax, ref->p, ref->i, ref->x, CUSPARSE_INDEX, CUSPARSE_INDEX, CUSPARSE_INDEX_BASE_ZERO, CUDA_R);
#endif
}
void cp2gpu(cusp& dest, const dspcell* srcc, int tocsr){
	dsp* src=dspcell2sp(srcc);
	dest=cusp(src, tocsr);
	dspfree(src);
}
void cp2gpu(cuspcell& dest, const dspcell* src, int tocsr){
	if(!src) return;
	if(!dest){
		dest=cuspcell(src->nx, src->ny);
	}
	for(int i=0; i<src->nx*src->ny; i++){
		dest[i]=cusp(src->p[i], tocsr);
	}
}
void cp2gpu(curmat &dest, const_anyarray src_, cudaStream_t stream){
	if(src_.c==0){
		dest.Zero(stream);
	}else if((src_.c->id&M_DBL)==M_DBL){
		cp2gpu(dest, src_.dm->p, src_.c->nx, src_.c->ny, stream);
	}else if((src_.c->id&M_FLT)==M_FLT){
		cp2gpu(dest, src_.sm->p, src_.c->nx, src_.c->ny, stream);
	}else{
		error("Invalid: id=%u\n", src_.c->id);
	}
}
void cp2gpu(cucmat &dest, const_anyarray src_, cudaStream_t stream){
	if(src_.c==0){
		dest.Zero(stream);
	}else if((src_.c->id&M_CMP)==M_CMP){
		cp2gpu(dest, src_.cm->p, src_.c->nx, src_.c->ny, stream);
	}else if((src_.c->id&M_ZMP)==M_ZMP){
		cp2gpu(dest, src_.zm->p, src_.c->nx, src_.c->ny, stream);
	}else{
		error("Invalid: id=%u\n", src_.c->id);
	}
}

__attribute__((weak)) int current_gpu(){
	int igpu;
	cudaGetDevice(&igpu);
	return igpu;
}
template <typename T, typename R, typename S>
void scale_add(T* p1, R alpha, const S* p2, R beta, long n){
	for(long i=0; i<n; i++){
		p1[i]=p1[i]*alpha+p2[i]*beta;
	}
}
template <>
void scale_add<double2, double, Comp>(double2* p1, double alpha, const Comp* p2, double beta, long n){
	for(long i=0; i<n; i++){
		p1[i].x=p1[i].x*alpha+p2[i].x*beta;
		p1[i].y=p1[i].y*alpha+p2[i].y*beta;
	}
}
template <>
void scale_add<float2, float, Comp>(float2* p1, float alpha, const Comp* p2, float beta, long n){
	for(long i=0; i<n; i++){
		p1[i].x=p1[i].x*alpha+p2[i].x*beta;
		p1[i].y=p1[i].y*alpha+p2[i].y*beta;
	}
}
#if CPU_SINGLE==0
template <>
void scale_add<float2, float, comp>(float2 *p1, float alpha, const comp *p2, float beta, long n){
	for(long i=0; i<n; i++){
		p1[i].x=p1[i].x*alpha+p2[i].x*beta;
		p1[i].y=p1[i].y*alpha+p2[i].y*beta;
	}
}
#endif
/**
   Convert device (Real) array and add to host double.
   dest = alpha * dest + beta *src;
*/
template <typename R, typename T, typename S>
static void add2cpu(T* restrict* dest, R alpha, const S* src, R beta, long n,
	cudaStream_t stream, pthread_mutex_t* mutex){
	S* tmp=0;
	tmp=(S*)malloc(n*sizeof(S));
	//For transfers from device memory to pageable host memory, the function
	//will return only once the copy has completed.
	DO(cudaMemcpyAsync(tmp, src, n*sizeof(S), D2H, stream));
	if(!*dest){
		*dest=(T*)malloc(sizeof(T)*n);
	}
	T* restrict p=*dest;
	CUDA_SYNC_STREAM;//just in case.
	if(mutex) LOCK(*mutex);
	scale_add(p, alpha, tmp, beta, n);
	if(mutex) UNLOCK(*mutex);
	free(tmp);
}
#define add2cpu_mat(D, T, C)						\
    void add2cpu(D##mat **out, T alpha, const NumArray<C, Gpu> &in, T beta, \
		 cudaStream_t stream, pthread_mutex_t *mutex){		\
	if(!in){							\
	    if(*out) D##scale(*out, alpha);				\
	    return;							\
	}								\
	if(!*out) {							\
	    *out=D##new(in.Nx(), in.Ny());				\
		if(in.keywords.length()) (*out)->keywords=strdup(in.keywords.c_str());\
	}else{								\
	    assert((*out)->nx*(*out)->ny==in.N());			\
	}								\
	add2cpu(&(*out)->p, alpha, in(), beta, in.N(), stream, mutex); \
    }

add2cpu_mat(s, float, Real)
add2cpu_mat(z, float, Comp)
#if CPU_SINGLE==0
add2cpu_mat(s, float, real)
add2cpu_mat(z, float, comp)
#endif
#if COMP_SINGLE==0
add2cpu_mat(d, real, Real)
add2cpu_mat(c, real, Comp)
#if CPU_SINGLE==0
add2cpu_mat(d, real, real)
add2cpu_mat(c, real, comp)
#endif
#endif
#define add2cpu_cell(D, T, C)				    \
    void add2cpu(D##cell **out, T alpha, const C &in, T beta,	\
		 cudaStream_t stream, pthread_mutex_t *mutex){		\
	if(!in){							\
	    if(*out) D##cellscale(*out, alpha);				\
	    return;							\
	}								\
	if(!*out) {							\
	    *out=D##cellnew(in.Nx(), in.Ny());				\
		if(in.keywords.length()) (*out)->keywords=strdup(in.keywords.c_str());\
	}else{								\
	    assert((*out)->nx*(*out)->ny==in.N());			\
	}								\
	for(int i=0; i<in.N(); i++){					\
	    add2cpu((*out)->p+i, alpha, in[i], beta, stream, mutex);	\
	}								\
    }

add2cpu_cell(d, real, curcell)
add2cpu_cell(c, real, cuccell)

add2cpu_cell(s, float, curcell)
add2cpu_cell(z, float, cuccell)

#define cp2cpu_same(dmat,dzero,dnew,T)				\
void cp2cpu(dmat **out, const NumArray<T, Gpu> &in, cudaStream_t stream){ \
	if(!in) {							\
	    if(*out) dzero(*out);					\
	    return;							\
	}								\
	if(!*out) {\
		*out=dnew(in.Nx(), in.Ny());				\
		if(in.keywords.length()) (*out)->keywords=strdup(in.keywords.c_str());\
	}else{\
		if(PN((*out))!=in.N()){\
			error("Dimension mismatch: %ldx%ld vs %ldx%ld\n", NX((*out)), NY((*out)), in.Nx(), in.Ny());\
		}\
	}\
	dmat *pout=*out;						\
	DO(cudaMemcpyAsync(P(pout), in(), in.N()*sizeof(T),D2H, stream));	\
	CUDA_SYNC_STREAM;\
	if(pout->keywords) free(pout->keywords);				\
	if(in.keywords.length()) pout->keywords=strdup(in.keywords.c_str());	\
}
#if COMP_SINGLE==0
	cp2cpu_same(dmat, dzero, dnew, double)
	cp2cpu_same(cmat, czero, cnew, double2)
#endif
	cp2cpu_same(smat, szero, snew, float)
	cp2cpu_same(zmat, zzero, znew, float2)
#if ! CUDA_DOUBLE
#if COMP_SINGLE==0
void cp2cpu(dmat** out, const curmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
void cp2cpu(cmat** out, const cucmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
#endif
#else
void cp2cpu(smat** out, const curmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
void cp2cpu(zmat** out, const cucmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
#endif
#define cp2cpu_cell(S, T)						\
void cp2cpu(S##cell **out, const NumCell<T, Gpu> &in, cudaStream_t stream){ \
	if(!in){							\
	    if(*out) S##cellzero(*out);					\
	    return;							\
	}								\
	if(!*out){\
		 *out=S##cellnew(in.Nx(), in.Ny());			\
		 if(in.keywords.length()) (*out)->keywords=strdup(in.keywords.c_str());\
	}\
	for(int i=0; i<in.N(); i++){					\
	    cp2cpu(&(*out)->p[i], in[i], stream);			\
	}								\
    }
cp2cpu_cell(s, Real)
cp2cpu_cell(d, Real)
cp2cpu_cell(c, Comp)
cp2cpu_cell(z, Comp)

void zfarr_push_scale(struct zfarr* ca, int i, const curmat& A, Real scale, cudaStream_t stream){
	X(mat)* tmp=NULL;
	if(scale==1) cp2cpu(&tmp, A, stream);
	else add2cpu(&tmp, 0, A, scale, stream);
	zfarr_push(ca, i, tmp);
	X(free)(tmp);
}

void zfarr_push_scale(struct zfarr *ca, int i, const cucmat &A, Real scale, cudaStream_t stream){
	XC(mat)* tmp=NULL;
	if(scale==1) cp2cpu(&tmp, A, stream);
	else add2cpu(&tmp, 0, A, scale, stream);
	zfarr_push(ca, i, tmp);
	XC(free)(tmp);
}

void zfarr_push_scale(struct zfarr *ca, int i, const curcell &A, Real scale, cudaStream_t stream){
	X(cell)* tmp=NULL;
	if(scale==1) cp2cpu(&tmp, A, stream);
	else add2cpu(&tmp, 0, A, scale, stream);
	zfarr_push(ca, i, tmp);
	X(cellfree)(tmp);
}

void zfarr_push_scale(struct zfarr *ca, int i, const cuccell &A, Real scale, cudaStream_t stream){
	XC(cell)* tmp=NULL;
	if(scale==1) cp2cpu(&tmp, A, stream);
	else add2cpu(&tmp, 0, A, scale, stream);
	zfarr_push(ca, i, tmp);
	XC(cellfree)(tmp);
}

void drawopdamp_gpu(const char* fig, loc_t* loc, const curmat& opd, cudaStream_t stream,
	const dmat* amp, real zlim,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(draw_current(fig, fn)){
		dmat* tmp=NULL;
		cp2cpu(&tmp, opd, stream);
		drawopdamp(fig, loc, tmp, amp, zlim, title, xlabel, ylabel, "%s", fn);
		dfree(tmp);
	}
}
void drawpsf_gpu(const char *fig, curmat &psf, int count, cudaStream_t stream, int zlog, real psfmin,
	const char *title, const char *xlabel, const char *ylabel,
	const char *format, ...){
	format2fn;
	if(draw_current(fig, fn)){
		dmat *psftemp=NULL;
		add2cpu(&psftemp, 0, psf, 1./count, stream);
		draw(fig, (plot_opts){ .image=psftemp, .zlim={psfmin,1}, .zlog=zlog }, title, xlabel, ylabel, "%s", fn);
		dfree(psftemp);
	}
}void curdraw_gpu(const char *fig, curmat &psf, int count, cudaStream_t stream, int zlog,
	const char* title, const char* xlabel, const char* ylabel,
	const char* format, ...){
	format2fn;
	if(draw_current(fig, fn)){
		dmat* psftemp=NULL;
		add2cpu(&psftemp, 0, psf, 1./count, stream);
		draw(fig, (plot_opts){.image=psftemp,.zlog=zlog}, title, xlabel, ylabel, "%s", fn);
		dfree(psftemp);
	}
}
void cucdraw_gpu(const char *fig, cucmat &psf, int count, cudaStream_t stream, int zlog,
	const char *title, const char *xlabel, const char *ylabel,
	const char *format, ...){
	format2fn;
	if(draw_current(fig, fn)){
		cmat *psftemp=NULL;
		cp2cpu(&psftemp, psf, stream);
		dmat *psfreal=NULL;
		cabs22d(&psfreal, 0, psftemp, 1./count);
		draw(fig, (plot_opts){.image=psfreal,.zlog=zlog}, title, xlabel, ylabel, "%s", fn);
		dfree(psfreal);
		cfree(psftemp);
	}
}
/**
   Free data if not referenced or reference is 1.
*/
#undef cudaFree
int mycudaFree(void* pp){
	if(!pp) return 0;
	int tofree=1;
	{
		lock_t tmp(cumemcache.mutex_hash);
		std::map<void*, int>::iterator it=cumemcache.memcount.find(pp);
		if(it!=cumemcache.memcount.end()){
			it->second--;
			tofree=!(it->second);
		}
	}
	if(tofree){
		return cudaFree(pp);
	} else{
		return 0;
	}
}
#undef cudaMalloc
int mycudaMalloc(void** p, size_t size){
	int ans=cudaMalloc(p, size);
	return ans;
}
