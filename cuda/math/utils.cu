/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include <pthread.h>
#include "utils.h"
#include "curmat.h"

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
/**1: allow memory deduplication. Useful during setup for const memory. Set to 1 during setup and 0 during simulation*/
int cuda_dedup=0;
cumemcache_t cumemcache;
/*
The generic one:
template<typename T, template<typename> class Dev>
template<typename S, template<typename> class Dev2>
void NumArray<T, Dev>::Copy(const NumArray<S, Dev2> &in, cudaStream_t stream){}
//Note that you cannot partially spatially a member function.
*/

/**
   Copy map_t to cumap_t. Allow multiple calling to override the data.

   @param[in,out] dest GPU map array. Allocated if empty; otherwise its
   parameters and data are overridden.
   @param[in] source CPU map array. Nothing is done if it is empty.
*/
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

/**
   Convert a host dsp array to GPU sparse array, optionally converted to CSR
   and/or transposed. Index and value arrays are deduplicated if cuda_dedup is
   set.

   @param[in] src_csc Source dsp in CSC. Nothing is done if NULL.
   @param[in] tocsr 0: Keep in CSC. 1: Convert to CSR.
   @param[in] transp 1: transpose.
*/
cusp::cusp(const dsp* src_csc, int tocsr, int transp){
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
	ref->nx=transp?src_csc->ny:src_csc->nx;//dimension of the stored matrix
	ref->ny=transp?src_csc->nx:src_csc->ny;
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
/**
   Concatenate a dspcell into a single sparse matrix and copy it to GPU.

   @param[out] dest GPU sparse matrix. Replaced.
   @param[in] srcc CPU sparse cell array.
   @param[in] tocsr 0: Keep in CSC. 1: Convert to CSR.
*/
void cp2gpu(cusp& dest, const dspcell* srcc, int tocsr){
	dsp* src=dspcell2sp(srcc);
	dest=cusp(src, tocsr);
	dspfree(src);
}
/**
   Copy a dspcell to GPU block by block. dest is allocated if empty.

   @param[in,out] dest GPU sparse cell array. Allocated if empty; each block is replaced.
   @param[in] src CPU sparse cell array. Nothing is done if NULL.
   @param[in] tocsr 0: Keep in CSC. 1: Convert to CSR.
*/
void cp2gpu(cuspcell& dest, const dspcell* src, int tocsr){
	if(!src) return;
	if(!dest){
		dest=cuspcell(src->nx, src->ny);
	}
	for(int i=0; i<src->nx*src->ny; i++){
		dest[i]=cusp(src->p[i], tocsr);
	}
}
/**
   Copy a real dmat or smat to GPU, converting to Real. dest is allocated if
   empty, otherwise its size must match. Zeros dest if src_ is NULL.

   @param[in,out] dest GPU array. Allocated if empty, otherwise copied into.
   @param[in] src_ CPU dmat or smat.
   @param[in] stream Stream to copy in.
*/
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
/**
   Copy a complex cmat or zmat to GPU, converting to Comp. dest is allocated
   if empty, otherwise its size must match. Zeros dest if src_ is NULL.

   @param[in,out] dest GPU array. Allocated if empty, otherwise copied into.
   @param[in] src_ CPU cmat or zmat.
   @param[in] stream Stream to copy in.
*/
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

/**
   Return the index of the current CUDA device. Weak symbol; overridden by
   cuda/sim/cudata.cu when linked.

   @return Index of the current CUDA device.
*/
__attribute__((weak)) int current_gpu(){
	int igpu;
	cudaGetDevice(&igpu);
	return igpu;
}
/**
   Scale and add host arrays: p1=p1*alpha+p2*beta. Specialized below for
   complex arrays of different precision.

   @param[in,out] p1 Array of length n.
   @param[in] alpha Scaling factor of p1.
   @param[in] p2 Array of length n.
   @param[in] beta Scaling factor of p2.
   @param[in] n Number of elements.
*/
template <typename T, typename R, typename S>
void scale_add(T* p1, R alpha, const S* p2, R beta, long n){
	for(long i=0; i<n; i++){
		p1[i]=p1[i]*alpha+p2[i]*beta;
	}
}
/**
   Specialization of scale_add() for double2 p1 and float2 p2.

   @param[in,out] p1 Array of length n.
   @param[in] alpha Scaling factor of p1.
   @param[in] p2 Array of length n.
   @param[in] beta Scaling factor of p2.
   @param[in] n Number of elements.
*/
template <>
void scale_add<double2, double, float2>(double2* p1, double alpha, const float2* p2, double beta, long n){
	for(long i=0; i<n; i++){
		p1[i].x=p1[i].x*alpha+p2[i].x*beta;
		p1[i].y=p1[i].y*alpha+p2[i].y*beta;
	}
}
/**
   Specialization of scale_add() for float2 p1 and double2 p2.

   @param[in,out] p1 Array of length n.
   @param[in] alpha Scaling factor of p1.
   @param[in] p2 Array of length n.
   @param[in] beta Scaling factor of p2.
   @param[in] n Number of elements.
*/
template <>
void scale_add<float2, float, double2>(float2* p1, float alpha, const double2* p2, float beta, long n){
	for(long i=0; i<n; i++){
		p1[i].x=p1[i].x*alpha+p2[i].x*beta;
		p1[i].y=p1[i].y*alpha+p2[i].y*beta;
	}
}
/**
   Copy device array to host, convert type and accumulate:
   dest = alpha * dest + beta *src. Synchronizes the stream.

   @param[in,out] dest Pointer to host array. Allocated (uninitialized) if NULL.
   @param[in] alpha Scaling factor of dest.
   @param[in] src Device array of length n.
   @param[in] beta Scaling factor of src.
   @param[in] n Number of elements.
   @param[in] stream Stream to copy in.
   @param[in,out] mutex If not NULL, locked during the accumulation.
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
/**
   Copy GPU array to CPU and accumulate: *out=*out*alpha+in*beta. *out is
   allocated if NULL; if in is empty, *out is only scaled by alpha. mutex, if
   not NULL, protects the accumulation. Synchronizes the stream.
*/
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
/*
	We provide add2cpu and cp2cpu for all CPU types.
*/
add2cpu_mat(d, real, Real)
add2cpu_mat(c, real, Comp)
add2cpu_mat(s, float, Real)
add2cpu_mat(z, float, Comp)
/**
   Cell version of add2cpu: *out=*out*alpha+in*beta for each block.
*/
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

/**
   Copy GPU array to CPU of the same precision. *out is allocated if NULL,
   otherwise its size must match; zeroed if in is empty. Keywords are copied.
   Synchronizes the stream.
*/
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
cp2cpu_same(dmat, dzero, dnew, real)
cp2cpu_same(cmat, czero, cnew, dcomplex)
cp2cpu_same(smat, szero, snew, float)
cp2cpu_same(zmat, zzero, znew, fcomplex)
#if CUDA_DOUBLE
#if CPU_SINGLE
/**
   Copy GPU array to CPU with precision conversion (via add2cpu).

   @param[in,out] out CPU array. Allocated if NULL, otherwise overwritten.
   @param[in] in GPU array.
   @param[in] stream Stream to copy in.
*/
void cp2cpu(dmat** out, const curmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
/**
   Copy GPU array to CPU with precision conversion (via add2cpu).

   @param[in,out] out CPU array. Allocated if NULL, otherwise overwritten.
   @param[in] in GPU array.
   @param[in] stream Stream to copy in.
*/
void cp2cpu(cmat** out, const cucmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
#endif
/**
   Copy GPU array to CPU with precision conversion (via add2cpu).

   @param[in,out] out CPU array. Allocated if NULL, otherwise overwritten.
   @param[in] in GPU array.
   @param[in] stream Stream to copy in.
*/
void cp2cpu(smat** out, const curmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
/**
   Copy GPU array to CPU with precision conversion (via add2cpu).

   @param[in,out] out CPU array. Allocated if NULL, otherwise overwritten.
   @param[in] in GPU array.
   @param[in] stream Stream to copy in.
*/
void cp2cpu(zmat** out, const cucmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
#else
#if CPU_SINGLE==0
/**
   Copy GPU array to CPU with precision conversion (via add2cpu).

   @param[in,out] out CPU array. Allocated if NULL, otherwise overwritten.
   @param[in] in GPU array.
   @param[in] stream Stream to copy in.
*/
void cp2cpu(dmat** out, const curmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
/**
   Copy GPU array to CPU with precision conversion (via add2cpu).

   @param[in,out] out CPU array. Allocated if NULL, otherwise overwritten.
   @param[in] in GPU array.
   @param[in] stream Stream to copy in.
*/
void cp2cpu(cmat** out, const cucmat& in, cudaStream_t stream){
	add2cpu(out, 0, in, 1, stream, 0);
}
#endif
#endif
/**
   Copy GPU cell array to CPU block by block. *out is allocated if NULL and
   zeroed if in is empty.
*/
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

/**
   Copy A to CPU, scale by scale, and save it as the i'th entry of ca.

   @param[in,out] ca File array to save to.
   @param[in] i Index of the entry.
   @param[in] A GPU array.
   @param[in] scale Scaling factor.
   @param[in] stream Stream to copy in.
*/
void zfarr_push_scale(struct zfarr* ca, int i, const curmat& A, Real scale, cudaStream_t stream){
	X(mat)* tmp=NULL;
	if(scale==1) cp2cpu(&tmp, A, stream);
	else add2cpu(&tmp, 0, A, scale, stream);
	zfarr_push(ca, i, tmp);
	X(free)(tmp);
}

/**
   Copy A to CPU, scale by scale, and save it as the i'th entry of ca.

   @param[in,out] ca File array to save to.
   @param[in] i Index of the entry.
   @param[in] A GPU array.
   @param[in] scale Scaling factor.
   @param[in] stream Stream to copy in.
*/
void zfarr_push_scale(struct zfarr *ca, int i, const cucmat &A, Real scale, cudaStream_t stream){
	XC(mat)* tmp=NULL;
	if(scale==1) cp2cpu(&tmp, A, stream);
	else add2cpu(&tmp, 0, A, scale, stream);
	zfarr_push(ca, i, tmp);
	XC(free)(tmp);
}

/**
   Copy A to CPU, scale by scale, and save it as the i'th entry of ca.

   @param[in,out] ca File array to save to.
   @param[in] i Index of the entry.
   @param[in] A GPU array.
   @param[in] scale Scaling factor.
   @param[in] stream Stream to copy in.
*/
void zfarr_push_scale(struct zfarr *ca, int i, const curcell &A, Real scale, cudaStream_t stream){
	X(cell)* tmp=NULL;
	if(scale==1) cp2cpu(&tmp, A, stream);
	else add2cpu(&tmp, 0, A, scale, stream);
	zfarr_push(ca, i, tmp);
	X(cellfree)(tmp);
}

/**
   Copy A to CPU, scale by scale, and save it as the i'th entry of ca.

   @param[in,out] ca File array to save to.
   @param[in] i Index of the entry.
   @param[in] A GPU array.
   @param[in] scale Scaling factor.
   @param[in] stream Stream to copy in.
*/
void zfarr_push_scale(struct zfarr *ca, int i, const cuccell &A, Real scale, cudaStream_t stream){
	XC(cell)* tmp=NULL;
	if(scale==1) cp2cpu(&tmp, A, stream);
	else add2cpu(&tmp, 0, A, scale, stream);
	zfarr_push(ca, i, tmp);
	XC(cellfree)(tmp);
}

/**
   Plot GPU opd defined on loc with amplitude amp, if the figure is currently
   drawn. See drawopdamp().

   @param[in] fig Figure category.
   @param[in] loc Coordinate of the points.
   @param[in] opd GPU OPD defined on loc.
   @param[in] stream Stream to copy in.
   @param[in] amp Amplitude defined on loc.
   @param[in] zlim Color range limit.
   @param[in] title Plot title.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the figure name, followed by arguments.
*/
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
/**
   Plot GPU psf divided by count as an image with color range [psfmin, 1], if
   the figure is currently drawn.

   @param[in] fig Figure category.
   @param[in] psf GPU PSF.
   @param[in] count Number of accumulations to normalize by.
   @param[in] stream Stream to copy in.
   @param[in] zlog Whether to plot in log scale.
   @param[in] psfmin Lower limit of the color range.
   @param[in] title Plot title.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the figure name, followed by arguments.
*/
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
}
/**
   Plot GPU real array psf divided by count as an image, if the figure is
   currently drawn.

   @param[in] fig Figure category.
   @param[in] psf GPU real array.
   @param[in] count Number of accumulations to normalize by.
   @param[in] stream Stream to copy in.
   @param[in] zlog Whether to plot in log scale.
   @param[in] title Plot title.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the figure name, followed by arguments.
*/
void curdraw_gpu(const char *fig, curmat &psf, int count, cudaStream_t stream, int zlog,
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
/**
   Plot abs(psf)^2/count of GPU complex array psf as an image, if the figure
   is currently drawn.

   @param[in] fig Figure category.
   @param[in] psf GPU complex array.
   @param[in] count Number of accumulations to normalize by.
   @param[in] stream Stream to copy in.
   @param[in] zlog Whether to plot in log scale.
   @param[in] title Plot title.
   @param[in] xlabel x axis label.
   @param[in] ylabel y axis label.
   @param[in] format printf style format of the figure name, followed by arguments.
*/
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
#undef cudaFree
/**
   Unreference deduplicated memory. Free data if not deduplicated or the
   reference count drops to 0.

   @param[in] pp Device memory. Nothing is done if NULL.
   @return CUDA error code of cudaFree, or 0 if not freed.
*/
int mycudaFree(void* pp){
	if(!pp) return 0;
	int tofree=1;
	{
		lock_t tmp(cumemcache.mutex_hash);
		std::map<void*, uint32_t[2]>::iterator it=cumemcache.memcount.find(pp);
		if(it!=cumemcache.memcount.end()){
			it->second[0]--;
			if((tofree=!(it->second[0]))){
				//dbg("memhash.erase %x\n", it->second[1]);
				cumemcache.memhash.erase(it->second[1]);
				cumemcache.memcount.erase(it);
			}
		}
	}
	if(tofree){
		return cudaFree(pp);
	} else{
		return 0;
	}
}
#undef cudaMalloc
/**
   Wrapper of cudaMalloc.

   @param[out] p Allocated device memory. The caller frees it with cudaFree.
   @param[in] size Number of bytes.
   @return CUDA error code.
*/
int mycudaMalloc(void** p, size_t size){
	int ans=cudaMalloc(p, size);
	return ans;
}
