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
#include "utils.h"
#include "curmat.h"
#include "cucmat.h"
#include <pthread.h>
const char *cufft_str[]={
    "success", 
    "invalid plan",
    "allocation failed",
    "",
    "invalid value",
    "internal errlr",
    "exec failed (error elsewhere caused cufft to fail)",
    "setup failed"
    "invalid size"
};
#ifndef I
#define I (__extension__ 1.0iF)
#endif
static cusparseMatDescr_t spdesc=NULL;
#if CUDA_VERSION < 4010
pthread_mutex_t cufft_mutex=PTHREAD_MUTEX_INITIALIZER;
#endif
int cuda_dedup=0;//1: allow memory deduplication. Useful during setup for const memory.
static __attribute((constructor)) void init(){
    DO(cusparseCreateMatDescr(&spdesc));
    cusparseSetMatType(spdesc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(spdesc, CUSPARSE_INDEX_BASE_ZERO);
}

/**
   Copy map_t to cumap_t. if type==1, use cudaArray, otherwise use Real
   array. Allow multiple calling to override the data.  */
void cp2gpu(cumap_t **dest0, map_t **source, int nps){
    if(nps==0) return;
    if(!*dest0){
	*dest0=new cumap_t[nps];
    }
    for(int ips=0; ips<nps; ips++){
	(*dest0)[ips].init(source[ips]);
    }
    cumap_t *dest=*dest0;
    for(int ips=0; ips<nps; ips++){
	cp2gpu(&dest[ips].p, (dmat*)source[ips]);
    }
}

/**
   Copy map_t to cumap_t. if type==1, use cudaArray, otherwise use Real
   array. Allow multiple calling to override the data.  */
void gpu2gpu(cumap_t **dest0, cumap_t *source, int nps){
    if(nps==0) return;
    if(!*dest0){
	*dest0=new cumap_t[nps];
	memcpy(*dest0, source, sizeof(cumap_t)*nps);
	for(int ips=0; ips<nps; ips++){
	    (*dest0)[ips].p=0;
	    (*dest0)[ips].cubic_cc=0;
	}	
    }
    cumap_t *dest=*dest0;
    for(int ips=0; ips<nps; ips++){
	if(source[ips].p){
	    gpu2gpu(&dest[ips].p, source[ips].p);
	}
	if(source[ips].cubic_cc){
	    gpu2gpu(&dest[ips].cubic_cc, source[ips].cubic_cc);
	}
    }
}
/*
  Convert a host dsp array to GPU sprase array.
*/
cusp::cusp(const dsp *src_csc, int tocsr)
    :p(NULL),i(NULL),x(NULL),nx(0),ny(0),nzmax(0),type(SP_CSC){
    if(!src_csc) return;
    dsp *src=const_cast<dsp*>(src_csc);
    if(tocsr){
	type=SP_CSR;
	src=sptrans(src_csc);
    }else{
	type=SP_CSC;
    }
    nx=src_csc->m;
    ny=src_csc->n;
    nzmax=src->nzmax;
    p=NULL; i=NULL; x=NULL;
    cp2gpu(&p, src->p, src->n+1, 1);
    cp2gpu(&i, src->i, src->nzmax, 1);
    cp2gpu(&x, src->x, src->nzmax, 1);
    if(tocsr){
	spfree(src);
    }
    nref=new int[1];
    nref[0]=1;
}
void cp2gpu(cusp **dest0, const spcell *srcc, int tocsr){
    dsp *src=spcell2sp(srcc);
    *dest0=new cusp(src, tocsr);
    spfree(src);
}
void cp2gpu(cuspcell **dest0, const spcell *src, int tocsr){
    if(!src) return;
    if(!*dest0){
	*dest0=cuspcellnew(src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	(*dest0)->p[i]=new cusp(src->p[i], tocsr);
    }
}
static const char *scsrmv_err[]={
    "Success",
    "Not initialized",
    "Allocation failed",
    "Invalid value",
    "Archtecture mismatch",
    "Mapping error",
    "Execution failed",
    "Internal error",
    "Matrix type not supported"
};

/*
  y=A*x where A is sparse. x, y are vectors. Slow for GS0.
*/

void cuspmul(Real *y, cusp *A, const Real *x, int ncolvec, char trans, Real alpha, cusparseHandle_t handle){
    cusparseOperation_t opr;
    int istrans=(trans=='t' || trans==1);
    if(A->type==SP_CSC){
	istrans=!istrans;
    }
    if(istrans){
	opr=CUSPARSE_OPERATION_TRANSPOSE;
    }else{
	opr=CUSPARSE_OPERATION_NON_TRANSPOSE;
    }
    int ncol=0, nrow=0;
    switch(A->type){
    case SP_CSR:
	nrow=A->nx; ncol=A->ny; break;
    case SP_CSC:
	nrow=A->ny; ncol=A->nx; break;
    default:
	error("Invalid type");
    }
    int status;
    Real one=1.f;
    if(ncolvec==1){
	status=CUSP(csrmv)(handle, opr,
			   nrow, ncol, A->nzmax, &alpha, spdesc,
			   A->x, A->p, A->i, x, &one, y);
    }else{
	int nlead=istrans?nrow:ncol;
	status=CUSP(csrmm)(handle, opr,
			   nrow, ncolvec, ncol, A->nzmax, &alpha, spdesc,
			   A->x, A->p, A->i, x, nlead, &one, y, nlead);
    }
    if(status!=0){
	error("cusparseScsrmv(m) failed with status '%s'\n", scsrmv_err[status]);
    }
}

/**
   Convert a source loc_t to device memory.
*/
void cp2gpu(Real (* restrict *dest)[2], const loc_t *src){
    Real (*tmp)[2]=(Real(*)[2])malloc(src->nloc*2*sizeof(Real));
    for(int iloc=0; iloc<src->nloc; iloc++){
	tmp[iloc][0]=(Real)src->locx[iloc];
	tmp[iloc][1]=(Real)src->locy[iloc];
    }
    cp2gpu((Real**)dest, (Real*)tmp, src->nloc*2, 1);
    free(tmp);
}

/**
   Convert dcell to curcell
*/
void cp2gpu(curcell *restrict *dest, const dcell *src){
    if(!src) {
	dzero(*dest);
	return;
    }
    if(!*dest) {
	long nc=src->nx*src->ny;
	long nx[nc];
	long ny[nc];
	for(long i=0; i<nc; i++){
	    if(src->p[i]){
		nx[i]=src->p[i]->nx;
		ny[i]=src->p[i]->ny;
	    }else{
		nx[i]=0;
		ny[i]=0;
	    }
	}
	*dest=curcellnew(src->nx, src->ny, nx, ny);
    }else if((*dest)->nx!=src->nx || (*dest)->ny!=src->ny){
	error("Mismatch: %ldx%ld vs %ldx%ld\n", 
	      (*dest)->nx, (*dest)->ny, src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	cp2gpu(&(*dest)->p[i], src->p[i]);
    }
}
/**
   Convert dcell to curcell
*/
void cp2gpu(cuccell *restrict *dest, const ccell *src){
    if(!src) {
	dzero(*dest);
	return;
    }
    if(!*dest) {
	long nc=src->nx*src->ny;
	long nx[nc];
	long ny[nc];
	for(long i=0; i<nc; i++){
	    if(src->p[i]){
		nx[i]=src->p[i]->nx;
		ny[i]=src->p[i]->ny;
	    }else{
		nx[i]=0;
		ny[i]=0;
	    }
	}
	*dest=cuccellnew(src->nx, src->ny, nx, ny);
    }else if((*dest)->nx!=src->nx || (*dest)->ny!=src->ny){
	error("Mismatch: %ldx%ld vs %ldx%ld\n", 
	      (*dest)->nx, (*dest)->ny, src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	cp2gpu(&(*dest)->p[i], src->p[i]);
    }
}

/*
  Write Real on gpu to file
*/ 
void gpu_write(const Real *p, int nx, int ny, const char *format, ...){
    format2fn;
    Real *tmp=(Real*)malloc(nx*ny*sizeof(Real));
    cudaMemcpy(tmp, p, nx*ny*sizeof(Real), cudaMemcpyDeviceToHost);
#if CUDA_DOUBLE
    writedbl(tmp,nx,ny,"%s",fn);
#else
    writeflt(tmp,nx,ny,"%s",fn);
#endif
    free(tmp);
}

/*
  Write Real on gpu to file
*/
void gpu_write(const Comp *p, int nx, int ny, const char *format, ...){
    format2fn;
    Comp *tmp=(Comp*)malloc(nx*ny*sizeof(Comp));
    cudaMemcpy(tmp, p, nx*ny*sizeof(Comp), cudaMemcpyDeviceToHost);
#if CUDA_DOUBLE
    writecmp((dcomplex*)tmp,nx,ny,"%s",fn);
#else
    writefcmp((fcomplex*)tmp,nx,ny,"%s",fn);
#endif
    free(tmp);
}
/*
  Write Real on gpu to file
*/
void gpu_write(const int *p, int nx, int ny, const char *format, ...){
    format2fn;
    int *tmp=(int*)malloc(nx*ny*sizeof(int));
    cudaMemcpy(tmp, p, nx*ny*sizeof(int), cudaMemcpyDeviceToHost);
    writeint(tmp,nx,ny,"%s",fn);
    free(tmp);
}
template <typename T, typename R, typename S>
    inline void scale_add(T *p1, R alpha, S *p2, R beta, long n){
    for(long i=0; i<n; i++){
	p1[i]=p1[i]*alpha+p2[i]*beta;
    }
}
template <>
inline void scale_add<dcomplex, double, Comp>(dcomplex *p1, double alpha, Comp *p2, double beta, long n){
    for(long i=0; i<n; i++){
	p1[i]=p1[i]*alpha+(p2[i].x+I*p2[i].y)*beta;
    }
}
template <>
inline void scale_add<fcomplex, float, Comp>(fcomplex *p1, float alpha, Comp *p2, float beta, long n){
    for(long i=0; i<n; i++){
	p1[i]=p1[i]*alpha+(p2[i].x+I*p2[i].y)*beta;
    }
}
/**
   Convert device (Real) array and add to host double.
   dest = alpha * dest + beta *src;
*/
template <typename R, typename T, typename S>
static void add2cpu(T * restrict *dest, R alpha, S *src, R beta, long n, 
		    cudaStream_t stream, pthread_mutex_t *mutex){
    extern int cuda_dedup;
    S *tmp=0;
    if(!cuda_dedup && cudata->memcache->count((void*)src)){
	tmp=(S*)(*cudata->memcache)[(void*)src];
    }else{
	tmp=(S*)malloc(n*sizeof(S));
	if(!cuda_dedup){
	    (*cudata->memcache)[(void*)src]=(void*)tmp;
	}
    }
    DO(cudaMemcpyAsync(tmp, src, n*sizeof(S), cudaMemcpyDeviceToHost, stream));
    if(!*dest){
	*dest=(T*)malloc(sizeof(T)*n);
    }
    T *restrict p=*dest;
    if(mutex) LOCK(*mutex);
    scale_add(p, alpha, tmp, beta, n);
    if(mutex) UNLOCK(*mutex);
    if(cuda_dedup) {
	free(tmp);
    }
}
#define add2cpu_mat(D, double, Comp)					\
void add2cpu(D##mat **out, double alpha, const cumat<Comp> *in, double beta,	\
	     cudaStream_t stream, pthread_mutex_t *mutex){		\
    if(!in){								\
	if(*out) D##scale(*out, alpha);					\
	return;								\
    }									\
    if(!*out) {								\
	*out=D##new(in->nx, in->ny);					\
    }else{								\
	assert((*out)->nx*(*out)->ny==in->nx*in->ny);			\
    }									\
    add2cpu(&(*out)->p, alpha, in->p, beta, in->nx*in->ny, stream, mutex);\
}
add2cpu_mat(s, float, Real)
add2cpu_mat(d, double,Real)
add2cpu_mat(z, float, Comp)
add2cpu_mat(c, double,Comp)

#define add2cpu_cell(D, double, curcell)				\
void add2cpu(D##cell **out, double alpha, const curcell *in, double beta, \
	     cudaStream_t stream, pthread_mutex_t *mutex){		\
    if(!in){								\
	if(*out) D##cellscale(*out, alpha);				\
	return;								\
    }									\
    if(!*out) {								\
	*out=D##cellnew(in->nx, in->ny);				\
    }else{								\
	assert((*out)->nx*(*out)->ny==in->nx*in->ny);			\
    }									\
    for(int i=0; i<in->nx*in->ny; i++){					\
	add2cpu(&(*out)->p[i], alpha, in->p[i], beta, stream, mutex);	\
    }									\
}
add2cpu_cell(d, double,curcell)
add2cpu_cell(s, float, curcell)
add2cpu_cell(c, double,cuccell)
add2cpu_cell(z, float, cuccell)
#define cp2cpu_same(dmat,dzero,dnew,double)				\
    void cp2cpu(dmat **out, const cumat<double> *in, cudaStream_t stream){ \
	if(!in) {							\
	if(*out) dzero(*out);						\
	return;								\
    }									\
    if(!*out) *out=dnew(in->nx, in->ny);				\
    DO(cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(double),	\
		       cudaMemcpyDeviceToHost, stream));		\
    if(in->header) (*out)->header=strdup(in->header);			\
    }

cp2cpu_same(dmat,dzero,dnew,double)
cp2cpu_same(cmat,czero,cnew,double2)
cp2cpu_same(smat,szero,snew,float)
cp2cpu_same(zmat,zzero,znew,float2)
#if CUDA_DOUBLE == 0
void cp2cpu(dmat **out, const curmat *in, cudaStream_t stream){
    add2cpu(out, 0, in, 1, stream, 0);
}
void cp2cpu(cmat **out, const cucmat *in, cudaStream_t stream){
    add2cpu(out, 0, in, 1, stream, 0);
}
#else
void cp2cpu(smat **out, const curmat *in, cudaStream_t stream){
    add2cpu(out, 0, in, 1, stream, 0);
}
void cp2cpu(zmat **out, const cucmat *in, cudaStream_t stream){
    add2cpu(out, 0, in, 1, stream, 0);
}
#endif
#define cp2cpu_cell(S, float)						\
    void cp2cpu(S##cell **out, const cucell<float> *in, cudaStream_t stream){	\
	if(!in){							\
	    if(*out) S##cellzero(*out);					\
	    return;							\
	}								\
	if(!*out) *out=S##cellnew(in->nx, in->ny);			\
	for(int i=0; i<in->nx*in->ny; i++){				\
	    cp2cpu(&(*out)->p[i], in->p[i], stream);			\
	}								\
    }
cp2cpu_cell(s, Real)
cp2cpu_cell(d, Real)
cp2cpu_cell(c, Comp)
cp2cpu_cell(z, Comp)

void cellarr_cur(struct cellarr *ca, int i, const curmat *A, cudaStream_t stream){
    X(mat) *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
#if CUDA_DOUBLE==1
    cellarr_dmat(ca, i, tmp);
#else
    cellarr_smat(ca, i, tmp);
#endif
    X(free)(tmp);
}

void cellarr_cuc(struct cellarr *ca, int i, const cucmat *A, cudaStream_t stream){
    C(mat) *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
#if CUDA_DOUBLE==1
    cellarr_cmat(ca, i, tmp);
#else
    cellarr_zmat(ca, i, tmp);
#endif
    C(free)(tmp);
}

void cellarr_curcell(struct cellarr *ca, int i, const curcell *A, cudaStream_t stream){
    X(cell) *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
#if CUDA_DOUBLE==1
    cellarr_dcell(ca, i, tmp);
#else
    cellarr_scell(ca, i, tmp);
#endif
    X(cellfree)(tmp);
}

void cellarr_cuccell(struct cellarr *ca, int i, const cuccell *A, cudaStream_t stream){
    C(cell) *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
#if CUDA_DOUBLE==1
    cellarr_ccell(ca, i, tmp);
#else
    cellarr_zcell(ca, i, tmp);
#endif
    C(cellfree)(tmp);
}
