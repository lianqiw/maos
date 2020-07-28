/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
void cp2gpu(cumapcell &dest, const mapcell *source){
    if(source->nx==0) return;
    if(!dest){
	dest=cumapcell(source->nx,1);
    }
    //initalize or override parameters.
    for(int ips=0; ips<source->nx; ips++){
	dest[ips]=(source->p[ips]);
    }
    for(int ips=0; ips<source->nx; ips++){
	cp2gpu(dest[ips].p, (const dmat*)source->p[ips]);
    }
}

/*
  Convert a host dsp array to GPU sprase array.g
*/
cusp::cusp(const dsp *src_csc, /**<Source dsp in CSC*/
	   int tocsr,        /**<0: Keep in CSC. 1: Convert to CSR */
	   int transp        /**<1: transpost*/
    )
    :p(NULL),i(NULL),x(NULL),nx(0),ny(0),nzmax(0),type(tocsr?SP_CSR:SP_CSC),nref(0){
    if(!src_csc) return;
    dsp *src_trans=0;
    const dsp *src=0;
    if(tocsr!=transp){
	src_trans=dsptrans(src_csc);
	src=src_trans;
    }else{
	src=src_csc;
    }
    nx=src_csc->nx;
    ny=src_csc->ny;
    nzmax=src->nzmax;
    cp2gpu(&p, src->p, src->ny+1, 1);
    cp2gpu(&i, src->i, src->nzmax, 1);
    cp2gpu(&x, src->x, src->nzmax, 1);
    nref=mymalloc(1, int);
    nref[0]=1;
    if(src_trans){
	dspfree(src_trans);
    }
}
void cp2gpu(cusp &dest, const dspcell *srcc, int tocsr){
    dsp *src=dspcell2sp(srcc);
    dest=cusp(src, tocsr);
    dspfree(src);
}
void cp2gpu(cuspcell &dest, const dspcell *src, int tocsr){
    if(!src) return;
    if(!dest){
	dest=cuspcell(src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	dest[i]=cusp(src->p[i], tocsr);
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

void cuspmul(Real *y, const cusp &A, const Real *x, int ncolvec, char trans, Real alpha, stream_t &stream){
    cusparseOperation_t opr;
    int istrans=(trans=='t' || trans==1);
    if(A.Type()==SP_CSC){
	istrans=!istrans;
    }
    if(istrans){
	opr=CUSPARSE_OPERATION_TRANSPOSE;
    }else{
	opr=CUSPARSE_OPERATION_NON_TRANSPOSE;
    }
    int ncol=0, nrow=0;
    switch(A.Type()){
    case SP_CSR:
	nrow=A.Nx(); ncol=A.Ny(); break;
    case SP_CSC:
	nrow=A.Ny(); ncol=A.Nx(); break;
    default:
	error("Invalid type: %d", A.Type());
    }
    int status;
    Real one=1.f;
    if(ncolvec==1){
	status=CUSP(csrmv)(stream.sparse(), opr,
			   nrow, ncol, A.Nzmax(), &alpha, spdesc,
			   A.Px(), A.Pp(), A.Pi(), x, &one, y);
    }else{
	int nlead=istrans?nrow:ncol;
	status=CUSP(csrmm)(stream.sparse(), opr,
			   nrow, ncolvec, ncol, A.Nzmax(), &alpha, spdesc,
			   A.Px(), A.Pp(), A.Pi(), x, nlead, &one, y, nlead);
    }
    if(status!=0){
	error("cusparseScsrmv(m) failed with status '%s'\n", scsrmv_err[status]);
    }
}

/**
   Convert a source loc_t to device memory. row vector is used.
*/
void cp2gpu(curmat &dest, const loc_t *src){
    Real2 *tmp=(Real2*)malloc(src->nloc*sizeof(Real2));
    for(int iloc=0; iloc<src->nloc; iloc++){
	tmp[iloc][0]=(Real)src->locx[iloc];
	tmp[iloc][1]=(Real)src->locy[iloc];
    }
    cp2gpu(dest, (Real*)tmp, 2, src->nloc);
    free(tmp);
}

/**
   Convert dcell to curcell
*/
void cp2gpu(curcell &dest, const dcell *src){
    if(!src) {
	dest.zero();
	return;
    }
    if(!dest) {
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
	dest=curcell(src->nx, src->ny, nx, ny);
    }else if(dest.Nx()!=src->nx || dest.Ny()!=src->ny){
	error("Mismatch: %ldx%ld vs %ldx%ld\n", 
	      dest.Nx(), dest.Ny(), src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	cp2gpu(dest[i], src->p[i]);
    }
}
/**
   Convert dcell to curcell
*/
void cp2gpu(cuccell &dest, const ccell *src){
    if(!src) {
	dest.zero();
	return;
    }
    if(!dest) {
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
	dest=cuccell(src->nx, src->ny, nx, ny);
    }else if(dest.Nx()!=src->nx || dest.Ny()!=src->ny){
	error("Mismatch: %ldx%ld vs %ldx%ld\n", 
	      dest.Nx(), dest.Ny(), src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	cp2gpu(dest[i], src->p[i]);
    }
}

/*
  Write Real on gpu to file
*/ 
void gpu_write(const Real *p, int nx, int ny, const char *format, ...){
    format2fn;
    Real *tmp=(Real*)malloc(nx*ny*sizeof(Real));
    cudaMemcpy(tmp, p, nx*ny*sizeof(Real), cudaMemcpyDeviceToHost);
    writearr(fn, 1, sizeof(Real), MCU_REAL, NULL, p, nx, ny);
    free(tmp);
}

/*
  Write Real on gpu to file
*/
void gpu_write(const Comp *p, int nx, int ny, const char *format, ...){
    format2fn;
    Comp *tmp=(Comp*)malloc(nx*ny*sizeof(Comp));
    cudaMemcpy(tmp, p, nx*ny*sizeof(Comp), cudaMemcpyDeviceToHost);
    writearr(fn, 1, sizeof(Comp), MCU_COMP, NULL, p, nx, ny);
    free(tmp);
}
/*
  Write Real on gpu to file
*/
void gpu_write(const int *p, int nx, int ny, const char *format, ...){
    format2fn;
    int *tmp=(int*)malloc(nx*ny*sizeof(int));
    cudaMemcpy(tmp, p, nx*ny*sizeof(int), cudaMemcpyDeviceToHost);
    writearr(fn, 1, sizeof(int), M_INT32, NULL, tmp, nx, ny);
    free(tmp);
}
template <typename T, typename R, typename S>
void scale_add(T *p1, R alpha, const S *p2, R beta, long n){
    for(long i=0; i<n; i++){
	p1[i]=p1[i]*alpha+p2[i]*beta;
    }
}
template <>
void scale_add<double2, double, Comp>(double2 *p1, double alpha, const Comp *p2, double beta, long n){
    for(long i=0; i<n; i++){
	p1[i].x=p1[i].x*alpha+p2[i].x*beta;
	p1[i].y=p1[i].y*alpha+p2[i].y*beta;
    }
}
template <>
 void scale_add<float2, float, Comp>(float2 *p1, float alpha, const Comp *p2, float beta, long n){
    for(long i=0; i<n; i++){
	p1[i].x=p1[i].x*alpha+p2[i].x*beta;
	p1[i].y=p1[i].y*alpha+p2[i].y*beta;
    }
}
/**
   Convert device (Real) array and add to host double.
   dest = alpha * dest + beta *src;
*/
template <typename R, typename T, typename S>
static void add2cpu(T * restrict *dest, R alpha, const S *src, R beta, long n, 
		    cudaStream_t stream, pthread_mutex_t *mutex){
    S *tmp=0;
    tmp=(S*)malloc(n*sizeof(S));
    CUDA_SYNC_STREAM;
    DO(cudaMemcpy(tmp, src, n*sizeof(S), cudaMemcpyDeviceToHost));
    if(!*dest){
	*dest=(T*)malloc(sizeof(T)*n);
    }
    T *restrict p=*dest;
    if(mutex) LOCK(*mutex);
    scale_add(p, alpha, tmp, beta, n);
    if(mutex) UNLOCK(*mutex);
    free(tmp);
}
#define add2cpu_mat(D, T, C)						\
    void add2cpu(D##mat **out, T alpha, const Array<C, Gpu> &in, T beta, \
		 cudaStream_t stream, pthread_mutex_t *mutex){		\
	if(!in){							\
	    if(*out) D##scale(*out, alpha);				\
	    return;							\
	}								\
	if(!*out) {							\
	    *out=D##new(in.Nx(), in.Ny());				\
	}else{								\
	    assert((*out)->nx*(*out)->ny==in.N());			\
	}								\
	add2cpu(&(*out)->p, alpha, in(), beta, in.N(), stream, mutex); \
    }

add2cpu_mat(s, float,Real)
add2cpu_mat(z, float,Comp)
#if COMP_SINGLE==0
add2cpu_mat(d, real, Real)
add2cpu_mat(c, real, Comp)
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
	}else{								\
	    assert((*out)->nx*(*out)->ny==in.N());			\
	}								\
	for(int i=0; i<in.N(); i++){					\
	    add2cpu((*out)->p+i, alpha, in[i], beta, stream, mutex);	\
	}								\
    }
add2cpu_cell(d, real, curcell)
add2cpu_cell(s, float,curcell)
add2cpu_cell(c, real, cuccell)
add2cpu_cell(z, float,cuccell)
#define cp2cpu_same(dmat,dzero,dnew,T)				\
    void cp2cpu(dmat **out, const Array<T, Gpu> &in, cudaStream_t stream){ \
	if(!in) {							\
	    if(*out) dzero(*out);					\
	    return;							\
	}								\
	if(!*out) *out=dnew(in.Nx(), in.Ny());				\
	dmat *pout=*out;						\
	CUDA_SYNC_STREAM;						\
	DO(cudaMemcpy(pout->p, in(), in.N()*sizeof(T),			\
		      cudaMemcpyDeviceToHost));				\
	if(pout->header) free(pout->header);				\
	if(in.header.length()) pout->header=strdup(in.header.c_str());	\
    }
#if COMP_SINGLE==0
cp2cpu_same(dmat,dzero,dnew,double)
cp2cpu_same(cmat,czero,cnew,double2)
#endif
cp2cpu_same(smat,szero,snew,float)
cp2cpu_same(zmat,zzero,znew,float2)
#if ! CUDA_DOUBLE 
#if COMP_SINGLE==0
void cp2cpu(dmat **out, const curmat &in, cudaStream_t stream){
    add2cpu(out, 0, in, 1, stream, 0);
}
void cp2cpu(cmat **out, const cucmat &in, cudaStream_t stream){
    add2cpu(out, 0, in, 1, stream, 0);
}
#endif
#else
void cp2cpu(smat **out, const curmat &in, cudaStream_t stream){
    add2cpu(out, 0, in, 1, stream, 0);
}
void cp2cpu(zmat **out, const cucmat &in, cudaStream_t stream){
    add2cpu(out, 0, in, 1, stream, 0);
}
#endif
#define cp2cpu_cell(S, T)						\
    void cp2cpu(S##cell **out, const Cell<T, Gpu> &in, cudaStream_t stream){ \
	if(!in){							\
	    if(*out) S##cellzero(*out);					\
	    return;							\
	}								\
	if(!*out) *out=S##cellnew(in.Nx(), in.Ny());			\
	for(int i=0; i<in.N(); i++){					\
	    cp2cpu(&(*out)->p[i], in[i], stream);			\
	}								\
    }
cp2cpu_cell(s, Real)
cp2cpu_cell(d, Real)
cp2cpu_cell(c, Comp)
cp2cpu_cell(z, Comp)

void zfarr_push(struct zfarr *ca, int i, const curmat &A, cudaStream_t stream){
    X(mat) *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    zfarr_push(ca, i, tmp);
    X(free)(tmp);
}

void zfarr_push(struct zfarr *ca, int i, const cucmat &A, cudaStream_t stream){
    XC(mat) *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    zfarr_push(ca, i, tmp);
    XC(free)(tmp);
}

void zfarr_push(struct zfarr *ca, int i, const curcell &A, cudaStream_t stream){
    X(cell) *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    zfarr_push(ca, i, tmp);
    X(cellfree)(tmp);
}

void zfarr_push(struct zfarr *ca, int i, const cuccell &A, cudaStream_t stream){
    XC(cell) *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    zfarr_push(ca, i, tmp);
    XC(cellfree)(tmp);
}

void drawopdamp_gpu(const char *fig, loc_t *loc, const curmat &opd, cudaStream_t stream, 
		    const real *amp, real *zlim,
		    const char *title, const char *xlabel, const char *ylabel,
		    const char* format,...){
    format2fn;
    if(draw_current(fig, fn)){
	dmat *tmp=NULL;
	cp2cpu(&tmp, opd, stream); 
	drawopdamp(fig, loc, tmp->p, amp, zlim, title, xlabel, ylabel, "%s", fn);
	dfree(tmp);
    }
}
void drawpsf_gpu(const char *fig, curmat &psf, int count, cudaStream_t stream, int plotpsf,
		  const char *title, const char *xlabel, const char *ylabel,
		  const char* format,...){
    format2fn;
    if(draw_current(fig, fn)){
	dmat *psftemp=NULL;
	cp2cpu(&psftemp, psf, stream);
	if(count!=1) dscale(psftemp, 1./count);
	if(plotpsf==2){
	    dcwlog10(psftemp);
	}
	ddraw(fig, psftemp, NULL, NULL, title, xlabel, ylabel, "%s", fn);
	dfree(psftemp);
    }
}
/**
   Free data if not referenced or reference is 1.
*/
#undef cudaFree
int mycudaFree(void *pp){
    if(!pp) return 0;
    int tofree=1;
    lock_t tmp(cuglobal->memmutex);
    std::map<void*, int>::iterator it=cuglobal->memcount.find(pp);
    if(it!=cuglobal->memcount.end()){
	it->second--;
	tofree=!(it->second);
    }
    if(tofree){
	return cudaFree(pp);
    }else{
	return 0;
    }
}
#undef cudaMalloc
int mycudaMalloc(void **p, size_t size){
    int ans=cudaMalloc(p, size);
    return ans;
}
