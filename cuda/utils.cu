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
#include <errno.h>
#include <pthread.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <cufft.h>
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
   Copy map_t to cumap_t. if type==1, use cudaArray, otherwise use float
   array. Allow multiple calling to override the data.  */
void cp2gpu(cumap_t **dest0, map_t **source, int nps){
    if(nps==0) return;
    if(!*dest0){
	*dest0=new cumap_t[nps];
	for(int ips=0; ips<nps; ips++){
	    (*dest0)[ips].init(source[ips]);
	}	
    }
    cumap_t *dest=*dest0;
    for(int ips=0; ips<nps; ips++){
	cp2gpu(&dest[ips].p, (dmat*)source[ips]);
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

//  Need debugging. result not correct.
/*void cusp::trans(){
  if(!p) return;
  int ncol, nrow;
  switch(type){
  case SP_CSR:
  nrow=nx;
  ncol=ny;
  break;
  case SP_CSC:
  nrow=ny;
  ncol=nx;
  break;
  default:
  nrow=0;ncol=0;
  error("Invalid format\n");
  }
  float *xnew;
  int *inew, *pnew;
  cudaMalloc(&xnew, nzmax*sizeof(float));
  cudaMalloc(&inew, nzmax*sizeof(int));
  cudaMalloc(&pnew, (ncol+1)*sizeof(int));
  stream_t stream;
  cusparseScsr2csc(stream, nrow, ncol,
  x, p, i, xnew, inew, pnew,
  1, CUSPARSE_INDEX_BASE_ZERO);
  cudaFree(x); x=xnew;
  cudaFree(i); i=inew;
  cudaFree(p); p=pnew;
  }*/
/*
  y=A*x where A is sparse. x, y are vectors. Slow for GS0.
*/

void cuspmul(float *y, cusp *A, const float *x, int ncolvec, char trans, float alpha, cusparseHandle_t handle){
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
    if(ncolvec==1){
	status=cusparseScsrmv(handle, opr,
			      nrow, ncol, alpha, spdesc,
			      A->x, A->p, A->i, x, 1.f, y);
    }else{
	int nlead=istrans?nrow:ncol;
	status=cusparseScsrmm(handle, opr,
			      nrow, ncolvec, ncol, alpha, spdesc,
			      A->x, A->p, A->i, x, nlead, 1.f, y, nlead);
    }
    if(status!=0){
	error("cusparseScsrmv(m) failed with status '%s'\n", scsrmv_err[status]);
    }
}

/**
   Convert a source loc_t to device memory.
*/
void cp2gpu(float (* restrict *dest)[2], const loc_t *src){
    float (*tmp)[2]=(float(*)[2])malloc(src->nloc*2*sizeof(float));
    for(int iloc=0; iloc<src->nloc; iloc++){
	tmp[iloc][0]=(float)src->locx[iloc];
	tmp[iloc][1]=(float)src->locy[iloc];
    }
    cp2gpu((float**)dest, (float*)tmp, src->nloc*2, 1);
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
/**
   Convert device (float) array and add to host double.
   dest = alpha * dest + beta *src;
*/
void cp2cpu(double * restrict *dest, double alpha, float *src, double beta, int n, 
	    cudaStream_t stream, pthread_mutex_t *mutex){
    CUDA_SYNC_STREAM;
    float *tmp=(float*)malloc(n*sizeof(float));
    DO(cudaMemcpy(tmp, src, n*sizeof(float), cudaMemcpyDeviceToHost));
    if(!*dest){
	*dest=(double*)malloc(sizeof(double)*n);
    }
    double *restrict p=*dest;
    if(mutex) LOCK(*mutex);
    for(int i=0; i<n; i++){
	p[i]=p[i]*alpha+beta*tmp[i];
    }
    if(mutex) UNLOCK(*mutex);
    free(tmp);
}
/*
  Write float on gpu to file
*/
void gpu_write(const float *p, int nx, int ny, const char *format, ...){
    format2fn;
    float *tmp=(float*)malloc(nx*ny*sizeof(float));
    cudaMemcpy(tmp, p, nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
    writeflt(tmp,nx,ny,"%s",fn);
    free(tmp);
}

/*
  Write float on gpu to file
*/
void gpu_write(const fcomplex *p, int nx, int ny, const char *format, ...){
    format2fn;
    fcomplex *tmp=(fcomplex*)malloc(nx*ny*sizeof(fcomplex));
    cudaMemcpy(tmp, p, nx*ny*sizeof(fcomplex), cudaMemcpyDeviceToHost);
    writefcmp((float complex*)tmp,nx,ny,"%s",fn);
    free(tmp);
}
/*
  Write float on gpu to file
*/
void gpu_write(const int *p, int nx, int ny, const char *format, ...){
    format2fn;
    int *tmp=(int*)malloc(nx*ny*sizeof(int));
    cudaMemcpy(tmp, p, nx*ny*sizeof(int), cudaMemcpyDeviceToHost);
    writeint(tmp,nx,ny,"%s",fn);
    free(tmp);
}

void cp2cpu(dmat **out, double alpha, const curmat *in, double beta, cudaStream_t stream, pthread_mutex_t *mutex){
    if(!in){
	if(*out) dzero(*out);
	return;
    }
    if(!*out) {
	*out=dnew(in->nx, in->ny);
    }else{
	assert((*out)->nx*(*out)->ny==in->nx*in->ny);
    }
    cp2cpu(&(*out)->p, alpha, in->p, beta, in->nx*in->ny, stream, mutex);
}
void cp2cpu(dcell **out, double alpha, const curcell *in, double beta, cudaStream_t stream, pthread_mutex_t *mutex){
    if(!in){
	if(*out) dcellzero(*out);
	return;
    }
    if(!*out) *out=dcellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	cp2cpu(&(*out)->p[i], alpha, in->p[i], beta, stream, mutex);
    }
}
void cp2cpu(smat **out, const curmat *in, cudaStream_t stream){
    if(!in) {
	if(*out) szero(*out);
	return;
    }
    if(!*out) *out=snew(in->nx, in->ny);
    DO(cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(float), cudaMemcpyDeviceToHost, stream));
    if(in->header) (*out)->header=strdup(in->header);
}


void cp2cpu(zmat **out, const cucmat *in, cudaStream_t stream){
    if(!in){
	if(*out) zzero(*out);
	return;
    }
    if(!*out) *out=znew(in->nx, in->ny);
    DO(cudaMemcpyAsync((*out)->p, in->p, in->nx*in->ny*sizeof(fcomplex), cudaMemcpyDeviceToHost, stream));
    if(in->header) (*out)->header=strdup(in->header);
}

void cp2cpu(scell **out, const curcell *in, cudaStream_t stream){
    if(!in){
	if(*out) scellzero(*out);
	return;
    }
    if(!*out) *out=scellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	cp2cpu(&(*out)->p[i], in->p[i], stream);
    }
}

void cp2cpu(zcell **out, const cuccell *in, cudaStream_t stream){
    if(!in){
	if(*out) zcellzero(*out);
	return;
    }
    if(!*out) *out=zcellnew(in->nx, in->ny);
    for(int i=0; i<in->nx*in->ny; i++){
	cp2cpu(&(*out)->p[i], in->p[i], stream);
    }
}

void cellarr_cur(struct cellarr *ca, int i, const curmat *A, cudaStream_t stream){
    smat *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_smat(ca, i, tmp);
    sfree(tmp);
}

void cellarr_cuc(struct cellarr *ca, int i, const cucmat *A, cudaStream_t stream){
    zmat *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_zmat(ca, i, tmp);
    zfree(tmp);
}

void cellarr_curcell(struct cellarr *ca, int i, const curcell *A, cudaStream_t stream){
    scell *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_scell(ca, i, tmp);
    scellfree(tmp);
}

void cellarr_cuccell(struct cellarr *ca, int i, const cuccell *A, cudaStream_t stream){
    zcell *tmp=NULL;
    cp2cpu(&tmp,A,stream);
    CUDA_SYNC_STREAM;
    cellarr_zcell(ca, i, tmp);
    zcellfree(tmp);
}
