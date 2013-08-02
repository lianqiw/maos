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

static __attribute((constructor)) void init(){
    DO(cusparseCreateMatDescr(&spdesc));
    cusparseSetMatType(spdesc, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(spdesc, CUSPARSE_INDEX_BASE_ZERO);
}

/**
   Convert double array to device memory (float)
*/

void cp2gpu(float * restrict *dest, const double *src, int n){
    if(!src) return;
    float *tmp=(float*)malloc(n*sizeof(float));
    for(int i=0; i<n; i++){
	tmp[i]=(float)src[i];
    }
    if(!*dest){
	DO(cudaMalloc((float**)dest, n*sizeof(float)));
    }
    DO(cudaMemcpy(*dest, tmp, n*sizeof(float),cudaMemcpyHostToDevice));
    free(tmp);
}
/**
   Convert double array to device memory (float)
*/

void cp2gpu(fcomplex * restrict *dest, const dcomplex *restrict src, int n){
    if(!src) return;
    fcomplex *tmp=(fcomplex*)malloc(n*sizeof(fcomplex));
    for(int i=0; i<n; i++){
	tmp[i]=(make_cuFloatComplex)(cuCreal(src[i]), cuCimag(src[i]));
    }
    if(!*dest){
	DO(cudaMalloc((fcomplex**)dest, n*sizeof(fcomplex)));
    }
    DO(cudaMemcpy(*dest, tmp, n*sizeof(fcomplex),cudaMemcpyHostToDevice));
    free(tmp);
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
    CUDA_SYNC_DEVICE;
}

/*
  Convert a host dsp array to GPU sprase array.
*/
void cp2gpu(cusp **dest0, const dsp *src_csc, int tocsr){
    if(!src_csc) return;
    if(!*dest0) *dest0=(cusp*)calloc(1, sizeof(cusp));
    cusp *dest=*dest0;
    dsp *src=const_cast<dsp*>(src_csc);
    if(tocsr){
	dest->type=SP_CSR;
	src=sptrans(src_csc);
    }else{
	dest->type=SP_CSC;
    }
    dest->nx=src_csc->m;
    dest->ny=src_csc->n;
    dest->nzmax=src->nzmax;
    dest->p=NULL; dest->i=NULL; dest->x=NULL;
    cp2gpu(&dest->p, src->p, src->n+1);
    cp2gpu(&dest->i, src->i, src->nzmax);
    cp2gpu(&dest->x, src->x, src->nzmax);
    if(tocsr){
	spfree(src);
    }
}
void cp2gpu(cusp **dest0, const spcell *srcc, int tocsr){
    dsp *src=spcell2sp(srcc);
    cp2gpu(dest0, src, tocsr);
    spfree(src);
}
void cp2gpu(cuspcell **dest0, const spcell *src, int tocsr){
    if(!src) return;
    if(!*dest0){
	*dest0=cuspcellnew(src->nx, src->ny);
    }
    for(int i=0; i<src->nx*src->ny; i++){
	cp2gpu(&(*dest0)->p[i], src->p[i], tocsr);
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
    int ncol, nrow;
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
    if(!*dest){
	DO(cudaMalloc((float**)dest, src->nloc*2*sizeof(float)));
    }
    DO(cudaMemcpy(*dest, tmp, src->nloc*2*sizeof(float),cudaMemcpyHostToDevice));
    free(tmp);
}

/**
   Convert dmat array to device memory.
*/
void cp2gpu(float * restrict *dest, const dmat *src){
    if(!src) return;
    cp2gpu(dest, src->p, src->nx*src->ny);
}
/**
   Convert dmat array to curmat
*/
void cp2gpu(curmat *restrict *dest, const dmat *src){
    if(!src){
	curzero(*dest);
	return;
    }
    if(!*dest){
	*dest=curnew(src->nx, src->ny);
    }else{
	assert(src->nx*src->ny==(*dest)->nx*(*dest)->ny);
    }
    cp2gpu(&(*dest)->p, src->p, src->nx*src->ny);
}
void cp2gpu(curmat *restrict *dest, const float *src, int nx, int ny, cudaStream_t stream){
    if(!src){
	curzero(*dest);
	return;
    }
    if(!*dest){
	*dest=curnew(nx, ny);
    }else{
	assert(nx*ny==(*dest)->nx*(*dest)->ny);
    }
    if(stream){
	DO(cudaMemcpyAsync((*dest)->p, src, nx*ny*sizeof(float),cudaMemcpyHostToDevice, stream));
    }else{
	DO(cudaMemcpy((*dest)->p, src, nx*ny*sizeof(float),cudaMemcpyHostToDevice));
    }
}
/*
  convert cmat to cucmat
*/
void cp2gpu(cucmat *restrict *dest, const cmat *src){
    if(!src){
	czero(*dest);
	return;
    }
    if(!*dest){
	*dest=cucnew(src->nx, src->ny);
    }else{
	assert(src->nx*src->ny==(*dest)->nx*(*dest)->ny);
    }
    cp2gpu(&(*dest)->p, (dcomplex*)src->p, (int)(src->nx*src->ny));
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
   Convert dmat array to device memory.
*/
void cp2gpu(fcomplex * restrict *dest, const cmat *src){
    if(src){
	cp2gpu(dest, (dcomplex*)src->p, src->nx*src->ny);
    }
}
/**
   Convert double array to device memory (float)
*/
void dbl2flt(float * restrict *dest, const double *src, int n){
    if(!src) return;
    if(!*dest){
	cudaMallocHost((float**)dest, n*sizeof(float));
    }
    for(int i=0; i<n; i++){
	(*dest)[i]=(float)src[i];
    }
}
/**
   Convert long array to device int
*/
void cp2gpu(int * restrict *dest, const long *src, int n){
    if(!src) return;
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    if(sizeof(long)==sizeof(int)){
	DO(cudaMemcpy(*dest, src, n*sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
    }else{
	int *tmp=(int*)malloc(sizeof(int)*n);
	for(int i=0; i<n; i++){
	    tmp[i]=(int)src[i];
	    if((long)tmp[i]!=src[i]){
		error("Overflow occured\n");
	    }
	}
	DO(cudaMemcpy(*dest, tmp, n*sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	free(tmp);
    }
}
void cp2gpu(int *restrict *dest, const int *src, int n){
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    DO(cudaMemcpy(*dest, src, sizeof(int)*n, cudaMemcpyHostToDevice));
}
/**
   Convert long array to device int
*/
void cp2gpu(int * restrict *dest, const spint *src, int n){
    if(!*dest){
	DO(cudaMalloc((int**)dest, n*sizeof(int)));
    }
    if(sizeof(spint)==sizeof(int)){
	DO(cudaMemcpy(*dest, src, n*sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
    }else{
	int *tmp=(int*)malloc(sizeof(int)*n);
	for(int i=0; i<n; i++){
	    tmp[i]=(int)src[i];
	    if((spint)tmp[i]!=src[i]){
		error("Overflow occured\n");
	    }
	}
	DO(cudaMemcpy(*dest, tmp, n*sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	free(tmp);
    }
}
/**
   Convert device (float) array and add to host double.
   dest = alpha * dest + beta *src;
*/
void cp2cpu(double * restrict *dest, double alpha, float *src, double beta, int n, 
	    cudaStream_t stream, pthread_mutex_t *mutex){
    float *tmp=(float*)malloc4async(n*sizeof(float));
    DO(cudaMemcpyAsync(tmp, src, n*sizeof(float), cudaMemcpyDeviceToHost, stream));
    if(!*dest){
	*dest=(double*)malloc(sizeof(double)*n);
    }
    double *restrict p=*dest;
    CUDA_SYNC_STREAM;
    if(mutex) LOCK(*mutex);
    for(int i=0; i<n; i++){
	p[i]=p[i]*alpha+beta*tmp[i];
    }
    if(mutex) UNLOCK(*mutex);
    free4async(tmp);
}
/*
  Write float on gpu to file
*/
void gpu_write(const float *p, int nx, int ny, const char *format, ...){
    format2fn;
    float *tmp=(float*)malloc(nx*ny*sizeof(float));
    cudaDeviceSynchronize();
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
    cudaDeviceSynchronize();
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
    cudaDeviceSynchronize();
    cudaMemcpy(tmp, p, nx*ny*sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
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
W01_T* gpu_get_W01(dsp *R_W0, dmat *R_W1){
    if(!R_W0 || !R_W1){
	error("R0, R1 must not be empty\n");
    }
    W01_T *W01=(W01_T*)calloc(1, sizeof(W01_T));
    cp2gpu(&W01->W1, R_W1);
    {
	/*W0 of partially illuminates subaps are stored as sparse matrix in
	  GPU. W0 of fully illuminated subaps are not.*/
	spint *pp=R_W0->p;
	spint *pi=R_W0->i;
	double *px=R_W0->x;
	dsp *W0new=spnew(R_W0->m, R_W0->n, R_W0->nzmax);
	spint *pp2=W0new->p;
	spint *pi2=W0new->i;
	double *px2=W0new->x;
	int *full;
	cudaMallocHost(&full, R_W0->n*sizeof(int));
	//#define W0_BW 1
	double W1max=dmax(R_W1);
	double thres=W1max*(1.f-1e-6);
	W01->W0v=(float)(W1max*4./9.);//max of W0 is 4/9 of max of W1. 
	info("W0v=%g\n", W01->W0v);
	int count=0;
	int count2=0;
	for(int ic=0; ic<R_W0->n; ic++){
	    pp2[ic]=count;
	    if(R_W1->p[ic]>thres){
		full[count2]=ic;
		count2++;
	    }else{
		int nv=pp[ic+1]-pp[ic];
		memcpy(pi2+count, pi+pp[ic], sizeof(spint)*nv);
		memcpy(px2+count, px+pp[ic], sizeof(double)*nv);
		count+=nv;
	    }
	}
	pp2[R_W0->n]=count;
	W0new->nzmax=count;
	cp2gpu(&W01->W0p, W0new, 1);
	cp2gpu(&W01->W0f, full, count2);
	W01->nW0f=count2;
	spfree(W0new);
	cudaFreeHost(full);
    }
    return W01;
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
