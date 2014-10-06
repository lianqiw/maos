#ifndef AOS_MEX_INTERFACE_H
#define AOS_MEX_INTERFACE_H
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <stdint.h>
#include <setjmp.h>
typedef uint16_t char16_t;
#include <mex.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "../lib/aos.h"
#ifdef __cplusplus
}
#endif
#ifndef INLINE
#define INLINE inline __attribute__((always_inline))
#endif
#define REFERENCE 0
/*
  2014-07-08:
  We reference data instead of copying.
 */
extern int donotquit;
INLINE mxArray *loc2mx(const loc_t*loc){
    if(!loc) return mxCreateDoubleMatrix(0,0,mxREAL);
    mxArray *out=mxCreateDoubleMatrix(loc->nloc, 2, mxREAL);
    memcpy(mxGetPr(out), loc->locx, loc->nloc*sizeof(double));
    memcpy(mxGetPr(out)+loc->nloc, loc->locy, loc->nloc*sizeof(double));
    return out;
}
INLINE mxArray *dsp2mx(const dsp*A){
    if(!A) return mxCreateSparse(0, 0, 0, mxREAL);
    mxArray *out=0;
    if(REFERENCE){
	out=mxCreateSparse(0, 0, 0, mxREAL);
	mxSetPr(out, A->x);
	mxSetJc(out, A->p);
	mxSetIr(out, A->i);
	mxSetM(out, A->nx);
	mxSetN(out, A->ny);
	mxSetNzmax(out, A->nzmax);
	if(A->nref) A->nref[0]++;
    }else{
	out=mxCreateSparse(A->m,A->n,A->nzmax,mxREAL);
	memcpy(mxGetIr(out),A->i,A->nzmax*sizeof(long));
	memcpy(mxGetJc(out),A->p,(A->n+1)*sizeof(long));
	memcpy(mxGetPr(out),A->x,A->nzmax*sizeof(double));
    }
    return out;
}
INLINE mxArray *csp2mx(const csp*A){
    if(!A) return mxCreateSparse(0, 0, 0, mxCOMPLEX);
    mxArray *out=0;
    out=mxCreateSparse(A->m,A->n,A->nzmax,mxCOMPLEX);
    memcpy(mxGetIr(out),A->i,A->nzmax*sizeof(long));
    memcpy(mxGetJc(out),A->p,(A->n+1)*sizeof(long));
    double *pr=mxGetPr(out);
    double *pi=mxGetPi(out);
    for(long i=0; i<A->nzmax; i++){
	pr[i]=creal(A->x[i]);
	pi[i]=cimag(A->x[i]);
    }
    return out;
}
INLINE mxArray *d2mx(const dmat *A){
    if(!A) return mxCreateDoubleMatrix(0,0,mxREAL);
    mxArray *out=0;
    if(REFERENCE && !A->mmap && A->nref){
	out=mxCreateDoubleMatrix(0,0,mxREAL);
	mxSetPr(out, A->p);
	mxSetM(out, A->nx);
	mxSetN(out, A->ny);
	if(A->nref) A->nref[0]++;
    }else{
	out=mxCreateDoubleMatrix(A->nx,A->ny,mxREAL);
	memcpy(mxGetPr(out),A->p,A->nx*A->ny*sizeof(double));
    }
    return out;
}
INLINE mxArray *c2mx(const cmat *A){
    if(!A) return mxCreateDoubleMatrix(0,0,mxCOMPLEX);
    mxArray *out=0;
    out=mxCreateDoubleMatrix(A->nx, A->ny, mxCOMPLEX);
    double *pr=mxGetPr(out);
    double *pi=mxGetPi(out);
    long i;
    for(i=0; i<A->nx*A->ny; i++){
	pr[i]=creal(A->p[i]);
	pi[i]=cimag(A->p[i]);
    }
    return out;
}
INLINE mxArray *lmat2mx(const lmat *A){
    if(!A) return mxCreateDoubleMatrix(0,0,mxREAL);
    mxArray *out;
    if(sizeof(long)==8){
	out=mxCreateNumericMatrix(A->nx, A->ny,mxINT64_CLASS,mxREAL);
    }else{
	out=mxCreateNumericMatrix(A->nx, A->ny,mxINT32_CLASS,mxREAL);
    }
    memcpy(mxGetPr(out),A->p,A->nx*A->ny*sizeof(long));
    return out;
}
INLINE mxArray *dcell2mx(const dcell *A){
    if(!A) return mxCreateCellMatrix(0,0);
    mxArray *out=mxCreateCellMatrix(A->nx,A->ny);
    for(int i=0; i<A->nx*A->ny; i++){
	if(A->p[i]){
	    mxSetCell(out, i, d2mx(A->p[i]));
	}
    }
    return out;
}
INLINE mxArray *ccell2mx(const ccell *A){
    if(!A) return mxCreateCellMatrix(0,0);
    mxArray *out=mxCreateCellMatrix(A->nx,A->ny);
    for(int i=0; i<A->nx*A->ny; i++){
	if(A->p[i]){
	    mxSetCell(out, i, c2mx(A->p[i]));
	}
    }
    return out;
}
mxArray *any2mx(const void *A_){
    mxArray *out=0;
    const cell *A=A_;
    long id=A?(A->id):0;
    switch(id){
    case MCC_ANY:
	out=mxCreateCellMatrix(A->nx, A->ny);
	for(int i=0; i<A->nx*A->ny; i++){
	    if(A->p[i]) mxSetCell(out, i, any2mx(A->p[i]));
	}
	break;
    case M_DBL:
	out=d2mx(A_);
	break;
    case M_CMP:
	out=c2mx(A_);
	break;
    case M_LOC64:
	out=loc2mx(A_);
	break;
    case M_DSP64:
	out=dsp2mx(A_);
	break;
    case M_INT64:
	out=lmat2mx(A_);
	break;
    case M_INT32:
	out=lmat2mx(A_);
	break;
    default:
	out=mxCreateCellMatrix(0,0);
    }
    return out;
}
INLINE mxArray *str2mx(const char *str){
    return mxCreateString(str);
}
INLINE dsp *mx2dsp(const mxArray *A){
    if(!mxIsDouble(A) || mxIsComplex(A)) error("Only double is supported\n");
    dsp *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=calloc(1, sizeof(dsp));
	out->m=mxGetM(A);
	out->n=mxGetN(A);
	out->p=mxGetJc(A);
	out->i=mxGetIr(A);
	out->x=mxGetPr(A);
	out->nzmax=mxGetNzmax(A);
    }
    return out;
}
INLINE loc_t *mx2loc(const mxArray *A){
    if(!mxIsDouble(A)) error("Only double is supported\n");
    loc_t *loc=calloc(1, sizeof(loc_t));
    loc->ref=1;
    loc->locx=mxGetPr(A);
    loc->nloc=mxGetM(A);
    loc->locy=loc->locx+loc->nloc;
    const double tol=1e-7;
    double dxd=INFINITY, dyd=INFINITY;
    for(long i=0; i<loc->nloc-1; i++){
	double dxi=fabs(loc->locx[i+1]-loc->locx[i]);
	if(dxi>tol && dxi+tol<dxd){
	    dxd=dxi;
	}
	double dyi=fabs(loc->locy[i+1]-loc->locy[i]);
	if(dyi>tol && dyi+tol<dyd){
	    dyd=dyi;
	}
    }
    loc->dx=dxd;
    loc->dy=dyd;
    return loc;
}

INLINE dmat *mx2d(const mxArray *A){
    if(!mxIsDouble(A)) error("Only double is supported\n");
    if(mxGetPi(A)){
	mexErrMsgTxt("A is complex");
    }
    if(mxGetIr(A)){
	mexErrMsgTxt("A is dsp");
    }
    dmat *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=dnew_ref(mxGetM(A), mxGetN(A), mxGetPr(A));
    }
    return out;
}
INLINE dmat *mx2dvec(const mxArray *A){
    dmat *out=mx2d(A);
    if(out->nx==1){
	out->nx=out->ny;
	out->ny=1;
    }else if(out->ny>1){
	fprintf(stderr, "Size is %ldx%ld\n", out->nx, out->ny);
	mexErrMsgTxt("Input is not a vector");
    }
    return out;
}
INLINE dcell *mx2dcell(const mxArray *A){
    if(!mxIsCell(A)){
	mexErrMsgTxt("A is not cell");
    }
    dcell *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=cellnew(mxGetM(A), mxGetN(A));
	for(int i=0; i<out->nx*out->ny; i++){
	    mxArray *Ai=mxGetCell(A, i);
	    out->p[i]=mx2d(Ai);
	}
    }
    return out;
}
static void *mx2any(const mxArray *A){
    if(!mxIsCell(A)){
	mexErrMsgTxt("A is not cell");
    }
    cell *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=cellnew(mxGetM(A), mxGetN(A));
	for(int i=0; i<out->nx*out->ny; i++){
	    mxArray *Ai=mxGetCell(A, i);
	    if(mxIsCell(Ai)){
		out->p[i]=mx2any(Ai);
	    }else{
		out->p[i]=(void*)mx2d(Ai);
	    }
	}
    }
    return out;
}
INLINE char *mx2str(const mxArray *A){
    int nlen=mxGetNumberOfElements(A)+1;
    char *fn=malloc(nlen);
    mxGetString(A, fn, nlen);
    return fn;
}

static void mex_signal_handler(int sig){
    if(sig){
	mexErrMsgTxt("Signal caught.\n");
    }else{
	info("signal 0 caught\n");
    }
}
static jmp_buf *exception_env;
static void mex_quitfun(const char *msg){
    if(exception_env){
	info2("longjmp\n");
	longjmp(*exception_env, 1);
    }else{
	info2("mexerrmsg\n");
	mexErrMsgTxt(msg);
    }
}
static void(*default_handler)(int)=NULL;
static void *calloc_mex(size_t nmemb, size_t size){
    void *p=mxCalloc(nmemb, size);
    mexMakeMemoryPersistent(p);
    return p;
}
static void *malloc_mex(size_t size){
    void *p=mxMalloc(size);
    mexMakeMemoryPersistent(p);
    return p;
}
static void *realloc_mex(void *p, size_t size){
    p=mxRealloc(p, size);
    mexMakeMemoryPersistent(p);
    return p;
}
static void free_mex(void*p){
    mxFree(p);
}
static __attribute__((constructor)) void init(){
    if(!default_handler){
	default_handler=signal(SIGTERM, mex_signal_handler);
    }
    quitfun=mex_quitfun;
    if(REFERENCE){
	CALLOC=calloc_mex;
	MALLOC=malloc_mex;
	REALLOC=realloc_mex;
	FREE=free_mex;
    }
}
static __attribute__((destructor)) void deinit(){
    fprintf(stderr, "mex unloaded\n");
    if(default_handler){
	signal(SIGTERM, default_handler);
    }else{
	signal(SIGTERM, SIG_DFL);
    }
}
#endif
