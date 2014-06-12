#ifndef AOS_MEX_INTERFACE_H
#define AOS_MEX_INTERFACE_H
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <stdint.h>
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
extern int donotquit;
INLINE mxArray *dsp2mx(const dsp*A){
    if(!A) return 0;
    mxArray *out=mxCreateSparse(A->m,A->n,A->nzmax,mxREAL);
    memcpy(mxGetIr(out),A->i,A->nzmax*sizeof(long));
    memcpy(mxGetJc(out),A->p,(A->n+1)*sizeof(long));
    memcpy(mxGetPr(out),A->x,A->nzmax*sizeof(double));
    return out;
}
INLINE dsp *mx2dsp(const mxArray *A){
    dsp *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=calloc(1, sizeof(dsp));
	out->m=mxGetM(A);
	out->n=mxGetN(A);
	out->p=mxGetJc(A);
	out->i=mxGetIr(A);
	out->x=mxGetPr(A);
	out->nzmax=mxGetNzmax(A);
	if(mxGetPi(A)){
	    mexErrMsgTxt("A is complex");
	}
    }
    return out;
}
INLINE mxArray *d2mx(const dmat *A){
    mxArray *out=0;
    if(A && A->nx && A->ny){
	out=mxCreateDoubleMatrix(A->nx,A->ny,mxREAL);
	memcpy(mxGetPr(out),A->p,A->nx*A->ny*sizeof(double));
    }
    return out;
}
INLINE mxArray *dcell2mx(const dcell *A){
    mxArray *out=mxCreateCellMatrix(A->nx,A->ny);
    for(int i=0; i<A->nx*A->ny; i++){
	if(A->p[i]){
	    mxSetCell(out, i, d2mx(A->p[i]));
	}
    }
    return out;
}
INLINE mxArray *c2mx(const cmat *A){
    mxArray *out=mxCreateDoubleMatrix(A->nx, A->ny, mxCOMPLEX);
    double *pr=mxGetPr(out);
    double *pi=mxGetPi(out);
    long i;
    for(i=0; i<A->nx*A->ny; i++){
	pr[i]=creal(A->p[i]);
	pi[i]=cimag(A->p[i]);
    }
    return out;
}

INLINE loc_t *mx2loc(const mxArray *A){
    loc_t *loc=calloc(1, sizeof(loc_t));
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
    if(mxGetPi(A)){
	mexErrMsgTxt("A is complex");
    }
    if(mxGetIr(A)){
	mexErrMsgTxt("A is dsp");
    }
    dmat *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=dnew_ref( mxGetM(A), mxGetN(A), mxGetPr(A));
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
	out=dcellnew(mxGetM(A), mxGetN(A));
	for(int i=0; i<out->nx*out->ny; i++){
	    mxArray *Ai=mxGetCell(A, i);
	    if(Ai && mxGetM(Ai) && mxGetN(Ai)){
		out->p[i]=dnew_ref( mxGetM(Ai), mxGetN(Ai), mxGetPr(Ai));
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
/*
static void mex_signal_handler(int sig){
    if(sig){
	mexErrMsgTxt("Signal caught.\n");
    }else{
	info("signal 0 caught\n");
    }
}
static void(*default_handler)(int)=NULL;
static __attribute__((constructor)) void init(){
    if(!default_handler){
	default_handler=signal(SIGTERM, mex_signal_handler);
    }
    quitfun=mexErrMsgTxt;
}
static __attribute__((destructor)) void deinit(){
    if(default_handler){
	signal(SIGTERM, default_handler);
    }else{
	signal(SIGTERM, SIG_DFL);
    }
    }*/
#endif
