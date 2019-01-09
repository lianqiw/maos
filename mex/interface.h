/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_MEX_INTERFACE_H
#define AOS_MEX_INTERFACE_H

#include <stdint.h>
#include <setjmp.h>
#if __GNUC__ && defined(__STDC_UTF_16__) && !defined(__cplusplus)
typedef int16_t char16_t;
#endif
#include <mex.h>
#include "../lib/aos.h"
#include "aolibmex.h"

#define REFERENCE 0
/*
  2014-07-08:
  We reference data instead of copying.
 */
extern int donotquit;
mxArray *loc2mx(const loc_t*loc){
    if(!loc) return mxCreateDoubleMatrix(0,0,mxREAL);
    mxArray *out=mxCreateDoubleMatrix(loc->nloc, 2, mxREAL);
    memcpy(mxGetPr(out), loc->locx, loc->nloc*sizeof(double));
    memcpy(mxGetPr(out)+loc->nloc, loc->locy, loc->nloc*sizeof(double));
    return out;
}
mxArray *dsp2mx(const dsp*A){
    if(!A) return mxCreateSparse(0, 0, 0, mxREAL);
    mxArray *out=0;
    if(REFERENCE){
	out=mxCreateSparse(0, 0, 0, mxREAL);
	mxSetPr(out, A->x);
	mxSetJc(out, (mwIndex*)A->p);
	mxSetIr(out, (mwIndex*)A->i);
	mxSetM(out, A->nx);
	mxSetN(out, A->ny);
	mxSetNzmax(out, A->nzmax);
	if(A->nref) A->nref[0]++;
    }else{
	out=mxCreateSparse(A->nx,A->ny,A->nzmax,mxREAL);
	memcpy(mxGetIr(out),A->i,A->nzmax*sizeof(long));
	memcpy(mxGetJc(out),A->p,(A->ny+1)*sizeof(long));
	memcpy(mxGetPr(out),A->x,A->nzmax*sizeof(double));
    }
    return out;
}
mxArray *csp2mx(const csp*A){
    if(!A) return mxCreateSparse(0, 0, 0, mxCOMPLEX);
    mxArray *out=0;
    out=mxCreateSparse(A->nx,A->ny,A->nzmax,mxCOMPLEX);
    memcpy(mxGetIr(out),A->i,A->nzmax*sizeof(long));
    memcpy(mxGetJc(out),A->p,(A->ny+1)*sizeof(long));
#if MX_HAS_INTERLEAVED_COMPLEX
    memcpy(mxGetData(out),A->x,A->nzmax*sizeof(dcomplex));
#else
    double *pr=mxGetPr(out);
    double *pi=mxGetPi(out);
    for(long i=0; i<A->nzmax; i++){
	pr[i]=creal(A->x[i]);
	pi[i]=cimag(A->x[i]);
    }
#endif
    return out;
}
mxArray *d2mx(const dmat *A){
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
mxArray *c2mx(const cmat *A){
    if(!A) return mxCreateDoubleMatrix(0,0,mxCOMPLEX);
    mxArray *out=0;
    out=mxCreateDoubleMatrix(A->nx, A->ny, mxCOMPLEX);
#if MX_HAS_INTERLEAVED_COMPLEX
    memcpy(mxGetData(out),A->p,A->nx*A->ny*sizeof(dcomplex));
#else
    double *pr=mxGetPr(out);
    double *pi=mxGetPi(out);
    long i;
    for(i=0; i<A->nx*A->ny; i++){
	pr[i]=creal(A->p[i]);
	pi[i]=cimag(A->p[i]);
    }
#endif
    return out;
}
lmat *d2l(const dmat *A){
    lmat *out=lnew(A->nx, A->ny);
    for(long i=0; i<A->nx*A->ny; i++){
	out->p[i]=(long)(A->p[i]);
    }
    return out;
}
dmat *l2d(const lmat *A){
    dmat *out=dnew(A->nx, A->ny);
    for(long i=0; i<A->nx*A->ny; i++){
	out->p[i]=(double)(A->p[i]);
    }
    return out;
}

mxArray *l2mx(const lmat *A){
    dmat *out=l2d(A);
    mxArray *B=d2mx(out);
    dfree(out);
    return B;
}
/*
mxArray *dcell2mx(const dcell *A){
    if(!A) return mxCreateCellMatrix(0,0);
    mxArray *out=mxCreateCellMatrix(A->nx,A->ny);
    for(int i=0; i<A->nx*A->ny; i++){
	if(A->p[i]){
	    mxSetCell(out, i, d2mx(A->p[i]));
	}
    }
    return out;
}
mxArray *ccell2mx(const ccell *A){
    if(!A) return mxCreateCellMatrix(0,0);
    mxArray *out=mxCreateCellMatrix(A->nx,A->ny);
    for(int i=0; i<A->nx*A->ny; i++){
	if(A->p[i]){
	    mxSetCell(out, i, c2mx(A->p[i]));
	}
    }
    return out;
    }*/
mxArray *any2mx(const void *A_){
    mxArray *out=0;
    const cell *A=(const cell*)A_;
    long id=A?(A->id):0;
    switch(id){
    case 0:
	break;
    case MCC_ANY:
	out=mxCreateCellMatrix(A->nx, A->ny);
	for(int i=0; i<A->nx*A->ny; i++){
	    if(A->p[i]) mxSetCell(out, i, any2mx(A->p[i]));
	}
	break;
    case M_DBL:
	out=d2mx((dmat*)A_);
	break;
    case M_CMP:
	out=c2mx((cmat*)A_);
	break;
    case M_LOC64:
	out=loc2mx((loc_t*)A_);
	break;
    case M_DSP64:
	out=dsp2mx((dsp*)A_);
	break;
    case M_INT64:
	out=l2mx((lmat*)A_);
	break;
    case M_INT32:
	out=l2mx((lmat*)A_);
	break;
    default:
	dbg("id=%ld is not handled.\n", id);
	out=mxCreateCellMatrix(0,0);
    }
    return out;
}
mxArray *str2mx(const char *str){
    return mxCreateString(str);
}
dsp *mx2dsp(const mxArray *A){
    if(!mxIsDouble(A) || mxIsComplex(A)) error("Only double is supported\n");
    dsp *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=(dsp*)calloc(1, sizeof(dsp));
	out->id=M_DSP64;
	out->nz=-1;
	out->nx=mxGetM(A);
	out->ny=mxGetN(A);
	out->p=(spint*)mxGetJc(A);
	out->i=(spint*)mxGetIr(A);
	out->x=mxGetPr(A);
	out->nzmax=mxGetNzmax(A);
    }
    return out;
}
#if MX_HAS_INTERLEAVED_COMPLEX
csp *mx2csp(const mxArray *A){
    if(!mxIsDouble(A) || !mxIsComplex(A)) error("Only double dcomplex is supported\n");
    csp *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=(csp*)calloc(1, sizeof(csp));
	out->id=M_CSP64;
	out->nz=-1;
	out->nx=mxGetM(A);
	out->ny=mxGetN(A);
	out->p=(spint*)mxGetJc(A);
	out->i=(spint*)mxGetIr(A);
	out->x=mxGetData(A);
	out->nzmax=mxGetNzmax(A);
    }
    return out;
}
#endif
loc_t *mx2loc(const mxArray *A){
    if(!mxIsDouble(A)) error("Only double is supported\n");
    loc_t *loc=(loc_t*)calloc(1, sizeof(loc_t));
    loc->nref=malloc(sizeof(int)); loc->nref[0]=1;
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

dmat *mx2d(const mxArray *A){
    if(!mxIsDouble(A)) error("Only double is supported\n");
    if(mxIsComplex(A)){
	mexErrMsgTxt("A is complex");
    }
    if(mxGetIr(A)){
	mexErrMsgTxt("A is sparse");
    }
    dmat *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=dnew_ref(mxGetM(A), mxGetN(A), mxGetPr(A));
    }
    return out;
}
cmat *mx2c(const mxArray *A){
#if MX_HAS_INTERLEAVED_COMPLEX
    if(!mxIsComplex(A)){
	mexErrMsgTxt("A is not complex");
    }
    if(mxGetIr(A)){
	mexErrMsgTxt("A is sparse");
    }
    cmat *out=0;
    if(A && mxGetM(A) && mxGetN(A)){
	out=cnew_ref(mxGetM(A), mxGetN(A), mxGetData(A));
    }
    return out;
#else
    mexErrMsgTxt("mx2c is not yet implemented\n");
    (void)A;
#endif
    return 0;
}
lmat *mx2l(const mxArray *A){
    dmat *B=mx2d(A);
    lmat *out=d2l(B);
    dfree(B);
    return out;
}
/*
dmat *mx2dvec(const mxArray *A){
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
dcell *mx2dcell(const mxArray *A){
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
    }*/
static void *mx2any(const mxArray *A, void*(*fun)(const mxArray*)){
    if(!A) return NULL;
    else if(!mxIsCell(A)){
	if(fun){
	    return fun(A);
	} else if(mxIsComplex(A)){
#if MX_HAS_INTERLEAVED_COMPLEX
	    if(mxIsSparse(A)){
		return mx2csp(A);
	    }else{
		return mx2c(A);
	    }
#else
	    error("Complex type not handled by mx2any\n");
	    return NULL;
#endif
	}else{
	    if(mxIsSparse(A)){
		return mx2dsp(A);
	    }else{
		return mx2d(A);
	    }
	}
    }else{
	cell *out=0;
	if(A && mxGetM(A) && mxGetN(A)){
	    out=cellnew(mxGetM(A), mxGetN(A));
	    for(int i=0; i<out->nx*out->ny; i++){
		mxArray *Ai=mxGetCell(A, i);
		out->p[i]=(cell*)mx2any(Ai, fun);
	    }
	}
	return out;
    }
}
static dcell *mx2dcell(const mxArray *A){
    return (dcell*)mx2any(A, (void*(*)(const mxArray*))mx2d);
}
static dcell *mx2ccell(const mxArray *A){
    return (dcell*)mx2any(A, (void*(*)(const mxArray*))mx2c);
}
static lcell *mx2lcell(const mxArray *A){
    return (lcell*)mx2any(A, (void*(*)(const mxArray*))mx2l);
}
static loccell *mx2loccell(const mxArray *A){
    return (loccell*)mx2any(A, (void*(*)(const mxArray*))mx2loc);
}
static dcell *mx2dspcell(const mxArray *A){
    return (dcell*)mx2any(A, (void*(*)(const mxArray*))mx2dsp);
}
static cell *mx2cell(const mxArray*A){
    return (cell*) mx2any(A, NULL);
}

static kalman_t *mx2kalman(const mxArray*A){
    kalman_t *kalman=(kalman_t*)calloc(1, sizeof(kalman_t));
    kalman->Ad=mx2d(mxGetField(A,0,"Ad"));
    kalman->Cd=mx2dcell(mxGetField(A,0,"Cd"));
    kalman->AdM=mx2d(mxGetField(A,0,"AdM"));
    kalman->FdM=mx2d(mxGetField(A,0,"FdM"));
    kalman->M=mx2dcell(mxGetField(A,0,"M"));
    kalman->P=mx2d(mxGetField(A,0,"P"));
    kalman->dthi=(double)mxGetScalar(mxGetField(A,0,"dthi"));
    kalman->dtrat=mx2d(mxGetField(A,0,"dtrat"));
    kalman->Gwfs=mx2dcell(mxGetField(A,0,"Gwfs"));
    kalman->Rwfs=mx2dcell(mxGetField(A,0,"Rwfs"));
    return kalman;
}
static mxArray* kalman2mx(kalman_t *kalman){
    const int nfield=12;
    const char *fieldnames[]={"Ad","Cd","AdM","FdM","Qn","Rn","M","P", "dthi", "dtrat", "Gwfs", "Rwfs"};
    mxArray *A=mxCreateStructMatrix(1,1,nfield,fieldnames);
    int pos=0;
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Ad));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Cd));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->AdM));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->FdM));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Qn));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Rn));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->M));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->P));
    mxSetFieldByNumber(A, 0, pos++, mxCreateDoubleScalar(kalman->dthi));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->dtrat));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Gwfs));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Rwfs));
    if(pos!=nfield){
	error("Invalid number of elements\n");
    }
    return A;
}
static mxArray *cn2est2mx(cn2est_t *cn2est){
    const int nfield=12;
    const char *fieldnames[]={"htrecon","wtrecon","r0m","ht","wt","r0","Pnk","iPnk","wtconvert","overlapi","cov2","cov1"};
    int pos=0;
    mxArray *A=mxCreateStructMatrix(1,1,nfield,fieldnames);
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->htrecon));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->wtrecon->p[0]));
    mxSetFieldByNumber(A, 0, pos++, mxCreateDoubleScalar(cn2est->r0m));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->ht));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->wt));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->r0));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->Pnk));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->iPnk));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->wtconvert));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->overlapi));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->cov2));
    mxSetFieldByNumber(A, 0, pos++, any2mx(cn2est->cov1));
    if(pos!=nfield){
	error("Invalid number of elements\n");
    }
    return A;

}
char *mx2str(const mxArray *A){
    int nlen=mxGetNumberOfElements(A)+1;
    char *fn=(char*)malloc(nlen);
    mxGetString(A, fn, nlen);
    return fn;
}
rand_t *mx2rand(const mxArray *A){
    int seed=(int)mxGetScalar(A);
    rand_t *out=(rand_t*)malloc(sizeof(rand_t));
    seed_rand(out, seed);
    return out;
}
static void mex_signal_handler(int sig){
    if(sig){
	mexErrMsgTxt("Signal caught.\n");
    }else{
	dbg("signal 0 caught\n");
    }
}
static jmp_buf *exception_env;
static void mex_quitfun(const char *msg){
    if(exception_env){
	info("error: %s\n", msg);
	longjmp(*exception_env, 1);
    }else{
	info("mexerrmsg\n");
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
	extern int mem_debug;
	mem_debug=1;
	extern void *(*calloc_custom)(size_t, size_t);
	extern void *(*malloc_custom)(size_t);
	extern void *(*realloc_custom)(void *, size_t);
	extern void  (*free_custom)(void *);
	calloc_custom=calloc_mex;
	malloc_custom=malloc_mex;
	realloc_custom=realloc_mex;
	free_custom=free_mex;
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
