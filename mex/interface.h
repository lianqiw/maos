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
#ifndef AOS_MEX_INTERFACE_H
#define AOS_MEX_INTERFACE_H
#include <signal.h>
#include <stdint.h>
#include <setjmp.h>
#if __GNUC__ && defined(__STDC_UTF_16__) && !defined(__cplusplus)
typedef int16_t char16_t;
#endif
#include <mex.h>
#include "../lib/aos.h"

/*
  2020-03-24: removed data reference from C to matlab capability. 
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
    out=mxCreateSparse(A->nx,A->ny,A->nzmax,mxREAL);
    memcpy(mxGetIr(out),A->pi,A->nzmax*sizeof(long));
    memcpy(mxGetJc(out),A->pp,(A->ny+1)*sizeof(long));
    memcpy(mxGetPr(out),A->px,A->nzmax*sizeof(double));
    return out;
}
mxArray *csp2mx(const csp*A){
    if(!A) return mxCreateSparse(0, 0, 0, mxCOMPLEX);
    mxArray *out=0;
    out=mxCreateSparse(A->nx,A->ny,A->nzmax,mxCOMPLEX);
    memcpy(mxGetIr(out),A->pi,A->nzmax*sizeof(long));
    memcpy(mxGetJc(out),A->pp,(A->ny+1)*sizeof(long));
#if MX_HAS_INTERLEAVED_COMPLEX
    memcpy(mxGetData(out),A->px,A->nzmax*sizeof(dcomplex));
#else
    double *pr=mxGetPr(out);
    double *pi=mxGetPi(out);
    for(long i=0; i<A->nzmax; i++){
	pr[i]=creal(A->px[i]);
	pi[i]=cimag(A->px[i]);
    }
#endif
    return out;
}
mxArray *d2mx(const dmat *A){
    if(!A) return mxCreateDoubleMatrix(0,0,mxREAL);
    mxArray *out=0;
    out=mxCreateDoubleMatrix(A->nx,A->ny,mxREAL);
    memcpy(mxGetPr(out),A->p,A->nx*A->ny*sizeof(double));
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
	out->nx=mxGetM(A);
	out->ny=mxGetN(A);
	out->pp=(spint*)mxGetJc(A);
	out->pi=(spint*)mxGetIr(A);
	out->px=mxGetPr(A);
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
	out->nx=mxGetM(A);
	out->ny=mxGetN(A);
	out->pp=(spint*)mxGetJc(A);
	out->pi=(spint*)mxGetIr(A);
	out->px=(comp*)mxGetData(A);
	out->nzmax=mxGetNzmax(A);
    }
    return out;
}
#endif
loc_t *mx2loc(const mxArray *A){
    if(!mxIsDouble(A)) error("Only double is supported\n");
    loc_t *loc=(loc_t*)calloc(1, sizeof(loc_t));
    //loc->nref=(int*)malloc(sizeof(int)); loc->nref[0]=0;
    loc->locx=mxGetPr(A);
    loc->nloc=mxGetM(A);
    loc->locy=loc->locx+loc->nloc;
    loc_dxdy(loc);
    /*const double tol=1e-7;
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
    loc->dy=dyd;*/
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
    if(A && mxGetNumberOfElements(A)){
	double *p=mxGetPr(A);
	//The user may supply 0 for empty matrix
	if(mxGetNumberOfElements(A)>1 || p[0]){
	    out=dnew_do(mxGetM(A), mxGetN(A), mxGetPr(A), 0);
	}
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
		out=cnew_do(mxGetM(A), mxGetN(A), (comp*)mxGetData(A), 0);
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
void *mx2any(const mxArray *A, void*(*fun)(const mxArray*)){
    if(!A) return NULL;
    else if(!mxIsCell(A)){
	if(fun){
	    warning("Requesting cell but input is not.\n");
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
dcell *mx2dcell(const mxArray *A){
    return (dcell*)mx2any(A, (void*(*)(const mxArray*))mx2d);
}
dcell *mx2ccell(const mxArray *A){
    return (dcell*)mx2any(A, (void*(*)(const mxArray*))mx2c);
}
lcell *mx2lcell(const mxArray *A){
    return (lcell*)mx2any(A, (void*(*)(const mxArray*))mx2l);
}
loccell *mx2loccell(const mxArray *A){
    return (loccell*)mx2any(A, (void*(*)(const mxArray*))mx2loc);
}
dcell *mx2dspcell(const mxArray *A){
    return (dcell*)mx2any(A, (void*(*)(const mxArray*))mx2dsp);
}
cell *mx2cell(const mxArray*A){
    return (cell*) mx2any(A, NULL);
}

static inline kalman_t *mx2kalman(const mxArray*A){
    kalman_t *kalman=(kalman_t*)calloc(1, sizeof(kalman_t));
	kalman->AdM=mx2d(mxGetField(A,0,"AdM"));
    kalman->BM=mx2d(mxGetField(A,0,"BM"));
    kalman->dthi=(double)mxGetScalar(mxGetField(A,0,"dthi"));
	kalman->dtrat_wfs=mx2l(mxGetField(A,0,"dtrat_wfs"));
	kalman->dtrats=mx2l(mxGetField(A,0,"dtrats"));
	kalman->Ad=mx2dcell(mxGetField(A,0,"Ad"));
    kalman->Cd=mx2dcell(mxGetField(A,0,"Cd"));
    kalman->Gwfs=mx2dcell(mxGetField(A,0,"Gwfs"));
    kalman->Cnn=mx2dcell(mxGetField(A,0,"Cnn"));
	kalman->M=mx2dcell(mxGetField(A,0,"M"));
    kalman->P=mx2dcell(mxGetField(A,0,"P"));
    return kalman;
}
static inline mxArray* kalman2mx(kalman_t *kalman){
    const int nfield=12;
    const char *fieldnames[]={"AdM","BM","dthi", "dtrat_wfs", "dtrats","Ad","Cd","Gwfs","Cnn","Rn","M","P"};
    mxArray *A=mxCreateStructMatrix(1,1,nfield,fieldnames);
    int pos=0;
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->AdM));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->BM));
    mxSetFieldByNumber(A, 0, pos++, mxCreateDoubleScalar(kalman->dthi));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->dtrat_wfs));
	mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->dtrats));
	mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Ad));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Cd));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Gwfs));
	mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Cnn));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->Rn));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->M));
    mxSetFieldByNumber(A, 0, pos++, any2mx(kalman->P));
    
    if(pos!=nfield){
	error("Invalid number of elements\n");
    }
    return A;
}
static inline mxArray *cn2est2mx(cn2est_t *cn2est){
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
static inline mxArray *dtf2mx(dtf_t *dtf){
    const int nfield=2;
    const char *fieldnames[]={"nominal","si"};
    mxArray *A=mxCreateStructMatrix(dtf->nwvl,1,nfield,fieldnames);
    for(int iwvl=0; iwvl<dtf->nwvl; iwvl++){
	int pos=0;
	mxSetFieldByNumber(A, iwvl, pos++, any2mx(dtf[iwvl].nominal));
	mxSetFieldByNumber(A, iwvl, pos++, any2mx(dtf[iwvl].si));
	if(pos!=nfield){
	    error("Invalid number of elements\n");
	}
    }
    return A;
}

static inline char *mx2str(const mxArray *A){
    int nlen=mxGetNumberOfElements(A)+1;
    char *fn=(char*)malloc(nlen);
    mxGetString(A, fn, nlen);
    return fn;
}
static inline rand_t *mx2rand(const mxArray *A){
    int seed=(int)mxGetScalar(A);
    rand_t *out=(rand_t*)malloc(sizeof(rand_t));
    seed_rand(out, seed);
    return out;
}


static inline void *calloc_mex(size_t nmemb, size_t size){
    void *p=mxCalloc(nmemb, size);
    mexMakeMemoryPersistent(p);
    return p;
}
static inline void *malloc_mex(size_t size){
    void *p=mxMalloc(size);
    mexMakeMemoryPersistent(p);
    return p;
}
static inline void *realloc_mex(void *p, size_t size){
    p=mxRealloc(p, size);
    mexMakeMemoryPersistent(p);
    return p;
}
static inline void free_mex(void*p){
    mxFree(p);
}
int mex_signal_handler(int sig){
    if(sig){
        mexErrMsgTxt("Signal caught.\n");
    } else{
        dbg("signal 0 caught\n");
    }
    return 1;
}
static __attribute__((constructor)) void init(){
    fprintf(stderr, "mex loaded\n");
    register_signal_handler(mex_signal_handler);
#if _OPENMP 
    NTHREAD=1;//MAOS with openmp does not play well with matlab. 
#endif
}
static __attribute__((destructor)) void deinit(){
    fprintf(stderr, "mex unloaded\n");
}
#endif
