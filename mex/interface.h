#ifndef AOS_MEX_INTERFACE_H
#define AOS_MEX_INTERFACE_H
#include "../lib/aos.h"
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif
#if __STDC_VERSION__ < 199901L /*C99*/
#define inline __inline
#endif

inline mxArray *dsp2mx(const dsp*A){
    mxArray *out=mxCreateSparse(A->m,A->n,A->nzmax,mxREAL);
    memcpy(mxGetIr(out),A->i,A->nzmax*sizeof(long));
    memcpy(mxGetJc(out),A->p,(A->n+1)*sizeof(long));
    memcpy(mxGetPr(out),A->x,A->nzmax*sizeof(double));
    return out;
}
inline dsp *mx2dsp(const mxArray *A){
    dsp *out=calloc(1, sizeof(dsp));
    out->p=mxGetJc(A);
    out->i=mxGetIr(A);
    out->x=mxGetPr(A);
    if(mxGetPi(A)){
	mexErrMsgTxt("A is complex");
    }
    return out;
}
inline mxArray *d2mx(const dmat *A){
    mxArray *out=mxCreateDoubleMatrix(A->nx,A->ny,mxREAL);
    memcpy(mxGetPr(out),A->p,A->nx*A->ny*sizeof(double));
    return out;
}
inline mxArray *c2mx(const cmat *A){
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

inline loc_t *mx2loc(const mxArray *A){
    loc_t *loc=calloc(1, sizeof(loc_t));
    loc->locx=mxGetPr(A);
    loc->nloc=mxGetM(A);
    loc->locy=loc->locx+loc->nloc;
    int found=0;
    int iloc;
    for(iloc=0; iloc<loc->nloc-2; iloc++){
	if(fabs(loc->locy[2+iloc]-loc->locy[iloc])<1.e-10 
	   && fabs(loc->locy[2+iloc]-loc->locy[1+iloc])<1.e-10
	   && fabs(loc->locx[2+iloc]+loc->locx[iloc]-loc->locx[1+iloc]*2)<1.e-10){
	    loc->dx=loc->locx[iloc+1]-loc->locx[iloc];
	    found=1;
	    break;
	}
    }
    if(!found){
	mexErrMsgTxt("Unable to determine dx");
    }
    return loc;
}

inline dmat *mx2d(const mxArray *A){
    if(mxGetPi(A)){
	mexErrMsgTxt("A is complex");
    }
    if(mxGetIr(A)){
	mexErrMsgTxt("A is dsp");
    }
    dmat *out=dnew_ref( mxGetM(A), mxGetN(A), mxGetPr(A));
    return out;
}
inline char *mx2str(const mxArray *A){
    int nlen=mxGetNumberOfElements(A)+1;
    char *fn=malloc(nlen);
    mxGetString(A, fn, nlen);
    return fn;
}

#endif
