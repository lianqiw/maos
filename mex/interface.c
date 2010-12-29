#include "interface.h"
mxArray *dsp2mx(const dsp*A){
    mxArray *out=mxCreateSparse(A->m,A->n,A->nzmax,mxREAL);
    memcpy(mxGetIr(out),A->i,A->nzmax*sizeof(long));
    memcpy(mxGetJc(out),A->p,(A->n+1)*sizeof(long));
    memcpy(mxGetPr(out),A->x,A->nzmax*sizeof(double));
    return out;
}
dsp *mx2dsp(const mxArray *A){
    dsp *out=calloc(1, sizeof(dsp));
    out->p=(long*)mxGetJc(A);
    out->i=(long*)mxGetIr(A);
    out->x=mxGetPr(A);
    if(mxGetPi(A)){
	mexErrMsgTxt("A is complex");
    }
    return out;
}
mxArray *d2mx(const dmat *A){
    mxArray *out=mxCreateDoubleMatrix(A->nx,A->ny,mxREAL);
    memcpy(mxGetPr(out),A->p,A->nx*A->ny*sizeof(double));
    return out;
}
mxArray *c2mx(const cmat *A){
    mxArray *out=mxCreateDoubleMatrix(A->nx, A->ny, mxCOMPLEX);
    double *pr=mxGetPr(out);
    double *pi=mxGetPi(out);
    for(long i=0; i<A->nx*A->ny; i++){
	pr[i]=creal(A->p[i]);
	pi[i]=cimag(A->p[i]);
    }
    return out;
}

loc_t *mx2loc(const mxArray *A){
    loc_t *loc=calloc(1, sizeof(loc_t));
    loc->locx=mxGetPr(A);
    loc->nloc=mxGetM(A);
    loc->locy=loc->locx+loc->nloc;
    if(fabs(loc->locx[2]+loc->locx[0]-loc->locx[1]*2)<1.e-10){
	loc->dx=loc->locx[1]-loc->locx[0];
    }else{
	mexErrMsgTxt("Unable to determine dx");
    }
    return loc;
}
dmat *mx2d(const mxArray *A){
    if(mxGetPi(A)){
	mexErrMsgTxt("A is complex");
    }
    if(mxGetIr(A)){
	mexErrMsgTxt("A is dsp");
    }
    dmat *out=dnew_ref(mxGetPr(A), mxGetM(A), mxGetN(A));
    return out;
}
char *mx2str(const mxArray *A){
    int nlen=mxGetNumberOfElements(A)+1;
    char *fn=malloc(nlen);
    mxGetString(A, fn, nlen);
    return fn;
}
