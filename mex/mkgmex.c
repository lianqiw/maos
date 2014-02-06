#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_XLOC,/*XLOC is where phase points are defined. XLOC==PLOC if on same plan*/
	P_DX,
	P_PLOC,/*PLOC is where AMP is defined.*/
	P_DP,
	P_AMP,
	P_SALOC,
	P_DSA,
	P_SCALE,
	P_DISPLACE,
	P_DOPARTIAL,
	P_TOT,
    };
    enum{
	PL_H,
	PL_TOT,
    };
    int do_partial;
    loc_t xloc,ploc,saloc;
    double *amp;
    double scale, *displace;
    if(nrhs!=P_TOT || nlhs!=PL_TOT){
	mexErrMsgTxt("Usage: G=mkgmex(xloc, dx, ploc, dp, amp, saloc, dsa, scale, displace[2], dopartial)\n");
    }
    if(mxGetNumberOfElements(prhs[P_DISPLACE])!=2){
	mexErrMsgTxt("Displace needs to be a two-vector.\n");
    }
#define LOC_FROM_MATLAB(A,B)						\
    {A.nloc=mxGetM(B); A.locx=mxGetPr(B); A.locy=A.locx+A.nloc;A.map=NULL;}
    LOC_FROM_MATLAB(xloc,prhs[P_XLOC]);
    LOC_FROM_MATLAB(ploc,prhs[P_PLOC]);
    LOC_FROM_MATLAB(saloc,prhs[P_SALOC]);
#undef LOC_FROM_MATLAB
    amp=mxGetPr(prhs[P_AMP]);
    xloc.dx=mxGetScalar(prhs[P_DX]);
    ploc.dx=mxGetScalar(prhs[P_DP]);
    saloc.dx=mxGetScalar(prhs[P_DSA]);
    
    /*scale and displace is between XLOC and PLOC->*/
    scale=mxGetScalar(prhs[P_SCALE]);
    displace=mxGetPr(prhs[P_DISPLACE]);
    do_partial=(int)mxGetScalar(prhs[P_DOPARTIAL]);
    dsp *GS0T=mkgt(&xloc, &ploc, amp, &saloc, 
		   1, scale, displace, do_partial);
    mxArray *GS0T2=mxCreateSparse(GS0T->m, GS0T->n, GS0T->nzmax, mxREAL);
    memcpy(mxGetPr(GS0T2), GS0T->x, GS0T->nzmax*sizeof(double));
    memcpy(mxGetIr(GS0T2), GS0T->i, GS0T->nzmax*sizeof(long));
    memcpy(mxGetJc(GS0T2), GS0T->p, (GS0T->n+1)*sizeof(long));
    spfree(GS0T);
    mexCallMATLAB(1,plhs,1, &GS0T2,"tanspose");
}
