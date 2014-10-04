#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_XLOC,/*XLOC is where phase points are defined. XLOC==PLOC if on same plan*/
	P_PLOC,/*PLOC is where AMP is defined.*/
	P_AMP,
	P_SALOC,
	P_SCALE,
	P_DISPX,
	P_DISPY,
	P_DOPARTIAL,
	P_TOT,
    };
    enum{
	PL_H,
	PL_TOT,
    };
    if(nrhs!=P_TOT){
	mexErrMsgTxt("Usage: G=mkgmex(xloc, ploc, amp, saloc, scale, dispx, dispy, dopartial)\n");
    }
    loc_t *xloc=mx2loc(prhs[P_XLOC]);
    loc_t *ploc=mx2loc(prhs[P_PLOC]);
    loc_t *saloc=mx2loc(prhs[P_SALOC]);
    double *amp=mxGetPr(prhs[P_AMP]);
    double dispx=mxGetScalar(prhs[P_DISPX]);
    double dispy=mxGetScalar(prhs[P_DISPY]);
    double scale=mxGetScalar(prhs[P_SCALE]);
    int do_partial=(int)mxGetScalar(prhs[P_DOPARTIAL]);
    dsp *GS0=mkg(xloc, ploc, amp, saloc, 1, scale, dispx, dispy, do_partial);
    plhs[0]=dsp2mx(GS0);
    dspfree(GS0);
    free(xloc);
    free(ploc);
    free(saloc);
}
