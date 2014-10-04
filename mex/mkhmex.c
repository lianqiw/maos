#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_LOCIN,/*XLOC is where phase points are defined. XLOC==PLOC if on same plan*/
	P_LOCOUT,
	P_AMPOUT,
	P_DISPX,
	P_DISPY,
	P_SCALE,
	P_CUBIC,
	P_CUBIC_IAC,
	P_TOT,
    };
    enum{
	PL_H,
	PL_TOT,
    };
    if(nrhs!=P_TOT){
	mexErrMsgTxt("Usage: H=mkgmex(locin, locout, ampout, dispx, dispy, scale, cubic, cubic_iac)\n");
    }
    loc_t *locin=mx2loc(prhs[P_LOCIN]);
    loc_t *locout=mx2loc(prhs[P_LOCOUT]);
    dmat *ampout=mx2d(prhs[P_AMPOUT]);
    double dispx=mxGetScalar(prhs[P_DISPX]);
    double dispy=mxGetScalar(prhs[P_DISPY]);
    double scale=mxGetScalar(prhs[P_SCALE]);
    int cubic=(int)mxGetScalar(prhs[P_CUBIC]);
    double cubic_iac=mxGetScalar(prhs[P_CUBIC_IAC]);
    if(ampout && ampout->nx!=locout->nloc){
	error("ampout invalid\n");
    }
    dsp *H=mkh(locin, locout, ampout?ampout->p:0, dispx, dispy, scale, cubic, cubic_iac);
    plhs[0]=dsp2mx(H);
    free(locin);
    free(locout);
    dfree(ampout);
    dspfree(H);
}
