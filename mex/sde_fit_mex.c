#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_PSD,
	P_COEFF,
	P_TMAX,
	P_MIN,
	P_MAX,
	P_DF,
	P_TOT,
    };
    enum{
	PL_COEFF,
	PL_TOT,
    };
    
    if(nrhs<P_MIN){
	mexErrMsgTxt("Usage: coeff=sde_fit_mex(psd, coeff0, tmax[, min, max, df])");
    }
    dmat *psd  = mx2d(prhs[P_PSD]);
    dmat *coeff0  = mx2d(prhs[P_COEFF]);
    double tmax  = (double)mxGetScalar(prhs[P_TMAX]);
    double min=0, max=INFINITY, df=0;
    if(nrhs>P_MAX){
	min  = (double)mxGetScalar(prhs[P_MIN]);
	max  = (double)mxGetScalar(prhs[P_MAX]);
    }
    if(nrhs > P_DF){
	df = (double)mxGetScalar(prhs[P_DF]);
    }
    dmat *coeff=sde_fit(psd, coeff0, tmax, min, max, df);
    plhs[PL_COEFF]=d2mx(coeff);
    dfree(coeff);
    dfree(psd);
    dfree(coeff0);
}
