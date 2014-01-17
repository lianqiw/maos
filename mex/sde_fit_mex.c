#include <mex.h>
#include <math.h>
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_PSD,
	P_COEFF,
	P_TMAX,
	P_MIN,
	P_MAX,
	P_TOT,
    };
    enum{
	PL_COEFF,
	PL_TOT,
    };
    
    if(nlhs!=PL_TOT || (nrhs!=P_TOT && nrhs!=P_MIN)){
	mexErrMsgTxt("Usage: coeff=sde_fit_mex(psd, coeff0, tmax[, min, max])");
    }
    dmat *psd  = mx2d(prhs[P_PSD]);
    dmat *coeff0  = mx2d(prhs[P_COEFF]);
    double tmax  = (double)mxGetScalar(prhs[P_TMAX]);
    double min=0, max=INFINITY;
    if(nrhs==P_TOT){
	min  = (double)mxGetScalar(prhs[P_MIN]);
	max  = (double)mxGetScalar(prhs[P_MAX]);
    }
    dmat *coeff=sde_fit(psd, coeff0, tmax, min, max);
    plhs[PL_COEFF]=d2mx(coeff);
    dfree(coeff);
    dfree(psd);
    dfree(coeff0);
}
