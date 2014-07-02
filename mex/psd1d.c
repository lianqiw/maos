/**
  Wrap of the function psd1d to mex routines.
*/
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_DATA,
	P_LSEG,
	P_TOT,
    };
    enum{
	PL_PSD,
	PL_TOT,
    }; 
    if(nrhs < P_TOT){
	mexErrMsgTxt("Usage: psd=psd1d(data, lseg)");
    }
    
    dmat *data=mx2d(prhs[P_DATA]);
    long lseg=(long)mxGetScalar(prhs[P_LSEG]);
    dmat *psd;
    if(nrhs > P_TOT){
	double dt=(double)mxGetScalar(prhs[P_TOT]);
	psd=psd1dt(data, lseg, dt);
    }else{
	psd=psd1d(data, lseg);
    }
    dfree(data);
    plhs[PL_PSD]=d2mx(psd);
    dfree(psd);
    dfree(data);
}
