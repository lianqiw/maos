/**
  Wrap of the function mk2dotf to mex routines.
*/
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_LOC,
	P_AMP,
	P_COV,
	P_NORM,
	P_TOT,
    };
    enum{
	PL_COV2D,
	PL_TOT,
    }; 
    if(nrhs !=P_TOT){
	mexErrMsgTxt("Usage: cov2d=mk2dcov(loc, amp, cov, normalization)");
    }
    loc_t *loc=mx2loc(prhs[P_LOC]);
    dmat *amp0=mx2d(prhs[P_AMP]);
    dmat *amp=ddup(amp0);
    dfree(amp0);
    dmat *cov=mx2d(prhs[P_COV]);
    int norm=(int)mxGetScalar(prhs[P_NORM]);
    dmat *cov2d=NULL;
    double *pamp=NULL;
    if(amp && amp->nx==loc->nloc){
	pamp=amp->p;
	normalize_max(pamp, loc->nloc, 1);
    }
    mk2dcov(&cov2d, loc, pamp, 0.5, cov, norm);
    plhs[PL_COV2D]=d2mx(cov2d);
    dfree(cov);
    dfree(cov2d);
    dfree(amp);
    free(loc);
}
