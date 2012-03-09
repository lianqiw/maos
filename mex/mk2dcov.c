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
	P_TOT,
    };
    enum{
	PL_COV2D,
	PL_TOT,
    }; 
    if(nlhs!=PL_TOT || nrhs !=P_TOT){
	mexErrMsgTxt("Usage: cov2d=mk2dcov(loc, amp, cov)");
    }
    loc_t *loc=mx2loc(prhs[P_LOC]);
    dmat *amp=mx2d(prhs[P_AMP]);
    dmat *cov=mx2d(prhs[P_COV]);
    dmat *cov2d=NULL;
    double *pamp=NULL;
    if(amp && amp->nx==loc->nloc){
	pamp=amp->p;
	normalize_max(pamp, loc->nloc, 1);
    }
    mk2dcov(&cov2d, loc, pamp, 0.5, cov);
    plhs[PL_COV2D]=d2mx(cov2d);
}
