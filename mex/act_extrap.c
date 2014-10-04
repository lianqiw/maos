/**
  Wrap of the function denc to mex routines.
*/
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_LOC,
	P_CPL,
	P_THRES,
	P_TOT,
    };
    enum{
	PL_H,
	PL_TOT,
    };
    if(nrhs !=P_TOT){
	mexErrMsgTxt("Usage: H=act_extrap(loc, cpl, thres)\n"
		     "loc: is nx2 coordinate of actuators\n"
		     "cpl: is nx1 coupling of each actuator\n"
		     "thres: when cpl below thres do extrapolation onto this act\n"
		     );
    }
    loc_t *loc=mx2loc(prhs[P_LOC]);
    dmat *cpl=mx2d(prhs[P_CPL]);
    double thres=(double)mxGetScalar(prhs[P_THRES]);
    if(cpl->nx==1 && cpl->ny>1){
	cpl->nx=cpl->ny;
	cpl->ny=1;
    }
    if(cpl->nx!=loc->nloc){
	error("loc and cpl does not match\n");
    }
    dsp *H=act_extrap_do(loc, cpl, thres);
    plhs[PL_H]=dsp2mx(H);
    dspfree(H);
    dfree(cpl);
    loc->locx=loc->locy=0;
    locfree(loc);
}
