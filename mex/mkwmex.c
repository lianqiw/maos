#include <mex.h>
#include <math.h>
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{/*input */
	P_LOC,
	P_R,
	P_TOT,
    };
    enum{/*output */
	PL_W0,
	PL_W1,
	PL_TOT,
    };
    if(nlhs!=PL_TOT || nrhs!=P_TOT){
	mexErrMsgTxt("Usage: [W0, W1]=mkwmex(loc, Radius)\n");
    }
    loc_t *loc=mx2loc(prhs[P_LOC]);
    double R=mxGetScalar(prhs[P_R]);
    dsp *W0;
    dmat *W1;
    mkw_circular(loc,0,0,R,&W0,&W1);
    plhs[PL_W0]=dsp2mx(W0);
    spfree(W0);
    plhs[PL_W1]=d2mx(W1);
    dfree(W1);
}
