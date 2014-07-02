#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_INPUT,
	P_KALMAN,
	P_TOT,
    };
    enum{
	PL_RES,
	PL_TOT,
    };
    
    if(nrhs!=P_TOT){
	mexErrMsgTxt("Usage: res=kalman_sim_mex(input, kalman)\n");
    }
    dmat *input=mx2d(prhs[P_INPUT]);
    kalman_t *kalman=calloc(1, sizeof(kalman_t));
    kalman->Ad=mx2d(mxGetField(prhs[P_KALMAN],0,"Ad"));
    kalman->Cd=mx2dcell(mxGetField(prhs[P_KALMAN],0,"Cd"));
    kalman->AdM=mx2d(mxGetField(prhs[P_KALMAN],0,"AdM"));
    kalman->FdM=mx2d(mxGetField(prhs[P_KALMAN],0,"FdM"));
    kalman->M=mx2dcell(mxGetField(prhs[P_KALMAN],0,"M"));
    kalman->P=mx2d(mxGetField(prhs[P_KALMAN],0,"P"));
    kalman->dthi=(double)mxGetScalar(mxGetField(prhs[P_KALMAN],0,"dthi"));
    kalman->dtrat=mx2d(mxGetField(prhs[P_KALMAN],0,"dtrat"));
    kalman->Gwfs=mx2dcell(mxGetField(prhs[P_KALMAN],0,"Gwfs"));
    kalman->Rwfs=mx2dcell(mxGetField(prhs[P_KALMAN],0,"Rwfs"));
    dmat *res=kalman_test(kalman, input);
    kalman_free(kalman);
    dfree(input);
    plhs[PL_RES]=d2mx(res);
    dfree(res);
}
