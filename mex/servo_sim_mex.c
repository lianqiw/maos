#include <mex.h>
#include <math.h>
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_INPUT,
	P_DT,
	P_DTRAT,
	P_SIGMAN,
	P_GAIN,
	P_TOT,
    };
    enum{
	PL_RES,
	PL_TOT,
    };
    
    if(nlhs!=PL_TOT || nrhs!=P_TOT){
	mexErrMsgTxt("Usage: gain=servo_test(input, dt, dtrat, sigman, gain);\n"
		     "input is input time series. Should match the PSD used for gain optimization\n"
		     "dt is the AO fundemental sampling period.\n"
		     "dtrat is ratio of sampling period of the WFS over dt.\n"
		     "sigman is the variance due to noise\n"
		     "gain is obtained from servo_optim\n");
    }
    dmat *input = mx2d(prhs[P_INPUT]);
    double dt   = mxGetScalar(prhs[P_DT]);
    int dtrat  = (int)mxGetScalar(prhs[P_DTRAT]);
    dmat *sigma2n = mx2d(prhs[P_SIGMAN]);/*m^2 */
    dmat *gain = mx2d(prhs[P_GAIN]);/*m^2 */
    dmat *res   = servo_test(input, dt, dtrat, sigma2n, gain);
    plhs[PL_RES]= d2mx(res);
    dfree(res);
    dfree(input);
    dfree(sigma2n);
    dgree(gain);
}
