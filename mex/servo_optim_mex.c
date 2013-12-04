#include <mex.h>
#include <math.h>
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_PSD,
	P_DT,
	P_DTRAT,
	P_SIGMAN,
	P_STYPE,
	P_TOT,
    };
    enum{
	PL_GAIN,
	PL_TOT,
    };

    if(nlhs!=PL_TOT || nrhs!=P_TOT){
	mexErrMsgTxt("Usage: gain=servo_optim(psd, dt, dtrat, sigman, servotype);\n"
		     "PSD psd should be in m^2/hz\n"
		     "dt is the AO fundemental sampling period.\n"
		     "dtrat is ratio of sampling period of the WFS over dt.\n"
		     "sigman is wavefront variance due to noise in m^2.\n"
		     "servotype is 1 or 2 for typeI, typeII controller.\n"
		     );
    }
    dmat *psd  = mx2d(prhs[P_PSD]);
    double dt  = mxGetScalar(prhs[P_DT]);
    long dtrat = (long)mxGetScalar(prhs[P_DTRAT]);
    double sigma = mxGetScalar(prhs[P_SIGMAN]);/*m^2 */
    dmat *sigma2 = dnew(1,1); 
    sigma2->p[0] = sigma;
    int servotype = (int)mxGetScalar(prhs[P_STYPE]);
    dcell *gain   = servo_optim(psd,dt,dtrat,M_PI/4,sigma2,servotype);
    plhs[PL_GAIN]  = d2mx(gain->p[0]);
    dfree(sigma2);
    dcellfree(gain);
}
