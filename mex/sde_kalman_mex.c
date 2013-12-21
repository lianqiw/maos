#include <mex.h>
#include <math.h>
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_COEFF,
	P_DTHI,
	P_DTRAT,
	P_GWFS,
	P_RWFS,
	P_PROJ,
	P_TOT,
    };
    enum{
	PL_KALMAN,
	PL_TOT,
    };
    
    if(nlhs!=PL_TOT || nrhs!=P_TOT){
	mexErrMsgTxt("Usage: kalman=sde_kalman_mex(coeff, dthi, dtrat, Gwfs, Rwfs)\n");
    }
    dmat *coeff  = mx2d(prhs[P_COEFF]);
    double dthi  = (double)mxGetScalar(prhs[P_DTHI]);
    double dtrat = (double)mxGetScalar(prhs[P_DTRAT]);
    dmat *Gwfs  = mx2d(prhs[P_GWFS]);
    dmat *Rwfs  = mx2d(prhs[P_RWFS]);
    dmat *Proj  = mx2d(prhs[P_PROJ]);
    kalman_t *kalman=sde_kalman(coeff, dthi, dtrat, Gwfs, Rwfs, Proj);
 
    int nfield=11;
    const char *fieldnames[]={"Ad","Fd","Cd","AdM","FdM","M","P", "dthi", "dtrat", "Gwfs", "Rwfs"};
    plhs[0]=mxCreateStructMatrix(1,1,nfield,fieldnames);
    mxSetFieldByNumber(plhs[0], 0, 0, d2mx(kalman->Ad));
    mxSetFieldByNumber(plhs[0], 0, 1, d2mx(kalman->Fd));
    mxSetFieldByNumber(plhs[0], 0, 2, d2mx(kalman->Cd));
    mxSetFieldByNumber(plhs[0], 0, 3, d2mx(kalman->AdM));
    mxSetFieldByNumber(plhs[0], 0, 4, d2mx(kalman->FdM));
    mxSetFieldByNumber(plhs[0], 0, 5, d2mx(kalman->M));
    mxSetFieldByNumber(plhs[0], 0, 6, d2mx(kalman->P));
    mxSetFieldByNumber(plhs[0], 0, 7, mxDuplicateArray(prhs[P_DTHI]));
    mxSetFieldByNumber(plhs[0], 0, 8, mxDuplicateArray(prhs[P_DTRAT]));
    mxSetFieldByNumber(plhs[0], 0, 9, mxDuplicateArray(prhs[P_GWFS]));
    mxSetFieldByNumber(plhs[0], 0, 10,mxDuplicateArray(prhs[P_RWFS]));
    
    kalman_free(kalman);
}
