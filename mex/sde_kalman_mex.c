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
    dmat *dtrat  = mx2d(prhs[P_DTRAT]);
    dcell *Gwfs  = mx2dcell(prhs[P_GWFS]);
    dcell *Rwfs  = mx2dcell(prhs[P_RWFS]);
    dmat *Proj   = mx2d(prhs[P_PROJ]);
    kalman_t *kalman=sde_kalman(coeff, dthi, dtrat, Gwfs, Rwfs, Proj);
    int nfield=12;
    const char *fieldnames[]={"Ad","Cd","AdM","FdM","Qn","Rn","M","P", "dthi", "dtrat", "Gwfs", "Rwfs"};
    plhs[0]=mxCreateStructMatrix(1,1,nfield,fieldnames);
    int pos=0;
    mxSetFieldByNumber(plhs[0], 0, pos++, d2mx(kalman->Ad));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(kalman->Cd));
    mxSetFieldByNumber(plhs[0], 0, pos++, d2mx(kalman->AdM));
    mxSetFieldByNumber(plhs[0], 0, pos++, d2mx(kalman->FdM));
    mxSetFieldByNumber(plhs[0], 0, pos++, d2mx(kalman->Qn));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(kalman->Rn));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(kalman->M));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(kalman->P));
    mxSetFieldByNumber(plhs[0], 0, pos++, mxDuplicateArray(prhs[P_DTHI]));
    mxSetFieldByNumber(plhs[0], 0, pos++, mxDuplicateArray(prhs[P_DTRAT]));
    mxSetFieldByNumber(plhs[0], 0, pos++, mxDuplicateArray(prhs[P_GWFS]));
    mxSetFieldByNumber(plhs[0], 0, pos++, mxDuplicateArray(prhs[P_RWFS]));
    
    kalman_free(kalman);
}
