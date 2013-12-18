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
	P_TOT,
    };
    enum{
	PL_AD,
	PL_FD,
	PL_CD,
	PL_ADM,
	PL_FDM,
	PL_M,
	PL_TOT,
    };
    
    if((nlhs!=PL_TOT && nlhs!=1) || nrhs!=P_TOT){
	mexErrMsgTxt("Usage: [Ad, Fd, Cd, AdM, FdM, M]=sde_kalman_mex(coeff, dthi, dtrat, Gwfs, Rwfs)\n"
		     "or kalman=sde_kalman_mex(coeff, dthi, dtrat, Gwfs, Rwfs)");
    }
    dmat *coeff  = mx2d(prhs[P_COEFF]);
    double dthi  = (double)mxGetScalar(prhs[P_DTHI]);
    double dtrat = (double)mxGetScalar(prhs[P_DTRAT]);
    dmat *Gwfs  = mx2d(prhs[P_GWFS]);
    dmat *Rwfs  = mx2d(prhs[P_RWFS]);
    kalman_t *kalman=sde_kalman(coeff, dthi, dtrat, Gwfs, Rwfs);
    if(nlhs==PL_TOT){
	plhs[PL_AD]=d2mx(kalman->Ad);
	plhs[PL_FD]=d2mx(kalman->Fd);
	plhs[PL_CD]=d2mx(kalman->Cd);
	plhs[PL_ADM]=d2mx(kalman->AdM);
	plhs[PL_FDM]=d2mx(kalman->FdM);
	plhs[PL_M]=d2mx(kalman->M);
    }else if(nlhs==1){//create a struct
	int nfield=7;
	const char *fieldnames[]={"Ad","Fd","Cd","AdM","FdM","M","P"};
	plhs[0]=mxCreateStructMatrix(1,1,nfield,fieldnames);
	mxSetFieldByNumber(plhs[0], 0, 0, d2mx(kalman->Ad));
	mxSetFieldByNumber(plhs[0], 0, 1, d2mx(kalman->Fd));
	mxSetFieldByNumber(plhs[0], 0, 2, d2mx(kalman->Cd));
	mxSetFieldByNumber(plhs[0], 0, 3, d2mx(kalman->AdM));
	mxSetFieldByNumber(plhs[0], 0, 4, d2mx(kalman->FdM));
	mxSetFieldByNumber(plhs[0], 0, 5, d2mx(kalman->M));
	mxSetFieldByNumber(plhs[0], 0, 6, d2mx(kalman->P));
    }
    kalman_free(kalman);
}
