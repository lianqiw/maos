/**
  Wrap of the function genotf to mex routines.
*/
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_WFSPAIR,
	P_WFSTHETA,
	P_SALOC,
	P_SAA,
	P_SAAT,
	P_HS,
	P_HT,
	P_HMAX,
	P_KEEPHT,
	P_L0,
	P_GRAD,
	P_TOT,
    };
    enum{
	PL_R0,
	PL_WT,
	PL_TOT,
    };
    if(nlhs!=PL_TOT || nrhs !=P_TOT){
	mexErrMsgTxt("Usage: [r0, wt]=cn2est(wfspair, wfstheta, saloc, saa, saat, hs, ht, hmax, keepht, L0, grad)");
    }
    dmat *wfspair=mx2d(prhs[P_WFSPAIR]);
    dmat *wfstheta=mx2d(prhs[P_WFSTHETA]);
    if(maxabs(wfstheta->p, wfstheta->nx*wfstheta->ny)>1){
	dscale(wfstheta, 1./206265);
    }
    loc_t *saloc=mx2loc(prhs[P_SALOC]);
    dmat *saa=mx2dvec(prhs[P_SAA]);
    double saat=mxGetScalar(prhs[P_SAAT]);
    double hs=mxGetScalar(prhs[P_HS]);
    dmat *htrecon=mx2dvec(prhs[P_HT]);
    double hmax=mxGetScalar(prhs[P_HMAX]);
    int keepht=(int)mxGetScalar(prhs[P_KEEPHT]);
    double L0=mxGetScalar(prhs[P_L0]);
    dcell *grad=mx2dcell(prhs[P_GRAD]);
    if(grad->nx==1 && grad->ny>1){
	grad->nx=grad->ny;
	grad->ny=1;
    }
    struct CN2EST_T *cn2est=cn2est_new(wfspair, wfstheta, saloc, saa, saat, hs, htrecon, hmax, keepht, L0);
    cn2est_push(cn2est, grad);
    cn2est_est(cn2est, 1, 0);
    plhs[PL_R0]=mxCreateDoubleScalar(cn2est->r0m);
    plhs[PL_WT]=d2mx(cn2est->wtrecon->p[0]);

    dfree(wfspair);
    dfree(wfstheta);
    free(saloc);
    dfree(saa);
    dfree(htrecon);
    dcellfree(grad);
}
