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
	P_KEEPHT,
	P_L0,
	P_GRAD,
	P_TOT,
    };
    enum{
	PL_RES,
	PL_TOT,
    };
    if(nrhs !=P_TOT){
	mexErrMsgTxt("Usage: res=cn2est(wfspair, wfstheta, saloc, saa, saat, hs, ht, keepht, L0, grad)");
    }
    dmat *wfspair=mx2d(prhs[P_WFSPAIR]);
    dmat *wfstheta=mx2d(prhs[P_WFSTHETA]);
    if(maxabs(wfstheta->p, wfstheta->nx*wfstheta->ny)>1){
	dscale(wfstheta, 1./206265);
    }
    loc_t *saloc=mx2loc(prhs[P_SALOC]);
    dmat *saa=mx2dvec(prhs[P_SAA]);
    double saat=mxGetScalar(prhs[P_SAAT]);
    dmat *hs=mx2dvec(prhs[P_HS]);
    dmat *htrecon=mx2dvec(prhs[P_HT]);
    int keepht=(int)mxGetScalar(prhs[P_KEEPHT]);
    double L0=mxGetScalar(prhs[P_L0]);
    dcell *grad=mx2dcell(prhs[P_GRAD]);
    if(grad->nx==1 && grad->ny>1){
	grad->nx=grad->ny;
	grad->ny=1;
    }
    struct CN2EST_T *cn2est=cn2est_new(wfspair, wfstheta, saloc, saa, saat, hs, htrecon, keepht, L0);
    cn2est_push(cn2est, grad);
    cn2est_est(cn2est, 1, 0);
#define nfield 11
    const char *fieldnames[nfield]={"htrecon","wtrecon","r0m","ht","wt","r0","Pnk","iPnk","overlap","cov2","cov1"};
    int pos=0;
    plhs[0]=mxCreateStructMatrix(1,1,nfield,fieldnames);
    mxSetFieldByNumber(plhs[0], 0, pos++, d2mx(cn2est->htrecon));
    mxSetFieldByNumber(plhs[0], 0, pos++, d2mx(cn2est->wtrecon->p[0]));
    mxSetFieldByNumber(plhs[0], 0, pos++, mxCreateDoubleScalar(cn2est->r0m));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(cn2est->ht));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(cn2est->wt));
    mxSetFieldByNumber(plhs[0], 0, pos++, d2mx(cn2est->r0));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(cn2est->Pnk));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(cn2est->iPnk));
    mxSetFieldByNumber(plhs[0], 0, pos++, d2mx(cn2est->overlap));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(cn2est->cov2));
    mxSetFieldByNumber(plhs[0], 0, pos++, dcell2mx(cn2est->cov1));

    dfree(wfspair);
    dfree(wfstheta);
    free(saloc);
    dfree(saa);
    dfree(htrecon);
    dcellfree(grad);
}
