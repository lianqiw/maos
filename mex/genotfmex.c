#define USE_MEM 0
#include "../lib/genotf.h"
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_LOC,
	P_AMP,
	P_OPDBIAS,
	P_AREA,
	P_THRES,
	P_WVL,
	P_DTHETA,
	P_R0,
	P_L0,
	P_NCOMPX,
	P_NCOMPY,
	P_NSA,
	P_PTTR,
	P_NTHREAD,
	P_TOT,
    };
    enum{
	PL_OTF,
	PL_TOT,
    };
    if(nlhs!=PL_TOT || nrhs !=P_TOT){
	mexErrMsgTxt("Usage: otf=genotfmex(loc, amp, opdbias, area, thres, wvl, dtheta, r0, l0, ncompx, ncompy, nsa, pttr, nthread)");
    }
    LOC_T *loc=mx2loc(prhs[P_LOC]);
    double *amp=mxGetPr(prhs[P_AMP]);
    double *opdbias=mxGetPr(prhs[P_OPDBIAS]);
    double *area=mxGetPr(prhs[P_AREA]);
    double thres=mxGetScalar(prhs[P_THRES]);
    double wvl=mxGetScalar(prhs[P_WVL]);
    double dtheta=mxGetScalar(prhs[P_DTHETA]);
    double r0=mxGetScalar(prhs[P_R0]);
    double l0=mxGetScalar(prhs[P_L0]);
    const int ncompx=mxGetScalar(prhs[P_NCOMPX]);
    const int ncompy=mxGetScalar(prhs[P_NCOMPY]);
    const int nsa=mxGetScalar(prhs[P_NSA]);
    const int pttr=mxGetScalar(prhs[P_PTTR]);
    const int nthread=mxGetScalar(prhs[P_NTHREAD]);
    ccell *otf=ccellnew(nsa,1);
    genotf(otf->p, loc, amp, opdbias, area, thres, wvl, dtheta, r0, l0, ncompx, ncompy,
	   nsa, pttr, nthread);
    mwSize dims[2];
    dims[0]=nsa;
    dims[1]=1;
    plhs[PL_OTF]=mxCreateCellArray(2,dims);
    for(int isa=0; isa<nsa; isa++){
	mxArray *iotf=c2mx(otf->p[isa]);
	mxSetCell(plhs[PL_OTF],isa,iotf);
    }
}
