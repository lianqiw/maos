/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
  Wrap of the function genotf to mex routines.
*/
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
	P_COV,
	P_R0,
	P_L0,
	P_NCOMPX,
	P_NCOMPY,
	P_NSA,
	P_PTTR,
	P_TOT,
    };
    enum{
	PL_OTF,
	PL_TOT,
    };
    if(nrhs !=P_TOT){
	mexErrMsgTxt("Usage: otf=genotfmex(loc, amp, opdbias, area, thres, wvl, "
		     "dtheta, cov, r0, l0, ncompx, ncompy, nsa, pttr)\n"
		     "loc is for one subaperture. amp is for all subapertures.");
    }
    loc_t *loc=mx2loc(prhs[P_LOC]);
    dmat *amp=mx2d(prhs[P_AMP]);
    dmat *opdbias=mx2d(prhs[P_OPDBIAS]);
    dmat *area=mx2d(prhs[P_AREA]);
    double thres=mxGetScalar(prhs[P_THRES]);
    double wvl=mxGetScalar(prhs[P_WVL]);
    double dtheta=mxGetScalar(prhs[P_DTHETA]);
    dmat *cov=prhs[P_COV]?mx2d(prhs[P_COV]):NULL;
    double r0=mxGetScalar(prhs[P_R0]);
    double l0=mxGetScalar(prhs[P_L0]);
    const int ncompx=mxGetScalar(prhs[P_NCOMPX]);
    const int ncompy=mxGetScalar(prhs[P_NCOMPY]);
    const int nsa=mxGetScalar(prhs[P_NSA]);
    const int pttr=mxGetScalar(prhs[P_PTTR]);
    ccell *otf=cellnew(nsa,1);
    genotf(otf->p, loc, amp, opdbias, area, thres, wvl, dtheta, cov, r0, l0, ncompx, ncompy, nsa, pttr);
    mwSize dims[2];
    dims[0]=nsa;
    dims[1]=1;
    if(nsa>1){
	plhs[PL_OTF]=mxCreateCellArray(2,dims);
	int isa;
	for(isa=0; isa<nsa; isa++){
	    mxArray *iotf=c2mx(otf->p[isa]);
	    mxSetCell(plhs[PL_OTF],isa,iotf);
	}
    }else{
	plhs[PL_OTF]=c2mx(otf->p[0]);
    }
    dfree(cov);
    ccellfree(otf);
}
