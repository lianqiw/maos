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
#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_LOCIN,/*XLOC is where phase points are defined. XLOC==PLOC if on same plan*/
	P_LOCOUT,
	P_AMPOUT,
	P_DISPX,
	P_DISPY,
	P_SCALE,
	P_CUBIC_IAC,
	P_TOT,
    };
    enum{
	PL_H,
	PL_TOT,
    };
    if(nrhs!=P_TOT){
	mexErrMsgTxt("Usage: H=mkgmex(locin, locout, ampout, dispx, dispy, scale, cubic_iac)\n");
    }
    loc_t *locin=mx2loc(prhs[P_LOCIN]);
    loc_t *locout=mx2loc(prhs[P_LOCOUT]);
    dmat *ampout=mx2d(prhs[P_AMPOUT]);
    double dispx=mxGetScalar(prhs[P_DISPX]);
    double dispy=mxGetScalar(prhs[P_DISPY]);
    double scale=mxGetScalar(prhs[P_SCALE]);
    double cubic_iac=mxGetScalar(prhs[P_CUBIC_IAC]);
    if(ampout && ampout->nx!=locout->nloc){
	error("ampout invalid\n");
    }
    dsp *H=mkh(locin, locout, ampout?ampout->p:0, dispx, dispy, scale, cubic_iac);
    plhs[0]=dsp2mx(H);
    free(locin);
    free(locout);
    dfree(ampout);
    dspfree(H);
}
