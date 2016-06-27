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
  Wrap of the function denc to mex routines.
*/
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{
	P_LOC,
	P_CPL,
	P_THRES,
	P_TOT,
    };
    enum{
	PL_H,
	PL_TOT,
    };
    (void)nlhs;
    if(nrhs !=P_TOT){
	mexErrMsgTxt("Usage: H=act_extrap(loc, cpl, thres)\n"
		     "loc: is nx2 coordinate of actuators\n"
		     "cpl: is nx1 coupling of each actuator\n"
		     "thres: when cpl below thres do extrapolation onto this act\n"
		     );
    }
    loc_t *loc=mx2loc(prhs[P_LOC]);
    dmat *cpl=mx2d(prhs[P_CPL]);
    double thres=(double)mxGetScalar(prhs[P_THRES]);
    if(cpl->nx==1 && cpl->ny>1){
	cpl->nx=cpl->ny;
	cpl->ny=1;
    }
    if(cpl->nx!=loc->nloc){
	error("loc and cpl does not match\n");
    }
    dsp *H=act_extrap_do(loc, cpl, thres);
    plhs[PL_H]=dsp2mx(H);
    dspfree(H);
    dfree(cpl);
    loc->locx=loc->locy=0;
    locfree(loc);
}
