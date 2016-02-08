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
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    enum{/*input */
	P_LOC,
	P_R,
	P_TOT,
    };
    enum{/*output */
	PL_W0,
	PL_W1,
	PL_TOT,
    };
    if(nlhs!=PL_TOT || nrhs!=P_TOT){
	mexErrMsgTxt("Usage: [W0, W1]=mkwmex(loc, Radius)\n");
    }
    loc_t *loc=mx2loc(prhs[P_LOC]);
    double R=mxGetScalar(prhs[P_R]);
    dsp *W0;
    dmat *W1;
    mkw_circular(loc,0,0,R,&W0,&W1);
    plhs[PL_W0]=dsp2mx(W0);
    dspfree(W0);
    plhs[PL_W1]=d2mx(W1);
    dfree(W1);
}
