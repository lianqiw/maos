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
    enum{
	P_INPUT,
	P_KALMAN,
	P_TOT,
    };
    enum{
	PL_RES,
	PL_TOT,
    };
    (void)nlhs;
    if(nrhs!=P_TOT){
	mexErrMsgTxt("Usage: res=kalman_sim_mex(input, kalman)\n");
    }
    dmat *input=mx2d(prhs[P_INPUT]);
    kalman_t *kalman=calloc(1, sizeof(kalman_t));
    kalman->Ad=mx2d(mxGetField(prhs[P_KALMAN],0,"Ad"));
    kalman->Cd=mx2dcell(mxGetField(prhs[P_KALMAN],0,"Cd"));
    kalman->AdM=mx2d(mxGetField(prhs[P_KALMAN],0,"AdM"));
    kalman->FdM=mx2d(mxGetField(prhs[P_KALMAN],0,"FdM"));
    kalman->M=mx2dcell(mxGetField(prhs[P_KALMAN],0,"M"));
    kalman->P=mx2d(mxGetField(prhs[P_KALMAN],0,"P"));
    kalman->dthi=(double)mxGetScalar(mxGetField(prhs[P_KALMAN],0,"dthi"));
    kalman->dtrat=mx2d(mxGetField(prhs[P_KALMAN],0,"dtrat"));
    kalman->Gwfs=mx2dcell(mxGetField(prhs[P_KALMAN],0,"Gwfs"));
    kalman->Rwfs=mx2dcell(mxGetField(prhs[P_KALMAN],0,"Rwfs"));
    dmat *res=kalman_test(kalman, input);
    kalman_free(kalman);
    dfree(input);
    plhs[PL_RES]=d2mx(res);
    dfree(res);
}
