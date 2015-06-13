/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
	P_XLOC,/*XLOC is where phase points are defined. XLOC==PLOC if on same plan*/
	P_PLOC,/*PLOC is where AMP is defined.*/
	P_AMP,
	P_SALOC,
	P_SCALE,
	P_DISPX,
	P_DISPY,
	P_DOPARTIAL,
	P_TOT,
    };
    enum{
	PL_H,
	PL_TOT,
    };
    if(nrhs!=P_TOT){
	mexErrMsgTxt("Usage: G=mkgmex(xloc, ploc, amp, saloc, scale, dispx, dispy, dopartial)\n");
    }
    loc_t *xloc=mx2loc(prhs[P_XLOC]);
    loc_t *ploc=mx2loc(prhs[P_PLOC]);
    loc_t *saloc=mx2loc(prhs[P_SALOC]);
    dmat *amp=mx2d(prhs[P_AMP]);
    double dispx=mxGetScalar(prhs[P_DISPX]);
    double dispy=mxGetScalar(prhs[P_DISPY]);
    double scale=mxGetScalar(prhs[P_SCALE]);
    int do_partial=(int)mxGetScalar(prhs[P_DOPARTIAL]);
    dsp *GS0=mkg(xloc, ploc, amp, saloc, 1, scale, dispx, dispy, do_partial);
    plhs[0]=dsp2mx(GS0);
    dspfree(GS0);
    free(xloc);
    free(ploc);
    free(saloc);
    dfree(amp);
}
