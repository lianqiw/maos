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
	P_NX,
	P_NY,
	P_CX,
	P_CY,
	P_R,
	P_TOT,
    };
    enum{/*output */
	PL_MAP,
	PL_TOT,
    };
    (void)nlhs;
    if(nrhs!=P_TOT){
	mexErrMsgTxt("Usage: [map]=mkcirmap(nx, ny, cx, cy, radius)\n");
    }
    int nx=(int)mxGetScalar(prhs[P_NX]);
    int ny=(int)mxGetScalar(prhs[P_NY]);
    double cx=mxGetScalar(prhs[P_CX])-1;//-1 to convert from matlab convention to C
    double cy=mxGetScalar(prhs[P_CY])-1;
    double R=mxGetScalar(prhs[P_R]);
    dmat *map=dnew(nx, ny);
    dcircle(map, cx, cy, 1, 1,  R, 1);
    plhs[PL_MAP]=d2mx(map);
    dfree(map);
}
