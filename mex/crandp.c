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
#include "interface.h"
#include "random.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    rand_t *strand;
    double *p,*pin;
    int nx,ny;
    int i;
    if(nrhs!=2){
	mexErrMsgTxt("Usage: val=crandp(stat, mean[])\n");
    }
    strand=(rand_t*)mxGetPr(prhs[0]);
    nx=mxGetM(prhs[1]);
    ny=mxGetN(prhs[1]);
    pin=mxGetPr(prhs[1]);
    plhs[0]=mxCreateDoubleMatrix(nx,ny,mxREAL);
    p=mxGetPr(plhs[0]);
    for(i=0; i<nx*ny; i++){
	p[i]=(double)randp(strand,pin[i]);
    }
}
