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
#include "random.h"
#include "interface.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    rand_t *strand;
    double *p;
    int nelem=1,M;
    int i;
    if(nrhs<1){ 
	mexErrMsgTxt("Usage: val=crandn(stat)\n");
    }
    strand=(rand_t*)mxGetPr(prhs[0]);
    if(nrhs==1){
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    }else if(nrhs>=2){
	M=mxGetNumberOfElements(prhs[1]);
	p=mxGetPr(prhs[1]);
	if(M==2){
	    plhs[0]=mxCreateDoubleMatrix((int)p[0],(int)p[1],mxREAL);
	}else if(M==1){
	    if(nrhs==2){
		plhs[0]=mxCreateDoubleMatrix((int)p[0],1,mxREAL);
	    }else if(nrhs==3){
		double *p2;
		p2=mxGetPr(prhs[2]);
		plhs[0]=mxCreateDoubleMatrix((int)p[0],(int)p2[0],mxREAL);
	    }else{
		mexErrMsgTxt("Unexpected input: accepts at most 3 inputs");
	    }
	}else{
	    mexErrMsgTxt("Unexpected input");
	}
    }
    nelem=mxGetNumberOfElements(plhs[0]);
    p=mxGetPr(plhs[0]);
    for(i=0; i<nelem; i++){
	p[i]=randn(strand);
    }
}
