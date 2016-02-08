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
/*
  Create a random stream with input seed
*/
#include "interface.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    unsigned int seed,nlem;
    rand_t *p;
    if(nrhs!=1){
	mexErrMsgTxt("Usage: stat=crandcreate(seed)");
    }
    seed=(unsigned int)mxGetScalar(prhs[0]);
    nlem=(unsigned int)ceil((double)sizeof(rand_t)
			    /(double)sizeof(mxINT32_CLASS));
    plhs[0]=mxCreateNumericMatrix(nlem,1,mxINT32_CLASS,mxREAL);
    p=(rand_t*)mxGetPr(plhs[0]);
    fprintf(stderr,"Seed=%u\n",seed);
    seed_rand(p,seed);
}
