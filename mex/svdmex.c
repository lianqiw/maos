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
	P_A,
	P_TOT,
    };
    enum{
	PL_U,
	PL_S,
	PL_V,
	PL_TOT,
    };
  if(nlhs!=PL_TOT || nrhs !=P_TOT){
      mexErrMsgTxt("Usage: [u,s,v]=svdmex(A)\n");
  }
  dmat *A=mx2d(prhs[P_A]);
  dmat *U=NULL, *Sdiag=NULL, *VT=NULL;
  dsvd(&U, &Sdiag, &VT, A);
  plhs[PL_U]=d2mx(U);
  plhs[PL_S]=d2mx(Sdiag);
  dmat *V=dtrans(VT);
  plhs[PL_V]=d2mx(V);
  dfree(U); dfree(Sdiag); dfree(VT); dfree(V); dfree(A); 
}
