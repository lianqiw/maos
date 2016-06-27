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
	P_PSF,
	P_DIAM,
	P_TYPE,
	P_TOT,
    };
    enum{
	PL_ENC,
	PL_TOT,
    };
    (void)nlhs;
    if(nrhs !=P_TOT){
	mexErrMsgTxt("Usage: enc=calcenc(psf, diam, type)"
		     "diam is The diameter for enclosed energy, or radius for azimuthal average"
		     "Type can be: "
		     "-1: azimuthal average, "
		     "0: within a square, "
		     "1: within a circle, "
		     "2: within a slit");
    }
    dmat *psf=mx2d(prhs[P_PSF]);
    dmat *diam=mx2d(prhs[P_DIAM]);
    if(diam->ny>1){
	if(diam->nx==1){
	    diam->nx=diam->ny;
	    diam->ny=1;
	}else{
	    error("diam has wrong size: %ldx%ld\n", diam->nx, diam->ny);
	}
    }
    int type=(int)mxGetScalar(prhs[P_TYPE]);
    static int inited=0;
    if(!inited){
	THREAD_POOL_INIT(NCPU); 
	inited=1;
    }
    dmat *enc=denc(psf, diam, type, NCPU);
    plhs[PL_ENC]=d2mx(enc);
    dfree(psf);
    dfree(diam);
    dfree(enc);
}
