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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "interface.h"
#ifdef MATLAB_MEX_FILE
#include <mex.h>
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    enum{
	P_SURF,
	P_X,
	P_Y,
	P_ALX,/*M3 tilt X . pi/2 is no tilt*/
	P_ALY,/*M3 tilt Y */
	P_THETAX,/*guide star offset*/
	P_THETAY,
	P_LOC,
	P_AMP,
	P_TOT,
    };
    enum{
	PL_OPD,
	PL_TOT,
    };
    if(P_TOT!=nrhs){
	mexErrMsgTxt("Usage: OPD=proj(surf, x, y, alx, aly, thetax, thetay, loc, amp)\n");
    }
    rmap_t *mapin=calloc(1, sizeof(rmap_t));
    mapin->p=mxGetPr(prhs[P_SURF]);
    mapin->nx=mxGetM(prhs[P_SURF]);
    mapin->ny=mxGetN(prhs[P_SURF]);
    double alx=mxGetScalar(prhs[P_ALX]);
    double aly=mxGetScalar(prhs[P_ALY]);
    double *X=mxGetPr(prhs[P_X]);
    double *Y=mxGetPr(prhs[P_Y]);
    mapin->dx=X[1]-X[0];
    mapin->dy=Y[1]-Y[0];
    mapin->ox=X[0];
    mapin->oy=Y[0];
    double bx=mxGetScalar(prhs[P_THETAX]);
    double by=mxGetScalar(prhs[P_THETAY]);;
    if (X[mapin->nx-1]-X[0]-(mapin->nx-1)*mapin->dx > 1.e-10)
	mexErrMsgTxt("X has non even spacing\n");
    if (Y[mapin->ny-1]-Y[0]-(mapin->ny-1)*mapin->dy > 1.e-10)
	mexErrMsgTxt("Y has non even spacing\n");
    
    double d_m3_f=20.;/*from m3 to focus*/
    double d_exitpupil_f=46.38661051;
    double d_exitpupil_m3=d_exitpupil_f-d_m3_f;
    double r_exitpupil=1.546220350;
    double r_pupil=15;
    loc_t* loc2=calloc(1, sizeof(loc_t));
    loc2->locx=mxGetPr(prhs[P_LOC]);
    loc2->nloc=mxGetM(prhs[P_LOC]);
    loc2->locy=loc2->locx+loc2->nloc;
    loc2->dx=loc2->locx[1]-loc2->locx[0];
    
    double *amp=mxGetPr(prhs[P_AMP]);
    plhs[PL_OPD]=mxCreateDoubleMatrix(loc2->nloc,1,mxREAL);
    double *phi2=mxGetPr(plhs[PL_OPD]);
    proj_rect_grid(mapin,alx,aly,
		   loc2,-r_exitpupil/r_pupil,r_exitpupil/r_pupil,
		   amp,phi2,-2,d_exitpupil_f,d_exitpupil_m3,bx,by);
    free(loc2);
}
#endif
