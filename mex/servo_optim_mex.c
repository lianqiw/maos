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
	P_PSD,
	P_DT,
	P_DTRAT,
	P_SIGMAN,
	P_STYPE,
	P_TOT,
    };
    enum{
	PL_GAIN,
	PL_EST,
	PL_TOT,
    };

    if(nrhs!=P_TOT){
	mexErrMsgTxt("Usage: [gain este]=servo_optim(psd, dt, dtrat, sigman, servotype);\n"
		     "PSD psd should be in m^2/hz\n"
		     "dt is the AO fundemental sampling period.\n"
		     "dtrat is ratio of sampling period of the WFS over dt.\n"
		     "sigman is wavefront variance due to noise in m^2.\n"
		     "servotype is 1 or 2 for typeI, typeII controller.\n"
		     );
    }
    dmat *psd  = mx2d(prhs[P_PSD]);
    double dt  = mxGetScalar(prhs[P_DT]);
    long dtrat = (long)mxGetScalar(prhs[P_DTRAT]);
    double sigma = mxGetScalar(prhs[P_SIGMAN]);/*m^2 */
    dmat *sigma2 = dnew(1,1); 
    sigma2->p[0] = sigma;
    int servotype = (int)mxGetScalar(prhs[P_STYPE]);
    dcell *res   = servo_optim(psd,dt,dtrat,M_PI/4,sigma2,servotype);
    int ng=res->p[0]->nx-2;
    dmat *gain=dnew(ng, 1);
    memcpy(gain->p, res->p[0]->p, ng*sizeof(double));
    dmat *est=dnew(1,1);
    est->p[0]=res->p[0]->p[ng]+res->p[0]->p[ng+1];
    plhs[PL_GAIN]  = d2mx(gain);
    if(nlhs>PL_EST){
	plhs[PL_EST] = d2mx(est);
    }
    dfree(sigma2);
    dcellfree(res);
    dfree(gain);
    dfree(est);
    dfree(psd);
}
