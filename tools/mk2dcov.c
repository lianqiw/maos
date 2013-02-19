/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
   Convert covariance matrices to OTF. Deprecated. Use the mex routine
   genotfmex.c in mex/ folder instead for easy access in MATLAB.
  */
#include "../lib/aos.h"
int main(int argc, char *argv[]){
    enum{
	P_EXE,
	P_RES,
	P_LOC,
	P_AMP,
	P_COV,
	P_NORM,
	P_TOT,
    };
    if(argc<P_TOT){
	info("Usage: %s res.bin loc.bin, amp.bin, cov.bin, normalize\n", argv[0]);
	exit(0);
    }
    
    loc_t *loc=locread("%s", argv[P_LOC]);
    dmat *amp=NULL;
    if(strcmp(argv[P_AMP],"NULL")) {
	amp=dread("%s",argv[P_AMP]);
	normalize_max(amp->p, amp->nx, 1);
    }
    dcell *cov=dcellread("%s",argv[P_COV]);
    dcell *cov2d=dcellnew(cov->nx, cov->ny);
    int norm=(int)strtol(argv[P_NORM], NULL, 10);
    for(int i=0; i<cov->nx*cov->ny; i++){
	mk2dcov(&cov2d->p[i], loc, amp?amp->p:NULL, 0.5, cov->p[i], norm);
    }
    dcellwrite(cov2d, "%s", argv[P_RES]); 
}
