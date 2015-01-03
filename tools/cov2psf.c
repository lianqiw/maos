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

/**
   \file cov2psf.c
   Convert covariance matrices to OTF. Deprecated. Use the mex routine
   genotfmex.c in mex/ folder instead for easy access in MATLAB.
  */
#include "../lib/aos.h"
int main(int argc, char *argv[]){
    if(argc<7){
	info("Usage: %s loc.bin amp.bin cov.bin ncomp pttr wvl1 wvl2 ...\n", argv[0]);
	exit(0);
    }
    enum{
	P_EXE,
	P_LOC,
	P_AMP,
	P_COV,
	P_NCOMP,
	P_PTTR,
	P_WVL
    };
    loc_t *loc=locread("%s",argv[P_LOC]);
    dmat *amp=NULL;
    if(strcmp(argv[P_AMP],"NULL")) {
	amp=dread("%s",argv[P_AMP]);
    }else{
	amp=dnew(loc->nloc,1);
	dset(amp,1);
    }
    dmat *cov=dread("%s",argv[P_COV]);
    long ncomp=strtol(argv[P_NCOMP], NULL, 10);
    long pttr=strtol(argv[P_PTTR], NULL, 10);
    cmat *otf=otf=cnew(ncomp, ncomp);
    dmat *psf=NULL;
    for(int iwvl=0; iwvl<argc-P_WVL; iwvl++){
	double wvl=strtod(argv[iwvl+P_WVL], NULL);
	double dtheta=wvl/(ncomp*loc->dx);
	genotf(&otf, loc, amp, NULL, NULL, 0, wvl, dtheta, cov, 0, 0, ncomp, ncomp, 1, pttr);
	writebin(otf, "%s_otf_%g.bin", argv[P_COV], wvl);
	cfftshift(otf);
	cfft2i(otf, 1);
	cfftshift(otf);
	creal2d(&psf, 0, otf, 1);
	writebin(psf, "%s_psf_%g.bin", argv[P_COV], wvl);
    }
    cfree(otf);
    dfree(psf);
    dfree(cov);
    dfree(amp);
    locfree(loc);
}
