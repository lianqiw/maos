/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "skyc.h"
#include "mtch.h"
/**
   \file skyc/mtch.c
   Routings regarding matched filter computation.
*/
/**
   Compute pixel intensity i0, and gradient of pixel intensity gx, gy from the
   PSF.
   \todo merge this file with maos/mtch.c and put in lib folder.
*/
void psf2i0gxgy(dmat *i0, dmat *gx, dmat *gy, dmat *psf, DTF_S *dtf, int deriv){
    cmat *otf=cnew(psf->nx, psf->ny);
    cmat *otfsave=0;
    ccpd(&otf, psf);//loaded psf has peak in corner 
    cfft2i(otf, -1);//turn to OTF, peak in corner.
    ccwm(otf, dtf->nominal);//apply pixel transfer function
    if(deriv){
	ccp(&otfsave, otf);//save for later
    }
    cfft2(otf, 1);//convert otf back to psf space
    dspmulcreal(i0->p, dtf->si, otf->p, 1);//sample psf to detectors.
    if(deriv){
	ccp(&otf, otfsave); //copy back otf
	//apply derivative.
	for(int iy=0; iy<otf->ny; iy++){
	    for(int ix=0; ix<otf->nx; ix++){
		P(otf,ix,iy)*=dtf->U->p[ix];
		P(otfsave,ix,iy)*=dtf->U->p[iy];
	    }
	}
	cfft2(otf, 1);
	cfft2(otfsave, 1);
	dspmulcreal(gx->p,dtf->si,otf->p,1);
	dspmulcreal(gy->p,dtf->si,otfsave->p,1);
	cfree(otfsave);
    }
    cfree(otf); 
}

/**
   Compute matched filter.
*/
void genmtch(dcell **mtche, dmat **sanea,
	     dcell *i0, dcell *gx, dcell *gy, real pixtheta, 
	     real rne, real bkgrnd, int cr){
    const long nsa=i0->nx;
    if(!*mtche){
	*mtche=dcellnew(nsa,1);
    }
    if(!*sanea){
	*sanea=dnew(nsa*2,1);
    }
    for(long isa=0; isa<nsa; isa++){    
	dmat *nea2=0;
	mtch(PP(*mtche, isa),&nea2, P(i0, isa), P(gx, isa), P(gy, isa), 0, 0, 0,
	     bkgrnd, bkgrnd, rne, pixtheta, pixtheta, 0, 0, cr);
	/*Drop coupling in x/y gradients. */
	(*sanea)->p[isa]=sqrt(nea2->p[0]);
	(*sanea)->p[isa+nsa]=sqrt(nea2->p[3]);
	dfree(nea2);
    }
}
