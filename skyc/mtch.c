/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
void psf2i0gxgy(dmat *i0, dmat *gx, dmat *gy, dmat *psf, DTF_S *dtf){
    cmat *otf=cnew(psf->nx, psf->ny);
    cmat *otfsave=cnew(psf->nx, psf->ny);
    cfft2plan(otf, 1);
    cfft2plan(otf, -1);
    cfft2plan(otfsave,1);
    ccpd(&otf, psf);/*loaded psf has peak in corner */
    cfft2i(otf, -1);/*turn to OTF, peak in corner. was 1, turn to -1 on 1/30/2013 */
    ccwm(otf, dtf->nominal);
    ccp(&otfsave, otf);
    cfft2(otf, 1);/*turn back. */
    spmulcreal(i0->p, dtf->si, otf->p, 1);
    ccp(&otf, otfsave);
    PCMAT(otf, potf);
    PCMAT(otfsave, potfsave);
    /*Now derivative */
    for(int iy=0; iy<otf->ny; iy++){
	for(int ix=0; ix<otf->nx; ix++){
	    potf[iy][ix]*=dtf->U->p[ix];
	    potfsave[iy][ix]*=dtf->U->p[iy];
	}
    }
    cfft2(otf, 1);//was 1, changed to -1 on 1/29/2013
    cfft2(otfsave, 1);//was 1, changed to -1 on 1/29/2013
    spmulcreal(gx->p,dtf->si,otf->p,1);
    spmulcreal(gy->p,dtf->si,otfsave->p,1);
    cfree(otf); cfree(otfsave);
}

/**
   shift without wraping i0 into i0x1 (+1) and i0x2 (-1)
*/
static void mki0shx(double *i0x1, double *i0x2, dmat *i0, double scale){
    int nx=i0->nx;
    double (*i0p)[nx]=(void*)i0->p;
    double (*i0x1p)[nx]=(void*)i0x1;
    double (*i0x2p)[nx]=(void*)i0x2;
    for(int iy=0; iy<i0->ny; iy++){
	for(int ix=0; ix<i0->nx-1; ix++){
	    i0x1p[iy][ix+1]=i0p[iy][ix]*scale;
	    i0x2p[iy][ix]=i0p[iy][ix+1]*scale;
	}
    }
}
/**
  shift without wraping i0 into i0y1 (+1) and i0y2 (-1)
*/
static void mki0shy(double *i0y1, double *i0y2, dmat *i0, double scale){
    int nx=i0->nx;
    double (*i0p)[nx]=(void*)i0->p;
    double (*i0y1p)[nx]=(void*)i0y1;
    double (*i0y2p)[nx]=(void*)i0y2;
    for(int iy=0; iy<i0->ny-1; iy++){
	for(int ix=0; ix<i0->nx; ix++){
	    i0y1p[iy+1][ix]=i0p[iy][ix]*scale;
	    i0y2p[iy][ix]=i0p[iy+1][ix]*scale;
	}
    }
}
/**
   Compute matched filter.
 */
void mtch(dcell **mtche, dmat **sanea,
	  dcell *i0, dcell *gx, dcell *gy, double pixtheta, 
	  double rne, double bkgrnd, int cr){
    double kp=1./pixtheta;
    const double psfvar=0;
    int nmod=3;
    if(cr){
	nmod=7;
    }
    const long nsa=i0->nx;
    const long nx=i0->p[0]->nx;
    const long ny=i0->p[0]->ny;
    const long npixtot=nx*ny;
    dmat *i0m=dnew(2,nmod);
    dmat *i0g=dnew(npixtot, nmod);
    PDMAT(i0m,pi0m);
    PDMAT(i0g,pi0g);
    dmat *wt=dnew(npixtot, 1);
    
    *mtche=dcellnew(nsa,1);
    *sanea=dnew(nsa*2,1);
    
    pi0m[0][0]=1;
    pi0m[1][1]=1;

    if(cr){
	pi0m[3][0]=1;
	pi0m[4][0]=-1;
	pi0m[5][1]=1;
	pi0m[6][1]=-1;
    }
    for(long isa=0; isa<nsa; isa++){    
	dzero(i0g);
	/*kp is here to ensure good conditioning */
	adddbl(pi0g[0], 1, gx->p[isa]->p, npixtot, 1, 0);
	adddbl(pi0g[1], 1, gy->p[isa]->p, npixtot, 1, 0);
	adddbl(pi0g[2], 1, i0->p[isa]->p, npixtot, kp, 0);
	if(cr){
	    mki0shx(pi0g[3], pi0g[4], i0->p[isa], kp);
	    mki0shy(pi0g[5], pi0g[6], i0->p[isa], kp);
	}
	for(long ipix=0; ipix<npixtot; ipix++){
	    wt->p[ipix]=1./(rne*rne+bkgrnd+i0->p[isa]->p[ipix]+psfvar*i0->p[isa]->p[ipix]);
	}
	
	dmat *tmp=dpinv(i0g, wt,NULL);
	dmm(&(*mtche)->p[isa],i0m,tmp,"nn",1);
	dfree(tmp);
	dcwpow(wt,-1);
	dmat *nea2=dtmcc((*mtche)->p[isa], wt);
	/*Drop coupling in x/y gradients. */
	(*sanea)->p[isa]=sqrt(nea2->p[0]);
	(*sanea)->p[isa+nsa]=sqrt(nea2->p[3]);
	dfree(nea2);
    }
    dfree(wt);
    dfree(i0g);
    dfree(i0m);
}
