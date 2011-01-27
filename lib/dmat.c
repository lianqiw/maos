/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#undef USE_COMPLEX
#include "mat.c"
#include "cell.c"
#include "matbin.c"

/**
   Compute the enclosed energy or azimuthal average of a.
*/
dmat *denc(dmat *psf, /**<The input array*/
	   dmat *dvec,/**<The diameter for enclosed energy, or radius for azimuthal average*/
	   int type   /**<The type. -1: azimuthal average, 0: within a square, 1: within a circle, 2: within a slit*/
	   ){
    double dmax=dvec->p[dvec->nx-1];
    double ncomp;
    if(type==-1){
	ncomp=nextpow2(dmax*2);
    }else{
	ncomp=nextpow2(dmax);
    }
    double ncomp_max=psf->nx>psf->ny?psf->nx:psf->ny;
    dmat *psfc;
    if(ncomp_max > ncomp){
	psfc=dnew(ncomp, ncomp);
	dembed(psfc, psf, 0);
    }else{
	ncomp=ncomp_max;
	psfc=dref(psf);
    }
    long ncomp2=ncomp*2;
    cmat *psf2=cnew(ncomp2, ncomp2);
    cfft2plan(psf2, -1);
    cembedd(psf2, psfc, 0);
    dfree(psfc);
    cfftshift(psf2);
    cfft2(psf2, -1);
    dmat *psf3=NULL;
    creal2d(&psf3, 0, psf2, 1);
    dscale(psf3, pow(ncomp2,-2));
    cfree(psf2);
    double dk=1./ncomp2;
    double pi2=2*M_PI;
    PDMAT(psf3, ppsf);

    dmat *enc=dnew(dvec->nx, 1);
    double *restrict dr=dvec->p;

    if(type==0){
	dmat *ksinc=dnew(dvec->nx, ncomp2);
	PDMAT(ksinc, pks);
	for(long iy=0; iy<ncomp2; iy++){
	    double ky=(iy<ncomp?iy:iy-ncomp2)*dk;
	    for(long ir=0; ir<dvec->nx; ir++){
		pks[iy][ir]=sinc(ky*dr[ir])*dr[ir];
	    }
	}
	for(long iy=0; iy<ncomp2; iy++){
	    for(long ix=0; ix<ncomp2; ix++){
		for(long ir=0; ir<dvec->nx; ir++){
		    double s=pks[iy][ir]*pks[ix][ir];
		    enc->p[ir]+=s*ppsf[iy][ix];
		}
	    }
	}
    }else{
	for(long iy=0; iy<ncomp2; iy++){
	    double ky=(iy<ncomp?iy:iy-ncomp2)*dk;
	    info("%ld of %ld for %ld\n", iy, ncomp2, dvec->nx);
	    for(long ix=0; ix<ncomp2; ix++){
		double kx=(ix<ncomp?ix:ix-ncomp2)*dk;
		switch(type){
		case -1: {//azimuthal average. dr is radius
		    double k=sqrt(kx*kx+ky*ky);
		    for(long ir=0; ir<dvec->nx; ir++){
			double s=j0(k*pi2*dr[ir]);
			enc->p[ir]+=s*ppsf[iy][ix];
		    }
		} break;
		case 0:
		    break;
		case 1: {//Encircled energy. dr is diameter
		    double k=sqrt(kx*kx+ky*ky);
		    for(long ir=0; ir<dvec->nx; ir++){
			const double r=dr[ir]*0.5;
			const double tmp=k*pi2*r;
			double s=j1(tmp)*r/k;
			if(!ix && !iy) s=pi2*r*r;//special case.
			enc->p[ir]+=s*ppsf[iy][ix];
		    }
		} break;
		case 2://Enstripped energe in a slit.
		    error("To implement: Do FFT only along 1-d\n");
		    break;
		default:
		    error("Not implemented\n");
		}
	    }
	}
    }
    dfree(psf3);
    return enc;
}
