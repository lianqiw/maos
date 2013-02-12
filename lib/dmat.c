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

#undef USE_COMPLEX
#include "mat.c"
#include "cell.c"
#include "matbin.c"
/*
  This file contains routines that are only useful for dmat (not cmat).
*/

typedef struct{
    dmat *enc; /**<Output*/
    dmat *dvec;/**<Radius wanted*/
    dmat *phat; /**<processed image.*/
    int type;  
}ENC_T;
void denc_thread(thread_t *pdata){
    ENC_T *data=pdata->data;
    const dmat *dvec=data->dvec;
    dmat *enc=data->enc;
    PDMAT(data->phat, ppsf);
    int type=data->type;
    const double *restrict dr=dvec->p;
    const long ncomp2=data->phat->nx;
    const long ncomp=ncomp2/2;
    const double dk=1./ncomp2;
    const double pi2=2*M_PI;
    if(type==0){
	dmat *ksinc=dnew(dvec->nx, ncomp2);
	PDMAT(ksinc, pks);
	/*Cache the data. */
	for(long iy=0; iy<ncomp2; iy++){
	    double ky=(iy<ncomp?iy:iy-ncomp2)*dk;
	    for(long ir=pdata->start; ir<pdata->end; ir++){
		pks[iy][ir]=sinc(ky*dr[ir])*dr[ir];
	    }
	}
	for(long iy=0; iy<ncomp2; iy++){
	    for(long ix=0; ix<ncomp2; ix++){
		for(long ir=pdata->start; ir<pdata->end; ir++){
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
		case -1: {/*azimuthal average. dr is radius */
		    double k=sqrt(kx*kx+ky*ky);
		    for(long ir=pdata->start; ir<pdata->end; ir++){
			double s=j0(k*pi2*dr[ir]);
			enc->p[ir]+=s*ppsf[iy][ix];
		    }
		} break;
		case 0:
		    break;
		case 1: {/*Encircled energy. dr is diameter */
		    double k=sqrt(kx*kx+ky*ky);
		    for(long ir=pdata->start; ir<pdata->end; ir++){
			const double r=dr[ir]*0.5;
			const double tmp=k*pi2*r;
			double s=j1(tmp)*r/k;
			if(!ix && !iy) s=pi2*r*r;/*special case. */
			enc->p[ir]+=s*ppsf[iy][ix];
		    }
		} break;
		case 2:/*Enstripped energe in a slit. */
		    error("To implement: Do FFT only along 1-d\n");
		    break;
		default:
		    error("Not implemented\n");
		}
	    }
	}
    }
}
/**
   Compute the enclosed energy or azimuthal average of a.
*/
dmat *denc(dmat *psf, /**<The input array*/
	   dmat *dvec,/**<The diameter for enclosed energy, or radius for azimuthal average*/
	   int type,  /**<The type. -1: azimuthal average, 0: within a square, 1: within a circle, 2: within a slit*/
	   int nthread
	   ){
    double rmax=dvec->p[dvec->nx-1];
    double ncomp;
    if(type==-1){
	ncomp=nextfftsize(rmax*2);
    }else{
	ncomp=nextfftsize(rmax);
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
    dmat *phat=NULL;
    creal2d(&phat, 0, psf2, 1);
    dscale(phat, pow(ncomp2,-2));
    cfree(psf2);
    dmat *enc=dnew(dvec->nx, 1);    
    ENC_T data={enc, dvec, phat, type};
    thread_t info[nthread];
    thread_prep(info, 0, dvec->nx, nthread, denc_thread, &data);
    CALL_THREAD(info, nthread, 0);
    dfree(phat);
    return enc;
}
