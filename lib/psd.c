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
#include "../sys/sys.h"
#include "dmat.h"
#include "cmat.h"
#include "fft.h"
/**
   Compute the PSD from a sequence.
*/
//#define W_J(i,N2) (1-pow((double)(i-N2)/(double)N2, 2))
#define W_J(i,N2) 1
dmat *psd1d(dmat *v, /**<[in] The data sequency*/
	   long lseg /**<[in] The length of overlapping segments*/
	   ){
    long lseg2=lseg>>1;
    if(v->nx==1){
	v->nx=v->ny;
	v->ny=1;
    }
    int nseg=(v->nx/lseg2-1)>>1; /*number of segments*/
    if(nseg<1){
	nseg=1;
	lseg=v->nx;
	lseg2=lseg>>1;
    }
    long ncol=v->ny;
    dmat *psd=dnew(lseg2+1, ncol);
    cmat *hat=cnew(lseg, 1);
    cfft2plan(hat, -1);
    for(long icol=0; icol<ncol; icol++){
	double *ppsd=psd->p+icol*(lseg2+1);
	for(int iseg=0; iseg<nseg; iseg++){
	    double* p=v->p+icol*v->nx+iseg*lseg2;
	    for(int ix=0; ix<lseg; ix++){
		hat->p[ix]=p[ix]*W_J(ix, lseg2);
	    }
	    cfft2(hat, -1);
	    ppsd[0]+=cabs2(hat->p[0]);
	    for(int ix=1; ix<lseg2; ix++){
		ppsd[ix]+=cabs2(hat->p[ix])+cabs2(hat->p[lseg-ix]);
	    }
	    ppsd[lseg2]+=cabs2(hat->p[lseg2]);
	}
    }
    double sumwt=0;
    for(int ix=0; ix<lseg; ix++){
	sumwt+=pow(W_J(ix, lseg2), 2);
    }
    sumwt*=lseg*nseg;
    dscale(psd, 1./sumwt);
    cfree(hat);
    return psd;
}
/**
  Wrap of psd1d to put the frequency along the first column.
*/
dmat *psd1dt(dmat *v, long lseg, double dt){
    dmat *psd=psd1d(v, lseg);
    dmat *psd2=dnew(psd->nx, psd->ny+1);
    int N=(psd->nx-1)*2;
    double df=1./(N*dt);
    for(int i=0; i<psd->nx; i++){
	psd2->p[i]=df*i;
    }
    dscale(psd, 1./df);//divide so the value is point, not integrated in a bin.
    memcpy(psd2->p+psd->nx, psd->p, psd->nx*psd->ny*sizeof(double));
    dfree(psd);
    return psd2;
}

/*Interpolate psd onto new f. We interpolate in log space which is more linear.*/
dmat *psdinterp1(const dmat *psdin, const dmat *fnew){
    dmat *f1=dsub(psdin, 0, 0, 0, 1);
    dmat *psd1=dsub(psdin, 0, 0, 1, 1);
    dmat *f2=ddup(fnew);
    dcwlog(psd1);
    dmat *psd2=dinterp1(f1, psd1, f2);
    dfree(f1); dfree(f2); dfree(psd1);
    dcwexp(psd2,1);
    return psd2;
}
