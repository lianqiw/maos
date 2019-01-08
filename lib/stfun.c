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
#include "../math/mathdef.h"
#include "stfun.h"

/**
   Data struct for structure function computation from digital data.
*/
struct stfun_t{
    long count;
    dmat *amp;
    cmat *hat0;
    cmat *hat1;
    cmat *hat2;
    cmat *hattot;
};

/**
   Initialize the stfun data.
*/
stfun_t *stfun_init(long nx, long ny, double *amp){
    stfun_t *A=mycalloc(1,struct stfun_t);
    A->count=0;
    A->hat0=cnew(nx*2, ny*2);
    A->hat1=cnew(nx*2, ny*2);
    A->hat2=cnew(nx*2, ny*2);
    A->hattot=cnew(nx*2, ny*2);

    //cfft2plan(A->hat0, -1);
    //cfft2plan(A->hat0, 1);
    //cfft2plan(A->hat1, -1);
    //cfft2plan(A->hat2, -1);
    //cfft2plan(A->hattot, 1);

    dmat *damp;
    if(amp){
	damp=dnew_ref(nx,ny,amp);
    }else{
	damp=dnew(nx, ny);
	dset(damp, 1);
    }
    cembedd(A->hat0,damp,0);
    A->amp=damp;
    cfft2(A->hat0, -1);
    return A;
}
void stfun_push(stfun_t *A, dmat *opd){
    A->count++;
    long ny=A->hat0->ny/2;/*maybe smaller than opd. */
    long nx=A->hat0->nx/2;
    long nx2=nx>>1;
    long ny2=ny>>1;
    dmat*  popd=opd;
    dmat* pamp=A->amp;
    cmat*  p1=A->hat1;
    cmat*  p2=A->hat2;
    cset(A->hat1, 0);
    cset(A->hat2, 0);
    for(long iy=0; iy<ny; iy++){
	for(long ix=0; ix<nx; ix++){
	    double o=IND(popd,ix,iy)*IND(pamp,ix,iy);
	    IND(p1,ix+nx2,iy+ny2)=o;
	    IND(p2,ix+nx2,iy+ny2)=o*o;
	}
    }
    cfft2(A->hat1, -1);
    cfft2(A->hat2, -1);
    for(long i=0; i<A->hat1->nx*A->hat1->ny; i++){
	/*real(t2*t0*)-t1*t1* */
	A->hattot->p[i]+=creal(A->hat2->p[i]*conj(A->hat0->p[i]))
	    -A->hat1->p[i]*conj(A->hat1->p[i]);
    }
}
dmat *stfun_finalize(stfun_t *A){
    cscale(A->hattot, 2./A->count);
    cfft2i(A->hattot, 1);
    cabs2toreal(A->hat0);
    cfft2i(A->hat0, 1);
    cfftshift(A->hattot);
    cfftshift(A->hat0);
    long nx=A->hat0->nx;
    long ny=A->hat0->ny;
    dmat *st=dnew(nx, ny);
    dmat* pst=st;
    cmat*  p1=A->hattot;
    cmat*  p2=A->hat0;
    for(long iy=1; iy<ny; iy++){/*skip first row/column where hat0 is 0. */
	for(long ix=1; ix<nx; ix++){
	    IND(pst,ix,iy)=creal(IND(p1,ix,iy)/IND(p2,ix,iy));
	}
    }
    cfree(A->hat0);
    cfree(A->hat1);
    cfree(A->hat2);
    cfree(A->hattot);
    dfree(A->amp);
    free(A);
    return st;
}
/**
   Generate the structure function of the phase of kolmogorov spectrum 
*/
dmat* stfun_kolmogorov(loc_t *loc, double r0){
    long nloc=loc->nloc;
    double *locx=loc->locx;
    double *locy=loc->locy;
    dmat *st=dnew(nloc,nloc);
    dmat* B=st;
    double coeff=6.88*pow(r0,-5./3.)*pow(0.5e-6/(2.*M_PI),2);
    for(int i=0; i<nloc; i++){
	for(int j=i; j<nloc; j++){
	    double rdiff2=pow(locx[i]-locx[j],2)+pow(locy[i]-locy[j],2);
	    IND(B,i,j)=IND(B,j,i)=coeff*pow(rdiff2,5./6.);
	}
    }
    return st;
}
