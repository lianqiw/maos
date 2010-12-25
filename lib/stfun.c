#include "type.h"
#include "stfun.h"
#include "dmat.h"
#include "cmat.h"
#include "cmat_extra.h"
#include "fft.h"

/**
   Data struct for structure function computation.
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
void stfun_init(stfun_t *A, long nx, long ny, double *amp){
    A->count=0;
    A->hat0=cnew(nx*2, ny*2);
    A->hat1=cnew(nx*2, ny*2);
    A->hat2=cnew(nx*2, ny*2);
    A->hattot=cnew(nx*2, ny*2);

    cfft2plan(A->hat0, -1);
    cfft2plan(A->hat0, 1);
    cfft2plan(A->hat1, -1);
    cfft2plan(A->hat2, -1);
    cfft2plan(A->hattot, 1);

    dmat *damp;
    if(amp){
	damp=dnew_data(amp,nx,ny);
    }else{
	damp=dnew(nx, ny);
	dset(damp, 1);
    }
    cembedd(A->hat0,damp,0);
    A->amp=damp;
    cfft2(A->hat0, -1);
}
void stfun_push(stfun_t *A, dmat *opd){
    A->count++;
    long ny=A->hat0->ny/2;//maybe smaller than opd.
    long nx=A->hat0->nx/2;
    long nx2=nx>>1;
    long ny2=ny>>1;
    PDMAT(opd, popd);
    PDMAT(A->amp,pamp);
    PCMAT(A->hat1, p1);
    PCMAT(A->hat2, p2);
    cset(A->hat1, 0);
    cset(A->hat2, 0);
    for(long iy=0; iy<ny; iy++){
	for(long ix=0; ix<nx; ix++){
	    double o=popd[iy][ix]*pamp[iy][ix];
	    p1[iy+ny2][ix+nx2]=o;
	    p2[iy+ny2][ix+nx2]=o*o;
	}
    }
    cfft2(A->hat1, -1);
    cfft2(A->hat2, -1);
    for(long i=0; i<A->hat1->nx*A->hat1->ny; i++){
	//real(t2*t0*)-t1*t1*
	A->hattot->p[i]+=creal(A->hat2->p[i]*conj(A->hat0->p[i]))
	    -A->hat1->p[i]*conj(A->hat1->p[i]);
    }
}
dmat *stfun_finalize(stfun_t *A){
    cscale(A->hattot, 2./A->count);
    cifft2(A->hattot, 1);
    cabs2toreal(A->hat0);
    cifft2(A->hat0, 1);
    cfftshift(A->hattot);
    cfftshift(A->hat0);
    long nx=A->hat0->nx;
    long ny=A->hat0->ny;
    dmat *st=dnew(nx, ny);
    PDMAT(st,pst);
    PCMAT(A->hattot, p1);
    PCMAT(A->hat0, p2);
    for(long iy=1; iy<ny; iy++){//skip first row/column where hat0 is 0.
	for(long ix=1; ix<nx; ix++){
	    pst[iy][ix]=creal(p1[iy][ix]/p2[iy][ix]);
	}
    }
    cfree(A->hat0);
    cfree(A->hat1);
    cfree(A->hat2);
    cfree(A->hattot);
    dfree(A->amp);
    return st;
}
