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
#include "fractal.h"
#include "turbulence.h"
/**
   Implementation of the fractal operation for atmospheric turbulence screen
   generation and reconstruction.

   This method does not produce right covariance statistics in interpolation
   because points that are not directly computed from one another. This renders
   the fractal screens to have lower high frequency components than FFT based screens.

*/

typedef struct vkcov_t{
    double r0;
    double L0;
    double dx;
    long n;
    long ninit;
    dmat *cov;
    dmat *C;/*the covariance matrix. */
    dmat *K;/*initial matrix for generating atmosphere. */
    dmat *KI;/*inverse of K. */
    struct vkcov_t *next;
}vkcov_t;
vkcov_t *head=NULL;
PNEW(mutex_cov);
static vkcov_t *vkcov_get(double r0, double L0, double dx, long n, long ninit){
    for(vkcov_t *p=head; p; p=p->next){
	if(fabs(p->r0-r0)<EPS && (fabs(p->L0-L0)<EPS || (!isfinite(p->L0) && !isfinite(L0)))
	   && fabs(p->dx-dx)<EPS && p->n == n && p->ninit==ninit){
	    return p;
	}
    }
    info("compute vkcov with r0=%g, L0=%g, dx=%g, n=%ld, ninit=%ld\n",
	 r0, L0, dx, n, ninit);
    return NULL;
}
void fractal_vkcov_free(){
    for(vkcov_t *p=head; p; p=head){ 
	head=p->next;
	dfree(p->cov);
	dfree(p->C);
	dfree(p->K);
	dfree(p->KI);
	free(p);
    }
}
static __attribute__((constructor)) void init(){
    register_deinit(fractal_vkcov_free, NULL);
}

/**
   Compute Von Karman covariance function at separations computed from the grid
   size nx and sampling dx, with Fried parameter of r0, and outerscale of L0.  
   ninit is the initial side size of the atm array to start with.
*/
static vkcov_t* vkcov_calc(double r0, double L0, double dx, long n, long ninit){
    if(L0>9000) L0=INFINITY;/*L0 bigger than 9000 is treated as infinity. */
    vkcov_t *node=vkcov_get(r0, L0, dx, n, ninit);
    if(node) return node;
    node=mycalloc(1,vkcov_t);
    node->r0=r0;
    node->L0=L0;
    node->dx=dx;
    node->n=n;
    node->ninit=ninit;
    long nroot=(long)round(log2((double)n-1));
    node->next=head;
    if(r0>=L0){
	error("Illegal parameter: r0=%g, L0=%g\n", r0, L0);
    }
    head=node;
    dmat *r=dnew(2, nroot+2);
    double sqrt2=sqrt(2);
    IND(r,0,0)=0;
    IND(r,1,0)=0;
    for(long i=0; i<=nroot; i++){
	long j=1<<i;
	IND(r,0,i+1)=j*dx;
	IND(r,1,i+1)=j*dx*sqrt2;
    }
    double D=(n-1)*dx;
    node->cov=turbcov(r, D*sqrt(2), r0, L0);
    dfree(r);
    dmat *rc0=dnew(ninit*ninit, ninit*ninit);
    dmat*  rc=rc0;
    double dx2=dx*(n-1)/(ninit-1);
    for(long j=0; j<ninit; j++){
	double y=dx2*j;
	for(long i=0; i<ninit; i++){
	    double x=dx2*i;
	    long k=i+j*ninit;
	    for(long j2=0; j2<ninit; j2++){
		double y2=dx2*j2;
		for(long i2=0; i2<ninit; i2++){
		    double x2=dx2*i2;
		    long k2=i2+j2*ninit;
		    IND(rc,k2,k)=sqrt((x-x2)*(x-x2)+(y-y2)*(y-y2));
		}
	    }
	}
    }
    node->C=turbcov(rc0, D*sqrt(2), r0, L0);
    dfree(rc0);
    dmat *u=NULL, *s=NULL, *v=NULL;
    dsvd(&u, &s, &v, node->C);
    dcwpow(s, 1./2.);
    node->K=ddup(u);
    dmuldiag(node->K, s);

    dcwpow(s, -1);
    dmuldiag(u, s);
    node->KI=dtrans(u);
    dfree(u);
    dfree(v);
    dfree(s);
    /*we have: K*K'==C */
    return node;
}

/**
   Apply forward fractal operator recursively.

   To confirm that fractal_trans is right, compare u'Av with v'A'u, where ' is transpose. 
   u,v are two vectors (reshape to 2-d arrays when applying A, or A')
*/
#define FRACTAL fractal_do
#define INVERSE 0
#define TRANSPOSE 0
#include "fractal_do.c"
#undef FRACTAL
#undef INVERSE
#undef TRANSPOSE

#define FRACTAL fractal_inv
#define INVERSE 1
#define TRANSPOSE 0
#include "fractal_do.c"
#undef FRACTAL
#undef INVERSE
#undef TRANSPOSE

#define FRACTAL fractal_trans
#define INVERSE 0
#define TRANSPOSE 1
#include "fractal_do.c"
#undef FRACTAL
#undef INVERSE
#undef TRANSPOSE


#define FRACTAL fractal_inv_trans
#define INVERSE 1
#define TRANSPOSE 1
#include "fractal_do.c"
#undef FRACTAL
#undef INVERSE
#undef TRANSPOSE

