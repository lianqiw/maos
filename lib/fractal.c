#include "fractal.h"
#include "nr/nr.h"
#include "sys/thread.h"
/**
   \file fractal.c

   Implementation of the fractal operation for atmospheric turbulence screen
   generation and reconstruction.

*/

/*
  Initialize some constants.
*/
static double vkcoeff;
static double vkcoeff0;
static void __attribute__((constructor)) init(void){
    vkcoeff=tgamma(11./6)/(pow(2, 5./6.) * pow(M_PI, 8./3.)) * pow(24./5*tgamma(6./5.), 5./6.)
	*pow(2*M_PI/0.5e-6, -2);
    //for cov(0)
    vkcoeff0=tgamma(11./6)*tgamma(5./6.)/ ( 2* pow(M_PI, 8./3.)) *pow(24./5*tgamma(6./5.), 5./6.)
	*pow(2*M_PI/0.5e-6, -2);
}
/**
   Compute Von Karman covariance function at separations computed from the grid
   size nx and sampling dx, with Fried parameter of r0, and outerscale of L0.  

   Return matrix:
   The first row is cov(i*dx) where i=0, 1, 2, 3, 4 ... n
   The second row is cov(i*sqrt(2)*dx) where i=0, 1, 2, 3, 4 ... n
   
*/
typedef struct vkcov_t{
    double r0;
    double L0;
    double dx;
    long n;
    dmat *cov;
    struct vkcov_t *next;
}vkcov_t;
vkcov_t *head=NULL;
PNEW(mutex_cov);
static dmat *vkcov_get(double r0, double L0, double dx, long n){
    for(vkcov_t *p=head; p; p=p->next){
	if(fabs(p->r0-r0)<EPS && (fabs(p->L0-L0)<EPS || (isinf(p->L0) && isinf(L0)))
	   && fabs(p->dx-dx)<EPS && p->n == n){
	    //info("found saved vkcov with r0=%g, L0=%g, dx=%g, n=%ld\n", r0, L0, dx, n);
	    return p->cov;
	}
    }
    info("not found saved vkcov with r0=%g, L0=%g, dx=%g, n=%ld\n",
	 r0, L0, dx, n);
    return NULL;
}
static dmat* vkcov_calc(double r0, double L0, double dx, long n){
    if(L0>9000) L0=INFINITY;//L0 bigger than 9000 is treated as infinity.
    dmat *cov=vkcov_get(r0, L0, dx, n);
    if(cov) return cov;
    vkcov_t *node=calloc(1, sizeof(vkcov_t));
    node->r0=r0;
    node->L0=L0;
    node->dx=dx;
    node->n=n;
    long nroot=(long)round(log2((double)n-1));
    node->cov=cov=dnew(2, nroot+2);
    node->next=head;
    if(r0>=L0){
	error("Illegal parameter: r0=%g, L0=%g\n", r0, L0);
    }
    head=node;
    PDMAT(cov, pcov);
    const double sqrt2=sqrt(2);
    if(isinf(L0)){//kolmogorov, compute from structure function
	const double power=5./3.;
	double coeff=6.88*pow(2*M_PI/0.5e-6, -2) * pow(r0, -power);
	double D=(n-1)*dx;
	double sigma2=0.5*coeff*pow(sqrt2*D, power);
	pcov[0][0]=pcov[0][1]=sigma2;
	for(long i=0; i<=nroot; i++){
	    long j=1<<i;
	    pcov[i+1][0]=sigma2-0.5*coeff*pow(j*dx, power);
	    pcov[i+1][1]=sigma2-0.5*coeff*pow(j*dx*sqrt2, power);
	}
    }else{//compute from Eq 12,16 in Roldolphe Conan's 2008 Paper.
	const double f0=1./L0;
	const double r0f0p=pow(r0*f0, -5./3.);
	double ri, rk, rip, rkp;
	double r2pif0;
	pcov[0][0]=pcov[0][1]=vkcoeff0*r0f0p;
	for(long i=0; i<=nroot; i++){
	    long j=1<<i;
	    r2pif0=(j*dx)*2*M_PI*f0;
	    bessik(r2pif0, 5./6., &ri, &rk, &rip, &rkp);
	    pcov[i+1][0]=vkcoeff*r0f0p*pow(r2pif0, 5./6.)*rk;
	
	    r2pif0=(j*dx*sqrt2)*2*M_PI*f0;
	    bessik(r2pif0, 5./6., &ri, &rk, &rip, &rkp);
	    pcov[i+1][1]=vkcoeff*r0f0p*pow(r2pif0, 5./6.)*rk;
	}
	if(nroot>0 && pcov[1][0]>pcov[0][0]){//when L0 is too big, some round of error happens
	    pcov[0][0]=pcov[1][0];
	}
    }
    return cov;
}

/**
   Apply forward fractal operator recursively.

   To confirm that fractal_trans is right, compare u'Av with v'A'u, where ' is transpose. 
   u,v are two vectors (reshape to 2-d arrays when applying A, or A')
*/
#define FRACTAL fractal
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

