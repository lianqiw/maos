#include "../lib/aos.h"
/*
  Test mrq min on psf.
 */
typedef struct info{
    cmat *otf;
    dmat *psf;
    dmat *psf2;//output;
    double rne2;
}info;
double funk(double *x, info *data){
    cmat *otf2=cnew(data->otf->nx, data->otf->ny);
    cfft2plan(otf2,1);
    ccp(&otf2, data->otf);
    ctilt(otf2, x[1], x[2], 0);
    cfft2(otf2, 1);
    creal2d(&data->psf2, 0, otf2, x[0]/(otf2->nx*otf2->ny));
    dfftshift(data->psf2);
    double sigma2=0;
    for(int i=0; i<otf2->nx*otf2->ny; i++){
	sigma2+=pow(data->psf2->p[i]-data->psf->p[i],2)/(data->psf2->p[i]+data->rne2);
    }
    //info("sigma2=%g. x=%g %g %g\n", sigma2, x[0], x[1], x[2]);
    return sigma2;
}

int main(){
    dcell *psfc=dcellread("pistat_seed1_sa2_x20_y40.bin.gz");
    dmat *psf=psfc->p[0];
    dscale(psf, 2000);
    cmat *otf;
    otf=cnew(psf->nx, psf->ny);
    cfft2plan(otf, -1);
    cfft2plan(otf, 1);
    cmat *otf2=NULL;
    otf2=cnew(psf->nx, psf->ny);
    cfft2plan(otf2, 1);

    cembedd(otf, psf, 0);
    cfft2(otf, -1);
   
    rand_t rstat;
    double rne=3;
    seed_rand(&rstat, 1);
    
    double xc=0;
    double yc=0.02;
    double scale=1.2;
    int ntot=100;
    dcell *psf2=dcellnew(ntot,1);
    for(int i=0; i<ntot; i++){
	ccp(&otf2, otf);
	ctilt(otf2, xc*i, yc*i, 0);
	cfft2(otf2, 1);
	creal2d(&psf2->p[i], 0, otf2, scale/(psf->nx*psf->ny));
	dfftshift(psf2->p[i]);
	
	for(int ix=0; ix<psf2->nx*psf2->ny; ix++){
	    psf2->p[i]->p[ix]=scale*randp(&rstat, psf2->p[i]->p[ix]) + randn(&rstat)*rne;
	}
    }
    dfftshift(psf);
    
 

    //  ddraw("psf", psf, NULL, NULL, "psf", "x", "y", "psf");
    //ddraw("psf", psf2, NULL, NULL, "psf2", "x", "y", "psf2");
    
    dmat *x=dnew(psf->nx, psf->ny);
    for(int i=0; i<x->nx*x->ny; i++){
	x->p[i]=i;
    }
    int nmod=3;
    double ftol=1e-9;

    //double pinit[4][3]={{1, 1.2, 1.8}, {0.8, 1.2, 2.5}, {1.2, 1.0, 1.0}, {1.3, 1.4, 1.0}};
    double pinit[4][3];
    double pinit0[4][3]={{0.8,1.5,2},{1.15,1.4,2},{1.22,1.6,2},{1.21,1.5,2.1}};
    dmat *yinit=dnew(nmod+1,1);
    info arg={otf,NULL,NULL,rne*rne};

    int nfunk;
    double *inits[nmod+1];
    for(int imod=0; imod<nmod+1; imod++){
	inits[imod]=pinit[imod];
    }
    dmat *out=dnew(ntot,3);
    int niter[ntot];
    for(int i=0; i<ntot; i++){
	arg.psf=psf2->p[i];
	memcpy(pinit, pinit0, sizeof(double)*12);
	for(int ic=0; ic<nmod+1; ic++){
	    yinit->p[ic]=funk(pinit[ic], &arg);
	}
	amoeba(inits, yinit->p, nmod, ftol, (double(*)(double *,void*))funk, &arg, &nfunk);
	info("%g %g %g\n", inits[0][0], inits[0][1], inits[0][2]);
	out->p[i+0]=inits[0][0];
	out->p[i+ntot]=inits[0][1];
	out->p[i+ntot*2]=inits[0][2];
	niter[i]=nfunk;
    }
    dwrite(out, "out");
    writeint(niter, ntot, 1, "niter");
    //funk(inits[0], &arg);
    //ddraw("psf", arg.psf2, NULL, NULL, "psf2", "x", "y", "psf2");
    //ddraw("psf", psf2, NULL, NULL, "psf x", "x", "y", "psf input");
}
