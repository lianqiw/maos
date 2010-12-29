#include "../lib/aos.h"
static void test_accuracy(){
    rand_t rstat;
    int seed=4;
    double r0=0.2;
    double dx=1./64;
    long N=1+64;
    long nx=N;
    long ny=N;
    seed_rand(&rstat, seed);
    map_t *atm=mapnew(nx, ny, dx, NULL);
    for(long i=0; i<nx*ny; i++){
	atm->p[i]=randn(&rstat);
    }
    map_t *atm2=mapnew(nx, ny, dx, NULL);
    for(long i=0; i<nx*ny; i++){
	atm2->p[i]=randn(&rstat);
    }
    sqmapwrite(atm, "atm_rand.bin");
    sqmapwrite(atm2, "atm2_rand.bin");
    fractal(atm->p, nx, ny, dx, r0);
    sqmapwrite(atm, "atm_frac.bin");
    fractal_inv(atm->p, nx, ny, dx, r0);
    sqmapwrite(atm, "atm_frac_inv.bin");

   
    fractal_trans(atm2->p, nx, ny, dx, r0);
    sqmapwrite(atm2, "atm2_frac_trans.bin");
    fractal_inv_trans(atm2->p, nx, ny, dx, r0);
    sqmapwrite(atm2, "atm2_frac_inv_trans.bin");
    
    /*
      atm2 is u, atm is v, now comparing u'Av against v'A'u. they should be the same.
      A^-1Av should be the same as v
      A'^-1A'u should be the same as u.
    */
    dmat *u=dread("atm2_rand.bin");
    dmat *Atu=dread("atm2_frac_trans.bin");
    dmat *Av=dread("atm_frac.bin");
    dmat *v=dread("atm_rand.bin");
    double dot1=dotdbl(u->p, Av->p, NULL, nx*ny);
    double dot2=dotdbl(v->p, Atu->p, NULL, nx*ny);
    info("dot1=%g, dot2=%g, diff=%g\n", dot1, dot2, dot1-dot2);

    dmat *u2=dread("atm2_frac_inv_trans.bin");
    dmat *v2=dread("atm_frac_inv.bin");
    double d1=ddiff(u, u2);
    double d2=ddiff(v, v2);
    info("d1=%g, d2=%g\n", d1, d2);
}
static void test_cov(){//not good
    rand_t rstat;
    int seed=4;
    double r0=0.2;
    double dx=1./64;
    long N=1+1024;
    long nx=N;
    long ny=N;
    long nframe=1;
    seed_rand(&rstat, seed);
    map_t *atm=mapnew(nx, ny, dx, NULL);
    cmat *atmhat=cnew((N+1)*3,(N+1)*3);
    dmat *atmhattot=dnew((N+1)*3,(N+1)*3);
    cfft2plan(atmhat,-1);
    cfft2plan(atmhat, 1);
    dset((dmat*)atm,1);
    cembedd(atmhat, (dmat*)atm, 0);
    cfft2(atmhat, -1);
    cabs22d(&atmhattot, 1, atmhat, 1);
    ccpd(&atmhat, atmhattot);
    cifft2(atmhat, 1);
    cfftshift(atmhat);
    dmat *denom=dnew((N+1)*3,(N+1)*3);
    dmat *cov=dnew((N+1)*3,(N+1)*3);
    creal2d(&denom, 0, atmhat, 1);
    dwrite(denom, "denom.bin");
    
    dzero(atmhattot);
    for(long i=0; i<nframe; i++){
	info("%ld of %ld\n", i, nframe);
	
	for(long j=0; j<nx*ny; j++){
	    atm->p[j]=randn(&rstat);
	}
	fractal(atm->p, nx, ny, dx, r0);
	//sqmapwrite(atm, "atm_%ld.bin", i);
	cembedd(atmhat, (dmat*)atm, 0);
	cfft2(atmhat, -1);
	cabs22d(&atmhattot, 1, atmhat, 1);

	if(i==0 || (i+1)%10==0){
	    dscale(atmhattot, 1./(i+1));
	    ccpd(&atmhat, atmhattot);
	    dwrite(atmhattot, "atm_psf_%ld.bin",i+1);
	    cifft2(atmhat, 1);
	    cfftshift(atmhat);
	    creal2d(&cov, 0, atmhat, 1);
	    for(long k=0; k<cov->nx*cov->ny; k++){
		cov->p[k]/=denom->p[k];
	    }
	    dwrite(cov, "atm_cov_%ld.bin",i+1);
	}
    }
}
static void test_corner(){//ok
    rand_t rstat;
    int seed=4;
    double r0=0.2;
    double dx=32;
    long N=1+1;
    long nx=N;
    long ny=N;
    long nframe=1000000;
    seed_rand(&rstat, seed);
    map_t *atm=mapnew(nx, ny, dx, NULL);
    dmat *vec=dnew_data(atm->p, N*N, 1);
    dmat *cov=NULL;
    for(long i=0; i<nframe; i++){
	info("%ld of %ld\n", i, nframe);
	for(long j=0; j<nx*ny; j++){
	    atm->p[j]=randn(&rstat);
	}
	fractal(atm->p, nx, ny, dx, r0);
	dmm(&cov, vec, vec, "nt", 1);
    }
    dscale(cov, 1./nframe);
    dwrite(cov,"cov.bin");
}
static void test_part(){
    rand_t rstat;
    int seed=4;
    double r0=0.2;
    double dx=1./64.;
    long N=1+2;
    long ofx=0;
    long ofy=0;
    long nx=N;
    long ny=N;
    long nframe=1000000;
    seed_rand(&rstat, seed);
    map_t *atm=mapnew(nx, ny, dx, NULL);
    dmat *vec=dnew(4,1);
    dmat *cov=NULL;
    PDMAT((dmat*)atm,pp);
    for(long i=0; i<nframe; i++){
	info("%ld of %ld\n", i, nframe);
	for(long j=0; j<nx*ny; j++){
	    atm->p[j]=randn(&rstat);
	}
	fractal(atm->p, nx, ny, dx, r0);
	vec->p[0]=pp[ofy+0][ofx+0];
	vec->p[1]=pp[ofy+0][ofx+1];
	vec->p[2]=pp[ofy+1][ofx+0];
	vec->p[3]=pp[ofy+1][ofx+1];
	dmm(&cov, vec, vec, "nt", 1);
    }
    dscale(cov, 1./nframe);
    dwrite(cov,"cov.bin");
}

static void test_atm(){
    rand_t rstat;
    int seed=4;
    double r0=0.2;
    double dx=1./64;
    long N=4096;
    long nx=N;
    long ny=N;
    long nframe=10;
    seed_rand(&rstat, seed);
    {
	map_t *atm=mapnew(nx+1, ny+1, dx, NULL);
	stfun_t *data=stfun_init(nx, ny, NULL);
	for(long i=0; i<nframe; i++){
	    for(long j=0; j<(nx+1)*(ny+1); j++){
		atm->p[j]=randn(&rstat);
	    }
	    fractal(atm->p, nx+1, ny+1, dx, r0);
	    stfun_push(data, (dmat*)atm);
	    if(i%100==0){
		info("%ld of %ld\n", i, nframe);
		ddraw("fractal", (dmat*)atm, NULL, "Atmosphere","x","y","%ld",i);
	    }
	}
	dmat *st=stfun_finalize(data);
	dwrite(st, "stfun_fractal.bin");
	ddraw("fractal", st, NULL, "Atmosphere","x","y","stfun");
    }
    {
	stfun_t *data=stfun_init(nx, ny, NULL);
	dmat *spect=vonkarman_spect(nx, ny, dx, r0, 100);
	cmat *atm=cnew(nx, ny);
	cfft2plan(atm, -1);
	dmat *atmr=dnew(atm->nx, atm->ny);
	dmat *atmi=dnew(atm->nx, atm->ny);
	spect->p[0]=0;
	for(long ii=0; ii<nframe; ii+=2){
	    for(long i=0; i<atm->nx*atm->ny; i++){
		atm->p[i]=(randn(&rstat)+I*randn(&rstat))*spect->p[i];
	    }
	    cfft2(atm, -1);
	    for(long i=0; i<atm->nx*atm->ny; i++){
		atmr->p[i]=creal(atm->p[i]);
		atmi->p[i]=cimag(atm->p[i]);
	    }
	    stfun_push(data, atmr);
	    stfun_push(data, atmi);
	    if(ii%100==0){
		info("%ld of %ld\n", ii, nframe);
		ddraw("fft", atmr, NULL, "Atmosphere","x","y","%ld",ii);
	    }
	}
	dmat *st=stfun_finalize(data);
	dwrite(st, "stfun_fft.bin");
	ddraw("fft", st, NULL, "Atmosphere","x","y","stfun");
    }
	
}
int main(){
    int ind=4;
    switch(ind){
    case 0:
	test_accuracy();
	break;
    case 1:
	test_cov();
	break;
    case 2:
	test_corner();
	break;
    case 3:
	test_part();
	break;
    case 4:
	test_atm();
    }
}
