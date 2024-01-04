/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


#include "../lib/aos.h"
real L0=INFINITY;
long ninit=2;
static void test_accuracy(){
    rand_t rstat;
    int seed=4;
    real r0=0.2;
    real dx=1./64.;
    long N=1+64;
    long nx=N;
    long ny=N;
    seed_rand(&rstat, seed);
    map_t *atm=mapnew(nx, ny, dx, dx);
    for(long i=0; i<nx*ny; i++){
	P(atm,i)=randn(&rstat);
    }
    map_t *atm2=mapnew(nx, ny, dx, dx);
    for(long i=0; i<nx*ny; i++){
	P(atm2,i)=randn(&rstat);
    }
    mapwrite(atm, "atm_rand.bin");
    mapwrite(atm2, "atm2_rand.bin");
    fractal_do((dmat*)atm, dx, r0,L0,ninit);
    mapwrite(atm, "atm_frac.bin");
    fractal_inv((dmat*)atm, dx, r0,L0,ninit);
    mapwrite(atm, "atm_frac_inv.bin");

   
    fractal_trans((dmat*)atm2, dx, r0,L0,ninit);
    mapwrite(atm2, "atm2_frac_trans.bin");
    fractal_inv_trans((dmat*)atm2, dx, r0,L0,ninit);
    mapwrite(atm2, "atm2_frac_inv_trans.bin");
    
    /*
      atm2 is u, atm is v, now comparing u'Av against v'A'u. they should be the same.
      A^-1Av should be the same as v
      A'^-1A'u should be the same as u.
    */
    dmat *u=dread("atm2_rand.bin");
    dmat *Atu=dread("atm2_frac_trans.bin");
    dmat *Av=dread("atm_frac.bin");
    dmat *v=dread("atm_rand.bin");
    real dot1=dvecdot(P(u), P(Av), NULL, nx*ny);
    real dot2=dvecdot(P(v), P(Atu), NULL, nx*ny);
    dbg("dot1=%g, dot2=%g, diff=%g\n", dot1, dot2, dot1-dot2);

    dmat *u2=dread("atm2_frac_inv_trans.bin");
    dmat *v2=dread("atm_frac_inv.bin");
    real d1=ddiff(u, u2);
    real d2=ddiff(v, v2);
    dbg("d1=%g, d2=%g\n", d1, d2);
}
static void test_cov(){/*not good */
    rand_t rstat;
    int seed=4;
    real r0=0.2;
    real dx=1./64;
    long N=1+1024;
    long nx=N;
    long ny=N;
    long nframe=1;
    seed_rand(&rstat, seed);
    map_t *atm=mapnew(nx, ny, dx,dx);
    cmat *atmhat=cnew((N+1)*3,(N+1)*3);
    dmat *atmhattot=dnew((N+1)*3,(N+1)*3);
    //cfft2plan(atmhat,-1);
    //cfft2plan(atmhat, 1);
    dset((dmat*)atm,1);
    cembedd(atmhat, (dmat*)atm, 0);
    cfft2(atmhat, -1);
    cabs22d(&atmhattot, 1, atmhat, 1);
    ccpd(&atmhat, atmhattot);
    cfft2i(atmhat, 1);
    cfftshift(atmhat);
    dmat *denom=dnew((N+1)*3,(N+1)*3);
    dmat *cov=dnew((N+1)*3,(N+1)*3);
    creal2d(&denom, 0, atmhat, 1);
    writebin(denom, "denom.bin");
    
    dzero(atmhattot);
    for(long i=0; i<nframe; i++){
	dbg("%ld of %ld\n", i, nframe);
	
	for(long j=0; j<nx*ny; j++){
	    P(atm,j)=randn(&rstat);
	}
	fractal_do((dmat*)atm, dx, r0,L0,ninit);
	/*mapwrite(atm, "atm_%ld.bin", i); */
	cembedd(atmhat, (dmat*)atm, 0);
	cfft2(atmhat, -1);
	cabs22d(&atmhattot, 1, atmhat, 1);

	if(i==0 || (i+1)%10==0){
	    dscale(atmhattot, 1./(i+1));
	    ccpd(&atmhat, atmhattot);
	    writebin(atmhattot, "atm_psf_%ld.bin",i+1);
	    cfft2i(atmhat, 1);
	    cfftshift(atmhat);
	    creal2d(&cov, 0, atmhat, 1);
	    for(long k=0; k<cov->nx*cov->ny; k++){
		P(cov,k)/=P(denom,k);
	    }
	    writebin(cov, "atm_cov_%ld.bin",i+1);
	}
    }
}
static void test_corner(){/*Compute the covariance of 4 corner points*/
    rand_t rstat;
    int seed=4;
    real r0=0.2;
    real dx=32;
    long N=1+1;
    long nx=N;
    long ny=N;
    long nframe=1000000;
    seed_rand(&rstat, seed);
    map_t *atm=mapnew(nx, ny, dx, dx);
    dmat *vec=dref_reshape((dmat*)atm, N*N, 1);
    dmat *cov=NULL;
    for(long i=0; i<nframe; i++){
	dbg("%ld of %ld\n", i, nframe);
	for(long j=0; j<nx*ny; j++){
	    P(atm,j)=randn(&rstat);
	}
	fractal_do((dmat*)atm, dx, r0,L0,ninit);
	dmm(&cov, 1, vec, vec, "nt", 1);
    }
    dscale(cov, 1./nframe);
    writebin(cov,"cov.bin");
}
static void test_part(){/**Compute the covariance of 4 points with various separation.*/
    rand_t rstat;
    int seed=4;
    real r0=0.2;
    real dx=1./64.;
    long N=1+2;
    long ofx=0;
    long ofy=0;
    long nx=N;
    long ny=N;
    long nframe=1000000;
    seed_rand(&rstat, seed);
    map_t *atm=mapnew(nx, ny, dx,dx);
    dmat *vec=dnew(4,1);
    dmat *cov=NULL;
    dmat* pp=(dmat*)atm;
    for(long i=0; i<nframe; i++){
	dbg("%ld of %ld\n", i, nframe);
	for(long j=0; j<nx*ny; j++){
	    P(atm,j)=randn(&rstat);
	}
	fractal_do((dmat*)atm, dx, r0,L0,ninit);
	P(vec,0)=P(pp,ofx+0,ofy+0);
	P(vec,1)=P(pp,ofx+1,ofy+0);
	P(vec,2)=P(pp,ofx+0,ofy+1);
	P(vec,3)=P(pp,ofx+1,ofy+1);
	dmm(&cov, 1, vec, vec, "nt", 1);
    }
    dscale(cov, 1./nframe);
    writebin(cov,"cov.bin");
}

static void test_stfun(){
    rand_t rstat;
    int seed=4;
    real r0=0.2;
    real dx=1./16;
    long N=32;
    long nx=N;
    long ny=N;
    long nframe=500;
    seed_rand(&rstat, seed);
    if(L0<9000){
	dmat *rr=dlinspace(0, N*dx, N);
	dmat *covvk=turbcov(rr, sqrt(2)*N*dx, r0, L0);
	writebin(covvk, "cov_vk");
	dfree(rr);
	dfree(covvk);
    }
    /*    return; */
    {
	map_t *atm=mapnew(nx+1, ny+1, dx, dx);
	stfun_t *data=stfun_init(nx, ny, NULL);
	zfarr *save=zfarr_init(nframe, 1, "fractal_atm.bin");
	for(long i=0; i<nframe; i++){
	    for(long j=0; j<(nx+1)*(ny+1); j++){
		P(atm,j)=randn(&rstat);
	    }
	    fractal_do((dmat*)atm, dx, r0,L0,ninit);
	    stfun_push(data, (dmat*)atm);
	    zfarr_push(save, i, (dmat*)atm);
	    if(i%100==0)
		dbg("%ld of %ld\n", i, nframe);
	}
	zfarr_close(save);
	dmat *st=stfun_finalize(data);
	writebin(st, "stfun_fractal.bin");
	ddraw("fractal", st, NULL,NULL, "Atmosphere","x","y","stfun");
    }
    /*exit(0); */
    {
	stfun_t *data=stfun_init(nx, ny, NULL);
	dmat *spect=turbpsd(nx, ny, dx, r0, 100, 0, 0.5);
	cmat *atm=cnew(nx, ny);
	//cfft2plan(atm, -1);
	dmat *atmr=dnew(atm->nx, atm->ny);
	dmat *atmi=dnew(atm->nx, atm->ny);
	P(spect,0)=0;
	for(long ii=0; ii<nframe; ii+=2){
	    for(long i=0; i<atm->nx*atm->ny; i++){
		P(atm,i)=COMPLEX(randn(&rstat), randn(&rstat))*P(spect,i);
	    }
	    cfft2(atm, -1);
	    for(long i=0; i<atm->nx*atm->ny; i++){
		P(atmr,i)=creal(P(atm,i));
		P(atmi,i)=cimag(P(atm,i));
	    }
	    stfun_push(data, atmr);
	    stfun_push(data, atmi);
	    if(ii%100==0)
		dbg("%ld of %ld\n", ii, nframe);
	}
	dmat *st=stfun_finalize(data);
	writebin(st, "stfun_fft.bin");
	ddraw("fft", st, NULL,NULL, "Atmosphere","x","y","stfun");
    }
	
}

static void test_psd(){
    rand_t rstat;
    int seed=4;
    real r0=0.2;
    real dx=1./64;
    long N=1024;
    long nx=N;
    long ny=N;
    long ratio=1;
    long xskip=nx*(ratio-1)/2;
    long yskip=ny*(ratio-1)/2;
    long nframe=512;
    seed_rand(&rstat, seed);
    if(1){
	map_t *atm=mapnew(nx+1, ny+1, dx,dx);
	cmat *hat=cnew(nx*ratio, ny*ratio);
	//cfft2plan(hat, -1);
	dmat *hattot=dnew(nx*ratio, ny*ratio);

	for(long i=0; i<nframe; i++){
	    info("%ld of %ld\n", i, nframe);
	    for(long j=0; j<(nx+1)*(ny+1); j++){
		P(atm,j)=randn(&rstat);
	    }
	    fractal_do((dmat*)atm, dx, r0,L0,ninit);
	    czero(hat);
	    for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
		    P(hat,ix+xskip,iy+yskip)=P(atm,ix,iy);
		}
	    }
	    cfftshift(hat);
	    cfft2i(hat, -1);
	    cabs22d(&hattot, 1, hat, 1);
	}
	dscale(hattot, 1./nframe);
	dfftshift(hattot);
	writebin(hattot, "PSD_fractal");
    }
    {
	dmat *spect=turbpsd(nx, ny, dx, r0, 100, 0, 0.5);
	writebin(spect, "spect");
	cmat *hat=cnew(nx*ratio, ny*ratio);
	//cfft2plan(hat, -1);
	dmat *hattot=dnew(nx*ratio, ny*ratio);
	cmat *atm=cnew(nx, ny);
	//cfft2plan(atm, -1);
	dmat *atmr=dnew(atm->nx, atm->ny);
	dmat *atmi=dnew(atm->nx, atm->ny);
	cmat*  phat=hat;
	dmat*  patmr=atmr;
	dmat*  patmi=atmi;
	for(long ii=0; ii<nframe; ii+=2){
	    info("%ld of %ld\n", ii, nframe);
	    for(long i=0; i<atm->nx*atm->ny; i++){
		P(atm,i)=COMPLEX(randn(&rstat), randn(&rstat))*P(spect,i);
	    }
	    cfft2(atm, -1);
	    for(long i=0; i<atm->nx*atm->ny; i++){
		P(atmr,i)=creal(P(atm,i));
		P(atmi,i)=cimag(P(atm,i));
	    }
	    czero(hat);
	    for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
		    P(phat,ix+xskip,iy+yskip)=P(patmr,ix,iy);
		}
	    }
	    cfftshift(hat);
	    cfft2i(hat, -1);
	    cabs22d(&hattot, 1, hat, 1);
	    czero(hat);
	    for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
		    P(phat,ix+xskip,iy+yskip)=P(patmi,ix,iy);
		}
	    }
	    cfftshift(hat);
	    cfft2i(hat, -1);
	    cabs22d(&hattot, 1, hat, 1);
	}
	dscale(hattot, 1./nframe);
	dfftshift(hattot);
	writebin(hattot, "PSD_fft");
    }
}
/*
  Compute cxx on atm to compare against L2, invpsd, fractal.
*/
static void test_cxx(){
    rand_t rstat;
    int seed=4;
    real r0=0.2;
    real dx=1./4;
    long N=16;
    long nx=N;
    long ny=N;
    long nframe=40960;
    seed_rand(&rstat, seed);
    {
	dmat *cxx=dnew(N*N,N*N);
	map_t *atm=mapnew(nx+1, ny+1, dx, dx);
	for(long i=0; i<nframe; i++){
	    dbg("%ld of %ld\n", i, nframe);
	    for(long j=0; j<(nx+1)*(ny+1); j++){
		P(atm,j)=randn(&rstat);
	    }
	    fractal_do((dmat*)atm, dx, r0, L0, ninit);
	    dmat *sec=dsub((dmat*)atm, 0, nx, 0, ny);
	    dmat *atmvec=dref_reshape(sec, nx*ny, 1);
	    dmm(&cxx,1, atmvec,atmvec,"nt",1);
	    dfree(atmvec);
	    dfree(sec);
	}
	dscale(cxx, 1./nframe);
	writebin(cxx, "cxx_fractal");
	dfree(cxx);
	mapfree(atm);
    }
    {
	dmat *cxx=dnew(N*N,N*N);
	dmat *spect=turbpsd(nx, ny, dx, r0, 100, 0, 0.5);
	P(spect,0)=P(spect,1);
	cmat *atm=cnew(nx, ny);
	//cfft2plan(atm, -1);
	dmat *atmr=dnew(nx*ny,1);
	dmat *atmi=dnew(nx*ny,1);
	for(long ii=0; ii<nframe; ii+=2){
	    dbg("%ld of %ld\n", ii, nframe);
	    for(long i=0; i<atm->nx*atm->ny; i++){
		P(atm,i)=COMPLEX(randn(&rstat), randn(&rstat))*P(spect,i);
	    }
	    cfft2(atm, -1);
	    for(long i=0; i<atm->nx*atm->ny; i++){
		P(atmr,i)=creal(P(atm,i));
		P(atmi,i)=cimag(P(atm,i));
	    }
	    dmm(&cxx,1, atmr,atmr,"nt",1);
	    dmm(&cxx,1, atmi,atmi,"nt",1);
	}
	dscale(cxx, 1./nframe);
	writebin(cxx, "cxx_fft");
	dfree(cxx);
	dfree(atmr);
	dfree(atmi);
	cfree(atm);
    }
    loc_t *loc=mksqloc_auto(16,16,1./4,1./4);
    locwrite(loc,"loc");
    dmat *B=stfun_kolmogorov(loc, r0);
    writebin(B, "B_theory");
}
int main(){
    int ind=0;
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
	test_stfun();
	break;
    case 5:
	test_psd();
	break;
    case 6:
	test_cxx();
    }
}
