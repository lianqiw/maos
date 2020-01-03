/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
TIC;
#if 0
static static void cmat_benchmark(){
    cmat *A=cnew(64,64);
    cmat *B=cnew(32,32);
    dmat *D=dnew(32,32);
    ccircle(B,36,36,16,1);
    for(int i=0; i<D->nx*D->ny; i++){
	D->p[i]=rand()*0.5e-6*2*M_PI;
    }
    tic;
    for(int i=0; i<10000;i++){
	cembed(A,B,0,C_FULL);
	/*initially: 0.008 ms in polaris. */
	/*after change interpolation: 0.009 ms. */
    }
    toc("cembed 0");
    tic;
    for(int i=0; i<10000;i++){
	cembed(A,B,-M_PI*0.25,C_FULL);
	/*3.3 ms in polaris. with switch. */
	/*after change interpolation: 0.098 ms in polaris */
    }
    toc("cembed -45");
    tic;
    tic;
    for(int i=0; i<10000;i++){
	cembed(A,B,-M_PI*0.25,C_REAL);
	/*convert xy to abs/phase, do interpolation, and change back. very slow. */
	/*takes 0.436ms. astatic void using this mode. set flag to 2. */
    }
    toc("cembed -45, abs/pi");
    tic;
    for(int i=0; i<10000;i++){
	cembed_wvf(A,D->p,D->p,D->nx,D->ny,0.5e-6,0);
	/*not usint LUT: */
	/*0.067 ms in polaris if D is zero. */
	/*0.274 ms in polaris if D is non zero; */
	/*using LUT: */
	/*0.049 ms if D is zero */
	/*0.066 ms if D is non zero. */
    }
    toc("cembed_wvf 0");
    tic;
    for(int i=0; i<10000;i++){
	cembed_wvf(A,D->p,D->p,D->nx,D->ny,0.5e-6,-45);
	/*0.127 ms in polaris. */
	/*0.138 ms after change interpolation. */
    }
    toc("cembed_wvf -45");
    /*cdraw("test_cmat",B,"B"); */
    cdraw("test_cmat",A,"A");
    tic;
    for(int i=0; i<10000;i++){
	cfftshift(A,C_REAL);/*0.014 ms in polaris. */
    }
    toc("cfftshift");
}
static static void cmat_correctness(static void){

    cmat *A=cnew(64,64);
    dmat *D=dnew(64,64);
    dmat *P=dnew(64,64);
    for(int iy=0;iy<P->ny; iy++){
	for(int ix=0; ix<P->nx; ix++){
	    P->p[ix+iy*P->nx]=ix*0.25;
	}
    }
    //cfft2plan(A,-1);
    dcircle(D, 40,32, 10, 1);
    ddraw("test_cmat",D,"Pupil");
    cembed_wvf(A,P->p,D->p,D->nx,D->ny,2,M_PI*0.25);
    cdraw("test_cmat",A,"Pupil");
    cfft2(A,-1);
    cdraw("test_cmat",A,"PSF");
    cfftshift(A,C_ABS2);
    cdraw("test_cmat",A,"A PSF shifted");
    writebin(A,"PSF");
    cmat *B=cnew(64,64);
    cembed(B,A,M_PI*0.75,C_FULL);
    writebin(B,"PSFrot");
    cdraw("test_cmat",B,"B embed");
    cembedscaleout(B,A,0.5,1,M_PI*0.75,C_FULL);
    cdraw("test_cmat",B,"B embedscale");
}
static static void test_sq2(){
    int N=4096;
    cmat *restrict C=cnew(N,N);
    cmat *restrict D=cnew(N,N);
    ccircle(C,N/2,N/2,N/4,1);
    //cfft2plan(C,-1);
    cfft2(C,-1);
    tic;
    for(int i=0; i<C->nx*C->ny; i++){
	D->p[i]=C->p[i]*conj(C->p[i]);
    }
    toc("*conj");
    tic;
    for(int i=0; i<C->nx*C->ny; i++){
	D->p[i]=C->p[i]*conj(C->p[i]);/*0.26 */
    }
    toc("*conj");
    tic;
    for(int i=0; i<C->nx*C->ny; i++){
	const comp tmp=C->p[i];
	D->p[i]=creal(tmp*conj(tmp));/*0.26 */
    }
    toc("*conj2");
    tic;
    for(int i=0; i<C->nx*C->ny; i++){
	const comp tmp=C->p[i];
	D->p[i]=creal(tmp)*creal(tmp)+cimag(tmp)*cimag(tmp);/*0.23 */
    }
    toc("real*real+imag*imag");
    tic;
    for(int i=0; i<C->nx*C->ny; i++){
	const comp tmp=C->p[i];
	D->p[i]=pow(creal(tmp),2)+pow(cimag(tmp),2);/*0.23 */
    }
    toc("pow real+pow imag");
    tic;
    for(int i=0; i<C->nx*C->ny; i++){
	D->p[i]=cabs2(C->p[i]);
    }
    toc("cabs2");
}
static static void test_cwm(){
    int N=4096;
    cmat *restrict C=cnew(N,N);
    cmat * D=cnew(N,N);
    ccircle(C,N/2,N/2,N/4,1);
    //cfft2plan(C,-1);
    cfft2(C,-1);
    ccp(&D,C);
    tic;
    for(int i=0; i<C->nx*C->ny; i++){
	D->p[i]*=conj(C->p[i]);
    }
    toc("*conj");
     tic;
    for(int i=0; i<C->nx*C->ny; i++){
	D->p[i]*=conj(C->p[i]);/*0.28 */
    }
    toc("*conj");
}

static void test_ctilt(){
    int N=64;
    cmat *C=cnew(N,N);
    ccircle(C,N/2,N/2,N/4,1);
    cdraw("test_cmat",C,"Cir");
    writebin(C,"C_psf");
    //cfft2plan(C,-1);
    //cfft2plan(C,1);
    cfft2(C,-1);
    /*cfftshift(C,C_FULL); */
    writebin(C,"C_otf");
    ctilt(C,0.1,0.1,0);
    writebin(C,"C_otf_tilt");
    /*cfftshift(C,C_FULL); */
    cfft2(C,1);
    writebin(C,"C_psf_shift");
    cscale(C,1./(N*N));
    cdraw("test_cmat",C,"shift");
}
#endif
static void bench_ccwm(void){
    int N=1024*4;
    rand_t strand;
    seed_rand(&strand,1);
    cmat *A=cnew(N,N);
    cmat *B=cnew(N,N);
    crandn(A,2,&strand);
    crandn(B,2,&strand);
    ccircle(B,N/2+N/4,N/2,1,1,N*1/4,1);
    ccircle(A,N/2+N/4,N/2,1,1,N*1/4,1);
    tic;
    ccwm(A,B);
    toc("ccwm");
    tic;
    ccwm(A,B);
    toc("ccwm");
    tic;
    cmat *D=cnew(N,1);
    crandn(D,3,&strand);
    cset(A,1);
    tic;
    ccwmcol(A,D);
    ccwm(A,B);
    toc("ccwm&ccwmcol");
    writebin(A,"A1.bin");
    tic;
    ccwmcol(A,D);
    ccwm(A,B);
    toc("ccwm&ccwmcol");
    cset(A,1);
    tic;
    ccwm3col(A,B,D,1,0,0);
    toc("ccwm3col");
    writebin(A,"A2.bin");
    tic;
    ccwm3col(A,B,D,1,0,0);
    toc("ccwm3col");
    cmat *E=cnew(N,N);
    czero(E);
    tic;
    ccwm(E,A);
    ccwm(E,B);
    toc("ccwm twice");
    writebin(E,"E.bin");
    tic;
    ccwm3(E,A,B,1,0,0);
    toc("ccwm3");
    writebin(E,"E2.bin");
}
/*TIC; */
#if 0
int test_ints(){
    int nopd=32;
    int npsf=64;
    real *opd=mycalloc(nopd*nopd,real);
    real *amp=mycalloc(nopd*nopd,real);
    comp *psf=mycalloc(npsf*npsf,comp);
    comp *psf2=mycalloc(npsf*npsf,comp);
    real wvkr=2*M_PI/0.5e-6;
    comp wvk=COMPLEX(0, wvkr);
    for(int i=0; i<nopd*nopd; i++){
	amp[i]=1;
	opd[i]=1;
	psf[i]=amp[i]*cexp(wvk*opd[i]);
	psf2[i]=amp[i]*cexp(wvk*opd[i]);
    }
    /*    tic; */
    for(int j=0; j<10000; j++){
	for(int i=0; i<nopd*nopd; i++){
	    psf[i]=amp[i]*cexp(wvk*opd[i]);
	}
    }
    /*takes 0.15 ms in debug mode, 0.12 ms in O3 mode */
    /*tic; */
    for(int j=0; j<10000; j++){
	for(int i=0; i<nopd*nopd; i++){
	    real junk=wvkr*opd[i];
	    psf[i]=amp[i]*COMPLEX(cos(junk), sin(junk));
	}
    }
    /*takes 0.9 ms in debug mode, 0.09ms in O3 mode. */
    /*tic; */
    for(int j=0; j<10000; j++){
	for(int i=0; i<nopd*nopd; i++){
	    real junk=cos(wvkr*opd[i]);
	    psf[i]=amp[i]*COMPLEX(junk, sqrt(1-junk*junk));
	}
    }
    /*takes 0.9 ms in debug mode, 0.09ms in O3 mode. */
 
    /*tic; */
    for(int i=0; i<10000; i++){
	embed_wvf(psf,npsf,opd,amp,nopd,0.5e-6,0.5);
    }
    toc("embed");
  
    fftw_plan plan_psf=fftw_plan_dft_2d
	(npsf,npsf,psf,psf,-1,FFTW_MEASURE);
    /*tic; */
    for(int i=0; i<10000; i++){
	fftw_execute(plan_psf);
    }
    /*takes 0.07ms todo 64x64 FFT. */
    /*tic; */
    for(int j=0; j<10000; j++){
	for(int i=0; i<npsf*npsf; i++){
	    psf[i]*=psf2[i];
	}
    }
    /*takes 0.08ms in debug mode, 0.03 in O3 mode. */
 
    for(int j=0; j<10000; j++){
	fftshift(psf,npsf,npsf);
    }
    toc("fftshift");/*0.014 ms in O3 */
    
    tic;
    for(int j=0; j<10000; j++){
	add_hf(psf,npsf,npsf);
    }
    toc("addhf");/*0.011 ms in O3 mode */
   
    }
#endif
int main(){
#if 0
    cmat_benchmark();/*passed */
    cmat_correctness();/*passed */
    test_sq2();/*passed */
    test_cwm();/*passed */
    comp a=COMPLEX(2342,3);
    dbg("a*conj(a)=%g\n",creal(a*conj(a)));
    dbg("abs2(a)=%g\n", cabs2(a));
    test_ctilt();
#endif
    /*bench_cembed(); */
    bench_ccwm();
}
