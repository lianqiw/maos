/*
  Copyright 2009-2026 Lianqi Wang

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
/*
  Unit tests for complex matrix routines (math/matcomp.c and the complex
  instances of matmath.c). Ported from test_cmat2.c, whose benchmarks and
  #if 0 blocks (cmat_benchmark, cmat_correctness, test_sq2, test_cwm,
  test_ctilt, test_ints, the commented ccwm3col/ccwm3 calls) are converted
  into correctness checks.
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-10
#endif
/*Assert two complex numbers agree to within tol.*/
#define assert_cnear(a, b, tol) assert_float_equal(cabs((a)-(b)), 0, tol)

/*Maximum absolute difference between two complex arrays of the same size.*/
static real cmaxdiff(const cmat *A, const cmat *B){
	assert_int_equal(NX(A), NX(B));
	assert_int_equal(NY(A), NY(B));
	real m=0;
	for(long i=0; i<PN(A); i++){
		real d=cabs(P(A, i)-P(B, i));
		if(d>m) m=d;
	}
	return m;
}
/*Location of the element with the largest modulus.*/
static void cargmax(long *ix, long *iy, const cmat *A){
	real m=-1;
	for(long y=0; y<NY(A); y++){
		for(long x=0; x<NX(A); x++){
			if(cabs(P(A, x, y))>m){
				m=cabs(P(A, x, y));
				*ix=x; *iy=y;
			}
		}
	}
}

/** cembed without rotation pads/crops around the FFT center so that padding then cropping is the identity. */
static void cmat_embed_pad(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	for(long n=6; n<=7; n++){
		cmat *B=cnew(n, n);
		crandn(B, 1, &rstat);
		cmat *A=cnew(2*n+2, 2*n+2);
		cembed(A, B, 0);
		const long off=(2*n+2)/2-n/2;
		real sumin=0;
		for(long y=0; y<NY(A); y++){
			for(long x=0; x<NX(A); x++){
				long bx=x-off, by=y-off;
				if(bx>=0&&bx<n&&by>=0&&by<n){
					assert_true(P(A, x, y)==P(B, bx, by));
				} else{
					assert_true(cabs(P(A, x, y))==0);
				}
				sumin+=cabs2(P(A, x, y));
			}
		}
		assert_float_equal(sumin, csumsq(B), TOL);
		cmat *B2=cnew(n, n);
		cembed(B2, A, 0);
		assert_float_equal(cmaxdiff(B, B2), 0, 0);
		cfree(A); cfree(B); cfree(B2);
	}
}

/** cembed with rotation: pi/2 permutes pixels exactly, 2*pi is the identity, and a centered disk is nearly invariant under -pi/4. */
static void cmat_embed_rotate(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long n=16;
	cmat *B=cnew(n, n);
	crandn(B, 1, &rstat);
	cmat *A=cnew(n, n);
	cembed(A, B, M_PI/2);
	//out(x,y)=in(y, -x) relative to the center (n/2, n/2)
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			long bx=iy, by=n-ix;
			if(bx>=1&&bx<=n-2&&by>=1&&by<=n-2){
				assert_cnear(P(A, ix, iy), P(B, bx, by), TOL);
			}
		}
	}
	cembed(A, B, 2*M_PI);
	for(long iy=1; iy<n-1; iy++){
		for(long ix=1; ix<n-1; ix++){
			assert_cnear(P(A, ix, iy), P(B, ix, iy), TOL);
		}
	}
	//port of cmat_benchmark: rotate a disk by -45 degree.
	cmat *C=cnew(64, 64);
	cmat *D=cnew(32, 32);
	ccircle(D, 16, 16, 1, 1, 8, 1);
	cembed(C, D, -M_PI*0.25);
	assert_float_equal(creal(csum(C)), creal(csum(D)), 0.01*creal(csum(D)));
	assert_float_equal(creal(P(C, 32, 32)), 1, TOL);
	assert_float_equal(creal(P(C, 32+12, 32)), 0, TOL);
	cfree(A); cfree(B); cfree(C); cfree(D);
}

/** cembed of a real valued complex array equals dembed of its real part and cembedd (C_REAL intent of cmat_benchmark). */
static void cmat_embed_real(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	dmat *d=dnew(20, 20);
	drandn(d, 1, &rstat);
	cmat *c=cnew(20, 20);
	for(long i=0; i<PN(d); i++) P(c, i)=P(d, i);
	const real thetas[]={0, -M_PI*0.25, 0.3};
	for(int it=0; it<3; it++){
		dmat *dout=dnew(30, 30);
		cmat *cout=cnew(30, 30);
		cmat *cout2=cnew(30, 30);
		dembed(dout, d, thetas[it]);
		cembed(cout, c, thetas[it]);
		cembedd(cout2, d, thetas[it]);
		dmat *re=NULL, *im=NULL;
		creal2d(&re, 0, cout, 1);
		cimag2d(&im, 0, cout, 1);
		assert_dmat_near(re, dout, TOL);
		assert_float_equal(dmaxabs(im), 0, 0);
		assert_float_equal(cmaxdiff(cout, cout2), 0, TOL);
		dfree(re); dfree(im); dfree(dout); cfree(cout); cfree(cout2);
	}
	dfree(d); cfree(c);
}

/** cembed_wvf builds amp*exp(2*pi*i*opd/wvl) centered in the array, with and without rotation. */
static void cmat_embed_wvf(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 3);
	const long nopd=32, npsf=64;
	const real wvl=0.5e-6;
	dmat *opd=dnew(nopd, nopd);
	dmat *amp=dnew(nopd, nopd);
	drandn(opd, 0.2e-6, &rstat);
	dcircle(amp, nopd/2, nopd/2, 1, 1, nopd/2-1, 1);
	cmat *A=cnew(npsf, npsf);
	cset(A, 5);//must be zeroed by cembed_wvf
	cembed_wvf(A, P(opd), P(amp), nopd, nopd, wvl, 0);
	const long skip=(npsf-nopd)/2;
	for(long iy=0; iy<npsf; iy++){
		for(long ix=0; ix<npsf; ix++){
			long ox=ix-skip, oy=iy-skip;
			comp expect=0;
			if(ox>=0&&ox<nopd&&oy>=0&&oy<nopd){
				real ph=2*M_PI/wvl*P(opd, ox, oy);
				expect=P(amp, ox, oy)*COMPLEX(cos(ph), sin(ph));
			}
			assert_cnear(P(A, ix, iy), expect, TOL);
		}
	}
	//Rotation by pi/2: opd/amp are interpolated at (y+nopd/2, -x+nopd/2).
	cembed_wvf(A, P(opd), P(amp), nopd, nopd, wvl, M_PI/2);
	long ncheck=0;
	for(long iy=0; iy<npsf; iy++){
		for(long ix=0; ix<npsf; ix++){
			long ox=iy-npsf/2+nopd/2, oy=-(ix-npsf/2)+nopd/2;
			if(ox>=1&&ox<=nopd-2&&oy>=1&&oy<=nopd-2){
				real ph=2*M_PI/wvl*P(opd, ox, oy);
				comp expect=P(amp, ox, oy)*COMPLEX(cos(ph), sin(ph));
				assert_cnear(P(A, ix, iy), expect, 1e-6);
				ncheck++;
			} else if(ox<-1||ox>nopd||oy<-1||oy>nopd){
				assert_float_equal(cabs(P(A, ix, iy)), 0, 0);
			}
		}
	}
	assert_int_equal(ncheck, (nopd-2)*(nopd-2));
	//-45 degree rotation of a uniform wavefront: unit modulus inside, zero far outside.
	dset(amp, 1);
	dset(opd, 0.1e-6);
	cembed_wvf(A, P(opd), P(amp), nopd, nopd, wvl, -M_PI/4);
	real ph=2*M_PI/wvl*0.1e-6;
	assert_cnear(P(A, npsf/2, npsf/2), COMPLEX(cos(ph), sin(ph)), TOL);
	assert_float_equal(cabs(P(A, npsf/2+nopd/2+2, npsf/2)), 1, TOL);//diagonal of the rotated square extends further
	assert_float_equal(cabs(P(A, 0, 0)), 0, 0);
	cfree(A); dfree(opd); dfree(amp);
}

/** Tilted wavefront gives a PSF peak at the expected offset; Parseval holds; cfftshift and cembed rotation move the peak (port of cmat_correctness). */
static void cmat_psf(void **state){
	(void)state;
	const long nopd=32, npsf=64;
	dmat *amp=dnew(nopd, nopd);
	dmat *opd=dnew(nopd, nopd);
	dcircle(amp, nopd/2, nopd/2, 1, 1, 12, 1);
	for(long iy=0; iy<nopd; iy++){
		for(long ix=0; ix<nopd; ix++){
			P(opd, ix, iy)=ix*0.25;//pi/4 phase per pixel at wvl=2
		}
	}
	cmat *A=cnew(npsf, npsf);
	cembed_wvf(A, P(opd), P(amp), nopd, nopd, 2, 0);
	const real e0=csumsq(A);
	cfft2(A, -1);
	assert_float_equal(csumsq(A), e0*npsf*npsf, TOL*e0*npsf*npsf);
	long px=-1, py=-1;
	cargmax(&px, &py, A);
	assert_int_equal(px, npsf/8);//pi/4 per pixel is 1/8 cycle per pixel.
	assert_int_equal(py, 0);
	assert_float_equal(cabs(P(A, px, py)), dsum(amp), TOL*dsum(amp));
	cfftshift(A);
	cargmax(&px, &py, A);
	assert_int_equal(px, npsf/2+npsf/8);
	assert_int_equal(py, npsf/2);
	//rotate the PSF intensity by pi/2 about the center.
	cabs2toreal(A, 1);
	cmat *B=cnew(npsf, npsf);
	cembed(B, A, M_PI/2);
	cargmax(&px, &py, B);
	assert_int_equal(px, npsf/2);
	assert_int_equal(py, npsf/2+npsf/8);
	//Rotating the wavefront by pi/4 rotates the PSF peak direction.
	cembed_wvf(A, P(opd), P(amp), nopd, nopd, 2, M_PI/4);
	cfft2(A, -1);
	cfftshift(A);
	cargmax(&px, &py, A);
	assert_int_equal(px, npsf/2+6);//8*cos(pi/4)=5.66
	assert_int_equal(py, npsf/2+6);
	cfree(A); cfree(B); dfree(amp); dfree(opd);
}

/** Different ways of computing |c|^2 agree (port of test_sq2 and the a*conj(a) check). */
static void cmat_abs2(void **state){
	(void)state;
	comp a=COMPLEX(2342, 3);
	assert_float_equal(creal(a*conj(a)), 2342.*2342.+9., 0);
	assert_float_equal(cabs2(a), 2342.*2342.+9., 0);
	const long N=64;
	cmat *C=cnew(N, N);
	ccircle(C, N/2, N/2, 1, 1, N/4, 1);
	cfft2(C, -1);
	dmat *D=NULL;
	cabs22d(&D, 0, C, 1);
	cmat *E=cdup(C);
	cabs2toreal(E, 1);
	for(long i=0; i<PN(C); i++){
		const comp tmp=P(C, i);
		const real v1=creal(tmp*conj(tmp));
		const real v2=creal(tmp)*creal(tmp)+cimag(tmp)*cimag(tmp);
		const real v3=pow(creal(tmp), 2)+pow(cimag(tmp), 2);
		const real v4=cabs2(tmp);
		const real tol=1e-12*(v2+1);
		assert_float_equal(v1, v2, tol);
		assert_float_equal(v3, v2, tol);
		assert_float_equal(v4, v2, tol);
		assert_float_equal(P(D, i), v2, tol);
		assert_float_equal(creal(P(E, i)), v2, tol);
		assert_float_equal(cimag(P(E, i)), 0, 0);
	}
	//scaled accumulation: D=2*D+3*|C|^2 = 5*|C|^2
	dmat *D0=ddup(D);
	cabs22d(&D, 2, C, 3);
	dscale(D0, 5);
	assert_dmat_near(D, D0, 1e-12*dmaxabs(D0));
	cfree(C); cfree(E); dfree(D); dfree(D0);
}

/** D.*=conj(C) with D=C gives |C|^2 as real numbers; ccwmc matches the explicit loop (port of test_cwm). */
static void cmat_cwmc(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 4);
	const long N=32;
	cmat *C=cnew(N, N);
	crandn(C, 1, &rstat);
	cmat *D=NULL;
	ccp(&D, C);
	ccwmc(D, C);
	for(long i=0; i<PN(C); i++){
		assert_float_equal(creal(P(D, i)), cabs2(P(C, i)), TOL);
		assert_float_equal(cimag(P(D, i)), 0, TOL);
	}
	cmat *B=cnew(N, N);
	crandn(B, 1, &rstat);
	ccp(&D, C);
	ccwmc(D, B);
	for(long i=0; i<PN(C); i++){
		assert_cnear(P(D, i), P(C, i)*conj(P(B, i)), TOL);
	}
	//NULL B is a no-op
	cmat *D2=cdup(D);
	ccwmc(D2, NULL);
	assert_float_equal(cmaxdiff(D, D2), 0, 0);
	//ccwmd multiplies by a real array
	dmat *w=dnew(N, N);
	drandn(w, 1, &rstat);
	ccp(&D, C);
	ccwmd(D, w);
	for(long i=0; i<PN(C); i++){
		assert_cnear(P(D, i), P(C, i)*P(w, i), TOL);
	}
	cfree(C); cfree(D); cfree(D2); cfree(B); dfree(w);
}

/** ctilt shifts the image: integer shifts are exact circular shifts, peak_corner=1 matches fftshift, fractional shifts move the centroid (port of test_ctilt). */
static void cmat_tilt(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 5);
	const long N=32;
	cmat *img=cnew(N, N);
	crandn(img, 1, &rstat);
	cmat *C=cdup(img);
	cfft2(C, -1);
	cmat *C1=cdup(C);
	//integer shift with peak in corner
	const long sx=3, sy=-5;
	ctilt(C, sx, sy, 0);
	cfft2i(C, 1);
	for(long y=0; y<N; y++){
		for(long x=0; x<N; x++){
			assert_cnear(P(C, (x+sx+N)%N, (y+sy+N)%N), P(img, x, y), TOL);
		}
	}
	//peak in center: equivalent after fftshift, for fractional shift too.
	cmat *C2=cdup(C1);
	ccp(&C, C1);
	ctilt(C, 0.3, -1.7, 0);
	cfftshift(C2);
	ctilt(C2, 0.3, -1.7, 1);
	cfftshift(C2);
	assert_float_equal(cmaxdiff(C, C2), 0, TOL*N*N);
	//ctilt2 with separate output equals in place ctilt
	cmat *C3=cnew(N, N);
	ctilt2(C3, C1, 0.3, -1.7, 0);
	assert_float_equal(cmaxdiff(C, C3), 0, TOL);
	//fractional shift of a smooth Gaussian moves its centroid and preserves energy.
	cmat *G=cnew(N, N);
	const real sigma=2;
	for(long y=0; y<N; y++){
		for(long x=0; x<N; x++){
			real r2=pow(x-N/2, 2)+pow(y-N/2, 2);
			P(G, x, y)=exp(-r2/(2*sigma*sigma));
		}
	}
	const real e0=csumsq(G);
	cfft2(G, -1);
	ctilt(G, 0.3, -0.45, 0);
	cfft2i(G, 1);
	assert_float_equal(csumsq(G), e0, TOL*e0);
	real sum=0, cx=0, cy=0;
	for(long y=0; y<N; y++){
		for(long x=0; x<N; x++){
			real v=creal(P(G, x, y));
			sum+=v; cx+=v*x; cy+=v*y;
			assert_float_equal(cimag(P(G, x, y)), 0, 1e-7);//only the tiny Nyquist component is not real
		}
	}
	assert_float_equal(cx/sum, N/2+0.3, 1e-6);
	assert_float_equal(cy/sum, N/2-0.45, 1e-6);
	cfree(img); cfree(C); cfree(C1); cfree(C2); cfree(C3); cfree(G);
}

/** ccwm, ccwmcol, ccwmrow match explicit loops; the disabled ccwm3col/ccwm3 combinations are reproduced by chaining; invalid sizes are ignored (port of bench_ccwm). */
static void cmat_cwm(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long N=64;
	cmat *A=cnew(N, N);
	cmat *B=cnew(N, N);
	crandn(A, 2, &rstat);
	crandn(B, 2, &rstat);
	ccircle(B, N/2+N/4, N/2, 1, 1, N/4, 1);
	ccircle(A, N/2+N/4, N/2, 1, 1, N/4, 1);
	cmat *A0=cdup(A);
	ccwm(A, B);
	for(long i=0; i<PN(A); i++){
		assert_cnear(P(A, i), P(A0, i)*P(B, i), TOL);
	}
	cmat *D=cnew(N, 1);
	crandn(D, 3, &rstat);
	cset(A, 1);
	ccwmcol(A, D);
	ccwm(A, B);
	for(long y=0; y<N; y++){
		for(long x=0; x<N; x++){
			assert_cnear(P(A, x, y), P(D, x)*P(B, x, y), TOL);
		}
	}
	//ccwm3col(A, W=B, D, 1, 0, 0): A=A.*W.*D per column
	ccp(&A, A0);
	ccwm(A, B);
	ccwmcol(A, D);
	for(long y=0; y<N; y++){
		for(long x=0; x<N; x++){
			assert_cnear(P(A, x, y), P(A0, x, y)*P(B, x, y)*P(D, x), 10*TOL);
		}
	}
	//ccwm3(E, A, B, 1, 0, 0) equals two ccwm.
	cmat *E=cnew(N, N);
	cset(E, 2);
	ccwm(E, A0);
	ccwm(E, B);
	for(long i=0; i<PN(E); i++){
		assert_cnear(P(E, i), 2.*P(A0, i)*P(B, i), 10*TOL);
	}
	//ccwmrow multiplies column iy by D(iy)
	ccp(&A, A0);
	ccwmrow(A, D);
	for(long y=0; y<N; y++){
		for(long x=0; x<N; x++){
			assert_cnear(P(A, x, y), P(A0, x, y)*P(D, y), TOL);
		}
	}
	//invalid sizes leave A unchanged.
	ccp(&A, A0);
	cmat *small=cnew(3, 3);
	cset(small, 7);
	ccwm(A, small);
	ccwmcol(A, small);
	ccwmrow(A, small);
	assert_float_equal(cmaxdiff(A, A0), 0, 0);
	cfree(A); cfree(A0); cfree(B); cfree(D); cfree(E); cfree(small);
}

/** ccircle fills about pi*r^2 with unit value inside and rejects invalid input. */
static void cmat_circle(void **state){
	(void)state;
	const long N=64;
	cmat *A=cnew(N, N);
	assert_int_equal(ccircle(A, N/2, N/2, 1, 1, N/4, 1), 0);
	const real area=M_PI*(N/4)*(N/4);
	assert_float_equal(creal(csum(A)), area, 0.005*area);
	assert_float_equal(cimag(csum(A)), 0, 0);
	assert_float_equal(creal(P(A, N/2, N/2)), 1, 0);
	assert_float_equal(creal(P(A, 0, 0)), 0, 0);
	//accumulates with complex value
	assert_int_equal(ccircle(A, N/2, N/2, 1, 1, N/4, COMPLEX(0, 2)), 0);
	assert_cnear(P(A, N/2, N/2), COMPLEX(1, 2), 0);
	//invalid input
	assert_int_equal(ccircle(A, N/2, N/2, 1, 1, -1, 1), -1);
	assert_int_equal(ccircle(NULL, N/2, N/2, 1, 1, 1, 1), -1);
	cfree(A);
}

/** ccpcorner copies the four corner quadrants with C_FULL, C_ABS2, C_REAL and C_ABS. */
static void cmat_cpcorner(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 6);
	const long nin=6, niny=4, nout=10, nouty=8;
	cmat *B=cnew(nin, niny);
	crandn(B, 1, &rstat);
	cmat *A=cnew(nout, nouty);
	const CEMBED flags[]={C_FULL, C_ABS2, C_REAL, C_ABS};
	for(int iflag=0; iflag<4; iflag++){
		cset(A, 9);
		ccpcorner(A, B, flags[iflag]);
		cmat *expect=cnew(nout, nouty);
		for(long y=0; y<niny; y++){
			long yo=y<niny/2?y:y-niny+nouty;
			for(long x=0; x<nin; x++){
				long xo=x<nin/2?x:x-nin+nout;
				comp v=P(B, x, y);
				switch(flags[iflag]){
				case C_FULL: P(expect, xo, yo)=v; break;
				case C_ABS2: P(expect, xo, yo)=cabs2(v); break;
				case C_REAL: P(expect, xo, yo)=creal(v); break;
				default: P(expect, xo, yo)=cabs(v); break;
				}
			}
		}
		assert_float_equal(cmaxdiff(A, expect), 0, TOL);
		cfree(expect);
	}
	//cropping the corners back recovers B.
	cmat *B2=cnew(nin, niny);
	ccpcorner(A, B, C_FULL);
	ccpcorner(B2, A, C_FULL);
	assert_float_equal(cmaxdiff(B, B2), 0, 0);
	cfree(A); cfree(B); cfree(B2);
}

/** Port of the disabled test_ints: equivalent ways of forming exp(i*k*opd), the PSF pipeline (embed, FFT, multiply, fftshift). */
static void cmat_ints(void **state){
	(void)state;
	const int nopd=32, npsf=64;
	const real wvl=0.5e-6;
	const real wvkr=2*M_PI/wvl;
	rand_t rstat;
	seed_rand(&rstat, 7);
	dmat *opd=dnew(nopd, nopd);
	dmat *amp=dnew(nopd, nopd);
	drandn(opd, 0.3e-6, &rstat);
	dset(amp, 1);
	for(int i=0; i<nopd*nopd; i++){
		const real junk=wvkr*P(opd, i);
		comp v1=P(amp, i)*cexp(COMPLEX(0, wvkr)*P(opd, i));
		comp v2=P(amp, i)*COMPLEX(cos(junk), sin(junk));
		assert_cnear(v1, v2, 1e-12);
		//cos with sqrt(1-cos^2) only agrees when sin>=0 (the sign problem noted in cembed_wvf).
		real cs=cos(junk);
		comp v3=P(amp, i)*COMPLEX(cs, sqrt(fmax(0., 1-cs*cs)));
		if(sin(junk)>=0){
			assert_cnear(v3, v2, 1e-7);
		} else{
			assert_cnear(v3, conj(v2), 1e-7);
		}
	}
	cmat *psf=cnew(npsf, npsf);
	cembed_wvf(psf, P(opd), P(amp), nopd, nopd, wvl, 0);
	cmat *psf2=cdup(psf);
	const real e0=csumsq(psf);
	assert_float_equal(e0, nopd*nopd, TOL);
	cfft2(psf, -1);
	//The PSF peak with a perfect wavefront is (number of pixels)^2; a random phase reduces it.
	assert_true(cmaxabs(psf)<nopd*nopd);
	assert_float_equal(csumsq(psf), e0*npsf*npsf, TOL*e0*npsf*npsf);
	//multiply then fftshift twice is the identity.
	cmat *psf3=cdup(psf);
	ccwm(psf, psf2);
	for(long i=0; i<PN(psf); i++){
		assert_cnear(P(psf, i), P(psf3, i)*P(psf2, i), TOL*npsf);
	}
	ccp(&psf3, psf);
	cfftshift(psf);
	cfftshift(psf);
	assert_float_equal(cmaxdiff(psf, psf3), 0, 0);
	cfree(psf); cfree(psf2); cfree(psf3); dfree(opd); dfree(amp);
}

UNITTEST_MAIN(
	cmocka_unit_test(cmat_embed_pad),
	cmocka_unit_test(cmat_embed_rotate),
	cmocka_unit_test(cmat_embed_real),
	cmocka_unit_test(cmat_embed_wvf),
	cmocka_unit_test(cmat_psf),
	cmocka_unit_test(cmat_abs2),
	cmocka_unit_test(cmat_cwmc),
	cmocka_unit_test(cmat_tilt),
	cmocka_unit_test(cmat_cwm),
	cmocka_unit_test(cmat_circle),
	cmocka_unit_test(cmat_cpcorner),
	cmocka_unit_test(cmat_ints)
)
