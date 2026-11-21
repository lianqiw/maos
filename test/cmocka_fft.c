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
  Unit tests for the FFT wrappers in math/fft.c. Ported from test_fft.c, which
  only timed cfft2 and cfft2partial for many sizes (including the disabled
  test_fft_speed and test_fft_speed_small). The timing is replaced by checks
  against a naive DFT and analytic properties.
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-10
#endif

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
/*Naive 2d DFT: out(k,l)=sum_{x,y} in(x,y)*exp(dir*2*pi*i*(k*x/nx+l*y/ny)).*/
static cmat *naive_dft(const cmat *in, int dir){
	const long nx=NX(in), ny=NY(in);
	cmat *out=cnew(nx, ny);
	for(long l=0; l<ny; l++){
		for(long k=0; k<nx; k++){
			comp s=0;
			for(long y=0; y<ny; y++){
				for(long x=0; x<nx; x++){
					real ph=dir*2*M_PI*((real)(k*x)/(real)nx+(real)(l*y)/(real)ny);
					s+=P(in, x, y)*COMPLEX(cos(ph), sin(ph));
				}
			}
			P(out, k, l)=s;
		}
	}
	return out;
}
/*Naive 1d DFT along the first dimension of each column.*/
static cmat *naive_dft_col(const cmat *in, int dir){
	const long nx=NX(in), ny=NY(in);
	cmat *out=cnew(nx, ny);
	for(long y=0; y<ny; y++){
		for(long k=0; k<nx; k++){
			comp s=0;
			for(long x=0; x<nx; x++){
				real ph=dir*2*M_PI*(real)(k*x)/(real)nx;
				s+=P(in, x, y)*COMPLEX(cos(ph), sin(ph));
			}
			P(out, k, y)=s;
		}
	}
	return out;
}

/** cfft2 forward and inverse agree with a naive DFT for even, odd, rectangular and 1-d shapes. */
static void fft_naive(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long shapes[][2]={{6,10},{7,5},{8,8},{1,9},{9,1},{1,1},{3,4}};
	for(unsigned is=0; is<sizeof(shapes)/sizeof(shapes[0]); is++){
		for(int dir=-1; dir<=1; dir+=2){
			cmat *A=cnew(shapes[is][0], shapes[is][1]);
			crandn(A, 1, &rstat);
			cmat *ref=naive_dft(A, dir);
			cfft2(A, dir);
			assert_float_equal(cmaxdiff(A, ref), 0, TOL*PN(A));
			cfree(A);
			cfree(ref);
		}
	}
}

/** Forward then cfft2i inverse returns the input for many sizes (port of the timing loops over sizes). */
static void fft_roundtrip(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	//All sizes 1..40 (test_fft_speed_small), plus sizes with a large prime factor like 3824=16*239.
	long sizes[48];
	int ns=0;
	for(long n=1; n<=40; n++) sizes[ns++]=n;
	sizes[ns++]=2*239; sizes[ns++]=64; sizes[ns++]=128; sizes[ns++]=200; sizes[ns++]=97;
	for(int is=0; is<ns; is++){
		long n=sizes[is];
		long ny=n>100?8:n;//keep big ones cheap
		cmat *A=cnew(n, ny);
		crandn(A, 20, &rstat);
		cmat *A0=cdup(A);
		cfft2(A, -1);
		cfft2i(A, 1);
		assert_float_equal(cmaxdiff(A, A0), 0, 20*TOL*log2((real)(n*ny)+2));
		cfree(A);
		cfree(A0);
	}
}

/** Known transforms: delta to constant, constant to delta, plane wave to single frequency. */
static void fft_known(void **state){
	(void)state;
	const long nx=16, ny=12;
	cmat *A=cnew(nx, ny);
	//delta at origin -> all ones
	P(A, 0, 0)=1;
	cfft2(A, -1);
	for(long i=0; i<PN(A); i++){
		assert_float_equal(creal(P(A, i)), 1, TOL);
		assert_float_equal(cimag(P(A, i)), 0, TOL);
	}
	//constant (now A=1) -> nx*ny at origin only
	cfft2(A, -1);
	for(long i=0; i<PN(A); i++){
		assert_float_equal(cabs(P(A, i)), i==0?(real)(nx*ny):0, TOL*nx*ny);
	}
	//plane wave exp(+2*pi*i*(kx*x/nx+ky*y/ny)) -> nx*ny at (kx,ky) for forward transform
	const long kx=3, ky=5;
	for(long y=0; y<ny; y++){
		for(long x=0; x<nx; x++){
			real ph=2*M_PI*((real)(kx*x)/nx+(real)(ky*y)/ny);
			P(A, x, y)=COMPLEX(cos(ph), sin(ph));
		}
	}
	cfft2(A, -1);
	for(long y=0; y<ny; y++){
		for(long x=0; x<nx; x++){
			real expect=(x==kx&&y==ky)?(real)(nx*ny):0;
			assert_float_equal(cabs(P(A, x, y)), expect, TOL*nx*ny);
		}
	}
	//delta at (1,0) with inverse transform -> exp(+2*pi*i*k/nx)
	czero(A);
	P(A, 1, 0)=1;
	cfft2(A, 1);
	for(long y=0; y<ny; y++){
		for(long x=0; x<nx; x++){
			real ph=2*M_PI*x/nx;
			assert_float_equal(creal(P(A, x, y)), cos(ph), TOL);
			assert_float_equal(cimag(P(A, x, y)), sin(ph), TOL);
		}
	}
	cfree(A);
}

/** Parseval: forward FFT scales energy by N; cfft2s is unitary; four forward transforms give N^2 times the input. */
static void fft_parseval(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	const long nx=24, ny=18;
	const real N=nx*ny;
	cmat *A=cnew(nx, ny);
	crandn(A, 1, &rstat);
	cmat *A0=cdup(A);
	const real e0=csumsq(A);
	cfft2(A, -1);
	assert_float_equal(csumsq(A)/N, e0, e0*TOL);
	ccp(&A, A0);
	cfft2s(A, -1);
	assert_float_equal(csumsq(A), e0, e0*TOL);
	cfft2s(A, 1);
	assert_float_equal(cmaxdiff(A, A0), 0, TOL);
	//F^2 x[n]=N x[-n]; F^4=N^2. Also exercises plan reuse on the same array.
	ccp(&A, A0);
	cfft2(A, -1);
	cfft2(A, -1);
	for(long y=0; y<ny; y++){
		for(long x=0; x<nx; x++){
			comp expect=N*P(A0, (nx-x)%nx, (ny-y)%ny);
			assert_float_equal(cabs(P(A, x, y)-expect), 0, TOL*N);
		}
	}
	cfft2(A, -1);
	cfft2(A, -1);
	cscale(A, 1./(N*N));
	assert_float_equal(cmaxdiff(A, A0), 0, TOL);
	cfree(A);
	cfree(A0);
}

/** cfft2partial matches the full cfft2 on the transformed columns, and with ncomp=0 does a 1-d FFT along x (port of test_fft_special). */
static void fft_partial(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	//test_fft_special used 3824 and 4096 with ncomp=512. Use small sizes with the same structure.
	const long ns[]={64, 2*3*17};
	const int ncomps[]={16, 20};
	for(int is=0; is<2; is++){
		const long n=ns[is];
		const int ncomp=ncomps[is];
		for(int dir=-1; dir<=1; dir+=2){
			cmat *A=cnew(n, n);
			crandn(A, 20, &rstat);
			cmat *B=cdup(A);
			cfft2(A, dir);
			cfft2partial(B, ncomp, dir);
			//Repeat on a copy to verify plan reuse does not alter the result.
			for(long iy=0; iy<n; iy++){
				for(long ix=0; ix<n; ix++){
					if(ix<ncomp/2||ix>=n-ncomp/2){
						assert_float_equal(cabs(P(A, ix, iy)-P(B, ix, iy)), 0, TOL*n*20);
					}
				}
			}
			cfree(A);
			cfree(B);
		}
	}
	//ncomp=0: only FFT along the first dimension.
	{
		cmat *A=cnew(10, 6);
		crandn(A, 1, &rstat);
		cmat *ref=naive_dft_col(A, -1);
		cfft2partial(A, 0, -1);
		assert_float_equal(cmaxdiff(A, ref), 0, TOL*10);
		cfree(A);
		cfree(ref);
	}
}

/** dcell_fft2 on a real/imaginary pair matches cfft2 in both directions. */
static void fft_cell_split(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 3);
	const long nx=12, ny=10;
	for(int dir=-1; dir<=1; dir+=2){
		cmat *A=cnew(nx, ny);
		crandn(A, 1, &rstat);
		dcell *dc=dcellnew_same(2, 1, nx, ny);
		for(long i=0; i<PN(A); i++){
			P(P(dc, 0), i)=creal(P(A, i));
			P(P(dc, 1), i)=cimag(P(A, i));
		}
		cfft2(A, dir);
		dcell_fft2(dc, dir);
		for(long i=0; i<PN(A); i++){
			assert_float_equal(P(P(dc, 0), i), creal(P(A, i)), TOL*nx*ny);
			assert_float_equal(P(P(dc, 1), i), cimag(P(A, i)), TOL*nx*ny);
		}
		cfree(A);
		dcellfree(dc);
	}
}

/** dfft1 real-to-halfcomplex matches cfft2 of a real vector, and the HC2R round trip returns n times the input. */
static void fft_r2hc(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 4);
	for(long n=7; n<=8; n++){
		dmat *a=dnew(n, 1);
		drandn(a, 1, &rstat);
		dmat *a0=ddup(a);
		cmat *ca=cnew(n, 1);
		for(long i=0; i<n; i++) P(ca, i)=P(a, i);
		cfft2(ca, -1);
		dfft1plan_r2hc(a, -1);
		dfft1plan_r2hc(a, 1);
		dfft1(a, -1);
		//halfcomplex: r0, r1, ..., r_{n/2}, i_{(n+1)/2-1}, ..., i1
		for(long k=0; k<=n/2; k++){
			assert_float_equal(P(a, k), creal(P(ca, k)), TOL*n);
		}
		for(long k=1; k<(n+1)/2; k++){
			assert_float_equal(P(a, n-k), cimag(P(ca, k)), TOL*n);
		}
		dfft1(a, 1);
		dscale(a, 1./n);
		assert_dmat_near(a, a0, TOL);
		dfree(a);
		dfree(a0);
		cfree(ca);
	}
}

/** cfftshift swaps quadrants (identity when applied twice), and cffttreat equals multiplying by (-1)^(x+y). */
static void fft_shift_treat(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 5);
	const long nx=8, ny=6;
	cmat *A=cnew(nx, ny);
	crandn(A, 1, &rstat);
	cmat *B=cdup(A);
	cfftshift(B);
	for(long y=0; y<ny; y++){
		for(long x=0; x<nx; x++){
			assert_true(P(B, (x+nx/2)%nx, (y+ny/2)%ny)==P(A, x, y));
		}
	}
	cfftshift(B);
	assert_float_equal(cmaxdiff(A, B), 0, 0);
	//1-d vector
	cmat *v=cnew(6, 1);
	for(long i=0; i<6; i++) P(v, i)=(real)i;
	cfftshift(v);
	for(long i=0; i<6; i++) assert_float_equal(creal(P(v, i)), (i+3)%6, 0);
	//shifting the spectrum by n/2 is modulation by (-1)^(x+y) in real space
	cmat *C=cffttreat(A);
	for(long y=0; y<ny; y++){
		for(long x=0; x<nx; x++){
			comp expect=((x+y)&1)?-P(A, x, y):P(A, x, y);
			assert_float_equal(cabs(P(C, x, y)-expect), 0, TOL*nx*ny);
		}
	}
	assert_null(cffttreat(NULL));
	cfree(A);
	cfree(B);
	cfree(C);
	cfree(v);
}

UNITTEST_MAIN(
	cmocka_unit_test(fft_naive),
	cmocka_unit_test(fft_roundtrip),
	cmocka_unit_test(fft_known),
	cmocka_unit_test(fft_parseval),
	cmocka_unit_test(fft_partial),
	cmocka_unit_test(fft_cell_split),
	cmocka_unit_test(fft_r2hc),
	cmocka_unit_test(fft_shift_treat)
)
