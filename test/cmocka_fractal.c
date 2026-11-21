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
  Tests of the fractal turbulence screen operators (fractal_do, fractal_inv,
  fractal_trans, fractal_inv_trans) and their statistics compared with theory
  and with FFT based screens. Ported from test_fractal.c, which ran one of
  test_accuracy, test_cov, test_corner, test_part, test_stfun, test_psd or
  test_cxx and wrote the results to files.
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-10
#endif
static const real r0=0.2;
static const int seed=4;
/*Structure function of Kolmogorov turbulence (m^2) at separation r.*/
static real kolmogorov(real r){
	return 6.88*pow(r0, -5./3.)*pow(0.5e-6/(2*M_PI), 2)*pow(r, 5./3.);
}
/*Distance from the origin.*/
static real dist(real x, real y){
	return sqrt(x*x+y*y);
}
/*Fill with random numbers and apply fractal_do*/
static void fractal_rand(dmat *atm, real dx, real L0, long ninit, rand_t *rstat){
	drandn(atm, 1, rstat);
	fractal_do(atm, dx, r0, L0, ninit);
}
/*Build the matrix A of the linear operator fractal_do for an n*n grid by applying it to unit vectors.
  The covariance of fractal screens is A*A'.*/
static dmat *fractal_matrix(long n, real dx, real L0, long ninit){
	dmat *A=dnew(n*n, n*n);
	dmat *atm=dnew(n, n);
	for(long j=0; j<n*n; j++){
		dzero(atm);
		P(atm, j)=1;
		fractal_do(atm, dx, r0, L0, ninit);
		memcpy(PCOL(A, j), P(atm), sizeof(real)*n*n);
	}
	dfree(atm);
	return A;
}
/*Theoretical covariance used by the fractal method for separation r on an n*n grid*/
static real fractal_cov_theory(real r, long n, real dx, real L0){
	dmat *rr=dnew(1, 1);
	P(rr, 0)=r;
	dmat *cov=turbcov(rr, (n-1)*dx*sqrt(2), r0, L0);
	real res=P(cov, 0);
	dfree(rr);
	dfree(cov);
	return res;
}
/*Generate a pair of FFT based periodic screens (real and imaginary part) using spect.*/
static void fft_screens(dmat *atmr, dmat *atmi, cmat *atm, const dmat *spect, rand_t *rstat){
	for(long i=0; i<PN(atm); i++){
		real ra=randn(rstat);
		real ia=randn(rstat);
		P(atm, i)=COMPLEX(ra, ia)*P(spect, i);
	}
	cfft2(atm, -1);
	for(long i=0; i<PN(atm); i++){
		P(atmr, i)=creal(P(atm, i));
		P(atmi, i)=cimag(P(atm, i));
	}
}

/** fractal_inv and fractal_inv_trans invert fractal_do and fractal_trans, and fractal_trans is the adjoint of fractal_do. (test_accuracy) */
static void test_accuracy(void **state){
	(void)state;
	const real dx=1./64.;
	const long N=65;
	const real L0s[2]={INFINITY, 30};
	const long ninits[2]={2, 3};
	rand_t rstat;
	seed_rand(&rstat, seed);
	for(int iL0=0; iL0<2; iL0++){
		for(int ii=0; ii<2; ii++){
			const real L0=L0s[iL0];
			const long ninit=ninits[ii];
			dmat *v=dnew(N, N);
			dmat *u=dnew(N, N);
			drandn(v, 1, &rstat);
			drandn(u, 1, &rstat);
			dmat *Av=ddup(v);
			fractal_do(Av, dx, r0, L0, ninit);
			dmat *Atu=ddup(u);
			fractal_trans(Atu, dx, r0, L0, ninit);
			//adjoint: u'(Av)=v'(A'u)
			real dot1=ddot(u, Av);
			real dot2=ddot(v, Atu);
			assert_float_equal(dot1, dot2, TOL*sqrt(dsumsq(u)*dsumsq(Av)));
			//inverse
			dmat *v2=ddup(Av);
			fractal_inv(v2, dx, r0, L0, ninit);
			assert_dmat_near(v2, v, TOL);
			dmat *u2=ddup(Atu);
			fractal_inv_trans(u2, dx, r0, L0, ninit);
			assert_dmat_near(u2, u, TOL);
			//inverse applied first
			dfree(v2);
			v2=ddup(v);
			fractal_inv(v2, dx, r0, L0, ninit);
			fractal_do(v2, dx, r0, L0, ninit);
			assert_dmat_near(v2, v, TOL);
			//linearity
			dmat *uv=ddup(u);
			dadd(&uv, 2, v, -3);
			fractal_do(uv, dx, r0, L0, ninit);
			dmat *Au=ddup(u);
			fractal_do(Au, dx, r0, L0, ninit);
			dadd(&Au, 2, Av, -3);
			assert_float_equal(ddiff(uv, Au), 0, TOL);
			//Output is an OPD: rms is of the order of r0 based turbulence
			assert_true(sqrt(dsumsq(Av)/PN(Av))<1e-5);
			dfree(uv); dfree(Au);
			dfree(v); dfree(u); dfree(Av); dfree(Atu); dfree(v2); dfree(u2);
		}
	}
}
/** On the initial coarse grid, the fractal screen has the exact theoretical covariance (exact A*A' and Monte Carlo). (test_corner) */
static void test_corner(void **state){
	(void)state;
	const real dx=32;
	const long N=2;
	for(int iL0=0; iL0<2; iL0++){
		const real L0=iL0?100:INFINITY;
		dmat *A=fractal_matrix(N, dx, L0, 2);
		dmat *C=NULL;
		dmm(&C, 0, A, A, "nt", 1);
		const real c0=fractal_cov_theory(0, N, dx, L0);
		for(long j=0; j<N*N; j++){
			for(long i=0; i<N*N; i++){
				real r=dx*dist((i%N)-(j%N), (i/N)-(j/N));
				assert_float_equal(P(C, i, j), fractal_cov_theory(r, N, dx, L0), TOL*c0);
			}
		}
		//Structure function between corners follows Kolmogorov theory
		if(isinf(L0)){
			assert_float_equal(P(C, 0, 0)+P(C, 1, 1)-2*P(C, 0, 1), kolmogorov(dx), TOL*kolmogorov(dx));
			assert_float_equal(P(C, 0, 0)+P(C, 3, 3)-2*P(C, 0, 3), kolmogorov(dx*sqrt(2)), TOL*kolmogorov(dx));
		}
		//Monte Carlo
		rand_t rstat;
		seed_rand(&rstat, seed);
		const long nframe=20000;
		dmat *atm=dnew(N, N);
		dmat *vec=dref_reshape(atm, N*N, 1);
		dmat *cov=NULL;
		for(long i=0; i<nframe; i++){
			fractal_rand(atm, dx, L0, 2, &rstat);
			dmm(&cov, 1, vec, vec, "nt", 1);
		}
		dscale(cov, 1./nframe);
		//statistical error is about sqrt(2/nframe)=1%
		assert_dmat_near(cov, C, 0.05*c0);
		dfree(cov); dfree(vec); dfree(atm);
		dfree(C); dfree(A);
	}
}
/** After one refinement step, the covariance of new points with their parents and their variance are exact. (test_part) */
static void test_part(void **state){
	(void)state;
	const real dx=1./64.;
	const long N=3;
	for(int iL0=0; iL0<2; iL0++){
		const real L0=iL0?30:INFINITY;
		dmat *A=fractal_matrix(N, dx, L0, 2);
		dmat *C=NULL;
		dmm(&C, 0, A, A, "nt", 1);
		const real c0=fractal_cov_theory(0, N, dx, L0);
		const real c1=fractal_cov_theory(dx, N, dx, L0);
		const real c2=fractal_cov_theory(dx*sqrt(2), N, dx, L0);
#define IND(ix,iy) ((ix)+(iy)*N)
		//all points have the same variance
		for(long i=0; i<N*N; i++){
			assert_float_equal(P(C, i, i), c0, TOL*c0);
		}
		//center with 4 corners (square configuration)
		const long corners[4]={IND(0, 0), IND(2, 0), IND(0, 2), IND(2, 2)};
		for(int ic=0; ic<4; ic++){
			assert_float_equal(P(C, IND(1, 1), corners[ic]), c2, TOL*c0);
		}
		//edge middle points with the two ends and the center (triangular configuration)
		assert_float_equal(P(C, IND(1, 0), IND(0, 0)), c1, TOL*c0);
		assert_float_equal(P(C, IND(1, 0), IND(2, 0)), c1, TOL*c0);
		assert_float_equal(P(C, IND(1, 0), IND(1, 1)), c1, TOL*c0);
		assert_float_equal(P(C, IND(0, 1), IND(0, 0)), c1, TOL*c0);
		assert_float_equal(P(C, IND(0, 1), IND(0, 2)), c1, TOL*c0);
		assert_float_equal(P(C, IND(0, 1), IND(1, 1)), c1, TOL*c0);
		assert_float_equal(P(C, IND(2, 1), IND(2, 2)), c1, TOL*c0);
		assert_float_equal(P(C, IND(1, 2), IND(0, 2)), c1, TOL*c0);
#undef IND
		//Monte Carlo of the 4 points in the lower left corner as in the old test
		rand_t rstat;
		seed_rand(&rstat, seed);
		const long nframe=20000;
		dmat *atm=dnew(N, N);
		dmat *vec=dnew(4, 1);
		dmat *cov=NULL;
		for(long i=0; i<nframe; i++){
			fractal_rand(atm, dx, L0, 2, &rstat);
			P(vec, 0)=P(atm, 0, 0);
			P(vec, 1)=P(atm, 1, 0);
			P(vec, 2)=P(atm, 0, 1);
			P(vec, 3)=P(atm, 1, 1);
			dmm(&cov, 1, vec, vec, "nt", 1);
		}
		dscale(cov, 1./nframe);
		const long ind[4]={0, 1, N, N+1};
		for(int j=0; j<4; j++){
			for(int i=0; i<4; i++){
				assert_float_equal(P(cov, i, j), P(C, ind[i], ind[j]), 0.05*c0);
			}
		}
		dfree(cov); dfree(vec); dfree(atm);
		dfree(C); dfree(A);
	}
}
/** Covariance of fractal screens is close to the theory: exact at the coarse grid and approximately elsewhere. (test_cov) */
static void test_cov(void **state){
	(void)state;
	const real dx=1./4.;
	const long N=17;
	const real L0=INFINITY;
	dmat *A=fractal_matrix(N, dx, L0, 2);
	dmat *C=NULL;
	dmm(&C, 0, A, A, "nt", 1);
	//Compare structure function with Kolmogorov theory for all pairs.
	real maxerr=0, sumerr=0;
	long count=0;
	for(long j=0; j<N*N; j++){
		for(long i=0; i<j; i++){
			real r=dx*dist((i%N)-(j%N), (i/N)-(j/N));
			real D=P(C, i, i)+P(C, j, j)-2*P(C, i, j);
			real err=D/kolmogorov(r)-1;
			if(fabs(err)>maxerr) maxerr=fabs(err);
			sumerr+=err;
			count++;
		}
	}
	//fractal method lacks some high frequency power, but agrees with theory reasonably well.
	assert_true(maxerr<0.3);
	assert_float_equal(sumerr/count, 0, 0.1);
	//the corners are on the initial grid and are exact.
	const long ic[4]={0, N-1, (N-1)*N, N*N-1};
	for(int j=0; j<4; j++){
		for(int i=0; i<j; i++){
			real r=dx*(N-1)*dist((ic[i]%N)!=(ic[j]%N), (ic[i]/N)!=(ic[j]/N));
			real D=P(C, ic[i], ic[i])+P(C, ic[j], ic[j])-2*P(C, ic[i], ic[j]);
			assert_float_equal(D, kolmogorov(r), TOL*D);
		}
	}
	//Covariance is symmetric positive definite
	dmat *U=NULL, *S=NULL, *VT=NULL;
	dsvd(&U, &S, &VT, C);
	assert_true(dmin(S)>0);
	dfree(U); dfree(S); dfree(VT);
	dfree(C);
	dfree(A);
}
/** Structure function computed by stfun from fractal and FFT screens agrees with Kolmogorov theory at small separations. (test_stfun) */
static void test_stfun(void **state){
	(void)state;
	const real dx=1./16;
	const long N=32;
	const long nframe=200;
	rand_t rstat;
	seed_rand(&rstat, seed);
	//stfun of a deterministic tilt is r^2*slope^2.
	{
		dmat *opd=dnew(N, N);
		for(long iy=0; iy<N; iy++){
			for(long ix=0; ix<N; ix++){
				P(opd, ix, iy)=0.3*ix*dx-0.4*iy*dx;
			}
		}
		dmat *st=stfun_batch(opd, NULL);
		assert_int_equal(NX(st), 2*N);
		for(long iy=-5; iy<=5; iy++){
			for(long ix=-5; ix<=5; ix++){
				real expect=pow((0.3*ix-0.4*iy)*dx, 2);
				assert_float_equal(P(st, N+ix, N+iy), expect, 1e-10);
			}
		}
		dfree(st);
		dfree(opd);
	}
	dmat *st[2]={NULL, NULL};
	{
		map_t *atm=mapnew(N+1, N+1, dx, dx);
		stfun_t *data=stfun_init(N, N, NULL);
		for(long i=0; i<nframe; i++){
			fractal_rand(atm->dmat, dx, INFINITY, 2, &rstat);
			stfun_push(data, atm->dmat);
		}
		st[0]=stfun_finalize(data);
		mapfree(atm);
	}
	{
		stfun_t *data=stfun_init(N, N, NULL);
		dmat *spect=turbpsd(N, N, dx, r0, 100, 0, 0.5);
		cmat *atm=cnew(N, N);
		dmat *atmr=dnew(N, N);
		dmat *atmi=dnew(N, N);
		for(long ii=0; ii<nframe; ii+=2){
			fft_screens(atmr, atmi, atm, spect, &rstat);
			stfun_push(data, atmr);
			stfun_push(data, atmi);
		}
		st[1]=stfun_finalize(data);
		dfree(spect); cfree(atm); dfree(atmr); dfree(atmi);
	}
	//Exact structure function of the periodic FFT screens: 2*(R(0)-R(r)) with R=Re(FFT(spect^2)).
	dmat *spect=turbpsd(N, N, dx, r0, 100, 0, 0.5);
	cmat *cov=cnew(N, N);
	for(long i=0; i<PN(spect); i++){
		P(cov, i)=pow(P(spect, i), 2);
	}
	cfft2(cov, -1);
	for(long ir=1; ir<=4; ir++){
		const real Dt=kolmogorov(ir*dx);
		//fractal screen: close to theory, but lacks some power at small separations.
		const real Fx=P(st[0], N+ir, N)/Dt;
		const real Fy=P(st[0], N, N+ir)/Dt;
		assert_true(Fx>0.7 && Fx<1.05);
		assert_true(Fy>0.7 && Fy<1.05);
		assert_float_equal(Fx, Fy, 0.05);//isotropy
		//FFT screen: agrees with its own spectrum (statistical error of a few percent)
		const real Dfft=2*creal(P(cov, 0, 0)-P(cov, ir, 0));
		assert_float_equal(P(st[1], N+ir, N)/Dfft, 1, 0.05);
		assert_float_equal(P(st[1], N, N+ir)/Dfft, 1, 0.05);
		//but a 2 m periodic screen without subharmonics lacks power compared to Kolmogorov theory
		assert_true(Dfft<0.7*Dt);
		assert_true(Fx>P(st[1], N+ir, N)/Dt);
	}
	cfree(cov);
	dfree(spect);
	//FFT screens lack large scale power relative to the fractal screens.
	assert_true(P(st[1], N+N/2, N)<P(st[0], N+N/2, N));
	dfree(st[0]);
	dfree(st[1]);
}
/*Copy the N*N lower left corner of screen to hat after removing piston, tip and
  tilt and applying a Hann window to suppress the spectral leakage of a non-periodic screen.*/
static void detrend_window(cmat *hat, const dmat *screen){
	const long N=NX(hat);
	real mean=0, tx=0, ty=0, norm=0;
	for(long iy=0; iy<N; iy++){
		for(long ix=0; ix<N; ix++){
			real x=ix-(N-1)*0.5, y=iy-(N-1)*0.5;
			mean+=P(screen, ix, iy);
			tx+=P(screen, ix, iy)*x;
			ty+=P(screen, ix, iy)*y;
			norm+=x*x;
		}
	}
	mean/=N*N; tx/=norm; ty/=norm;//x and y are orthogonal on the full grid
	for(long iy=0; iy<N; iy++){
		real wy=0.5-0.5*cos(2*M_PI*(iy+0.5)/N);
		for(long ix=0; ix<N; ix++){
			real wx=0.5-0.5*cos(2*M_PI*(ix+0.5)/N);
			real x=ix-(N-1)*0.5, y=iy-(N-1)*0.5;
			P(hat, ix, iy)=(P(screen, ix, iy)-mean-tx*x-ty*y)*wx*wy;
		}
	}
}
/** PSD of FFT screens equals the input spectrum, and fractal screen PSD is comparable to that of FFT screens at the same size. (test_psd) */
static void test_psd(void **state){
	(void)state;
	const real dx=1./16;
	const long N=32;
	const long nframe=64;
	rand_t rstat;
	seed_rand(&rstat, seed);
	dmat *spect=turbpsd(N, N, dx, r0, 100, 0, 0.5);
	//FFT based screen: |FFT(screen)|^2=N^4*spect^2 on average.
	{
		cmat *atm=cnew(N, N);
		cmat *hat=cnew(N, N);
		dmat *atmr=dnew(N, N);
		dmat *atmi=dnew(N, N);
		dmat *psd=dnew(N, N);
		for(long ii=0; ii<nframe; ii+=2){
			fft_screens(atmr, atmi, atm, spect, &rstat);
			for(int ir=0; ir<2; ir++){
				cembedd(hat, ir?atmi:atmr, 0);
				cfft2(hat, -1);
				for(long i=0; i<PN(hat); i++){
					P(psd, i)+=pow(cabs(P(hat, i)), 2);
				}
			}
		}
		real ratio=0;
		for(long i=1; i<PN(psd); i++){
			ratio+=P(psd, i)/(nframe*pow(N, 4)*pow(P(spect, i), 2));
		}
		ratio/=(PN(psd)-1);
		assert_float_equal(ratio, 1, 0.03);
		assert_float_equal(P(psd, 0), 0, 1e-20);//no piston
		cfree(atm); cfree(hat); dfree(atmr); dfree(atmi); dfree(psd);
	}
	//Fractal screen versus FFT screen (both of size 2N, cut to N) in radial bins.
	{
		const long nbin=8;
		dmat *bins[2]={dnew(nbin, 1), dnew(nbin, 1)};
		map_t *frac=mapnew(2*N+1, 2*N+1, dx, dx);
		dmat *spect2=turbpsd(2*N, 2*N, dx, r0, 100, 0, 0.5);
		cmat *atm=cnew(2*N, 2*N);
		dmat *atmr=dnew(2*N, 2*N);
		dmat *atmi=dnew(2*N, 2*N);
		cmat *hat=cnew(N, N);
		for(long ii=0; ii<nframe; ii++){
			for(int is=0; is<2; is++){
				dmat *screen=NULL;
				if(is==0){
					fractal_rand(frac->dmat, dx, INFINITY, 2, &rstat);
					screen=frac->dmat;
				}else{
					if(ii%2==0) fft_screens(atmr, atmi, atm, spect2, &rstat);
					screen=ii%2==0?atmr:atmi;
				}
				detrend_window(hat, screen);
				cfft2(hat, -1);
				for(long iy=0; iy<N; iy++){
					long fy=iy<N/2?iy:N-iy;
					for(long ix=0; ix<N; ix++){
						long fx=ix<N/2?ix:N-ix;
						long ib=(long)(sqrt(fx*fx+fy*fy)/(N/2)*nbin);
						if(ib>0&&ib<nbin){
							P(bins[is], ib)+=pow(cabs(P(hat, ix, iy)), 2);
						}
					}
				}
			}
		}
		for(long ib=1; ib<nbin; ib++){
			real ratio=P(bins[0], ib)/P(bins[1], ib);
			assert_true(ratio>0.4&&ratio<3);
		}
		//Relative to FFT screens which have no power above the Nyquist frequency, fractal screens
		//contain more high frequency power (aliased) after removing tip/tilt.
		assert_true(P(bins[0], nbin-1)/P(bins[1], nbin-1)>P(bins[0], 1)/P(bins[1], 1));
		dfree(bins[0]); dfree(bins[1]);
		mapfree(frac); dfree(spect2); cfree(atm); dfree(atmr); dfree(atmi); cfree(hat);
	}
	dfree(spect);
}
/** Sample covariance of fractal and FFT screens agree with the exact covariance, and both are consistent with Kolmogorov structure function. (test_cxx) */
static void test_cxx(void **state){
	(void)state;
	const real dx=1./4;
	const long N=8;
	const long nframe=4000;
	rand_t rstat;
	seed_rand(&rstat, seed);
	loc_t *loc=mksqloc_auto(N, N, dx, dx);
	dmat *B=stfun_kolmogorov(loc, r0);
	real Dfrac=0;//ratio of fractal screen structure function to theory at the smallest separation
	//fractal screens
	{
		dmat *A=fractal_matrix(N+1, dx, INFINITY, 2);
		dmat *Cfull=NULL;
		dmm(&Cfull, 0, A, A, "nt", 1);
		//extract the N*N subgrid covariance
		dmat *Cex=dnew(N*N, N*N);
		for(long j=0; j<N*N; j++){
			long jj=(j%N)+(j/N)*(N+1);
			for(long i=0; i<N*N; i++){
				long ii=(i%N)+(i/N)*(N+1);
				P(Cex, i, j)=P(Cfull, ii, jj);
			}
		}
		dmat *cxx=dnew(N*N, N*N);
		map_t *atm=mapnew(N+1, N+1, dx, dx);
		dmat *atmvec=dnew(N*N, 1);
		for(long i=0; i<nframe; i++){
			fractal_rand(atm->dmat, dx, INFINITY, 2, &rstat);
			for(long k=0; k<N*N; k++){
				P(atmvec, k)=P(atm, k%N, k/N);
			}
			dmm(&cxx, 1, atmvec, atmvec, "nt", 1);
		}
		dscale(cxx, 1./nframe);
		const real c0=P(Cex, 0, 0);
		assert_dmat_near(cxx, Cex, 0.1*c0);
		//structure function vs theory
		Dfrac=(P(Cex, 0, 0)+P(Cex, 1, 1)-2*P(Cex, 0, 1))/P(B, 0, 1);
		for(long j=0; j<N*N; j++){
			for(long i=0; i<j; i++){
				real D=P(Cex, i, i)+P(Cex, j, j)-2*P(Cex, i, j);
				assert_float_equal(D/P(B, i, j), 1, 0.3);
			}
		}
		dfree(atmvec); mapfree(atm); dfree(cxx); dfree(Cex); dfree(Cfull); dfree(A);
	}
	//FFT screens
	{
		dmat *spect=turbpsd(N, N, dx, r0, 100, 0, 0.5);
		P(spect, 0)=P(spect, 1);//keep some piston power as in the old test
		//Exact covariance: Re(FFT(spect^2)) at the separation (periodic)
		cmat *cov=cnew(N, N);
		for(long i=0; i<PN(spect); i++){
			P(cov, i)=pow(P(spect, i), 2);
		}
		cfft2(cov, -1);
		dmat *Cex=dnew(N*N, N*N);
		for(long j=0; j<N*N; j++){
			for(long i=0; i<N*N; i++){
				long sx=((i%N)-(j%N)+N)%N;
				long sy=((i/N)-(j/N)+N)%N;
				P(Cex, i, j)=creal(P(cov, sx, sy));
			}
		}
		dmat *cxx=dnew(N*N, N*N);
		cmat *atm=cnew(N, N);
		dmat *atmr=dnew(N*N, 1);
		dmat *atmi=dnew(N*N, 1);
		for(long ii=0; ii<nframe; ii+=2){
			fft_screens(atmr, atmi, atm, spect, &rstat);
			dmm(&cxx, 1, atmr, atmr, "nt", 1);
			dmm(&cxx, 1, atmi, atmi, "nt", 1);
		}
		dscale(cxx, 1./nframe);
		const real c0=P(Cex, 0, 0);
		assert_dmat_near(cxx, Cex, 0.1*c0);
		//structure function at the smallest separation is close to theory
		real D=P(Cex, 0, 0)+P(Cex, 1, 1)-2*P(Cex, 0, 1);
		assert_true(D/P(B, 0, 1)<0.6);//small periodic FFT screen lacks power
		assert_true(D/P(B, 0, 1)<Dfrac);//fractal screen is closer to theory
		assert_true(D/P(B, 0, 1)>0.2);
		dfree(cxx); cfree(atm); dfree(atmr); dfree(atmi); dfree(Cex); cfree(cov); dfree(spect);
	}
	dfree(B);
	locfree(loc);
}

UNITTEST_MAIN(
	cmocka_unit_test(test_accuracy),
	cmocka_unit_test(test_corner),
	cmocka_unit_test(test_part),
	cmocka_unit_test(test_cov),
	cmocka_unit_test(test_stfun),
	cmocka_unit_test(test_psd),
	cmocka_unit_test(test_cxx)
)
