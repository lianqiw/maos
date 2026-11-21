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
#include "unittest.h"
/*
  Tests of the random number generators (port of test_rand.c). The statistical
  tolerances are 6 sigma of the sample estimate so that a correct generator
  essentially never fails, while a broken one (wrong scaling, bias) does.
*/
#define NSAMPLE 200000

/** The same seed reproduces the same sequence for every generator; different seeds differ. */
static void rand_reproducible(void **state){
	(void)state;
	rand_t a, b;
	seed_rand(&a, 1);
	seed_rand(&b, 1);
	for(int i=0; i<1000; i++){
		assert_true(lrand(&a)==lrand(&b));
		assert_true(randu(&a)==randu(&b));
		assert_true(randn(&a)==randn(&b));
		assert_int_equal(randp(&a, 30), randp(&b, 30));
	}
	//lrandu combines two consecutive 32 bit values, also across state refreshes.
	seed_rand(&a, 1);
	seed_rand(&b, 1);
	for(int i=0; i<5000; i++){
		real v=lrandu(&a);
		assert_true(v>=0&&v<1);
		real u1=randu(&b);
		real u2=randu(&b);
		assert_float_equal(v, u1+u2*(1./4294967296.), 1e-15);
	}
	//Reseeding restarts the sequence (test_number in the old test).
	seed_rand(&a, 1);
	real first=randn(&a);
	seed_rand(&a, 1);
	assert_true(first==randn(&a));
	//Seed a new stream with an integer drawn from another one.
	seed_rand(&a, 1);
	unsigned long i1=lrand(&a);
	unsigned long i2=lrand(&a);
	assert_true(i1!=i2);
	seed_rand(&b, i2);
	rand_t c;
	seed_rand(&c, i2);
	assert_true(randn(&b)==randn(&c));
	//Different seeds give different sequences.
	seed_rand(&a, 1);
	seed_rand(&b, 2);
	int ndiff=0;
	for(int i=0; i<100; i++){
		if(lrand(&a)!=lrand(&b)) ndiff++;
	}
	assert_true(ndiff>95);
	//The matrix fillers use the same stream as the scalar ones.
	seed_rand(&a, 3);
	seed_rand(&b, 3);
	dmat *A=dnew(5, 4);
	drandu(A, 2, &a);
	for(long i=0; i<PN(A); i++){
		assert_float_equal(P(A, i), 2*randu(&b), 1e-6);
	}
	drandn(A, 3, &a);
	for(long i=0; i<PN(A); i++){
		assert_float_equal(P(A, i), 3*randn(&b), 1e-5);
	}
	dfree(A);
}
/** lrand returns 32 bit integers with uniformly distributed low and high bits. */
static void rand_lrand(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	long nodd=0, nhigh=0;
	for(int i=0; i<NSAMPLE; i++){
		unsigned long v=lrand(&rstat);
		assert_true(v<=0xffffffffUL);
		nodd+=v&1;
		nhigh+=(v>>31)&1;
	}
	real tol=6*sqrt(0.25/NSAMPLE);
	assert_float_equal((real)nodd/NSAMPLE, 0.5, tol);
	assert_float_equal((real)nhigh/NSAMPLE, 0.5, tol);
}
/**
   Check that x has the given mean and variance (6 sigma tolerance for N
   samples given the fourth central moment m4).
*/
static void check_mean_var(const real *x, long n, real mean, real var, real m4){
	real sum=0, sum2=0;
	for(long i=0; i<n; i++) sum+=x[i];
	real m=sum/n;
	for(long i=0; i<n; i++) sum2+=(x[i]-m)*(x[i]-m);
	real v=sum2/(n-1);
	assert_float_equal(m, mean, 6*sqrt(var/n));
	assert_float_equal(v, var, 6*sqrt((m4-var*var)/n));
}
/** randu is in [0,1) with mean 1/2, variance 1/12, and a flat histogram. */
static void rand_randu(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	dmat *p=dnew(NSAMPLE, 1);
	{
		const int nbin=10;
		long hist[10]={0};
		for(long i=0; i<NSAMPLE; i++){
			real v=randu(&rstat);
			assert_true(v>=0&&v<1);
			P(p, i)=v;
			hist[(int)(v*nbin)]++;
		}
		check_mean_var(P(p), NSAMPLE, 0.5, 1./12., 1./80.);
		//Each bin count is binomial(N, 0.1).
		real ex=(real)NSAMPLE/nbin;
		for(int ib=0; ib<nbin; ib++){
			assert_float_equal(hist[ib], ex, 6*sqrt(ex*0.9));
		}
		//Adjacent samples are uncorrelated: E[(x_i-0.5)(x_{i+1}-0.5)]=0 with std (1/12)/sqrt(N).
		real cor=0;
		for(long i=0; i+1<NSAMPLE; i++){
			cor+=(P(p, i)-0.5)*(P(p, i+1)-0.5);
		}
		cor/=NSAMPLE-1;
		assert_float_equal(cor, 0, 6./12./sqrt(NSAMPLE));
	}
	//drandu scales by max
	drandu(p, 10, &rstat);
	assert_true(dmin(p)>=0);
	assert_true(dmax(p)<10);
	check_mean_var(P(p), NSAMPLE, 5, 100./12., 10000./80.);
	dfree(p);
	//NULL input is ignored.
	drandu(NULL, 1, &rstat);
}
/** randn is standard normal: mean, variance, kurtosis, and the 1 and 3 sigma tail fractions. */
static void rand_randn(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	dmat *p=dnew(NSAMPLE, 1);
	drandn(p, 1, &rstat);
	check_mean_var(P(p), NSAMPLE, 0, 1, 3);
	long n1=0, n3=0;
	real m3=0, m4=0;
	for(long i=0; i<NSAMPLE; i++){
		real v=P(p, i);
		if(fabs(v)<1) n1++;
		if(fabs(v)>3) n3++;
		m3+=v*v*v;
		m4+=v*v*v*v;
	}
	m3/=NSAMPLE; m4/=NSAMPLE;
	assert_float_equal(m3, 0, 6*sqrt(15./NSAMPLE));//var of x^3 is 15
	assert_float_equal(m4, 3, 6*sqrt(96./NSAMPLE));//var of x^4 is 105-9
	real p1=erf(1/sqrt(2.));
	real p3=erfc(3/sqrt(2.));
	assert_float_equal((real)n1/NSAMPLE, p1, 6*sqrt(p1*(1-p1)/NSAMPLE));
	assert_float_equal((real)n3/NSAMPLE, p3, 6*sqrt(p3*(1-p3)/NSAMPLE));
	//Symmetry: as many positive as negative.
	long npos=0;
	for(long i=0; i<NSAMPLE; i++) npos+=P(p, i)>0;
	assert_float_equal((real)npos/NSAMPLE, 0.5, 6*sqrt(0.25/NSAMPLE));
	//drandn scales by sigma
	drandn(p, 10, &rstat);
	check_mean_var(P(p), NSAMPLE, 0, 100, 3e4);
	dfree(p);
	drandn(NULL, 1, &rstat);
}
/** randp is Poisson: mean and variance equal xm for the direct (xm<=200) and Gaussian (xm>200) paths. */
static void rand_randp(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	const long n=NSAMPLE/4;
	dmat *p=dnew(n, 1);
	const real xms[]={0.5, 3, 30, 1000};
	for(unsigned int im=0; im<sizeof(xms)/sizeof(xms[0]); im++){
		real xm=xms[im];
		long nzero=0;
		for(long i=0; i<n; i++){
			long v=randp(&rstat, xm);
			assert_true(v>=0);
			if(!v) nzero++;
			P(p, i)=v;
		}
		//Poisson: variance=xm, fourth central moment=xm(1+3xm)
		check_mean_var(P(p), n, xm, xm, xm*(1+3*xm));
		if(xm<=30){
			real p0=exp(-xm);
			assert_float_equal((real)nzero/n, p0, 6*sqrt(p0*(1-p0)/n)+1./n);
		}
	}
	//Zero or negative mean gives zero.
	assert_int_equal(randp(&rstat, 0), 0);
	assert_int_equal(randp(&rstat, -1), 0);
	dfree(p);
}

UNITTEST_MAIN(
	cmocka_unit_test(rand_reproducible),
	cmocka_unit_test(rand_lrand),
	cmocka_unit_test(rand_randu),
	cmocka_unit_test(rand_randn),
	cmocka_unit_test(rand_randp),
)
