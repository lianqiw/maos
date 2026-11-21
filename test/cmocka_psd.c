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

void test_psd2d(void **state){
	(void)state;
	dmat *screen=(dmat *)genatm_simple(0.186, 0, -11./3., 1./64., 128, 1);
	dmat *psd=psd2d_aniso(screen, 1./64.);
	real inte1=sqrt(psd_inte2(psd))*1e9;
	//writebin(psd, "psd_%ld", nx);
	dfree(psd);
	dmat *extra=NULL;
	psd=psd2d(&extra, screen, 1./64.);
	real inte2=sqrt(psd_inte2(psd))*1e9;
	dfree(psd);
	real rms1=sqrt(dsumsq(screen)/PN(screen))*1e9;
	info("size is %ld rms is %g, psdint is %g, psdint 2 is %g, recon r0=%g, slope=%g\n",
		NX(screen), rms1, inte1, inte2, P(extra, 0), P(extra, 1));
	assert_float_equal(rms1, 169.648, 1e-3);
	assert_float_equal(inte1, 56.4647, 1e-3);
	assert_float_equal(inte2, 216.578, 1e-3);
	assert_float_equal(P(extra,0), 0.128266, 1e-3);
	assert_float_equal(P(extra, 1), -3.68282, 1e-3);
	dfree(extra);
	dfree(screen);
}

void test_stfun(void **state){
	(void)state;
	dmat *screen=(dmat *)genatm_simple(0.186, 0, -11./3., 1./64., 128, 1);
	dmat *st=stfun_batch(screen, NULL);
	info("sum=%g, max=%g\n", dsum(st), dmax(st));
	assert_float_equal(dsum(st)*1e9, 3.01473, 1e-3);
	assert_float_equal(dmax(st)*1e13, 1.15708, 1e-3);
	//dshow(st, "st");
	dfree(st);
	dfree(screen);
}
/*Two column PSD [f, psd] on a given frequency grid.*/
static dmat *psd_table(const dmat *f, real (*fun)(real)){
	dmat *psd=dnew(NX(f), 2);
	for(long i=0; i<NX(f); i++){
		P(psd, i, 0)=P(f, i);
		P(psd, i, 1)=fun(P(f, i));
	}
	return psd;
}
static real fun_exp(real f){ return exp(-f*0.1); }
static real fun_lin(real f){ return 2+3*f; }
static real fun_pow(real f){ return 1e-6*pow(f, -2); }
static real fun_one(real f){ (void)f; return 1; }

/** psd1d satisfies Parseval's theorem (sum of PSD equals the mean square,
	averaged over segments) and psd1dt puts frequency in the first column. */
static void test_psd1d(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long n=1024;
	dmat *v=dnew(n, 2);
	drandn(v, 1, &rstat);
	dmat *psd=psd1d(v, 1);
	assert_int_equal(NX(psd), n/2+1);
	assert_int_equal(NY(psd), 2);
	for(int ic=0; ic<2; ic++){
		real ms=0;
		for(long i=0; i<n; i++) ms+=pow(P(v, i, ic), 2);
		ms/=n;
		real sum=0;
		for(long i=0; i<NX(psd); i++) sum+=P(psd, i, ic);
		assert_float_equal(sum, ms, 1e-10*ms);
	}
	dfree(psd);
	//3 segments of half overlapping windows of length 2*n/4
	psd=psd1d(v, 3);
	assert_int_equal(NX(psd), n/4+1);
	real ms=0;
	for(int iseg=0; iseg<3; iseg++){
		for(long i=0; i<n/2; i++) ms+=pow(P(v, iseg*n/4+i, 0), 2);
	}
	ms/=(3*n/2);
	real sum=0;
	for(long i=0; i<NX(psd); i++) sum+=P(psd, i, 0);
	assert_float_equal(sum, ms, 1e-10*ms);
	dfree(psd);
	//row vector input is treated as a single sequence. nseg<1 is treated as 1.
	dmat *vr=dnew(1, n);
	for(long i=0; i<n; i++) P(vr, i)=P(v, i, 0);
	psd=psd1d(vr, 0);
	dmat *psd0=psd1d(v, 1);
	assert_int_equal(NY(psd), 1);
	for(long i=0; i<NX(psd); i++){
		assert_float_equal(P(psd, i), P(psd0, i, 0), 1e-12);
	}
	dfree(psd0);
	dfree(psd);
	//psd1dt: white noise of unit variance has one sided PSD density 2*dt
	const real dt=1e-3;
	dmat *psdt=psd1dt(v, 4, dt);
	assert_int_equal(NY(psdt), 3);
	long nf=NX(psdt);
	real df=P(psdt, 1, 0);
	assert_float_equal(df, 1./((nf-1)*2*dt), 1e-12);
	assert_float_equal(P(psdt, nf-1, 0), 0.5/dt, 1e-9);
	real mean=0;
	for(long i=1; i<nf-1; i++) mean+=P(psdt, i, 1)+P(psdt, i, 2);
	mean/=2*(nf-2);
	assert_float_equal(mean, 2*dt, 0.1*2*dt);
	dfree(psdt);
	dfree(vr);
	dfree(v);
}

/** psd_interp1 is exact for functions linear in the interpolation space and
	returns a tiny value outside of the range. */
static void test_psd_interp1(void **state){
	(void)state;
	dmat *f=dlinspace(1, 1, 100);
	dmat *fnew=dlinspace(1.25, 0.5, 150);//1.25 to 75.75
	dmat *fout=dlinspace(200, 1, 2);
	dmat *psd=psd_table(f, fun_lin);
	dmat *res=psd_interp1(psd, fnew, 0);
	assert_int_equal(NX(res), 150);
	for(long i=0; i<NX(fnew); i++){
		assert_float_equal(P(res, i), fun_lin(P(fnew, i)), 1e-10);
	}
	dfree(res);
	dfree(psd);
	psd=psd_table(f, fun_exp);
	res=psd_interp1(psd, fnew, 1);
	for(long i=0; i<NX(fnew); i++){
		assert_float_equal(P(res, i), fun_exp(P(fnew, i)), 1e-10*fun_exp(P(fnew, i)));
	}
	dfree(res);
	res=psd_interp1(psd, fout, 1);
	assert_float_equal(P(res, 0), 1e-40, 1e-50);
	dfree(res);
	res=psd_interp1(psd, fout, 0);
	assert_float_equal(P(res, 1), 1e-40, 1e-50);
	dfree(res);
	dfree(psd);
	dfree(fout);
	dfree(fnew);
	dfree(f);
}

/** psd_vibid finds the vibration peaks injected on a slightly rippled flat background. */
static void test_psd_vibid(void **state){
	(void)state;
	/*psd_vibid() starts its second order low pass filter from 0, only detects
	  peaks above 1 Hz, and uses the first difference as the initial slope
	  estimate. It needs a fine grid below 1 Hz and some ripple in the
	  background to work properly.*/
	const long nf=10000;
	dmat *psd=dnew(nf, 2);//0.01 to 100 Hz
	const real fv[2]={23.0, 87.0};
	const real amp=1e-15, wd=0.3;
	for(int ipeak=0; ipeak<2; ipeak++){
		for(long i=0; i<nf; i++){
			P(psd, i, 0)=0.01*(i+1);
			P(psd, i, 1)=1e-20*(1+0.1*sin(i*1.7));
			if(ipeak){
				for(int iv=0; iv<2; iv++){
					real x=(P(psd, i, 0)-fv[iv])/wd;
					P(psd, i, 1)+=amp*exp(-x*x);
				}
			}
		}
		dmat *vib=psd_vibid(psd);
		assert_int_equal(NX(vib), 4);
		if(!ipeak){//no peak on the background
			assert_int_equal(NY(vib), 0);
		}else{
			assert_int_equal(NY(vib), 2);
			for(int iv=0; iv<2; iv++){
				assert_float_equal(P(vib, 0, iv), fv[iv], 0.01);//centroid
				//integrated power of the gaussian peak
				assert_float_equal(P(vib, 1, iv), amp*wd*sqrt(M_PI), 0.01*amp*wd*sqrt(M_PI));
				//start and end index bracket the peak
				assert_true(P(psd, (long)P(vib, 2, iv), 0)<fv[iv]-wd);
				assert_true(P(psd, (long)P(vib, 3, iv), 0)>fv[iv]+wd);
			}
		}
		dfree(vib);
	}
	//a peak that extends to the end of the frequency range is reported
	for(long i=0; i<nf; i++){
		P(psd, i, 1)=1e-20*(1+0.1*sin(i*1.7));
		real x=(P(psd, i, 0)-100)/wd;
		P(psd, i, 1)+=amp*exp(-x*x);
	}
	dmat *vib=psd_vibid(psd);
	assert_int_equal(NY(vib), 1);
	assert_int_equal((long)P(vib, 3, 0), nf-1);
	assert_true(P(vib, 0, 0)>99);
	dfree(vib);
	//too short input returns an empty result without out of range access
	dmat *psd2=dsub(psd, 0, 1, 0, 2);
	vib=psd_vibid(psd2);
	assert_int_equal(PN(vib), 0);
	dfree(vib);
	dfree(psd2);
	dfree(psd);
}

/** psd_t2s and psd_s2t are inverse of each other and transform the power law as expected. */
static void test_psd_t2s(void **state){
	(void)state;
	dmat *psdt=dnew(1, 4);
	P(psdt, 0)=-8./3.; P(psdt, 1)=1e-6; P(psdt, 2)=0.1; P(psdt, 3)=100;
	const real v=12;
	dmat *psds=psd_t2s(psdt, v);
	assert_int_equal(NY(psds), 4);
	assert_float_equal(P(psds, 0), -11./3., 1e-12);
	assert_float_equal(P(psds, 2), 0.1/v, 1e-12);
	assert_float_equal(P(psds, 3), 100/v, 1e-12);
	//beta2=beta/B(1/2,(n-1)/2)*v^(2+alpha2) with n=11/3
	real n=11./3.;
	real bfun=tgamma(0.5)*tgamma((n-1)*0.5)/tgamma(n*0.5);
	assert_float_equal(P(psds, 1), 1e-6/bfun*pow(v, 2-n), 1e-12*P(psds, 1));
	dmat *psdt2=psd_s2t(psds, v);
	assert_dmat_near(psdt2, psdt, 1e-12*1e-6+1e-12);
	assert_float_equal(P(psdt2, 1), 1e-6, 1e-15);
	dfree(psdt2);
	dfree(psds);
	//two element form
	dmat *ps=dnew(1, 2);
	P(ps, 0)=-11./3.; P(ps, 1)=2;
	dmat *pt=psd_s2t(ps, 5);
	assert_int_equal(NY(pt), 2);
	dmat *ps2=psd_t2s(pt, 5);
	assert_dmat_near(ps2, ps, 1e-12);
	dfree(ps2);
	dfree(pt);
	dfree(ps);
	dfree(psdt);
}

/** psd_inte and psd_inte2 integrate on linear and logarithmic grids. */
static void test_psd_inte(void **state){
	(void)state;
	dmat *f=dlinspace(1, 0.5, 101);//1 to 51
	dmat *psd=psd_table(f, fun_one);
	//rectangle rule on linear grid: n*dnu
	assert_float_equal(psd_inte2(psd), 101*0.5, 1e-10);
	assert_float_equal(psd_inte(P(f), PCOL(psd, 1), NX(f)), 101*0.5, 1e-10);
	dfree(psd);
	dfree(f);
	//log grid: sum psd*nu*dlognu ~ integral
	f=dlogspace(-1, 2, 3001);
	psd=psd_table(f, fun_pow);
	real ref=1e-6*(1./0.1-1./100);
	assert_float_equal(psd_inte2(psd), ref, 2e-3*ref);
	dfree(psd);
	//irregular grid returns 0
	P(f, 1)=P(f, 0)*1.5;
	psd=psd_table(f, fun_one);
	assert_float_equal(psd_inte2(psd), 0, 0);
	dfree(psd);
	dfree(f);
}

/** psd_reverse_cumu: first entry is sqrt of the total integral, monotonic decreasing. */
static void test_psd_reverse_cumu(void **state){
	(void)state;
	//linear grid with flat PSD: cumu(i)=sqrt((n-i)*df)
	dmat *f=dlinspace(1, 2, 50);
	dmat *psd=dnew(50, 3);
	for(long i=0; i<50; i++){
		P(psd, i, 0)=P(f, i);
		P(psd, i, 1)=1;
		P(psd, i, 2)=4;
	}
	dmat *cumu=psd_reverse_cumu(psd, 10);
	assert_int_equal(NX(cumu), 50);
	assert_int_equal(NY(cumu), 3);
	for(long i=0; i<50; i++){
		assert_float_equal(P(cumu, i, 0), P(f, i), 0);
		assert_float_equal(P(cumu, i, 1), 10*sqrt((50-i)*2.), 1e-10);
		assert_float_equal(P(cumu, i, 2), 2*P(cumu, i, 1), 1e-10);
	}
	dfree(cumu);
	dfree(psd);
	dfree(f);
	//log grid
	f=dlogspace(-1, 2, 1001);
	psd=psd_table(f, fun_pow);
	cumu=psd_reverse_cumu(psd, 1);
	real tot=psd_inte2(psd);
	assert_float_equal(P(cumu, 0, 1), sqrt(tot), 1e-10*sqrt(tot));
	for(long i=1; i<NX(cumu); i++){
		assert_true(P(cumu, i, 1)<P(cumu, i-1, 1));
	}
	dfree(cumu);
	dfree(psd);
	dfree(f);
}

/** psd2ts produces time series whose variance matches the integral of the PSD
	for both the power law and the tabulated form. psd2ts2 is psd2ts with a seed. */
static void test_psd2ts(void **state){
	(void)state;
	const real dt=1e-3;
	const int nstep=10000;
	//tabulated flat PSD up to 100 Hz
	dmat *f=dlinspace(0, 1, 101);
	dmat *psd=psd_table(f, fun_one);
	real var=0;
	const int nseed=4;
	for(int seed=1; seed<=nseed; seed++){
		dmat *ts=psd2ts2(psd, seed, dt, nstep);
		assert_int_equal(NX(ts), nstep);
		assert_int_equal(NY(ts), 1);
		var+=dsumsq(ts)/nstep;
		dfree(ts);
	}
	var/=nseed;
	assert_float_equal(var, 100, 5);
	//same seed gives same result
	rand_t rstat;
	seed_rand(&rstat, 3);
	dmat *ts1=psd2ts(psd, &rstat, dt, 1000);
	dmat *ts2=psd2ts2(psd, 3, dt, 1000);
	assert_dmat_near(ts1, ts2, 0);
	dfree(ts1);
	dfree(ts2);
	dfree(psd);
	//power law form: flat (alpha=0) between fmin and fmax
	dmat *pl=dnew(4, 1);
	P(pl, 0)=0; P(pl, 1)=2; P(pl, 2)=10; P(pl, 3)=60;
	var=0;
	for(int seed=1; seed<=nseed; seed++){
		dmat *ts=psd2ts2(pl, seed, dt, nstep);
		var+=dsumsq(ts)/nstep;
		dfree(ts);
	}
	var/=nseed;
	assert_float_equal(var, 2*50, 6);
	//two element form spans the full frequency range up to 1/dt
	dresize(pl, 2, 1);
	P(pl, 1)=1e-4;
	dmat *ts=psd2ts2(pl, 1, dt, 4096);
	var=dsumsq(ts)/4096;
	assert_float_equal(var, 1e-4/dt, 0.1*1e-4/dt);
	dfree(ts);
	dfree(pl);
	dfree(f);
}

/** add_psd adds PSDs on the same grid exactly and interpolates otherwise. */
static void test_add_psd(void **state){
	(void)state;
	dmat *f=dlinspace(1, 1, 50);
	dmat *p1=psd_table(f, fun_lin);
	dmat *p2=psd_table(f, fun_one);
	dmat *sum=add_psd(p1, p2, 3);
	for(long i=0; i<50; i++){
		assert_float_equal(P(sum, i, 0), P(f, i), 0);
		assert_float_equal(P(sum, i, 1), fun_lin(P(f, i))+3, 1e-12);
	}
	dfree(sum);
	//different grid: linear interpolation of a linear function is exact
	dmat *f2=dlinspace(0.5, 0.5, 120);
	dmat *p3=psd_table(f2, fun_lin);
	sum=add_psd(p2, p3, 1);
	for(long i=0; i<50; i++){
		assert_float_equal(P(sum, i, 1), 1+fun_lin(P(f, i)), 1e-10);
	}
	dfree(sum);
	//add_psd2 to NULL makes a scaled copy
	dmat *out=NULL;
	add_psd2(&out, p1, 2);
	for(long i=0; i<50; i++){
		assert_float_equal(P(out, i, 0), P(f, i), 0);
		assert_float_equal(P(out, i, 1), 2*fun_lin(P(f, i)), 1e-12);
	}
	dfree(out);
	add_psd2(&out, p1, 1);
	assert_dmat_near(out, p1, 0);
	add_psd2(&out, p1, -1);
	for(long i=0; i<50; i++){
		assert_float_equal(P(out, i, 1), 0, 1e-12);
	}
	dfree(out);
	dfree(p3);
	dfree(f2);
	dfree(p2);
	dfree(p1);
	dfree(f);
}

/** psd_select sums the selected modes, scales and removes the DC component. */
static void test_psd_select(void **state){
	(void)state;
	dmat *psd=dnew(10, 4);
	for(long i=0; i<10; i++){
		P(psd, i, 0)=i;//f=0 at first
		for(int j=1; j<4; j++){
			P(psd, i, j)=j*10+i;
		}
	}
	//single mode, DC removed
	dmat *s=psd_select(psd, 1, 1, 0, 2);
	assert_int_equal(NX(s), 9);
	assert_int_equal(NY(s), 2);
	for(long i=0; i<9; i++){
		assert_float_equal(P(s, i, 0), i+1, 0);
		assert_float_equal(P(s, i, 1), 2*(20+i+1), 0);
	}
	dfree(s);
	//two modes, keep DC. nm<=0 means 1.
	s=psd_select(psd, 0, 2, 1, 1);
	assert_int_equal(NX(s), 10);
	assert_float_equal(P(s, 0, 1), 10+20, 0);
	dfree(s);
	s=psd_select(psd, 2, 0, 1, 1);
	assert_float_equal(P(s, 3, 1), 33, 0);
	dfree(s);
	//number of modes is clipped. all modes with im=-1.
	s=psd_select(psd, 1, 10, 1, 1);
	assert_float_equal(P(s, 0, 1), 20+30, 0);
	dfree(s);
	s=psd_select(psd, -1, 0, 1, 1);
	assert_float_equal(P(s, 5, 1), 15+25+35, 0);
	dfree(s);
	//no DC to remove if first frequency is not 0
	P(psd, 0, 0)=0.5;
	s=psd_select(psd, -1, 0, 0, 1);
	assert_int_equal(NX(s), 10);
	dfree(s);
	dfree(psd);
}

/** psd2d on a cell of screens averages; empty input returns NULL; non-square screens are rejected. */
static void test_psd2d_cell(void **state){
	(void)state;
	map_t *s1=genatm_simple(0.2, 0, -11./3., 1./32., 64, 1);
	map_t *s2=genatm_simple(0.2, 0, -11./3., 1./32., 64, 2);
	dcell *sc=dcellnew(2, 1);
	P(sc, 0)=dref(s1->dmat);
	P(sc, 1)=dref(s2->dmat);
	dmat *pa=psd2d(NULL, s1->dmat, 1./32.);
	dmat *pb=psd2d(NULL, s2->dmat, 1./32.);
	dmat *extra=NULL;
	dmat *pc=psd2d(&extra, sc, 1./32.);
	assert_int_equal(NX(pc), NX(pa));
	//The power law correction factor depends on the fitted slope, so compare the shape loosely.
	for(long i=1; i<NX(pc); i++){
		assert_float_equal(P(pc, i, 0), P(pa, i, 0), 1e-12);
	}
	//fitted slope is close to -11/3 and r0 close to 0.2
	assert_float_equal(P(extra, 1), -11./3., 0.5);
	assert_float_equal(P(extra, 0), 0.2, 0.1);
	real ratio=psd_inte2(pc)/(0.5*(psd_inte2(pa)+psd_inte2(pb)));
	assert_float_equal(ratio, 1, 0.3);
	dfree(extra);
	dfree(pa); dfree(pb); dfree(pc);
	assert_null(psd2d(NULL, (dmat *)NULL, 1));
	//non-square input
	dmat *ns=dnew(4, 8);
	signal_caught=0;
	assert_null(psd2d_aniso(ns, 1));
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	dfree(ns);
	dcellfree(sc);
	mapfree(s1);
	mapfree(s2);
}

UNITTEST_MAIN(
	cmocka_unit_test(test_psd2d),
	cmocka_unit_test(test_stfun),
	cmocka_unit_test(test_psd1d),
	cmocka_unit_test(test_psd_interp1),
	cmocka_unit_test(test_psd_vibid),
	cmocka_unit_test(test_psd_t2s),
	cmocka_unit_test(test_psd_inte),
	cmocka_unit_test(test_psd_reverse_cumu),
	cmocka_unit_test(test_psd2ts),
	cmocka_unit_test(test_add_psd),
	cmocka_unit_test(test_psd_select),
	cmocka_unit_test(test_psd2d_cell)
)
