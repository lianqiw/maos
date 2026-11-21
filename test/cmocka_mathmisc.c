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
  Unit tests for math/mathmisc.c (integer helpers, golden section search,
  readvec/readspint type conversion) and the vector routines of math/vec.c for
  all numeric types. Data for readvec is passed through a pipe, not a file.
*/
#include <unistd.h>
#include "unittest.h"

/** factorial() computes n1*(n1+1)*...*n2 and is 1 for an empty range. */
static void mathmisc_factorial(void **state){
	(void)state;
	assert_float_equal((double)factorial(1, 5), 120, 0);
	assert_float_equal((double)factorial(3, 5), 60, 0);
	assert_float_equal((double)factorial(7, 7), 7, 0);
	assert_float_equal((double)factorial(5, 4), 1, 0);//empty range
	assert_float_equal((double)factorial(0, 3), 0, 0);//includes 0
	//20! fits in long; 25! exceeds long but not long double
	long f20=1;
	for(long i=2; i<=20; i++) f20*=i;
	assert_float_equal((double)factorial(1, 20), (double)f20, 0);
	assert_float_equal((double)(factorial(1, 25)/factorial(1, 24)), 25, 1e-9);
	//Overflow is reported by error() but the infinite value is still returned.
	assert_true(isinf((double)factorial(1, 5000)));
}

/** calc_gcd() and calc_lcm() for assorted values, including negative ones. */
static void mathmisc_gcd_lcm(void **state){
	(void)state;
	assert_int_equal(calc_gcd(12, 18), 6);
	assert_int_equal(calc_gcd(18, 12), 6);
	assert_int_equal(calc_gcd(7, 7), 7);
	assert_int_equal(calc_gcd(17, 5), 1);
	assert_int_equal(calc_gcd(0, 5), 5);
	assert_int_equal(calc_gcd(5, 0), 5);
	assert_int_equal(calc_gcd(1071, 462), 21);
	assert_int_equal(calc_lcm(4, 6), 12);
	assert_int_equal(calc_lcm(6, 4), 12);
	assert_int_equal(calc_lcm(7, 13), 91);
	assert_int_equal(calc_lcm(0, 0), 0);
	assert_int_equal(calc_lcm(5, 5), 5);
	//negative input uses absolute values
	assert_int_equal(calc_gcd(-12, 18), 6);
	assert_int_equal(calc_gcd(12, -18), 6);
	assert_int_equal(calc_gcd(-12, -18), 6);
	assert_int_equal(calc_gcd(-7, 0), 7);
	assert_int_equal(calc_lcm(-4, 6), 12);
	//brute force comparison
	for(long a=1; a<30; a++){
		for(long b=1; b<30; b++){
			long g=1;
			for(long k=1; k<=a&&k<=b; k++){
				if(a%k==0&&b%k==0) g=k;
			}
			assert_int_equal(calc_gcd(a, b), g);
			assert_int_equal(calc_lcm(a, b), a*b/g);
		}
	}
}

/** invperm() inverts a permutation, maxminlong() finds extremes. */
static void mathmisc_perm_maxmin(void **state){
	(void)state;
	long p[6]={3, 0, 5, 1, 4, 2};
	long *pi=invperm(p, 6);
	assert_non_null(pi);
	for(long i=0; i<6; i++){
		assert_int_equal(pi[p[i]], i);
		assert_int_equal(p[pi[i]], i);
	}
	free(pi);
	assert_null(invperm(NULL, 0));

	long v[7]={4, -9, 12, 0, 12, -9, 3};
	long mx=0, mn=0;
	maxminlong(v, 7, &mx, &mn);
	assert_int_equal(mx, 12);
	assert_int_equal(mn, -9);
	maxminlong(v, 7, NULL, &mn);//optional output
	assert_int_equal(mn, -9);
	mx=mn=100;
	maxminlong(NULL, 7, &mx, &mn);//no-op
	assert_int_equal(mx, 100);
	assert_int_equal(mn, 100);
	maxminlong(v, 0, &mx, &mn);//empty: no-op, does not read v[0]
	assert_int_equal(mx, 100);
	assert_int_equal(mn, 100);
}

/** nextpow2(), nextfftsize() and mylog2(). */
static void mathmisc_sizes(void **state){
	(void)state;
	assert_int_equal(nextpow2(1), 1);
	assert_int_equal(nextpow2(2), 2);
	assert_int_equal(nextpow2(3), 4);
	assert_int_equal(nextpow2(64), 64);
	assert_int_equal(nextpow2(65), 128);
	assert_int_equal(nextpow2(1000000), 1048576);
	for(unsigned long m=0; m<40; m++){
		assert_int_equal(mylog2(1UL<<m), m);
	}
	//nextfftsize: smallest even number >=n with only factors 2,3,5,7.
	for(long n=0; n<400; n++){
		long expect=n<2?2:n;
		for(;; expect++){
			if(expect&1) continue;
			long r=expect;
			const int rad[4]={2, 3, 5, 7};
			for(int i=0; i<4; i++) while(r%rad[i]==0) r/=rad[i];
			if(r==1) break;
		}
		assert_int_equal(nextfftsize(n), expect);
	}
	assert_int_equal(nextfftsize(22), 24);
	assert_int_equal(nextfftsize(98), 98);
	assert_int_equal(nextfftsize(1022), 1024);
}

static real parabola(void *param, real x){
	real x0=*(real *)param;
	return (x-x0)*(x-x0)+1;
}
static real absval(void *param, real x){
	(void)param;
	return fabs(x-0.25);
}
/** golden_section_search() finds the minimum of unimodal functions. */
static void mathmisc_golden(void **state){
	(void)state;
	real x0=1.234;
	real xm=golden_section_search(parabola, &x0, -3, 7, 1e-8);
	assert_float_equal(xm, x0, 1e-5);
	x0=-2.5;
	xm=golden_section_search(parabola, &x0, -10, 1, 1e-8);
	assert_float_equal(xm, x0, 1e-5);
	xm=golden_section_search(absval, NULL, 0, 2, 1e-8);
	assert_float_equal(xm, 0.25, 1e-6);
}

/*Write raw bytes into a pipe and open the read end as file_t.*/
static file_t *pipe_with(const void *data, size_t len){
	int fds[2];
	if(pipe(fds)) return NULL;
	if(len&&write(fds[1], data, len)!=(ssize_t)len){
		close(fds[0]); close(fds[1]);
		return NULL;
	}
	close(fds[1]);
	return zfdopen_as(fds[0], "r", ".bin");
}
/*Fill buf with values from vals as type magic.*/
static size_t fill_as(void *buf, uint32_t magic, const long *vals, int n){
	switch(magic){
	case M_DBL: for(int i=0; i<n; i++) ((double *)buf)[i]=(double)vals[i]; return n*sizeof(double);
	case M_FLT: for(int i=0; i<n; i++) ((float *)buf)[i]=(float)vals[i]; return n*sizeof(float);
	case M_INT32: for(int i=0; i<n; i++) ((int *)buf)[i]=(int)vals[i]; return n*sizeof(int);
	case M_INT64: for(int i=0; i<n; i++) ((int64_t *)buf)[i]=(int64_t)vals[i]; return n*sizeof(int64_t);
	case M_CMP: for(int i=0; i<n; i++) ((dcomplex *)buf)[i]=COMPLEX((double)vals[i], 0.); return n*sizeof(dcomplex);
	case M_ZMP: for(int i=0; i<n; i++) ((fcomplex *)buf)[i]=FCOMPLEX((float)vals[i], 0.f); return n*sizeof(fcomplex);
	default: return 0;
	}
}
/*Get element i of buf of type magic as double.*/
static double get_as(const void *buf, uint32_t magic, int i){
	switch(magic){
	case M_DBL: return ((const double *)buf)[i];
	case M_FLT: return ((const float *)buf)[i];
	case M_INT32: return ((const int *)buf)[i];
	case M_INT64: return (double)((const int64_t *)buf)[i];
	case M_CMP: return creal(((const dcomplex *)buf)[i]);
	case M_ZMP: return crealf(((const fcomplex *)buf)[i]);
	default: return NAN;
	}
}
/** readvec() converts between all supported data types. */
static void mathmisc_readvec(void **state){
	(void)state;
	const uint32_t types[]={M_DBL, M_FLT, M_INT32, M_INT64, M_CMP, M_ZMP};
	const int nt=sizeof(types)/sizeof(types[0]);
	const long vals[5]={3, -7, 0, 1024, 42};
	const int n=5;
	for(int ifile=0; ifile<nt; ifile++){
		for(int ip=0; ip<nt; ip++){
			uint32_t mf=types[ifile], mp=types[ip];
			char src[5*16], dst[5*16];
			memset(dst, 0, sizeof(dst));
			size_t len=fill_as(src, mf, vals, n);
			size_t size=(mp==M_DBL||mp==M_INT64||mp==M_ZMP)?8:(mp==M_CMP?16:4);
			file_t *fp=pipe_with(src, len);
			assert_non_null(fp);
			int ans=readvec(dst, mp, mf, size, n, fp);
			zfclose(fp);
			int fcomp=(mf==M_CMP||mf==M_ZMP);
			int pcomp=(mp==M_CMP||mp==M_ZMP);
			if(fcomp&&!pcomp){//complex to real is not supported
				assert_int_not_equal(ans, 0);
			}else{
				assert_int_equal(ans, 0);
				for(int i=0; i<n; i++){
					assert_float_equal(get_as(dst, mp, i), vals[i], 0);
				}
			}
		}
	}
	//invalid arguments
	double d[2];
	assert_int_equal(readvec(NULL, M_DBL, M_DBL, 8, 2, NULL), -1);
	file_t *fp=pipe_with(NULL, 0);
	assert_int_equal(readvec(d, M_DBL, M_DBL, 8, 0, fp), 0);//nothing to read
	assert_int_not_equal(readvec(d, M_DBL, M_DBL, 8, 2, fp), 0);//premature end of file
	zfclose(fp);
}
/** readspint() reads a dense array and converts it to spint. */
static void mathmisc_readspint(void **state){
	(void)state;
	int fds[2];
	assert_int_equal(pipe(fds), 0);
	file_t *fw=zfdopen_as(fds[1], "w", ".bin");
	dmat *A=dnew(3, 2);
	for(long i=0; i<PN(A); i++) P(A, i)=i*3-4;
	writedata(fw, A, 0);//small enough to fit in the pipe buffer
	zfclose(fw);
	file_t *fr=zfdopen_as(fds[0], "r", ".bin");
	long nx=0, ny=0;
	spint *p=readspint(fr, &nx, &ny);
	zfclose(fr);
	assert_non_null(p);
	assert_int_equal(nx, 3);
	assert_int_equal(ny, 2);
	for(long i=0; i<PN(A); i++){
		assert_int_equal(p[i], (spint)P(A, i));
	}
	free(p);
	dfree(A);
}

/*Check the vec* functions of one type. V is the value of index i.*/
#define TEST_VEC(X, T, R, V, ISCOMP, TOL)				\
	{								\
		const long n=7;						\
		T p[7], q[7]; R w[7];					\
		T s0=0; T sq=0; T sw=0; T sqw=0; T spw=0; T spqw=0; R amax=(R)ISCOMP(V(0)), amin=amax; \
		for(long i=0; i<n; i++){				\
			p[i]=V(i); q[i]=V(n-i); w[i]=(R)(i+1);		\
			s0+=p[i]; sq+=q[i]; sw+=(T)w[i]; sqw+=q[i]*w[i]; spw+=p[i]*w[i]; spqw+=p[i]*q[i]*w[i]; \
			R a=(R)ISCOMP(p[i]);				\
			if(a>amax){ amax=a; } if(a<amin){ amin=a; }	\
		}							\
		T spq=0; for(long i=0; i<n; i++) spq+=p[i]*q[i];	\
		assert_float_equal(creal(X##vecsum(p, n)), creal(s0), TOL); \
		assert_float_equal(creal(X##vecsum(NULL, 0)), 0, 0);	\
		assert_float_equal(creal(X##vecdot(p, q, w, n)), creal(spqw), TOL); \
		assert_float_equal(creal(X##vecdot(NULL, q, w, n)), creal(sqw), TOL); \
		assert_float_equal(creal(X##vecdot(p, NULL, w, n)), creal(spw), TOL); \
		assert_float_equal(creal(X##vecdot(NULL, NULL, w, n)), creal(sw), TOL); \
		assert_float_equal(creal(X##vecdot(p, q, NULL, n)), creal(spq), TOL); \
		assert_float_equal(creal(X##vecdot(NULL, q, NULL, n)), creal(sq), TOL); \
		assert_float_equal(creal(X##vecdot(p, NULL, NULL, n)), creal(s0), TOL); \
		assert_float_equal(creal(X##vecdot(NULL, NULL, NULL, n)), n, 0); \
		R mx=0, mn=0;						\
		X##vecmaxmin(p, n, &mx, &mn);				\
		assert_float_equal(mx, amax, TOL);			\
		assert_float_equal(mn, amin, TOL);			\
		X##vecmaxmin(p, n, NULL, NULL);				\
		R ma=X##vecmaxabs(p, n);				\
		assert_float_equal(ma, fabs((double)amax)>fabs((double)amin)?fabs((double)amax):fabs((double)amin), TOL); \
		assert_float_equal(X##vecmaxabs(NULL, 0), 0, 0);	\
	}

#define VD(i) ((double)(i)*1.5-4.)
#define VS(i) ((float)(i)*1.5f-4.f)
#define VL(i) ((long)(i)*3-10)
#define VC(i) COMPLEX((double)(i)-2., 0.5*(double)(i))
#define VZ(i) FCOMPLEX((float)(i)-2.f, 0.5f*(float)(i))
#define ID_R(a) (a)
#define ABS_C(a) cabs(a)
#define ABS_Z(a) cabsf(a)
/** vecsum/vecdot/vecmaxmin/vecmaxabs for d, s, c, z, l types. */
static void mathmisc_vec(void **state){
	(void)state;
	TEST_VEC(d, real, real, VD, ID_R, 1e-10);
	TEST_VEC(s, float, float, VS, ID_R, 1e-4);
	TEST_VEC(l, long, long, VL, ID_R, 0);
	TEST_VEC(c, comp, real, VC, ABS_C, 1e-10);
	TEST_VEC(z, fcomplex, float, VZ, ABS_Z, 1e-4);
}
/** vecmaxmin skips NaN except for the first element; vecnormalize_* rescale. */
static void mathmisc_vec_normalize(void **state){
	(void)state;
	real p[5]={1, NAN, -3, 2, NAN};
	real mx, mn;
	dvecmaxmin(p, 5, &mx, &mn);
	assert_float_equal(mx, 2, 0);
	assert_float_equal(mn, -3, 0);
	float ps[4]={2, -6, 1, 7};
	svecnormalize_sum(ps, 4, 10);
	assert_float_equal(ps[0]+ps[1]+ps[2]+ps[3], 10, 1e-5);
	assert_float_equal(ps[1], -15, 1e-4);
	svecnormalize_sumabs(ps, 4, 2);
	assert_float_equal(fabs(ps[0])+fabs(ps[1])+fabs(ps[2])+fabs(ps[3]), 2, 1e-5);
	assert_float_equal(ps[1], -0.75, 1e-5);
	real pd[3]={1, 2, 5};
	dvecnormalize_sum(pd, 3, 1);
	assert_float_equal(pd[2], 0.625, 1e-12);
	dvecnormalize_sumabs(pd, 3, 16);
	assert_float_equal(pd[0], 2, 1e-12);
	comp pc[2]={COMPLEX(3., 4.), COMPLEX(0., -5.)};
	cvecnormalize_sumabs(pc, 2, 1);
	assert_float_equal(cabs(pc[0]), 0.5, 1e-12);
	assert_float_equal(cimag(pc[1]), -0.5, 1e-12);
	cvecnormalize_sum(pc, 2, 2);
	assert_float_equal(creal(pc[0]+pc[1]), 2, 1e-12);
	assert_float_equal(cimag(pc[0]+pc[1]), 0, 1e-12);
	long pl[3]={2, 4, 6};
	lvecnormalize_sum(pl, 3, 24);//integer scale factor 2
	assert_int_equal(pl[2], 12);
	//NULL with zero length is a no-op
	dvecnormalize_sum(NULL, 0, 1);
	dvecnormalize_sumabs(NULL, 0, 1);
	//sinc and mysqrt macros
	assert_float_equal(sinc(0.), 1, 0);
	assert_float_equal(sinc(0.5), 2/M_PI, 1e-15);
	assert_float_equal(sinc(2.), 0, 1e-15);
	real m4=-4;
	assert_float_equal(mysqrt(m4), -2, 0);
	assert_float_equal(mysqrt(9.), 3, 0);
}
/**
   dbessik() and sbessik() against known values of the modified Bessel functions.
   x<2 uses the Temme series (Chebyshev expansion of Gamma1 and Gamma2 in
   beschb), x>=2 uses the continued fraction.
*/
static void mathmisc_bessik(void **state){
	(void)state;
	//nu=1/2 has closed forms.
	const double xs[]={0.05, 0.3, 1.0, 1.7, 2.5, 6.0};
	for(int i=0; i<6; i++){
		const double x=xs[i];
		double ri, rk, rip, rkp;
		dbessik(x, 0.5, &ri, &rk, &rip, &rkp);
		const double rk0=sqrt(M_PI/(2*x))*exp(-x);
		const double ri0=sqrt(2/(M_PI*x))*sinh(x);
		//the float precision Chebyshev coefficients of Gamma1/2 gave errors of 1e-9.
		assert_float_equal(rk, rk0, 1e-12*rk0);
		assert_float_equal(ri, ri0, 1e-12*ri0);
		assert_float_equal(rkp, -rk0*(1+0.5/x), 1e-12*rk0*(1+0.5/x));
		//Wronskian: I*K'-I'*K=-1/x
		assert_float_equal(ri*rkp-rip*rk, -1./x, 1e-14/x);
	}
	//Reference values from numerical integration (K) and power series (I).
	const struct{
		double nu, x, rk, ri;
	}refs[]={
		{5./6., 0.05, 12.098978810177696, 0.04916627304756478},
		{5./6., 0.3, 2.4207665691474385, 0.22146092075008514},
		{5./6., 1.0, 0.5405300878883402, 0.6816661681007694},
		{5./6., 1.7, 0.1949587123768023, 1.3440171825012768},
		{5./6., 2.5, 0.07017009013303842, 2.7168036645992255},
		{5./6., 6.0, 0.0013125989438328982, 63.080049983827344},
		{1./3., 0.05, 3.991017706867745, 0.32759861981214433},
		{1./3., 0.3, 1.5091129245822228, 0.6050965184152769},
		{1./3., 1.0, 0.43843063344154654, 1.0646313978895294},
		{1./3., 1.7, 0.16992206277162125, 1.7313127900072365},
	};
	for(unsigned int i=0; i<sizeof(refs)/sizeof(refs[0]); i++){
		double ri, rk, rip, rkp;
		dbessik(refs[i].x, refs[i].nu, &ri, &rk, &rip, &rkp);
		assert_float_equal(rk, refs[i].rk, 1e-12*refs[i].rk);
		assert_float_equal(ri, refs[i].ri, 1e-12*refs[i].ri);
		assert_float_equal(ri*rkp-rip*rk, -1./refs[i].x, 1e-14/refs[i].x);
		float fi, fk, fip, fkp;
		sbessik((float)refs[i].x, (float)refs[i].nu, &fi, &fk, &fip, &fkp);
		assert_float_equal(fk, refs[i].rk, 1e-5*refs[i].rk);
		assert_float_equal(fi, refs[i].ri, 1e-5*refs[i].ri);
	}
}

UNITTEST_MAIN(
	cmocka_unit_test(mathmisc_factorial),
	cmocka_unit_test(mathmisc_gcd_lcm),
	cmocka_unit_test(mathmisc_perm_maxmin),
	cmocka_unit_test(mathmisc_sizes),
	cmocka_unit_test(mathmisc_golden),
	cmocka_unit_test(mathmisc_readvec),
	cmocka_unit_test(mathmisc_readspint),
	cmocka_unit_test(mathmisc_vec),
	cmocka_unit_test(mathmisc_vec_normalize),
	cmocka_unit_test(mathmisc_bessik),
)
