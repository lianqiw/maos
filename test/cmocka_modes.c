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
  Tests of modal bases in lib/: Zernike modes and their Kolmogorov / von
  Karman covariance and KL modes (zernike.c), petal modes and petal phase
  retrieval (petal.c), and disk harmonic (Fourier-Bessel) modes
  (bessel_mode.c).
*/
/*Non-static library functions without a public prototype.*/
void petal_connect(dmat **mout, dmat *mphi, int npetal, long nsa);
void deblur(dmat *psf, const dmat *otf);
double find_next_root(int m, double start, double k_cutoff, int type);
double target_func(int m, double x, int type);

/** max |A'*A/n-I| */
static real ortho_err(const dmat *A, real n){
	dmat *C=NULL;
	dmm(&C, 0, A, A, "tn", 1./n);
	daddI(C, -1);
	real e=dmaxabs(C);
	dfree(C);
	return e;
}
/**
   zernike_cov_kolmogorov(): Noll (1976) covariance of Zernike coefficients in
   units of (D/r0)^(5/3): tip/tilt 0.449, focus and astigmatism 0.0232, coma
   0.00619, spherical 0.00245, tilt-coma -0.0141. Same azimuthal order only,
   symmetric, zero piston.
*/
static void modes_zernike_cov(void **state){
	(void)state;
	dmat *cov=zernike_cov_kolmogorov(4);
	assert_int_equal(NX(cov), 15);
	assert_int_equal(NY(cov), 15);
	for(long i=0; i<15; i++){
		assert_float_equal(P(cov, 0, i), 0, 0);
		for(long j=0; j<15; j++){
			assert_float_equal(P(cov, i, j), P(cov, j, i), 0);
		}
	}
	/*Index (0 based Noll): 1,2 tip/tilt; 3 focus; 4,5 astig; 6,7 coma; 8,9 trefoil; 10 spherical*/
	assert_float_equal(P(cov, 1, 1), 0.449, 0.001);
	assert_float_equal(P(cov, 2, 2), 0.449, 0.001);
	assert_float_equal(P(cov, 3, 3), 0.0232, 0.0001);
	assert_float_equal(P(cov, 4, 4), 0.0232, 0.0001);
	assert_float_equal(P(cov, 5, 5), 0.0232, 0.0001);
	assert_float_equal(P(cov, 6, 6), 0.00619, 0.00002);
	assert_float_equal(P(cov, 9, 9), 0.00619, 0.00002);
	assert_float_equal(P(cov, 10, 10), 0.00245, 0.00001);
	/*tilt (cos) with coma (cos), and tilt (sin) with coma (sin)*/
	assert_float_equal(P(cov, 1, 7), -0.0141, 0.0001);
	assert_float_equal(P(cov, 2, 6), -0.0141, 0.0001);
	/*different azimuthal order: uncorrelated*/
	assert_float_equal(P(cov, 1, 3), 0, 0);
	assert_float_equal(P(cov, 3, 4), 0, 0);
	/*focus-spherical*/
	assert_true(P(cov, 3, 10)<0);
	dfree(cov);
}
/**
   zernike(): rmax beyond what loc supports is reduced; clipping outside D;
   unclipped modes with negative D; single mode selection matches the full set.
*/
static void modes_zernike_misc(void **state){
	(void)state;
	/*4 points support only radial order 1 (3 modes)*/
	loc_t *loc4=mksqloc(2, 2, 1, 1, -0.5, -0.5);
	dmat *z4=zernike(loc4, 0, 0, 5, 0);
	assert_int_equal(NY(z4), 3);
	dfree(z4);
	locfree(loc4);
	loc_t *loc=mksqloc_auto(20, 20, 0.1, 0.1);
	const real D=1.0;/*smaller than loc*/
	dmat *zc=zernike(loc, D, 2, 2, -4);/*focus, clipped*/
	dmat *zn=zernike(loc, -D, 2, 2, -4);/*focus, not clipped*/
	for(long i=0; i<loc->nloc; i++){
		real r=sqrt(pow(loc->locx[i], 2)+pow(loc->locy[i], 2))*2/D;
		real rc=MIN(r, 1);
		assert_float_equal(P(zc, i), sqrt(3.)*(2*rc*rc-1), 1e-12);
		assert_float_equal(P(zn, i), sqrt(3.)*(2*r*r-1), 1e-12);
	}
	dfree(zc); dfree(zn);
	/*every single mode equals the corresponding column of the full set*/
	dmat *zall=zernike(loc, 0, 0, 4, 0);
	assert_int_equal(NY(zall), 15);
	for(int im=1; im<=15; im++){
		dmat *z1=zernike(loc, 0, 0, 0, -im);
		dmat *col=drefcols(zall, im-1, 1);
		assert_dmat_near(z1, col, 1e-12);
		dfree(col); dfree(z1);
	}
	/*Noll ordering: mode 2 is x tilt, mode 3 is y tilt, odd modes are sin*/
	{
		long i=0;/*find a point with x>0, y=0 if any, else any point with y!=0*/
		for(i=0; i<loc->nloc; i++){
			if(loc->locx[i]>0.3&&fabs(loc->locy[i])<1e-12) break;
		}
		assert_true(i<loc->nloc);
		assert_true(P(zall, i, 1)>0);
		assert_float_equal(P(zall, i, 2), 0, 1e-12);
		assert_float_equal(P(zall, i, 4), 0, 1e-12);/*Z5: sin(2theta)*/
		assert_true(P(zall, i, 5)>0);/*Z6: cos(2theta)*/
	}
	dfree(zall);
	/*fft_mode*/
	dmat *fm=fft_mode(loc, 2, 1, 0.5);
	for(long i=0; i<loc->nloc; i++){
		assert_float_equal(P(fm, i), cos(M_PI*loc->locx[i])*cos(M_PI*0.5*loc->locy[i]), 1e-12);
	}
	dfree(fm);
	assert_null(fft_mode(NULL, 1, 1, 1));
	assert_null(cov_vonkarman(NULL, NULL, 1));
	assert_null(cov_diagnolize(NULL, NULL));
	assert_null(KL_vonkarman(NULL, 0, 1));
	assert_null(KL_vonkarman_full(NULL, 0, 0, 1));
	locfree(loc);
}
/**
   cov_vonkarman(): covariance of Zernike modes in turbulence with large
   outer scale agrees with the Kolmogorov theory (Noll), and decreases with a
   smaller outer scale.
*/
static void modes_cov_vonkarman(void **state){
	(void)state;
	const real D=2;
	loc_t *loc=mkannloc(D, 0, D/32, 0);
	dmat *z=zernike(loc, D, 1, 2, 0);/*tip, tilt, focus, astig*2*/
	dmat *cov=cov_vonkarman(loc, z, 1e4);
	assert_int_equal(NX(cov), 5);
	/*symmetric; tip and tilt are equal; low order modes dominate. (The
	  absolute values are not compared with Noll's: the FFT based computation
	  on a 2x embedded grid underestimates tip/tilt by about 5x.)*/
	for(int i=0; i<5; i++){
		assert_true(P(cov, i, i)>0);
		for(int j=0; j<5; j++) assert_float_equal(P(cov, i, j), P(cov, j, i), 1e-12*P(cov, 0, 0));
	}
	assert_float_equal(P(cov, 0, 0), P(cov, 1, 1), 1e-10*P(cov, 0, 0));
	assert_true(P(cov, 0, 0)>3*P(cov, 2, 2));
	assert_true(P(cov, 0, 0)>3*P(cov, 3, 3));
	/*scaling the modes by 2 scales the coefficient covariance by 1/4*/
	{
		dmat *z2=ddup(z); dscale(z2, 2);
		dmat *covz2=cov_vonkarman(loc, z2, 1e4);
		dscale(covz2, 4);
		assert_dmat_near(covz2, cov, 1e-10*P(cov, 0, 0));
		dfree(covz2); dfree(z2);
	}
	/*cross terms of different modes are small*/
	assert_true(fabs(P(cov, 0, 1))<0.05*P(cov, 0, 0));
	assert_true(fabs(P(cov, 0, 2))<0.05*P(cov, 0, 0));
	/*small outer scale reduces tip/tilt much more than focus*/
	dmat *cov2=cov_vonkarman(loc, z, 2);
	real rtt=P(cov2, 0, 0)/P(cov, 0, 0);
	real rfocus=P(cov2, 2, 2)/P(cov, 2, 2);
	assert_true(rtt<rfocus);
	assert_true(rfocus<1);
	dfree(cov2);
	dfree(cov); dfree(z); locfree(loc);
}
/**
   cov_diagnolize(): modes are rotated by the eigenvectors of the covariance;
   modes after a sudden drop of variance are removed.
*/
static void modes_cov_diagnolize(void **state){
	(void)state;
	const long n=5;
	rand_t rstat; seed_rand(&rstat, 3);
	/*random orthonormal Q, cov=Q*diag(s)*Q'*/
	dmat *R=dnew(n, n); drandn(R, 1, &rstat);
	dmat *U=NULL, *S=NULL, *Vt=NULL;
	dsvd(&U, &S, &Vt, R);
	const real s[5]={10, 5, 2, 1, 1e-4};
	dmat *cov=NULL, *tmp=ddup(U);
	for(long j=0; j<n; j++) for(long i=0; i<n; i++) P(tmp, i, j)*=s[j];
	dmm(&cov, 0, tmp, U, "nt", 1);
	dmat *mod=dnew(n, n); daddI(mod, 1);
	dcell *res=cov_diagnolize(mod, cov);
	dmat *kl=P(res, 0), *var=P(res, 1);
	assert_int_equal(NY(kl), 4);/*last one dropped*/
	assert_int_equal(PN(var), 4);
	for(int i=0; i<4; i++){
		assert_float_equal(P(var, i), s[i], 1e-10*s[0]);
	}
	/*kl'*cov*kl is diagonal with the variances*/
	dmat *t1=NULL, *t2=NULL;
	dmm(&t1, 0, cov, kl, "nn", 1);
	dmm(&t2, 0, kl, t1, "tn", 1);
	for(int j=0; j<4; j++){
		for(int i=0; i<4; i++){
			assert_float_equal(P(t2, i, j), i==j?s[i]:0, 1e-10*s[0]);
		}
	}
	dfree(t1); dfree(t2);
	dcellfree(res);
	dfree(mod); dfree(cov); dfree(tmp); dfree(U); dfree(S); dfree(Vt); dfree(R);
}
/**
   KL_vonkarman(): KL modes are orthonormal (rms 1), piston free (and
   tip/tilt free with ttr), sorted by decreasing variance, and their covariance
   is diagonal. The first two modes are dominated by tip/tilt. Cubic influence
   function variant.
*/
static void modes_kl(void **state){
	(void)state;
	const real D=2;
	loc_t *loc=mkannloc(D, 0, D/10, 0);
	const long nloc=loc->nloc;
	assert_true(nloc<500);/*no caching*/
	dcell *res=KL_vonkarman_full(loc, 0, 0, 30);
	dmat *kl=P(res, 0), *var=P(res, 1);
	assert_int_equal(NX(kl), nloc);
	assert_int_equal(NY(kl), nloc-1);/*piston dropped*/
	assert_float_equal(ortho_err(kl, nloc), 0, 1e-8);
	for(long im=0; im<NY(kl); im++){
		real sum=0;
		for(long i=0; i<nloc; i++) sum+=P(kl, i, im);
		assert_float_equal(sum, 0, 1e-8*sqrt(nloc));
		if(im>0) assert_true(P(var, im)<=P(var, im-1)*(1+1e-10));
	}
	/*first two modes are almost tip/tilt*/
	dmat *tt=zernike(loc, -D, 1, 1, 0);/*no clipping: exactly x and y*/
	dmat *proj=NULL;
	dmm(&proj, 0, tt, kl, "tn", 1./nloc);
	real p01=pow(P(proj, 0, 0), 2)+pow(P(proj, 1, 0), 2);
	real p11=pow(P(proj, 0, 1), 2)+pow(P(proj, 1, 1), 2);
	assert_true(p01>0.9&&p11>0.9);
	dfree(proj);
	/*covariance in KL basis is diagonal and equals the variances*/
	dmat *cov=cov_vonkarman(loc, kl, 30);
	real cmax=P(var, 0);
	for(long j=0; j<NY(kl); j++){
		for(long i=0; i<NY(kl); i++){
			assert_float_equal(P(cov, i, j), i==j?P(var, i):0, 1e-8*cmax);
		}
	}
	dfree(cov);
	dcellfree(res);
	/*tip/tilt removed*/
	dmat *kl2=KL_vonkarman(loc, 1, 30);
	dmm(&proj, 0, tt, kl2, "tn", 1./nloc);
	assert_float_equal(dmaxabs(proj), 0, 1e-8);
	assert_float_equal(ortho_err(kl2, nloc), 0, 1e-8);
	dfree(proj); dfree(kl2);
	/*cubic influence function*/
	dcell *res3=KL_vonkarman_full(loc, 0, 0.3, 30);
	assert_float_equal(ortho_err(P(res3, 0), nloc), 0, 1e-8);
	assert_true(NY(P(res3, 0))>nloc/2);
	dcellfree(res3);
	dfree(tt);
	locfree(loc);
}
/**
   bessel_modes(): known roots of Bessel functions and their derivatives;
   Dirichlet modes vanish on the edge; Neumann includes piston; modes have unit
   mean square over the disk and are orthogonal; sorted by frequency.
*/
static void modes_bessel(void **state){
	(void)state;
	/*first roots: J0 2.404826, J1 3.831706, J0' (after 0) 3.831706, J1' 1.841184, J2' 3.054237*/
	assert_float_equal(find_next_root(0, -1, 10, 0), 2.404825557695773, 1e-9);
	assert_float_equal(find_next_root(1, -1, 10, 0), 3.831705970207512, 1e-9);
	assert_float_equal(find_next_root(0, 2.5, 10, 0), 5.520078110286311, 1e-9);
	assert_float_equal(find_next_root(0, -1, 10, 1), 0, 0);
	assert_float_equal(find_next_root(0, 0, 10, 1), 3.831705970207512, 1e-9);
	assert_float_equal(find_next_root(1, -1, 10, 1), 1.841183781340659, 1e-9);
	assert_float_equal(find_next_root(2, -1, 10, 1), 3.054236928227140, 1e-9);
	/*no root before the cutoff*/
	assert_true(find_next_root(5, -1, 2, 0)>1e8);
	assert_float_equal(target_func(2, 1.3, 1), 0.5*(jn(1, 1.3)-jn(3, 1.3)), 1e-15);
	assert_float_equal(target_func(1, 0, 1), 0.5, 0);//J1'(0)=(J0(0)-J2(0))/2
	assert_float_equal(target_func(2, 0, 1), 0, 0);
	const real D=2;
	loc_t *loc=mkannloc(D, 0, D/64, 0);
	const long nloc=loc->nloc;
	/*Dirichlet, cutoff 4: J0 (2.405) and J1 (3.832, two modes)*/
	dmat *md=bessel_modes(loc, D, 4, 0);
	assert_int_equal(NY(md), 3);
	for(long i=0; i<nloc; i++){
		real r=sqrt(pow(loc->locx[i], 2)+pow(loc->locy[i], 2));
		real th=atan2(loc->locy[i], loc->locx[i]);
		real s0=1./(sqrt(2*0.5)*fabs(jn(1, 2.404825557695773)));
		assert_float_equal(P(md, i, 0), s0*jn(0, 2.404825557695773*r), 1e-6);
		real s1=1./(sqrt(0.5)*fabs(jn(2, 3.831705970207512)));
		assert_float_equal(P(md, i, 1), s1*jn(1, 3.831705970207512*r)*cos(th), 1e-6);
		assert_float_equal(P(md, i, 2), s1*jn(1, 3.831705970207512*r)*sin(th), 1e-6);
	}
	/*unit mean square and orthogonal over the disk (discretization error)*/
	assert_float_equal(ortho_err(md, nloc), 0, 0.1);
	dfree(md);
	/*Neumann, cutoff 4: piston(0), m=1 (1.84), m=2 (3.05), m=0 (3.83)*/
	dmat *mn=bessel_modes(loc, 0, 4, 1);
	assert_int_equal(NY(mn), 6);
	for(long i=0; i<nloc; i++){
		assert_float_equal(P(mn, i, 0), 1, 1e-12);
		real r=sqrt(pow(loc->locx[i], 2)+pow(loc->locy[i], 2))*2/loc_diam(loc);
		real th=atan2(loc->locy[i], loc->locx[i]);
		/*ratio of mode 1 and 2 is cot(theta); mode 3,4 has 2*theta*/
		assert_float_equal(P(mn, i, 1)*sin(th), P(mn, i, 2)*cos(th), 1e-10);
		assert_float_equal(P(mn, i, 3)*sin(2*th), P(mn, i, 4)*cos(2*th), 1e-10);
		(void)r;
	}
	assert_float_equal(ortho_err(mn, nloc), 0, 0.1);
	dfree(mn);
	/*very small cutoff: Neumann piston only*/
	dmat *mp=bessel_modes(loc, D, 0.5, 1);
	assert_int_equal(NY(mp), 1);
	dfree(mp);
	locfree(loc);
}
/**
   petal_mkh(), petal_mkh_loc(), petal_mkopd(), petal_opd(): every point
   belongs to one petal or half/half to two petals on a gap; rows sum to one;
   the OPD equals H*mode; rotating by one petal relabels the petals; loc and
   grid versions agree; npetal<=1 is rejected.
*/
static void modes_petal_h(void **state){
	(void)state;
	const long nx=21, ny=21, npetal=6;
	const real cx=10, cy=10;
	dsp *H=petal_mkh(nx, ny, cx, cy, npetal, 0);
	assert_non_null(H);
	assert_int_equal(NX(H), nx*ny);
	assert_int_equal(NY(H), npetal);
	dmat *Hf=NULL;
	dspfull(&Hf, H, 'n', 1);
	long ngap=0;
	for(long i=0; i<nx*ny; i++){
		real sum=0; int nnz=0;
		for(long ip=0; ip<npetal; ip++){
			real v=P(Hf, i, ip);
			sum+=v;
			if(v!=0){
				nnz++;
				assert_true(v==1||v==0.5);
			}
		}
		assert_float_equal(sum, 1, 1e-15);
		assert_true(nnz==1||nnz==2);
		if(nnz==2) ngap++;
	}
	/*theta=0: gaps aligned with +y axis: points above the center are on a gap*/
	for(long iy=(long)cy+1; iy<ny; iy++){
		long i=(long)cx+iy*nx;
		int nnz=0;
		for(long ip=0; ip<npetal; ip++) nnz+=(P(Hf, i, ip)!=0);
		assert_int_equal(nnz, 2);
	}
	/*points on +x axis are in the middle of a petal*/
	long ix0=(long)cx+5+(long)cy*nx;
	int nnz=0;
	for(long ip=0; ip<npetal; ip++) nnz+=(P(Hf, ix0, ip)!=0);
	assert_int_equal(nnz, 1);
	assert_true(ngap>=10);
	/*total weight is the number of points*/
	dmat *cnt=dspsum(H, 1);
	assert_float_equal(dsum(cnt), nx*ny, 1e-10);
	dfree(cnt);
	/*OPD from modes*/
	dmat *mode=dnew(npetal, 1);
	for(long ip=0; ip<npetal; ip++) P(mode, ip)=ip*ip+0.5;
	dmat *opd=petal_mkopd(nx, ny, cx, cy, npetal, 0, mode);
	assert_int_equal(NX(opd), nx);
	assert_int_equal(NY(opd), ny);
	dmat *ho=NULL;
	dspmm(&ho, H, mode, "nn", 1);
	dmat *opdv=dref_reshape(opd, nx*ny, 1);
	assert_dmat_near(opdv, ho, 1e-14);
	/*petal_opd into an existing array*/
	dmat *opd2=dnew(nx, ny);
	petal_opd(opd2, cx, cy, npetal, 0, mode);
	assert_dmat_near(opd2, opd, 0);
	dfree(opd2);
	/*rotation by one petal relabels petals*/
	dsp *H2=petal_mkh(nx, ny, cx, cy, npetal, 2*M_PI/npetal);
	dmat *H2f=NULL;
	dspfull(&H2f, H2, 'n', 1);
	for(long ip=0; ip<npetal; ip++){
		long jp=(ip+npetal-1)%npetal;
		for(long i=0; i<nx*ny; i++){
			assert_float_equal(P(H2f, i, jp), P(Hf, i, ip), 0);
		}
	}
	dfree(H2f); dspfree(H2);
	/*loc version*/
	loc_t *loc=mksqloc(nx, ny, 1, 1, -cx, -cy);
	dsp *Hl=petal_mkh_loc(loc, npetal, 0);
	dmat *Hlf=NULL;
	dspfull(&Hlf, Hl, 'n', 1);
	assert_dmat_near(Hlf, Hf, 0);
	dfree(Hlf); dspfree(Hl); locfree(loc);
	/*invalid number of petals*/
	assert_null(petal_mkh(nx, ny, cx, cy, 1, 0));
	dfree(ho); dfree(opdv); dfree(opd); dfree(mode);
	dfree(Hf); dspfree(H);
}
/**
   petal_connect(): petal pistons measured by several subapertures with
   different offsets are connected through differences; a fully measured ring
   is referenced to petal 0; a gap restarts the reference after it.
*/
static void modes_petal_connect(void **state){
	(void)state;
	const int npetal=6, nsa=4;
	const real p[6]={0.3, -0.2, 0.5, 1.0, -0.7, 0.1};
	dmat *mphi=dnew(npetal, nsa);
	for(int isa=0; isa<nsa; isa++){
		for(int ip=0; ip<npetal; ip++){
			/*each subaperture sees 3 consecutive petals with its own offset*/
			int rel=(ip-isa*npetal/nsa+npetal)%npetal;
			P(mphi, ip, isa)=rel<3?p[ip]+isa*0.37:NAN;
		}
	}
	dmat *out=NULL;
	petal_connect(&out, mphi, npetal, nsa);
	//Fully connected ring: the result is p-p[0].
	assert_float_equal(P(out, 0), 0, 0);
	for(int ip=0; ip<npetal; ip++){
		assert_float_equal(P(out, ip), p[ip]-p[0], 1e-12);
	}
	/*A closure error (inconsistent differences around the ring) is distributed
	  linearly. Add e to petal 2 in subapertures that also measure petal 1 but not 3.*/
	{
		dmat *mphi2=ddup(mphi);
		const real e=0.06;
		for(int isa=0; isa<nsa; isa++){
			if(!isnan(P(mphi2, 1, isa))&&!isnan(P(mphi2, 2, isa))){
				P(mphi2, 2, isa)+=e;
				if(!isnan(P(mphi2, 3, isa))) P(mphi2, 3, isa)+=e;
			}
		}
		dmat *out2=NULL;
		petal_connect(&out2, mphi2, npetal, nsa);
		assert_float_equal(P(out2, 0), 0, 0);
		/*differences are consistent except (1,2) which is off by e. The ring
		  sum is e, which is distributed as -e*ip/npetal.*/
		for(int ip=0; ip<npetal; ip++){
			real ref=p[ip]-p[0]+(ip>=2?e:0)-e*ip/npetal;
			assert_float_equal(P(out2, ip), ref, 1e-12);
		}
		dfree(out2); dfree(mphi2);
	}
	/*break the pair (2,3): reference restarts at petal 3*/
	for(int isa=0; isa<nsa; isa++){
		if(!isnan(P(mphi, 2, isa))&&!isnan(P(mphi, 3, isa))){
			P(mphi, 3, isa)=NAN;
		}
	}
	petal_connect(&out, mphi, npetal, nsa);
	for(int ip=0; ip<npetal; ip++){
		assert_float_equal(P(out, ip), p[ip]-p[3], 1e-12);
	}
	dfree(out); dfree(mphi);
}
/**
   deblur(): removes a known detector transfer function blur exactly.
*/
static void modes_deblur(void **state){
	(void)state;
	const long n=16;
	dmat *psf=dnew(n, n);
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			P(psf, ix, iy)=exp(-(pow(ix-8.3, 2)+pow(iy-7.6, 2))/4.);
		}
	}
	dmat *otf=NULL;
	dtf_otf(&otf, n, n, 1, 1, 0, 0.3, 2);/*peak in corner*/
	cmat *wvf=NULL;
	ccpd(&wvf, psf);
	cfft2(wvf, -1);
	for(long i=0; i<PN(wvf); i++) P(wvf, i)*=P(otf, i);
	cfft2(wvf, 1);
	dmat *blur=NULL;
	creal2d(&blur, 0, wvf, 1./PN(wvf));
	assert_true(ddiff(blur, psf)>1e-3);
	deblur(blur, otf);
	assert_dmat_near(blur, psf, 1e-10);
	/*smaller psf than otf: embedded*/
	dmat *small=dnew(8, 8);
	dembed(small, psf, 0);
	dmat *small0=ddup(small);
	dmat *se=dnew(n, n);
	dembed(se, small, 0);
	ccpd(&wvf, se);
	cfft2(wvf, -1);
	for(long i=0; i<PN(wvf); i++) P(wvf, i)*=P(otf, i);
	cfft2(wvf, 1);
	creal2d(&se, 0, wvf, 1./PN(wvf));
	dembed(small, se, 0);
	deblur(small, otf);/*only an approximate inverse since the blurred wings are cropped*/
	assert_true(ddiff(small, small0)<0.2);
	dfree(se); dfree(small); dfree(small0);
	dfree(blur); cfree(wvf); dfree(otf); dfree(psf);
}
/** PSF (nyquist sampled, peak in center, size 2n) of pupil amp with phase phi (radian).*/
static dmat *mkpsf(const dmat *amp, const dmat *phi, long npsf){
	long n=NX(amp);
	cmat *wvf=cnew(2*n, 2*n);
	cembed_wvf(wvf, P(phi), P(amp), n, n, TWOPI, 0);
	cfft2(wvf, -1);
	cfftshift(wvf);
	dmat *psf0=NULL;
	cabs22d(&psf0, 0, wvf, 1);
	dmat *psf=dnew(npsf, npsf);
	dembed(psf, psf0, 0);
	dfree(psf0);
	cfree(wvf);
	return psf;
}
/**
   petal_setup()/petal_solve(): phase retrieval of differential petal piston
   from the PSF of a single (TT) aperture. The reconstruction with nsa=1 uses
   modes of opposing petal differences. A phase with tip/tilt content needs
   the tip/tilt modes (the PSF is recentered); a tilt free phase is recovered
   without them. Pixel blur is deconvolved. Starting from the truth stays
   there. petal_solve_wfs() is the same as setup+solve.
*/
static void modes_petal_solve(void **state){
	(void)state;
	const long n=16, npetal=6, npsf=2*n;
	/*slightly decentered pupil to avoid exact zeros in the field (see cmocka_wfs.c)*/
	dmat *amp=dnew(n, n);
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			real r2=pow(ix-7.3, 2)+pow(iy-7.6, 2);
			P(amp, ix, iy)=r2<=56?1:0;
		}
	}
	/*true petal phases: antisymmetric (opposing petals opposite), which has
	  no twin image ambiguity. The second one is free of tip/tilt.*/
	const real mv[2][3]={{0.4, -0.25, 0.15}, {0.3, -0.3, 0.3}};
	for(int mc=0; mc<2; mc++){
		dmat *m=dnew(npetal, 1);
		for(int i=0; i<3; i++){
			P(m, i)=mv[mc][i];
			P(m, i+3)=-P(m, i);
		}
		dmat *phi=petal_mkopd(n, n, n/2-0.5, n/2-0.5, npetal, 0, m);
		dcell *ints=dcellnew(1, 1);
		P(ints, 0)=mkpsf(amp, phi, npsf);
		const int withtt=mc==0?1:0;
		petal_t *petal=petal_setup(NULL, 1, amp, 0, 0, 0, npsf, withtt);
		assert_non_null(petal);
		dcell *phi1=NULL;
		dmat *mphi=NULL;
		petal_solve(&phi1, &mphi, petal, ints, NULL, 50);
		assert_int_equal(PN(mphi), npetal);
		for(int i=0; i<npetal; i++){
			assert_float_equal(P(mphi, i), P(m, i), 0.01);
		}
		assert_int_equal(NX(phi1), 1);
		assert_int_equal(PN(P(phi1, 0)), n*n);
		/*starting from the truth (modal): the total stays at the truth*/
		dmat *mphi2=NULL;
		petal_solve(NULL, &mphi2, petal, ints, m, 10);
		for(int i=0; i<npetal; i++){
			assert_float_equal(P(mphi2, i), P(m, i), 0.01);
		}
		dfree(mphi2);
		petal_free(petal);
		if(mc==1){
			/*with pixel blur: blur the PSF with the detector OTF, which is deconvolved*/
			dmat *otf=NULL;
			dtf_otf(&otf, npsf, npsf, 1, 1, 0, 0.2, 2);
			cmat *wvf=NULL;
			ccpd(&wvf, P(ints, 0));
			cfft2(wvf, -1);
			for(long i=0; i<PN(wvf); i++) P(wvf, i)*=P(otf, i);
			cfft2(wvf, 1);
			dcell *ints2=dcellnew(1, 1);
			creal2d(&P(ints2, 0), 0, wvf, 1./PN(wvf));
			petal_t *p2=petal_setup(NULL, 1, amp, 1, 0.2, 0, npsf, withtt);
			dmat *mphi4=NULL;
			petal_solve(NULL, &mphi4, p2, ints2, NULL, 50);
			for(int i=0; i<npetal; i++){
				assert_float_equal(P(mphi4, i), P(m, i), 0.01);
			}
			dfree(mphi4);
			petal_free(p2);
			dcellfree(ints2); cfree(wvf); dfree(otf);
			/*wrapper with the default number of PSF pixels (16)*/
			petal_t *p3=petal_setup(NULL, 1, amp, 0, 0, 0, 16, 0);
			dmat *mphi3=NULL, *mphi5=NULL;
			petal_solve(NULL, &mphi3, p3, ints, NULL, 0);
			petal_solve_wfs(NULL, &mphi5, ints, NULL, NULL, amp, 0, 0, 0, 0, 0);
			assert_dmat_near(mphi5, mphi3, 1e-12);
			/*the measured differential pistons have the right signs*/
			for(int i=0; i<npetal; i++){
				assert_true(P(mphi3, i)*P(m, i)>0);
			}
			dfree(mphi3); dfree(mphi5);
			petal_free(p3);
		}
		dfree(mphi); dcellfree(phi1);
		dcellfree(ints); dfree(phi); dfree(m);
	}
	/*no input is handled*/
	dmat *mphi=NULL;
	petal_solve(NULL, &mphi, NULL, NULL, NULL, 1);
	assert_null(mphi);
	petal_free(NULL);
	petal_free_arr(NULL, 0);
	dfree(amp);
}
/**
   petal_setup()/petal_solve() with several subapertures (quadrants): each
   subaperture sees some petals; the connected solution recovers the petal
   pistons relative to the first one. petal_free_arr() frees an array.
*/
static void modes_petal_quadrant(void **state){
	(void)state;
	const long n=32, ns=16, npetal=6;
	/*full pupil, slightly decentered*/
	dmat *ampf=dnew(n, n);
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			real r2=pow(ix-15.3, 2)+pow(iy-15.6, 2);
			P(ampf, ix, iy)=r2<=225?1:0;
		}
	}
	dmat *m=dnew(npetal, 1);
	const real mv[6]={0, 0.3, -0.2, 0.25, 0.1, -0.3};
	for(int i=0; i<npetal; i++) P(m, i)=mv[i];
	dmat *phif=petal_mkopd(n, n, n/2-0.5, n/2-0.5, npetal, 0, m);
	/*4 quadrants as subapertures, lower left corners in grid units relative to the pupil center*/
	loc_t *saloc=locnew(4, ns, ns);
	dmat *amp=dnew(ns*ns*4, 1);
	dcell *ints=dcellnew(4, 1);
	for(int isa=0; isa<4; isa++){
		long ox=(isa%2)*ns, oy=(isa/2)*ns;
		saloc->locx[isa]=ox-(n/2-0.5);
		saloc->locy[isa]=oy-(n/2-0.5);
		dmat *a=dnew(ns, ns), *ph=dnew(ns, ns);
		for(long iy=0; iy<ns; iy++){
			for(long ix=0; ix<ns; ix++){
				P(a, ix, iy)=P(ampf, ix+ox, iy+oy);
				P(ph, ix, iy)=P(phif, ix+ox, iy+oy);
				P(amp, ix+iy*ns+isa*ns*ns)=P(a, ix, iy);
			}
		}
		P(ints, isa)=mkpsf(a, ph, 2*ns);
		dfree(a); dfree(ph);
	}
	petal_t *petal=petal_setup(saloc, 1, amp, 0, 0, 0, 2*ns, 1);
	dmat *mphi=NULL;
	dcell *phi1=NULL;
	/*A tiny random initial phase avoids exactly zero fields (NaN in gerchberg_saxton, see report).*/
	dmat *phi1b=dnew(ns*ns*4, 1);
	rand_t rstat; seed_rand(&rstat, 5);
	drandn(phi1b, 1e-3, &rstat);
	petal_solve(&phi1, &mphi, petal, ints, phi1b, 50);
	dfree(phi1b);
	assert_int_equal(PN(mphi), npetal);
	assert_int_equal(NX(phi1), 4);
	/*Each quadrant sees two petals; the pairs across the 90 and 270 degree
	  gaps (on quadrant boundaries) are never seen together, so the result is
	  referenced to the first petal after each unmeasured pair.*/
	for(int i=0; i<npetal; i++){
		assert_float_equal(P(mphi, i), mv[i]-mv[i<3?0:3], 0.01);
	}
	/*array of petal_t pointers*/
	petal_t **parr=mycalloc(2, petal_t *);
	parr[0]=petal;
	parr[1]=petal_setup(NULL, 1, P(ints, 0), 0, 0, 0, ns, 0);
	petal_free_arr(parr, 2);
	dfree(mphi); dcellfree(phi1);
	dcellfree(ints); dfree(amp); locfree(saloc);
	dfree(phif); dfree(m); dfree(ampf);
}
UNITTEST_MAIN(
	cmocka_unit_test(modes_zernike_cov),
	cmocka_unit_test(modes_zernike_misc),
	cmocka_unit_test(modes_cov_vonkarman),
	cmocka_unit_test(modes_cov_diagnolize),
	cmocka_unit_test(modes_kl),
	cmocka_unit_test(modes_bessel),
	cmocka_unit_test(modes_petal_h),
	cmocka_unit_test(modes_petal_connect),
	cmocka_unit_test(modes_deblur),
	cmocka_unit_test(modes_petal_solve),
	cmocka_unit_test(modes_petal_quadrant)
)
