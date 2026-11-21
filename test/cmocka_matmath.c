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
void test_embed(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	for(long nx=6; nx<8; nx++){
		for(double ang=0; ang<M_PI/2; ang+=M_PI/3){
			dmat *in=dnew(nx,nx);
			drandn(in, 1, &rstat);
			cmat *cin=cnew(nx,nx);
			cembedd(cin, in, 0);
			//dshow(in, "in");
			for(long i=0; i<10; i++){
				dmat *out=dnew(i, i);
				dembed(out, in, ang);
				//dshow(out, "out");
				if(i>0){
					assert_float_equal(P(in, nx>>1, nx>>1),P(out, i>>1, i>>1),1e-15);
				}
				dfree(out);
				cmat *cout=cnew(i, i);
				cembedd(cout, in, ang);
				if(i>0){
					assert_float_equal(P(in, nx>>1, nx>>1), creal(P(cout, i>>1, i>>1)), 1e-15);
				}
				cmat *cout2=cnew(i, i);
				cembed(cout2, cin, ang);
				assert_float_equal(csumdiffsq(cout, cout2),0, 1e-12);
				cfree(cout);
				cfree(cout2);
			}
			cfree(cin);
			dfree(in);
		}
	}
}
void test_denc(void **state){
	(void)state;
	dmat *screen=(dmat *)genatm_simple(0.186, 0, -11./3., 1./64., 129, 1);
	/*dmat *screen=dnew(129,129);
	rand_t rstat;
	seed_rand(&rstat, 1);
	drandn(screen,1,&rstat);*/
	dmat *xx, *enc;
	real ss1=0, ss2=0, ss3=0;
	xx=dlinspace(0, 1, 60); enc=denc(screen, xx, -1, 10); ss1=P(enc, 10); dfree(enc); dfree(xx);
	xx=dlinspace(0, 1, 70); enc=denc(screen, xx, -1, 10); ss2=P(enc, 10); dfree(enc); dfree(xx);
	dmat *screen2=dnew(128, 128);
	dembed(screen2, screen, 0);
	xx=dlinspace(0, 1, 70); enc=denc(screen2, xx, -1, 10); ss3=P(enc, 10); dfree(enc); dfree(xx);
	//info("ss=%g,%g,%g\n", ss1, ss2, ss3);
	assert_float_equal(ss1, ss2, 1e-9);
	assert_float_equal(ss1, ss3, 1e-9);
	dfree(screen2);
	dfree(screen);
}
void test_dcircle(void **state){
	(void)state;
	//test dcircle, dcog, dshift2center
	long nn=10;
	dmat *A=dnew(nn, nn);
	real r=3;
	real val=2;
	real cx=3;
	real cy=5;
	assert_int_equal(dcircle(NULL, cx, cy, 1, 1, r, val), -1);
	assert_int_equal(dcircle(A, cx, cy, 1, 1, -1, val), -1);
	assert_int_equal(dcircle(A, cx, cy, 1, 1, r, val), 0);
	//dshow(A, "A");
	real r2=sqrt(dsum(A)/val/M_PI);
	assert_float_equal(r, r2, 0.05);
	//test dcog 
	real grad[2];
	dcog(grad, A, 0, 0, 0, 0, 0, NULL);
	assert_float_equal(grad[0], cx-(nn-1)*0.5, 0.1);
	assert_float_equal(grad[1], cy-(nn-1)*0.5, 0.1);
	//test shift2center
	dshift2center(A, 0, 0);
	dcog(grad, A, 0, 0, 0, 0, 0, NULL);
	assert_float_equal(grad[0], 0, 0.1);
	assert_float_equal(grad[1], 0, 0.1);
	dshift2center(A, 0.1, -0.7);
	dcog(grad, A, 0, 0, 0, 0, 0, NULL);
	assert_float_equal(grad[0], 0.1, 0.1);
	assert_float_equal(grad[1], -0.7, 0.1);

	r2=sqrt(dsum(A)/M_PI);
	assert_float_equal(r+1, r2, 0.5);
	info("r=%g, r2=%g, grad=%g %g\n", r, r2, grad[0], grad[1]);
	dfree(A);
}

#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-10
#endif
static dmat *randmat(long nx, long ny, rand_t *rstat){
	dmat *A=dnew(nx, ny);
	drandn(A, 1, rstat);
	return A;
}
/** norm, std, isnan, normalize_sum and normalize_sumabs for several types. */
static void test_norm_std(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *A=randmat(5, 4, &rstat);
	real ss=0, s=0;
	for(long i=0; i<PN(A); i++){ ss+=P(A, i)*P(A, i); s+=P(A, i); }
	assert_float_equal(dnorm(A), sqrt(ss/4), TOL);
	real var=0, mean=s/20;
	for(long i=0; i<PN(A); i++) var+=pow(P(A, i)-mean, 2);
	assert_float_equal(dstd(A), sqrt(var/19), TOL);
	assert_float_equal(dmean(A), mean, TOL);
	assert_float_equal(dmean(NULL), 0, 0);
	assert_float_equal(dnorm(NULL), 0, 0);
	assert_float_equal(dstd(NULL), 0, 0);
	assert_false(disnan(A));
	P(A, 7)=NAN;
	assert_true(disnan(A));
	assert_false(disnan(NULL));
	P(A, 7)=1;
	dnormalize_sum(A, 3);
	assert_float_equal(dsum(A), 3, TOL);
	dnormalize_sumabs(A, 2);
	assert_float_equal(dsumabs(A), 2, TOL);
	dnormalize_sum(NULL, 1);
	dnormalize_sumabs(NULL, 1);
	//complex: std uses |sum|^2
	cmat *C=cnew(3, 3);
	crandn(C, 1, &rstat);
	comp cs=csum(C);
	real css=csumsq(C);
	assert_float_equal(cstd(C), sqrt((css-cabs2(cs)/9)/8), TOL);
	assert_float_equal(cnorm(C), sqrt(css/3), TOL);
	cnormalize_sumabs(C, 1);
	assert_float_equal(csumabs(C), 1, TOL);
	P(C, 2)=COMPLEX(NAN, 0.);
	assert_true(cisnan(C));
	//single precision
	smat *S=snew(4, 1);
	srandn(S, 1, &rstat);
	snormalize_sum(S, 5);
	assert_float_equal(ssum(S), 5, 1e-4);
	real sm=ssum(S)/4., sv=0;
	for(int i=0; i<4; i++) sv+=pow(P(S, i)-sm, 2);
	assert_float_equal(sstd(S), sqrt(sv/3), 1e-4);
	sfree(S);
	cfree(C);
	dfree(A);
}
/** dot, wdot, cwdiv, mulvec, mulvec3, mcc, tmcc and imcc against naive loops. */
static void test_dot_mcc(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	dmat *A=randmat(6, 3, &rstat);
	dmat *B=randmat(6, 3, &rstat);
	real ref=0;
	for(long i=0; i<PN(A); i++) ref+=P(A, i)*P(B, i);
	assert_float_equal(ddot(A, B), ref, TOL);
	dmat *Bv=randmat(18, 1, &rstat);
	assert_float_equal(ddot(A, Bv), 0, 0);//mismatch
	//empty or NULL input gives 0 so that sums over blocks are not corrupted.
	assert_float_equal(ddot(A, NULL), 0, 0);
	assert_float_equal(ddot(NULL, NULL), 0, 0);
	dmat *E=dnew(0, 0);
	assert_float_equal(ddot(E, E), 0, 0);
	dfree(E);
	//wdot: a'*W*b
	dmat *W=randmat(6, 3, &rstat);
	real a[6], b[3];
	for(int i=0; i<6; i++) a[i]=i-2.5;
	for(int i=0; i<3; i++) b[i]=1.5*i+1;
	ref=0;
	for(int j=0; j<3; j++) for(int i=0; i<6; i++) ref+=a[i]*P(W, i, j)*b[j];
	assert_float_equal(dwdot(a, W, b), ref, TOL);
	assert_float_equal(dwdot(NULL, W, b), 0, 0);
	//cwdiv with NaN replacement
	dmat *N=dnew(3, 1), *D=dnew(3, 1);
	P(N, 0)=6; P(N, 1)=0; P(N, 2)=-1;
	P(D, 0)=3; P(D, 1)=0; P(D, 2)=4;
	dcwdiv(N, D, 7);
	assert_float_equal(P(N, 0), 2, 0);
	assert_float_equal(P(N, 1), 7, 0);
	assert_float_equal(P(N, 2), -0.25, 0);
	dcwdiv(N, A, 0);//mismatch: no change
	assert_float_equal(P(N, 0), 2, 0);
	//mulvec: y+=A*x*alpha
	real x[3]={1, -2, 0.5}, y[6]={0}, yref[6]={0};
	for(int i=0; i<6; i++){
		y[i]=yref[i]=i;
		for(int j=0; j<3; j++) yref[i]+=P(A, i, j)*x[j]*2;
	}
	dmulvec(y, A, x, 2);
	for(int i=0; i<6; i++) assert_float_equal(y[i], yref[i], TOL);
	dmulvec(NULL, A, x, 2);
	//mulvec3
	dmat *A3=randmat(3, 3, &rstat);
	real y3[3];
	dmulvec3(y3, A3, x);
	for(int i=0; i<3; i++){
		real r=0;
		for(int j=0; j<3; j++) r+=P(A3, i, j)*x[j];
		assert_float_equal(y3[i], r, TOL);
	}
	y3[0]=99;
	dmulvec3(y3, A, x);//wrong size: no change
	assert_float_equal(y3[0], 99, 0);
	//mcc=A'*diag(wt)*A, tmcc=A*diag(wt)*A'
	dmat *wt=dnew(6, 1);
	for(int i=0; i<6; i++) P(wt, i)=0.5+i;
	dmat *mcc=dmcc(A, wt);
	dmat *wt3=dnew(3, 1);
	for(int i=0; i<3; i++) P(wt3, i)=1+i;
	dmat *tmcc=dtmcc(A, wt3);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			real r=0;
			for(int k=0; k<6; k++) r+=P(A, k, i)*P(wt, k)*P(A, k, j);
			assert_float_equal(P(mcc, i, j), r, TOL);
		}
	}
	for(int i=0; i<6; i++){
		for(int j=0; j<6; j++){
			real r=0;
			for(int k=0; k<3; k++) r+=P(A, i, k)*P(wt3, k)*P(A, j, k);
			assert_float_equal(P(tmcc, i, j), r, TOL);
		}
	}
	assert_null(dmcc(A, wt3));
	assert_null(dtmcc(A, wt));
	dmat *imcc=dimcc(A, wt);
	dmat *Id=NULL;
	dmm(&Id, 0, imcc, mcc, "nn", 1);
	for(int i=0; i<3; i++) for(int j=0; j<3; j++) assert_float_equal(P(Id, i, j), i==j, 1e-8);
	//complex dot is without conjugation
	cmat *Ca=cnew(2, 1), *Cb=cnew(2, 1);
	P(Ca, 0)=COMPLEX(1., 2.); P(Ca, 1)=COMPLEX(0., 1.);
	P(Cb, 0)=COMPLEX(3., -1.); P(Cb, 1)=COMPLEX(2., 0.);
	comp cd=cdot(Ca, Cb);
	assert_float_equal(creal(cd), 5, 0);
	assert_float_equal(cimag(cd), 7, 0);
	comp ca[2]={COMPLEX(1., 1.), COMPLEX(2., 0.)};
	comp cb[1]={COMPLEX(0., 1.)};
	cd=cwdot(ca, Ca, cb);
	assert_float_equal(creal(cd), -5, 1e-14);
	assert_float_equal(cimag(cd), -1, 1e-14);
	cfree(Ca); cfree(Cb);
	dfree(Id); dfree(imcc); dfree(mcc); dfree(tmcc); dfree(wt); dfree(wt3);
	dfree(A3); dfree(N); dfree(D); dfree(W); dfree(Bv); dfree(A); dfree(B);
}
/** rectangle, rotvec, rotvect, rotvecnn, corr and para3. */
static void test_geometry(void **state){
	(void)state;
	dmat *A=dnew(40, 40);
	assert_int_equal(drectangle(A, 20, 20, 6, 3, 0, 1), 0);
	assert_float_equal(dsum(A), 4*6*3, 1e-10);
	assert_float_equal(P(A, 20, 20), 1, 0);
	assert_float_equal(P(A, 27, 20), 0, 0);
	assert_float_equal(P(A, 20, 24), 0, 0);
	dzero(A);
	drectangle(A, 20.5, 20, 6, 3, 0, 1);//half pixel offset
	assert_float_equal(dsum(A), 4*6*3, 1e-3*4*6*3);
	dzero(A);
	drectangle(A, 19.5, 20.3, 8, 4, M_PI/6, 2);//rotated
	//The summed value equals the area times val for rotated rectangles.
	assert_float_equal(dsum(A), 2*4*8*4, 2e-3*2*4*8*4);
	assert_float_equal(P(A, 19+6, 20+3), 2, 0);//inside along the rotated axis
	assert_float_equal(P(A, 19+6, 20-3), 0, 0);
	{//more random rotated rectangles
		rand_t rs;
		seed_rand(&rs, 7);
		for(int i=0; i<20; i++){
			real rx=2+8*randu(&rs), ry=2+8*randu(&rs), th=randu(&rs)*M_PI;
			dzero(A);
			drectangle(A, 19.5+randu(&rs), 19.5+randu(&rs), rx, ry, th, 1);
			assert_float_equal(dsum(A), 4*rx*ry, 1e-2*4*rx*ry);
		}
	}
	assert_int_equal(drectangle(A, 0, 0, -1, 1, 0, 1), -1);
	assert_int_equal(drectangle(NULL, 0, 0, 1, 1, 0, 1), -1);
	dfree(A);
	//rotvec rotates each row (x,y)
	dmat *V=dnew(3, 2);
	P(V, 0, 0)=1; P(V, 0, 1)=0;
	P(V, 1, 0)=0; P(V, 1, 1)=2;
	P(V, 2, 0)=1; P(V, 2, 1)=1;
	drotvec(V, M_PI/2);
	assert_float_equal(P(V, 0, 0), 0, 1e-15);
	assert_float_equal(P(V, 0, 1), 1, 1e-15);
	assert_float_equal(P(V, 1, 0), -2, 1e-15);
	assert_float_equal(P(V, 2, 0), -1, 1e-15);
	assert_float_equal(P(V, 2, 1), 1, 1e-15);
	dmat *Vt=dtrans(V);
	drotvect(Vt, -M_PI/2);//back
	assert_float_equal(P(Vt, 0, 0), 1, 1e-15);
	assert_float_equal(P(Vt, 1, 1), 2, 1e-15);
	drotvect(V, 1);//wrong shape: no change
	drotvec(Vt, 1);
	assert_float_equal(P(V, 1, 0), -2, 1e-15);
	dfree(V); dfree(Vt);
	//rotvecnn: B=R*A*R'
	dmat *M=dnew(2, 2);
	P(M, 0, 0)=3; P(M, 1, 0)=1; P(M, 0, 1)=1; P(M, 1, 1)=2;
	dmat *Mr=NULL;
	real th=0.4;
	drotvecnn(&Mr, M, th);
	dmat *R=dnew(2, 2);
	P(R, 0, 0)=cos(th); P(R, 0, 1)=-sin(th); P(R, 1, 0)=sin(th); P(R, 1, 1)=cos(th);
	dmat *RM=NULL, *RMR=NULL;
	dmm(&RM, 0, R, M, "nn", 1);
	dmm(&RMR, 0, RM, R, "nt", 1);
	assert_dmat_near(Mr, RMR, TOL);
	dmat *M3=dnew(3, 3);
	drotvecnn(&Mr, M3, th);//wrong size: no change
	assert_dmat_near(Mr, RMR, TOL);
	dfree(M3); dfree(M); dfree(Mr); dfree(R); dfree(RM); dfree(RMR);

	//corr and para3 find the shift between two images
	const int n=16;
	dmat *g1=dnew(n, n), *g2=dnew(n, n);
	for(int iy=0; iy<n; iy++){
		for(int ix=0; ix<n; ix++){
			P(g1, ix, iy)=exp(-(pow(ix-7.5, 2)+pow(iy-7.5, 2))/4);
			P(g2, ix, iy)=exp(-(pow(ix-7.5-2, 2)+pow(iy-7.5+1, 2))/4);
		}
	}
	dmat *cc=NULL;
	dcorr(&cc, g2, g1);
	assert_int_equal(NX(cc), 2*n-1);
	//direct check of one entry
	real r=0;
	for(int iy=0; iy<n; iy++) for(int ix=0; ix<n; ix++){
		int jx=ix-2, jy=iy-1;
		if(jx>=0&&jx<n&&jy>=0&&jy<n) r+=P(g2, ix, iy)*P(g1, jx, jy);
	}
	assert_float_equal(P(cc, n-1+2, n-1+1), r, TOL);
	real grad[2];
	dpara3(grad, cc);
	assert_float_equal(grad[0], 2, 0.05);
	assert_float_equal(grad[1], -1, 0.05);
	dcorr(&cc, g2, g1);//reuse output
	dmat *g3=dnew(3, 3);
	dcorr(&cc, g2, g3);//mismatch: no change
	dfree(g3);
	dfree(cc); dfree(g1); dfree(g2);
}
/** cog with mask, flux, threshold and shift2center. */
static void test_cog(void **state){
	(void)state;
	dmat *A=dnew(10, 8);
	P(A, 6, 2)=3; P(A, 7, 2)=1;
	real grad[2];
	dcog(grad, A, 0, 0, 0, 0, 0, NULL);
	assert_float_equal(grad[0], (6*3+7)/4.-4.5, 1e-14);
	assert_float_equal(grad[1], 2-3.5, 1e-14);
	//offset of the reference
	dcog(grad, A, 1, -1, 0, 0, 0, NULL);
	assert_float_equal(grad[0], (6*3+7)/4.-5.5, 1e-14);
	assert_float_equal(grad[1], 2-2.5, 1e-14);
	//mask selects pixels; threshold is ignored
	dmat *mask=dnew(10, 8);
	P(mask, 7, 2)=1;
	dcog(grad, A, 0, 0, 100, 0, 0, mask);
	assert_float_equal(grad[0], 7-4.5, 1e-14);
	//flux overrides the denominator
	dcog(grad, A, 0, 0, 0, 0, 8, NULL);
	assert_float_equal(grad[0], (6*3+7)/8.-4.5, 1e-14);
	//background and threshold
	dadds(A, 0.5);
	dcog(grad, A, 0, 0, 0.6, 0.5, 0, NULL);
	assert_float_equal(grad[0], (6*3+7)/4.-4.5, 1e-12);
	dcog(grad, A, 0, 0, 10, 0.5, 0, NULL);//all below threshold
	assert_float_equal(grad[0], 0, 0);
	assert_float_equal(grad[1], 0, 0);
	dfree(mask);
	dfree(A);
	//shift2center on a complex array
	cmat *C=cnew(16, 16);
	for(int iy=0; iy<16; iy++) for(int ix=0; ix<16; ix++) P(C, ix, iy)=COMPLEX(exp(-(pow(ix-5, 2)+pow(iy-6, 2))/4.), 0.);
	cshift2center(C, 0, 0);
	real gc[2];
	ccog(gc, C, 0, 0, 0, 0, 0, NULL);
	assert_float_equal(gc[0], 0, 0.15);
	assert_float_equal(gc[1], 0, 0.15);
	cfree(C);
}
/** gramschmidt orthonormalizes columns with and without weighting. */
static void test_gramschmidt(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 3);
	dmat *M=randmat(20, 5, &rstat);
	dmat *M0=ddup(M);
	real amp[20];
	for(int i=0; i<20; i++) amp[i]=0.5+(i%3);
	dgramschmidt(M, amp);
	real wsum=dvecsum(amp, 20);
	for(int i=0; i<5; i++){
		for(int j=0; j<5; j++){
			real d=dvecdot(PCOL(M, i), PCOL(M, j), amp, 20)/wsum;
			assert_float_equal(d, i==j?1:0, 1e-10);
		}
	}
	//first column keeps its direction
	real ratio=P(M, 0, 0)/P(M0, 0, 0);
	for(int i=0; i<20; i++) assert_float_equal(P(M, i, 0), P(M0, i, 0)*ratio, 1e-10);
	//dependent column is detected and left as zero
	dcp(&M, M0);
	for(int i=0; i<20; i++) P(M, i, 3)=P(M, i, 0)*2-P(M, i, 1);
	dgramschmidt(M, NULL);
	dmat *M3=drefcols(M, 3, 1);
	assert_float_equal(dnorm(M3), 0, 1e-6);
	dfree(M3);
	for(int i=0; i<5; i++){
		if(i==3) continue;
		real d=dvecdot(PCOL(M, i), PCOL(M, 4), NULL, 20)/20;
		assert_float_equal(d, i==4?1:0, 1e-10);
	}
	//complex
	cmat *C=cnew(10, 3);
	crandn(C, 1, &rstat);
	cgramschmidt(C, NULL);
	for(int i=0; i<3; i++){//orthonormal with the conjugate inner product
		for(int j=0; j<3; j++){
			comp d=0;
			for(int k=0; k<10; k++) d+=conj(P(C, k, i))*P(C, k, j);
			assert_float_equal(cabs(d/10.-(i==j?1.:0.)), 0, 1e-10);
		}
	}
	cfree(C);
	dfree(M0);
	dfree(M);
}
/** cwpow_thres, polyval, addI, add, addcol, add_relax, adds, logspace, linspace, cwlog, cwlog10. */
static void test_elementwise(void **state){
	(void)state;
	dmat *A=dnew(4, 1);
	P(A, 0)=4; P(A, 1)=-0.1; P(A, 2)=0.5; P(A, 3)=9;
	dcwpow_thres(A, 0.5, 0.05);//thres=0.45
	assert_float_equal(P(A, 0), 2, 0);
	assert_float_equal(P(A, 1), 0, 0);
	assert_float_equal(P(A, 2), sqrt(0.5), 1e-15);
	assert_float_equal(P(A, 3), 3, 0);
	//polyval: p(0)*x^2+p(1)*x+p(2) with p given as a row vector
	dmat *x=dlinspace(-1, 0.5, 5);
	dmat *p=dnew(1, 3);
	P(p, 0)=2; P(p, 1)=-1; P(p, 2)=0.5;
	dmat *y=ddup(x);
	dpolyval(y, p);
	assert_int_equal(NX(p), 3);
	for(int i=0; i<5; i++){
		real xi=P(x, i);
		assert_float_equal(P(y, i), 2*xi*xi-xi+0.5, 1e-14);
	}
	//addI on a non-square matrix
	dmat *B=dnew(2, 3);
	daddI(B, 2);
	assert_float_equal(dsum(B), 4, 0);
	assert_float_equal(P(B, 1, 1), 2, 0);
	daddI(NULL, 1);
	//add with only scaling of B
	dmat *B2=dnew(2, 3);
	dadd(&B, 3, B2, 0);
	assert_float_equal(P(B, 1, 1), 6, 0);
	dmat *Bc=ddup(B);
	dadd(&B, 1, B2, 0);//ac=0 and bc=1: no change
	assert_dmat_near(B, Bc, 0);
	dmat *C=NULL;
	dadd(&C, 5, B, 2);//C is created: bc ignored
	assert_float_equal(P(C, 1, 1), 12, 0);
	dadd(&C, 0, B, -1);//assign
	assert_float_equal(P(C, 1, 1), -6, 0);
	dadd(&C, 0, B, 1);//copy
	assert_float_equal(P(C, 1, 1), 6, 0);
	//addcol
	dmat *v=dnew(2, 1);
	P(v, 0)=1; P(v, 1)=2;
	daddcol(B, 2, 2, v, 3);//B(:,2)=2*B(:,2)+3*v
	assert_float_equal(P(B, 0, 2), 3, 0);
	assert_float_equal(P(B, 1, 2), 6, 0);
	daddcol(B, 0, 0, v, 1);//copy
	assert_float_equal(P(B, 0, 0), 1, 0);
	assert_float_equal(P(B, 1, 0), 2, 0);
	daddcol(B, 1, 0, v, -2);//assign scaled
	assert_float_equal(P(B, 1, 1), -4, 0);
	daddcol(B, 1, 0, v, 0);//no-op
	assert_float_equal(P(B, 1, 1), -4, 0);
	//add_relax adds the overlapping part only
	dmat *D=dnew(3, 2);
	dset(D, 1);
	dadd_relax(&D, 2, B, 1);
	assert_float_equal(P(D, 0, 0), 2+1, 0);
	assert_float_equal(P(D, 1, 1), 2-4, 0);
	assert_float_equal(P(D, 2, 1), 1, 0);//outside of B: unchanged
	dmat *E=NULL;
	dadd_relax(&E, 2, B, 3);
	assert_float_equal(P(E, 1, 1), -12, 0);
	dadd_relax(&E, 2, NULL, 3);//no-op
	assert_float_equal(P(E, 1, 1), -12, 0);
	//adds
	dadds(E, 1.5);
	assert_float_equal(P(E, 1, 1), -10.5, 0);
	dadds(NULL, 1);
	//logspace and linspace
	dmat *ls=dlogspace(-1, 2, 4);
	assert_float_equal(P(ls, 0), 0.1, 1e-15);
	assert_float_equal(P(ls, 3), 100, 1e-12);
	dcwlog10(ls);
	for(int i=0; i<4; i++) assert_float_equal(P(ls, i), i-1, 1e-14);
	dmat *ex=dlinspace(0, 1, 3);
	dcwexp(ex, 1);
	dcwlog(ex);
	for(int i=0; i<3; i++) assert_float_equal(P(ex, i), i, 1e-14);
	smat *sl=slogspace(0, 1, 2);
	assert_float_equal(P(sl, 1), 10, 1e-5);
	scwlog(sl);
	assert_float_equal(P(sl, 1), log(10), 1e-5);
	sfree(sl);
	dfree(ls); dfree(ex);
	dfree(E); dfree(D); dfree(B2); dfree(v); dfree(C); dfree(Bc); dfree(B);
	dfree(y); dfree(p); dfree(x); dfree(A);
}
/** interp1 with linear, logarithmic and arbitrary spacing, and with two column input. */
static void test_interp1(void **state){
	(void)state;
	//log spacing
	dmat *xin=dlogspace(0, 2, 11);
	dmat *yin=dnew(11, 2);
	for(int i=0; i<11; i++){
		P(yin, i, 0)=log10(P(xin, i))*2+1;//linear in log(x)
		P(yin, i, 1)=-log10(P(xin, i));
	}
	dmat *xnew=dnew(5, 1);
	P(xnew, 0)=0.5; P(xnew, 1)=1; P(xnew, 2)=5; P(xnew, 3)=100; P(xnew, 4)=200;
	dmat *ynew=dinterp1(xin, yin, xnew, NAN);
	assert_int_equal(NY(ynew), 2);
	assert_float_equal(P(ynew, 0, 0), 1, 1e-12);//below range: first value
	assert_float_equal(P(ynew, 1, 0), 1, 1e-12);
	assert_float_equal(P(ynew, 2, 0), log10(5)*2+1, 1e-12);
	assert_float_equal(P(ynew, 3, 0), 5, 1e-12);
	assert_float_equal(P(ynew, 4, 0), 5, 1e-12);//above range: last value
	assert_float_equal(P(ynew, 2, 1), -log10(5), 1e-12);
	dfree(ynew);
	ynew=dinterp1log(xin, yin, xnew, -1);
	assert_float_equal(P(ynew, 0, 0), -1, 0);
	assert_float_equal(P(ynew, 4, 1), -1, 0);
	assert_float_equal(P(ynew, 3, 1), -2, 1e-12);
	dfree(ynew);
	//arbitrary spacing
	dmat *xa=dnew(4, 1), *ya=dnew(4, 1);
	P(xa, 0)=0; P(xa, 1)=1; P(xa, 2)=3; P(xa, 3)=7;
	for(int i=0; i<4; i++) P(ya, i)=2*P(xa, i)-1;
	dmat *xn=dnew(5, 1);
	P(xn, 0)=-1; P(xn, 1)=0.5; P(xn, 2)=2; P(xn, 3)=6; P(xn, 4)=9;
	ynew=dinterp1(xa, ya, xn, NAN);//extrapolate linearly
	for(int i=0; i<5; i++) assert_float_equal(P(ynew, i), 2*P(xn, i)-1, 1e-12);
	dfree(ynew);
	ynew=dinterp1(xa, ya, xn, 100);//default value outside
	assert_float_equal(P(ynew, 0), 100, 0);
	assert_float_equal(P(ynew, 4), 100, 0);
	assert_float_equal(P(ynew, 2), 3, 1e-12);
	dfree(ynew);
	//two column input: second column is the value
	dmat *xy=dcat(xa, ya, 2);
	ynew=dinterp1(xy, NULL, xn, NAN);
	assert_float_equal(P(ynew, 3), 11, 1e-12);
	dfree(ynew);
	//several columns: columns 2 and later are the values; linear spacing uses the first column of xin
	dmat *xyy=dcat(xy, ya, 2);
	ynew=dinterp1(xyy, NULL, xn, NAN);
	assert_int_equal(NY(ynew), 2);
	assert_float_equal(P(ynew, 3, 1), 11, 1e-12);
	dfree(ynew);
	dmat *xl2=dlinspace(0, 1, 4);
	dmat *xly=dcat(xl2, xa, 2);
	ynew=dinterp1linear(xly, NULL, xn, NAN);
	assert_float_equal(P(ynew, 2), 3, 1e-12);//xa at x=2
	dfree(ynew);
	assert_null(dinterp1linear(xl2, NULL, xn, 0));//no values
	ynew=dinterp1(xa, ya, xy, 0);//only the first column of xnew is used
	assert_int_equal(NX(ynew), 4); assert_int_equal(NY(ynew), 1);
	assert_float_equal(P(ynew, 2), 5, 1e-12);
	dfree(ynew);
	dmat *xp=dnew(1, 1);
	assert_null(dinterp1(xp, xp, xn, 0));//single point; xp is not freed
	assert_float_equal(P(xp, 0), 0, 0);
	dfree(xp);
	dfree(xly); dfree(xl2); dfree(xyy);
	//linear spacing with single precision
	smat *sx=slinspace(0, 0.5, 5), *sy=slinspace(1, 1, 5), *sn=snew(2, 1);
	P(sn, 0)=0.25; P(sn, 1)=1.75;
	smat *syn=sinterp1(sx, sy, sn, 0);
	assert_float_equal(P(syn, 0), 1.5, 1e-6);
	assert_float_equal(P(syn, 1), 4.5, 1e-6);
	sfree(syn); sfree(sx); sfree(sy); sfree(sn);
	//monotonically decreasing coordinates: linear, logarithmic, and arbitrary spacing
	{
		dmat *xd=dlinspace(4, -0.5, 9);//4, 3.5, ..., 0
		dmat *yd=dnew(9, 1);
		for(int i=0; i<9; i++) P(yd, i)=2*P(xd, i)+1;
		dmat *xq=dnew(4, 1);
		P(xq, 0)=1.25; P(xq, 1)=3.9; P(xq, 2)=-1; P(xq, 3)=5;
		ynew=dinterp1linear(xd, yd, xq, -7);
		assert_non_null(ynew);
		assert_float_equal(P(ynew, 0), 3.5, 1e-12);
		assert_float_equal(P(ynew, 1), 8.8, 1e-12);
		assert_float_equal(P(ynew, 2), -7, 0);//outside
		assert_float_equal(P(ynew, 3), -7, 0);
		dfree(ynew);
		ynew=dinterp1(xd, yd, xq, NAN);//nearest end value outside
		assert_float_equal(P(ynew, 0), 3.5, 1e-12);
		assert_float_equal(P(ynew, 2), 1, 1e-12);
		assert_float_equal(P(ynew, 3), 9, 1e-12);
		dfree(ynew);
		//log spacing, decreasing
		dmat *xl=dnew(11, 1), *yl=dnew(11, 1);
		for(int i=0; i<11; i++){
			P(xl, i)=P(xin, 10-i);
			P(yl, i)=log10(P(xl, i))*2+1;
		}
		ynew=dinterp1log(xl, yl, xnew, -1);
		assert_non_null(ynew);
		assert_float_equal(P(ynew, 0), -1, 0);
		assert_float_equal(P(ynew, 2), log10(5)*2+1, 1e-12);
		assert_float_equal(P(ynew, 3), 5, 1e-12);
		dfree(ynew);
		ynew=dinterp1(xl, yl, xnew, NAN);
		assert_float_equal(P(ynew, 2), log10(5)*2+1, 1e-12);
		dfree(ynew);
		//arbitrary spacing, decreasing
		dmat *xr=dnew(4, 1), *yr=dnew(4, 1);
		for(int i=0; i<4; i++){
			P(xr, i)=P(xa, 3-i);
			P(yr, i)=P(ya, 3-i);
		}
		ynew=dinterp1(xr, yr, xn, NAN);//extrapolate linearly
		for(int i=0; i<5; i++) assert_float_equal(P(ynew, i), 2*P(xn, i)-1, 1e-12);
		dfree(ynew);
		ynew=dinterp1(xr, yr, xn, 100);
		assert_float_equal(P(ynew, 0), 100, 0);
		assert_float_equal(P(ynew, 4), 100, 0);
		assert_float_equal(P(ynew, 2), 3, 1e-12);
		dfree(ynew);
		//invalid input returns NULL
		assert_null(dinterp1linear(xr, yr, xn, 0));//not linearly spaced
		assert_null(dinterp1log(xd, yd, xq, 0));//not log spaced
		assert_null(dinterp1(NULL, yd, xq, 0));
		dmat *x1=dnew(1, 1);
		assert_null(dinterp1linear(x1, x1, xq, 0));//single point
		dfree(x1);
		dfree(xr); dfree(yr); dfree(xl); dfree(yl); dfree(xd); dfree(yd); dfree(xq);
	}
	dfree(xy); dfree(xn); dfree(xa); dfree(ya);
	dfree(xnew); dfree(yin); dfree(xin);
}
/** blend merges B into the center of A with a soft edge. */
static void test_blend(void **state){
	(void)state;
	dmat *A=dnew(10, 10);
	dmat *B=dnew(6, 6);
	dset(B, 1);
	dblend(A, B, 2);
	//center of B is copied, edges are weighted
	assert_float_equal(P(A, 5, 5), 1, 0);
	assert_float_equal(P(A, 2, 5), 0, 0);//weight 0 at the edge
	assert_float_equal(P(A, 3, 5), 1, 0);
	assert_float_equal(P(A, 0, 0), 0, 0);
	//B larger than A: center part of B is used
	dmat *C=dnew(4, 4);
	dmat *D=dnew(8, 8);
	for(int i=0; i<64; i++) P(D, i)=i;
	dblend(C, D, 2);
	assert_float_equal(P(C, 1, 1), P(D, 3, 3), 0);
	assert_float_equal(P(C, 2, 2), P(D, 4, 4), 0);
	assert_float_equal(P(C, 3, 3), 0, 0);
	assert_float_equal(P(C, 0, 1), 0, 0);
	dfree(C); dfree(D);
	dfree(A); dfree(B);
}
/** 1-d and 2-d cubic splines reproduce the data at nodes and linear functions everywhere. */
static void test_spline(void **state){
	(void)state;
	dmat *x=dlinspace(0, 0.25, 17);
	dmat *y=dnew(17, 1);
	for(int i=0; i<17; i++) P(y, i)=sin(P(x, i));
	dmat *ys=dspline(x, y, x);
	assert_dmat_near(ys, y, 1e-12);
	dfree(ys);
	dmat *xn=dlinspace(0.1, 0.37, 10);
	ys=dspline(x, y, xn);
	for(int i=0; i<10; i++) assert_float_equal(P(ys, i), sin(P(xn, i)), P(xn, i)<0.25?2e-3:3e-4);
	dfree(ys);
	//two column input with a linear function
	dmat *xy=dnew(17, 2);
	for(int i=0; i<17; i++){ P(xy, i, 0)=P(x, i); P(xy, i, 1)=3*P(x, i)-2; }
	dmat *coeff=dspline_prep(xy, NULL);
	ys=dspline_eval(coeff, x, xn);
	for(int i=0; i<10; i++) assert_float_equal(P(ys, i), 3*P(xn, i)-2, 1e-12);
	dfree(ys); dfree(coeff); dfree(xy);
	//2-d spline of a bilinear function
	dmat *gx=dlinspace(-1, 0.5, 6), *gy=dlinspace(0, 0.25, 5);
	dmat *z=dnew(6, 5);
	for(int iy=0; iy<5; iy++) for(int ix=0; ix<6; ix++) P(z, ix, iy)=1+2*P(gx, ix)-P(gy, iy)+P(gx, ix)*P(gy, iy);
	dcell *bc=dbspline_prep(gx, gy, z);
	dmat *px=dnew(3, 1), *py=dnew(3, 1);
	P(px, 0)=-1; P(py, 0)=0;
	P(px, 1)=0.3; P(py, 1)=0.6;
	P(px, 2)=1.4; P(py, 2)=0.95;
	dmat *pz=dbspline_eval(bc, gx, gy, px, py);
	for(int i=0; i<3; i++){
		real xx=P(px, i), yy=P(py, i);
		assert_float_equal(P(pz, i), 1+2*xx-yy+xx*yy, 1e-10);
	}
	dfree(pz); dfree(px); dfree(py); dcellfree(bc);
	dfree(z); dfree(gx); dfree(gy);
	dfree(xn); dfree(x); dfree(y);
}
/** fwhm, gauss_fit, gauss_sigma and fwhm_gauss of an elliptical Gaussian. */
static void test_gauss(void **state){
	(void)state;
	const int n=64;
	const real sa=4, sb=2.5, th=0.3, amp=2, cx=31.2, cy=32.4;
	dmat *A=dnew(n, n);
	const real ct=cos(th), st=sin(th);
	for(int iy=0; iy<n; iy++){
		for(int ix=0; ix<n; ix++){
			real x=ix-cx, y=iy-cy;
			real u=x*ct+y*st, v=-x*st+y*ct;
			P(A, ix, iy)=amp*exp(-0.5*(u*u/(sa*sa)+v*v/(sb*sb)));
		}
	}
	real mr, ma, mb, angle, a;
	dmat *img=NULL;
	dgauss_fit(&mr, &ma, &mb, &angle, &a, &img, A, 0);
	assert_float_equal(ma, sa, 1e-3);
	assert_float_equal(mb, sb, 1e-3);
	assert_float_equal(mr, sqrt(sa*sb), 1e-3);
	assert_float_equal(angle, th, 1e-3);
	assert_float_equal(a, amp, 1e-3);
	assert_dmat_near(img, A, 1e-3);
	assert_float_equal(dgauss_sigma(A, 0), mr, 1e-12);
	//fwhm_gauss uses a threshold of 1% which truncates the wings slightly
	real mr1=dgauss_sigma(A, 0.01);
	assert_true(mr1<mr && mr1>0.95*mr);
	assert_float_equal(dfwhm_gauss(A), mr1*2*1.17741, 1e-12);
	dfree(img);
	//empty image
	dmat *Z=dnew(4, 4);
	dgauss_fit(&mr, &ma, &mb, &angle, &a, NULL, Z, 0);
	assert_float_equal(mr, 0, 0);
	assert_float_equal(a, 0, 0);
	dfree(Z);
	//fwhm counts pixels above half maximum: exact for a circular Gaussian
	dmat *G=dnew(n, n);
	const real sg=5;
	for(int iy=0; iy<n; iy++) for(int ix=0; ix<n; ix++) P(G, ix, iy)=exp(-0.5*(pow(ix-32, 2)+pow(iy-32, 2))/(sg*sg));
	assert_float_equal(dfwhm(G), 2.35482*sg, 0.3);
	assert_float_equal(dfwhm(NULL), 0, 0);
	dfree(G);
	dfree(A);
}
/** enc: ensquared, encircled, enslitted energy and azimuthal average of a uniform patch. */
static void test_enc_types(void **state){
	(void)state;
	//uniform 1/(32*32) image: enclosed energy is proportional to the area
	const int n=32;
	dmat *psf=dnew(n, n);
	dset(psf, 1./(n*n));
	dmat *dvec=dnew(3, 1);
	P(dvec, 0)=4; P(dvec, 1)=8; P(dvec, 2)=16;
	dmat *e0=denc(psf, dvec, 0, 2);
	dmat *e1=denc(psf, dvec, 1, 1);
	dmat *e2=denc(psf, dvec, 2, 2);
	dmat *em=denc(psf, dvec, -1, 2);
	for(int i=0; i<3; i++){
		real d=P(dvec, i);
		assert_float_equal(P(e0, i), d*d/(n*n), 0.02);
		assert_float_equal(P(e1, i), M_PI*d*d/4/(n*n), 0.02);
		assert_float_equal(P(e2, i), d/n, 0.02);
		assert_float_equal(P(em, i), 1./(n*n), 0.2/(n*n));//radius within the patch
	}
	assert_null(denc(psf, dvec, 3, 1));//invalid type
	//larger psf is cropped to the needed size
	dmat *big=dnew(128, 128);
	dembed(big, psf, 0);
	dmat *e0b=denc(big, dvec, 0, 0);
	assert_dmat_near(e0b, e0, 1e-3);
	dfree(e0b); dfree(big);
	dfree(e0); dfree(e1); dfree(e2); dfree(em);
	dfree(dvec); dfree(psf);
}
/** trapz with shared, per column and no abscissa. */
static void test_trapz(void **state){
	(void)state;
	dmat *x=dlinspace(0, 0.25, 9);//0..2
	dmat *y=dnew(9, 2);
	for(int i=0; i<9; i++){
		P(y, i, 0)=P(x, i);//integral 2
		P(y, i, 1)=1;//integral 2
	}
	assert_float_equal(dtrapz(x, y), 4, 1e-14);
	dmat *x2=dnew(9, 2);
	for(int i=0; i<9; i++){ P(x2, i, 0)=P(x, i); P(x2, i, 1)=-2*P(x, i); }//abs of spacing is used
	assert_float_equal(dtrapz(x2, y), 2+4, 1e-14);
	assert_float_equal(dtrapz(NULL, y), 16, 1e-14);//unit spacing
	assert_float_equal(dtrapz(x, NULL), 0, 0);
	cmat *cy=cnew(3, 1);
	P(cy, 0)=COMPLEX(0., 1.); P(cy, 1)=COMPLEX(1., 1.); P(cy, 2)=COMPLEX(2., 1.);
	comp ct=ctrapz(NULL, cy);
	assert_float_equal(creal(ct), 2, 1e-14);
	assert_float_equal(cimag(ct), 2, 1e-14);
	cfree(cy);
	dfree(x2); dfree(x); dfree(y);
}

UNITTEST_MAIN(
	cmocka_unit_test(test_embed),
	cmocka_unit_test(test_denc),
	cmocka_unit_test(test_dcircle),
	cmocka_unit_test(test_norm_std),
	cmocka_unit_test(test_dot_mcc),
	cmocka_unit_test(test_geometry),
	cmocka_unit_test(test_cog),
	cmocka_unit_test(test_gramschmidt),
	cmocka_unit_test(test_elementwise),
	cmocka_unit_test(test_interp1),
	cmocka_unit_test(test_blend),
	cmocka_unit_test(test_spline),
	cmocka_unit_test(test_gauss),
	cmocka_unit_test(test_enc_types),
	cmocka_unit_test(test_trapz),
)
