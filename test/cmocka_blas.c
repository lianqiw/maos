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

static void blas_dpinv(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *A=dnew(10, 4);
	drandn(A, 1, &rstat);
	dmat *w=dnew(10, 1);
	dset(w, 2);
	dmat *Ap=dpinv(A, w);
	dmat *ApA=NULL;
	dmm(&ApA, 0, Ap, A, "nn", 1);
	dsp *spw=dspnewdiag(w->nx, P(w), 1);
	dmat *Ap2=dpinv(A, spw);
	assert_float_equal(dsum(Ap), 0.675211, 1e-6);
	assert_float_equal(ddiff(Ap, Ap2), 0, 1e-6);

	dcell *At=dcellnew_same(2,2,5,2);
	dcell *Ac=NULL;
	d2cell(&Ac, A, At);
	dcell *Acp=dcellpinv(Ac, NULL);
	dmat *Acp2=dcell2m(Acp);
	assert_float_equal(ddiff(Acp2, Ap2), 0, 1e-6);

	dfree(A); dfree(w); dfree(Ap); dfree(ApA); dspfree(spw); dfree(Ap2);
	dcellfree(At);
	dcellfree(Ac);
	dcellfree(Acp);
	dfree(Acp2);
}

#define ID(a) (a)
/*Naive op(A)(i,k) for trans character t*/
#define OPA(A, t, i, k, CONJF) ((t)=='n'?P(A, i, k):((t)=='t'?P(A, k, i):CONJF(P(A, k, i))))
/*Compare X(mm) for all transposes against a naive triple loop. A is m*k, B is k*n before transposing.*/
#define TEST_MM(X, T, CONJF, ABSF, TOL, NTRANS)				\
	{								\
		rand_t rstat;						\
		seed_rand(&rstat, 5);					\
		const char *trans[]={"nn", "nt", "tn", "tt", "cn", "nc", "cc", "ct", "tc"}; \
		for(int it=0; it<NTRANS; it++){				\
			const char *tr=trans[it];			\
			const long m=5, k=4, n=3;			\
			X##mat *A=tr[0]=='n'?X##new(m, k):X##new(k, m);	\
			X##mat *B=tr[1]=='n'?X##new(k, n):X##new(n, k);	\
			X##randn(A, 1, &rstat);				\
			X##randn(B, 1, &rstat);				\
			X##mat *C=X##new(m, n);				\
			X##randn(C, 1, &rstat);				\
			X##mat *C0=X##dup(C);				\
			X##mm(&C, (T)0.5, A, B, tr, (T)2);		\
			for(long j=0; j<n; j++){			\
				for(long i=0; i<m; i++){		\
					T ref=P(C0, i, j)*(T)0.5;	\
					for(long l=0; l<k; l++){	\
						ref+=(T)2*OPA(A, tr[0], i, l, CONJF)*OPA(B, tr[1], l, j, CONJF); \
					}				\
					assert_float_equal(ABSF(P(C, i, j)-ref), 0, TOL); \
				}					\
			}						\
			/*mm3: op(A)*op(B)*op(C) with C=identity*/	\
			X##mat *Id=X##new(n, n);			\
			X##addI(Id, (T)1);				\
			char tr3[4]={tr[0], tr[1], 'n', 0};		\
			X##mat *D=NULL;					\
			X##mm3(&D, (T)0, A, B, Id, tr3, (T)1);		\
			X##mat *D2=NULL;				\
			X##mm(&D2, (T)0, A, B, tr, (T)1);		\
			for(long i=0; i<PN(D); i++){			\
				assert_float_equal(ABSF(P(D, i)-P(D2, i)), 0, TOL); \
			}						\
			X##free(D); X##free(D2); X##free(Id);		\
			X##free(A); X##free(B); X##free(C); X##free(C0); \
		}							\
		X##mat *E=NULL;						\
		X##mm(&E, (T)0, NULL, NULL, "nn", (T)1);/*no-op*/	\
		assert_null(E);						\
	}
static float abssf(float a){ return fabsf(a); }
static float abszf(fcomplex a){ return cabsf(a); }
static real absc(comp a){ return cabs(a); }
static fcomplex conjzf(fcomplex a){ return conjf(a); }
static comp conjc(comp a){ return conj(a); }
/** X(mm) and X(mm3) for d, s, c, z types and all transposes against naive products. */
static void blas_mm(void **state){
	(void)state;
	TEST_MM(d, real, ID, fabs, 1e-12, 4);
	TEST_MM(s, float, ID, abssf, 1e-4, 4);
	TEST_MM(c, comp, conjc, absc, 1e-12, 9);
	TEST_MM(z, fcomplex, conjzf, abszf, 1e-4, 9);
}
/*Assert A*B is identity.*/
static void assert_inverse(const dmat *A, const dmat *B, real tol){
	dmat *C=NULL;
	dmm(&C, 0, A, B, "nn", 1);
	for(long j=0; j<NY(C); j++){
		for(long i=0; i<NX(C); i++){
			assert_float_equal(P(C, i, j), i==j?1:0, tol);
		}
	}
	dfree(C);
}
/** inv, invspd, inv_inplace, invspd_inplace and chol for several types. */
static void blas_inverse(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	dmat *R=dnew(6, 6);
	drandn(R, 1, &rstat);
	dmat *S=NULL;
	dmm(&S, 0, R, R, "tn", 1);
	daddI(S, 0.5);
	dmat *Si=dinvspd(S);
	assert_inverse(S, Si, 1e-10);
	dmat *Ri=dinv(R);
	assert_inverse(R, Ri, 1e-9);
	assert_null(dinv(NULL));
	assert_null(dinvspd(NULL));
	dinv_inplace(NULL);
	dinvspd_inplace(NULL);
	{//empty matrices are no-op
		dmat *E=dnew(0, 0);
		dinvspd_inplace(E);
		dinv_inplace(E);
		assert_int_equal(PN(E), 0);
		dfree(E);
	}
	{//failure raises an error without writing any file and without touching non-square input
		dmat *E=dnew(3, 3);//singular / not positive definite
		dinvspd_inplace(E);
		dinv_inplace(E);
		dmat *N=dnew(3, 2);
		P(N, 0)=1;
		dinvspd_inplace(N);
		dinv_inplace(N);
		assert_float_equal(P(N, 0), 1, 0);
		assert_null(dchol(N));
		dfree(N);
		dfree(E);
	}
	//cholesky: S=L*L', L lower triangular
	dmat *L=dchol(S);
	for(long j=0; j<6; j++) for(long i=0; i<j; i++) assert_float_equal(P(L, i, j), 0, 0);
	dmat *LL=NULL;
	dmm(&LL, 0, L, L, "nt", 1);
	assert_dmat_near(LL, S, 1e-10);
	dmat *one=dnew(1, 1);
	P(one, 0)=9;
	dmat *L1=dchol(one);
	assert_float_equal(P(L1, 0), 3, 0);
	assert_null(dchol(NULL));
	dfree(L1); dfree(one); dfree(LL); dfree(L);
	//complex
	cmat *C=cnew(4, 4);
	crandn(C, 1, &rstat);
	cmat *Ci=cinv(C);
	cmat *CC=NULL;
	cmm(&CC, 0, C, Ci, "nn", 1);
	for(long j=0; j<4; j++) for(long i=0; i<4; i++) assert_float_equal(cabs(P(CC, i, j)-(i==j?1.:0.)), 0, 1e-10);
	//hermitian positive definite
	cmat *H=NULL;
	cmm(&H, 0, C, C, "cn", 1);
	caddI(H, 1);
	cmat *Hi=cinvspd(H);
	cmm(&CC, 0, H, Hi, "nn", 1);
	for(long j=0; j<4; j++) for(long i=0; i<4; i++) assert_float_equal(cabs(P(CC, i, j)-(i==j?1.:0.)), 0, 1e-10);
	cmat *Hl=cchol(H);
	cmm(&CC, 0, Hl, Hl, "nc", 1);
	for(long i=0; i<16; i++) assert_float_equal(cabs(P(CC, i)-P(H, i)), 0, 1e-10);
	cfree(Hl); cfree(Hi); cfree(H); cfree(CC); cfree(Ci); cfree(C);
	//single precision
	smat *Ss=snew(3, 3);
	P(Ss, 0, 0)=4; P(Ss, 1, 1)=2; P(Ss, 2, 2)=1; P(Ss, 1, 0)=P(Ss, 0, 1)=1;
	smat *Ssi=sinvspd(Ss);
	smat *SS=NULL;
	smm(&SS, 0, Ss, Ssi, "nn", 1);
	for(long j=0; j<3; j++) for(long i=0; i<3; i++) assert_float_equal(P(SS, i, j), i==j?1:0, 1e-5);
	sinv_inplace(Ss);
	for(long i=0; i<9; i++) assert_float_equal(P(Ss, i), P(Ssi, i), 1e-5);
	sfree(SS); sfree(Ssi); sfree(Ss);
	dfree(Ri); dfree(Si); dfree(S); dfree(R);
}
/*U*diag(S)*VT.*/
#define TEST_SVD(X, XR, T, ABSF, TOL, M, N)				\
	{								\
		rand_t rstat;						\
		seed_rand(&rstat, 7);					\
		X##mat *A=X##new(M, N);					\
		X##randn(A, 1, &rstat);					\
		X##mat *U=NULL, *VT=NULL; XR##mat *Sd=NULL;		\
		X##svd(&U, &Sd, &VT, A);				\
		long ns=M<N?M:N;					\
		assert_int_equal(NX(Sd), ns);				\
		assert_int_equal(NY(U), ns);				\
		assert_int_equal(NX(VT), ns);				\
		for(long i=1; i<ns; i++) assert_true(P(Sd, i)<=P(Sd, i-1)); \
		X##mat *US=X##dup(U);					\
		X##muldiag(US, Sd);					\
		X##mat *A2=NULL;					\
		X##mm(&A2, (T)0, US, VT, "nn", (T)1);			\
		for(long i=0; i<PN(A); i++) assert_float_equal(ABSF(P(A2, i)-P(A, i)), 0, TOL); \
		X##mat *U2=NULL, *VT2=NULL; XR##mat *Sd2=NULL;		\
		X##svd_cache(&U2, &Sd2, &VT2, A);			\
		for(long i=0; i<ns; i++) assert_float_equal(P(Sd2, i), P(Sd, i), TOL); \
		X##free(U2); X##free(VT2); XR##free(Sd2);		\
		X##svd(NULL, &Sd2, NULL, A);/*only singular values*/	\
		for(long i=0; i<ns; i++) assert_float_equal(P(Sd2, i), P(Sd, i), TOL); \
		XR##free(Sd2);						\
		X##free(US); X##free(A2);				\
		X##free(U); X##free(VT); XR##free(Sd); X##free(A);	\
	}
/** svd for d, s, c, z types and tall/wide matrices reconstructs the input. */
static void blas_svd(void **state){
	(void)state;
	TEST_SVD(d, d, real, fabs, 1e-10, 6, 4);
	TEST_SVD(d, d, real, fabs, 1e-10, 3, 5);
	TEST_SVD(s, s, float, abssf, 1e-4, 5, 5);
	TEST_SVD(c, d, comp, absc, 1e-10, 4, 6);
	TEST_SVD(z, s, fcomplex, abszf, 1e-4, 5, 3);
	//empty input
	dmat *U=dnew(1, 1), *VT=dnew(1, 1), *Sd=dnew(1, 1);
	dmat *U0=U, *VT0=VT, *Sd0=Sd;
	dmat *E=dnew(0, 0);
	dsvd(&U, &Sd, &VT, E);
	assert_null(U);
	assert_null(VT);
	assert_null(Sd);
	dfree(U0); dfree(VT0); dfree(Sd0); dfree(E);
}
/** svd_pow2 drops small singular values according to both thresholds. */
static void blas_svd_pow(void **state){
	(void)state;
	//A=Q*diag(s)*Q' with known singular values
	rand_t rstat;
	seed_rand(&rstat, 9);
	const int n=4;
	dmat *Q=dnew(n, n);
	drandn(Q, 1, &rstat);
	dgramschmidt(Q, NULL);
	dscale(Q, 1./sqrt(n));//orthonormal columns
	const real s[4]={1, 0.1, 1e-3, 1e-9};
	dmat *A=NULL;
	dmat *QS=ddup(Q);
	dmat *sv=dnew(n, 1);
	for(int i=0; i<n; i++) P(sv, i)=s[i];
	dmuldiag(QS, sv);
	dmm(&A, 0, QS, Q, "nt", 1);
	//reference with given power for the kept modes
	for(int nkeep=1; nkeep<=4; nkeep++){
		real t1=0, t2=0;
		switch(nkeep){
		case 1: t1=0; t2=0.5; break;//relative to previous
		case 2: t1=1e-2; t2=0; break;
		case 3: t1=1e-7; t2=1e-3; break;
		case 4: t1=0; t2=0; break;
		}
		dmat *B=ddup(A);
		dsvd_pow2(B, -0.5, t1, t2);
		dmat *ref=NULL;
		for(int i=0; i<n; i++) P(sv, i)=i<nkeep?pow(s[i], -0.5):0;
		dcp(&QS, Q);
		dmuldiag(QS, sv);
		dmm(&ref, 0, QS, Q, "nt", 1);
		assert_dmat_near(B, ref, 1e-6*pow(s[nkeep-1], -0.5));
		dfree(ref);
		dfree(B);
	}
	dmat *B=ddup(A);
	dsvd_pow(B, 1);//default thresholds keep the first 3 modes
	for(int i=0; i<n; i++) P(sv, i)=i<3?s[i]:0;
	dcp(&QS, Q);
	dmuldiag(QS, sv);
	dmat *ref=NULL;
	dmm(&ref, 0, QS, Q, "nt", 1);
	assert_dmat_near(B, ref, 1e-10);
	dsvd_pow2(NULL, 1, 0, 0);
	//complex
	cmat *C=cnew(3, 3);
	crandn(C, 1, &rstat);
	cmat *Ci=cdup(C);
	csvd_pow(Ci, -1);
	cmat *CC=NULL;
	cmm(&CC, 0, C, Ci, "nn", 1);
	for(long j=0; j<3; j++) for(long i=0; i<3; i++) assert_float_equal(cabs(P(CC, i, j)-(i==j?1.:0.)), 0, 1e-10);
	cfree(CC); cfree(Ci); cfree(C);
	dfree(ref); dfree(B);
	dfree(sv); dfree(QS); dfree(A); dfree(Q);
}
/** expm of nilpotent, diagonal and rotation generators. */
static void blas_expm(void **state){
	(void)state;
	dmat *N=dnew(2, 2);
	P(N, 0, 1)=1;
	dmat *E=NULL;
	dexpm(&E, 0, N, 3);//I+3N
	assert_float_equal(P(E, 0, 0), 1, 1e-12);
	assert_float_equal(P(E, 0, 1), 3, 1e-12);
	assert_float_equal(P(E, 1, 0), 0, 1e-12);
	assert_float_equal(P(E, 1, 1), 1, 1e-12);
	dmat *D=dnew(3, 3);
	P(D, 0, 0)=1; P(D, 1, 1)=-2; P(D, 2, 2)=5;//large: needs scaling and squaring
	dmat *ED=NULL;
	dexpm(&ED, 0, D, 1);
	assert_float_equal(P(ED, 0, 0), exp(1), 1e-10);
	assert_float_equal(P(ED, 1, 1), exp(-2), 1e-10);
	assert_float_equal(P(ED, 2, 2), exp(5), 1e-8*exp(5));
	assert_float_equal(P(ED, 0, 1), 0, 1e-12);
	//rotation and accumulation
	dmat *G=dnew(2, 2);
	P(G, 1, 0)=1; P(G, 0, 1)=-1;
	const real th=2.5;
	dexpm(&E, 2, G, th);//E=2*E+R(th)
	assert_float_equal(P(E, 0, 0), 2+cos(th), 1e-10);
	assert_float_equal(P(E, 1, 0), sin(th), 1e-10);
	assert_float_equal(P(E, 0, 1), 6-sin(th), 1e-10);
	assert_float_equal(P(E, 1, 1), 2+cos(th), 1e-10);
	dfree(G); dfree(ED); dfree(D); dfree(E); dfree(N);
}
/** pinv2 with dense matrix weighting and thresholds. */
static void blas_pinv2(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 3);
	dmat *A=dnew(8, 3);
	drandn(A, 1, &rstat);
	dmat *W=dnew(8, 8);
	for(int i=0; i<8; i++) P(W, i, i)=1+i;
	P(W, 1, 0)=P(W, 0, 1)=0.5;
	dmat *Ap=dpinv(A, W);
	//Ap*A=I
	assert_inverse(Ap, A, 1e-10);
	//Ap=inv(A'WA)A'W
	dmat *AtW=NULL, *AtWA=NULL, *ref=NULL;
	dmm(&AtW, 0, A, W, "tn", 1);
	dmm(&AtWA, 0, AtW, A, "nn", 1);
	dinv_inplace(AtWA);
	dmm(&ref, 0, AtWA, AtW, "nn", 1);
	assert_dmat_near(Ap, ref, 1e-10);
	assert_null(dpinv(NULL, W));
	//rank deficient: duplicated column is handled by the threshold
	dmat *A2=dnew(8, 3);
	for(int i=0; i<8; i++){ P(A2, i, 0)=P(A, i, 0); P(A2, i, 1)=P(A, i, 1); P(A2, i, 2)=P(A, i, 0); }
	dmat *A2p=dpinv2(A2, NULL, 1e-7, 0);
	dmat *P2=NULL;
	dmm(&P2, 0, A2, A2p, "nn", 1);//projection onto the range of A2
	dmat *P2A=NULL;
	dmm(&P2A, 0, P2, A2, "nn", 1);
	assert_dmat_near(P2A, A2, 1e-8);
	assert_false(disnan(A2p));
	dfree(P2A); dfree(P2); dfree(A2p); dfree(A2);
	dfree(AtW); dfree(AtWA); dfree(ref); dfree(Ap); dfree(W); dfree(A);
}

UNITTEST_MAIN(
	cmocka_unit_test(blas_dpinv),
	cmocka_unit_test(blas_mm),
	cmocka_unit_test(blas_inverse),
	cmocka_unit_test(blas_svd),
	cmocka_unit_test(blas_svd_pow),
	cmocka_unit_test(blas_expm),
	cmocka_unit_test(blas_pinv2),
)
