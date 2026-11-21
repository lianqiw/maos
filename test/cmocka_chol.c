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
  Tests of the cholmod based sparse Cholesky factorization (port of
  test_chol.c). The test matrix is a 2-D Laplacian plus identity, which is
  symmetric positive definite.
*/
#if CPU_SINGLE
#define RTOL 1e-4
#else
#define RTOL 1e-10
#endif
#define NGRID 20 //grid size; the matrix is NGRID^2 square

/** 2-D 5-point Laplacian on an n*n grid plus alpha*I, as a sparse matrix. */
static dsp *mklaplacian(long n, real alpha){
	long N=n*n;
	dsp *A=dspnew(N, N, 5*N);
	long count=0;
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			long ic=ix+iy*n;
			A->pp[ic]=count;
			//rows in increasing order
			if(iy>0){ A->pi[count]=ic-n; A->px[count]=-1; count++; }
			if(ix>0){ A->pi[count]=ic-1; A->px[count]=-1; count++; }
			A->pi[count]=ic; A->px[count]=4+alpha; count++;
			if(ix+1<n){ A->pi[count]=ic+1; A->px[count]=-1; count++; }
			if(iy+1<n){ A->pi[count]=ic+n; A->px[count]=-1; count++; }
		}
	}
	A->pp[N]=count;
	dspsetnzmax(A, count);
	return A;
}
/** Relative residual |A*x-y|/|y| for each column, maximum over columns. */
static real residual(const dsp *A, const dmat *x, const dmat *y){
	dmat *r=ddup(y);
	dspmm(&r, A, x, "nn", -1);
	real maxres=0;
	for(long ic=0; ic<NY(y); ic++){
		real rn=0, yn=0;
		for(long ix=0; ix<NX(y); ix++){
			rn+=P(r, ix, ic)*P(r, ix, ic);
			yn+=P(y, ix, ic)*P(y, ix, ic);
		}
		real res=sqrt(rn/yn);
		if(res>maxres) maxres=res;
	}
	dfree(r);
	return maxres;
}
/** Check that L*L'=A(p,p), i.e. (L*L')(i,j)=A(p[i],p[j]). */
static void check_factor(const dsp *L, const long *p, const dsp *A){
	long N=NX(A);
	dmat *FL=NULL, *FA=NULL, *LLt=NULL;
	dspfull(&FL, L, 'n', 1);
	dspfull(&FA, A, 'n', 1);
	dmm(&LLt, 0, FL, FL, "nt", 1);
	real maxd=0;
	for(long j=0; j<N; j++){
		for(long i=0; i<N; i++){
			real d=fabs(P(LLt, i, j)-P(FA, p[i], p[j]));
			if(d>maxd) maxd=d;
		}
	}
	assert_float_equal(maxd, 0, RTOL*10);
	dfree(FL); dfree(FA); dfree(LLt);
}

/** chol_factorize + chol_solve with the native cholmod factor: small residual, multi-column solve matches column solves, reuse of x. */
static void chol_cholmod(void **state){
	(void)state;
	dsp *A=mklaplacian(NGRID, 1);
	long N=NX(A);
	spchol *R=chol_factorize(A);
	assert_non_null(R);
	assert_non_null(R->L);
	assert_null(R->Cl);
	assert_null(R->Cu);
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *y=dnew(N, 1);
	drandn(y, 1, &rstat);
	dmat *x=NULL;
	chol_solve(&x, R, y);
	assert_int_equal(NX(x), N);
	assert_int_equal(NY(x), 1);
	assert_true(residual(A, x, y)<RTOL);
	//Solving again into the existing x gives the same answer.
	dmat *x0=ddup(x);
	chol_solve(&x, R, y);
	assert_dmat_near(x, x0, 0);
	//A known solution is recovered.
	dmat *xt=dnew(N, 1);
	drandn(xt, 1, &rstat);
	dmat *yt=dnew(N, 1);
	dspmv(yt, A, xt, 'n', 1);
	dmat *xs=NULL;
	chol_solve(&xs, R, yt);
	assert_dmat_near(xs, xt, RTOL*10);
	//Multiple right hand sides (threaded path) match single column solves.
	const long nrhs=5;
	dmat *Y=dnew(N, nrhs);
	drandn(Y, 1, &rstat);
	dmat *X=NULL;
	chol_solve(&X, R, Y);
	assert_int_equal(NY(X), nrhs);
	assert_true(residual(A, X, Y)<RTOL);
	for(long ic=0; ic<nrhs; ic++){
		dmat *yc=dsub(Y, 0, N, ic, 1);
		dmat *xc=NULL;
		chol_solve(&xc, R, yc);
		dmat *Xc=dsub(X, 0, N, ic, 1);
		assert_dmat_near(xc, Xc, RTOL);
		dfree(xc); dfree(yc); dfree(Xc);
	}
	dfree(X); dfree(Y);
	dfree(xs); dfree(xt); dfree(yt);
	dfree(x0); dfree(x); dfree(y);
	chol_free(R);
	assert_null(R);
	dspfree(A);
}

/** Only the upper triangular part of A is used: factorizing triu(A) gives the same solution. */
static void chol_upper_input(void **state){
	(void)state;
	dsp *A=mklaplacian(NGRID, 0.5);
	long N=NX(A);
	dmat *FA=NULL;
	dspfull(&FA, A, 'n', 1);
	for(long iy=0; iy<N; iy++){
		for(long ix=iy+1; ix<N; ix++){
			P(FA, ix, iy)=0;
		}
	}
	dsp *Au=d2sp(FA, 0);
	assert_int_equal(dspcheck(Au), 2);
	spchol *R1=chol_factorize(A);
	spchol *R2=chol_factorize(Au);
	rand_t rstat;
	seed_rand(&rstat, 2);
	dmat *y=dnew(N, 1);
	drandn(y, 1, &rstat);
	dmat *x1=NULL, *x2=NULL;
	chol_solve(&x1, R1, y);
	chol_solve(&x2, R2, y);
	assert_dmat_near(x1, x2, RTOL);
	assert_true(residual(A, x2, y)<RTOL);
	dfree(x1); dfree(x2); dfree(y);
	chol_free(R1); chol_free(R2);
	dspfree(Au); dfree(FA);
	dspfree(A);
}

/**
   chol_convert to the lower factor Cl (keep=0) and then to the upper factor
   Cu=Cl': the factors reproduce A(p,p) and all three solve paths agree (old
   test_chol main).
*/
static void chol_convert_paths(void **state){
	(void)state;
	dsp *A=mklaplacian(NGRID, 1);
	long N=NX(A);
	spchol *R=chol_factorize(A);
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *y=dnew(N, 3);
	drandn(y, 1, &rstat);
	dmat *x=NULL, *x2=NULL, *x3=NULL;
	chol_solve(&x, R, y);//cholmod

	chol_convert(R, 0);//convert to Cl and free the cholmod factor
	assert_null(R->L);
	assert_null(R->c);
	assert_non_null(R->Cl);
	assert_non_null(R->Cp);
	assert_int_equal(NX(R->Cl), N);
	assert_int_equal(NY(R->Cl), N);
	assert_int_equal(dspcheck(R->Cl), 1);//lower triangular
	//Cp is a permutation
	long *perm=mycalloc(N, long);
	int *seen=mycalloc(N, int);
	for(long i=0; i<N; i++){
		perm[i]=R->Cp[i];
		assert_true(perm[i]>=0&&perm[i]<N);
		assert_false(seen[perm[i]]);
		seen[perm[i]]=1;
	}
	free(seen);
	check_factor(R->Cl, perm, A);
	//converting again is a no-op
	dsp *Cl=R->Cl;
	chol_convert(R, 0);
	assert_ptr_equal(R->Cl, Cl);
	chol_solve(&x2, R, y);//lower
	assert_dmat_near(x2, x, RTOL);
	//single column through the lower path
	dmat *y1=dsub(y, 0, N, 1, 1);
	dmat *x21=NULL;
	chol_solve(&x21, R, y1);
	dmat *x1=dsub(x, 0, N, 1, 1);
	assert_dmat_near(x21, x1, RTOL);
	dfree(x21);

	R->Cu=dsptrans(R->Cl);
	dspfree(R->Cl);
	assert_int_equal(dspcheck(R->Cu), 2);//upper triangular
	chol_solve(&x3, R, y);//upper
	assert_dmat_near(x3, x, RTOL);
	assert_true(residual(A, x3, y)<RTOL);
	dmat *x31=NULL;
	chol_solve(&x31, R, y1);
	assert_dmat_near(x31, x1, RTOL);
	dfree(x31); dfree(x1); dfree(y1);
	free(perm);
	dfree(x); dfree(x2); dfree(x3); dfree(y);
	chol_free(R);
	dspfree(A);
}

/** chol_convert with keep=1 keeps the cholmod factor (used for solving) and chol_factorize2 returns Cl and the permutation. */
static void chol_keep(void **state){
	(void)state;
	dsp *A=mklaplacian(NGRID, 2);
	long N=NX(A);
	spchol *R=chol_factorize(A);
	rand_t rstat;
	seed_rand(&rstat, 3);
	dmat *y=dnew(N, 1);
	drandn(y, 1, &rstat);
	dmat *x=NULL, *x2=NULL;
	chol_solve(&x, R, y);
	chol_convert(R, 1);
	assert_non_null(R->L);
	assert_non_null(R->Cl);
	assert_non_null(R->Cp);
	chol_solve(&x2, R, y);
	assert_dmat_near(x2, x, RTOL);

	lmat *Cp=NULL;
	dsp *Cl=chol_factorize2(&Cp, A);
	assert_non_null(Cl);
	assert_int_equal(PN(Cp), N);
	assert_int_equal(dspcheck(Cl), 1);
	//Same ordering and factor as the converted one.
	for(long i=0; i<N; i++){
		assert_int_equal(P(Cp, i), R->Cp[i]);
	}
	dmat *F1=NULL, *F2=NULL;
	dspfull(&F1, Cl, 'n', 1);
	dspfull(&F2, R->Cl, 'n', 1);
	assert_dmat_near(F1, F2, RTOL);
	check_factor(Cl, P(Cp), A);
	//Without Cp output
	dsp *Cl2=chol_factorize2(NULL, A);
	dmat *F3=NULL;
	dspfull(&F3, Cl2, 'n', 1);
	assert_dmat_near(F3, F1, 0);
	dfree(F1); dfree(F2); dfree(F3);
	dspfree(Cl); dspfree(Cl2); lfree(Cp);
	dfree(x); dfree(x2); dfree(y);
	chol_free(R);
	dspfree(A);
}

/**
   Solving for columns of the identity gives columns of A^-1, which satisfy
   A*X=I and the symmetry of A^-1. chol_spsolve() with a sparse right hand
   side gives the same result.
*/
static void chol_inverse(void **state){
	(void)state;
	dsp *A=mklaplacian(NGRID, 1);
	long N=NX(A);
	spchol *R=chol_factorize(A);
	const long nc=4;
	const long cols[4]={0, 7, N/2, N-1};
	dmat *Yd=dnew(N, nc);
	for(long ic=0; ic<nc; ic++){
		P(Yd, cols[ic], ic)=1;
	}
	dmat *Xd=NULL;
	chol_solve(&Xd, R, Yd);
	assert_true(residual(A, Xd, Yd)<RTOL);
	//A^-1 is symmetric: X(cols[i], j)==X(cols[j], i)
	for(long i=0; i<nc; i++){
		for(long j=0; j<nc; j++){
			assert_float_equal(P(Xd, cols[i], j), P(Xd, cols[j], i), RTOL);
		}
		//A is a diagonally dominant M-matrix, so A^-1 is positive with a dominant diagonal.
		for(long ix=0; ix<N; ix++){
			assert_true(P(Xd, ix, i)>0);
			assert_true(P(Xd, ix, i)<=P(Xd, cols[i], i));
		}
	}
	//sparse right hand side
	dsp *Ys=d2sp(Yd, 0);
	dsp *Xs=chol_spsolve(R, Ys);
	assert_non_null(Xs);
	dmat *Xs2=NULL;
	dspfull(&Xs2, Xs, 'n', 1);
	assert_dmat_near(Xs2, Xd, RTOL);
	dfree(Xs2);
	dspfree(Xs);
	dspfree(Ys);
	dfree(Xd); dfree(Yd);
	chol_free(R);
	dspfree(A);
}

/** NULL input: chol_factorize and chol_factorize2 return NULL; chol_convert and chol_free accept NULL. */
static void chol_null(void **state){
	(void)state;
	spchol *R=chol_factorize(NULL);
	assert_null(R);
	chol_convert(NULL, 0);
	chol_free(R);
	assert_null(R);
	assert_null(chol_factorize2(NULL, NULL));
}

UNITTEST_MAIN(
	cmocka_unit_test(chol_cholmod),
	cmocka_unit_test(chol_upper_input),
	cmocka_unit_test(chol_convert_paths),
	cmocka_unit_test(chol_keep),
	cmocka_unit_test(chol_inverse),
	cmocka_unit_test(chol_null),
)
