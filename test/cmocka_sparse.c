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
  Tests of sparse matrix routines (port of test_sparse.c). Each sparse
  operation is compared against the same operation on the dense equivalent
  obtained with dspfull.
*/
#define TOL (100*EPS)

/** Dense copy of sparse A. */
static dmat *full(const dsp *A){
	dmat *F=NULL;
	dspfull(&F, A, 'n', 1);
	return F;
}
/** Random sparse matrix with values in [-0.5, 0.5). */
static dsp *randsp(long nx, long ny, real fill, rand_t *rstat){
	dsp *A=dspnewrandu(nx, ny, 1, fill, rstat);
	for(long i=0; i<A->pp[A->ny]; i++){
		A->px[i]-=0.5;
	}
	return A;
}
/** Assert sparse A equals dense B. */
static void assert_sp_dense(const dsp *A, const dmat *B, real tol){
	dmat *F=full(A);
	assert_dmat_near(F, B, tol);
	dfree(F);
}

/** Whether row indices are strictly increasing in each column. */
static int sorted(const dsp *A){
	for(long ic=0; ic<NY(A); ic++){
		for(long i=A->pp[ic]+1; i<A->pp[ic+1]; i++){
			if(A->pi[i]<=A->pi[i-1]) return 0;
		}
	}
	return 1;
}
/** Construction: random fill factor, diagonal create/extract, d2sp thresholding, spcheck triangle detection, dup. */
static void sparse_create(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long nx=200, ny=100;
	const real fill=0.05;
	dsp *A=dspnewrandu(nx, ny, 2, fill, &rstat);
	assert_int_equal(NX(A), nx);
	assert_int_equal(NY(A), ny);
	long nz=A->pp[ny];
	//number of nonzeros is binomial(nx*ny, fill)
	assert_float_equal(nz, nx*ny*fill, 6*sqrt(nx*ny*fill*(1-fill)));
	for(long i=0; i<nz; i++){
		assert_true(A->px[i]>=0&&A->px[i]<2);
		assert_true(A->pi[i]>=0&&A->pi[i]<nx);
	}
	assert_int_equal(dspcheck(A), 0);//general matrix
	assert_true(sorted(A));
	//dup is an exact copy
	dsp *B=dspdup(A);
	dmat *FA=full(A);
	assert_sp_dense(B, FA, 0);
	//d2sp of the full matrix recovers A; the threshold drops small values.
	dsp *C=d2sp(FA, 0);
	assert_int_equal(C->pp[ny], nz);
	assert_sp_dense(C, FA, 0);
	dspfree(C);
	C=d2sp(FA, 1);
	long nbig=0;
	for(long i=0; i<PN(FA); i++){
		if(fabs(P(FA, i))>1) nbig++; else P(FA, i)=0;
	}
	assert_int_equal(C->pp[ny], nbig);
	assert_sp_dense(C, FA, 0);
	dspfree(C);
	dfree(FA);
	dspfree(A);
	dspfree(B);

	//diagonal matrix
	const long n=10;
	real w[10];
	for(long i=0; i<n; i++) w[i]=i+1;
	A=dspnewdiag(n, w, 2);
	assert_int_equal(dspcheck(A), 3);
	dmat *d=dspdiag(A);
	for(long i=0; i<n; i++){
		assert_float_equal(P(d, i), 2*(i+1), 0);
	}
	assert_float_equal(dsptrace(A, 1), 2*n*(n+1)/2, TOL);
	assert_float_equal(dsptrace(A, 2), 4*n*(n+1)*(2*n+1)/6, TOL*1000);
	dfree(d);
	dspfree(A);
	A=dspnewdiag(n, NULL, 3);
	d=dspdiag(A);
	for(long i=0; i<n; i++){
		assert_float_equal(P(d, i), 3, 0);
	}
	dfree(d);
	dspfree(A);

	//triangular matrices
	dmat *L=dnew(n, n);
	for(long iy=0; iy<n; iy++){
		for(long ix=iy; ix<n; ix++){
			P(L, ix, iy)=ix+iy+1;
		}
	}
	A=d2sp(L, 0);
	assert_int_equal(dspcheck(A), 1);
	B=dsptrans(A);
	assert_int_equal(dspcheck(B), 2);
	dspfree(A);
	dspfree(B);
	dfree(L);

	//NULL input
	assert_null(dspdup(NULL));
	assert_null(dsptrans(NULL));
	assert_null(dspdiag(NULL));
	assert_null(dspaddI(NULL, 1));
	assert_null(dspsum(NULL, 1));
	assert_null(dspsumabs(NULL, 1));
	assert_null(d2sp(NULL, 0));
	assert_null(dspcat(NULL, NULL, 1));
	assert_float_equal(dsptrace(NULL, 1), 0, 0);
	dmat *F=NULL;
	dspfull(&F, NULL, 'n', 1);
	assert_null(F);
	assert_int_equal(dspcheck(NULL), 3);
}

/** dspmv: y+=alpha*op(A)*x for 'n' and 't' against dense, accumulation, and the adjoint identity (old test_spmul). */
static void sparse_spmv(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long nx=200, ny=100;
	dsp *A0=randsp(nx, ny, 0.05, &rstat);
	dsp *A=dspaddI(A0, 1);
	dmat *F=full(A);
	dmat *F0=full(A0);
	//spaddI adds 1 to the first min(nx,ny) diagonal elements only.
	for(long i=0; i<ny; i++) P(F0, i, i)+=1;
	assert_dmat_near(F, F0, TOL);
	dfree(F0);
	dmat *x=dnew(ny, 1);
	drandn(x, 1, &rstat);
	dmat *y=dnew(nx, 1);
	dmat *yref=NULL;
	dmm(&yref, 0, F, x, "nn", 1);
	dspmv(y, A, x, 'n', 1);
	assert_dmat_near(y, yref, TOL);
	//accumulates: y now has 1+2=3 times A*x
	dspmv(y, A, x, 'n', 2);
	dscale(yref, 3);
	assert_dmat_near(y, yref, TOL);
	//transpose
	drandn(y, 1, &rstat);
	dmat *x2=dnew(ny, 1);
	dmat *x2ref=NULL;
	dmm(&x2ref, 0, F, y, "tn", 1);
	dspmv(x2, A, y, 't', 1);
	assert_dmat_near(x2, x2ref, TOL);
	//adjoint identity: y'*(A*x)=(A'*y)'*x
	dmat *Ax=dnew(nx, 1);
	dspmv(Ax, A, x, 'n', 1);
	assert_float_equal(ddot(y, Ax), ddot(x2, x), TOL*nx);
	assert_float_equal(dspwdot(y, A, x), ddot(y, Ax), TOL*nx);
	//spwdot with NULL A is a plain dot product
	assert_float_equal(dspwdot(x, NULL, x2), ddot(x, x2), TOL);
	assert_float_equal(dspwdot(NULL, A, x), 0, 0);
	//NULL A does nothing
	dmat *Ax2=ddup(Ax);
	dspmv(Ax, NULL, x, 'n', 1);
	assert_dmat_near(Ax, Ax2, 0);
	dfree(Ax2);
	dfree(Ax);
	dfree(x); dfree(y); dfree(yref); dfree(x2); dfree(x2ref);
	dfree(F);
	dspfree(A);
	dspfree(A0);
}

/** dspmm and dmulsp for all transpose combinations against dmm on the dense matrices. */
static void sparse_spmm(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	const long n1=30, n2=20, nc=7;
	dsp *A=randsp(n1, n2, 0.2, &rstat);
	dmat *F=full(A);
	const char *trans[]={"nn", "tn", "nt", "tt"};
	for(int it=0; it<4; it++){
		const char *tr=trans[it];
		//op(A)*op(x)
		long nax=tr[0]=='n'?n2:n1;//inner dimension
		dmat *x=tr[1]=='n'?dnew(nax, nc):dnew(nc, nax);
		drandn(x, 1, &rstat);
		dmat *y=NULL, *yref=NULL;
		dspmm(&y, A, x, tr, 2);
		dmm(&yref, 0, F, x, tr, 2);
		assert_dmat_near(y, yref, TOL);
		//accumulates
		dspmm(&y, A, x, tr, -1);
		dmm(&yref, 1, F, x, tr, -1);
		assert_dmat_near(y, yref, TOL);
		dfree(y); dfree(yref);
		//op(x)*op(A): x now on the left
		long nbx=tr[1]=='n'?n1:n2;//inner dimension
		dmat *z=tr[0]=='n'?dnew(nc, nbx):dnew(nbx, nc);
		drandn(z, 1, &rstat);
		dmulsp(&y, z, A, tr, 3);
		dmm(&yref, 0, z, F, tr, 3);
		assert_dmat_near(y, yref, TOL);
		dfree(y); dfree(yref);
		dfree(x); dfree(z);
	}
	//dcellmm dispatches sparse*dense, dense*sparse, and sparse*sparse into dense.
	dmat *x=dnew(n2, nc);
	drandn(x, 1, &rstat);
	dmat *y=NULL, *yref=NULL;
	dcellmm(&y, A, x, "nn", 1);
	dmm(&yref, 0, F, x, "nn", 1);
	assert_dmat_near(y, yref, TOL);
	dfree(y); dfree(yref);
	dmat *AtA=dnew(n2, n2);
	dmat *AtAref=NULL;
	dcellmm(&AtA, A, A, "tn", 1);
	dmm(&AtAref, 0, F, F, "tn", 1);
	assert_dmat_near(AtA, AtAref, TOL);
	dfree(AtA); dfree(AtAref);
	dfree(x);
	dfree(F);
	dspfree(A);
}

/** dspmulsp/dspmulsp2 for all transpose combinations; multiplying by a diagonal matrix scales columns (old disabled test_spmul). */
static void sparse_spmulsp(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long n=10;
	//A*W with W diagonal equals A with column i scaled by w[i].
	dsp *A=dspnewrandu(n, n, 10, 0.2, &rstat);
	real w[10];
	for(long i=0; i<n; i++) w[i]=2+i;
	dsp *W=dspnewdiag(n, w, 1);
	dsp *AW=dspmulsp(A, W, "nn");
	dsp *AW2=dspdup(A);
	dspmuldiag(AW2, w, 1);
	dsp *diff=dspadd2(AW, 1, AW2, -1);
	dmat *Fd=full(diff);
	assert_float_equal(dmaxabs(Fd), 0, TOL*10);
	dfree(Fd);
	dspfree(diff);
	//same with spscaley
	dsp *AW3=dspdup(A);
	dmat *wv=dnew(n, 1);
	for(long i=0; i<n; i++) P(wv, i)=w[i];
	dspscaley(AW3, wv);
	dmat *FAW=full(AW);
	assert_sp_dense(AW3, FAW, TOL*10);
	//W*A scales rows: spscalex
	dsp *WA=dspmulsp(W, A, "nn");
	dsp *WA2=dspdup(A);
	dspscalex(WA2, wv);
	dmat *FWA=full(WA);
	assert_sp_dense(WA2, FWA, TOL*10);
	dfree(FWA); dfree(FAW); dfree(wv);
	dspfree(WA); dspfree(WA2); dspfree(AW3);
	dspfree(AW); dspfree(AW2); dspfree(W); dspfree(A);

	//General products against dense.
	const long n1=25, n2=15, n3=20;
	const char *trans[]={"nn", "tn", "nt", "tt"};
	for(int it=0; it<4; it++){
		const char *tr=trans[it];
		dsp *B1=tr[0]=='n'?randsp(n1, n2, 0.2, &rstat):randsp(n2, n1, 0.2, &rstat);
		dsp *B2=tr[1]=='n'?randsp(n2, n3, 0.2, &rstat):randsp(n3, n2, 0.2, &rstat);
		dmat *F1=full(B1), *F2=full(B2);
		dmat *Cref=NULL;
		dmm(&Cref, 0, F1, F2, tr, 1);
		dsp *C=dspmulsp(B1, B2, tr);
		assert_int_equal(NX(C), n1);
		assert_int_equal(NY(C), n3);
		assert_sp_dense(C, Cref, TOL);
		//spmulsp2 accumulates with scale
		dspmulsp2(&C, B1, B2, tr, -3);
		dscale(Cref, -2);
		assert_sp_dense(C, Cref, TOL*10);
		//dcellmm with sparse operands and NULL output gives sparse result
		dsp *C2=NULL;
		dcellmm(&C2, B1, B2, tr, 2);
		assert_int_equal(C2->id, M_DSP);
		dscale(Cref, -1);
		assert_sp_dense(C2, Cref, TOL*10);
		dspfree(C2);
		dspfree(C);
		dfree(Cref); dfree(F1); dfree(F2);
		dspfree(B1); dspfree(B2);
	}
	//NULL operands do nothing
	dsp *C=NULL;
	dspmulsp2(&C, NULL, NULL, "nn", 1);
	assert_null(C);
}

/** Transpose: full(A')==full(A)', (A')'==A, and spcelltrans transposes both block layout and blocks. */
static void sparse_trans(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 3);
	dsp *A=randsp(40, 25, 0.1, &rstat);
	dmat *F=full(A);
	dsp *At=dsptrans(A);
	dmat *Ft=dtrans(F);
	assert_sp_dense(At, Ft, 0);
	dmat *Ft2=NULL;
	dspfull(&Ft2, A, 't', 1);
	assert_dmat_near(Ft2, Ft, 0);
	dsp *Att=dsptrans(At);
	assert_sp_dense(Att, F, 0);
	assert_true(sorted(At));//transpose produces sorted rows
	assert_true(sorted(Att));
	dfree(Ft2);

	dspcell *Ac=dspcellnew(2, 3);
	for(long i=0; i<PN(Ac); i++){
		P(Ac, i)=randsp(5+i%2, 4+i/2, 0.3, &rstat);
	}
	dspcell *Act=dspcelltrans(Ac);
	assert_int_equal(NX(Act), 3);
	assert_int_equal(NY(Act), 2);
	for(long iy=0; iy<NY(Ac); iy++){
		for(long ix=0; ix<NX(Ac); ix++){
			dmat *Fi=NULL;
			dspfull(&Fi, P(Ac, ix, iy), 't', 1);
			assert_sp_dense(P(Act, iy, ix), Fi, 0);
			dfree(Fi);
		}
	}
	assert_null(dspcelltrans(NULL));
	dspcellfree(Act);
	dspcellfree(Ac);
	dfree(F); dfree(Ft);
	dspfree(A); dspfree(At); dspfree(Att);
}

/** dspsum/dspsumabs along both dimensions against dense sums (old test_spsum). */
static void sparse_sum(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 4);
	const long nx=50, ny=30;
	dsp *A=randsp(nx, ny, 0.1, &rstat);
	dmat *F=full(A);
	dmat *s1=dspsum(A, 1);
	dmat *s2=dspsum(A, 2);
	dmat *a1=dspsumabs(A, 1);
	dmat *a2=dspsumabs(A, 2);
	assert_int_equal(NX(s1), 1); assert_int_equal(NY(s1), ny);
	assert_int_equal(NX(s2), nx); assert_int_equal(NY(s2), 1);
	assert_int_equal(NX(a1), 1); assert_int_equal(NY(a1), ny);
	assert_int_equal(NX(a2), nx); assert_int_equal(NY(a2), 1);
	for(long iy=0; iy<ny; iy++){
		real s=0, a=0;
		for(long ix=0; ix<nx; ix++){
			s+=P(F, ix, iy); a+=fabs(P(F, ix, iy));
		}
		assert_float_equal(P(s1, iy), s, TOL);
		assert_float_equal(P(a1, iy), a, TOL);
	}
	for(long ix=0; ix<nx; ix++){
		real s=0, a=0;
		for(long iy=0; iy<ny; iy++){
			s+=P(F, ix, iy); a+=fabs(P(F, ix, iy));
		}
		assert_float_equal(P(s2, ix), s, TOL);
		assert_float_equal(P(a2, ix), a, TOL);
	}
	assert_float_equal(dsum(s1), dsum(F), TOL*10);
	assert_float_equal(dsum(s2), dsum(F), TOL*10);
	assert_float_equal(dsum(a1), dsumabs(F), TOL*10);
	assert_float_equal(dsum(a2), dsumabs(F), TOL*10);
	//abs sums dominate plain sums
	for(long i=0; i<ny; i++) assert_true(P(a1, i)>=fabs(P(s1, i))-TOL);
	dfree(s1); dfree(s2); dfree(a1); dfree(a2);
	//trace of a square block
	dsp *S=randsp(20, 20, 0.3, &rstat);
	dmat *FS=full(S);
	real tr=0, tr2=0;
	for(long i=0; i<20; i++){ tr+=P(FS, i, i); tr2+=P(FS, i, i)*P(FS, i, i); }
	assert_float_equal(dsptrace(S, 1), tr, TOL);
	assert_float_equal(dsptrace(S, 2), tr2, TOL);
	dfree(FS);
	dspfree(S);
	dfree(F);
	dspfree(A);
}

/** spadd2/spadd/spaddI/spscale/spdroptol/spsym/spcat/spdropemptycol against dense equivalents. */
static void sparse_add(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 5);
	const long nx=30, ny=20;
	dsp *A=randsp(nx, ny, 0.1, &rstat);
	dsp *B=randsp(nx, ny, 0.1, &rstat);
	dmat *FA=full(A), *FB=full(B);
	dmat *ref=NULL;
	dadd(&ref, 0, FA, 2);
	dadd(&ref, 1, FB, -3);
	dsp *C=dspadd2(A, 2, B, -3);
	assert_sp_dense(C, ref, TOL);
	//spadd2 (CSparse cs_add) may leave row indices unsorted; spsort fixes the order without changing values.
	dspsort(C);
	assert_true(sorted(C));
	assert_sp_dense(C, ref, TOL);
	dspfree(C);
	//spadd into NULL is a scaled copy
	dspadd(&C, 5, B, -3);
	dadd(&ref, 0, FB, -3);
	assert_sp_dense(C, ref, TOL);
	//spadd accumulates A0*alpha+B*beta
	dspadd(&C, 1, A, 2);
	dadd(&ref, 1, FA, 2);
	assert_sp_dense(C, ref, TOL);
	//NULL B does nothing
	dspadd(&C, 1, NULL, 2);
	assert_sp_dense(C, ref, TOL);
	//A-A has no nonzeros left
	dsp *Z=dspadd2(A, 1, A, -1);
	assert_int_equal(Z->pp[NY(Z)], 0);
	dspfree(Z);
	//spscale
	dspscale(C, 0.5);
	dscale(ref, 0.5);
	assert_sp_dense(C, ref, TOL);
	dspfree(C);

	//spaddI on a matrix with missing diagonal elements, rectangular wide and tall
	for(int it=0; it<2; it++){
		dsp *S=it?randsp(8, 12, 0.2, &rstat):randsp(12, 8, 0.2, &rstat);
		dmat *FS=full(S);
		dsp *SI=dspaddI(S, 2);
		for(long i=0; i<MIN(NX(FS), NY(FS)); i++) P(FS, i, i)+=2;
		assert_sp_dense(SI, FS, TOL);
		assert_true(sorted(SI));
		dfree(FS);
		dspfree(SI);
		dspfree(S);
	}

	//spdroptol drops elements not above thres*max|A|
	C=dspdup(A);
	real maxa=dmaxabs(FA);
	dspdroptol(C, 0.5);
	dmat *FAd=ddup(FA);
	for(long i=0; i<PN(FAd); i++){
		if(fabs(P(FAd, i))<=0.5*maxa) P(FAd, i)=0;
	}
	assert_sp_dense(C, FAd, 0);
	dfree(FAd);
	dspfree(C);

	//spsym gives (S+S')/2
	dsp *S=randsp(15, 15, 0.2, &rstat);
	dmat *FS=full(S);
	dmat *FSt=dtrans(FS);
	dadd(&FS, 0.5, FSt, 0.5);
	dspsym(&S);
	assert_sp_dense(S, FS, TOL);
	dsp *St=dsptrans(S);
	assert_sp_dense(St, FS, TOL);
	dspfree(St);
	dspfree(S);
	dfree(FS); dfree(FSt);

	//spcat side by side
	dsp *A2=randsp(nx, 7, 0.2, &rstat);
	dmat *FA2=full(A2);
	C=dspcat(A, A2, 1);
	assert_int_equal(NX(C), nx);
	assert_int_equal(NY(C), ny+7);
	dmat *Fcat=dnew(nx, ny+7);
	for(long iy=0; iy<ny+7; iy++){
		for(long ix=0; ix<nx; ix++){
			P(Fcat, ix, iy)=iy<ny?P(FA, ix, iy):P(FA2, ix, iy-ny);
		}
	}
	assert_sp_dense(C, Fcat, 0);
	dspfree(C);
	//spcat with a NULL operand is a copy
	C=dspcat(A, NULL, 1);
	assert_sp_dense(C, FA, 0);
	dspfree(C);
	dfree(Fcat);
	dfree(FA2);
	dspfree(A2);

	//spdropemptycol removes empty columns and returns the column map
	dmat *D=dnew(5, 6);
	P(D, 1, 0)=1; P(D, 3, 2)=2; P(D, 0, 5)=3; P(D, 4, 5)=4;
	C=d2sp(D, 0);
	lmat *map=dspdropemptycol(C);
	assert_int_equal(NY(C), 3);
	assert_int_equal(PN(map), 3);
	assert_int_equal(P(map, 0), 0);
	assert_int_equal(P(map, 1), 2);
	assert_int_equal(P(map, 2), 5);
	dmat *FC=full(C);
	for(long iy=0; iy<3; iy++){
		for(long ix=0; ix<5; ix++){
			assert_float_equal(P(FC, ix, iy), P(D, ix, P(map, iy)), 0);
		}
	}
	dfree(FC);
	lfree(map);
	dspfree(C);
	dfree(D);

	dfree(ref); dfree(FA); dfree(FB);
	dspfree(A); dspfree(B);
}

/**
   Sparse cell operations (old disabled benchmark): A=G0'*G0 through dcellmm,
   and x'*A*x computed by spcellwdot, by dcellmm+dcelldot, and as |G0*x|^2
   all agree. spcell2sp and spcellfull agree.
*/
static void sparse_cell(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 6);
	const long nr[2]={40, 30};//rows of each block row
	const long nc[2]={20, 25};//cols of each block column
	dspcell *G0=dspcellnew(2, 2);
	for(int iy=0; iy<2; iy++){
		for(int ix=0; ix<2; ix++){
			P(G0, ix, iy)=randsp(nr[ix], nc[iy], 0.1, &rstat);
		}
	}
	dspcell *G0T=dspcelltrans(G0);
	dspcell *A=NULL;
	dcellmm(&A, G0T, G0, "nn", 1);
	assert_int_equal(NX(A), 2);
	assert_int_equal(NY(A), 2);
	dspcell *A2=NULL;
	dcellmm(&A2, G0, G0, "tn", 1);
	dcell *x=dcellnew(2, 1);
	for(int i=0; i<2; i++){
		P(x, i)=dnew(nc[i], 1);
		drandn(P(x, i), 1, &rstat);
	}
	real res=dspcellwdot(x, A, x);
	real res1=dspcellwdot(x, A2, x);
	dcell *y=NULL;
	dcellmm(&y, A, x, "nn", 1);
	real res2=dcelldot(x, y);
	dcell *z=NULL;
	dcellmm(&z, G0, x, "nn", 1);
	real res3=dcelldot(z, z);
	assert_true(res3>0);
	assert_float_equal(res, res3, TOL*res3*10);
	assert_float_equal(res1, res3, TOL*res3*10);
	assert_float_equal(res2, res3, TOL*res3*10);
	//NULL A means identity
	assert_float_equal(dspcellwdot(x, NULL, x), dcelldot(x, x), TOL);
	//G0'*(G0*x) through the transpose flag equals G0T*z
	dcell *u=NULL, *u2=NULL;
	dcellmm(&u, G0, z, "tn", 1);
	dcellmm(&u2, G0T, z, "nn", 1);
	for(int i=0; i<2; i++){
		assert_dmat_near(P(u, i), P(u2, i), TOL*100);
	}
	//spcell2sp concatenates the blocks; it agrees with spcellfull
	dsp *Gs=dspcell2sp(G0);
	assert_int_equal(NX(Gs), nr[0]+nr[1]);
	assert_int_equal(NY(Gs), nc[0]+nc[1]);
	dcell *Gf=NULL;
	dspcellfull(&Gf, G0, 'n', 1);
	dmat *Gm=dcell2m(Gf);
	assert_sp_dense(Gs, Gm, 0);
	//spcellfull with transpose
	dcell *Gft=NULL;
	dspcellfull(&Gft, G0, 't', 1);
	dmat *Gmt=dcell2m(Gft);
	dmat *Gmt2=dtrans(Gm);
	assert_dmat_near(Gmt, Gmt2, 0);
	//dcellscale/spcellscale and spcellsym
	dspcellscale(A, 2);
	assert_float_equal(dspcellwdot(x, A, x), 2*res3, TOL*res3*20);
	dfree(Gmt); dfree(Gmt2); dcellfree(Gft);
	dfree(Gm); dcellfree(Gf); dspfree(Gs);
	dcellfree(u); dcellfree(u2);
	dcellfree(x); dcellfree(y); dcellfree(z);
	dspcellfree(A); dspcellfree(A2);
	dspcellfree(G0T); dspcellfree(G0);
}

/*Element (i,j) of op(M) where op is n: none, N: conjugate, t: transpose, c: conjugate transpose.*/
static comp opc(const cmat *M, char t, long i, long j){
	switch(t){
	case 'n': return P(M, i, j);
	case 'N': return conj(P(M, i, j));
	case 't': return P(M, j, i);
	default: return conj(P(M, j, i));
	}
}
static int istrans(char t){
	return t=='t'||t=='c';
}
/*Random complex dense matrix.*/
static cmat *crand(long nx, long ny, rand_t *rstat){
	cmat *A=cnew(nx, ny);
	crandn(A, 1, rstat);
	return A;
}
/** Complex cspmm and cmulsp with all combinations of n, N, t, c against a naive product. */
static void sparse_complex_mm(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 21);
	const long m=5, k=4, n=3;
	const char ops[4]={'n', 'N', 't', 'c'};
	const comp alpha=COMPLEX(0.5, -1.5);
	for(int ia=0; ia<4; ia++){
		for(int ix=0; ix<4; ix++){
			const char ta=ops[ia], tx=ops[ix];
			const char tr[3]={ta, tx, 0};
			//y=op(A)*op(x): A is m*k after op, x is k*n after op
			csp *A=istrans(ta)?cspnewrandu(k, m, COMPLEX(1., 1.), 0.6, &rstat):cspnewrandu(m, k, COMPLEX(1., 1.), 0.6, &rstat);
			cmat *x=istrans(tx)?crand(n, k, &rstat):crand(k, n, &rstat);
			cmat *Af=NULL;
			cspfull(&Af, A, 'n', 1);
			cmat *y=crand(m, n, &rstat);
			cmat *y0=cdup(y);
			cspmm(&y, A, x, tr, alpha);
			for(long j=0; j<n; j++){
				for(long i=0; i<m; i++){
					comp ref=P(y0, i, j);
					for(long l=0; l<k; l++) ref+=alpha*opc(Af, ta, i, l)*opc(x, tx, l, j);
					assert_float_equal(cabs(P(y, i, j)-ref), 0, 1e-12);
				}
			}
			cfree(y); cfree(y0); cfree(x); cfree(Af); cspfree(A);
			if(ta=='N'||tx=='N') continue;//mulsp only supports n, t, c
			//y=op(x)*op(A): x is m*k after op, A is k*n after op
			A=istrans(tx)?cspnewrandu(n, k, COMPLEX(1., -1.), 0.6, &rstat):cspnewrandu(k, n, COMPLEX(1., -1.), 0.6, &rstat);
			x=istrans(ta)?crand(k, m, &rstat):crand(m, k, &rstat);
			cspfull(&Af, A, 'n', 1);
			y=NULL;
			cmulsp(&y, x, A, tr, alpha);
			assert_int_equal(NX(y), m);
			assert_int_equal(NY(y), n);
			for(long j=0; j<n; j++){
				for(long i=0; i<m; i++){
					comp ref=0;
					for(long l=0; l<k; l++) ref+=alpha*opc(x, ta, i, l)*opc(Af, tx, l, j);
					assert_float_equal(cabs(P(y, i, j)-ref), 0, 1e-12);
				}
			}
			cfree(y); cfree(x); cfree(Af); cspfree(A);
		}
	}
	//spmv with conjugate transpose and a row vector x is treated as a vector
	csp *A=cspnewrandu(4, 3, COMPLEX(1., 1.), 0.7, &rstat);
	cmat *Af=NULL;
	cspfull(&Af, A, 'n', 1);
	cmat *x=crand(4, 1, &rstat);
	cmat *y=cnew(3, 1);
	cspmv(y, A, x, 'c', 2);
	for(long i=0; i<3; i++){
		comp ref=0;
		for(long l=0; l<4; l++) ref+=2.*conj(P(Af, l, i))*P(x, l);
		assert_float_equal(cabs(P(y, i)-ref), 0, 1e-12);
	}
	cmat *xr=crand(1, 3, &rstat);//row vector
	cmat *yr=NULL;
	cspmm(&yr, A, xr, "nn", 1);
	assert_int_equal(NX(yr), 4);
	assert_int_equal(NY(yr), 1);
	for(long i=0; i<4; i++){
		comp ref=0;
		for(long l=0; l<3; l++) ref+=P(Af, i, l)*P(xr, l);
		assert_float_equal(cabs(P(yr, i)-ref), 0, 1e-12);
	}
	cfree(yr);
	cmat *xc=crand(1, 4, &rstat);//row vector for transpose
	cspmm(&yr, A, xc, "tt", 1);
	assert_int_equal(NX(yr), 3);
	cfree(yr); cfree(xc); cfree(xr);
	//sparse times sparse. csptrans is the conjugate transpose.
	csp *B=cspnewrandu(4, 2, COMPLEX(1., 2.), 0.7, &rstat);
	cmat *Bf=NULL;
	cspfull(&Bf, B, 'n', 1);
	csp *Ah=csptrans(A);
	csp *C=cspmulsp(Ah, B, "nn");
	cspfree(Ah);
	cmat *Cf=NULL;
	cspfull(&Cf, C, 'n', 1);
	for(long j=0; j<2; j++){
		for(long i=0; i<3; i++){
			comp ref=0;
			for(long l=0; l<4; l++) ref+=conj(P(Af, l, i))*P(Bf, l, j);
			assert_float_equal(cabs(P(Cf, i, j)-ref), 0, 1e-12);
		}
	}
	cfree(Cf); cspfree(C);
	//cspmulsp with all n/t/c combinations against the dense reference: op(L)*op(R) is 3x2.
	{
		const char ts[3]={'n', 't', 'c'};
		for(int ia=0; ia<3; ia++){
			for(int ib=0; ib<3; ib++){
				const char tr[3]={ts[ia], ts[ib], 0};
				csp *L=istrans(tr[0])?cspnewrandu(4, 3, COMPLEX(1., 2.), 0.7, &rstat):cspnewrandu(3, 4, COMPLEX(1., 2.), 0.7, &rstat);
				csp *R=istrans(tr[1])?cspnewrandu(2, 4, COMPLEX(-1., 1.), 0.7, &rstat):cspnewrandu(4, 2, COMPLEX(-1., 1.), 0.7, &rstat);
				cmat *Lf=NULL, *Rf=NULL, *Of=NULL;
				cspfull(&Lf, L, 'n', 1);
				cspfull(&Rf, R, 'n', 1);
				csp *O=cspmulsp(L, R, tr);
				cspfull(&Of, O, 'n', 1);
				assert_int_equal(NX(Of), 3);
				assert_int_equal(NY(Of), 2);
				for(long j=0; j<2; j++){
					for(long i=0; i<3; i++){
						comp ref=0;
						for(long l=0; l<4; l++) ref+=opc(Lf, tr[0], i, l)*opc(Rf, tr[1], l, j);
						assert_float_equal(cabs(P(Of, i, j)-ref), 0, 1e-12);
					}
				}
				cfree(Lf); cfree(Rf); cfree(Of); cspfree(L); cspfree(R); cspfree(O);
			}
		}
	}
	//A*conj(E') with scaling, accumulated twice
	csp *E=cspnewrandu(2, 3, COMPLEX(-1., 1.), 0.7, &rstat);
	cmat *Ef=NULL;
	cspfull(&Ef, E, 'n', 1);
	csp *Eh=csptrans(E);
	csp *C2=NULL;
	cspmulsp2(&C2, A, Eh, "nn", 2);
	cspmulsp2(&C2, A, Eh, "nn", 2);
	cspfree(Eh);
	cspfull(&Cf, C2, 'n', 1);
	for(long j=0; j<2; j++){
		for(long i=0; i<4; i++){
			comp ref=0;
			for(long l=0; l<3; l++) ref+=4.*P(Af, i, l)*conj(P(Ef, j, l));
			assert_float_equal(cabs(P(Cf, i, j)-ref), 0, 1e-12);
		}
	}
	cfree(Cf); cfree(Ef); cspfree(E); cspfree(C2);
	//cspconj
	csp *Ac=cspdup(A);
	cspconj(Ac);
	for(long i=0; i<A->pp[A->ny]; i++) assert_float_equal(cabs(Ac->px[i]-conj(A->px[i])), 0, 0);
	cspfree(Ac);
	//trace of complex sparse uses the magnitude
	csp *D=cspnewdiag(3, NULL, COMPLEX(3., 4.));
	assert_float_equal(csptrace(D, 1), 15, 1e-12);
	assert_float_equal(csptrace(D, 2), 75, 1e-10);
	cspfree(D);
	cfree(Bf); cspfree(B);
	cfree(y); cfree(x); cfree(Af); cspfree(A);
}
/** spmulcreal, spref, spnew2, spclean, spfull of dense input and spconvolvop with many points. */
static void sparse_misc2(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 22);
	dsp *A=randsp(5, 4, 0.5, &rstat);
	dmat *Af=full(A);
	//y+=alpha*A*real(x)
	comp x[4];
	for(int i=0; i<4; i++) x[i]=COMPLEX((double)i-1.5, 10.*i);
	real y[5]={1, 1, 1, 1, 1};
	dspmulcreal(y, A, x, 2);
	for(int i=0; i<5; i++){
		real ref=1;
		for(int l=0; l<4; l++) ref+=2*P(Af, i, l)*creal(x[l]);
		assert_float_equal(y[i], ref, TOL);
	}
	dspmulcreal(NULL, A, x, 2);
	//spref shares data and counts references
	dsp *Ar=dspref(A);
	assert_ptr_equal(Ar->px, A->px);
	dspfree(A);
	dmat *Arf=full(Ar);
	assert_dmat_near(Arf, Af, 0);
	dfree(Arf);
	assert_null(dspref(NULL));
	//spnew2 has the same shape and capacity
	dsp *B=dspnew2(Ar);
	assert_int_equal(NX(B), 5);
	assert_int_equal(NY(B), 4);
	assert_int_equal(B->nzmax, Ar->pp[Ar->ny]);
	dspfree(B);
	assert_null(dspnew2(NULL));
	//spclean drops explicit zeros
	dsp *Cl=dspdup(Ar);
	long nz=Cl->pp[Cl->ny];
	Cl->px[0]=0;
	dspclean(Cl);
	assert_int_equal(Cl->pp[Cl->ny], nz-1);
	dspclean(NULL);
	dspfree(Cl);
	//spfull also accepts dense input
	dmat *out=NULL;
	dspfull(&out, (dsp *)Af, 'n', 2);
	assert_float_equal(dsum(out), 2*dsum(Af), TOL);
	dmat *outt=NULL;
	dspfull(&outt, (dsp *)Af, 't', 1);
	assert_int_equal(NX(outt), 4);
	assert_float_equal(P(outt, 2, 3), P(Af, 3, 2), 0);
	dfree(out); dfree(outt);
	//spconvolvop with more than 10 coupled points
	dmat *K=dnew(4, 4);
	for(int i=0; i<12; i++) P(K, i)=i+1;
	dsp *Cv=dspconvolvop(K);
	assert_int_equal(NX(Cv), 16);
	assert_int_equal(Cv->pp[16], 16*12);
	dmat *Cvf=full(Cv);
	//C(i,j)=K(i-j) with periodic wrapping
	for(long jy=0; jy<4; jy++) for(long jx=0; jx<4; jx++){
		long j=jx+jy*4;
		for(long iy=0; iy<4; iy++) for(long ix=0; ix<4; ix++){
			long i=ix+iy*4;
			assert_float_equal(P(Cvf, i, j), P(K, (ix-jx+4)%4, (iy-jy+4)%4), 0);
		}
	}
	dfree(Cvf); dspfree(Cv); dfree(K);
	//spnewrandu that needs to grow its storage
	dsp *R=dspnewrandu(50, 50, 1, 0.01, &rstat);
	assert_true(R->pp[R->ny]>0);
	dspfree(R);
	R=dspnewrandu(10, 10, 1, 2, &rstat);//fill is clipped to 1
	assert_int_equal(R->pp[R->ny], 100);
	dspfree(R);
	//reverse row permutation
	long perm[5]={4, 2, 0, 1, 3};
	dsp *Pr=dspperm(Ar, 1, perm, NULL);
	dmat *Prf=full(Pr);
	for(long j=0; j<4; j++) for(long i=0; i<5; i++) assert_float_equal(P(Prf, perm[i], j), P(Af, i, j), 0);
	dfree(Prf); dspfree(Pr);
	dspfree(Ar);
	dfree(Af);
}

UNITTEST_MAIN(
	cmocka_unit_test(sparse_create),
	cmocka_unit_test(sparse_spmv),
	cmocka_unit_test(sparse_spmm),
	cmocka_unit_test(sparse_spmulsp),
	cmocka_unit_test(sparse_trans),
	cmocka_unit_test(sparse_sum),
	cmocka_unit_test(sparse_add),
	cmocka_unit_test(sparse_cell),
	cmocka_unit_test(sparse_complex_mm),
	cmocka_unit_test(sparse_misc2),
)
