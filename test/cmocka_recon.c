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
  Tests of the reconstruction building blocks in lib/: preconditioned
  conjugate gradient (pcg.c), sparse plus low rank operators and solvers
  (muv.c), laplacian regularization (laplacian.c), and actuator slaving,
  stuck/floating actuator handling and extrapolation (slaving.c).
*/
#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-9
#endif
static const long NB1=8, NB2=6;/*block sizes*/
/**
   Dense SPD and diagonally dominant matrix of size n with seed.
*/
static dmat *mkspd(long n, int seed){
	rand_t rstat;
	seed_rand(&rstat, seed);
	dmat *R=dnew(n, n);
	drandn(R, 1, &rstat);
	dmat *A=NULL;
	dmm(&A, 0, R, R, "nt", 0.2);
	for(long i=0; i<n; i++){
		P(A, i, i)+=n;
	}
	dfree(R);
	return A;
}
static dmat *mkrand(long nx, long ny, int seed){
	rand_t rstat;
	seed_rand(&rstat, seed);
	dmat *A=dnew(nx, ny);
	drandn(A, 1, &rstat);
	return A;
}
/**
   Split a full matrix into a 2x2 block dcell (dense) with block sizes NB1, NB2.
*/
static dcell *split_dense(const dmat *A){
	dcell *out=dcellnew(2, 2);
	long n[2]={NB1, NB2};
	long o[2]={0, NB1};
	for(int iy=0; iy<2; iy++){
		for(int ix=0; ix<2; ix++){
			P(out, ix, iy)=dsub(A, o[ix], n[ix], o[iy], n[iy]);
		}
	}
	return out;
}
/**
   Split a full matrix into a 2x2 block dspcell (sparse).
*/
static dspcell *split_sparse(const dmat *A){
	dcell *tmp=split_dense(A);
	dspcell *out=dspcellnew(2, 2);
	for(int i=0; i<4; i++){
		P(out, i)=d2sp(P(tmp, i), 0);
	}
	dcellfree(tmp);
	return out;
}
/**
   Split a vector into a 2x1 dcell.
*/
static dcell *split_vec(const dmat *x){
	dcell *out=dcellnew(2, 1);
	P(out, 0)=dsub(x, 0, NB1, 0, 1);
	P(out, 1)=dsub(x, NB1, NB2, 0, 1);
	return out;
}
/**
   The test problem: M is SPD, U is a single low rank column, and V=-U so that
   the full operator M-U*V'=M+U*U' remains SPD.
*/
typedef struct{
	dmat *M;   /**<full sparse part*/
	dmat *U;   /**<full U*/
	dmat *Afull;/**<M-UV'*/
	dmat *b;   /**<right hand side*/
	dmat *xref;/**<Afull\b*/
}problem_t;
static void problem_init(problem_t *p){
	long n=NB1+NB2;
	p->M=mkspd(n, 1);
	p->U=mkrand(n, 1, 2);
	dscale(p->U, 0.5);
	p->Afull=ddup(p->M);
	dmm(&p->Afull, 1, p->U, p->U, "nt", 1);
	p->b=mkrand(n, 1, 3);
	dmat *Ai=dinv(p->Afull);
	p->xref=NULL;
	dmm(&p->xref, 0, Ai, p->b, "nn", 1);
	dfree(Ai);
}
static void problem_free(problem_t *p){
	dfree(p->M); dfree(p->U); dfree(p->Afull); dfree(p->b); dfree(p->xref);
}
/**
   Fill muv_t with the problem. sparse: M as dspcell, otherwise dense dcell.
*/
static void problem_muv(muv_t *A, const problem_t *p, int sparse, int lowrank){
	memset(A, 0, sizeof(muv_t));
	if(lowrank){
		A->M=sparse?(cell *)split_sparse(p->M):(cell *)split_dense(p->M);
		A->U=split_vec(p->U);
		dmat *V=ddup(p->U); dscale(V, -1);
		A->V=split_vec(V);
		dfree(V);
	} else{
		A->M=sparse?(cell *)split_sparse(p->Afull):(cell *)split_dense(p->Afull);
	}
}
static void assert_cell_near_vec(const dcell *x, const dmat *ref, real tol){
	dmat *xm=dcell2m(x);
	reshape(xm, PN(xm), 1);
	assert_dmat_near(xm, ref, tol*(1+dmaxabs(ref)));
	dfree(xm);
}
/** Diagonal (Jacobi) preconditioner. M is a dcell of the inverse diagonal of each block. */
static void jacobi_pre(dcell **xout, const void *M, const dcell *xin){
	const dcell *D=(const dcell *)M;
	if(!*xout) *xout=dcellnew2(xin);
	for(int ib=0; ib<NX(xin); ib++){
		for(long i=0; i<PN(P(xin, ib)); i++){
			P(P(*xout, ib), i)=P(P(D, ib), i)*P(P(xin, ib), i);
		}
	}
}
/** Function form of the operator. A is a dense full matrix. */
static void fun_mul(dcell **xout, const void *A, const dcell *xin, const real alpha){
	const dmat *Af=(const dmat *)A;
	dmat *x=dcell2m(xin); reshape(x, PN(x), 1);
	dmat *y=NULL;
	dmm(&y, 0, Af, x, "nn", alpha);
	if(!*xout) *xout=dcellnew2(xin);
	dmat *yo=dcell2m(*xout); reshape(yo, PN(yo), 1);
	dadd(&yo, 1, y, 1);
	d2cell(xout, yo, xin);
	dfree(x); dfree(y); dfree(yo);
}
static void fun_mult(dcell **xout, const void *A, const dcell *xin, const real alpha){
	dmat *At=dtrans((const dmat *)A);
	fun_mul(xout, At, xin, alpha);
	dfree(At);
}

/**
   muv() applies (M-UV')*x*alpha and accumulates; muv_trans() applies
   (M'-VU')*x*alpha for a nonsymmetric M with independent U and V.
*/
static void recon_muv_apply(void **state){
	(void)state;
	long n=NB1+NB2;
	dmat *M=mkrand(n, n, 11);
	dmat *U=mkrand(n, 2, 12);
	dmat *V=mkrand(n, 2, 13);
	dmat *x=mkrand(n, 1, 14);
	dmat *Afull=ddup(M);
	dmm(&Afull, 1, U, V, "nt", -1);
	muv_t A={0};
	A.M=(cell *)split_dense(M);
	dcell *Uc=dcellnew(2, 1), *Vc=dcellnew(2, 1);
	P(Uc, 0)=dsub(U, 0, NB1, 0, 2); P(Uc, 1)=dsub(U, NB1, NB2, 0, 2);
	P(Vc, 0)=dsub(V, 0, NB1, 0, 2); P(Vc, 1)=dsub(V, NB1, NB2, 0, 2);
	A.U=Uc; A.V=Vc;
	dcell *xc=split_vec(x);
	dmat *yref=NULL;
	dmm(&yref, 0, Afull, x, "nn", 2.);
	dcell *y=NULL;
	muv(&y, &A, xc, 2.);
	assert_cell_near_vec(y, yref, TOL);
	/*accumulates*/
	muv(&y, &A, xc, 2.);
	dscale(yref, 2);
	assert_cell_near_vec(y, yref, TOL);
	dcellfree(y);
	/*transpose*/
	dfree(yref);
	dmm(&yref, 0, Afull, x, "tn", -1.);
	muv_trans(&y, &A, xc, -1.);
	assert_cell_near_vec(y, yref, TOL);
	dcellfree(y);
	/*muv_trans with NULL input does nothing*/
	muv_trans(&y, &A, NULL, 1);
	assert_null(y);
	/*block application with muv_ib: output block xb gets A(xb,yb)*x(yb). The
	  argument has the layout of the private MUV_IB_T struct in muv.c.*/
	struct{const muv_t *A; int xb; int yb;} B={&A, 0, 1};
	dcell *x1=dcellnew2(xc);
	dcp(&P(x1, 1), P(xc, 1));/*only block 1 nonzero*/
	dcell *yfull=NULL;
	muv(&yfull, &A, x1, 1);
	muv_ib(&y, &B, xc, 1);
	assert_non_null(y);
	assert_dmat_near(P(y, 0), P(yfull, 0), TOL*10);
	assert_null(P(y, 1));
	dcellfree(y);
	/*xb=-1 falls back to full muv*/
	B.xb=-1;
	muv_ib(&y, &B, xc, 1);
	dfree(yref);
	dmm(&yref, 0, Afull, x, "nn", 1.);
	assert_cell_near_vec(y, yref, TOL);
	dcellfree(y);
	muv_free(&A);
	assert_null(A.M);
	assert_null(A.U);
	dcellfree(x1); dcellfree(yfull); dcellfree(xc);
	dfree(M); dfree(U); dfree(V); dfree(x); dfree(Afull); dfree(yref);
}
/**
   muv() and muv_trans() with function operators (Mfun, Mtfun) and the error
   when muv_trans() has no transpose function.
*/
static void recon_muv_fun(void **state){
	(void)state;
	long n=NB1+NB2;
	dmat *Af=mkrand(n, n, 21);
	dmat *x=mkrand(n, 1, 22);
	dcell *xc=split_vec(x);
	muv_t A={0};
	A.Mfun=fun_mul;
	A.Mdata=Af;
	dcell *y=NULL;
	muv(&y, &A, xc, 3);
	dmat *yref=NULL;
	dmm(&yref, 0, Af, x, "nn", 3);
	assert_cell_near_vec(y, yref, TOL);
	dcellfree(y);
	/*no Mtfun: error() and nothing done*/
	muv_trans(&y, &A, xc, 1);
	assert_null(y);
	A.Mtfun=fun_mult;
	muv_trans(&y, &A, xc, 1);
	dmm(&yref, 0, Af, x, "tn", 1);
	assert_cell_near_vec(y, yref, TOL);
	dcellfree(y); dcellfree(xc); dfree(x); dfree(Af); dfree(yref);
}
/**
   pcg() solves an SPD system to machine precision in n iterations, with and
   without a preconditioner, honors warm restart, handles a zero right hand
   side, and returns the residual after maxiter.
*/
static void recon_pcg(void **state){
	(void)state;
	problem_t p={0};
	problem_init(&p);
	const long n=NB1+NB2;
	muv_t A;
	problem_muv(&A, &p, 0, 1);
	dcell *b=split_vec(p.b);
	dcell *x=NULL;
	real res=pcg(&x, muv, &A, NULL, NULL, b, 0, n);
	assert_cell_near_vec(x, p.xref, 1e3*TOL);
	assert_true(res<1e-6);
	/*Jacobi preconditioner*/
	dcell *D=dcellnew(2, 1);
	for(int ib=0; ib<2; ib++){
		long nb=ib==0?NB1:NB2;
		long off=ib==0?0:NB1;
		P(D, ib)=dnew(nb, 1);
		for(long i=0; i<nb; i++){
			P(P(D, ib), i)=1./P(p.Afull, i+off, i+off);
		}
	}
	dcell *x2=NULL;
	pcg(&x2, muv, &A, jacobi_pre, D, b, 0, n);
	assert_cell_near_vec(x2, p.xref, 1e3*TOL);
	/*few iterations: residual is reported and decreases*/
	dcell *x3=NULL;
	real res1=pcg(&x3, muv, &A, NULL, NULL, b, 0, 1);
	dcellzero(x3);
	real res3=pcg(&x3, muv, &A, NULL, NULL, b, 0, 3);
	assert_true(res1>0&&res1<1);
	assert_true(res3<res1);
	/*warm restart from the converged solution stays converged after 1 iteration*/
	dcellcp(&x3, x);
	pcg(&x3, muv, &A, NULL, NULL, b, 1, 1);
	assert_cell_near_vec(x3, p.xref, 1e3*TOL);
	/*cold start zeroes the input*/
	pcg(&x3, muv, &A, NULL, NULL, b, 0, 1);
	dmat *x3m=dcell2m(x3);
	dmat *xm=dcell2m(x);
	assert_true(ddiff(x3m, xm)>1e-6);
	dfree(x3m); dfree(xm);
	/*zero rhs: terminates early with zero solution and zero residual*/
	dcell *b0=dcellnew2(b);
	dcell *x0=NULL;
	real r0=pcg(&x0, muv, &A, NULL, NULL, b0, 0, 5);
	assert_float_equal(r0, 0, 0);
	assert_float_equal(dcellnorm(x0), 0, 0);
	dcellfree(b0); dcellfree(x0);
	dcellfree(x); dcellfree(x2); dcellfree(x3); dcellfree(D); dcellfree(b);
	muv_free(&A);
	problem_free(&p);
}
/**
   Direct solvers: Cholesky (with low rank Woodbury correction) and SVD agree with
   the dense solution; muv_solve() with alg 0, 1, 2 and with a right hand side
   operator; in place solve; muv_direct_spsolve() computes A^-1 X A^-T.
*/
static void recon_muv_direct(void **state){
	(void)state;
	problem_t p={0};
	problem_init(&p);
	const long n=NB1+NB2;
	dcell *b=split_vec(p.b);
	for(int alg=0; alg<3; alg++){
		muv_t A;
		problem_muv(&A, &p, 1, 1);
		A.alg=alg;
		A.maxit=n;
		if(alg==0) muv_direct_prep(&A, 0);
		if(alg==2) muv_direct_prep(&A, 1);
		if(alg==0){
			assert_non_null(A.C);
			assert_non_null(A.Up);
			assert_non_null(A.Vp);
		} else if(alg==2){
			assert_non_null(A.MI);
			assert_null(A.U);/*absorbed into MI*/
		}
		dcell *x=NULL;
		muv_solve(&x, &A, NULL, b);
		assert_cell_near_vec(x, p.xref, 1e3*TOL);
		dcellfree(x);
		muv_free(&A);
	}
	/*single block dcell and in-place mat solve*/
	{
		muv_t A={0};
		A.M=(cell *)dspcellnew(1, 1);
		P(A.M, 0)=(cell *)d2sp(p.Afull, 0);
		muv_direct_prep(&A, 0);
		dcell *b1=dcellnew(1, 1); P(b1, 0)=ddup(p.b);
		dcell *x=NULL;
		muv_direct_solve(&x, &A, b1);
		assert_dmat_near(P(x, 0), p.xref, 1e3*TOL);
		dmat *bb=ddup(p.b);
		muv_direct_solve_mat(&bb, &A, bb);
		assert_dmat_near(bb, p.xref, 1e3*TOL);
		/*NULL input does nothing*/
		dcell *xn=NULL;
		muv_direct_solve(&xn, &A, NULL);
		assert_null(xn);
		/*right hand side operator R=2*I: x=A\(2b)*/
		muv_t R={0};
		R.M=(cell *)dspcellnew(1, 1);
		dmat *two=dnew(n, 1); dset(two, 2);
		P(R.M, 0)=(cell *)dspnewdiag(n, P(two), 1);
		dfree(two);
		A.alg=0;
		dcellfree(x);
		muv_solve(&x, &A, &R, b1);
		dmat *x2=ddup(p.xref); dscale(x2, 2);
		assert_dmat_near(P(x, 0), x2, 1e3*TOL);
		dfree(x2);
		/*spsolve: A^-1*X*A^-T with chol (sparse output) and svd (dense output)*/
		dmat *X=mkrand(n, n, 31);
		for(long i=0; i<PN(X); i++){
			if(fabs(P(X, i))<1) P(X, i)=0;
		}
		dsp *Xs=d2sp(X, 0);
		dmat *Ai=dinv(p.Afull);
		dmat *tmp=NULL, *ref=NULL;
		dmm(&tmp, 0, Ai, X, "nn", 1);
		dmm(&ref, 0, tmp, Ai, "nt", 1);
		dsp *outs=(dsp *)muv_direct_spsolve(&A, Xs);
		dmat *outd=NULL;
		dspfull(&outd, outs, 'n', 1);
		assert_dmat_near(outd, ref, 1e3*TOL);
		dspfree(outs); dfree(outd);
		muv_direct_prep(&A, 1e-12);/*svd with explicit threshold*/
		assert_null(A.C);
		dmat *outm=(dmat *)muv_direct_spsolve(&A, Xs);
		assert_dmat_near(outm, ref, 1e3*TOL);
		dmat *bb2=ddup(p.b);
		muv_direct_solve_mat(&bb2, &A, bb2);/*in place with MI*/
		assert_dmat_near(bb2, p.xref, 1e3*TOL);
		dfree(bb2);
		dfree(outm); dfree(X); dspfree(Xs); dfree(Ai); dfree(tmp); dfree(ref);
		dfree(bb); dcellfree(x); dcellfree(b1);
		muv_free(&A);
		muv_free(&R);
	}
	/*invalid algorithm is rejected*/
	{
		muv_t A;
		problem_muv(&A, &p, 0, 0);
		A.alg=5;
		dcell *x=NULL;
		muv_solve(&x, &A, NULL, b);
		assert_null(x);
		muv_free(&A);
	}
	dcellfree(b);
	problem_free(&p);
}
/**
   Block Gauss-Seidel solver converges to the dense solution for a diagonally
   dominant matrix with Cholesky, CG, or SVD block solvers, with and without
   low rank terms.
*/
static void recon_muv_bgs(void **state){
	(void)state;
	problem_t p={0};
	problem_init(&p);
	dcell *b=split_vec(p.b);
	for(int alg=0; alg<3; alg++){
		for(int lowrank=0; lowrank<2; lowrank++){
			muv_t A;
			problem_muv(&A, &p, alg!=2, lowrank);
			A.alg=alg;
			A.bgs=30;
			A.maxit=20;
			A.warm=0;
			if(alg==0) muv_direct_diag_prep(&A, 0);
			if(alg==2) muv_direct_diag_prep(&A, 1);
			if(alg!=1){
				assert_int_equal(A.nb, 2);
				if(lowrank){
					assert_non_null(A.UpB);
				}
			}
			dcell *x=NULL;
			muv_solve(&x, &A, NULL, b);
			assert_cell_near_vec(x, p.xref, 1e-6);
			/*block solve agrees with dense block inverse*/
			if(alg!=1&&!lowrank){
				dmat *blk=dsub(p.Afull, 0, NB1, 0, NB1);
				dmat *bi=dinv(blk);
				dmat *rhs=mkrand(NB1, 1, 41);
				dmat *y=NULL, *yref=NULL;
				muv_direct_diag_solve(&y, &A, rhs, 0);
				dmm(&yref, 0, bi, rhs, "nn", 1);
				assert_dmat_near(y, yref, 1e3*TOL);
				muv_direct_diag_solve(&rhs, &A, rhs, 0);/*in place*/
				assert_dmat_near(rhs, yref, 1e3*TOL);
				dmat *yn=NULL;
				muv_direct_diag_solve(&yn, &A, NULL, 0);
				assert_null(yn);
				dfree(blk); dfree(bi); dfree(rhs); dfree(y); dfree(yref);
			}
			dcellfree(x);
			muv_free(&A);
		}
	}
	dcellfree(b);
	problem_free(&p);
}
/**
   Laplacian on a periodic map: coefficient formula, constant is in the null
   space, a cosine is an eigenfunction with known eigenvalue, and the sparse
   matrix agrees with direct application.
*/
static void recon_laplacian_map(void **state){
	(void)state;
	const real r0=0.2, wt=0.5, dx=1./16;
	const real cf=laplacian_coef(r0, wt, dx);
	const real cfref=(2.*M_PI/0.5e-6)*sqrt(pow(r0/dx, 5./3.)/3.44/wt)*pow(2.-pow(2, -1./3.)-pow(2, -1./6.), -0.5);
	assert_float_equal(cf, cfref, 1e-12*cfref);
	/*scales as weight^-1/2*/
	assert_float_equal(laplacian_coef(r0, 4*wt, dx), cf*0.5, 1e-12*cf);
	const long nx=8, ny=6;
	dmat *opd=dnew(nx, ny);
	dmat *out=dnew(nx, ny);
	/*constant*/
	dset(opd, 3);
	apply_laplacian_map(out, opd, dx, r0, wt);
	assert_float_equal(dmaxabs(out), 0, 1e-12*cf);
	/*cosine along x with k=2 periods, and along y with 1 period*/
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			P(opd, ix, iy)=cos(2*M_PI*2*ix/nx)+cos(2*M_PI*iy/ny);
		}
	}
	apply_laplacian_map(out, opd, dx, r0, wt);
	const real ex=0.5*(cos(2*M_PI*2/nx)-1)*cf;
	const real ey=0.5*(cos(2*M_PI/ny)-1)*cf;
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			real ref=ex*cos(2*M_PI*2*ix/nx)+ey*cos(2*M_PI*iy/ny);
			assert_float_equal(P(out, ix, iy), ref, 1e-10*cf);
		}
	}
	/*matrix form*/
	dsp *L=mklaplacian_map(nx, ny, dx, r0, wt);
	assert_int_equal(NX(L), nx*ny);
	assert_int_equal(NY(L), nx*ny);
	assert_int_equal(L->pp[NY(L)], nx*ny*5);
	dmat *v=dref_reshape(opd, nx*ny, 1);
	dmat *y=NULL;
	dspmm(&y, L, v, "nn", 1);
	dmat *outv=dref_reshape(out, nx*ny, 1);
	assert_dmat_near(y, outv, 1e-10*cf);
	/*symmetric*/
	dsp *Lt=dsptrans(L);
	dsp *diff=dspadd2(L, 1, Lt, -1);
	dmat *fd=NULL;
	dspfull(&fd, diff, 'n', 1);
	assert_float_equal(dmaxabs(fd), 0, 1e-12*cf);
	dfree(fd); dspfree(diff); dspfree(Lt);
	dfree(y); dfree(v); dfree(outv); dspfree(L); dfree(opd); dfree(out);
}
/**
   Laplacian on a loc grid: interior stencil, partial stencils with their
   coefficients on edges and corners, zero row sums, linear functions vanish in
   the interior, and degenerate points get zero weight.
*/
static void recon_laplacian_loc(void **state){
	(void)state;
	const real r0=0.3, wt=1, dx=0.5;
	const long n=5;
	loc_t *loc=mksqloc(n, n, dx, dx, 0, 0);
	dsp *L=mklaplacian_loc(loc, r0, wt);
	dmat *Lf=NULL;
	dspfull(&Lf, L, 'n', 1);
	const real cf=laplacian_coef(r0, wt, dx);
	const real cf4=(2.*M_PI/0.5e-6)*sqrt(pow(r0/dx, 5./3.)/3.44/wt)*pow(2.-pow(2, -4./3.)-pow(2, -1./6.), -0.5);
	const real cf3=(2.*M_PI/0.5e-6)*sqrt(pow(r0/dx, 5./3.)/3.44/wt)*pow(2.-pow(2, -1./6.), -0.5);
#define IND(ix,iy) ((ix)+(iy)*n)
	assert_float_equal(loc->locx[IND(3, 2)], 3*dx, 1e-12);
	assert_float_equal(loc->locy[IND(3, 2)], 2*dx, 1e-12);
	/*interior point*/
	const long ic=IND(2, 2);
	assert_float_equal(P(Lf, ic, ic), -cf, 1e-10*cf);
	assert_float_equal(P(Lf, ic, IND(1, 2)), 0.25*cf, 1e-10*cf);
	assert_float_equal(P(Lf, ic, IND(3, 2)), 0.25*cf, 1e-10*cf);
	assert_float_equal(P(Lf, ic, IND(2, 1)), 0.25*cf, 1e-10*cf);
	assert_float_equal(P(Lf, ic, IND(2, 3)), 0.25*cf, 1e-10*cf);
	/*bottom edge point: missing down neighbor, up neighbor doubled*/
	const long ie=IND(2, 0);
	assert_float_equal(P(Lf, ie, ie), -cf4, 1e-10*cf);
	assert_float_equal(P(Lf, ie, IND(2, 1)), 0.5*cf4, 1e-10*cf);
	assert_float_equal(P(Lf, ie, IND(1, 0)), 0.25*cf4, 1e-10*cf);
	/*corner point*/
	const long io=IND(0, 0);
	assert_float_equal(P(Lf, io, io), -cf3, 1e-10*cf);
	assert_float_equal(P(Lf, io, IND(1, 0)), 0.5*cf3, 1e-10*cf);
	assert_float_equal(P(Lf, io, IND(0, 1)), 0.5*cf3, 1e-10*cf);
	/*row sums are zero; linear function vanishes in the interior*/
	dmat *ones=dnew(n*n, 1); dset(ones, 1);
	dmat *y=NULL;
	dspmm(&y, L, ones, "nn", 1);
	assert_float_equal(dmaxabs(y), 0, 1e-10*cf);
	dmat *lin=dnew(n*n, 1);
	for(long i=0; i<n*n; i++){
		P(lin, i)=2*loc->locx[i]-3*loc->locy[i];
	}
	dzero(y);
	dspmm(&y, L, lin, "nn", 1);
	for(long iy=1; iy<n-1; iy++){
		for(long ix=1; ix<n-1; ix++){
			assert_float_equal(P(y, IND(ix, iy)), 0, 1e-8*cf);
		}
	}
	/*edge does not vanish for a gradient normal to it*/
	assert_true(fabs(P(y, IND(2, 0)))>1e-3*cf);
#undef IND
	dfree(ones); dfree(lin); dfree(y); dfree(Lf); dspfree(L); locfree(loc);
	/*two isolated adjacent points: one neighbor only, weight is zero*/
	loc_t *loc2=locnew(2, dx, dx);
	loc2->locx[0]=0; loc2->locy[0]=0;
	loc2->locx[1]=dx; loc2->locy[1]=0;
	dsp *L2=mklaplacian_loc(loc2, r0, wt);
	assert_int_equal(NX(L2), 2);
	dmat *L2f=NULL;
	dspfull(&L2f, L2, 'n', 1);
	assert_float_equal(dmaxabs(L2f), 0, 0);
	dfree(L2f); dspfree(L2); locfree(loc2);
}
/** Square actuator grid of n*n with unit spacing centered at origin. */
static loccell *mkaloc(long n){
	loccell *aloc=loccellnew(1, 1);
	P(aloc, 0)=mksqloc(n, n, 1, 1, -(n-1)*0.5, -(n-1)*0.5);
	return aloc;
}
/** Coupling: border actuators get val, interior get 1. */
static dcell *mkcpl_border(long n, real val){
	dcell *c=dcellnew(1, 1);
	P(c, 0)=dnew(n*n, 1);
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			P(P(c, 0), ix+iy*n)=(ix==0||iy==0||ix==n-1||iy==n-1)?val:1;
		}
	}
	return c;
}
static dmat *spmulvec(const dsp *A, const dmat *x){
	dmat *y=NULL;
	dspmm(&y, A, x, "nn", 1);
	return y;
}
static void assert_sp_symmetric(const dsp *A, real tol){
	dsp *At=dsptrans(A);
	dsp *diff=dspadd2(A, 1, At, -1);
	dmat *fd=NULL;
	dspfull(&fd, diff, 'n', 1);
	assert_float_equal(dmaxabs(fd), 0, tol);
	dfree(fd); dspfree(diff); dspfree(At);
}
/**
   genactcpl(): coupling is normalized column sum of |HA| (or HA'*W1), and an
   actuator fully coupled in any direction is clipped to 1.
*/
static void recon_genactcpl(void **state){
	(void)state;
	dmat *H=dnew(3, 4);
	P(H, 0, 0)=1; P(H, 2, 0)=1;
	P(H, 0, 2)=-2; P(H, 1, 2)=1;
	P(H, 2, 3)=0.5;
	dspcell *HA=dspcellnew(1, 1);
	P(HA, 0)=d2sp(H, 0);
	dcell *cpl=genactcpl(HA, NULL);
	const real ref1[4]={2./3, 0, 1, 1./6};
	for(int i=0; i<4; i++){
		assert_float_equal(P(P(cpl, 0), i), ref1[i], 1e-14);
	}
	dcellfree(cpl);
	/*with weighting: HA'*W1*/
	dmat *W1=dnew(3, 1); dset(W1, 1);
	cpl=genactcpl(HA, W1);
	const real ref2[4]={2./1, 0, -1./1, 0.5/1};/*column sum [2,0,-1,0.5], max 2*/
	for(int i=0; i<4; i++){
		assert_float_equal(P(P(cpl, 0), i), ref2[i]*0.5, 1e-14);
	}
	dcellfree(cpl);
	/*two fit directions: threshold is max/2; above it clipped to 1*/
	dspcell *HA2=dspcellnew(2, 1);
	P(HA2, 0)=d2sp(H, 0);
	P(HA2, 1)=d2sp(H, 0);
	cpl=genactcpl(HA2, NULL);
	const real ref3[4]={1, 0, 1, 1./3};
	for(int i=0; i<4; i++){
		assert_float_equal(P(P(cpl, 0), i), ref3[i], 1e-14);
	}
	dcellfree(cpl);
	dspcellfree(HA2); dspcellfree(HA); dfree(W1); dfree(H);
}
/**
   slaving(): the regularization S=H'H is symmetric, scaled by sclsq, has
   constants in its null space (mode 1, 3, 4), touches only inactive
   actuators and their neighbors, constrains stuck actuators individually, and
   in mode 2 averages each active island.
*/
static void recon_slaving(void **state){
	(void)state;
	const long n=7;
	const long nact=n*n;
	loccell *aloc=mkaloc(n);
	dmat *ones=dnew(nact, 1); dset(ones, 1);
	for(int mode=1; mode<=4; mode++){
		if(mode==2) continue;
		dcell *cpl=mkcpl_border(n, 0.05);
		dspcell *S=slaving(aloc, cpl, NULL, NULL, 0.5, 1, mode);
		assert_non_null(S);
		dsp *S0=P(S, 0, 0);
		assert_int_equal(NX(S0), nact);
		assert_sp_symmetric(S0, 1e-12);
		dmat *y=spmulvec(S0, ones);
		assert_float_equal(dmaxabs(y), 0, 1e-12);
		dfree(y);
		/*center actuator is far from inactive ones: not constrained*/
		dmat *Sf=NULL;
		dspfull(&Sf, S0, 'n', 1);
		const long ic=3+3*n;
		for(long i=0; i<nact; i++){
			assert_float_equal(P(Sf, ic, i), 0, 0);
		}
		if(mode==1){
			/*every border actuator is constrained*/
			for(long i=0; i<n; i++){
				assert_true(P(Sf, i, i)>0);
			}
			/*sclsq scaling*/
			dspcell *S2=slaving(aloc, cpl, NULL, NULL, 0.5, 2.5, mode);
			dmat *S2f=NULL;
			dspfull(&S2f, P(S2, 0, 0), 'n', 1);
			dscale(Sf, 2.5);
			assert_dmat_near(S2f, Sf, 1e-12);
			dfree(S2f);
			dspcellfree(S2);
		}
		dfree(Sf);
		dspcellfree(S);
		dcellfree(cpl);
	}
	/*mode 1 with a stuck interior actuator and a floating one*/
	{
		dcell *cpl=mkcpl_border(n, 0.05);
		lcell *stuck=lcellnew(1, 1);
		P(stuck, 0)=lnew(nact, 1);
		const long is=3+3*n;
		P(P(stuck, 0), is)=1;
		lcell *flt=lcellnew(1, 1);
		P(flt, 0)=lnew(nact, 1);
		const long ifl=2+3*n;
		P(P(flt, 0), ifl)=1;
		dspcell *S=slaving(aloc, cpl, stuck, flt, 0.5, 1, 1);
		/*actcpl is updated: stuck to 1, float to 0*/
		assert_float_equal(P(P(cpl, 0), is), 1, 0);
		assert_float_equal(P(P(cpl, 0), ifl), 0, 0);
		dmat *Sf=NULL;
		dspfull(&Sf, P(S, 0, 0), 'n', 1);
		/*stuck actuator only constrained by itself*/
		assert_float_equal(P(Sf, is, is), 1, 1e-14);
		for(long i=0; i<nact; i++){
			if(i!=is) assert_float_equal(P(Sf, is, i), 0, 0);
		}
		/*floating actuator is slaved*/
		assert_true(P(Sf, ifl, ifl)>0);
		dfree(Sf);
		dspcellfree(S);
		lcellfree(stuck); lcellfree(flt); dcellfree(cpl);
	}
	/*mode 2: two separate active islands (columns 0-1 and 5-6)*/
	{
		dcell *cpl=dcellnew(1, 1);
		P(cpl, 0)=dnew(nact, 1);
		for(long iy=0; iy<n; iy++){
			for(long ix=0; ix<n; ix++){
				P(P(cpl, 0), ix+iy*n)=(ix<2||ix>4)?1:0.01;
			}
		}
		dspcell *S=slaving(aloc, cpl, NULL, NULL, 0.5, 1, 2);
		dsp *S0=P(S, 0, 0);
		assert_sp_symmetric(S0, 1e-12);
		dmat *y=spmulvec(S0, ones);
		for(long iy=0; iy<n; iy++){
			for(long ix=0; ix<n; ix++){
				real ref=(ix<2||ix>4)?1./(2*n):0;
				assert_float_equal(P(y, ix+iy*n), ref, 1e-12);
			}
		}
		/*a differential piston between islands is not penalized more than common*/
		dmat *v=dnew(nact, 1);
		for(long iy=0; iy<n; iy++){
			for(long ix=0; ix<n; ix++){
				P(v, ix+iy*n)=ix<2?1:(ix>4?-1:0);
			}
		}
		assert_float_equal(dspwdot(v, S0, v), dspwdot(ones, S0, ones), 1e-12);
		dfree(v);
		dfree(y);
		dspcellfree(S);
		dcellfree(cpl);
	}
	/*missing actcplc is an error*/
	assert_null(slaving(aloc, NULL, NULL, NULL, 0.5, 1, 1));
	dfree(ones);
	cellfree(aloc);
}
/**
   act_stuck(), act_zero(), act_stuck_cmd(): zero the columns/rows/entries of
   stuck or dead actuators in dense, sparse, and vector forms; invalid input
   is rejected.
*/
static void recon_act_stuck(void **state){
	(void)state;
	const long n=3, nact=n*n, nfit=4;
	loccell *aloc=mkaloc(n);
	lcell *stuck=lcellnew(1, 1);
	P(stuck, 0)=lnew(nact, 1);
	P(P(stuck, 0), 4)=1;
	P(P(stuck, 0), 7)=1;
	dmat *H=mkrand(nfit, nact, 51);
	/*dense interaction matrix*/
	dcell *HA=dcellnew(1, 1);
	P(HA, 0)=ddup(H);
	act_stuck(aloc, HA, stuck);
	for(long ia=0; ia<nact; ia++){
		for(long i=0; i<nfit; i++){
			real ref=(ia==4||ia==7)?0:P(H, i, ia);
			assert_float_equal(P(P(HA, 0), i, ia), ref, 0);
		}
	}
	/*coupling vector*/
	dcell *cpl=dcellnew(1, 1);
	P(cpl, 0)=dnew(nact, 1); dset(P(cpl, 0), 1);
	act_stuck(aloc, cpl, stuck);
	assert_float_equal(dsum(P(cpl, 0)), nact-2, 0);
	assert_float_equal(P(P(cpl, 0), 4), 0, 0);
	/*sparse*/
	dspcell *HS=dspcellnew(1, 1);
	P(HS, 0)=d2sp(H, 0);
	act_stuck(aloc, HS, stuck);
	dmat *HSf=NULL;
	dspfull(&HSf, P(HS, 0), 'n', 1);
	assert_dmat_near(HSf, P(HA, 0), 0);
	dfree(HSf);
	/*NULL stuck: nothing done*/
	act_stuck(aloc, HS, NULL);
	/*invalid format and invalid size*/
	dcell *bad=dcellnew(2, 2);
	act_stuck(aloc, bad, stuck);
	dcellfree(bad);
	dcell *bad2=dcellnew(1, 1);
	P(bad2, 0)=dnew(2, 3);
	act_stuck(aloc, bad2, stuck);
	assert_float_equal(dsum(P(bad2, 0)), 0, 0);
	dcellfree(bad2);
	/*act_zero: rows of dead actuators in mode vectors*/
	dcell *HB=dcellnew(1, 2);
	P(HB, 0, 0)=mkrand(nact, 3, 52);
	P(HB, 0, 1)=mkrand(nact, 2, 53);
	dcell *HB0=dcelldup(HB);
	act_zero(aloc, HB, stuck);
	for(int im=0; im<2; im++){
		for(long iy=0; iy<NY(P(HB, 0, im)); iy++){
			for(long ia=0; ia<nact; ia++){
				real ref=(ia==4||ia==7)?0:P(P(HB0, 0, im), ia, iy);
				assert_float_equal(P(P(HB, 0, im), ia, iy), ref, 0);
			}
		}
	}
	act_zero(aloc, HB, NULL);
	/*act_stuck_cmd: zero commands*/
	dcell *adm=dcellnew(1, 1);
	P(adm, 0)=dnew(nact, 2); dset(P(adm, 0), 3);
	act_stuck_cmd(aloc, adm, stuck);
	assert_float_equal(dsum(P(adm, 0)), 3*2*(nact-2), 0);
	assert_float_equal(P(P(adm, 0), 7, 1), 0, 0);
	act_stuck_cmd(aloc, adm, NULL);
	dcellfree(adm); dcellfree(HB); dcellfree(HB0);
	dspcellfree(HS); dcellfree(cpl); dcellfree(HA); dfree(H);
	lcellfree(stuck);
	cellfree(aloc);
}
/**
   act_float(): the influence of a floating actuator is moved equally to its 4
   neighbors (preserving the total response to a piston), for sparse and dense
   matrices consistently. act_float_interp(): floating actuators take the
   average of their neighbors.
*/
static void recon_act_float(void **state){
	(void)state;
	const long n=5, nact=n*n, nfit=4;
	loccell *aloc=mkaloc(n);
	lcell *flt=lcellnew(1, 1);
	P(flt, 0)=lnew(nact, 1);
	const long ic=2+2*n;/*center: 4 neighbors*/
	const long io=0;    /*corner: 2 neighbors*/
	P(P(flt, 0), ic)=1;
	P(P(flt, 0), io)=1;
	dmat *H=mkrand(nfit, nact, 61);
	dmat *ones=dnew(nact, 1); dset(ones, 1);
	dmat *y0=NULL;
	dmm(&y0, 0, H, ones, "nn", 1);
	/*sparse*/
	dspcell *HA=dspcellnew(1, 1);
	P(HA, 0)=d2sp(H, 0);
	act_float(aloc, &HA, NULL, flt);
	dmat *HAf=NULL;
	dspfull(&HAf, P(HA, 0), 'n', 1);
	dmat *y1=NULL;
	dmm(&y1, 0, HAf, ones, "nn", 1);
	assert_dmat_near(y1, y0, 1e-12);
	for(long i=0; i<nfit; i++){
		assert_float_equal(P(HAf, i, ic), 0, 0);
		assert_float_equal(P(HAf, i, io), 0, 0);
		/*right neighbor of center gets 1/4 of the center column*/
		assert_float_equal(P(HAf, i, ic+1), P(H, i, ic+1)+0.25*P(H, i, ic), 1e-12);
		/*corner neighbors get 1/2 each*/
		assert_float_equal(P(HAf, i, io+1), P(H, i, io+1)+0.5*P(H, i, io), 1e-12);
		assert_float_equal(P(HAf, i, io+n), P(H, i, io+n)+0.5*P(H, i, io), 1e-12);
	}
	/*dense gives the same result*/
	dcell *HB=dcellnew(1, 1);
	P(HB, 0)=ddup(H);
	act_float(aloc, NULL, HB, flt);
	assert_dmat_near(P(HB, 0), HAf, 1e-12);
	/*nothing to do*/
	act_float(aloc, NULL, NULL, flt);
	act_float(aloc, &HA, NULL, NULL);
	/*interpolator*/
	dspcell *Hi=act_float_interp(aloc, flt);
	dmat *lin=dnew(nact, 1);
	for(long i=0; i<nact; i++){
		P(lin, i)=2*P(aloc, 0)->locx[i]+5*P(aloc, 0)->locy[i]+1;
	}
	dmat *y=spmulvec(P(Hi, 0, 0), lin);
	for(long i=0; i<nact; i++){
		if(i==io){
			assert_float_equal(P(y, i), 0.5*(P(lin, 1)+P(lin, n)), 1e-12);
		} else{/*exact for linear function at interior and identity elsewhere*/
			assert_float_equal(P(y, i), P(lin, i), 1e-12);
		}
	}
	dfree(y);
	dspcellfree(Hi);
	/*without floating actuators it is identity*/
	Hi=act_float_interp(aloc, NULL);
	y=spmulvec(P(Hi, 0, 0), lin);
	assert_dmat_near(y, lin, 0);
	dfree(y); dspcellfree(Hi);
	dfree(lin); dfree(HAf); dfree(y1); dfree(y0); dfree(ones); dfree(H);
	dcellfree(HB); dspcellfree(HA); lcellfree(flt); cellfree(aloc);
}
/**
   act_extrap(): active actuators are unchanged, inactive ones are extrapolated
   from active neighbors (iteratively for corners); constants are preserved.
   With low order removal (lor) piston/tip/tilt are reproduced exactly.
*/
static void recon_act_extrap(void **state){
	(void)state;
	const long n=6, nact=n*n;
	loccell *aloc=mkaloc(n);
	dcell *cpl=mkcpl_border(n, 0);
	rand_t rstat; seed_rand(&rstat, 71);
	dmat *x=dnew(nact, 1); drandn(x, 1, &rstat);
	dmat *ones=dnew(nact, 1); dset(ones, 1);
	for(int lor=0; lor<2; lor++){
		dspcell *E=act_extrap(aloc, cpl, 0.5, lor);
		dsp *E0=P(E, 0, 0);
		dmat *y=spmulvec(E0, ones);
		if(!lor){/*constant is preserved*/
			assert_dmat_near(y, ones, 1e-10);
		}
		dfree(y);
		y=spmulvec(E0, x);
		for(long iy=0; iy<n; iy++){
			for(long ix=0; ix<n; ix++){
				long i=ix+iy*n;
				int border=(ix==0||iy==0||ix==n-1||iy==n-1);
				if(!border){
					assert_float_equal(P(y, i), P(x, i), 1e-10);
				} else if(!lor){
					int corner=(ix==0||ix==n-1)&&(iy==0||iy==n-1);
					if(!corner){
						/*only one active neighbor (the inner one)*/
						long jx=ix==0?1:(ix==n-1?n-2:ix);
						long jy=iy==0?1:(iy==n-1?n-2:iy);
						assert_float_equal(P(y, i), P(x, jx+jy*n), 1e-10);
					}
				}
			}
		}
		dfree(y);
		if(lor){
			/*the low order (zernike radial order 1 to 3) modes are reproduced exactly*/
			dmat *Mz=zernike(P(aloc, 0), loc_diam(P(aloc, 0)), 1, 3, 0);
			assert_int_equal(NY(Mz), 9);
			dmat *EM=NULL;
			dspmm(&EM, E0, Mz, "nn", 1);
			assert_dmat_near(EM, Mz, 1e-9);
			dfree(EM); dfree(Mz);
		}
		dspcellfree(E);
	}
	dfree(x); dfree(ones); dcellfree(cpl); cellfree(aloc);
}
UNITTEST_MAIN(
	cmocka_unit_test(recon_muv_apply),
	cmocka_unit_test(recon_muv_fun),
	cmocka_unit_test(recon_pcg),
	cmocka_unit_test(recon_muv_direct),
	cmocka_unit_test(recon_muv_bgs),
	cmocka_unit_test(recon_laplacian_map),
	cmocka_unit_test(recon_laplacian_loc),
	cmocka_unit_test(recon_genactcpl),
	cmocka_unit_test(recon_slaving),
	cmocka_unit_test(recon_act_stuck),
	cmocka_unit_test(recon_act_float),
	cmocka_unit_test(recon_act_extrap)
)
