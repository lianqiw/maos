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
  Unit tests for assorted matrix operations and the control related routines in
  lib/. Ported from test_mat.c, where main() exited in the first call
  (test_servo) so that none of the other tests ever ran, and most of them read
  or wrote files. All inputs are now synthesized in memory.
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-10
#endif
/*Assert that A is (close to) the identity matrix.*/
static void assert_identity(const dmat *A, real tol){
	assert_non_null(A);
	assert_int_equal(NX(A), NY(A));
	for(long iy=0; iy<NY(A); iy++){
		for(long ix=0; ix<NX(A); ix++){
			assert_float_equal(P(A, ix, iy), ix==iy?1:0, tol);
		}
	}
}

/** dclip limits values to [min,max], returns the number clipped, and handles NULL and infinite bounds (test_clip). */
static void matops_clip(void **state){
	(void)state;
	const int N=64;
	dmat *A=dnew(N, N);
	rand_t rstat;
	seed_rand(&rstat, 1);
	drandn(A, 1, &rstat);
	dmat *A0=ddup(A);
	long nexpect=0;
	for(long i=0; i<PN(A); i++){
		if(P(A, i)>0.5||P(A, i)<-0.5) nexpect++;
	}
	assert_true(nexpect>0);
	assert_int_equal(dclip(A, -0.5, 0.5), nexpect);
	for(long i=0; i<PN(A); i++){
		real expect=P(A0, i)>0.5?0.5:(P(A0, i)<-0.5?-0.5:P(A0, i));
		assert_float_equal(P(A, i), expect, 0);
	}
	//Clipping again does nothing.
	assert_int_equal(dclip(A, -0.5, 0.5), 0);
	//Infinite bounds and NULL are no-ops.
	assert_int_equal(dclip(A0, -INFINITY, INFINITY), 0);
	assert_int_equal(dclip(NULL, -1, 1), 0);
	//One sided limit.
	dcp(&A, A0);
	assert_true(dclip(A, -INFINITY, 0)>0);
	assert_float_equal(dmax(A), 0, 0);
	dfree(A);
	dfree(A0);
}

/** dhistfill accumulates one histogram column per element, with out of range values in the end bins (test_hist). */
static void matops_histfill(void **state){
	(void)state;
	const int N=8, nbin=100;
	dmat *A=dnew(N, N);
	dmat *B=NULL;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const int nrep=100;
	for(int i=0; i<nrep; i++){
		drandn(A, 1, &rstat);
		dhistfill(&B, A, 0, 0.1, nbin);
	}
	assert_int_equal(NX(B), nbin);
	assert_int_equal(NY(B), N*N);
	for(long i=0; i<NY(B); i++){
		real sum=0;
		for(long ib=0; ib<nbin; ib++) sum+=P(B, ib, i);
		assert_float_equal(sum, nrep, 0);
	}
	//Known bin placement
	dmat *v=dnew(4, 1);
	P(v, 0)=0; P(v, 1)=0.21; P(v, 2)=-100; P(v, 3)=100;
	dmat *h=NULL;
	dhistfill(&h, v, 0, 0.1, 10);
	assert_float_equal(P(h, 5, 0), 1, 0);//center maps to n/2
	assert_float_equal(P(h, 7, 1), 1, 0);//0.21/0.1 rounds to 2
	assert_float_equal(P(h, 0, 2), 1, 0);//below range goes to first bin
	assert_float_equal(P(h, 9, 3), 1, 0);//above range goes to last bin
	assert_float_equal(dsum(h), 4, 0);
	//NULL input is a no-op
	dhistfill(&h, NULL, 0, 0.1, 10);
	assert_float_equal(dsum(h), 4, 0);
	dfree(A); dfree(B); dfree(v); dfree(h);
}

/** dcellpinv of a block diagonal TT matrix with sparse weighting equals the blockwise dpinv and the dense formula (test_dcellpinv with synthesized TT/saneai). */
static void matops_cellpinv(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const int nwfs=3, nsa=10;
	dcell *TT=dcellnew(nwfs, nwfs);
	dspcell *saneai=dspcellnew(nwfs, nwfs);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		dmat *tt=P(TT, iwfs, iwfs)=dnew(nsa*2, 2);
		drandn(tt, 0.1, &rstat);
		for(int isa=0; isa<nsa; isa++){
			P(tt, isa, 0)+=1;
			P(tt, isa+nsa, 1)+=1;
		}
		dmat *wt=dnew(nsa*2, 1);
		for(long i=0; i<PN(wt); i++) P(wt, i)=1+randu(&rstat);
		P(saneai, iwfs, iwfs)=dspnewdiag(nsa*2, P(wt), 1);
		dfree(wt);
	}
	dcell *PTT=dcellpinv(TT, saneai);
	assert_non_null(PTT);
	assert_int_equal(NX(PTT), nwfs);
	assert_int_equal(NY(PTT), nwfs);
	for(int iy=0; iy<nwfs; iy++){
		for(int ix=0; ix<nwfs; ix++){
			if(ix!=iy){
				assert_true(!P(PTT, ix, iy)||dmaxabs(P(PTT, ix, iy))==0);
				continue;
			}
			dmat *ptt=dpinv(P(TT, ix, ix), P(saneai, ix, ix));
			assert_dmat_near(P(PTT, ix, ix), ptt, TOL);
			//left inverse: PTT*TT=I
			dmat *I2=NULL;
			dmm(&I2, 0, ptt, P(TT, ix, ix), "nn", 1);
			assert_identity(I2, TOL);
			//dense formula inv(A'WA)A'W
			dmat *W=NULL;
			dspfull(&W, P(saneai, ix, ix), 'n', 1);
			dmat *AtW=NULL, *AtWA=NULL, *ref=NULL;
			dmm(&AtW, 0, P(TT, ix, ix), W, "tn", 1);
			dmm(&AtWA, 0, AtW, P(TT, ix, ix), "nn", 1);
			dmat *inv=dinv(AtWA);
			dmm(&ref, 0, inv, AtW, "nn", 1);
			assert_dmat_near(ptt, ref, TOL);
			dfree(ptt); dfree(I2); dfree(W); dfree(AtW); dfree(AtWA); dfree(inv); dfree(ref);
		}
	}
	assert_null(dcellpinv(NULL, saneai));
	assert_null(dpinv(NULL, NULL));
	dcellfree(TT);
	dspcellfree(saneai);
	dcellfree(PTT);
}

/** dcellreduce merges along rows (dim=1) or columns (dim=2) consistently with dcell2m, keeping empty blocks zero (test_dcellcat). */
static void matops_cellreduce(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long nxs[2]={2, 3}, nys[3]={1, 2, 4};
	dcell *A=dcellnew(2, 3);
	for(int iy=0; iy<3; iy++){
		for(int ix=0; ix<2; ix++){
			if(ix==0&&iy==1) continue;//leave one block empty
			P(A, ix, iy)=dnew(nxs[ix], nys[iy]);
			drandn(P(A, ix, iy), 1, &rstat);
		}
	}
	dmat *full=dcell2m(A);
	assert_int_equal(NX(full), 5);
	assert_int_equal(NY(full), 7);
	dcell *B=dcellreduce(A, 1);
	assert_int_equal(NX(B), 1);
	assert_int_equal(NY(B), 3);
	long col0=0;
	for(int iy=0; iy<3; iy++){
		assert_int_equal(NX(P(B, iy)), 5);
		assert_int_equal(NY(P(B, iy)), nys[iy]);
		for(long jy=0; jy<nys[iy]; jy++){
			for(long jx=0; jx<5; jx++){
				assert_float_equal(P(P(B, iy), jx, jy), P(full, jx, col0+jy), 0);
			}
		}
		col0+=nys[iy];
	}
	dcell *C=dcellreduce(A, 2);
	assert_int_equal(NX(C), 2);
	assert_int_equal(NY(C), 1);
	long row0=0;
	for(int ix=0; ix<2; ix++){
		assert_int_equal(NX(P(C, ix)), nxs[ix]);
		assert_int_equal(NY(P(C, ix)), 7);
		for(long jy=0; jy<7; jy++){
			for(long jx=0; jx<nxs[ix]; jx++){
				assert_float_equal(P(P(C, ix), jx, jy), P(full, row0+jx, jy), 0);
			}
		}
		row0+=nxs[ix];
	}
	//Merging both ways gives the full matrix.
	dcell *D=dcellreduce(B, 2);
	assert_dmat_near(P(D, 0), full, 0);
	//Invalid input.
	assert_null(dcellreduce(A, 3));
	assert_null(dcellreduce(NULL, 1));
	dcellfree(A); dcellfree(B); dcellfree(C); dcellfree(D); dfree(full);
}

/** dcellnew3/dcellnew_file/dnew_file allocate zeroed, contiguous blocks of the requested shapes without files when the name starts with '-' (test_save, test_async). */
static void matops_cellnew(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	long nnx[6]={2, 3, 4, 5, 6, 0};
	long nny[6]={3, 4, 2, 5, 6, 3};
	dcell *b=dcellnew3(2, 3, nnx, nny);
	assert_int_equal(NX(b), 2);
	assert_int_equal(NY(b), 3);
	assert_non_null(b->m);
	long tot=0;
	for(int i=0; i<6; i++){
		assert_int_equal(PN(P(b, i)), nnx[i]*nny[i]);
		if(nnx[i]>0&&nny[i]>0){
			assert_int_equal(NX(P(b, i)), nnx[i]);
			assert_int_equal(NY(P(b, i)), nny[i]);
			assert_ptr_equal(P(P(b, i)), P(b->m)+tot);//contiguous
		}
		tot+=nnx[i]*nny[i];
	}
	assert_int_equal(PN(b->m), tot);
	assert_float_equal(dsumsq(b->m), 0, 0);
	//writing into blocks writes into the backing vector
	for(int i=0; i<5; i++) drandn(P(b, i), 1, &rstat);
	real ss=0;
	for(int i=0; i<5; i++) ss+=dsumsq(P(b, i));
	assert_float_equal(dsumsq(b->m), ss, 1e-12*ss);
	dcellfree(b);
	//same blocks
	dcell *s=dcellnew_same(2, 2, 3, 4);
	for(int i=0; i<4; i++){
		assert_int_equal(NX(P(s, i)), 3);
		assert_int_equal(NY(P(s, i)), 4);
	}
	dcellfree(s);
	//File backed creation with saving disabled by a leading '-'.
	dmat *xx=dnew_file(10, 10, "test_mat", "-test_mat.bin");
	assert_non_null(xx);
	assert_null(xx->fp);
	assert_string_equal(xx->keywords, "test_mat");
	long sizes[3]={10, 0, 10};
	dcell *xc=dcellnew_file(3, 1, sizes, sizes, "test0", "-test_cell.bin");
	assert_non_null(xc);
	assert_null(xc->fp);
	assert_string_equal(xc->keywords, "test0");
	for(int i=0; i<PN(xc); i++){
		dmat *x=P(xc, i);
		for(long iy=0; iy<NY(x); iy++){
			for(long ix=0; ix<NX(x); ix++){
				P(x, ix, iy)=ix+iy*10;
			}
		}
	}
	dcp(&xx, P(xc, 0));
	for(long i=0; i<100; i++) assert_float_equal(P(xx, i), i, 0);
	assert_float_equal(dsum(P(xc, 2)), 4950, 0);
	assert_int_equal(PN(P(xc, 1)), 0);
	dfree(xx);
	dcellfree(xc);
	//Zero size returns NULL
	assert_null(dnew_file(0, 3, NULL, "-x"));
	assert_null(dcellnew_file(0, 1, sizes, sizes, NULL, "-x"));
}

/** 1-d cubic spline: exact at nodes, reproduces linear functions everywhere and quadratics in the interior; all interfaces agree (test_spline). */
static void matops_spline(void **state){
	(void)state;
	const long nx=20;
	dmat *x=dnew(nx, 1);
	dmat *y=dnew(nx, 1);
	dmat *x2=dnew(nx*100, 1);
	for(int i=0; i<nx; i++){
		P(x, i)=i;
		P(y, i)=sin(i*M_PI/2.);
	}
	for(int i=0; i<nx*100; i++){
		P(x2, i)=0.01*i*(nx-1)/nx;//stay within [0, nx-1]
	}
	dmat *coeff=dspline_prep(x, y);
	assert_int_equal(NX(coeff), 4);
	dmat *y2=dspline_eval(coeff, x, x2);
	dmat *y3=dspline(x, y, x2);
	assert_dmat_near(y2, y3, 0);
	//nodes
	dmat *yn=dspline_eval(coeff, x, x);
	assert_dmat_near(yn, y, TOL);
	//two column input form
	dmat *xy=dnew(nx, 2);
	for(int i=0; i<nx; i++){
		P(xy, i, 0)=P(x, i);
		P(xy, i, 1)=P(y, i);
	}
	dmat *y4=dspline(xy, NULL, x2);
	assert_dmat_near(y4, y2, 0);
	//smooth function is well approximated in the interior
	for(int i=0; i<nx; i++){
		P(x, i)=i*0.5+1;
		P(y, i)=sin(P(x, i)*0.6);
	}
	dmat *y5=dspline(x, y, x2);//x2 is partially outside: just check it is finite
	for(long i=0; i<PN(y5); i++) assert_false(isnan(P(y5, i)));
	dmat *xs=dlinspace(2, 0.013, 500);
	dmat *ys=dspline(x, y, xs);
	for(long i=0; i<PN(xs); i++){
		assert_float_equal(P(ys, i), sin(P(xs, i)*0.6), 2e-3);
	}
	//linear function: exact everywhere, quadratic: exact in interior intervals
	dmat *yl=dnew(nx, 1), *yq=dnew(nx, 1);
	for(int i=0; i<nx; i++){
		P(yl, i)=3*P(x, i)-2;
		P(yq, i)=P(x, i)*P(x, i)-P(x, i);
	}
	dmat *xe=dlinspace(1, 0.01, (long)((nx-1)*0.5/0.01));//[1, 10.5)
	dmat *zl=dspline(x, yl, xe);
	dmat *zq=dspline(x, yq, xe);
	for(long i=0; i<PN(xe); i++){
		real xv=P(xe, i);
		assert_float_equal(P(zl, i), 3*xv-2, 1e-10);
		if(xv>=1.5&&xv<=10.0){
			assert_float_equal(P(zq, i), xv*xv-xv, 1e-10);
		}
	}
	dfree(x); dfree(y); dfree(x2); dfree(coeff); dfree(y2); dfree(y3); dfree(yn); dfree(xy); dfree(y4);
	dfree(y5); dfree(xs); dfree(ys); dfree(yl); dfree(yq); dfree(xe); dfree(zl); dfree(zq);
}

/** 2-d bicubic spline: exact at grid nodes, reproduces bilinear functions, approximates a smooth function (test_spline_2d). */
static void matops_bspline(void **state){
	(void)state;
	const long nx=20, ny=30;
	dmat *x=dnew(nx, 1);
	dmat *y=dnew(ny, 1);
	dmat *z=dnew(nx, ny);
	dmat *zb=dnew(nx, ny);
	for(int i=0; i<ny; i++) P(y, i)=i*1.2+10;
	for(int i=0; i<nx; i++) P(x, i)=i*2+10;
	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			P(z, ix, iy)=sin(ix*M_PI/14.2)*sin(iy*M_PI/12.3);
			P(zb, ix, iy)=1+0.5*P(x, ix)-0.3*P(y, iy)+0.01*P(x, ix)*P(y, iy);
		}
	}
	dcell *coeff=dbspline_prep(x, y, z);
	dcell *coeffb=dbspline_prep(x, y, zb);
	assert_int_equal(NX(coeff), nx);
	assert_int_equal(NY(coeff), ny);
	//grid nodes
	dmat *xg=dnew(nx, ny), *yg=dnew(nx, ny);
	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			P(xg, ix, iy)=P(x, ix);
			P(yg, ix, iy)=P(y, iy);
		}
	}
	dmat *zg=dbspline_eval(coeff, x, y, xg, yg);
	assert_dmat_near(zg, z, 1e-10);
	//Fine grid within the domain (old test used 11x oversampling with extrapolation).
	const long nxnew=nx*5, nynew=ny*5;
	dmat *xnew=dnew(nxnew, nynew);
	dmat *ynew=dnew(nxnew, nynew);
	for(int iy=0; iy<nynew; iy++){
		for(int ix=0; ix<nxnew; ix++){
			P(xnew, ix, iy)=10+ix*(2.*(nx-1))/nxnew;
			P(ynew, ix, iy)=10+iy*(1.2*(ny-1))/nynew;
		}
	}
	dmat *z22=dbspline_eval(coeff, x, y, xnew, ynew);
	dmat *zbb=dbspline_eval(coeffb, x, y, xnew, ynew);
	for(long i=0; i<PN(xnew); i++){
		real xv=P(xnew, i), yv=P(ynew, i);
		real fx=(xv-10)/2, fy=(yv-10)/1.2;//grid index units
		assert_float_equal(P(z22, i), sin(fx*M_PI/14.2)*sin(fy*M_PI/12.3), 1e-2);
		assert_float_equal(P(zbb, i), 1+0.5*xv-0.3*yv+0.01*xv*yv, 1e-9);
	}
	dfree(x); dfree(y); dfree(z); dfree(zb); dfree(xg); dfree(yg); dfree(zg);
	dfree(xnew); dfree(ynew); dfree(z22); dfree(zbb);
	dcellfree(coeff); dcellfree(coeffb);
}

/** dsvd decomposes A=U*S*VT with orthogonal U, V and sorted S; dsvd_pow computes inverse, square root and pseudo inverse (test_svd with a synthesized sparse input). */
static void matops_svd(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long n=8;
	dmat *A=dnew(n, n);
	drandn(A, 1, &rstat);
	//svd reconstruction
	dmat *U=NULL, *S=NULL, *VT=NULL;
	dsvd(&U, &S, &VT, A);
	assert_int_equal(PN(S), n);
	for(long i=1; i<n; i++) assert_true(P(S, i)<=P(S, i-1));
	dmat *US=ddup(U);
	dmuldiag(US, S);
	dmat *A2=NULL;
	dmm(&A2, 0, US, VT, "nn", 1);
	assert_dmat_near(A2, A, TOL);
	dmat *UtU=NULL, *VVt=NULL;
	dmm(&UtU, 0, U, U, "tn", 1);
	dmm(&VVt, 0, VT, VT, "nt", 1);
	assert_identity(UtU, TOL);
	assert_identity(VVt, TOL);
	//inverse of a general matrix
	dmat *Ai=ddup(A);
	dsvd_pow(Ai, -1);
	dmat *I1=NULL;
	dmm(&I1, 0, A, Ai, "nn", 1);
	assert_identity(I1, 1e-8);
	//inverting twice (the old test called dsvd_pow twice on a matrix read from "SVD")
	dsvd_pow(Ai, -1);
	assert_dmat_near(Ai, A, 1e-8);
	//sparse input converted to dense, as in the old test.
	dsp *As=d2sp(A, 0);
	dmat *Af=NULL;
	dspfull(&Af, As, 'n', 1);
	assert_dmat_near(Af, A, 0);
	//square root of an SPD matrix
	dmat *SPD=NULL;
	dmm(&SPD, 0, A, A, "nt", 1);
	daddI(SPD, 1);
	dmat *R=ddup(SPD);
	dsvd_pow(R, 0.5);
	dmat *R2=NULL;
	dmm(&R2, 0, R, R, "nn", 1);
	assert_dmat_near(R2, SPD, 1e-8);
	//pseudo inverse of a rank deficient symmetric matrix: truncated eigenvalues
	dmat *V=dnew(n, 2);
	drandn(V, 1, &rstat);
	dmat *L=NULL;
	dmm(&L, 0, V, V, "nt", 1);//rank 2
	dmat *Li=ddup(L);
	dsvd_pow(Li, -1);
	dmat *T1=NULL, *T2=NULL;
	dmm3(&T1, 0, L, Li, L, "nnn", 1);
	dmm3(&T2, 0, Li, L, Li, "nnn", 1);
	assert_dmat_near(T1, L, 1e-8);
	assert_dmat_near(T2, Li, 1e-8);
	//empty input is a no-op
	dmat *E=dnew(0, 0);
	dsvd_pow(E, -1);
	assert_int_equal(PN(E), 0);
	dfree(A); dfree(U); dfree(S); dfree(VT); dfree(US); dfree(A2); dfree(UtU); dfree(VVt); dfree(Ai); dfree(I1);
	dspfree(As); dfree(Af); dfree(SPD); dfree(R); dfree(R2); dfree(V); dfree(L); dfree(Li); dfree(T1); dfree(T2); dfree(E);
}

/** csvd_pow inverts a complex matrix and csvd reconstructs it (test_svd2 with a synthesized input). */
static void matops_csvd(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	const long n=6;
	cmat *A=cnew(n, n);
	crandn(A, 1, &rstat);
	cmat *Ai=cdup(A);
	csvd_pow(Ai, -1);
	cmat *I1=NULL;
	cmm(&I1, 0, A, Ai, "nn", 1);
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			assert_float_equal(cabs(P(I1, ix, iy)-(ix==iy?1.:0.)), 0, 1e-8);
		}
	}
	cmat *U=NULL, *VT=NULL;
	dmat *S=NULL;
	csvd(&U, &S, &VT, A);
	cmat *US=cdup(U);
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			P(US, ix, iy)*=P(S, iy);
		}
	}
	cmat *A2=NULL;
	cmm(&A2, 0, US, VT, "nn", 1);
	for(long i=0; i<PN(A); i++){
		assert_float_equal(cabs(P(A2, i)-P(A, i)), 0, TOL);
	}
	cfree(A); cfree(Ai); cfree(I1); cfree(U); cfree(VT); dfree(S); cfree(US); cfree(A2);
}

/** psd1d satisfies Parseval (sum of PSD equals mean square), locates a sinusoid, and psd1dt adds the frequency axis (test_psd1d). */
static void matops_psd1d(void **state){
	(void)state;
	const long n=100;
	dmat *tmp=dnew(n, 1);
	for(int i=0; i<n; i++){
		P(tmp, i)=sin(i/10.);
	}
	dmat *psd=psd1d(tmp, 1);
	assert_int_equal(NX(psd), n/2+1);
	assert_int_equal(NY(psd), 1);
	assert_float_equal(dsum(psd), dsumsq(tmp)/n, 1e-12);
	//pure tone at bin 5
	for(int i=0; i<n; i++){
		P(tmp, i)=sin(2*M_PI*5*i/n);
	}
	dfree(psd);
	psd=psd1d(tmp, 1);
	for(long i=0; i<NX(psd); i++){
		assert_float_equal(P(psd, i), i==5?0.5:0, 1e-12);
	}
	//overlapping segments: average of the mean squares of each segment
	dfree(psd);
	const long nseg=3;
	psd=psd1d(tmp, nseg);
	const long lseg2=n/(nseg+1);
	real ms=0;
	for(long iseg=0; iseg<nseg; iseg++){
		for(long i=0; i<2*lseg2; i++) ms+=pow(P(tmp, iseg*lseg2+i), 2);
	}
	ms/=2*lseg2*nseg;
	assert_int_equal(NX(psd), lseg2+1);
	assert_float_equal(dsum(psd), ms, 1e-12);
	//row vector input is the same as column vector
	dmat *row=dtrans(tmp);
	dmat *psd2=psd1d(row, nseg);
	assert_dmat_near(psd2, psd, 0);
	//with time step: frequency column and density
	const real dt=0.01;
	dmat *psdt=psd1dt(tmp, 1, dt);
	assert_int_equal(NY(psdt), 2);
	const real df=1./(n*dt);
	assert_float_equal(P(psdt, 5, 0), 5*df, 1e-12);
	assert_float_equal(P(psdt, 5, 1), 0.5/df, 1e-10);
	dfree(tmp); dfree(psd); dfree(row); dfree(psd2); dfree(psdt);
}

/** Second order harmonic oscillator filter: unit DC gain, low pass response, linearity, input shape handling (test_sho). */
static void matops_sho(void **state){
	(void)state;
	const real f0=200, zeta=0.9;
	const real dt=1./64000.;
	const long n=6400;//0.1 s
	dmat *x=dnew(n, 3);
	for(long i=0; i<n; i++){
		P(x, i, 0)=1;//step
		P(x, i, 1)=2;//scaled step
		P(x, i, 2)=sin(2*M_PI*2000*i*dt);//10x above resonance
	}
	dmat *y=sho_filter(x, dt, f0, zeta, 0);
	assert_int_equal(NX(y), n);
	assert_int_equal(NY(y), 3);
	assert_true(P(y, 0, 0)<0.01);
	assert_float_equal(P(y, n-1, 0), 1, 1e-4);
	assert_true(dmax(y)<=2*1.01);//overshoot is ~0.15% for zeta=0.9
	for(long i=0; i<n; i++) assert_float_equal(P(y, i, 1), 2*P(y, i, 0), 1e-12);
	real amp=0;
	for(long i=n/2; i<n; i++) amp=fmax(amp, fabs(P(y, i, 2)));
	assert_true(amp<0.02);//|H|~1/r^2=0.01
	//low frequency sinusoid at the slow sampling of the old test (input2 at 800 Hz) passes nearly unchanged
	const real dt2=1./800.;
	const long n2=800;
	dmat *x2=dnew(1, n2);//row vector
	for(long i=0; i<n2; i++) P(x2, 0, i)=sin(2*M_PI*5*i*dt2);
	dmat *y2=sho_filter(x2, dt2, f0, zeta, 0);
	assert_int_equal(NX(y2), n2);
	assert_int_equal(NY(y2), 1);
	assert_int_equal(NX(x2), 1);//input is not modified
	real amp2=0;
	for(long i=n2/2; i<n2; i++) amp2=fmax(amp2, fabs(P(y2, i)));
	assert_float_equal(amp2, 1, 0.02);
	//averaging mode also converges to the step value
	dmat *x3=dnew(100, 1);
	dset(x3, 1);
	dmat *y3=sho_filter(x3, dt2, f0, zeta, 1);
	assert_float_equal(P(y3, NX(y3)-1, 0), 1, 1e-4);
	//zero input gives zero output
	dmat *x0=dnew(10, 1);
	dmat *y0=sho_filter(x0, dt, f0, zeta, 0);
	assert_float_equal(dmaxabs(y0), 0, 0);
	dfree(x); dfree(y); dfree(x2); dfree(y2); dfree(x3); dfree(y3); dfree(x0); dfree(y0);
}

/** dcellmm handles cell*cell, mat*cell and mat*mat with cell or matrix output and matches dmm (test_mm). */
static void matops_cellmm(void **state){
	(void)state;
	dcell *A=dcellnew(1, 1);
	dcell *B=dcellnew(1, 1);
	P(A, 0)=dnew(3, 4);
	dset(P(A, 0), 1);
	P(B, 0)=dnew(4, 5);
	dset(P(B, 0), 1);
	dcell *C=NULL;
	dcellmm(&C, A, B, "nn", 1);
	assert_non_null(C);
	assert_int_equal(NX(C), 1);
	assert_int_equal(NY(C), 1);
	assert_int_equal(NX(P(C, 0)), 3);
	assert_int_equal(NY(P(C, 0)), 5);
	assert_float_equal(dmax(P(C, 0)), 4, 0);
	assert_float_equal(dmin(P(C, 0)), 4, 0);
	dcellmm(&C, P(A, 0), B, "nn", 1);//accumulates
	assert_float_equal(dmin(P(C, 0)), 8, 0);
	dcellmm(&C, P(A, 0), P(B, 0), "nn", 1);
	assert_float_equal(dmin(P(C, 0)), 12, 0);
	assert_float_equal(dmax(P(C, 0)), 12, 0);
	dmat *C2=NULL;
	dcellmm(&C2, P(A, 0), B, "nn", 1);
	assert_non_null(C2);
	assert_int_equal(NX(C2), 3);
	assert_int_equal(NY(C2), 5);
	assert_float_equal(dmin(C2), 4, 0);
	dcellmm(&C2, P(A, 0), P(B, 0), "nn", 2);
	assert_float_equal(dmin(C2), 12, 0);
	assert_float_equal(dmax(C2), 12, 0);
	//random blocks with transposes against dmm
	rand_t rstat;
	seed_rand(&rstat, 1);
	dcell *X=dcellnew(2, 2), *Y=dcellnew(2, 1);
	for(int i=0; i<4; i++){ P(X, i)=dnew(3, 3); drandn(P(X, i), 1, &rstat); }
	for(int i=0; i<2; i++){ P(Y, i)=dnew(3, 2); drandn(P(Y, i), 1, &rstat); }
	dcell *Z=NULL;
	dcellmm(&Z, X, Y, "tn", 0.5);
	dmat *Xm=dcell2m(X), *Ym=dcell2m(Y), *Zm=NULL;
	dmm(&Zm, 0, Xm, Ym, "tn", 0.5);
	dmat *Zc=dcell2m(Z);
	assert_dmat_near(Zc, Zm, TOL);
	//NULL input is a no-op
	dcellmm(&Z, NULL, Y, "nn", 1);
	dmat *Zc2=dcell2m(Z);
	assert_dmat_near(Zc2, Zm, TOL);
	dcellfree(A); dcellfree(B); dcellfree(C); dfree(C2);
	dcellfree(X); dcellfree(Y); dcellfree(Z); dfree(Xm); dfree(Ym); dfree(Zm); dfree(Zc); dfree(Zc2);
}

/** dexpm matches known matrix exponentials (zero, diagonal, nilpotent, rotation), satisfies exp(A)exp(-A)=I and accumulates with alpha (test_expm). */
static void matops_expm(void **state){
	(void)state;
	dmat *B=NULL;
	dmat *Z=dnew(3, 3);
	dexpm(&B, 0, Z, 1);
	assert_identity(B, TOL);
	//diagonal
	dmat *D=dnew(3, 3);
	P(D, 0, 0)=1; P(D, 1, 1)=-2; P(D, 2, 2)=0.5;
	dexpm(&B, 0, D, 1.5);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			assert_float_equal(P(B, i, j), i==j?exp(1.5*P(D, i, i)):0, 1e-9*exp(1.5));
		}
	}
	dfree(B);//the output must match the input size when not NULL.
	//nilpotent: exp([0 1;0 0])=[1 1;0 1]
	dmat *N=dnew(2, 2);
	P(N, 0, 1)=1;
	dexpm(&B, 0, N, 1);
	assert_float_equal(P(B, 0, 0), 1, TOL);
	assert_float_equal(P(B, 0, 1), 1, TOL);
	assert_float_equal(P(B, 1, 0), 0, TOL);
	assert_float_equal(P(B, 1, 1), 1, TOL);
	//rotation generator
	dmat *G=dnew(2, 2);
	P(G, 0, 1)=-1; P(G, 1, 0)=1;
	const real t=2.3;
	dexpm(&B, 0, G, t);
	assert_float_equal(P(B, 0, 0), cos(t), 1e-9);
	assert_float_equal(P(B, 0, 1), -sin(t), 1e-9);
	assert_float_equal(P(B, 1, 0), sin(t), 1e-9);
	assert_float_equal(P(B, 1, 1), cos(t), 1e-9);
	//random: exp(A)*exp(-A)=I
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *A=dnew(5, 5);
	drandn(A, 1, &rstat);
	dmat *E1=NULL, *E2=NULL, *EE=NULL;
	dexpm(&E1, 0, A, 1);
	dexpm(&E2, 0, A, -1);
	dmm(&EE, 0, E1, E2, "nn", 1);
	assert_identity(EE, 1e-9);
	//accumulation out=2*out+exp(A)
	dmat *E3=ddup(E1);
	dexpm(&E3, 2, A, 1);
	dscale(E1, 3);
	assert_dmat_near(E3, E1, 1e-9*dmaxabs(E1));
	dfree(B); dfree(Z); dfree(D); dfree(N); dfree(G); dfree(A); dfree(E1); dfree(E2); dfree(EE); dfree(E3);
}

/** reccati converges to the solution of the discrete algebraic Riccati equation (the disabled test_reccati with synthesized input). */
static void matops_reccati(void **state){
	(void)state;
	//scalar case a=0.9, q=r=c=1: P=a^2*P*r/(P+r)+q reduces to P^2-0.81P-1=0
	dmat *A=dnew(1, 1), *Q=dnew(1, 1), *C=dnew(1, 1), *Rn=dnew(1, 1);
	P(A, 0)=0.9; P(Q, 0)=1; P(C, 0)=1; P(Rn, 0)=1;
	dmat *K=NULL, *Pi=NULL;
	real diff=reccati(&K, &Pi, A, Q, C, Rn);
	assert_true(diff<1e-10);
	const real Pexp=(0.81+sqrt(0.81*0.81+4))/2;
	assert_float_equal(P(Pi, 0), Pexp, 1e-9);
	assert_float_equal(P(K, 0), Pexp/(Pexp+1), 1e-9);
	dfree(A); dfree(Q); dfree(C); dfree(Rn); dfree(K); dfree(Pi);
	//matrix case: check the DARE residual P=APA'-APC'(CPC'+R)^-1CPA'+Q
	A=dnew(2, 2);
	P(A, 0, 0)=0.9; P(A, 0, 1)=0.2; P(A, 1, 0)=-0.1; P(A, 1, 1)=0.7;
	Q=dnew(2, 2); P(Q, 0, 0)=1; P(Q, 1, 1)=0.5; P(Q, 0, 1)=P(Q, 1, 0)=0.1;
	C=dnew(1, 2); P(C, 0, 0)=1; P(C, 0, 1)=0.5;
	Rn=dnew(1, 1); P(Rn, 0)=0.3;
	diff=reccati(&K, &Pi, A, Q, C, Rn);
	assert_true(diff<1e-10);
	dmat *APAt=NULL, *APCt=NULL, *S=NULL, *res=NULL, *tmp=NULL;
	dmm3(&APAt, 0, A, Pi, A, "nnt", 1);
	dmm3(&APCt, 0, A, Pi, C, "nnt", 1);
	dmm3(&S, 0, C, Pi, C, "nnt", 1);
	dadd(&S, 1, Rn, 1);
	dinv_inplace(S);
	dcp(&res, APAt);
	dmm(&tmp, 0, APCt, S, "nn", 1);
	dmm(&res, 1, tmp, APCt, "nt", -1);
	dadd(&res, 1, Q, 1);
	assert_dmat_near(res, Pi, 1e-9);
	//Pi is symmetric positive definite
	assert_float_equal(P(Pi, 0, 1), P(Pi, 1, 0), 1e-9);
	assert_true(P(Pi, 0, 0)>0&&P(Pi, 0, 0)*P(Pi, 1, 1)-P(Pi, 0, 1)*P(Pi, 1, 0)>0);
	//Kinf=P*C'*(CPC'+R)^-1
	dmat *Kref=NULL;
	dmm3(&Kref, 0, Pi, C, S, "ntn", 1);
	assert_dmat_near(K, Kref, 1e-9);
	dfree(A); dfree(Q); dfree(C); dfree(Rn); dfree(K); dfree(Pi);
	dfree(APAt); dfree(APCt); dfree(S); dfree(res); dfree(tmp); dfree(Kref);
}

/*SDE coefficients [c1=2*zeta*omega0, c2=omega0^2, sigma] of a damped oscillator.*/
static dmat *sde_coeff(real f0, real zeta, real sigma){
	dmat *coeff=dnew(3, 1);
	const real omega0=2*M_PI*f0;
	P(coeff, 0)=2*zeta*omega0;
	P(coeff, 1)=omega0*omega0;
	P(coeff, 2)=sigma;
	return coeff;
}
/*Two column [f, psd] of the SDE.*/
static dmat *sde_psd_table(const dmat *coeff, real df, real fmax){
	long nf=(long)round(fmax/df);
	dmat *f=dlinspace(df, df, nf);
	dmat *psd=sde_psd2(f, coeff);
	dmat *out=dnew(nf, 2);
	for(long i=0; i<nf; i++){
		P(out, i, 0)=P(f, i);
		P(out, i, 1)=P(psd, i);
	}
	dfree(f); dfree(psd);
	return out;
}

/** sde_fit recovers the coefficients of a PSD generated by sde_psd2 from a perturbed initial guess (test_sde with synthesized PSD). */
static void matops_sde_fit(void **state){
	(void)state;
	dmat *coeff0=sde_coeff(10, 0.3, 1);
	dmat *psd=sde_psd_table(coeff0, 0.01, 50);
	dmat *coeff=ddup(coeff0);
	P(coeff, 0)*=1.5;
	P(coeff, 1)*=0.8;
	real diff=sde_fit(&coeff, psd, 0, 0);
	assert_true(diff<0.01);
	assert_int_equal(NX(coeff), 3);
	for(int i=0; i<3; i++){
		assert_float_equal(P(coeff, i)/P(coeff0, i), 1, 0.02);
	}
	dfree(coeff0); dfree(psd); dfree(coeff);
}

/** sde_kalman builds a consistent filter (state propagation is expm(Ac*dt), covariance is SPD) and kalman_test reduces the residual of a synthesized time series (test_kalman). */
static void matops_kalman(void **state){
	(void)state;
	const real dthi=1./800;
	dmat *coeff=sde_coeff(5, 0.5, 1);
	lmat *dtrat_wfs=lnew(1, 1);
	P(dtrat_wfs, 0)=1;
	dcell *Gwfs=dcellnew(1, 1);
	P(Gwfs, 0)=dnew(2, 1);
	dset(P(Gwfs, 0), 1);//two identical measurements of the mode
	dcell *Cnn=dcellnew(1, 1);
	P(Cnn, 0)=dnew(2, 2);
	daddI(P(Cnn, 0), 1e-9);//measurement noise well below the signal
	kalman_t *k=sde_kalman(coeff, dthi, dtrat_wfs, NULL, Gwfs, Cnn, NULL);
	assert_non_null(k);
	assert_int_equal(PN(k->dtrats), 1);
	assert_int_equal(P(k->dtrats, 0), 1);
	//continuous state propagation [-c1 -c2; 1 0]
	dmat *Ac=dnew(2, 2);
	P(Ac, 0, 0)=-P(coeff, 0); P(Ac, 0, 1)=-P(coeff, 1); P(Ac, 1, 0)=1;
	dmat *AdM=NULL;
	dexpm(&AdM, 0, Ac, dthi);
	assert_dmat_near(k->AdM, AdM, 1e-10);
	assert_dmat_near(P(k->Ad, 0), AdM, 1e-10);
	assert_int_equal(NX(P(k->P, 0)), 1);
	assert_true(P(P(k->P, 0), 0)>0);
	assert_int_equal(NX(P(k->M, 0)), 2);//nstate
	assert_int_equal(NY(P(k->M, 0)), 2);//nmeas
	assert_int_equal(NX(P(k->Cd, 0)), 2);
	//least squares reconstructor of [1;1] with equal noise is [0.5 0.5]
	dmat *Rlsq=dcell2m(P(k->Rlsq, 0));
	assert_float_equal(P(Rlsq, 0), 0.5, 1e-10);
	assert_float_equal(P(Rlsq, 1), 0.5, 1e-10);
	//time series following the SDE PSD
	dmat *psd=sde_psd_table(coeff, 0.05, 400);
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *ts=psd2ts(psd, &rstat, dthi, 4096);
	const real rms_in=sqrt(dsumsq(ts)/PN(ts));
	assert_true(rms_in>0);
	real rms_res[2];
	for(int flag=0; flag<2; flag++){//0: LQG, 1: integrator
		dmat *res=kalman_test(k, NULL, ts, flag);
		assert_int_equal(PN(res), PN(ts));
		rms_res[flag]=sqrt(dsumsq(res)/PN(res));
		assert_true(rms_res[flag]<0.2*rms_in);
		dfree(res);
	}
	assert_true(rms_res[0]<rms_res[1]);//LQG beats the integrator on its own model
	kalman_free(k);
	//Two modes and two WFS at different rates (1 and 4): the fast WFS sees mode 0,
	//the slow WFS sees mode 1, which is directly controlled by the slow loop.
	dmat *coeff2=dnew(3, 2);
	for(int i=0; i<3; i++){
		P(coeff2, i, 0)=P(coeff, i);
		P(coeff2, i, 1)=P(coeff, i);
	}
	P(coeff2, 1, 1)*=0.25;//mode 1 is slower (2.5 Hz)
	lmat *dtrat2=lnew(2, 1);
	P(dtrat2, 0)=1; P(dtrat2, 1)=4;
	lmat *mdirect=lnew(2, 1);
	P(mdirect, 1)=1;
	dcell *Gwfs2=dcellnew(2, 1);
	P(Gwfs2, 0)=dnew(1, 2); P(P(Gwfs2, 0), 0, 0)=1;
	P(Gwfs2, 1)=dnew(1, 2); P(P(Gwfs2, 1), 0, 1)=1;
	dcell *Cnn2=dcellnew(2, 2);
	for(int i=0; i<4; i++){
		P(Cnn2, i)=dnew(1, 1);
		daddI(P(Cnn2, i), 1e-9);
	}
	k=sde_kalman(coeff2, dthi, dtrat2, mdirect, Gwfs2, Cnn2, NULL);
	assert_non_null(k);
	assert_int_equal(PN(k->dtrats), 2);
	assert_int_equal(P(k->dtrats, 1), 4);
	assert_int_equal(NX(P(k->Ad, 1)), 4);//nmod*order states
	assert_int_equal(NX(P(k->Cd, 1)), 2);//both WFS measure at the slow rate
	assert_int_equal(NX(P(k->Cd, 0)), 2);
	assert_int_equal(NX(P(k->P, 1)), 2);
	//kalman_test does not modify its input.
	assert_int_equal(NX(ts), 4096);
	assert_int_equal(NY(ts), 1);
	dmat *ts2=dnew(2, PN(ts));
	for(long i=0; i<PN(ts); i++){
		P(ts2, 0, i)=P(ts, i);
		P(ts2, 1, i)=P(ts, PN(ts)-1-i);
	}
	dmat *res=kalman_test(k, NULL, ts2, 0);
	real rms2=sqrt(dsumsq(res)/PN(res));
	assert_true(rms2<0.5*rms_in);
	dfree(res); dfree(ts2); dfree(coeff2);
	kalman_free(k);
	lfree(dtrat2); lfree(mdirect); dcellfree(Gwfs2); dcellfree(Cnn2);
	kalman_free(NULL);
	dfree(coeff); lfree(dtrat_wfs); dcellfree(Gwfs); dcellfree(Cnn); dfree(Ac); dfree(AdM); dfree(Rlsq);
	dfree(psd); dfree(ts);
}

/*Reference integrator gain for phase margin pm from the analytic open loop
  transfer function used by servo_calc_init: phase=-pi/2-w*T/2-w*dt*(1+al).*/
static real servo_gain_ref(real dt, long dtrat, real al, real pm){
	const real T=dt*dtrat;
	const real w=(M_PI/2-pm)/(T/2+dt*(1+al));
	const real x=w*T/2;
	real mag=1./(2*sin(x));
	if(dtrat>1) mag*=pow(sin(x)/x, 2);
	return 1./mag;
}
/** servo_optim_margin matches the analytic gain for a given phase margin, decreases with latency and margin, and is reduced by an SHO (test_servo). */
static void matops_servo(void **state){
	(void)state;
	const real dt=1/600.;
	real g=servo_optim_margin(dt, 1, 1, M_PI/4, 0, 0);
	assert_float_equal(g, 2*sin(M_PI/20), 1e-3*g);//closed form for dtrat=1, al=1, pm=pi/4
	assert_float_equal(g, servo_gain_ref(dt, 1, 1, M_PI/4), 1e-3*g);
	//default margin is pi/4
	assert_float_equal(servo_optim_margin(dt, 1, 1, 0, 0, 0), g, 0);
	real g2=servo_optim_margin(dt, 4, 0.5, M_PI/3, 0, 0);
	assert_float_equal(g2, servo_gain_ref(dt, 4, 0.5, M_PI/3), 2e-3*g2);
	//more latency or margin requires lower gain
	assert_true(servo_optim_margin(dt, 1, 2, M_PI/4, 0, 0)<g);
	assert_true(servo_optim_margin(dt, 1, 1, M_PI/3, 0, 0)<g);
	//a DM with a SHO response adds phase lag
	real gs=servo_optim_margin(dt, 1, 1, M_PI/4, 100, 0.5);
	assert_true(gs>0&&gs<g);
}

UNITTEST_MAIN(
	cmocka_unit_test(matops_clip),
	cmocka_unit_test(matops_histfill),
	cmocka_unit_test(matops_cellpinv),
	cmocka_unit_test(matops_cellreduce),
	cmocka_unit_test(matops_cellnew),
	cmocka_unit_test(matops_spline),
	cmocka_unit_test(matops_bspline),
	cmocka_unit_test(matops_svd),
	cmocka_unit_test(matops_csvd),
	cmocka_unit_test(matops_psd1d),
	cmocka_unit_test(matops_sho),
	cmocka_unit_test(matops_cellmm),
	cmocka_unit_test(matops_expm),
	cmocka_unit_test(matops_reccati),
	cmocka_unit_test(matops_sde_fit),
	cmocka_unit_test(matops_kalman),
	cmocka_unit_test(matops_servo)
)
