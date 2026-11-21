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
  Unit tests for the generic cell functions in math/cell.c and the cell
  variants of the type generic routines in math/mat.c, math/matmath.c,
  math/blas.c, math/spmm.c and math/sp.c. Block results are compared against
  the same operation on the assembled dense matrix (X(cell2m)).
*/
#include <unistd.h>
#include <sys/socket.h>
#include "unittest.h"

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
/*Create a nx*ny dcell with block sizes from mx and my, filled with random numbers.*/
static dcell *randcell(long nx, long ny, const long *mx, const long *my, rand_t *rstat){
	dcell *A=dcellnew(nx, ny);
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			P(A, ix, iy)=randmat(mx[ix], my[iy], rstat);
		}
	}
	return A;
}
/*Assert that cell A assembles to dense matrix B.*/
static void assert_cell_dense(const dcell *A, const dmat *B, real tol){
	dmat *Am=dcell2m(A);
	assert_dmat_near(Am, B, tol);
	dfree(Am);
}

/** cellnew, cellinit, cellinit2, cell_is_diag and celldim. */
static void cell_basic(void **state){
	(void)state;
	cell *A=cellnew(-1, 3);
	assert_non_null(A);
	assert_int_equal(NX(A), 0);
	assert_int_equal(NY(A), 3);
	assert_null(P(A));
	assert_true(iscell(A));
	cellfree(A);
	assert_null(A);
	assert_false(iscell(NULL));

	//cellinit creates or checks
	cell *B=NULL;
	cellinit(&B, 2, 3);
	assert_non_null(B);
	assert_int_equal(NX(B), 2);
	assert_int_equal(NY(B), 3);
	cell *B0=B;
	cellinit(&B, 2, 3);//same size: kept
	assert_ptr_equal(B, B0);
	cellfree(B);

	//cellinit2 with dcell reference creates a dcell with continuous memory
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long mx[2]={3, 2}, my[2]={4, 1};
	dcell *R=randcell(2, 2, mx, my, &rstat);
	dfree(P(R, 1, 0));//make one block empty
	dcell *C=NULL;
	cellinit2(&C, R);
	assert_non_null(C);
	assert_non_null(C->m);
	assert_int_equal(NX(C), 2);
	assert_int_equal(NY(C), 2);
	assert_int_equal(NX(C, 0, 1), 3);
	assert_int_equal(NY(C, 0, 1), 1);
	assert_int_equal(PN(C, 1, 0), 0);//empty placeholder
	assert_int_equal(PN(C->m), 12+3+2);
	dcell *C0=C;
	cellinit2(&C, R);//exists: only checked
	assert_ptr_equal(C, C0);
	dcellfree(C);
	//cellinit2 with ccell reference
	ccell *Rc=ccellnew(2, 1);
	P(Rc, 0)=cnew(2, 2);
	P(Rc, 1)=cnew(3, 1);
	ccell *Cc=NULL;
	cellinit2(&Cc, Rc);
	assert_non_null(Cc->m);
	assert_int_equal(PN(Cc->m), 7);
	assert_int_equal(Cc->m->id, M_COMP);
	ccellfree(Cc);
	ccellfree(Rc);
	//cellinit2 with lcell reference gives a generic cell
	lcell *Rl=lcellnew(1, 2);
	P(Rl, 0)=lnew(2, 2);
	cell *Cl=NULL;
	cellinit2(&Cl, Rl);
	assert_non_null(Cl);
	assert_int_equal(NX(Cl), 1);
	assert_int_equal(NY(Cl), 2);
	assert_null(P(Cl, 0));
	cellfree(Cl);
	lcellfree(Rl);

	//celldim
	long nx, ny, *nxs, *nys;
	celldim(R, &nx, &ny, &nxs, &nys);
	assert_int_equal(nx, 5);
	assert_int_equal(ny, 5);
	assert_int_equal(nxs[0], 3);
	assert_int_equal(nxs[1], 2);
	assert_int_equal(nys[0], 4);
	assert_int_equal(nys[1], 1);
	free(nxs); free(nys);
	celldim((dcell *)NULL, &nx, &ny, &nxs, &nys);
	assert_int_equal(nx, 0);
	assert_int_equal(ny, 0);
	assert_null(nxs);
	assert_null(nys);
	//cell_is_diag
	assert_false(cell_is_diag(R));
	dfree(P(R, 0, 1));
	assert_true(cell_is_diag(R));
	assert_false(cell_is_diag(P(R, 0)));//not a cell
	dcellfree(R);
}

/** cellref references every supported type, cellhash, cellconvert and cellfree of nested cells. */
static void cell_ref_hash(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	cell *A=cellnew(3, 4);
	dmat *d=randmat(3, 2, &rstat);
	P(A, 0)=(cell *)d;
	cmat *c=cnew(2, 2); crandn(c, 1, &rstat); P(A, 1)=(cell *)c;
	smat *s=snew(4, 1); srandn(s, 1, &rstat); P(A, 2)=(cell *)s;
	zmat *z=znew(1, 3); zrandn(z, 1, &rstat); P(A, 3)=(cell *)z;
	lmat *l=lnew(2, 1); P(l, 1)=7; P(A, 4)=(cell *)l;
	dsp *ds=dspnewrandu(4, 4, 1, 0.5, &rstat); P(A, 5)=(cell *)ds;
	csp *cs=cspnewrandu(3, 3, COMPLEX(1., 1.), 0.5, &rstat); P(A, 6)=(cell *)cs;
	ssp *ss=sspnewrandu(3, 5, 1, 0.5, &rstat); P(A, 7)=(cell *)ss;
	zsp *zs=zspnewrandu(2, 2, FCOMPLEX(1.f, 0.f), 1, &rstat); P(A, 8)=(cell *)zs;
	cell *sub=cellnew(1, 2);//nested cell
	P(sub, 0)=(cell *)randmat(2, 2, &rstat);
	P(A, 9)=sub;
	A->keywords=strdup("a=1");
	A->m=(cell *)dnew(2, 2);

	cell *B=cellref(A);
	assert_non_null(B);
	assert_int_equal(NX(B), 3);
	assert_int_equal(NY(B), 4);
	assert_string_equal(B->keywords, "a=1");
	assert_ptr_equal(P((dmat *)B->m), P((dmat *)A->m));
	assert_ptr_equal(P((dmat *)P(B, 0)), P(d));
	assert_ptr_equal(P((cmat *)P(B, 1)), P(c));
	assert_ptr_equal(P((smat *)P(B, 2)), P(s));
	assert_ptr_equal(P((zmat *)P(B, 3)), P(z));
	assert_ptr_equal(P((lmat *)P(B, 4)), P(l));
	assert_ptr_equal(((dsp *)P(B, 5))->px, ds->px);
	assert_ptr_equal(((csp *)P(B, 6))->px, cs->px);
	assert_ptr_equal(((ssp *)P(B, 7))->px, ss->px);
	assert_ptr_equal(((zsp *)P(B, 8))->px, zs->px);
	assert_true(iscell(P(B, 9)));
	assert_ptr_not_equal(P(B, 9), sub);
	assert_ptr_equal(P((dmat *)P(P(B, 9), 0)), P((dmat *)P(sub, 0)));
	assert_null(P(B, 10));
	assert_int_equal(mem_nref(d->mem), 2);
	assert_null(cellref((cell *)NULL));
	dmat *dr=(dmat *)cellref(d);//non-cell is referenced
	assert_ptr_equal(P(dr), P(d));
	dfree(dr);

	//cellhash chains the hash of elements in order; unknown types are skipped.
	uint32_t key=cellhash(A, 5);
	uint32_t k2=dhash(d, 5);
	k2=chash(c, k2);
	k2=shash(s, k2);
	k2=zhash(z, k2);
	k2=dhash((dmat *)P(sub, 0), k2);
	assert_int_equal(key, k2);
	assert_int_equal(cellhash(B, 5), key);
	assert_int_equal(cellhash((cell *)NULL, 9), 9);
	cellfree(B);
	assert_int_equal(mem_nref(d->mem), 1);
	cellfree(A);

	//cellconvert converts nested elements in place
	cell *C=cellnew(2, 1);
	P(C, 0)=(cell *)randmat(3, 3, &rstat);
	P(C, 1)=cellnew(1, 1);
	P(P(C, 1), 0)=(cell *)dnew(2, 2);
	mapcell *M=(mapcell *)cellconvert(C, map_convert_cell);
	assert_ptr_equal(M, C);
	assert_int_equal(P(M, 0)->id, M_MAP);
	assert_int_equal(P(P(C, 1), 0)->id, M_MAP);
	assert_null(cellconvert(NULL, map_convert_cell));
	cellfree(M);
}

/** cellresize preserves overlapping blocks and frees the rest. */
static void cell_resize(void **state){
	(void)state;
	dcell *A=dcellnew(3, 2);
	for(long i=0; i<PN(A); i++){
		P(A, i)=dnew(1, 1);
		P(P(A, i), 0)=(real)i;
	}
	//Change ny with the same nx: linear realloc
	cellresize(A, 3, 3);
	assert_int_equal(NY(A), 3);
	for(long i=0; i<6; i++) assert_float_equal(P(P(A, i), 0), i, 0);
	for(long i=6; i<9; i++) assert_null(P(A, i));
	cellresize(A, 3, 1);//shrink: frees the dropped blocks
	assert_int_equal(PN(A), 3);
	for(long i=0; i<3; i++) assert_float_equal(P(P(A, i), 0), i, 0);
	//column vector growing in nx keeps the layout
	cellresize(A, 5, 1);
	assert_int_equal(NX(A), 5);
	assert_float_equal(P(P(A, 2), 0), 2, 0);
	assert_null(P(A, 4));
	//General reshape: blocks are copied by 2-d index
	P(A, 4)=dnew(1, 1); P(P(A, 4), 0)=4;
	cellresize(A, 2, 2);
	assert_int_equal(NX(A), 2);
	assert_int_equal(NY(A), 2);
	assert_float_equal(P(P(A, 0, 0), 0), 0, 0);
	assert_float_equal(P(P(A, 1, 0), 0), 1, 0);
	assert_null(P(A, 0, 1));
	assert_null(P(A, 1, 1));
	cellresize(A, 3, 4);
	assert_float_equal(P(P(A, 1, 0), 0), 1, 0);
	assert_null(P(A, 2, 3));
	dcellresize(A, 0, 0);
	assert_int_equal(PN(A), 0);
	dcellfree(A);
	//resizing a cell with m set is refused
	dcell *B=dcellnew_same(2, 1, 3, 3);
	cellresize(B, 4, 1);
	assert_int_equal(NX(B), 2);
	dcellfree(B);
}

/** Cell constructors of mat.c for several types: cell_cast, cellnew2/3/_same/_file, cellref, cellset, celltrans. */
static void cell_constructors(void **state){
	(void)state;
	long nnx[4]={2, 3, 2, 3};
	long nny[4]={1, 1, 2, 2};
	dcell *A=dcellnew3(2, 2, nnx, nny);
	assert_non_null(A->m);
	assert_int_equal(PN(A->m), 2+3+4+6);
	for(int i=0; i<4; i++){
		assert_int_equal(NX(A, i), nnx[i]);
		assert_int_equal(NY(A, i), nny[i]);
	}
	assert_ptr_equal(P(P(A, 1)), P(A->m)+2);
	dcell *A2=dcellnew3(2, 1, nnx, NULL);//NULL nny means 1 column
	assert_int_equal(NY(A2, 1), 1);
	dcellfree(A2);
	dcellset(A, 2);
	assert_float_equal(dsum(A->m), 2*15, 0);
	assert_float_equal(dcellsum(A), 2*15, 0);
	dmat *se=dcellsum_each(A);
	for(int i=0; i<4; i++) assert_float_equal(P(se, i), 2*nnx[i]*nny[i], 0);
	dfree(se);
	assert_null(dcellsum_each(NULL));
	assert_float_equal(dcellsum(NULL), 0, 0);
	//cellref shares data and m
	A->keywords=strdup("k=2");
	dcell *Ar=dcellref(A);
	assert_ptr_equal(P(Ar->m), P(A->m));
	assert_ptr_equal(P(P(Ar, 3)), P(P(A, 3)));
	assert_string_equal(Ar->keywords, "k=2");
	assert_null(dcellref(NULL));
	dcellfree(Ar);
	//celltrans transposes blocks and layout
	P(P(A, 1, 1), 2, 1)=-5;
	dcell *At=dcelltrans(A);
	assert_int_equal(NX(At, 0, 1), nny[1]);
	assert_int_equal(NY(At, 0, 1), nnx[1]);
	assert_int_equal(NX(At, 1, 1), nny[3]);
	assert_float_equal(P(P(At, 1, 1), 1, 2), -5, 0);
	assert_null(dcelltrans(NULL));
	dmat *Am=dcell2m(A), *Atm=dcell2m(At), *Amt=dtrans(Am);
	assert_dmat_near(Atm, Amt, 0);
	dfree(Am); dfree(Atm); dfree(Amt);
	dcellfree(At);
	//cellnew2 copies the layout; empty blocks become placeholders
	dfree(P(A, 2));
	dcell *B=dcellnew2(A);
	assert_int_equal(PN(B->m), 2+3+6);
	assert_non_null(P(B, 2));
	assert_int_equal(PN(B, 2), 0);
	assert_float_equal(dsum(B->m), 0, 0);
	assert_null(dcellnew2(NULL));
	dcellfree(B);
	//cell_cast
	assert_ptr_equal(dcell_cast((cell *)A), A);
	cell *mixed=cellnew(2, 1);
	P(mixed, 0)=(cell *)dnew(1, 1);
	P(mixed, 1)=(cell *)cnew(1, 1);
	assert_null(dcell_cast(mixed));
	assert_null(ccell_cast(mixed));
	assert_null(dcell_cast((cell *)P(mixed, 0)));//not a cell
	cellfree(mixed);
	dcellfree(A);
	//cellnew_file without a file name
	dcell *F=dcellnew_file(2, 1, (long *)-3, (long *)-2, "x=1", NULL);
	assert_non_null(F);
	assert_string_equal(F->keywords, "x=1");
	assert_int_equal(PN(F->m), 12);
	assert_null(F->fp);
	dcellfree(F);
	F=dcellnew_file(2, 1, (long *)-3, (long *)-2, NULL, "-skip");
	assert_null(F->fp);
	dcellfree(F);
	assert_null(dcellnew_file(0, 1, (long *)-3, (long *)-2, NULL, NULL));
	ccell *Fc=ccellnewsame_file(1, 3, 2, 2, "y=2", NULL);
	assert_non_null(Fc);
	assert_string_equal(Fc->keywords, "y=2");
	assert_int_equal(PN(Fc->m), 12);
	ccellfree(Fc);
	Fc=ccellnewsame_file(1, 1, 2, 2, NULL, "-%d", 1);
	assert_null(Fc->keywords);
	ccellfree(Fc);
	assert_null(ccellnewsame_file(1, 0, 2, 2, NULL, NULL));
	//other types
	scell *S=scellnew_same(2, 2, 2, 1);
	scellset(S, 1.5f);
	assert_float_equal(scellsum(S), 12, 1e-6);
	scell *St=scelltrans(S);
	assert_int_equal(NX(St, 0), 1);
	scellfree(St);
	scellfree(S);
	zcell *Z=zcellnew_same(1, 2, 2, 2);
	zcellset(Z, FCOMPLEX(1.f, 2.f));
	assert_float_equal(cimagf(zcellsum(Z)), 16, 1e-6);
	zcellfree(Z);
	lcell *L=lcellnew_same(2, 1, 3, 1);
	lcellset(L, 3);
	assert_int_equal(lcellsum(L), 18);
	lcell *L2=lcellref(L);
	assert_int_equal(lsum(P(L2, 1)), 9);
	lcellfree(L2);
	lcellfree(L);
	dcellset(NULL, 1);//no-op
}

/** cellreduce, celldropempty, 2cellref, 2cell, cellsub, cell_col. */
static void cell_reshape(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 3);
	const long mx[3]={2, 3, 1}, my[2]={2, 3};
	dcell *A=randcell(3, 2, mx, my, &rstat);
	dmat *Am=dcell2m(A);
	assert_int_equal(NX(Am), 6);
	assert_int_equal(NY(Am), 5);
	//reduce along rows: 1x2 cell
	dcell *R1=dcellreduce(A, 1);
	assert_int_equal(NX(R1), 1);
	assert_int_equal(NY(R1), 2);
	assert_cell_dense(R1, Am, 0);
	dcell *R2=dcellreduce(A, 2);
	assert_int_equal(NX(R2), 3);
	assert_int_equal(NY(R2), 1);
	assert_cell_dense(R2, Am, 0);
	assert_null(dcellreduce(NULL, 1));
	dcellfree(R1); dcellfree(R2);
	//An empty block still keeps the layout as zeros
	dcell *Ae=dcelldup(A);
	dzero(P(Ae, 1, 0));
	dcell *Ae2=dcelldup(Ae);
	dfree(P(Ae2, 1, 0));
	dcell *Re=dcellreduce(Ae2, 1);
	dmat *Aem=dcell2m(Ae);
	assert_cell_dense(Re, Aem, 0);
	dcellfree(Re);
	Re=dcellreduce(Ae2, 2);
	assert_cell_dense(Re, Aem, 0);
	dcellfree(Re);
	dfree(Aem);
	dcellfree(Ae);

	//2cell splits a dense matrix with the layout of A
	dcell *B=NULL;
	d2cell(&B, Am, A);
	assert_true(dcelldiff(A, B)<TOL);
	assert_float_equal(dcelldiff(A, B), 0, 0);
	dzero(P(B, 0));
	d2cell(&B, Am, NULL);//B is used as reference
	assert_float_equal(dcelldiff(A, B), 0, 0);
	dcellfree(B);
	//2cell with an empty block in the reference allocates it
	dcell *C=NULL;
	dcellcp(&C, Ae2);
	dfree(P(C, 1, 0));
	dcell *C2=dcelldup(Ae2);
	dmat *C2m=dcell2m(C2);
	d2cell(&C, C2m, NULL);
	assert_non_null(P(C, 1, 0));
	dfree(C2m);
	dcellfree(C);
	dcellfree(C2);
	dcellfree(Ae2);

	//2cellref references a vector
	dmat *v=dlinspace(0, 1, 10);
	long dims[3]={2, 5, 3};
	dcell *vr=d2cellref(v, dims, 3);
	assert_int_equal(NX(vr, 1), 5);
	assert_ptr_equal(P(P(vr, 2)), P(v)+7);
	assert_float_equal(P(P(vr, 1), 0), 2, 0);
	dcellfree(vr);
	dfree(v);

	//cellsub references a sub range
	dcell *S=dcellsub(A, 1, 2, 1, 0);
	assert_int_equal(NX(S), 2);
	assert_int_equal(NY(S), 1);
	assert_ptr_equal(P(P(S, 0, 0)), P(P(A, 1, 1)));
	dcellfree(S);
	S=dcellsub(A, 0, 0, 0, 0);//all
	assert_int_equal(PN(S), 6);
	dcellfree(S);
	assert_null(dcellsub(NULL, 0, 1, 0, 1));

	//celldropempty
	dcell *D=dcelldup(A);
	dmat *Dm=D->m;//contiguous storage (may be NULL) is transferred to the new container
	real *pD20=P(P(D, 2, 0));
	dfree(P(D, 1, 0)); dfree(P(D, 1, 1)); P(D, 1, 1)=dnew(0, 0);//empty row 1 (NULL and 0x0 blocks)
	dcelldropempty(&D, 1);
	assert_int_equal(NX(D), 2);
	assert_int_equal(NX(D, 1, 0), 1);
	assert_ptr_equal(D->m, Dm);
	assert_ptr_equal(P(P(D, 1, 0)), pD20);
	dfree(P(D, 0, 1)); dfree(P(D, 1, 1));//empty column 1
	dcelldropempty(&D, 2);
	assert_int_equal(NY(D), 1);
	assert_int_equal(NX(D), 2);
	dcelldropempty(&D, 1);//nothing to drop
	assert_int_equal(NX(D), 2);
	dfree(P(D, 0)); dfree(P(D, 1));
	dcell *D2=dcelldup(D);
	dcelldropempty(&D, 1);//everything is dropped
	assert_null(D);
	dcelldropempty(&D2, 2);
	assert_null(D2);
	dcelldropempty(NULL, 1);

	//cell_col with and without m
	dcell *G=dcellnew_same(3, 2, 2, 2);
	for(long i=0; i<PN(G->m); i++) P(G->m, i)=(real)i;
	dmat *g=dcell_col(G, 1);
	assert_int_equal(NX(g), 4);
	assert_int_equal(NY(g), 3);
	assert_float_equal(P(g, 0, 0), 12, 0);
	assert_float_equal(P(g, 3, 2), 23, 0);
	dfree(g);
	dcell *G2=dcelldup(G);
	dfree(G2->m);//force copying block by block
	g=dcell_col(G2, 1);
	assert_float_equal(P(g, 0, 0), 12, 0);
	assert_float_equal(P(g, 3, 2), 23, 0);
	dfree(g);
	assert_null(dcell_col(NULL, 0));
	dcellfree(G2);
	dcellfree(G);
	dfree(Am);
	dcellfree(A);
}

/** cellcat, cellcat2 and cellcat_each. */
static void cell_cat(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 4);
	const long mx[2]={2, 3}, my[2]={1, 2};
	dcell *A=randcell(2, 2, mx, my, &rstat);
	dcell *B=randcell(2, 2, mx, my, &rstat);
	dmat *Am=dcell2m(A), *Bm=dcell2m(B);
	dcell *C=dcellcat(A, B, 1);
	assert_int_equal(NX(C), 4);
	dmat *ABx=dcat(Am, Bm, 1);
	assert_cell_dense(C, ABx, 0);
	dcellfree(C);
	C=dcellcat(A, B, 2);
	assert_int_equal(NY(C), 4);
	dmat *ABy=dcat(Am, Bm, 2);
	assert_cell_dense(C, ABy, 0);
	dcellfree(C);
	C=dcellcat(NULL, B, 1);
	assert_float_equal(dcelldiff(C, B), 0, 0);
	dcellfree(C);
	assert_null(dcellcat(NULL, NULL, 1));
	//cellcat2 appends in place
	C=NULL;
	dcellcat2(&C, A, 1);//copies when C is NULL
	dcellcat2(&C, B, 2);
	dcellcat2(&C, NULL, 2);//no-op
	assert_int_equal(NX(C), 2);
	assert_int_equal(NY(C), 4);
	assert_cell_dense(C, ABy, 0);
	dcellfree(C);
	C=dcellnew_same(2, 1, 2, 1);//m is freed with a warning
	dcellcat2(&C, B, 2);
	assert_null(C->m);
	assert_int_equal(NY(C), 3);
	assert_float_equal(P(P(C, 1, 2), 0), P(P(B, 1, 1), 0), 0);
	dcellfree(C);
	//cellcat_each concatenates each block
	dcell *E=dcellcat_each(A, B, 1);
	assert_int_equal(NX(E, 1, 1), 6);
	dmat *e11=dcat(P(A, 1, 1), P(B, 1, 1), 1);
	assert_dmat_near(P(E, 1, 1), e11, 0);
	dfree(e11);
	dcellfree(E);
	E=dcellcat_each(NULL, B, 1);
	assert_float_equal(dcelldiff(E, B), 0, 0);
	dcellfree(E);
	E=dcellcat_each(A, NULL, 1);
	assert_float_equal(dcelldiff(E, A), 0, 0);
	dcellfree(E);
	assert_null(dcellcat_each(NULL, NULL, 2));
	dfree(Am); dfree(Bm); dfree(ABx); dfree(ABy);
	dcellfree(A); dcellfree(B);
}

/*D is clipped to [-0.5,0.5].*/
static int dcellmax_abs_le(const dcell *D){
	for(long i=0; i<PN(D); i++){
		if(dmaxabs(P(D, i))>0.5) return 0;
	}
	return 1;
}
/** Element wise cell math: cellnorm, celldot, cellcwm, celldropzero, celldiff, cellclip, cellcwpow, cellscale, cellzero. */
static void cell_elementwise(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 5);
	const long mx[3]={4, 2, 3}, my[1]={1};
	dcell *A=randcell(3, 1, mx, my, &rstat);
	dcell *B=randcell(3, 1, mx, my, &rstat);
	dmat *Am=dcell2m(A), *Bm=dcell2m(B);
	assert_float_equal(dcellnorm(A), dnorm(P(A, 0))+dnorm(P(A, 1))+dnorm(P(A, 2)), TOL);
	assert_float_equal(dcelldot(A, B), ddot(Am, Bm), TOL);
	assert_float_equal(dcelldot(A, NULL), 0, 0);
	{//empty (NULL) blocks contribute 0
		dcell *A2=dcelldup(A), *B2=dcelldup(B);
		dfree(P(A2, 1)); dfree(P(B2, 1));
		real ref=ddot(P(A, 0), P(B, 0))+ddot(P(A, 2), P(B, 2));
		assert_float_equal(dcelldot(A2, B2), ref, TOL);
		dcell *B3=dcellnew(2, 1);
		assert_float_equal(dcelldot(A, B3), 0, 0);//mismatch: error and return 0
		dcellfree(A2); dcellfree(B2); dcellfree(B3);
	}
	//cwm
	dcell *C=dcelldup(B);
	dcellcwm(C, A);
	dmat *Cm=ddup(Bm);
	dcwm(Cm, Am);
	assert_cell_dense(C, Cm, TOL);
	dfree(Cm);
	//celldiff is zero for identical cells and grows with difference
	assert_float_equal(dcelldiff(A, A), 0, 0);
	real dd=dcelldiff(A, B);
	assert_true(dd>0.1);
	//cellclip
	dcell *D=dcelldup(A);
	long nexpect=0;
	for(long i=0; i<PN(Am); i++) if(P(Am, i)>0.5||P(Am, i)<-0.5) nexpect++;
	assert_int_equal(dcellclip(D, -0.5, 0.5), nexpect);
	assert_true(dcellmax_abs_le(D));
	assert_int_equal(dcellclip(D, -INFINITY, INFINITY), 0);
	assert_int_equal(dcellclip(NULL, -1, 1), 0);
	//cellcwpow
	dcellcp(&D, A);
	dcellcwpow(D, 2);
	dmat *A2=ddup(Am);
	dcwm(A2, Am);
	assert_cell_dense(D, A2, TOL);
	dcellcwpow(NULL, 2);
	dfree(A2);
	//cellscale and cellzero
	dcellscale(D, -2);
	assert_float_equal(dcellsum(D), -2*dsumsq(Am), TOL*10);
	dcellzero(D);
	assert_float_equal(dcellnorm(D), 0, 0);
	dcellscale((dcell *)NULL, 2);
	//celldropzero
	dset(P(D, 1), 1e-5);
	P(P(D, 2), 1)=1;
	dcelldropzero(D, 1e-3);
	assert_null(P(D, 0));
	assert_null(P(D, 1));
	assert_non_null(P(D, 2));
	//cellcp from NULL zeros the destination
	dcellcp(&C, NULL);
	assert_float_equal(dcellnorm(C), 0, 0);
	dcellfree(C); dcellfree(D);
	dfree(Am); dfree(Bm);
	dcellfree(A); dcellfree(B);
}

/** celladd handles dense, sparse, single cell and nested cases. */
static void cell_add(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 6);
	//cell + cell allocates the output
	const long mx[2]={3, 2}, my[2]={2, 3};
	dcell *A=randcell(2, 2, mx, my, &rstat);
	dcell *B=randcell(2, 2, mx, my, &rstat);
	dmat *Am=dcell2m(A), *Bm=dcell2m(B);
	dcell *C=NULL;
	dcelladd(&C, 1, A, 2);
	dmat *ref=NULL;
	dadd(&ref, 0, Am, 2);
	assert_cell_dense(C, ref, TOL);
	dcelladd(&C, 3, B, -1);
	dadd(&ref, 3, Bm, -1);
	assert_cell_dense(C, ref, TOL);
	dcelladd(&C, 0.5, (dcell *)NULL, 1);//only scale
	dscale(ref, 0.5);
	assert_cell_dense(C, ref, TOL);
	dcellfree(C);
	//single element cell to matrix and vice versa
	dcell *A1=dcellnew(1, 1);
	P(A1, 0)=randmat(3, 3, &rstat);
	dmat *m=dnew(3, 3);
	dcelladd(&m, 1, A1, 2);
	dmat *m2=NULL;
	dadd(&m2, 0, P(A1, 0), 2);
	assert_dmat_near(m, m2, TOL);
	dcell *A1c=dcellnew(1, 1);
	dcelladd(&A1c, 1, m, 1);//matrix added to the first element
	assert_dmat_near(P(A1c, 0), m2, TOL);
	dcellfree(A1c);
	//sparse added to dense
	dsp *S=dspnewrandu(3, 3, 1, 0.4, &rstat);
	dmat *Sf=NULL;
	dspfull(&Sf, S, 'n', 1);
	dcelladd(&m, 2, S, 3);
	dadd(&m2, 2, Sf, 3);
	assert_dmat_near(m, m2, TOL);
	//sparse added to sparse
	dsp *S2=dspnewrandu(3, 3, 1, 0.4, &rstat);
	dmat *S2f=NULL;
	dspfull(&S2f, S2, 'n', 1);
	dsp *S3=dspdup(S2);
	dcelladd(&S3, -1, S, 0.5);
	dmat *S3f=NULL;
	dspfull(&S3f, S3, 'n', 1);
	dadd(&S2f, -1, Sf, 0.5);
	assert_dmat_near(S3f, S2f, TOL);
	//sparse cell plus dense cell
	dspcell *Sc=dspcellnew(2, 1);
	P(Sc, 0)=dspref(S);
	P(Sc, 1)=dspref(S2);
	dcell *Dc=dcellnew_same(2, 1, 3, 3);
	dcelladd(&Dc, 1, Sc, 1);
	dmat *Dm=dcell2m(Dc), *Scm=dcell2m(Sc);
	assert_dmat_near(Dm, Scm, 0);
	dfree(Dm); dfree(Scm);
	//cellscale with sparse elements
	dcellscale(Sc, 2);
	dmat *Sff=NULL;
	dspfull(&Sff, P(Sc, 0), 'n', 1);
	dscale(Sf, 2);
	assert_dmat_near(Sff, Sf, TOL);
	dfree(Sff);

	dfree(m); dfree(m2); dfree(Sf); dfree(S2f); dfree(S3f);
	dspfree(S); dspfree(S2); dspfree(S3);
	dspcellfree(Sc); dcellfree(Dc);
	dcellfree(A1);
	dfree(Am); dfree(Bm); dfree(ref);
	dcellfree(A); dcellfree(B);
}

/** celladdI with dense and sparse diagonal blocks. */
static void cell_addI(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 7);
	cell *A=cellnew(3, 3);
	P(A, 0, 0)=(cell *)randmat(3, 3, &rstat);
	P(A, 1, 1)=(cell *)dspnewrandu(2, 2, 1, 0.5, &rstat);
	P(A, 1, 0)=(cell *)randmat(2, 3, &rstat);
	//A(2,2) is empty
	P(A, 0, 1)=(cell *)dnew(3, 2);
	dmat *before=dcell2m(A);
	dcelladdI(A, 2.5);
	dmat *after=dcell2m(A);
	for(long i=0; i<5; i++) P(before, i, i)+=2.5;
	assert_dmat_near(after, before, TOL);
	assert_int_equal(P(A, 1, 1)->id, M_DSP);
	dfree(before); dfree(after);
	dcelladdI((cell *)NULL, 1);
	cellfree(A);
}

/*Reference dense product op(A)*op(B) with trans.*/
static void check_cellmm(const cell *A, const cell *B, const char *trans, real alpha){
	dcell *C=NULL;
	dcellmm(&C, A, B, trans, alpha);
	dmat *Am=dcell2m(A), *Bm=dcell2m(B);
	dmat *ref=NULL;
	dmm(&ref, 0, Am, Bm, trans, alpha);
	assert_cell_dense(C, ref, TOL);
	//accumulate
	dcellmm(&C, A, B, trans, alpha);
	dscale(ref, 2);
	assert_cell_dense(C, ref, TOL);
	dfree(ref); dfree(Am); dfree(Bm);
	dcellfree(C);
}
/** cellmm on cells of dense and sparse blocks for all transposes against the dense product. */
static void cell_mm(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 8);
	const long ma[2]={3, 2}, na[3]={2, 4, 1};
	const long nb[2]={3, 2};
	dcell *A=randcell(2, 3, ma, na, &rstat);//5x7
	dcell *B=randcell(3, 2, na, nb, &rstat);//7x5
	dcell *At=dcelltrans(A);
	dcell *Bt=dcelltrans(B);
	check_cellmm((cell *)A, (cell *)B, "nn", 1);
	check_cellmm((cell *)At, (cell *)B, "tn", 2);
	check_cellmm((cell *)A, (cell *)Bt, "nt", -1);
	check_cellmm((cell *)At, (cell *)Bt, "tt", 0.5);
	//cell2 returns a new cell
	dcell *C=dcellmm2(A, B, "nn");
	dcell *C2=NULL;
	dcellmm(&C2, A, B, "nn", 1);
	assert_float_equal(dcelldiff(C, C2), 0, 0);
	dcellfree(C); dcellfree(C2);

	//sparse blocks: A is a sparse cell
	dspcell *As=dspcellnew(2, 2);
	P(As, 0, 0)=dspnewrandu(3, 3, 1, 0.5, &rstat);
	P(As, 1, 0)=dspnewrandu(2, 3, 1, 0.5, &rstat);
	P(As, 0, 1)=dspnewrandu(3, 2, 1, 0.5, &rstat);
	P(As, 1, 1)=dspnewrandu(2, 2, 1, 0.5, &rstat);
	const long mb[2]={3, 2}, nb2[1]={4};
	dcell *Bv=randcell(2, 1, mb, nb2, &rstat);
	check_cellmm((cell *)As, (cell *)Bv, "nn", 1);
	dcell *Bvt=dcelltrans(Bv);
	check_cellmm((cell *)Bvt, (cell *)As, "nn", 1);//dense times sparse
	check_cellmm((cell *)Bvt, (cell *)As, "nt", 1);
	check_cellmm((cell *)As, (cell *)Bv, "tn", 1);
	//sparse times sparse into sparse result
	dspcell *Ss=NULL;
	dcellmm(&Ss, As, As, "nn", 1);
	assert_int_equal(P(Ss, 0, 0)->id, M_DSP);
	dmat *Ssm=dcell2m(Ss), *Asm=dcell2m(As), *ref=NULL;
	dmm(&ref, 0, Asm, Asm, "nn", 1);
	assert_dmat_near(Ssm, ref, TOL);
	dfree(Ssm); dfree(ref);
	dspcellfree(Ss);
	//sparse times sparse into dense result
	dmat *Sdm=dnew(3, 2);
	dcellmm(&Sdm, P(As, 0, 0), P(As, 0, 1), "tn", 2);
	dmat *a00=NULL, *a01=NULL;
	dspfull(&a00, P(As, 0, 0), 'n', 1);
	dspfull(&a01, P(As, 0, 1), 'n', 1);
	dmm(&ref, 0, a00, a01, "tn", 2);
	assert_dmat_near(Sdm, ref, TOL);
	dfree(Sdm); dfree(ref); dfree(a00); dfree(a01);
	//single element cell times matrix and mm_cell
	dcell *A1=dcellnew(1, 1);
	P(A1, 0)=randmat(3, 4, &rstat);
	dmat *x=randmat(4, 2, &rstat);
	dmat *y=NULL, *yref=NULL;
	dcellmm(&y, A1, x, "nn", 1);
	dmm(&yref, 0, P(A1, 0), x, "nn", 1);
	assert_dmat_near(y, yref, TOL);
	dmm_cell(&y, (cell *)A1, x, "nn", -1);
	assert_float_equal(dnorm(y), 0, TOL);
	dmm_cell(NULL, (cell *)A1, x, "nn", -1);//no-op
	//cell output with a single element, allocated or not
	dcell *yc=dcellnew(1, 1);
	P(yc, 0)=dnew(3, 2);
	dcellmm(&yc, P(A1, 0), x, "nn", 1);
	assert_dmat_near(P(yc, 0), yref, TOL);
	dfree(P(yc, 0));
	dcellmm(&yc, P(A1, 0), x, "nn", 1);//NULL element is allocated
	assert_non_null(P(yc, 0));
	assert_dmat_near(P(yc, 0), yref, TOL);
	dcellfree(yc);
	//Unsupported shape combinations raise an error and return instead of recursing forever
	yc=dcellnew(2, 1);
	dcellmm(&yc, P(A1, 0), x, "nn", 1);//cell output with 2 elements for matrix product
	assert_null(P(yc, 0));
	assert_null(P(yc, 1));
	dcellfree(yc);
	dcell *A2=dcellnew(2, 1);
	P(A2, 0)=randmat(3, 4, &rstat);
	P(A2, 1)=randmat(3, 4, &rstat);
	dmat *ybad=NULL;
	dcellmm(&ybad, A2, x, "nn", 1);//multi-element cell times matrix
	assert_null(ybad);
	dcellfree(A2);
	//dimension mismatch is an error and leaves the output untouched
	dmat *xbad=randmat(5, 2, &rstat);
	dcellmm(&ybad, A1, xbad, "nn", 1);
	assert_null(ybad);
	dfree(xbad);
	dfree(x); dfree(y); dfree(yref);
	dcellfree(A1);
	dcellmm(&y, (dcell *)NULL, A, "nn", 1);//no-op
	assert_null(y);
	dfree(Asm);
	dcellfree(Bvt); dcellfree(Bv);
	dspcellfree(As);
	dcellfree(A); dcellfree(B); dcellfree(At); dcellfree(Bt);
}
/** cellmm on complex and single precision cells. */
static void cell_mm_types(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 9);
	ccell *A=ccellnew_same(2, 1, 3, 2);
	ccell *B=ccellnew_same(2, 1, 3, 2);
	crandn(A->m, 1, &rstat);
	crandn(B->m, 1, &rstat);
	ccell *C=NULL;
	ccellmm(&C, A, B, "cn", 1);//2x2 blocks: sum of A(i)'*B(i)
	cmat *ref=NULL;
	for(int i=0; i<2; i++){
		cmm(&ref, 1, P(A, i), P(B, i), "cn", 1);
	}
	assert_int_equal(NX(C), 1);
	assert_int_equal(NY(C), 1);
	for(long i=0; i<PN(ref); i++){
		assert_float_equal(cabs(P(P(C, 0), i)-P(ref, i)), 0, TOL);
	}
	cfree(ref);
	ccellfree(C);
	scell *As=scellnew_same(1, 2, 2, 2);
	scell *Bs=scellnew_same(2, 1, 2, 3);
	srandn(As->m, 1, &rstat);
	srandn(Bs->m, 1, &rstat);
	scell *Cs=scellmm2(As, Bs, "nn");
	smat *rs=NULL;
	smm(&rs, 0, P(As, 0), P(Bs, 0), "nn", 1);
	smm(&rs, 1, P(As, 1), P(Bs, 1), "nn", 1);
	for(long i=0; i<PN(rs); i++){
		assert_float_equal(P(P(Cs, 0), i), P(rs, i), 1e-5);
	}
	sfree(rs);
	scellfree(Cs);
	scellfree(As); scellfree(Bs);
	ccellfree(A); ccellfree(B);
}

/** cellinv, cellinvspd, cellinvspd_each, cellpinv and cellsvd_pow against dense results. */
static void cell_inverse(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 10);
	//SPD block matrix
	dmat *R=randmat(6, 6, &rstat);
	dmat *S=NULL;
	dmm(&S, 0, R, R, "tn", 1);
	daddI(S, 1);
	dcell *Sc=dcellnew_same(2, 2, 3, 3);
	d2cell(&Sc, S, NULL);
	dmat *Si=dinv(S);
	dcell *Ci=dcellinvspd(Sc);
	assert_cell_dense(Ci, Si, 1e-8);
	dcellfree(Ci);
	Ci=dcellinv(Sc);
	assert_cell_dense(Ci, Si, 1e-8);
	dcellfree(Ci);
	assert_null(dcellinvspd(NULL));
	//cellsvd_pow on a full block matrix
	dcell *Sp=dcelldup(Sc);
	dcellsvd_pow(Sp, -1);
	assert_cell_dense(Sp, Si, 1e-8);
	dcellfree(Sp);
	//block diagonal: each block is processed independently
	dcell *Sd=dcelldup(Sc);
	dfree(P(Sd, 1, 0)); dfree(P(Sd, 0, 1));
	dcell *Sdp=dcelldup(Sd);
	dcellsvd_pow(Sdp, -1);
	for(int i=0; i<2; i++){
		dmat *inv=dinv(P(Sd, i, i));
		assert_dmat_near(P(Sdp, i, i), inv, 1e-8);
		dfree(inv);
	}
	assert_int_equal(PN(Sdp, 1, 0), 0);
	//cellinvspd_each inverts every block. Empty (NULL or 0x0) blocks are skipped.
	dcell *Se=dcellnew(4, 1);
	P(Se, 0)=ddup(P(Sd, 0, 0));
	P(Se, 1)=ddup(P(Sd, 1, 1));
	P(Se, 3)=dnew(0, 0);
	dcell *Sde=dcellinvspd_each(Se);
	for(int i=0; i<2; i++){
		dmat *inv=dinv(P(Se, i));
		assert_dmat_near(P(Sde, i), inv, 1e-8);
		dfree(inv);
	}
	assert_int_equal(PN(Sde, 2), 0);
	assert_int_equal(PN(Sde, 3), 0);
	assert_null(dcellinvspd_each(NULL));
	dcellfree(Se); dcellfree(Sde); dcellfree(Sdp); dcellfree(Sd);
	//pseudo inverse of a tall block matrix
	dmat *T=randmat(8, 3, &rstat);
	dcell *Tc=dcellnew(2, 1);
	P(Tc, 0)=dsub(T, 0, 5, 0, 3);
	P(Tc, 1)=dsub(T, 5, 3, 0, 3);
	dcell *Tp=dcellpinv(Tc, NULL);
	dmat *Tpm=dcell2m(Tp);
	dmat *Tpref=dpinv(T, NULL);
	assert_dmat_near(Tpm, Tpref, 1e-8);
	dfree(Tpm); dfree(Tpref); dcellfree(Tp);
	//with diagonal weighting given per block
	dcell *W=dcellnew(2, 1);
	P(W, 0)=dnew(5, 1); dset(P(W, 0), 2);
	P(W, 1)=dnew(3, 1); dset(P(W, 1), 0.5);
	dmat *Wm=dcell2m(W);
	Tp=dcellpinv(Tc, W);
	Tpm=dcell2m(Tp);
	Tpref=dpinv(T, Wm);
	assert_dmat_near(Tpm, Tpref, 1e-8);
	dfree(Tpm); dfree(Tpref); dcellfree(Tp);
	//with a full weighting matrix as a 2x2 cell
	dmat *Wf=dnew(8, 8);
	for(int i=0; i<8; i++) P(Wf, i, i)=P(Wm, i);
	dcell *Wc=dcellnew(2, 2);
	P(Wc, 0, 0)=dsub(Wf, 0, 5, 0, 5);
	P(Wc, 1, 0)=dsub(Wf, 5, 3, 0, 5);
	P(Wc, 0, 1)=dsub(Wf, 0, 5, 5, 3);
	P(Wc, 1, 1)=dsub(Wf, 5, 3, 5, 3);
	Tp=dcellpinv(Tc, Wc);
	Tpm=dcell2m(Tp);
	Tpref=dpinv(T, Wf);
	assert_dmat_near(Tpm, Tpref, 1e-8);
	dfree(Tpm); dfree(Tpref); dcellfree(Tp);
	assert_null(dcellpinv(NULL, NULL));
	dcellfree(Wc); dfree(Wf);
	dfree(Wm); dcellfree(W);
	dcellfree(Tc); dfree(T);
	dfree(R); dfree(S); dfree(Si); dcellfree(Sc);
}

/** spcell operations: ref, cast, trans, full, 2sp, sym, sort, droptol, scale and wdot. */
static void cell_spcell(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 11);
	dspcell *A=dspcellnew(2, 2);
	const long n[2]={4, 3};
	for(int iy=0; iy<2; iy++){
		for(int ix=0; ix<2; ix++){
			P(A, ix, iy)=dspnewrandu(n[ix], n[iy], 1, 0.5, &rstat);
		}
	}
	dmat *Am=dcell2m(A);
	//ref shares data
	dspcell *Ar=dspcellref(A);
	assert_ptr_equal(P(Ar, 1)->px, P(A, 1)->px);
	dspcellfree(Ar);
	assert_null(dspcellref(NULL));
	//cast
	assert_ptr_equal(dspcell_cast((cell *)A), A);
	cell *mixed=cellnew(2, 1);
	P(mixed, 0)=(cell *)dnew(2, 2);
	assert_null(dspcell_cast(mixed));
	assert_null(dspcell_cast(NULL));
	cellfree(mixed);
	//trans
	dspcell *At=dspcelltrans(A);
	dmat *Atm=dcell2m(At), *Amt=dtrans(Am);
	assert_dmat_near(Atm, Amt, 0);
	dfree(Atm);
	//full
	dcell *Af=NULL;
	dspcellfull(&Af, A, 'n', 2);
	dmat *Afm=dcell2m(Af);
	dmat *Am2=ddup(Am);
	dscale(Am2, 2);
	assert_dmat_near(Afm, Am2, TOL);
	dfree(Afm); dcellfree(Af);
	dspcellfull(&Af, A, 't', 1);
	Afm=dcell2m(Af);
	assert_dmat_near(Afm, Amt, TOL);
	dfree(Afm); dcellfree(Af);
	//spcell2sp, including a dense block
	dsp *As=dspcell2sp(A);
	dmat *Asf=NULL;
	dspfull(&Asf, As, 'n', 1);
	assert_dmat_near(Asf, Am, 0);
	dfree(Asf); dspfree(As);
	cell *Amix=cellnew(2, 1);
	P(Amix, 0)=(cell *)dspref(P(A, 0, 0));
	P(Amix, 1)=(cell *)randmat(3, 4, &rstat);
	As=dspcell2sp((dspcell *)Amix);
	dspfull(&Asf, As, 'n', 1);
	dmat *Amixm=dcell2m(Amix);
	assert_dmat_near(Asf, Amixm, 0);
	dfree(Amixm); dfree(Asf); dspfree(As);
	cellfree(Amix);
	//sym: (A+A')/2
	dspcell *Sy=dspcellref(A);
	dspcellfree(Sy);
	Sy=dspcellnew(2, 2);
	for(int i=0; i<4; i++) P(Sy, i)=dspdup(P(A, i));
	dspcellsym(&Sy);
	dmat *Sym=dcell2m(Sy);
	dmat *ref=NULL;
	dadd(&ref, 0, Am, 0.5);
	dadd(&ref, 1, Amt, 0.5);
	assert_dmat_near(Sym, ref, TOL);
	dspcellsort(Sy);//sorting does not change the values
	dmat *Sym2=dcell2m(Sy);
	assert_dmat_near(Sym2, Sym, 0);
	dfree(Sym2); dfree(Sym); dfree(ref);
	//droptol removes values small relative to the maximum of each block
	long nz0=P(Sy, 1)->pp[NY(P(Sy, 1))];
	P(Sy, 1)->px[0]*=1e-20;
	dspcellscale(Sy, 3);
	dspcelldroptol(Sy, 1e-10);
	assert_int_equal(P(Sy, 1)->pp[NY(P(Sy, 1))], nz0-1);
	dspcellfree(Sy);
	//wdot: y'*A*x
	dcell *x=dcellnew(2, 1), *y=dcellnew(2, 1);
	for(int i=0; i<2; i++){
		P(x, i)=randmat(n[i], 1, &rstat);
		P(y, i)=randmat(n[i], 1, &rstat);
	}
	dmat *xm=dcell2m(x), *ym=dcell2m(y), *Ax=NULL;
	dmm(&Ax, 0, Am, xm, "nn", 1);
	assert_float_equal(dspcellwdot(y, A, x), ddot(ym, Ax), TOL);
	assert_float_equal(dspcellwdot(y, NULL, x), dcelldot(y, x), TOL);
	dfree(xm); dfree(ym); dfree(Ax);
	dcellfree(x); dcellfree(y);
	dfree(Amt); dfree(Am2);
	dspcellfree(At);
	dfree(Am);
	dspcellfree(A);
}

/** maprot rotates matrices and cells of matrices. */
static void cell_maprot(void **state){
	(void)state;
	dmat *A=dnew(9, 9);
	P(A, 6, 4)=1;//point at (+2, 0) from center (4,4)
	dcell *Ac=dcellnew(2, 1);
	P(Ac, 0)=ddup(A);
	P(Ac, 1)=ddup(A);
	dmaprot(A, M_PI/2);
	assert_float_equal(P(A, 4, 6), 1, 1e-10);//rotated to (0,+2)
	assert_float_equal(dsum(A), 1, 1e-10);
	dmaprot(Ac, M_PI);
	assert_float_equal(P(P(Ac, 1), 2, 4), 1, 1e-10);
	dmaprot(Ac, 0);//no-op
	assert_float_equal(P(P(Ac, 0), 2, 4), 1, 1e-10);
	dfree(A);
	dcellfree(Ac);
}

/*Write each object to the write end of a pipe and read back from the read end with readbin_fp.*/
static cell *pipe_roundtrip(const cell **A, int n, M_ID id, int level){
	int fds[2];
	if(pipe(fds)) return NULL;
	file_t *fw=zfdopen_as(fds[1], "w", ".bin");
	for(int i=0; i<n; i++){
		writedata(fw, A[i], 0);//small enough to fit in the pipe buffer
	}
	zfclose(fw);
	file_t *fr=zfdopen_as(fds[0], "r", ".bin");
	cell *B=readbin_fp(fr, id, level);
	zfclose(fr);
	return B;
}
/** writesock/readsock and readbin_fp with scanning, single element and type conversion. */
static void cell_io(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 12);
	//all element types through a socket pair
	cell *A=cellnew(4, 2);
	zmat *z=znew(2, 3); zrandn(z, 1, &rstat); P(A, 0)=(cell *)z;
	ssp *ss=sspnewrandu(4, 3, 1, 0.5, &rstat); P(A, 1)=(cell *)ss;
	csp *cs=cspnewrandu(3, 3, COMPLEX(1., 2.), 0.6, &rstat); P(A, 2)=(cell *)cs;
	zsp *zs=zspnewrandu(2, 4, FCOMPLEX(1.f, -1.f), 0.6, &rstat); P(A, 3)=(cell *)zs;
	lmat *l=lnew(3, 1); P(l, 2)=-5; P(A, 4)=(cell *)l;
	dmat *d=randmat(2, 2, &rstat); P(A, 5)=(cell *)d;
	//A(6) is empty
	A->keywords=strdup("hello=1");
	int sv[2];
	assert_int_equal(socketpair(AF_UNIX, SOCK_STREAM, 0, sv), 0);
	writesock(A, sv[0]);
	close(sv[0]);//end of data
	cell *B=readsock(sv[1]);
	close(sv[1]);
	assert_non_null(B);
	assert_true(iscell(B));
	assert_int_equal(NX(B), 4);
	assert_int_equal(NY(B), 2);
	assert_non_null(B->keywords);
	assert_non_null(strstr(B->keywords, "hello=1"));
	zmat *z2=(zmat *)P(B, 0);
	assert_int_equal(z2->id, M_ZMP);
	for(long i=0; i<PN(z); i++) assert_float_equal(cabsf(P(z2, i)-P(z, i)), 0, 0);
	smat *f1=NULL, *f2=NULL;
	sspfull(&f1, ss, 'n', 1);
	sspfull(&f2, (ssp *)P(B, 1), 'n', 1);
	for(long i=0; i<PN(f1); i++) assert_float_equal(P(f1, i), P(f2, i), 0);
	sfree(f1); sfree(f2);
	cmat *c1=NULL, *c2=NULL;
	cspfull(&c1, cs, 'n', 1);
	cspfull(&c2, (csp *)P(B, 2), 'n', 1);
	for(long i=0; i<PN(c1); i++) assert_float_equal(cabs(P(c1, i)-P(c2, i)), 0, 0);
	cfree(c1); cfree(c2);
	zmat *zz1=NULL, *zz2=NULL;
	zspfull(&zz1, zs, 'n', 1);
	zspfull(&zz2, (zsp *)P(B, 3), 'n', 1);
	for(long i=0; i<PN(zz1); i++) assert_float_equal(cabsf(P(zz1, i)-P(zz2, i)), 0, 0);
	zfree(zz1); zfree(zz2);
	assert_int_equal(P((lmat *)P(B, 4), 2), -5);
	assert_dmat_near((dmat *)P(B, 5), d, 0);
	assert_true(isempty(P(B, 6)));
	cellfree(B);
	cellfree(A);

	//several arrays written one after another are read as a cell (scan)
	const int nw=12;//more than the initial scan buffer
	const cell *arr[12];
	for(int i=0; i<nw; i++){
		dmat *t=dnew(1, 1);
		P(t, 0)=i;
		arr[i]=(cell *)t;
	}
	B=pipe_roundtrip(arr, nw, M_0, 1);
	assert_non_null(B);
	assert_int_equal(NX(B), nw);
	for(int i=0; i<nw; i++) assert_float_equal(P((dmat *)P(B, i), 0), i, 0);
	cellfree(B);
	//automatic mode with a single array is not wrapped
	B=pipe_roundtrip(arr, 1, M_0, -1);
	assert_false(iscell(B));
	assert_float_equal(P((dmat *)B, 0), 0, 0);
	cellfree(B);
	//request a matrix but find a cell: the first element is used
	dcell *C=dcellnew(2, 1);
	P(C, 0)=dref((dmat *)arr[3]);
	P(C, 1)=dref((dmat *)arr[4]);
	const cell *one[1]={(cell *)C};
	B=pipe_roundtrip(one, 1, M_0, 0);
	assert_false(iscell(B));
	assert_float_equal(P((dmat *)B, 0), 3, 0);
	cellfree(B);
	//type conversion on read: double to float
	B=pipe_roundtrip(one, 1, M_FLT, 1);
	assert_int_equal(P(B, 1)->id, M_FLT);
	assert_float_equal(P((smat *)P(B, 1), 0), 4, 0);
	cellfree(B);
	//unsupported requested type
	B=pipe_roundtrip(arr, 1, M_INT32, 0);
	assert_null(B);
	//empty cell read as matrix is NULL
	dcell *E=dcellnew(0, 0);
	const cell *empty[1]={(cell *)E};
	B=pipe_roundtrip(empty, 1, M_0, 0);
	assert_null(B);
	dcellfree(E);
	dcellfree(C);
	for(int i=0; i<nw; i++){
		cell *t=(cell *)arr[i];
		cellfree(t);
	}
	//invalid combinations are refused
	int fds[2];
	assert_int_equal(pipe(fds), 0);
	file_t *fw=zfdopen_as(fds[1], "w", ".bin");
	dmat *x=dnew(2, 2);
	writedata(fw, x, 2);//fp and ncol together: error, nothing written
	writecell_async(x, 0);//ncol=0: nothing to do
	zfclose(fw);
	char buf[8];
	assert_int_equal(read(fds[0], buf, sizeof(buf)), 0);
	close(fds[0]);
	dfree(x);
}

/*Write a 3x2 sparse matrix with 32 bit index in format magic to a pipe: [1 0; 0 2; 3 0]*/
static file_t *pipe_sp32(uint32_t magic){
	int fds[2];
	if(pipe(fds)) return NULL;
	file_t *fw=zfdopen_as(fds[1], "w", ".bin");
	header_t header={magic, 3, 2, NULL};
	write_header(&header, fw);
	uint64_t nzmax=3;
	int pp[3]={0, 2, 3};
	int pi[3]={0, 2, 1};
	zfwrite(&nzmax, sizeof(uint64_t), 1, fw);
	zfwrite(pp, sizeof(int), 3, fw);
	zfwrite(pi, sizeof(int), 3, fw);
	const double v[3]={1, 3, 2};
	switch(magic){
	case M_DSP32: zfwrite(v, sizeof(double), 3, fw); break;
	case M_SSP32: {float vf[3]={1, 3, 2}; zfwrite(vf, sizeof(float), 3, fw);} break;
	case M_CSP32: {dcomplex vc[3]={COMPLEX(1., 0.), COMPLEX(3., 0.), COMPLEX(2., 0.)}; zfwrite(vc, sizeof(dcomplex), 3, fw);} break;
	case M_ZSP32: {fcomplex vz[3]={FCOMPLEX(1.f, 0.f), FCOMPLEX(3.f, 0.f), FCOMPLEX(2.f, 0.f)}; zfwrite(vz, sizeof(fcomplex), 3, fw);} break;
	}
	zfclose(fw);
	return zfdopen_as(fds[0], "r", ".bin");
}
/** Sparse and dense readdata with header read from file and index/precision conversion. */
static void cell_io_convert(void **state){
	(void)state;
	const uint32_t magics[4]={M_DSP32, M_SSP32, M_CSP32, M_ZSP32};
	for(int im=0; im<4; im++){
		file_t *fr=pipe_sp32(magics[im]);
		dmat *full=NULL;
		if(im<2){
			dsp *A=dspreaddata(fr, NULL);
			dspfull(&full, A, 'n', 1);
			dspfree(A);
		}else{
			csp *A=cspreaddata(fr, NULL);
			cmat *cf=NULL;
			cspfull(&cf, A, 'n', 1);
			creal2d(&full, 0, cf, 1);
			cfree(cf);
			cspfree(A);
		}
		zfclose(fr);
		assert_non_null(full);
		assert_float_equal(P(full, 0, 0), 1, 0);
		assert_float_equal(P(full, 2, 0), 3, 0);
		assert_float_equal(P(full, 1, 1), 2, 0);
		assert_float_equal(dsum(full), 6, 0);
		dfree(full);
	}
	//single precision sparse read as double and complex single read as complex double
	rand_t rstat;
	seed_rand(&rstat, 13);
	ssp *ss=sspnewrandu(4, 4, 1, 0.5, &rstat);
	zsp *zs=zspnewrandu(3, 3, FCOMPLEX(1.f, 1.f), 0.5, &rstat);
	csp *cs=cspnewrandu(3, 3, COMPLEX(1., 1.), 0.5, &rstat);
	const cell *arr[1]={(cell *)ss};
	cell *B=pipe_roundtrip(arr, 1, M_DSP, 0);
	assert_int_equal(B->id, M_DSP);
	dmat *bf=NULL;
	dspfull(&bf, (dsp *)B, 'n', 1);
	smat *sf=NULL;
	sspfull(&sf, ss, 'n', 1);
	for(long i=0; i<PN(sf); i++) assert_float_equal(P(bf, i), P(sf, i), 0);
	dfree(bf); sfree(sf); cellfree(B);
	arr[0]=(cell *)zs;
	B=pipe_roundtrip(arr, 1, M_CSP, 0);
	assert_int_equal(B->id, M_CSP);
	cmat *cf=NULL;
	cspfull(&cf, (csp *)B, 'n', 1);
	zmat *zf=NULL;
	zspfull(&zf, zs, 'n', 1);
	for(long i=0; i<PN(zf); i++) assert_float_equal(cabs(P(cf, i)-(comp)P(zf, i)), 0, 0);
	cfree(cf); zfree(zf); cellfree(B);
	arr[0]=(cell *)cs;
	B=pipe_roundtrip(arr, 1, M_ZSP, 0);
	assert_int_equal(B->id, M_ZSP);
	cellfree(B);
	//dense readdata with the header read from the file
	int fds[2];
	assert_int_equal(pipe(fds), 0);
	file_t *fw=zfdopen_as(fds[1], "w", ".bin");
	smat *x=snew(2, 3);
	for(int i=0; i<6; i++) P(x, i)=(float)i*0.5f;
	swritedata(fw, x, 0);
	zfclose(fw);
	file_t *fr=zfdopen_as(fds[0], "r", ".bin");
	dmat *xd=dreaddata(fr, NULL);
	zfclose(fr);
	assert_int_equal(NX(xd), 2);
	assert_int_equal(NY(xd), 3);
	for(int i=0; i<6; i++) assert_float_equal(P(xd, i), P(x, i), 0);
	dwritedata(NULL, NULL, -1);//invalid usage: ignored
	dfree(xd); sfree(x);
	sspfree(ss); zspfree(zs); cspfree(cs);
}

UNITTEST_MAIN(
	cmocka_unit_test(cell_basic),
	cmocka_unit_test(cell_ref_hash),
	cmocka_unit_test(cell_resize),
	cmocka_unit_test(cell_constructors),
	cmocka_unit_test(cell_reshape),
	cmocka_unit_test(cell_cat),
	cmocka_unit_test(cell_elementwise),
	cmocka_unit_test(cell_add),
	cmocka_unit_test(cell_addI),
	cmocka_unit_test(cell_mm),
	cmocka_unit_test(cell_mm_types),
	cmocka_unit_test(cell_inverse),
	cmocka_unit_test(cell_spcell),
	cmocka_unit_test(cell_maprot),
	cmocka_unit_test(cell_io),
	cmocka_unit_test(cell_io_convert),
)
