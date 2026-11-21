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

void mat_basic(void **state){
    (void)state;
    dmat *b=dnew(3, 3);
    rand_t rstat;
    seed_rand(&rstat, 1);
    drandn(b, 1, &rstat);
	if(sizeof(b->p[0])==sizeof(double)){
    	assert_int_equal(dhash(b, 1), 0x190394e1);
	}else{
		assert_int_equal(dhash(b, 1), 0x814f5609);
	}
    assert_float_equal(dsum(b), 3.748561, 0.000001);
    assert_float_equal(dtrace(b), 0.820539, 0.000001);
    assert_float_equal(dmaxabs(b), 1.852614, 1e-6);
    assert_float_equal(dsumabs(b), 6.988739, 1e-6);
    //test dnew(), dset(), dref(), dfree()
    assert_int_equal(mem_nref(b->mem), 1);
    dmat *a=dref(b);
    assert_int_equal(mem_nref(b->mem), 2);
    assert_ptr_equal(P(a), P(b));
    dfree(a);
    assert_null(a);
    assert_int_equal(mem_nref(b->mem), 1);
    mem_t *pmem=mem_ref(mem_new(NULL));
    assert_non_null(pmem);
    a=dnew_do(3, 3, NULL, pmem);//invalid usage.
    assert_ptr_not_equal(a->mem, pmem);
    mem_unref(&pmem);
    dfree(a);
    real *p=mycalloc(9, real);
    a=dnew_do(3, 3, p, NULL);//a is weak pointer to p
    //dset
    dset(a, 2);
    assert_float_equal(dsum(a), 2*PN(a), 0.000001);
    assert_null(a->mem);
    assert_int_equal(mem_nref(a->mem), 0);
    assert_non_null(P(a));
    assert_ptr_equal(P(a), p);
    assert_int_equal(dresize(a,4,4), -1);//resize weak pointer is invalid
    dfree(a);
    free(p);
    assert_null(dref(a));
    //dinit and dresize
    assert_int_equal(dinit(&a, 3, 3), 0);//create
    assert_non_null(a);
    assert_int_equal(dinit(&a, 3, 3), 0);//recreate does nothing
    assert_non_null(a);
    assert_int_equal(mem_nref(a->mem), 1);
    assert_int_equal(dinit(&a, 0, 0), -1);//must fail
    assert_int_equal(mem_nref(a->mem), 1);
    assert_int_equal(NX(a), 3);
    assert_int_equal(NY(a), 3);
    assert_int_equal(dinit(&a, 4, 0), -1);//must fail
    assert_int_equal(mem_nref(a->mem), 1);
    assert_int_equal(NX(a), 3);
    assert_int_equal(NY(a), 3);
    assert_int_equal(dresize(a, 4, 5), 0);//resize
    assert_int_equal(NX(a), 4);
    assert_int_equal(NY(a), 5);
    assert_int_equal(mem_nref(a->mem), 1);
    assert_int_equal(dresize(a, 6, -1), 0);//resize with dimension -1 perserves
    assert_int_equal(NX(a), 6);
    assert_int_equal(NY(a), 5);
    assert_int_equal(mem_nref(a->mem), 1);
    assert_int_equal(dresize(a, 6, -1), 0);//resize with dimension -1 perserves
    assert_int_equal(NX(a), 6);
    assert_int_equal(NY(a), 5);
    assert_int_equal(mem_nref(a->mem), 1);
    assert_int_equal(dresize(a, -1, -1), 0);//resize with dimension -1 perserves
    assert_int_equal(NX(a), 6);
    assert_int_equal(NY(a), 5);
	assert_int_equal(dresize(a, 0, 0), 0);//resize with both dimension 0 frees
    assert_int_equal(NX(a), 0);
    assert_int_equal(NY(a), 0);
	assert_null(P(a));
    assert_int_equal(mem_nref(a->mem), 0);
    dfree(a);
    assert_int_equal(dresize(a, 6, 0), 0);//resize NULL with dimension 0 is ok
    assert_null(a);
    assert_int_equal(dresize(a, 6, 5), -1);//resize NULL with valid dimension is not valid
    a=dnew(0,0);//create empty dimension data is valid
    assert_non_null(a);
    assert_null(P(a));
    assert_int_equal(dresize(a, 6, 5), 0);//resize empty dimension with valid dimension is ok
    assert_int_equal(NX(a), 6);
    assert_int_equal(NY(a), 5);
    assert_int_equal(mem_nref(a->mem), 1);
    dfree(a);
    cell* ac=cellnew(0,0);
    assert_int_equal(dresize((dmat*)ac, 6, 5), -1);//resize cell to dmat is not allowed due to size incompatibility.
    /*assert_int_equal(NX(ac), 6);
    assert_int_equal(NY(ac), 5);
    assert_int_equal(ac->id, M_REAL);*/
    cellfree(ac);
    ac=cellnew(1,1);
    assert_int_equal(dresize((dmat*)ac, 6, 5), -1);//resize cell to dmat is not valid
    cellfree(ac);

    a=dref(b);
    assert_int_equal(dresize(a,5,3),-1);//resize referenced array is invalid
    assert_int_equal(dresize(a, 3, 3), -1);//resize referenced array is invalid even with the same dimension
    dfree(a);
    a=dnew(3, 0);
    dinit(&a, 3, 3);
    assert_int_equal(PN(a), 9);
    dfree(a);

    //dref_reshape
    assert_non_null(b);
    dmat *c=dref_reshape(b, 9, 1);
    assert_int_equal(mem_nref(b->mem), 2);
    assert_ptr_equal(&P(c, 8), &P(b, 2, 2));
    dfree(c);
    assert_ptr_equal(dref_reshape(b, 6, 1), NULL);//should fail
    assert_ptr_equal(dref_reshape(b, 6, 0), NULL);//should fail
    assert_ptr_not_equal(c=dref_reshape(b, 9, 0), NULL);//should not fail
    assert_int_equal(NX(c), 9);
    assert_int_equal(NY(c), 1);
    assert_int_equal(mem_nref(b->mem), 2);
    dfree(c);
    assert_ptr_not_equal(c=dref_reshape(b, 0, 1), NULL);//should not fail
    assert_int_equal(NX(c), 9);
    assert_int_equal(NY(c), 1);
    dfree(c);

    //drefcols
    dmat *d=drefcols(b, 1, 2);
    assert_int_equal(b->nx, d->nx);
    assert_int_equal(d->ny, 2);
    assert_ptr_equal(&P(b, 0, 1), &P(d, 0, 0));
    reshape(d, 2, 3);
    assert_int_equal(NX(d), 2);
    assert_int_equal(NY(d), 3);
    dfree(d);

    //dsub
    dmat *e=dsub(b, 1, 2, 2, 1);//copies data
    assert_int_equal(mem_nref(e->mem), 1);
    assert_ptr_not_equal(&P(b, 1, 2), &P(e, 0, 0));
    assert_int_equal(e->nx, 2);
    assert_int_equal(e->ny, 1);
    assert_int_equal(P(b, 1, 2), P(e, 0, 0));
    dfree(e);
    e=dsub(b, 1, 2, 2, 2);//not enough data
    assert_null(e);
    d=dsub(b, 1, 0, 2, 0);
    assert_int_equal(NX(d), 2);
    assert_int_equal(NY(d), 1);
    assert_float_equal(P(d, 1), P(b, 2, 2), 0);
    dfree(d);
    // test dcat
    dmat *f=ddup(b);
    assert_int_equal(mem_nref(f->mem), 1);
    assert_int_equal(mem_nref(b->mem), 1);
    dmat *f2=dcat(b, f, 0);//invalid usage
    assert_null(f2);
    f2=dcat(b, f, 1);//along x dimension
    assert_int_equal(f2->nx, 6);
    assert_int_equal(f2->ny, 3);
    dfree(f2);
    f2=dcat(b, f, 2);//along y dimension
    assert_int_equal(f2->nx, 3);
    assert_int_equal(f2->ny, 6);
    dfree(f2);
    dmat *fv=dref_reshape(f, 0, 1);//into a vector

    f2=dcat(b, fv, 1);//mismatch
    assert_null(f2);
    f2=dcat(b, fv, 2);//mismatch
    assert_null(f2);
    dfree(fv);
    dfree(f);
    f2=dcat(b, NULL, 1);
    assert_float_equal(P(f2,0,0), P(b,0,0), 1e-6);
    dfree(f2);
    f2=dcat( NULL, b, 1);
    assert_float_equal(P(f2, 0, 0), P(b, 0, 0), 1e-6);
    dfree(f2);
    f2=dcat(NULL, NULL, 1);
    assert_null(f2);
    //zerocols
    dmat *g=ddup(b);
    assert_int_equal(dzerocol(g, 3), -1);
    assert_int_equal(dzerocol(g, 1), 0);
    assert_float_equal(P(g, 0, 1), 0, 0);
    assert_float_equal(P(g, 1, 1), 0, 0);
    assert_float_equal(P(g, 2, 1), 0, 0);
    dcp(&g, NULL);//copy null is zero
    assert_float_equal(dsum(g), 0, 0);
    //trans
    dmat *h=dtrans(g);
    dshow(g, "c");
    assert_float_equal(P(h, 1, 0), 0, 0);
    assert_float_equal(P(h, 1, 1), 0, 0);
    assert_float_equal(P(h, 1, 2), 0, 0);
    assert_float_equal(dsum(h), dsum(g), 0);

    dsort(h, 1);
    dmat *hc=drefcols(h, 0, 1);
    assert_float_equal(P(h, 0, 0), dmin(hc), 0);
    assert_float_equal(P(h, 2, 0), dmax(hc), 0);
    dsort(h, 0);
    assert_float_equal(P(h, 0, 0), dmax(hc), 0);
    assert_float_equal(P(h, 2, 0), dmin(hc), 0);
    dfree(hc);
    dfree(h);

    dzero(g);
    assert_float_equal(dsum(g), 0, 0);
    dfree(g);

    dmat *k=ddup(b);
    assert_int_equal(dflip(k, 0),0);//flip up-down (1) and left-right (2) or both (0)
    assert_int_equal(dflip(k, -1), -1);//invalid index
    assert_int_equal(dflip(NULL, 0),0);//empty array
    assert_int_equal(dflip(NULL, -1), 0);//empty array
    assert_float_equal(P(k, 2, 2), P(b, 0, 0), 0);
    dfree(k);

    //cpcorner2center
    dcpcorner2center(k, b);
    assert_null(k);
    dcp(&k,b);
    dcpcorner2center(k, b);
    assert_float_equal(P(k,1,1), P(b,0,0),1e-6);
    
    assert_int_equal(dshift(&k, b, 1, 1), 0);
    assert_float_equal(P(k, 1, 1), P(b, 0, 0), 1e-6);
    dfree(k);
    
    dfree(b);
    assert_null(b);
}

void mat_cell(void **state){
	(void)state;
    dmat *in=dlinspace(0,1,16);
    long dims[2]={4,12};
    dcell *in2vec=d2cellref(in, dims, 2);
    assert_non_null(in2vec->m);
    assert_int_equal(NX(in2vec),2);
    assert_int_equal(NX(P(in2vec,1)),dims[1]);
    assert_float_equal(dcellsum(in2vec), dsum(in), 1e-12);
    dcell *inc=dcellnew_same(2, 2, 2, 2);assert_non_null(inc->m);
    reshape(in, 4, 4);
    dcell *out=0;
    d2cell(&out, in, inc);
    assert_float_equal(dcellsum(out), dsum(in), 1e-12);
    assert_int_equal(NX(out), NX(inc));
    assert_int_equal(NY(out), NY(inc));
    for(long i=0; i<PN(out); i++){
        assert_int_equal(NX(out, i), 2);
        assert_int_equal(NY(out, i), 2);
    }
    assert_float_equal(P(P(out, 1, 0), 0, 1), P(in, 2, 1), 0);
    dcellfree(out);
    dcellfree(inc);
    dcellfree(in2vec);
    dfree(in);
}
/** dinit errors, zerocol, cols, mean and maxmin. */
void mat_init_cols(void **state){
	(void)state;
	assert_int_equal(dinit(NULL, 2, 2), -1);
	cmat *c=cnew(2, 2);
	dmat *cd=(dmat *)c;
	assert_int_equal(dinit(&cd, 2, 2), -1);//type mismatch
	cfree(c);
	assert_int_equal(dzerocol(NULL, 0), -1);
	//cols overrides a stack struct without reference counting
	dmat *A=dnew(3, 4);
	for(long i=0; i<PN(A); i++) P(A, i)=(real)i;
	dmat col;
	memset(&col, 0, sizeof(col));
	dmat *pc=&col;
	dcols(pc, A, 1, 2);
	assert_int_equal(NX(pc), 3);
	assert_int_equal(NY(pc), 2);
	assert_ptr_equal(P(pc), PCOL(A, 1));
	assert_int_equal(mem_nref(A->mem), 1);
	dcols(pc, A, 3, 1);//again: the old content is released
	assert_float_equal(P(pc, 0), 9, 0);
	dcols(NULL, A, 0, 1);//no-op
	//mean and maxmin
	assert_float_equal(dmean(A), 5.5, 1e-15);
	real mx=0, mn=0;
	dmaxmin(A, &mx, &mn);
	assert_float_equal(mx, 11, 0);
	assert_float_equal(mn, 0, 0);
	dmaxmin(NULL, &mx, &mn);
	assert_float_equal(mx, 11, 0);
	cmat *C=cnew(2, 1);
	P(C, 0)=COMPLEX(3., 4.); P(C, 1)=COMPLEX(0., 1.);
	cmaxmin(C, &mx, &mn);//magnitude
	assert_float_equal(mx, 5, 1e-15);
	assert_float_equal(mn, 1, 1e-15);
	assert_float_equal(creal(cmean(C)), 1.5, 1e-15);
	assert_float_equal(cimag(cmean(C)), 2.5, 1e-15);
	lmat *L=lnew(4, 1);
	P(L, 0)=3; P(L, 1)=-7; P(L, 2)=12; P(L, 3)=0;
	long lmx=0, lmn=0;
	lmaxmin(L, &lmx, &lmn);
	assert_int_equal(lmx, 12);
	assert_int_equal(lmn, -7);
	assert_int_equal(lmean(L), 2);
	assert_int_equal(lmaxabs(L), 12);
	lfree(L);
	cfree(C);
	dfree(A);
}
/** vecperm/vecpermi are inverse permutations; negative entries conjugate. */
void mat_perm(void **state){
	(void)state;
	dmat *in=dlinspace(10, 1, 5);
	long perm[5]={2, 4, 1, 3, -0};//entry 0 maps to in(0)
	dmat *out=dnew(5, 1);
	dvecperm(out, in, perm);
	for(int i=0; i<5; i++) assert_float_equal(P(out, i), P(in, perm[i]>0?perm[i]:-perm[i]), 0);
	dmat *back=dnew(5, 1);
	dvecpermi(back, out, perm);
	assert_dmat_near(back, in, 0);
	dmat *bad=dnew(4, 1);
	dvecperm(bad, in, perm);//mismatch: no-op
	assert_float_equal(dsum(bad), 0, 0);
	dfree(bad);
	//complex with conjugation
	cmat *cin=cnew(3, 1);
	for(int i=0; i<3; i++) P(cin, i)=COMPLEX((double)i, (double)(i+1));
	long cperm[3]={-2, 1, 0};
	cmat *cout=cnew(3, 1);
	cvecperm(cout, cin, cperm);
	assert_float_equal(creal(P(cout, 0)), 2, 0);
	assert_float_equal(cimag(P(cout, 0)), -3, 0);
	assert_float_equal(cimag(P(cout, 1)), 2, 0);
	assert_float_equal(cimag(P(cout, 2)), -1, 0);//0 is treated as negative
	cmat *cback=cnew(3, 1);
	cvecpermi(cback, cout, cperm);
	for(int i=0; i<3; i++) assert_float_equal(cabs(P(cback, i)-P(cin, i)), 0, 0);
	cfree(cin); cfree(cout); cfree(cback);
	dfree(in); dfree(out); dfree(back);
}
/** flip along each axis for several types. */
void mat_flip(void **state){
	(void)state;
	dmat *A=dnew(3, 4);
	for(long i=0; i<PN(A); i++) P(A, i)=(real)i;
	dmat *B=ddup(A);
	assert_int_equal(dflip(B, 1), 0);
	for(long iy=0; iy<4; iy++) for(long ix=0; ix<3; ix++) assert_float_equal(P(B, ix, iy), P(A, 2-ix, iy), 0);
	dcp(&B, A);
	assert_int_equal(dflip(B, 2), 0);
	for(long iy=0; iy<4; iy++) for(long ix=0; ix<3; ix++) assert_float_equal(P(B, ix, iy), P(A, ix, 3-iy), 0);
	dflip(B, 1);
	dflip(B, 1);
	dflip(B, 2);
	assert_dmat_near(B, A, 0);
	assert_int_equal(dflip(B, 3), -1);
	lmat *L=lnew(2, 3);
	for(long i=0; i<6; i++) P(L, i)=i;
	lflip(L, 2);
	assert_int_equal(P(L, 0, 0), 4);
	assert_int_equal(P(L, 1, 2), 1);
	lfree(L);
	zmat *Z=znew(3, 1);
	P(Z, 0)=FCOMPLEX(1.f, 1.f);
	zflip(Z, 1);
	assert_float_equal(crealf(P(Z, 2)), 1, 0);
	zfree(Z);
	dfree(A); dfree(B);
}
/** unique, sort for real and integer types. */
void mat_unique(void **state){
	(void)state;
	dmat *A=dnew(7, 1);
	real v[7]={3, 1, 3, -2, 1, 5, 3};
	memcpy(P(A), v, sizeof(v));
	dmat *u=dunique(A, 1);
	assert_int_equal(NX(u), 4);
	assert_float_equal(P(u, 0), -2, 0);
	assert_float_equal(P(u, 1), 1, 0);
	assert_float_equal(P(u, 2), 3, 0);
	assert_float_equal(P(u, 3), 5, 0);
	dfree(u);
	u=dunique(A, 0);
	assert_float_equal(P(u, 0), 5, 0);
	assert_float_equal(P(u, 3), -2, 0);
	dfree(u);
	lmat *L=lnew(5, 1);
	P(L, 0)=4; P(L, 1)=4; P(L, 2)=-1; P(L, 3)=9; P(L, 4)=-1;
	lmat *lu=lunique(L, 1);
	assert_int_equal(NX(lu), 3);
	assert_int_equal(P(lu, 0), -1);
	assert_int_equal(P(lu, 2), 9);
	lsort(L, 0);
	assert_int_equal(P(L, 0), 9);
	assert_int_equal(P(L, 4), -1);
	smat *S=snew(3, 2);
	for(int i=0; i<6; i++) P(S, i)=(float)((i*5)%6);
	ssort(S, 1);//each column
	assert_float_equal(P(S, 0, 0), 0, 0);
	assert_float_equal(P(S, 2, 0), 5, 0);
	assert_float_equal(P(S, 0, 1), 1, 0);
	sfree(S);
	lfree(lu); lfree(L);
	dfree(A);
}
/** fftshift, cpcorner2center and shift for even, odd and vector shapes. */
void mat_fftshift(void **state){
	(void)state;
	dmat *A=dnew(4, 4);
	for(long i=0; i<16; i++) P(A, i)=(real)i;
	dmat *B=ddup(A);
	dfftshift(B);
	for(long iy=0; iy<4; iy++) for(long ix=0; ix<4; ix++) assert_float_equal(P(B, ix, iy), P(A, (ix+2)%4, (iy+2)%4), 0);
	dfftshift(B);
	assert_dmat_near(B, A, 0);
	//vector
	cmat *v=cnew(6, 1);
	for(long i=0; i<6; i++) P(v, i)=COMPLEX((double)i, 0.);
	cfftshift(v);
	assert_float_equal(creal(P(v, 0)), 3, 0);
	assert_float_equal(creal(P(v, 5)), 2, 0);
	cfree(v);
	//odd size vector: the last element is not moved (with a warning)
	dmat *o=dlinspace(0, 1, 5);
	dfftshift(o);
	assert_float_equal(P(o, 0), 2, 0);
	assert_float_equal(P(o, 2), 0, 0);
	assert_float_equal(P(o, 4), 4, 0);
	dfree(o);
	//cpcorner2center into a larger array zero pads
	dmat *C=dnew(6, 6);
	dset(C, 9);
	assert_int_equal(dcpcorner2center(C, A), 0);
	assert_float_equal(P(C, 3, 3), P(A, 0, 0), 0);
	assert_float_equal(P(C, 2, 2), P(A, 3, 3), 0);
	assert_float_equal(P(C, 0, 0), 0, 0);
	assert_float_equal(dsum(C), dsum(A), 0);
	assert_int_equal(dcpcorner2center(C, C), -1);
	//shift by zero copies; shifting into itself is refused
	dmat *D=NULL;
	assert_int_equal(dshift(&D, A, 0, 4), 0);
	assert_dmat_near(D, A, 0);
	assert_int_equal(dshift(&D, A, -1, 1), 0);
	assert_float_equal(P(D, 3, 1), P(A, 0, 0), 0);
	dmat *Dr=dref(D);
	assert_int_equal(dshift(&Dr, D, 1, 0), -1);
	dfree(Dr);
	assert_int_equal(dshift(&D, NULL, 1, 1), -1);
	//show and show_diag only print
	dshow(A, "A");
	dshow_diag(A, "A");
	dmat *big=dnew(20, 20);
	dshow_diag(big, "big");
	dmat *row=dnew(30, 1);
	dshow(row, "row");
	dshow(NULL, "null");
	dshow_diag(NULL, "null");
	cshow_diag((cmat *)NULL, "null");
	dfree(row); dfree(big);
	dfree(D); dfree(C); dfree(A); dfree(B);
}
/** 2cell with an empty reference and cell_cast of non-cells. */
void mat_2cell_empty(void **state){
	(void)state;
	dcell *ref=dcellnew(1, 1);
	dcell *B=NULL;
	d2cell(&B, NULL, ref);
	assert_non_null(B);
	assert_int_equal(NX(B), 1);
	assert_null(P(B, 0));
	dcellfree(B);
	dcellfree(ref);
}

UNITTEST_MAIN(
	cmocka_unit_test(mat_basic),
	cmocka_unit_test(mat_cell),
	cmocka_unit_test(mat_init_cols),
	cmocka_unit_test(mat_perm),
	cmocka_unit_test(mat_flip),
	cmocka_unit_test(mat_unique),
	cmocka_unit_test(mat_fftshift),
	cmocka_unit_test(mat_2cell_empty),
)
