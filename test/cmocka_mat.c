/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>

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


#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdint.h>
#include <cmocka.h>
#include "../lib/aos.h"

static void mat_basic(void **state){
    (void)state;
    dmat *b=dnew(3, 3);
    rand_t rstat;
    seed_rand(&rstat, 1);
    drandn(b, 1, &rstat);
    assert_int_equal(dhash(b, 1), 0x190394e1);
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
    mem_t *pmem=mem_new(NULL);
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
    assert_int_equal(dresize(a, 6, 0), 0);//resize with dimension 0 perserves
    assert_int_equal(NX(a), 6);
    assert_int_equal(NY(a), 5);
    assert_int_equal(mem_nref(a->mem), 1);
    assert_int_equal(dresize(a, 6, 0), 0);//resize with dimension 0 perserves
    assert_int_equal(NX(a), 6);
    assert_int_equal(NY(a), 5);
    assert_int_equal(mem_nref(a->mem), 1);
    assert_int_equal(dresize(a, 0, 0), 0);//resize with dimension 0 perserves
    assert_int_equal(NX(a), 6);
    assert_int_equal(NY(a), 5);
    assert_int_equal(mem_nref(a->mem), 1);
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
    assert_int_equal(dresize((dmat*)ac, 6, 5), 0);//resize empty dimension cell to dmat is ok
    assert_int_equal(NX(ac), 6);
    assert_int_equal(NY(ac), 5);
    assert_int_equal(ac->id, M_REAL);
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
    
    dfree(b);
    assert_null(b);
}
void print_mat(dmat *out){
    for(long iy=0; iy<NY(out); iy++){
        for(long ix=0; ix<NX(out); ix++){
            print_message("%4.0f ", P(out, ix, iy));
        }
        print_message("\n");
    }
}
void set_mat(dmat *out){
    int ny2=NY(out)>>1;
    int nx2=NX(out)>>1;
    for(long iy=0; iy<NY(out); iy++){
        for(long ix=0; ix<NX(out); ix++){
            P(out, ix, iy)=(ix-nx2)+(iy-ny2)*10;
        }
    }
}
void mat_embed_n(int nin, int nout){
    dmat *in=dnew(nin, nin);
    dmat *out=dnew(nout, nout);
    set_mat(in);
    dembed(out, in, 0);
   //print_mat(in);
    //print_mat(out);
    int nin2=nin>>1;
    int nout2=nout>>1;
    for(int io=-1; io<2; io++){
        for(int jo=-1; jo<2; jo++){
            assert_float_equal(P(in, nin2+io, nin2+jo), P(out, nout2+io, nout2+jo), 1e-16);
        }
    }
    dzero(out);
    dembed(out, in, M_PI/2);
    //print_mat(in);
    //rint_mat(out);
    assert_float_equal(P(in, nin2+1, nin2), P(out, nout2, nout2+1), 1e-16);

}
//test the embedding routines
static void mat_embed(void **state){
    (void)state;
    mat_embed_n(4, 8);
    mat_embed_n(4, 7);
    mat_embed_n(5, 8);
    mat_embed_n(5, 7);
    mat_embed_n(8, 4);
    mat_embed_n(7, 4);
    mat_embed_n(8, 5);
    mat_embed_n(7, 5);
}
static void mat_fresnel_prop(void **state){
    (void)state;
    if(zfexist("wvf0")){
        cmat *wvf0=cread("wvf0");
        cmat *wvf1=0;
        real dxout=0;
        fresnel_prop(&wvf1, &dxout, wvf0, 1./64, 0.589e-6, 100, 1, 1);
    }
}

int main(void){
    register_signal_handler(dummy_signal_handler);//captures error().
    //register_malloc(_test_malloc, _test_calloc, _test_realloc, _test_free);
    (void)dummy_signal_handler;
    (void)mat_basic;
    (void)mat_embed;
    (void)mat_fresnel_prop;
#if 0
    mat_basic(NULL);
    mat_embed(NULL);
    mat_fresnel_prop(NULL);
#else
    LOG_LEVEL=-4;
    const struct CMUnitTest tests[]={
        cmocka_unit_test(mat_basic),
        cmocka_unit_test(mat_embed),
        cmocka_unit_test(mat_fresnel_prop),
    };
    //(void)tests;
    return cmocka_run_group_tests(tests, NULL, NULL);
#endif
}
