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

static void mat_basic(void** state){
    (void)state;
    dmat* b=NULL;
    {
        dmat* a=dnew(3, 3);
        assert_int_equal(mem_nref(a->mem), 1);
        b=dref(a);
        P(b, 2, 2)=1;
        assert_float_equal(P(a, 2, 2), 1, 0);
        dset(a, 2);
        assert_float_equal(P(b, 2, 2), 2, 0);
        assert_int_equal(mem_nref(b->mem), 2);
        dfree(a);
        assert_null(a);
        assert_int_equal(mem_nref(b->mem), 1);
    }
    
    {
        dmat* c=dref_reshape(b, 9, 1);
        assert_int_equal(mem_nref(b->mem), 2);
        assert_float_equal(P(c, 8), 2, 0);
        dfree(c);
    }
    {
        dmat* d=drefcols(b, 1, 2);
        assert_int_equal(b->nx, d->nx);
        assert_int_equal(d->ny, 2);
        assert_ptr_equal(&P(b, 0, 1), &P(d, 0, 0));
        dresize(d, 2, 3);
        dfree(d);
    }
    {
        dmat* e=dsub(b, 1, 2, 2, 1);
        assert_int_equal(mem_nref(e->mem), 1);
        assert_int_equal(e->nx, 2);
        assert_int_equal(e->ny, 1);
        assert_int_equal(P(b, 1, 2), P(e,0, 0));
        dfree(e);
        e=dsub(b, 1, 2, 2, 2);
        assert_null(e);
    }
    {// test dcat
        dmat* f=ddup(b);
        assert_int_equal(mem_nref(f->mem), 1);
        assert_int_equal(mem_nref(b->mem), 1);
        dmat* f2=dcat(b, f, 0);
        assert_null(f2);
        f2=dcat(b, f, 1);
        assert_int_equal(f2->nx, 6);
        assert_int_equal(f2->ny, 3);
        dfree(f2);
        f2=dcat(b, f, 2);
        assert_int_equal(f2->nx, 3);
        assert_int_equal(f2->ny, 6);
        dfree(f2);
        dmat* fv=dref_reshape(f, 6, 1);
        f2=dcat(b, fv, 1);
        assert_null(f2);
        f2=dcat(b, fv, 2);
        assert_null(f2);
        dfree(fv);
        dfree(f2);
    }
    {
        dmat* c=ddup(b);
        dzerocol(c, 1);
        assert_float_equal(P(c, 0, 1), 0, 0); 
        assert_float_equal(P(c, 1, 1), 0, 0);
        assert_float_equal(P(c, 2, 1), 0, 0);
        assert_float_equal(dsum(c), 12, 0);
        
        dmat* d=dtrans(c);
        dshow(c,"c");
        assert_float_equal(P(d, 1, 0), 0, 0);
        assert_float_equal(P(d, 1, 1), 0, 0);
        assert_float_equal(P(d, 1, 2), 0, 0);
        assert_float_equal(dsum(d), 12, 0);
        
        dsort(d, 1);
        assert_float_equal(P(d, 0, 0), 0, 0);
        assert_float_equal(P(d, 2, 0), 2, 0);
        dsort(d, 0);
        assert_float_equal(P(d, 0, 0), 2, 0);
        assert_float_equal(P(d, 2, 0), 0, 0);
        
        dzero(c);
        assert_float_equal(dsum(c), 0, 0);
        assert_float_equal(dsum(b), 18, 0);
        dcp(&c, b);
        assert_float_equal(dsum(c), 18, 0);
        assert_float_equal(dtrace(c), 6, 0);
    }
    {
        dmat *c=ddup(b);
        dflip(c,0);
        assert_float_equal(P(c,2,2),P(b,0,0), 0);
        dfree(c);
    }
    {
        uint32_t d=dhash(b, 1);
        assert_int_equal(d, 2756278489);
    }
    dfree(b);
    assert_null(b);
}
void mat_embed_n(int nin, int nout){
    dmat *in=dnew(nin, nin);
    dmat *out=dnew(nout, nout);
    int nin2=nin>>1;
    int nout2=nout>>1;
    P(in, nin2)=1;
    P(in, nin2+1)=2;
    dembed(out, in, 0);
    for(long ix=0; ix<nout*nout; ix++){
        print_message("%g ", P(out,ix));
    }
    assert_float_equal(P(in, nin2), P(out, nout2), 1e-16);
}
//test the embedding routines
static void mat_embed(void** state){
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
static void mat_fresnel_prop(void** state){
    (void)state;
    if(zfexist("wvf0")){
        cmat *wvf0=cread("wvf0");
        cmat *wvf1=0;
        real dxout=0;
        fresnel_prop(&wvf1, &dxout, wvf0, 1./64, 0.589e-6, 100, 1, 1);
    }
}
static int dummy_signal_handler(int sig){
    info("Signal=%d caught, will ignore.\n", sig);
    return 1;
}
int main(void){
    register_signal_handler(dummy_signal_handler);
    LOG_LEVEL=-4;
    const struct CMUnitTest tests[]={
        cmocka_unit_test(mat_basic),
        cmocka_unit_test(mat_embed),
        cmocka_unit_test(mat_fresnel_prop),
    };
    (void)tests;
    return cmocka_run_group_tests(tests, NULL, NULL);
}
