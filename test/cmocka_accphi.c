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
  Tests of ray tracing (prop) between map/loc/pts/locstat inputs and outputs and
  the ray tracing matrices (mkh, mkht, mkh_cubic). Ported from test_accphi.c.
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 2e-5
#else
#define TOL 1e-10
#endif

/*Common geometry: a screen at height ht seen from a source at hs (cone effect
  scale) with output on a square grid of subapertures at the ground.*/
typedef struct{
	real dx;      //sampling of screen and output
	real ht;      //screen height
	real hs;      //source range
	real scale;   //1-ht/hs
	real shiftx, shifty;
	map_t *screen;//input screen
	loc_t *locin; //loc of the screen grid points (same ordering as screen)
	pts_t *pts;   //output as pts
	loc_t *locout;//output as loc, with stat and map
	loc_t *locout2;//output as loc, without stat
	map_t *mapout;//output as map (covers locout->map)
}geom_t;

static void geom_init(geom_t *g){
	memset(g, 0, sizeof(*g));
	const real D=4;   //output diameter
	const real D2=6;  //screen size
	const real dsa=0.5;
	g->dx=1./16.;
	g->scale=0.7;
	g->ht=10000;
	g->hs=g->ht/(1-g->scale);
	g->shiftx=0.01;
	g->shifty=0.05;
	long nx=(long)(D2/g->dx);
	g->screen=mapnew(nx, nx, g->dx, g->dx);
	g->screen->ht=g->ht;
	g->locin=mksqloc(nx, nx, g->dx, g->dx, g->screen->ox, g->screen->oy);
	g->locin->ht=g->ht;
	long nsa1=(long)(D/dsa);
	loc_t *saloc=mksqloc_auto(nsa1, nsa1, dsa, dsa);
	long nxsa=(long)(dsa/g->dx);
	g->pts=ptsnew(saloc->nloc, dsa, dsa, nxsa, nxsa, g->dx, g->dx);
	for(long isa=0; isa<saloc->nloc; isa++){
		g->pts->origx[isa]=saloc->locx[isa];
		g->pts->origy[isa]=saloc->locy[isa];
	}
	locfree(saloc);
	g->locout=pts2loc(g->pts);
	loc_create_stat(g->locout);
	loc_create_map(g->locout);
	g->locout2=pts2loc(g->pts);
	g->mapout=mapnew2(g->locout->map);
}
static void geom_free(geom_t *g){
	mapfree(g->screen);
	locfree(g->locin);
	ptsfree(g->pts);
	locfree(g->locout);
	locfree(g->locout2);
	mapfree(g->mapout);
}
/*Fill the screen with a smooth periodic function.*/
static void fill_sin(map_t *screen){
	const long nx=NX(screen), ny=NY(screen);
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			P(screen, ix, iy)=sin((real)ix/nx*2*M_PI)*sin((real)iy/ny*2*M_PI)+0.3*cos((real)(ix+2*iy)/nx*2*M_PI);
		}
	}
}
/*Fill the screen with a*x+b*y+c.*/
static void fill_linear(map_t *screen, real a, real b, real c){
	for(long iy=0; iy<NY(screen); iy++){
		for(long ix=0; ix<NX(screen); ix++){
			P(screen, ix, iy)=a*(screen->ox+ix*screen->dx)+b*(screen->oy+iy*screen->dy)+c;
		}
	}
}
/*propdata with common geometry filled in.*/
static propdata_t pd_base(const geom_t *g){
	propdata_t pd;
	memset(&pd, 0, sizeof(pd));
	pd.alpha=1;
	pd.hs=g->hs;
	pd.shiftx=g->shiftx;
	pd.shifty=g->shifty;
	return pd;
}
/*Propagate from the screen to all output types and return them in a dcell (all on locout ordering).
  out[0]: loc (prop_grid), out[1]: pts, out[2]: locstat, out[3]: map, out[4]: nongrid loc, out[5]: mkh*/
static dcell *prop_all(geom_t *g, real alpha, int wrap){
	const long nloc=g->locout->nloc;
	dcell *out=dcellnew_same(6, 1, nloc, 1);
	propdata_t pd;
	//map to loc
	pd=pd_base(g); pd.mapin=g->screen; pd.locout=g->locout2; pd.phiout=P(P(out, 0)); pd.alpha=alpha; pd.wrap=wrap; prop(&pd);
	//map to pts
	pd=pd_base(g); pd.mapin=g->screen; pd.ptsout=g->pts; pd.phiout=P(P(out, 1)); pd.alpha=alpha; pd.wrap=wrap; prop(&pd);
	//map to locstat
	pd=pd_base(g); pd.mapin=g->screen; pd.ostat=g->locout->stat; pd.phiout=P(P(out, 2)); pd.alpha=alpha; pd.wrap=wrap; prop(&pd);
	//map to map
	dzero(g->mapout->dmat);
	pd=pd_base(g); pd.mapin=g->screen; pd.mapout=g->mapout; pd.alpha=alpha; pd.wrap=wrap; prop(&pd);
	loc_extract(P(out, 3), g->locout, g->mapout);
	//loc to loc.
	pd=pd_base(g); pd.locin=g->locin; pd.phiin=P(g->screen); pd.locout=g->locout; pd.phiout=P(P(out, 4)); pd.alpha=alpha; prop(&pd);
	//ray tracing matrix
	pd=pd_base(g); pd.locin=g->locin; pd.locout=g->locout;
	dsp *H=g->locin->iac?mkh_cubic(&pd, g->locin->iac):mkh(&pd);
	dspmv(P(out, 5), H, g->screen->dmat, 'n', alpha);
	dspfree(H);
	return out;
}
/*Same as prop_all but with cubic influence function.
  out[0]: loc (prop_grid_cubic), out[1]: pts, out[2]: locstat, out[3]: map, out[4]: nongrid loc, out[5]: mkh_cubic*/
static dcell *prop_all_cubic(geom_t *g, real alpha, real iac){
	g->screen->iac=iac;
	g->locin->iac=iac;
	dcell *out=prop_all(g, alpha, 0);
	g->screen->iac=0;
	g->locin->iac=0;
	return out;
}
/*Assert that all entries in a dcell agree with the first one.*/
static void assert_all_agree(const dcell *out, real tol){
	for(long i=1; i<PN(out); i++){
		assert_dmat_near(P(out, 0), P(out, i), tol);
	}
}

/** All linear ray tracing methods (map to loc/pts/locstat/map, loc to loc, mkh) agree. */
static void test_linear_agree(void **state){
	(void)state;
	geom_t g; geom_init(&g);
	fill_sin(g.screen);
	dcell *out=prop_all(&g, 1, 0);
	assert_all_agree(out, TOL);
	//the result is not trivial
	assert_true(dmaxabs(P(out, 0))>0.1);
	cellfree(out);
	geom_free(&g);
}
/** All cubic ray tracing methods agree and differ from linear interpolation by a small amount. */
static void test_cubic_agree(void **state){
	(void)state;
	geom_t g; geom_init(&g);
	fill_sin(g.screen);
	dcell *outc=prop_all_cubic(&g, 1, 0.3);
	assert_all_agree(outc, TOL);
	dcell *outl=prop_all(&g, 1, 0);
	real diff=ddiff(P(outc, 0), P(outl, 0));
	assert_true(diff>0);//cubic is different
	assert_dmat_near(P(outc, 0), P(outl, 0), 1e-2);//but close for a smooth screen
	cellfree(outc);
	cellfree(outl);
	geom_free(&g);
}
/** Output is accumulated: alpha=1 followed by alpha=-1 gives zero, results scale with alpha and alpha=0 defaults to 1. */
static void test_alpha(void **state){
	(void)state;
	geom_t g; geom_init(&g);
	fill_sin(g.screen);
	dcell *out1=prop_all(&g, 1, 0);
	dcell *out2=prop_all(&g, 2.5, 0);
	dcell *out0=prop_all(&g, 0, 0);
	for(long i=0; i<PN(out1); i++){
		dmat *tmp=ddup(P(out1, i));
		dscale(tmp, 2.5);
		assert_dmat_near(tmp, P(out2, i), TOL);
		dfree(tmp);
		if(i<5){//prop() treats alpha=0 as 1 (out0[5] is from dspmv)
			assert_dmat_near(P(out0, i), P(out1, i), TOL);
		}
	}
	//accumulation
	const long nloc=g.locout->nloc;
	dmat *phi=dnew(nloc, 1);
	propdata_t pd=pd_base(&g);
	pd.mapin=g.screen; pd.locout=g.locout2; pd.phiout=P(phi);
	prop(&pd);
	pd.alpha=-1;
	prop(&pd);
	assert_float_equal(dmaxabs(phi), 0, TOL);
	//pts and stat
	pd=pd_base(&g); pd.mapin=g.screen; pd.ptsout=g.pts; pd.phiout=P(phi); prop(&pd);
	pd.alpha=-1; prop(&pd);
	pd=pd_base(&g); pd.mapin=g.screen; pd.ostat=g.locout->stat; pd.phiout=P(phi); prop(&pd);
	pd.alpha=-1; prop(&pd);
	assert_float_equal(dmaxabs(phi), 0, TOL);
	//accumulate into existing value
	dset(phi, 3);
	pd=pd_base(&g); pd.mapin=g.screen; pd.locout=g.locout2; pd.phiout=P(phi); prop(&pd);
	dadds(phi, -3);
	assert_dmat_near(phi, P(out1, 0), TOL);
	dfree(phi);
	cellfree(out0);
	cellfree(out1);
	cellfree(out2);
	geom_free(&g);
}
/** A constant screen propagates to the constant times alpha for linear and cubic methods. */
static void test_constant(void **state){
	(void)state;
	geom_t g; geom_init(&g);
	dset(g.screen->dmat, 1.7);
	const long nloc=g.locout->nloc;
	dmat *expect=dnew(nloc, 1);
	dset(expect, 1.7*-0.6);
	dcell *out=prop_all(&g, -0.6, 0);
	for(long i=0; i<PN(out); i++){
		assert_dmat_near(P(out, i), expect, TOL);
	}
	cellfree(out);
	out=prop_all_cubic(&g, -0.6, 0.3);
	for(long i=0; i<PN(out); i++){
		assert_dmat_near(P(out, i), expect, TOL);
	}
	cellfree(out);
	//The rows of the ray tracing matrices sum to 1.
	propdata_t pd=pd_base(&g); pd.locin=g.locin; pd.locout=g.locout;
	dsp *H=mkh(&pd);
	dmat *rs=dspsum(H, 2);
	dmat *one=dnew(nloc, 1); dset(one, 1);
	assert_dmat_near(rs, one, TOL);
	dfree(rs); dspfree(H);
	H=mkh_cubic(&pd, 0.3);
	rs=dspsum(H, 2);
	assert_dmat_near(rs, one, TOL);
	dfree(rs); dspfree(H);
	dfree(one);
	dfree(expect);
	geom_free(&g);
}
/** A tilted plane is reproduced exactly with the cone scaling, shift and angular offset (ht*theta). */
static void test_linear_screen(void **state){
	(void)state;
	geom_t g; geom_init(&g);
	const real a=0.3, b=-0.7, c=0.2;
	fill_linear(g.screen, a, b, c);
	const real thetax=2e-5, thetay=-1e-5;
	const long nloc=g.locout->nloc;
	dmat *expect=dnew(nloc, 1);
	for(long i=0; i<nloc; i++){
		real x=g.locout->locx[i]*g.scale+g.ht*thetax+g.shiftx;
		real y=g.locout->locy[i]*g.scale+g.ht*thetay+g.shifty;
		P(expect, i)=a*x+b*y+c;
	}
	for(int icubic=0; icubic<2; icubic++){
		g.screen->iac=g.locin->iac=icubic?0.3:0;
		dmat *phi=dnew(nloc, 1);
		propdata_t pd=pd_base(&g);
		pd.thetax=thetax; pd.thetay=thetay;
		pd.mapin=g.screen; pd.locout=g.locout2; pd.phiout=P(phi);
		prop(&pd);
		assert_dmat_near(phi, expect, TOL);
		dzero(phi);
		pd=pd_base(&g);
		pd.thetax=thetax; pd.thetay=thetay;
		pd.locin=g.locin; pd.phiin=P(g.screen); pd.locout=g.locout; pd.phiout=P(phi);
		prop(&pd);
		assert_dmat_near(phi, expect, TOL);
		dfree(phi);
	}
	g.screen->iac=g.locin->iac=0;
	//no cone effect when hs is infinity (hs=0 defaults to infinity)
	dmat *phi=dnew(nloc, 1);
	propdata_t pd=pd_base(&g);
	pd.hs=0;
	pd.mapin=g.screen; pd.locout=g.locout2; pd.phiout=P(phi);
	prop(&pd);
	for(long i=0; i<nloc; i++){
		P(expect, i)=a*(g.locout->locx[i]+g.shiftx)+b*(g.locout->locy[i]+g.shifty)+c;
	}
	assert_dmat_near(phi, expect, TOL);
	//rotation for loc to loc (rotates the coordinate)
	const real rot=0.3;
	dzero(phi);
	pd=pd_base(&g);
	pd.locin=g.locin; pd.phiin=P(g.screen); pd.locout=g.locout; pd.phiout=P(phi); pd.rot=rot;
	prop(&pd);
	for(long i=0; i<nloc; i++){
		real x0=g.locout->locx[i], y0=g.locout->locy[i];
		real x=(x0*cos(rot)+y0*sin(rot))*g.scale+g.shiftx;
		real y=(-x0*sin(rot)+y0*cos(rot))*g.scale+g.shifty;
		P(expect, i)=a*x+b*y+c;
	}
	assert_dmat_near(phi, expect, TOL);
	dfree(phi);
	dfree(expect);
	geom_free(&g);
}
/** wrap=1 makes the screen periodic; wrap=0 leaves out-of-screen points untouched. */
static void test_wrap(void **state){
	(void)state;
	geom_t g; geom_init(&g);
	fill_sin(g.screen);
	dcell *out1=prop_all(&g, 1, 1);
	//shift by a full period
	g.shiftx+=NX(g.screen)*g.dx;
	g.shifty-=2*NY(g.screen)*g.dx;
	dcell *out2=prop_all(&g, 1, 1);
	for(long i=0; i<4; i++){//only grid based methods support wrap
		assert_dmat_near(P(out1, i), P(out2, i), TOL);
	}
	cellfree(out2);
	//partial shift: points that falls off the screen are wrapped to the other side.
	g.shiftx=g.shiftx-NX(g.screen)*g.dx+2.5;
	g.shifty=g.shifty+2*NY(g.screen)*g.dx;
	out2=prop_all(&g, 1, 1);
	for(long i=1; i<4; i++){
		assert_dmat_near(P(out2, 0), P(out2, i), TOL);
	}
	assert_true(ddiff(P(out1, 0), P(out2, 0))>0.1);
	cellfree(out2);
	//equivalent shift by -period is identical
	g.shiftx-=NX(g.screen)*g.dx;
	dcell *out3=prop_all(&g, 1, 1);
	out2=prop_all(&g, 1, 1);
	for(long i=0; i<4; i++){
		assert_dmat_near(P(out3, i), P(out2, i), TOL);
	}
	cellfree(out3);
	cellfree(out2);
	cellfree(out1);
	//Without wrap, points outside of the screen are not touched.
	g.shiftx=100;
	const long nloc=g.locout->nloc;
	dmat *phi=dnew(nloc, 1);
	dset(phi, 7);
	dmat *expect=ddup(phi);
	propdata_t pd=pd_base(&g); pd.mapin=g.screen; pd.locout=g.locout2; pd.phiout=P(phi); prop(&pd);
	pd=pd_base(&g); pd.mapin=g.screen; pd.ptsout=g.pts; pd.phiout=P(phi); prop(&pd);
	pd=pd_base(&g); pd.mapin=g.screen; pd.ostat=g.locout->stat; pd.phiout=P(phi); prop(&pd);
	pd=pd_base(&g); pd.locin=g.locin; pd.phiin=P(g.screen); pd.locout=g.locout; pd.phiout=P(phi); prop(&pd);
	assert_dmat_near(phi, expect, 0);
	dfree(phi);
	dfree(expect);
	geom_free(&g);
}
/** The transpose propagations and mkht are adjoints of the forward operators: y'(Hx)=(H'y)'x. */
static void test_transpose(void **state){
	(void)state;
	geom_t g; geom_init(&g);
	rand_t rstat; seed_rand(&rstat, 1);
	drandn(g.screen->dmat, 1, &rstat);
	const long nloc=g.locout->nloc;
	dmat *y=dnew(nloc, 1);
	drandn(y, 1, &rstat);
	//map to locstat
	dmat *hx=dnew(nloc, 1);
	propdata_t pd=pd_base(&g); pd.mapin=g.screen; pd.ostat=g.locout->stat; pd.phiout=P(hx); prop(&pd);
	map_t *hty=mapnew2(g.screen);
	pd=pd_base(&g); pd.mapin=hty; pd.ostat=g.locout->stat; pd.phiout=P(y); pd.transpose=1; prop(&pd);
	real d1=ddot(y, hx);
	real d2=ddot(hty->dmat, g.screen->dmat);
	assert_float_equal(d1, d2, TOL*fabs(d1));
	//map to map
	map_t *ymap=mapnew2(g.mapout);
	drandn(ymap->dmat, 1, &rstat);
	dzero(g.mapout->dmat);
	pd=pd_base(&g); pd.mapin=g.screen; pd.mapout=g.mapout; prop(&pd);
	dzero(hty->dmat);
	pd=pd_base(&g); pd.mapin=hty; pd.mapout=ymap; pd.transpose=1; prop(&pd);
	d1=ddot(ymap->dmat, g.mapout->dmat);
	d2=ddot(hty->dmat, g.screen->dmat);
	assert_float_equal(d1, d2, TOL*fabs(d1));
	//the transpose of the forward ray tracing matrix equals mkht.
	for(int icubic=0; icubic<2; icubic++){
		g.locin->iac=icubic?0.3:0;
		pd=pd_base(&g); pd.locin=g.locin; pd.locout=g.locout;
		dsp *H=mkh(&pd);
		dsp *Ht=mkht(&pd);
		assert_int_equal(NX(H), nloc);
		assert_int_equal(NY(H), g.locin->nloc);
		assert_int_equal(NX(Ht), g.locin->nloc);
		assert_int_equal(NY(Ht), nloc);
		dsp *Htt=dsptrans(Ht);
		dsp *Hd=dspadd2(H, 1, Htt, -1);
		dmat *sa=dspsumabs(Hd, 1);
		assert_float_equal(dmaxabs(sa), 0, TOL);
		dfree(sa);
		//and the matrix is consistent with prop
		dmat *hx2=dnew(nloc, 1);
		dspmv(hx2, H, g.screen->dmat, 'n', 1);
		dzero(hx);
		g.screen->iac=g.locin->iac;
		pd=pd_base(&g); pd.mapin=g.screen; pd.locout=g.locout; pd.phiout=P(hx); prop(&pd);
		g.screen->iac=0;
		assert_dmat_near(hx, hx2, TOL);
		dfree(hx2);
		dspfree(Hd);
		dspfree(Htt);
		dspfree(Ht);
		dspfree(H);
	}
	g.locin->iac=0;
	mapfree(ymap);
	mapfree(hty);
	dfree(hx);
	dfree(y);
	geom_free(&g);
}
/** A screen above the source is skipped: prop leaves the output untouched and mkh returns NULL. */
static void test_above_source(void **state){
	(void)state;
	geom_t g; geom_init(&g);
	fill_sin(g.screen);
	const long nloc=g.locout->nloc;
	dmat *phi=dnew(nloc, 1);
	propdata_t pd=pd_base(&g);
	pd.hs=g.ht/2;
	pd.mapin=g.screen; pd.locout=g.locout2; pd.phiout=P(phi);
	prop(&pd);
	assert_float_equal(dmaxabs(phi), 0, 0);
	pd=pd_base(&g);
	pd.hs=g.ht/2;
	pd.locin=g.locin; pd.locout=g.locout;
	dsp *H=mkh(&pd);
	assert_null(H);
	dfree(phi);
	geom_free(&g);
}

UNITTEST_MAIN(
	cmocka_unit_test(test_linear_agree),
	cmocka_unit_test(test_cubic_agree),
	cmocka_unit_test(test_alpha),
	cmocka_unit_test(test_constant),
	cmocka_unit_test(test_linear_screen),
	cmocka_unit_test(test_wrap),
	cmocka_unit_test(test_transpose),
	cmocka_unit_test(test_above_source)
)
