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
  Tests of projection of a tilted rectangular grid surface (like M3) onto the
  pupil (proj_rect_grid, m3proj). Ported from test_proj.c, which read the M3
  surface and aperture from files. The surface is now synthesized in memory.
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 1e-5
#else
#define TOL 1e-10
#endif

/*Geometry from the old test (TMT M3).*/
static const real d_m3_f=20.;/*from m3 to focus */
static const real d_exitpupil_f=46.38661051;
static const real r_exitpupil=1.546220350;
static const real r_pupil=15;
#define D_EXITPUPIL_M3 (d_exitpupil_f-d_m3_f)

/*Create a rectangular map with different x/y spacing, centered at origin.*/
static rmap_t *rmap_make(long nx, long ny, real dx, real dy){
	rmap_t *map=rmap_convert(dnew(nx, ny));
	map->dx=dx;
	map->dy=dy;
	map->ox=-nx/2*dx;
	map->oy=-ny/2*dy;
	return map;
}
/*Fill map with a*X+b*Y+c*/
static void rmap_linear(rmap_t *map, real a, real b, real c){
	for(long iy=0; iy<NY(map); iy++){
		for(long ix=0; ix<NX(map); ix++){
			P(map, ix, iy)=a*(map->ox+ix*map->dx)+b*(map->oy+iy*map->dy)+c;
		}
	}
}
/*Expected obliquity factor and surface coordinate of a pupil point following the
  geometry of a surface tilted by thetax, thetay with the beam angle betax, betay.*/
static void proj_expect(real *cosang, real *X, real *Y, real x, real y, real thetax, real thetay,
	real ratiox, real ratioy, real hs, real ht, real betax, real betay){
	real offx=hs*betax;
	real offy=hs*betay;
	if(ratiox<0) offx=-offx;
	if(ratioy<0) offy=-offy;
	const real px=x*ratiox+offx;
	const real py=y*ratioy+offy;
	//surface normal and ray direction
	const real nx=-cos(thetax), ny=-cos(thetay), nz=-sqrt(1-nx*nx-ny*ny);
	*cosang=(px*nx+py*ny-hs*nz)/sqrt((px*px+py*py+hs*hs)*(nx*nx+ny*ny+nz*nz));
	//intersection of ray with the tilted surface by law of sines.
	const real alx=atan2(px, hs);
	const real aly=atan2(py, hs);
	*X=(hs-ht)*sin(alx)/sin(thetax-alx)-offx/sin(thetax);
	*Y=(hs-ht)*sin(aly)/sin(thetay-aly)-offy/sin(thetay);
}

/** A constant surface projects to the constant times the scaling and the obliquity factor. */
static void test_constant(void **state){
	(void)state;
	rmap_t *map=rmap_make(120, 110, 0.02, 0.015);
	dset(map->dmat, 1.3);
	loc_t *loc=mksqloc_auto(32, 32, 1, 1);
	dmat *phi=dnew(loc->nloc, 1);
	const real thetax=M_PI*0.75, thetay=M_PI*0.5, sc=-2;
	const real ratio=r_exitpupil/r_pupil;
	proj_rect_grid(P(phi), map, thetax, thetay, loc, -ratio, ratio, NULL, sc, d_exitpupil_f, D_EXITPUPIL_M3, 0, 0);
	for(long i=0; i<loc->nloc; i++){
		real cosang, X, Y;
		proj_expect(&cosang, &X, &Y, loc->locx[i], loc->locy[i], thetax, thetay, -ratio, ratio, d_exitpupil_f, D_EXITPUPIL_M3, 0, 0);
		assert_float_equal(P(phi, i), sc*1.3*cosang, TOL);
		//for 45 degree tilt the obliquity is close to cos(45)
		assert_float_equal(cosang, cos(M_PI/4), 0.03);
	}
	//Without tilting, obliquity factor is just the cosine of the ray angle.
	dzero(phi);
	proj_rect_grid(P(phi), map, M_PI/2, M_PI/2, loc, -ratio, ratio, NULL, 1, d_exitpupil_f, D_EXITPUPIL_M3, 0, 0);
	for(long i=0; i<loc->nloc; i++){
		real x=loc->locx[i]*ratio, y=loc->locy[i]*ratio;
		real cosang=d_exitpupil_f/sqrt(x*x+y*y+d_exitpupil_f*d_exitpupil_f);
		assert_float_equal(P(phi, i), 1.3*cosang, TOL);
	}
	dfree(phi);
	locfree(loc);
	rmapfree(map);
}
/** Without tilt, a linear surface projects exactly to a linear function scaled by (1-ht/hs) and offset by the beam angle. */
static void test_linear_untilted(void **state){
	(void)state;
	rmap_t *map=rmap_make(120, 110, 0.02, 0.015);
	const real a=0.4, b=-0.9, c=0.1;
	rmap_linear(map, a, b, c);
	loc_t *loc=mksqloc_auto(32, 32, 1, 1);
	dmat *phi=dnew(loc->nloc, 1);
	const real ratio=r_exitpupil/r_pupil;
	const real hs=d_exitpupil_f, ht=D_EXITPUPIL_M3;
	const real betax=1e-3, betay=-2e-3;
	proj_rect_grid(P(phi), map, M_PI/2, M_PI/2, loc, ratio, ratio, NULL, 1, hs, ht, betax, betay);
	for(long i=0; i<loc->nloc; i++){
		real x=loc->locx[i]*ratio, y=loc->locy[i]*ratio;
		real cosang=hs/sqrt(pow(x+hs*betax, 2)+pow(y+hs*betay, 2)+hs*hs);
		//similar triangles: surface coordinate of the ray
		real X=x*(1-ht/hs)-betax*ht;
		real Y=y*(1-ht/hs)-betay*ht;
		assert_float_equal(P(phi, i), cosang*(a*X+b*Y+c), TOL);
	}
	dfree(phi);
	locfree(loc);
	rmapfree(map);
}
/** With tilt, a linear surface is projected following the tilted surface geometry. */
static void test_linear_tilted(void **state){
	(void)state;
	rmap_t *map=rmap_make(120, 110, 0.02, 0.015);
	const real a=0.4, b=-0.9, c=0.1;
	rmap_linear(map, a, b, c);
	loc_t *loc=mksqloc_auto(32, 32, 1, 1);
	dmat *phi=dnew(loc->nloc, 1);
	const real ratio=r_exitpupil/r_pupil;
	const real hs=d_exitpupil_f, ht=D_EXITPUPIL_M3;
	const real betax=-1e-3, betay=2e-3;
	const real sc=-2;
	for(int itilt=0; itilt<2; itilt++){
		const real thetax=itilt?M_PI*0.5:M_PI*0.75;
		const real thetay=itilt?M_PI*0.6:M_PI*0.5;
		dzero(phi);
		proj_rect_grid(P(phi), map, thetax, thetay, loc, -ratio, ratio, NULL, sc, hs, ht, betax, betay);
		real Xmin=INFINITY, Xmax=-INFINITY;
		for(long i=0; i<loc->nloc; i++){
			real cosang, X, Y;
			proj_expect(&cosang, &X, &Y, loc->locx[i], loc->locy[i], thetax, thetay, -ratio, ratio, hs, ht, betax, betay);
			assert_float_equal(P(phi, i), sc*cosang*(a*X+b*Y+c), TOL);
			if(X<Xmin) Xmin=X;
			if(X>Xmax) Xmax=X;
		}
		//The footprint is stretched by 1/sin(theta) along the tilt direction.
		if(!itilt){
			real span0=(32-1)*ratio*(1-ht/hs);
			assert_true((Xmax-Xmin)>span0*1.3);
		}
	}
	dfree(phi);
	locfree(loc);
	rmapfree(map);
}
/** Output is accumulated, scales with sc, and points with zero amplitude or outside the surface are untouched. */
static void test_amp_and_bounds(void **state){
	(void)state;
	rmap_t *map=rmap_make(120, 110, 0.02, 0.015);
	rand_t rstat; seed_rand(&rstat, 2);
	drandn(map->dmat, 1, &rstat);
	loc_t *loc=mksqloc_auto(32, 32, 1, 1);
	const long nloc=loc->nloc;
	const real ratio=r_exitpupil/r_pupil;
	const real hs=d_exitpupil_f, ht=D_EXITPUPIL_M3;
	dmat *phi1=dnew(nloc, 1);
	proj_rect_grid(P(phi1), map, M_PI*0.75, M_PI*0.5, loc, -ratio, ratio, NULL, 1, hs, ht, 0, 0);
	//all points are inside the surface and get nonzero values.
	for(long i=0; i<nloc; i++){
		assert_true(P(phi1, i)!=0);
	}
	//accumulation with scaling
	dmat *phi2=ddup(phi1);
	proj_rect_grid(P(phi2), map, M_PI*0.75, M_PI*0.5, loc, -ratio, ratio, NULL, -3, hs, ht, 0, 0);
	dmat *expect=ddup(phi1);
	dscale(expect, -2);
	assert_dmat_near(phi2, expect, TOL);
	//amplitude mask: points with zero amplitude are skipped
	dmat *amp=dnew(nloc, 1);
	for(long i=0; i<nloc; i++){
		P(amp, i)=(RSS(loc->locx[i], loc->locy[i])<12)?1:0;
	}
	dset(phi2, 5);
	proj_rect_grid(P(phi2), map, M_PI*0.75, M_PI*0.5, loc, -ratio, ratio, P(amp), 1, hs, ht, 0, 0);
	for(long i=0; i<nloc; i++){
		assert_float_equal(P(phi2, i), P(amp, i)?(5+P(phi1, i)):5, TOL);
	}
	//A surface smaller than the beam footprint: points outside are untouched
	rmap_t *small=rmap_make(40, 40, 0.02, 0.015);
	dset(small->dmat, 1);
	dset(phi2, 5);
	proj_rect_grid(P(phi2), small, M_PI/2, M_PI/2, loc, ratio, ratio, NULL, 1, hs, ht, 0, 0);
	long ninside=0, noutside=0;
	for(long i=0; i<nloc; i++){
		real X=loc->locx[i]*ratio*(1-ht/hs);
		real Y=loc->locy[i]*ratio*(1-ht/hs);
		//strictly inside or outside with margin of one pixel
		int inside=X>small->ox+small->dx&&X<small->ox+(NX(small)-2)*small->dx
			&&Y>small->oy+small->dy&&Y<small->oy+(NY(small)-2)*small->dy;
		int outside=X<small->ox-small->dx||X>small->ox+NX(small)*small->dx
			||Y<small->oy-small->dy||Y>small->oy+NY(small)*small->dy;
		if(inside){
			ninside++;
			assert_float_equal(P(phi2, i), 5+1*hs/sqrt(pow(loc->locx[i]*ratio, 2)+pow(loc->locy[i]*ratio, 2)+hs*hs), TOL);
		}else if(outside){
			noutside++;
			assert_float_equal(P(phi2, i), 5, 0);
		}
	}
	assert_true(ninside>0);
	assert_true(noutside>0);
	//Too much tilt is rejected and output is untouched
	dset(phi2, 5);
	proj_rect_grid(P(phi2), map, 0.04, M_PI/2, loc, ratio, ratio, NULL, 1, hs, ht, 0, 0);
	assert_float_equal(dmaxabs(phi2), 5, 0);
	assert_float_equal(dmin(phi2), 5, 0);
	rmapfree(small);
	dfree(amp);
	dfree(expect);
	dfree(phi1);
	dfree(phi2);
	locfree(loc);
	rmapfree(map);
}
/** m3proj() with the M3 parameters in the map header equals the direct call to proj_rect_grid used by the old test. */
static void test_m3proj(void **state){
	(void)state;
	rmap_t *map=rmap_make(120, 110, 0.02, 0.015);
	rand_t rstat; seed_rand(&rstat, 3);
	drandn(map->dmat, 1, &rstat);
	map->txdeg=135;
	map->tydeg=90;
	map->fexit=d_exitpupil_f;
	map->fsurf=d_m3_f;
	map->ftel=d_exitpupil_f*r_pupil/r_exitpupil;
	loc_t *loc=mksqloc_auto(32, 32, 1, 1);
	const long nloc=loc->nloc;
	dmat *phi1=dnew(nloc, 1);
	dmat *phi2=dnew(nloc, 1);
	const real thx=1e-5, thy=-2e-5;//guide star direction
	m3proj(map, phi1, loc, thx, thy, INFINITY);
	const real bx=thx*map->ftel/map->fexit;
	const real by=thy*map->ftel/map->fexit;
	proj_rect_grid(P(phi2), map, M_PI*0.75, M_PI*0.5,
		loc, -r_exitpupil/r_pupil, r_exitpupil/r_pupil,
		NULL, -2, d_exitpupil_f, D_EXITPUPIL_M3, bx, by);
	assert_dmat_near(phi1, phi2, TOL);
	assert_true(dmaxabs(phi1)>0.1);
	//m3proj2 takes the parameters from the keywords string (rmap_make_keywords round trip)
	rmap_make_keywords((cell*)map);
	char *keywords=strdup(map->keywords);
	dmat *data=ddup(map->dmat);
	dmat *phi3=m3proj2(data, keywords, loc, thx, thy, INFINITY);
	assert_dmat_near(phi1, phi3, 1e3*TOL);
	dfree(phi3);
	dfree(data);
	free(keywords);
	dfree(phi1);
	dfree(phi2);
	locfree(loc);
	rmapfree(map);
}

UNITTEST_MAIN(
	cmocka_unit_test(test_constant),
	cmocka_unit_test(test_linear_untilted),
	cmocka_unit_test(test_linear_tilted),
	cmocka_unit_test(test_amp_and_bounds),
	cmocka_unit_test(test_m3proj)
)
