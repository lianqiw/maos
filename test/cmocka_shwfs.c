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
  Test the Shack-Hartmann WFS subaperture geometry, geometric gradients, and
  CoG (lib/shwfs.c) and the laser launch telescope (lib/llt.c).
*/
#include "unittest.h"

/** Subaperture grid for different orders, central obscuration, and misregistration. */
static void shwfs_saloc_test(void **state){
	(void)state;
	const real dsa=0.5;
	loc_t *saloc=shwfs_saloc(1, dsa, 0, 0, 0);//one subaperture centered on the aperture
	assert_int_equal(saloc->nloc, 1);
	assert_float_equal(saloc->locx[0], -0.25, 1e-15);
	assert_float_equal(saloc->locy[0], -0.25, 1e-15);
	locfree(saloc);

	saloc=shwfs_saloc(2, dsa, 0, 0, 0);//2x2 (quadrant) subapertures
	assert_int_equal(saloc->nloc, 4);
	assert_float_equal(dmin(saloc->dmat), -dsa, 1e-15);
	assert_float_equal(dmax(saloc->dmat), 0, 1e-15);
	locfree(saloc);

	saloc=shwfs_saloc(4, dsa, 0, 0, 0);//all 4x4 overlap with the circular aperture
	assert_int_equal(saloc->nloc, 16);
	assert_float_equal(saloc->dx, dsa, 0);
	real sumx=0, sumy=0;
	for(long i=0; i<saloc->nloc; i++){
		sumx+=saloc->locx[i]+dsa*0.5;
		sumy+=saloc->locy[i]+dsa*0.5;
	}
	assert_float_equal(sumx, 0, 1e-12);//symmetric about the center
	assert_float_equal(sumy, 0, 1e-12);
	locfree(saloc);

	saloc=shwfs_saloc(4, dsa, 2, 0, 0);//the central 2x2 are within the obscuration
	assert_int_equal(saloc->nloc, 12);
	for(long i=0; i<saloc->nloc; i++){
		assert_true(fabs(saloc->locx[i]+dsa*0.5)>dsa||fabs(saloc->locy[i]+dsa*0.5)>dsa);
	}
	locfree(saloc);

	saloc=shwfs_saloc(4, dsa, 0, 0.25, 0.1);//misregistration keeps more subapertures
	assert_true(saloc->nloc>16);
	locfree(saloc);
}

/** Points in each subaperture with the sampling adjusted to divide the subaperture size. */
static void shwfs_pts_test(void **state){
	(void)state;
	const real dsa=0.5;
	loc_t *saloc=shwfs_saloc(2, dsa, 0, 0, 0);
	pts_t *pts=shwfs_pts(saloc, dsa, 1./64.);
	assert_int_equal(pts->nsa, 4);
	assert_int_equal(pts->nxsa, 32);
	assert_int_equal(pts->nysa, 32);
	assert_float_equal(pts->dx, 1./64., 1e-15);
	assert_float_equal(pts->dsa, dsa, 0);
	for(long i=0; i<pts->nsa; i++){
		assert_float_equal(pts->origx[i], saloc->locx[i]+pts->dx*0.5, 1e-15);
		assert_float_equal(pts->origy[i], saloc->locy[i]+pts->dx*0.5, 1e-15);
	}
	ptsfree(pts);
	pts=shwfs_pts(saloc, dsa, 0.03);//adjusted to dsa/16
	assert_int_equal(pts->nxsa, 16);
	assert_float_equal(pts->dx, dsa/16, 1e-15);
	ptsfree(pts);
	locfree(saloc);
}

/** The lenslet fill factor map keeps the requested fraction of light. */
static void shwfs_safill_test(void **state){
	(void)state;
	const int nxsa=16;
	dmat *ampi=shwfs_safill(nxsa, 1);
	assert_float_equal(dsum(ampi), nxsa*nxsa, 1e-12);
	assert_float_equal(dmin(ampi), 1, 0);
	dfree(ampi);
	ampi=shwfs_safill(nxsa, 0.8);//the outermost ring is partially transmitting
	assert_float_equal(dsum(ampi), nxsa*nxsa*0.8, 1e-10);
	assert_true(P(ampi, 0, 0)>0 && P(ampi, 0, 0)<1);
	assert_float_equal(P(ampi, 1, 1), 1, 0);
	dfree(ampi);
	ampi=shwfs_safill(nxsa, 0.5);//two rings are off and the third is partially transmitting
	assert_float_equal(dsum(ampi), nxsa*nxsa*0.5, 1e-10);
	assert_float_equal(P(ampi, 1, 0), 0, 0);
	assert_float_equal(P(ampi, 1, 1), 0, 0);
	assert_true(P(ampi, 2, 2)>0 && P(ampi, 2, 2)<1);
	assert_float_equal(P(ampi, nxsa/2, nxsa/2), 1, 0);
	assert_true(dmax(ampi)<=1);
	dfree(ampi);
}

/** Subaperture area and its maximum or minimum over WFS. */
static void shwfs_saa_test(void **state){
	(void)state;
	dmat *amp=dnew(8, 1);
	real val[8]={1, 1, 1, 1, 0, 1, 0, 1};
	for(int i=0; i<8; i++) P(amp, i)=val[i];
	dmat *saa=shwfs_saa(amp, 4);
	assert_int_equal(NX(saa), 2);
	assert_float_equal(P(saa, 0), 1, 1e-15);
	assert_float_equal(P(saa, 1), 0.5, 1e-15);

	dcell *saas=dcellnew(2, 1);
	P(saas, 0)=dref(saa);
	P(saas, 1)=dnew(2, 1);
	P(P(saas, 1), 0)=0.2;
	P(P(saas, 1), 1)=0.9;
	dmat *smax=shwfs_saa_max(saas, 1);
	dmat *smin=shwfs_saa_max(saas, -1);
	assert_float_equal(P(smax, 0), 1, 0);
	assert_float_equal(P(smax, 1), 0.9, 0);
	assert_float_equal(P(smin, 0), 0.2, 0);
	assert_float_equal(P(smin, 1), 0.5, 0);
	dfree(smax);
	dfree(smin);
	dcellfree(saas);
	saas=dcellnew(1, 1);
	P(saas, 0)=dref(saa);
	smax=shwfs_saa_max(saas, 1);//a reference for a single WFS
	assert_ptr_equal(P(smax), P(saa));
	dfree(smax);
	dcellfree(saas);
	assert_null(shwfs_saa_max(NULL, 1));
	dfree(saa);
	dfree(amp);
}

/** Isolated subapertures are turned off. */
static void shwfs_sa_isolated_test(void **state){
	(void)state;
	const real dsa=0.5;
	loc_t *saloc=shwfs_saloc(4, dsa, 0, 0, 0);//4x4
	pts_t *pts=shwfs_pts(saloc, dsa, dsa/4);
	dmat *saa=dnew(pts->nsa, 1);
	dset(saa, 1);
	shwfs_sa_isolated(pts, saa, 0.5);//all valid: unchanged
	assert_float_equal(dsum(saa), pts->nsa, 0);
	assert_null(pts->map);//the map is freed

	dzero(saa);
	P(saa, 5)=1;//a single valid subaperture
	shwfs_sa_isolated(pts, saa, 0.5);
	assert_float_equal(dsum(saa), 0, 0);

	dzero(saa);
	P(saa, 5)=1;//two neighbors have only one valid neighbor each
	P(saa, 6)=1;
	shwfs_sa_isolated(pts, saa, 0.5);
	assert_float_equal(dsum(saa), 0, 0);

	dzero(saa);//a 2x2 block: each has three valid neighbors
	P(saa, 5)=1; P(saa, 6)=1; P(saa, 9)=1; P(saa, 10)=1;
	shwfs_sa_isolated(pts, saa, 0.5);
	assert_float_equal(dsum(saa), 4, 0);
	dfree(saa);
	ptsfree(pts);
	locfree(saloc);

	saloc=shwfs_saloc(2, dsa, 0, 0, 0);//nothing is done with 4 or fewer subapertures
	pts=shwfs_pts(saloc, dsa, dsa/4);
	saa=dnew(pts->nsa, 1);
	P(saa, 0)=1;
	shwfs_sa_isolated(pts, saa, 0.5);
	assert_float_equal(dsum(saa), 1, 0);
	dfree(saa);
	ptsfree(pts);
	locfree(saloc);
}

/** The geometric gradients of a tilted OPD equal the tilt in every subaperture for both the average gradient and the zernike tilt. */
static void shwfs_grad_test(void **state){
	(void)state;
	const real dsa=0.5;
	const real tx=1e-6, ty=-2e-6;//tilt in radian
	loc_t *saloc=shwfs_saloc(4, dsa, 0, 0, 0);
	pts_t *pts=shwfs_pts(saloc, dsa, dsa/8);
	loc_t *loc=pts2loc(pts);
	const long nsa=pts->nsa;
	const long npsa=pts->nxsa*pts->nysa;
	dcell *amp=dcellnew(1, 1);
	P(amp, 0)=dnew(loc->nloc, 1);
	dset(P(amp, 0), 1);
	dcell *saa=dcellnew(1, 1);
	P(saa, 0)=shwfs_saa(P(amp, 0), npsa);
	dmat *opd=dnew(loc->nloc, 1);
	for(long i=0; i<loc->nloc; i++){
		P(opd, i)=tx*loc->locx[i]+ty*loc->locy[i];
	}
	dspcell *GS0=shwfs_mkg(loc, amp, saloc, saa, 0.3);
	assert_int_equal(NX(GS0), 1);
	dccell *saimcc=shwfs_saimcc(pts, amp);
	assert_int_equal(NX(saimcc), 1);
	assert_int_equal(NX(P(saimcc, 0)), nsa);
	for(int iz=0; iz<2; iz++){
		dmat *grad=NULL;
		shwfs_grad_geom(&grad, pts, opd, P(amp, 0), iz?NULL:P(GS0, 0), iz?P(saimcc, 0):NULL);
		assert_int_equal(NX(grad), nsa*2);
		for(long isa=0; isa<nsa; isa++){
			assert_float_equal(P(grad, isa), tx, 1e-12);
			assert_float_equal(P(grad, isa+nsa), ty, 1e-12);
		}
		shwfs_grad_geom(&grad, pts, opd, P(amp, 0), iz?NULL:P(GS0, 0), iz?P(saimcc, 0):NULL);//accumulates
		assert_float_equal(P(grad, 0), 2*tx, 1e-12);
		dfree(grad);
	}
	cellfree(saimcc);
	dspcellfree(GS0);
	dfree(opd);
	dcellfree(saa);
	dcellfree(amp);
	locfree(loc);
	ptsfree(pts);
	locfree(saloc);
}

/** The LLT pupil, its Gaussian amplitude, focus, and tip/tilt/focus removal. */
static void llt_test(void **state){
	(void)state;
	const real lltd=0.3, dsa=0.5, dx=1./64.;
	pts_t *pts=llt_pts(lltd, dsa, dx);//as large as the WFS subaperture
	assert_int_equal(pts->nsa, 1);
	assert_int_equal(pts->nxsa, 32);
	assert_float_equal(pts->origx[0], (dx-dsa)*0.5, 1e-15);
	assert_float_equal(pts->origy[0], (dx-dsa)*0.5, 1e-15);
	ptsfree(pts);
	pts=llt_pts(1, dsa, dx);//larger LLT
	assert_int_equal(pts->nxsa, 64);
	ptsfree(pts);

	pts=llt_pts(lltd, dsa, dx);
	const long nxsa=pts->nxsa;
	dmat *amp=llt_amp(pts, lltd, 1, 0, 0);
	assert_int_equal(NX(amp), nxsa*nxsa);
	real sum2=ddot(amp, amp);
	assert_true(sum2>0.5 && sum2<1);//clipped part of a beam normalized to 1
	assert_float_equal(P(amp, 0), 0, 0);//corner is outside of the LLT
	assert_float_equal(P(amp, nxsa/2*nxsa+nxsa/2-1), P(amp, (nxsa/2-1)*nxsa+nxsa/2), 1e-15);//symmetric
	dmat *amp2=llt_amp(pts, lltd, 1, 0.05, 0);//misregistration clips the beam
	assert_true(ddot(amp2, amp2)<sum2);
	dfree(amp2);

	loc_t *loc=mksqloc(nxsa, nxsa, dx, dx, pts->origx[0], pts->origy[0]);
	dcell *ncpa=dcellnew(1, 1);
	P(ncpa, 0)=dnew(nxsa, nxsa);
	llt_add_focus(ncpa, loc, amp, lltd, 100);
	real var=0, piston=0;
	long count=0;
	for(long i=0; i<loc->nloc; i++){
		if(P(amp, i)>0){
			count++;
			var+=P(P(ncpa, 0), i)*P(P(ncpa, 0), i);
			piston+=P(P(ncpa, 0), i);
		}
	}
	assert_float_equal(piston/count, 0, 1e-20);
	assert_float_equal(sqrt(var/count), 100e-9, 1e-15);//rms focus in m

	dmat *pttf=NULL, *proj=NULL;
	llt_remove_ttf(ncpa, loc, amp, lltd, 2, &pttf, &proj);//focus is removed
	assert_int_equal(NY(pttf), 4);
	assert_int_equal(NX(proj), 4);
	real resmax=0;//residual within the beam
	for(long i=0; i<loc->nloc; i++){
		if(P(amp, i)>0 && fabs(P(P(ncpa, 0), i))>resmax) resmax=fabs(P(P(ncpa, 0), i));
	}
	assert_float_equal(resmax, 0, 1e-15);
	dfree(pttf);
	dfree(proj);
	dzero(P(ncpa, 0));
	llt_add_focus(ncpa, loc, amp, lltd, 100);
	llt_remove_ttf(ncpa, loc, amp, lltd, 1, NULL, NULL);//focus is kept
	resmax=0;
	for(long i=0; i<loc->nloc; i++){
		if(P(amp, i)>0 && fabs(P(P(ncpa, 0), i))>resmax) resmax=fabs(P(P(ncpa, 0), i));
	}
	assert_true(resmax>1e-8);
	dcellfree(ncpa);
	locfree(loc);
	dfree(amp);
	ptsfree(pts);
}

/** CoG mask, gradient offset, and noise equivalent angle. */
static void shwfs_cog_test(void **state){
	(void)state;
	const long nsa=2, npix=8;
	dcell *i0=dcellnew_same(nsa, 1, npix, npix);
	for(long isa=0; isa<nsa; isa++){
		for(long iy=0; iy<npix; iy++){
			for(long ix=0; ix<npix; ix++){//a spot shifted by isa pixels along x
				real rx=ix-3.5-isa, ry=iy-3.5;
				P(P(i0, isa), ix, iy)=100*exp(-(rx*rx+ry*ry)/2);
			}
		}
	}
	dcell *mask=shwfs_cog_mask(i0, 50);
	assert_int_equal(NX(mask), nsa);
	assert_float_equal(dmax(P(mask, 0)), 1, 0);
	assert_float_equal(P(P(mask, 0), 0), P(P(i0, 0), 0)/50, 1e-15);
	
	const real pixtheta=1e-6;
	dmat *gradoff=dnew(nsa*2, 1);
	shwfs_cog_gradoff(gradoff, P(i0), NULL, nsa, 0, 0, pixtheta, 2*pixtheta);
	for(long isa=0; isa<nsa; isa++){
		real g[2];
		dcog(g, P(i0, isa), 0, 0, 0, 0, 0, NULL);
		assert_float_equal(P(gradoff, isa), g[0]*pixtheta, 1e-20);
		assert_float_equal(P(gradoff, isa+nsa), g[1]*2*pixtheta, 1e-20);
	}
	assert_float_equal(P(gradoff, nsa+0), 0, 1e-18);//no shift along y
	assert_true(P(gradoff, 1)>P(gradoff, 0));//the second spot is shifted along x
	shwfs_cog_gradoff(gradoff, P(i0), P(mask), nsa, 0, 0, pixtheta, 2*pixtheta);//accumulates
	assert_true(P(gradoff, 1)>0);
	dfree(gradoff);

	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *nea=shwfs_cog_nea(P(i0), NULL, NULL, NULL, NULL, nsa, 0, 0, 200, &rstat, 10, 0, 3, pixtheta, pixtheta);
	assert_int_equal(NX(nea), nsa);
	assert_int_equal(NY(nea), 3);
	assert_true(P(nea, 0, 0)>0 && P(nea, 0, 1)>0);
	real srot[2]={M_PI/2, M_PI/2};
	seed_rand(&rstat, 1);//same noise, rotated by 90 degree: xx and yy are swapped
	dmat *nea2=shwfs_cog_nea(P(i0), NULL, NULL, NULL, srot, nsa, 0, 0, 200, &rstat, 10, 0, 3, pixtheta, pixtheta);
	for(long isa=0; isa<nsa; isa++){
		assert_float_equal(P(nea2, isa, 0), P(nea, isa, 1), 1e-12*P(nea, isa, 1));
		assert_float_equal(P(nea2, isa, 1), P(nea, isa, 0), 1e-12*P(nea, isa, 0));
	}
	dfree(nea);
	dfree(nea2);
	dcellfree(mask);
	dcellfree(i0);
}

UNITTEST_MAIN(
	cmocka_unit_test(shwfs_saloc_test),
	cmocka_unit_test(shwfs_pts_test),
	cmocka_unit_test(shwfs_safill_test),
	cmocka_unit_test(shwfs_saa_test),
	cmocka_unit_test(shwfs_sa_isolated_test),
	cmocka_unit_test(shwfs_grad_test),
	cmocka_unit_test(llt_test),
	cmocka_unit_test(shwfs_cog_test),
)
