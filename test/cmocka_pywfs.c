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
  Tests of the pyramid wavefront sensor (lib/pywfs.c) using the library API
  with a small configuration. pywfs_mkg() is not tested since it caches the
  interaction matrix in a file; the interaction matrix is built here by
  poking modes with pywfs_ints() and pywfs_grad().
*/
static const real D=1;
static const int ORDER=8;
static const real WVL=0.8e-6;
static const real SIGLEV=100;
/** Configuration with the essential parameters; others take default values. */
static pywfs_cfg_t *mkcfg(int nside){
	pywfs_cfg_t *cfg=mycalloc(1, pywfs_cfg_t);
	cfg->order=ORDER;
	cfg->siglev=SIGLEV;
	cfg->wvl=dnew(1, 1);
	P(cfg->wvl, 0)=WVL;
	cfg->D=D;
	cfg->dx=D/(ORDER*4);
	cfg->nside=nside;
	cfg->modulpos=8;/*fewer positions to keep the test fast*/
	return cfg;
}
/** Gradients of OPD (on locfft->loc). */
static dmat *pygrad(const pywfs_t *py, const dmat *opd){
	dmat *ints=NULL, *grad=NULL;
	pywfs_ints(&ints, py, opd, py->cfg->siglev);
	pywfs_grad(&grad, py, ints);
	dfree(ints);
	return grad;
}
/** OPD of tilt (radian) on loc. */
static dmat *mktilt(const loc_t *loc, real tx, real ty){
	dmat *opd=dnew(loc->nloc, 1);
	real ptt[3]={0, tx, ty};
	loc_add_ptt(opd, ptt, loc);
	return opd;
}
static real vmean(const real *p, long n){
	real s=0;
	for(long i=0; i<n; i++) s+=p[i];
	return s/n;
}
/**
   pywfs_new() with defaults: dimensions, filled defaults, subaperture
   illumination normalized, optical gain calibrated so the T/T response
   averages to one; intensities are non-negative and sum to the signal level
   within the detected fraction; gradients of a flat wavefront are zero.
*/
static void pywfs_setup(void **state){
	(void)state;
	pywfs_cfg_t *cfg=mkcfg(0);
	cfg->modulpos=0;
	pywfs_t *py=pywfs_new(cfg, NULL, NULL, NULL);
	assert_non_null(py);
	/*defaults*/
	assert_int_equal(cfg->nside, 4);
	assert_int_equal(cfg->modulpos, 32);
	assert_int_equal(cfg->modulring, 1);
	assert_int_equal(cfg->sigmatch, 2);
	assert_float_equal(cfg->modulate, 5*WVL/D, 1e-20);
	assert_float_equal(cfg->dsa, D/ORDER, 1e-15);
	assert_float_equal(cfg->poke, 1e-7, 0);
	assert_float_equal(cfg->pixblur, 0.3, 0);
	assert_true(isinf(cfg->hs));
	assert_float_equal(P(cfg->wvlwts, 0), 1, 0);
	const long nsa=py->saloc->nloc;
	assert_int_equal(nsa, ORDER*ORDER);
	assert_int_equal(NX(py->si), 4);
	assert_int_equal(NX(P(py->si, 0)), nsa);
	assert_float_equal(dsum(py->saa), nsa, 1e-8);
	assert_true(py->gain>0);
	assert_int_equal(NX(py->GTT), nsa*2);
	assert_int_equal(NY(py->GTT), 2);
	/*the T/T response averages to one radian per radian*/
	assert_float_equal(dsum(py->GTT)/(2*nsa), 1, 1e-10);
	assert_float_equal(vmean(PCOL(py->GTT, 0), nsa), 1, 0.05);
	assert_float_equal(vmean(PCOL(py->GTT, 1)+nsa, nsa), 1, 0.05);
	assert_float_equal(vmean(PCOL(py->GTT, 0)+nsa, nsa), 0, 0.05);
	/*intensity of a flat wavefront*/
	const loc_t *loc=py->locfft->loc;
	dmat *opd=dnew(loc->nloc, 1);
	dmat *ints=NULL;
	pywfs_ints(&ints, py, opd, SIGLEV);
	assert_int_equal(NX(ints), nsa);
	assert_int_equal(NY(ints), 4);
	real isum=dsum(ints);
	assert_true(dmin(ints)>-1e-8*isum);
	assert_true(isum<=nsa*SIGLEV*1.0001&&isum>0.5*nsa*SIGLEV);
	/*the four pupil images are equal for a flat wavefront (by symmetry)*/
	real s0=dvecsum(PCOL(ints, 0), nsa);
	for(int i=1; i<4; i++){
		assert_float_equal(dvecsum(PCOL(ints, i), nsa), s0, 0.02*s0);
	}
	/*accumulates; scales with signal level*/
	pywfs_ints(&ints, py, opd, SIGLEV);
	assert_float_equal(dsum(ints), 2*isum, 1e-8*isum);
	/*gradient of a flat wavefront is zero after the reference subtraction*/
	dmat *grad=pygrad(py, opd);
	assert_int_equal(PN(grad), nsa*2);
	assert_float_equal(dmaxabs(grad), 0, 1e-12);
	/*gradoff is the gradient before subtraction: nonzero but small*/
	assert_non_null(py->gradoff);
	dfree(grad); dfree(ints); dfree(opd);
	pywfs_free(py);
	pycfg_free(cfg);
	/*missing required parameters*/
	pywfs_free(NULL);
	pycfg_free(NULL);
}
/**
   pywfs_new() with dim subapertures removed (saat>0) derives the intensity of
   the remaining subapertures from that of all subapertures: the gradient of a
   flat wavefront is still zero. With opdbias, gradoff is the gradient of
   opdbias, while the optical gain and subaperture areas are unchanged.
*/
static void pywfs_gradoff(void **state){
	(void)state;
	loc_t *loc=mkannloc(D, 0, D/(ORDER*4), 0);
	pywfs_cfg_t *cfg=mkcfg(4);
	cfg->saat=0.5;
	pywfs_t *py=pywfs_new(cfg, loc, NULL, NULL);
	const long nsa=py->saloc->nloc;
	assert_true(nsa<ORDER*ORDER);/*some subapertures are removed*/
	dmat *opd0=dnew(loc->nloc, 1);
	dmat *grad=pygrad(py, opd0);
	assert_float_equal(dmaxabs(grad), 0, 1e-12*dmaxabs(py->gradoff));
	dfree(grad);
	/*opdbias: astigmatism of 20 nm*/
	dmat *opd=dnew(loc->nloc, 1);
	for(long i=0; i<loc->nloc; i++){
		P(opd, i)=20e-9*(loc->locx[i]*loc->locx[i]-loc->locy[i]*loc->locy[i])/(D*D*0.25);
	}
	pywfs_cfg_t *cfg2=mkcfg(4);
	cfg2->saat=0.5;
	pywfs_t *py2=pywfs_new(cfg2, loc, NULL, opd);
	assert_int_equal(py2->saloc->nloc, nsa);
	assert_float_equal(py2->gain, py->gain, 0);
	assert_dmat_near(py2->saa, py->saa, 0);
	assert_dmat_near(py2->GTT, py->GTT, 0);
	/*the gradient of opdbias is zero after the reference subtraction*/
	grad=pygrad(py2, opd);
	assert_float_equal(dmaxabs(grad), 0, 1e-12*dmaxabs(py2->gradoff));
	/*and equals the difference of the offsets with the flat reference*/
	dmat *grad1=pygrad(py, opd);
	dmat *diff=ddup(py2->gradoff);
	dadd(&diff, 1, py->gradoff, -1);
	assert_true(dmaxabs(diff)>0);
	assert_dmat_near(grad1, diff, 1e-10*dmaxabs(diff));
	dfree(grad); dfree(grad1); dfree(diff); dfree(opd); dfree(opd0);
	pywfs_free(py); pywfs_free(py2);
	pycfg_free(cfg); pycfg_free(cfg2);
	locfree(loc);
}
/**
   Gradients respond linearly to a small tilt with unit gain in the tilt
   direction and no cross talk; reversing the tilt reverses the gradient.
   A modal interaction matrix built by poking Zernike modes has the expected
   dimensions and its tip/tilt columns agree with GTT.
*/
static void pywfs_linear(void **state){
	(void)state;
	pywfs_cfg_t *cfg=mkcfg(4);
	pywfs_t *py=pywfs_new(cfg, NULL, NULL, NULL);
	const loc_t *loc=py->locfft->loc;
	const long nsa=py->saloc->nloc;
	const real a=0.002*AS2RAD;
	dmat *opd=mktilt(loc, a, 0);
	dmat *g1=pygrad(py, opd);
	assert_float_equal(vmean(P(g1), nsa)/a, 1, 0.05);
	assert_float_equal(vmean(P(g1)+nsa, nsa)/a, 0, 0.05);
	dfree(opd);
	opd=mktilt(loc, 2*a, 0);
	dmat *g2=pygrad(py, opd);
	dmat *g2h=ddup(g2); dscale(g2h, 0.5);
	assert_dmat_near(g2h, g1, 0.05*dmaxabs(g1));
	dfree(opd); dfree(g2h);
	opd=mktilt(loc, -a, 0);
	dmat *gm=pygrad(py, opd);
	dscale(gm, -1);
	assert_dmat_near(gm, g1, 0.05*dmaxabs(g1));
	dfree(opd);
	/*y tilt*/
	opd=mktilt(loc, 0, a);
	dmat *gy=pygrad(py, opd);
	assert_float_equal(vmean(P(gy)+nsa, nsa)/a, 1, 0.05);
	assert_float_equal(vmean(P(gy), nsa)/a, 0, 0.05);
	dfree(opd);
	/*interaction matrix by poking zernike modes (radial order 1 to 3)*/
	dmat *mod=zernike(loc, D, 1, 3, 0);
	const long nmod=NY(mod);
	assert_int_equal(nmod, 9);
	dmat *G=dnew(nsa*2, nmod);
	const real poke=20e-9;
	for(long im=0; im<nmod; im++){
		dmat *opdi=NULL;
		dmat *col=drefcols(mod, im, 1);
		dadd(&opdi, 0, col, poke);
		dmat *gi=pygrad(py, opdi);
		for(long i=0; i<nsa*2; i++) P(G, i, im)=P(gi, i)/poke;
		dfree(gi); dfree(opdi); dfree(col);
	}
	assert_int_equal(NX(G), nsa*2);
	assert_int_equal(NY(G), nmod);
	/*Z2 = 2x/R (R=D/2): a unit coefficient is a tilt of 2/R*/
	const real ts=2/(D*0.5);
	dmat *gtt=drefcols(py->GTT, 0, 1);
	dmat *gz=drefcols(G, 0, 1);
	dmat *gtts=ddup(gtt); dscale(gtts, ts);
	assert_dmat_near(gz, gtts, 0.1*dmaxabs(gtts));
	dfree(gtt); dfree(gz); dfree(gtts);
	/*all modes are sensed: G'G is well conditioned*/
	dmat *GG=NULL;
	dmm(&GG, 0, G, G, "tn", 1);
	dmat *U=NULL, *S=NULL, *Vt=NULL;
	dsvd(&U, &S, &Vt, GG);
	assert_true(P(S, nmod-1)>1e-4*P(S, 0));
	dfree(U); dfree(S); dfree(Vt); dfree(GG);
	dfree(G); dfree(mod);
	dfree(g1); dfree(g2); dfree(gm); dfree(gy);
	pywfs_free(py);
	pycfg_free(cfg);
}
/**
   Variants: raw intensities as gradients, 3-sided pyramid, roof (2 sides),
   Zernike sensor (1 side), individual and no
   signal level matching, subaperture threshold, pupil shift, field stop,
   pyramid edge defects, a single modulation position, multiple rings, and
   pupil elongation. Flat wavefront always gives zero gradients and a tilt is
   sensed with the right sign.
*/
static void pywfs_variants(void **state){
	(void)state;
	const real a=0.002*AS2RAD;
	for(int iv=0; iv<11; iv++){
		pywfs_cfg_t *cfg=mkcfg(iv==1?3:(iv==9?2:(iv==10?1:4)));
		if(iv==1) cfg->modulpos=0;/*default: 30 positions, multiple of 3*/
		switch(iv){
		case 0: cfg->raw=1; break;
		case 2: cfg->sigmatch=1; break;
		case 3: cfg->sigmatch=0; break;
		case 4: cfg->saat=0.3; break;
		case 5:
			cfg->psx=dnew(4, 1); cfg->psy=dnew(4, 1);
			P(cfg->psx, 1)=0.2; P(cfg->psy, 2)=-0.1; P(cfg->psx, 3)=1;/*integer part is dropped*/
			break;
		case 6: cfg->fieldstop=3*AS2RAD; cfg->flate=0.001*AS2RAD; cfg->flatv=0.002*AS2RAD; break;
		case 7: cfg->modulring=2; cfg->modulpos_i=0; break;
		case 8: cfg->pupelong=0.1; break;
		}
		pywfs_t *py=pywfs_new(cfg, NULL, NULL, NULL);
		const loc_t *loc=py->locfft->loc;
		const long nsa=py->saloc->nloc;
		const int ng=pywfs_ng(cfg);
		assert_true(py->gain>0);
		dmat *opd=dnew(loc->nloc, 1);
		dmat *g0=pygrad(py, opd);
		assert_int_equal(PN(g0), nsa*ng);
		assert_float_equal(dmaxabs(g0), 0, 1e-12);
		dfree(opd);
		opd=mktilt(loc, a, 0);
		dmat *g=pygrad(py, opd);
		if(iv==9){
			/*roof: two pupils; light moves from one to the other*/
			assert_int_equal(ng, 2);
			real s0=vmean(P(g), nsa), s1=vmean(P(g)+nsa, nsa);
			assert_true(fabs(s0)>0);
			assert_float_equal(s0, -s1, 0.05*fabs(s0));
		} else if(iv==10){
			/*zernike sensor: one pupil image*/
			assert_int_equal(ng, 1);
			assert_true(dmaxabs(g)>0);
		} else if(iv==0){
			/*raw: normalized intensities of the 4 pupils; x tilt moves light
			  from pupils 1,3 to 0,2 (or the reverse), not between top/bottom*/
			assert_int_equal(ng, 4);
			real s[4];
			for(int i=0; i<4; i++) s[i]=vmean(P(g)+i*nsa, nsa);
			assert_true(fabs(s[0]+s[2]-s[1]-s[3])>10*fabs(s[0]+s[1]-s[2]-s[3]));
		} else{
			assert_int_equal(ng, 2);
			real gxm=vmean(P(g), nsa)/a;
			real gym=vmean(P(g)+nsa, nsa)/a;
			assert_true(gxm>0.5&&gxm<1.5);
			assert_true(fabs(gym)<0.2);
		}
		if(iv==4){
			assert_true(nsa<ORDER*ORDER);/*corners dropped*/
			assert_float_equal(dsum(py->saa), nsa, 1e-8);
		}
		if(iv==5){
			assert_non_null(py->pupilshift);
			assert_float_equal(P(py->pupilshift, 3, 0), 0, 0);
			assert_float_equal(P(py->pupilshift, 1, 0), 0.2, 1e-15);
		}
		if(iv==8){
			assert_non_null(py->msaloc);
			assert_int_equal(NX(py->msaloc), 4);
		}
		dfree(g); dfree(g0); dfree(opd);
		pywfs_free(py);
		pycfg_free(cfg);
	}
}
/**
   pywfs_simu() equals pywfs_new()+pywfs_ints()+pywfs_grad() with the same
   configuration, including an amplitude map and a user supplied loc.
   pywfs_gain_calibrate() keeps the gain near its calibrated value for a
   wavefront within the linear range.
*/
static void pywfs_simu_gain(void **state){
	(void)state;
	loc_t *loc=mkannloc(D, 0, D/(ORDER*4), 0);
	dmat *amp=dnew(loc->nloc, 1); dset(amp, 1);
	dmat *opd=mktilt(loc, 0.001*AS2RAD, -0.0005*AS2RAD);
	/*reference*/
	pywfs_cfg_t *cfg=mkcfg(4);
	cfg->D=0; cfg->dx=0;/*taken from loc*/
	pywfs_t *py=pywfs_new(cfg, loc, amp, NULL);
	assert_float_equal(cfg->dx, loc->dx, 0);
	assert_float_equal(cfg->D, loc_diam(loc), 0);
	assert_ptr_equal(py->loc->locx, loc->locx);
	dmat *gref=pygrad(py, opd);
	/*wrapper with a configuration*/
	pywfs_cfg_t *cfg2=mkcfg(4);
	cfg2->D=0; cfg2->dx=0;
	dmat *ints=NULL, *grad=NULL;
	//siglev is optional when pycfg is set: pycfg->siglev is used.
	pywfs_simu(&ints, &grad, cfg2, 0, NULL, 0, loc, amp, opd);
	assert_dmat_near(grad, gref, 1e-12*dmaxabs(gref));
	assert_int_equal(NY(ints), 4);
	real isum=dsum(ints);
	dfree(ints); dfree(grad);
	//a different siglev argument does not override pycfg->siglev.
	pywfs_simu(&ints, &grad, cfg2, 0, NULL, SIGLEV*10, loc, amp, opd);
	assert_float_equal(dsum(ints), isum, 1e-10*isum);
	dfree(ints); dfree(grad);
	/*wrapper without configuration and without ints output*/
	dmat *wvl=dnew(1, 1); P(wvl, 0)=WVL;
	pywfs_simu(NULL, &grad, NULL, ORDER, wvl, SIGLEV, loc, amp, opd);
	assert_int_equal(PN(grad), ORDER*ORDER*2);
	/*default modulation (32 positions): average gradient equals the tilt*/
	assert_float_equal(vmean(P(grad), ORDER*ORDER)/(0.001*AS2RAD), 1, 0.05);
	assert_float_equal(vmean(P(grad)+ORDER*ORDER, ORDER*ORDER)/(-0.0005*AS2RAD), 1, 0.05);
	dfree(grad); dfree(wvl);
	pycfg_free(cfg2);
	/*gain calibration (r0=0: no fitting error added)*/
	real gain0=py->gain;
	pywfs_gain_calibrate(py, gref, 0);
	assert_true(isfinite(py->gain));
	assert_true(py->gain>0.5*gain0&&py->gain<2*gain0);
	/*no debug: pywfs_test returns immediately*/
	pywfs_test(py);
	pywfs_test(NULL);
	dfree(gref);
	pywfs_free(py);
	pycfg_free(cfg);
	dfree(opd); dfree(amp); locfree(loc);
}
UNITTEST_MAIN(
	cmocka_unit_test(pywfs_setup),
	cmocka_unit_test(pywfs_gradoff),
	cmocka_unit_test(pywfs_linear),
	cmocka_unit_test(pywfs_variants),
	cmocka_unit_test(pywfs_simu_gain)
)
