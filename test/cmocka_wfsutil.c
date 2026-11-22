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
  Test WFS utilities moved from maos to lib: NEA conversions (lib/nea.c),
  dithering (lib/dither.c), the pupil amplitude (mkamp), and the parabolic peak
  of a correlation (lib/libmisc.c).
*/
#include "unittest.h"

/** NEA: Cholesky factor, its product, inverse, and the sparse matrix. */
static void nea_test(void **state){
	(void)state;
	const long nsa=3;
	dmat *nea=dnew(nsa, 3);
	for(long isa=0; isa<nsa; isa++){
		P(nea, isa, 0)=4+isa;//xx
		P(nea, isa, 1)=2+isa;//yy
		P(nea, isa, 2)=isa*0.5;//xy
	}
	dmat *chol=NULL;
	nea_chol(&chol, nea, 2);
	dmat *nea2=NULL;
	nea_mm(&nea2, chol, 2);//LL'
	assert_dmat_near(nea2, nea, 1e-14);
	dmat *inv=NULL;
	nea_inv(&inv, nea, 2, 1);
	for(long isa=0; isa<nsa; isa++){//inv*nea=identity
		real a=P(nea, isa, 0), b=P(nea, isa, 1), c=P(nea, isa, 2);
		real ia=P(inv, isa, 0), ib=P(inv, isa, 1), ic=P(inv, isa, 2);
		assert_float_equal(ia*a+ic*c, 1, 1e-14);
		assert_float_equal(ia*c+ic*b, 0, 1e-14);
		assert_float_equal(ic*c+ib*b, 1, 1e-14);
	}
	nea_inv(&inv, nea, 2, 2);//scaled
	assert_float_equal(P(inv, 0, 0), 2*P(nea, 0, 1)/(P(nea, 0, 0)*P(nea, 0, 1)), 1e-14);

	dsp *sp=nea2sp(nea, 1, 1, 2);
	dmat *full=NULL;
	dspfull(&full, sp, 'n', 1);
	assert_int_equal(NX(full), nsa*2);
	for(long isa=0; isa<nsa; isa++){
		assert_float_equal(P(full, isa, isa), P(nea, isa, 0), 0);
		assert_float_equal(P(full, isa+nsa, isa+nsa), P(nea, isa, 1), 0);
		assert_float_equal(P(full, isa, isa+nsa), P(nea, isa, 2), 0);
		assert_float_equal(P(full, isa+nsa, isa), P(nea, isa, 2), 0);
	}
	dfree(full);
	dspfree(sp);
	sp=nea2sp(nea, 0, 0, 2);//diagonal only
	assert_int_equal(sp->nzmax, nsa*2);
	dspfree(sp);

	dmat *col=dnew(nsa*2, 1);//a single column is reshaped
	nea_check(col, nsa, 2);
	assert_int_equal(NX(col), nsa);
	assert_int_equal(NY(col), 2);
	dfree(col);
	dfree(inv);
	dfree(nea2);
	dfree(chol);
	dfree(nea);
}

/** Dithering: position of the dither and the estimated amplitude. */
static void dither_test(void **state){
	(void)state;
	const int npoint=4;
	real cs, ss;
	for(int isim=0; isim<8; isim++){//no delay
		dither_position(&cs, &ss, 0, 1, npoint, isim, 0);
		real angle=M_PI*0.5+isim*2*M_PI/npoint;
		assert_float_equal(cs, cos(angle), 1e-14);
		assert_float_equal(ss, sin(angle), 1e-14);
	}
	const real amp=0.3;
	const long nframe=npoint*10;
	dmat *signal=dnew(1, nframe);//one mode
	for(long i=0; i<nframe; i++){
		P(signal, i)=amp*cos(M_PI*0.5+i*2*M_PI/npoint)+1e-3*i;//with a drift
	}
	assert_float_equal(dither_amp_calc(NULL, signal, 1, npoint, 1, 0), amp, 1e-12);
	dmat *tt=dnew(2, nframe);//tip and tilt
	for(long i=0; i<nframe; i++){
		real angle=M_PI*0.5+i*2*M_PI/npoint;
		P(tt, 0, i)=amp*cos(angle);
		P(tt, 1, i)=amp*sin(angle);
	}
	assert_float_equal(dither_amp_calc(NULL, tt, 1, npoint, 0, 1), amp, 1e-12);
	dmat *res=NULL;
	dither_amp_calc(&res, tt, 1, npoint, 0, 0);
	assert_int_equal(NX(res), 2);
	assert_float_equal(P(res, 0), amp, 1e-12);
	assert_float_equal(P(res, 1), amp, 1e-12);
	dfree(res);
	dfree(tt);
	dfree(signal);
}

/** Pupil amplitude of an annular aperture without a map. */
static void mkamp_test(void **state){
	(void)state;
	loc_t *loc=mksqloc_auto(64, 64, 1./16, 1./16);
	dmat *amp=mkamp(loc, NULL, 0, 0, 3, 1);
	for(long i=0; i<loc->nloc; i++){
		real r=sqrt(pow(loc->locx[i], 2)+pow(loc->locy[i], 2));
		if(r<0.4||r>1.6){
			assert_float_equal(P(amp, i), 0, 0);
		} else if(r>0.6&&r<1.4){
			assert_float_equal(P(amp, i), 1, 1e-12);
		}
	}
	real area=dsum(amp)/(16*16);
	assert_float_equal(area, M_PI*(1.5*1.5-0.5*0.5), 0.05);
	dfree(amp);
	locfree(loc);
}

/** The peak of a parabola is found exactly. */
static void parabolic_peak_test(void **state){
	(void)state;
	dmat *corr=dnew(7, 1);
	for(long i=0; i<7; i++){
		P(corr, i)=10-pow(i-3.3, 2);
	}
	assert_float_equal(parabolic_peak_1d(corr), 0.3, 1e-12);//relative to the center
	dmat *corr2=dnew(2, 1);
	assert_float_equal(parabolic_peak_1d(corr2), 0, 0);//too short
	dfree(corr2);
	dmat *c2=dnew(7, 7);//symmetric: no shift
	for(long iy=0; iy<7; iy++){
		for(long ix=0; ix<7; ix++){
			P(c2, ix, iy)=10-pow(ix-3, 2)-pow(iy-3, 2);
		}
	}
	real grad[2]={1,1};
	parabolic_peak_sum(grad, c2, 5);
	assert_float_equal(grad[0], 0, 1e-12);
	assert_float_equal(grad[1], 0, 1e-12);
	dfree(c2);
	dfree(corr);
}

UNITTEST_MAIN(
	cmocka_unit_test(nea_test),
	cmocka_unit_test(dither_test),
	cmocka_unit_test(mkamp_test),
	cmocka_unit_test(parabolic_peak_test),
)
