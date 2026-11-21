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

typedef struct{
	loc_t *loc;
	dmat *amp;
}pupil_t;
/**
   Create a circular pupil of diameter D sampled at dx.
*/
static pupil_t pupil_new(real D, real dx){
	pupil_t pp;
	pp.loc=mkannloc(D, 0, dx, 0.5);
	pp.amp=dnew(pp.loc->nloc, 1);
	dset(pp.amp, 1);
	return pp;
}
static void pupil_free(pupil_t *pp){
	locfree(pp->loc);
	dfree(pp->amp);
}
/** Index of the maximum of abs(A). */
static long cmaxloc(const cmat *A){
	long imax=0;
	real vmax=-1;
	for(long i=0; i<PN(A); i++){
		real v=cabs(P(A, i));
		if(v>vmax){
			vmax=v;
			imax=i;
		}
	}
	return imax;
}
/** Diffraction limited PSF has peak (Strehl) 1 at the center, or total energy 1 with sum2one. */
static void locfft_strehl(void **state){
	(void)state;
	pupil_t pp=pupil_new(1, 1./32);
	dmat *wvl=dnew(1, 1);
	P(wvl, 0)=1e-6;
	locfft_t *locfft=locfft_init(pp.loc, pp.amp, wvl, NULL, 2, 0);
	long nembed=P(locfft->nembed, 0);
	assert_true(nembed>=2*32);
	dmat *opd=dnew(pp.loc->nloc, 1);
	ccell *psf=NULL;
	locfft_psf(&psf, locfft, opd, NULL, 0);
	cmat *psf0=P(psf, 0);
	assert_int_equal(NX(psf0), nembed);
	assert_int_equal(cmaxloc(psf0), nembed/2*(1+nembed));
	assert_float_equal(cabs(P(psf0, nembed/2, nembed/2)), 1, 1e-10);
	//piston does not change the Strehl
	dset(opd, 0.3e-6);
	locfft_psf(&psf, locfft, opd, NULL, 0);
	assert_float_equal(cabs(P(P(psf, 0), nembed/2, nembed/2)), 1, 1e-10);
	//sum2one normalizes the total energy.
	locfft_psf(&psf, locfft, opd, NULL, 1);
	assert_float_equal(csumsq(P(psf, 0)), 1, 1e-10);
	//psfsize=1 only computes the Strehl.
	lmat *psfsize=lnew(1, 1);
	P(psfsize, 0)=1;
	ccell *psf1=NULL;
	locfft_psf(&psf1, locfft, opd, psfsize, 0);
	assert_int_equal(PN(P(psf1, 0)), 1);
	assert_float_equal(cabs(P(P(psf1, 0), 0)), 1, 1e-10);
	ccellfree(psf1);
	lfree(psfsize);
	ccellfree(psf);
	dfree(opd);
	locfft_free(locfft);
	dfree(wvl);
	pupil_free(&pp);
}
/** Small random aberration reduces the Strehl as exp(-sigma^2) (Marechal approximation). */
static void locfft_marechal(void **state){
	(void)state;
	pupil_t pp=pupil_new(1, 1./32);
	dmat *wvl=dnew(1, 1);
	P(wvl, 0)=1e-6;
	locfft_t *locfft=locfft_init(pp.loc, pp.amp, wvl, NULL, 2, 0);
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat *opd=dnew(pp.loc->nloc, 1);
	drandn(opd, 50e-9, &rstat);
	real rms=sqrt(dsumsq(opd)/PN(opd)-pow(dsum(opd)/PN(opd), 2));
	real sigma=TWOPI*rms/P(wvl, 0);
	ccell *psf=NULL;
	locfft_psf(&psf, locfft, opd, NULL, 0);
	long nembed=P(locfft->nembed, 0);
	real strehl=cabs2(P(P(psf, 0), nembed/2, nembed/2));
	assert_float_equal(strehl, exp(-sigma*sigma), 0.01);
	assert_true(strehl<1);
	ccellfree(psf);
	dfree(opd);
	locfft_free(locfft);
	dfree(wvl);
	pupil_free(&pp);
}
/** A tilt of k*lambda/(nembed*dx) moves the PSF peak by k pixels. */
static void locfft_tilt(void **state){
	(void)state;
	pupil_t pp=pupil_new(1, 1./32);
	dmat *wvl=dnew(1, 1);
	P(wvl, 0)=1e-6;
	locfft_t *locfft=locfft_init(pp.loc, pp.amp, wvl, NULL, 2, 0);
	long nembed=P(locfft->nembed, 0);
	real dtheta=P(wvl, 0)/(nembed*pp.loc->dx);//PSF sampling in radian
	dmat *opd=dnew(pp.loc->nloc, 1);
	const int kx=3, ky=-2;
	for(long i=0; i<pp.loc->nloc; i++){
		P(opd, i)=pp.loc->locx[i]*kx*dtheta+pp.loc->locy[i]*ky*dtheta;
	}
	ccell *psf=NULL;
	locfft_psf(&psf, locfft, opd, NULL, 0);
	long imax=cmaxloc(P(psf, 0));
	long ix=imax%nembed, iy=imax/nembed;
	//the sign follows the FFT convention: check the magnitude and that x and y are opposite.
	assert_int_equal(labs(ix-nembed/2), abs(kx));
	assert_int_equal(labs(iy-nembed/2), abs(ky));
	assert_true((ix-nembed/2)*(iy-nembed/2)<0);
	ccellfree(psf);
	dfree(opd);
	locfft_free(locfft);
	dfree(wvl);
	pupil_free(&pp);
}
/** Cropped PSF (including the partial FFT path for small crops) equals the center of the full PSF, for each wavelength. */
static void locfft_crop(void **state){
	(void)state;
	pupil_t pp=pupil_new(2, 1./64);
	dmat *wvl=dnew(2, 1);
	P(wvl, 0)=1e-6;
	P(wvl, 1)=2e-6;
	lmat *fftsize=lnew(2, 1);
	P(fftsize, 0)=P(fftsize, 1)=512;
	locfft_t *locfft=locfft_init(pp.loc, pp.amp, wvl, fftsize, 2, 0);
	assert_int_equal(P(locfft->nembed, 0), 512);
	assert_ptr_equal(P(P(locfft->embed, 1)), P(P(locfft->embed, 0)));//same size shares the embedding
	rand_t rstat;
	seed_rand(&rstat, 2);
	dmat *opd=dnew(pp.loc->nloc, 1);
	drandn(opd, 30e-9, &rstat);
	ccell *full=NULL;
	locfft_psf(&full, locfft, opd, NULL, 0);
	const long ncrop[2]={32, 400};//32: partial FFT is used. 400: full FFT then crop.
	for(int ic=0; ic<2; ic++){
		lmat *psfsize=lnew(2, 1);
		P(psfsize, 0)=P(psfsize, 1)=ncrop[ic];
		ccell *crop=NULL;
		locfft_psf(&crop, locfft, opd, psfsize, 0);
		for(int iwvl=0; iwvl<2; iwvl++){
			cmat *c=P(crop, iwvl);
			cmat *f=P(full, iwvl);
			assert_int_equal(NX(c), ncrop[ic]);
			long off=(NX(f)-NX(c))/2;
			real maxdiff=0;
			for(long iy=0; iy<NY(c); iy++){
				for(long ix=0; ix<NX(c); ix++){
					maxdiff=MAX(maxdiff, cabs(P(c, ix, iy)-P(f, ix+off, iy+off)));
				}
			}
			assert_float_equal(maxdiff, 0, 1e-10);
		}
		ccellfree(crop);
		lfree(psfsize);
	}
	//Longer wavelength has higher Strehl
	long n2=512/2;
	assert_true(cabs(P(P(full, 1), n2, n2))>cabs(P(P(full, 0), n2, n2)));
	ccellfree(full);
	dfree(opd);
	locfft_free(locfft);
	lfree(fftsize);
	dfree(wvl);
	pupil_free(&pp);
}
/** A wide field stop keeps a smooth OPD, while a narrow one removes high spatial frequencies. */
static void locfft_fieldstop_test(void **state){
	(void)state;
	pupil_t pp=pupil_new(1, 1./32);
	dmat *wvl=dnew(1, 1);
	P(wvl, 0)=1e-6;
	dmat *wvlwts=dnew(1, 1);
	P(wvlwts, 0)=1;
	dmat *opd=dnew(pp.loc->nloc, 1);
	for(long i=0; i<pp.loc->nloc; i++){//low order (focus like) plus a high frequency ripple
		real x=pp.loc->locx[i], y=pp.loc->locy[i];
		P(opd, i)=20e-9*(x*x+y*y)*4+5e-9*cos(TWOPI*x/(4*pp.loc->dx));
	}
	locfft_t *wide=locfft_init(pp.loc, pp.amp, wvl, NULL, 2, 1);//1 radian: no effect
	dmat *opd1=ddup(opd);
	locfft_fieldstop(wide, opd1, wvlwts);
	assert_dmat_near(opd, opd1, 1e-12);
	locfft_t *narrow=locfft_init(pp.loc, pp.amp, wvl, NULL, 2, 4*P(wvl, 0));//4 lambda/D
	dmat *opd2=ddup(opd);
	locfft_fieldstop(narrow, opd2, wvlwts);
	//the ripple is removed: project the OPD onto the ripple pattern.
	real rr=0, c0=0, c2=0;
	for(long i=0; i<pp.loc->nloc; i++){
		real rip=cos(TWOPI*pp.loc->locx[i]/(4*pp.loc->dx));
		rr+=rip*rip;
		c0+=rip*P(opd, i);
		c2+=rip*P(opd2, i);
	}
	assert_float_equal(c0/rr, 5e-9, 1e-9);
	assert_float_equal(c2/rr, 0, 0.5e-9);
	dfree(opd1);
	dfree(opd2);
	locfft_free(wide);
	locfft_free(narrow);
	dfree(opd);
	dfree(wvlwts);
	dfree(wvl);
	pupil_free(&pp);
}
/** Create a circular aperture complex field of nx*nx. */
static cmat *wvf_circle(long nx){
	cmat *wvf=cnew(nx, nx);
	ccircle(wvf, nx/2, nx/2, 1, 1, nx/4, 1);
	return wvf;
}
/** Automatic selection uses angular spectrum (output sampling equals input) in the near field and the FFT method otherwise. */
static void fresnel_auto(void **state){
	(void)state;
	const long nx=64;
	const real dx=1./64, wvl=0.5e-6;
	cmat *in=wvf_circle(nx);
	cmat *out=NULL, *out2=NULL;
	real dxout=0, dxout2=0;
	real z=1e4;//dp=z*wvl/(dx*nx)<dx: angular spectrum
	fresnel_prop(&out, &dxout, in, dx, wvl, z, 1, -1);
	assert_float_equal(dxout, dx, 0);
	fresnel_prop(&out2, &dxout2, in, dx, wvl, z, 1, 2);
	assert_float_equal(creal(cdiff(out, out2)), 0, 0);
	z=1e5;//dp>dx: FFT method
	fresnel_prop(&out, &dxout, in, dx, wvl, z, 1, -1);
	assert_float_equal(dxout, z*wvl/(dx*nx), 1e-15);
	fresnel_prop(&out2, &dxout2, in, dx, wvl, z, 1, 1);
	assert_float_equal(creal(cdiff(out, out2)), 0, 0);
	cfree(in);
	cfree(out);
	cfree(out2);
}
/** FFT propagation sets the output sampling to z*wvl/(dx*nx), scales the energy by nx*ny, and scale=0 normalizes the output. */
static void fresnel_fft(void **state){
	(void)state;
	const long nx=32;
	const real dx=1./32, wvl=0.5e-6, z=1e5;
	cmat *in=wvf_circle(nx);
	cmat *out=NULL;
	real dxout=0;
	fresnel_prop(&out, &dxout, in, dx, wvl, z, 1, 1);
	assert_float_equal(dxout, z*wvl/(dx*nx), 1e-15);
	assert_float_equal(csumsq(out), csumsq(in)*PN(in), 1e-9*csumsq(in)*PN(in));
	//The far field of a symmetric aperture is symmetric and peaks at the center.
	real vc=cabs(P(out, nx/2, nx/2));
	for(long i=0; i<PN(out); i++){
		assert_true(cabs(P(out, i))<=vc*(1+1e-12));
	}
	assert_float_equal(cabs(P(out, nx/2+3, nx/2)), cabs(P(out, nx/2-3, nx/2)), 1e-9*vc);
	assert_float_equal(cabs(P(out, nx/2, nx/2+5)), cabs(P(out, nx/2+5, nx/2)), 1e-9*vc);
	//scale=0 determines the scale from the input to normalize the energy.
	fresnel_prop(&out, &dxout, in, dx, wvl, z, 0, 1);
	assert_float_equal(csumsq(out), 1, 1e-9);
	cfree(in);
	cfree(out);
}
/** Angular spectrum propagation scales the energy by nx*ny like the FFT method, and propagating back recovers the input. */
static void fresnel_angular(void **state){
	(void)state;
	const long nx=64;
	const real dx=1./64, wvl=0.5e-6, z=1e4;
	cmat *in=wvf_circle(nx);
	cmat *out=NULL, *back=NULL;
	real dxout=0, dxback=0;
	fresnel_prop(&out, &dxout, in, dx, wvl, z, 1, 2);
	assert_float_equal(dxout, dx, 0);
	assert_float_equal(csumsq(out), csumsq(in)*PN(in), 1e-9*csumsq(in)*PN(in));
	fresnel_prop(&back, &dxback, out, dxout, wvl, -z, 1./PN(in), 2);
	assert_float_equal(creal(cdiff(back, in)), 0, 1e-9*sqrt(csumsq(in)));
	cfree(in);
	cfree(out);
	cfree(back);
}
/** The FFT and angular spectrum methods agree when their output sampling is the same. */
static void fresnel_agree(void **state){
	(void)state;
	const long nx=64;
	const real dx=0.5, wvl=0.5e-6;
	const real z=dx*dx*nx/wvl;//output sampling of the FFT method equals dx
	cmat *in=wvf_circle(nx);
	cmat *out1=NULL, *out2=NULL;
	real dx1=0, dx2=0;
	fresnel_prop(&out1, &dx1, in, dx, wvl, z, 1, 1);
	fresnel_prop(&out2, &dx2, in, dx, wvl, z, 1, 2);
	assert_float_equal(dx1, dx2, 1e-12);
	//compare the intensities as the methods may differ by a global phase.
	real maxdiff=0, maxval=0;
	for(long i=0; i<PN(in); i++){
		maxdiff=MAX(maxdiff, fabs(cabs2(P(out1, i))-cabs2(P(out2, i))));
		maxval=MAX(maxval, cabs2(P(out1, i)));
	}
	assert_float_equal(maxdiff/maxval, 0, 1e-9);
	cfree(in);
	cfree(out1);
	cfree(out2);
}
/** Brute force propagation agrees with the FFT method in the far field. */
static void fresnel_brute(void **state){
	(void)state;
	const long nx=16;
	const real dx=1./16, wvl=0.5e-6, z=1e5;
	cmat *in=wvf_circle(nx);
	cmat *out1=NULL, *out0=cnew(nx, nx);
	real dx1=0;
	fresnel_prop(&out1, &dx1, in, dx, wvl, z, 1, 1);
	real dx0=dx1;
	fresnel_prop(&out0, &dx0, in, dx, wvl, z, 1, 0);
	assert_float_equal(dx0, dx1, 0);
	real maxdiff=0, maxval=0;
	for(long i=0; i<PN(in); i++){
		maxdiff=MAX(maxdiff, fabs(cabs(P(out1, i))-cabs(P(out0, i))));
		maxval=MAX(maxval, cabs(P(out1, i)));
	}
	assert_float_equal(maxdiff/maxval, 0, 1e-3);
	cfree(in);
	cfree(out0);
	cfree(out1);
}

UNITTEST_MAIN(
	cmocka_unit_test(locfft_strehl),
	cmocka_unit_test(locfft_marechal),
	cmocka_unit_test(locfft_tilt),
	cmocka_unit_test(locfft_crop),
	cmocka_unit_test(locfft_fieldstop_test),
	cmocka_unit_test(fresnel_auto),
	cmocka_unit_test(fresnel_fft),
	cmocka_unit_test(fresnel_angular),
	cmocka_unit_test(fresnel_agree),
	cmocka_unit_test(fresnel_brute),
)
