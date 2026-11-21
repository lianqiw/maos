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
  Tests of the Shack-Hartmann WFS building blocks in lib/: average gradient
  operator (mkg.c), Z-tilt operator (mkz.c), detector and elongation transfer
  functions (mkdtf.c), subaperture OTF (genotf.c), short exposure pixel
  intensities (gensei.c), matched filter (mtch.c), Gerchberg-Saxton phase
  retrieval (phase.c), and the cumulative reconstructor (cure.c).
*/
#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-9
#endif
/*Common SHWFS geometry: 4x4 subapertures of 0.5 m sampled at 1/16 m.*/
static const real DSA=0.5;
static const real DXP=1./16;
static const long NSA1=4;/*subapertures along each dimension*/

/** OPD a*x+b*y+c*x^2 on loc. */
static dmat *mkopd(const loc_t *loc, real a, real b, real c){
	dmat *o=dnew(loc->nloc, 1);
	for(long i=0; i<loc->nloc; i++){
		real x=loc->locx[i], y=loc->locy[i];
		P(o, i)=a*x+b*y+c*x*x;
	}
	return o;
}
/** Pupil grid covering the subapertures with a margin of 2 points.*/
static loc_t *mkploc(void){
	long n=(long)round(NSA1*DSA/DXP)+4;
	return mksqloc(n, n, DXP, DXP, -NSA1*DSA*0.5-2*DXP, -NSA1*DSA*0.5-2*DXP);
}
static loc_t *mksaloc(void){
	return mksqloc(NSA1, NSA1, DSA, DSA, -NSA1*DSA*0.5, -NSA1*DSA*0.5);
}
static dmat *spmulvec(const dsp *A, const dmat *x){
	dmat *y=NULL;
	dspmm(&y, A, x, "nn", 1);
	return y;
}
/** Check that g contains slope a for x and b for y for all nsa subapertures (skipping those with skip[isa]). */
static void check_slope(const dmat *g, long nsa, real a, real b, real tol, const int *skip){
	assert_int_equal(PN(g), nsa*2);
	for(long isa=0; isa<nsa; isa++){
		if(skip&&skip[isa]) continue;
		assert_float_equal(P(g, isa), a, tol);
		assert_float_equal(P(g, isa+nsa), b, tol);
	}
}
/**
   mkg(): the average gradient of a linear phase equals its slope and of a
   piston is zero, for the same and different (shifted, coarser) OPD grids,
   with and without amplitude, partial subapertures, and subaperture area
   threshold. The average gradient of x^2 equals twice the subaperture center.
   mkgt() is the transpose.
*/
static void wfs_mkg(void **state){
	(void)state;
	loc_t *ploc=mkploc();
	loc_t *saloc=mksaloc();
	const long nsa=saloc->nloc;
	dsp *G=mkg(ploc, ploc, NULL, saloc, NULL, 0, 1, 0, 0, 0);
	assert_int_equal(NX(G), nsa*2);
	assert_int_equal(NY(G), ploc->nloc);
	dmat *opd=mkopd(ploc, 0.3, -0.2, 0);
	dmat *g=spmulvec(G, opd);
	check_slope(g, nsa, 0.3, -0.2, TOL, NULL);
	dfree(g); dfree(opd);
	/*piston*/
	opd=dnew(ploc->nloc, 1); dset(opd, 5);
	g=spmulvec(G, opd);
	assert_float_equal(dmaxabs(g), 0, TOL);
	dfree(g); dfree(opd);
	/*x^2: average gradient over the subaperture is 2*center*/
	opd=mkopd(ploc, 0, 0, 1);
	g=spmulvec(G, opd);
	for(long isa=0; isa<nsa; isa++){
		assert_float_equal(P(g, isa), 2*(saloc->locx[isa]+DSA*0.5), 1e-6);
		assert_float_equal(P(g, isa+nsa), 0, TOL);
	}
	dfree(g); dfree(opd);
	/*transpose*/
	dsp *Gt=mkgt(ploc, ploc, NULL, saloc, NULL, 0, 1, 0, 0, 0);
	dsp *Gtt=dsptrans(Gt);
	dsp *diff=dspadd2(G, 1, Gtt, -1);
	dmat *fd=NULL;
	dspfull(&fd, diff, 'n', 1);
	assert_float_equal(dmaxabs(fd), 0, TOL);
	dfree(fd); dspfree(diff); dspfree(Gtt); dspfree(Gt);
	/*uniform amplitude and do_partial give the same operator for full subapertures*/
	dmat *amp=dnew(ploc->nloc, 1); dset(amp, 1);
	dsp *G2=mkg(ploc, ploc, amp, saloc, NULL, 0, 1, 0, 0, 1);
	opd=mkopd(ploc, -0.7, 0.4, 0);
	g=spmulvec(G2, opd);
	check_slope(g, nsa, -0.7, 0.4, TOL, NULL);
	dfree(g); dspfree(G2);
	/*circular amplitude: partially illuminated subapertures still measure the slope*/
	dmat *ampc=dnew(ploc->nloc, 1);
	for(long i=0; i<ploc->nloc; i++){
		real r=sqrt(pow(ploc->locx[i], 2)+pow(ploc->locy[i], 2));
		P(ampc, i)=r<1.05?1:0;
	}
	G2=mkg(ploc, ploc, ampc, saloc, NULL, 0, 1, 0, 0, 0);
	g=spmulvec(G2, opd);
	check_slope(g, nsa, -0.7, 0.4, 1e-6, NULL);
	dfree(g); dspfree(G2);
	/*shifted subapertures (dispx/dispy) on the same grid*/
	G2=mkg(ploc, ploc, NULL, saloc, NULL, 0, 1, DXP*0.5, -DXP, 0);
	g=spmulvec(G2, opd);
	check_slope(g, nsa, -0.7, 0.4, 1e-6, NULL);
	dfree(g); dspfree(G2);
	dfree(opd);
	/*subaperture area threshold drops the weak subaperture*/
	dmat *saa=dnew(nsa, 1); dset(saa, 1);
	P(saa, 3)=0.1;
	G2=mkg(ploc, ploc, NULL, saloc, saa, 0.5, 1, 0, 0, 0);
	dmat *G2f=NULL;
	dspfull(&G2f, G2, 'n', 1);
	for(long i=0; i<ploc->nloc; i++){
		assert_float_equal(P(G2f, 3, i), 0, 0);
		assert_float_equal(P(G2f, 3+nsa, i), 0, 0);
	}
	dfree(G2f); dspfree(G2); dfree(saa);
	/*different OPD grid (coarser, shifted, larger): bilinear interpolation of a linear phase*/
	loc_t *xloc=mksqloc(24, 24, 0.125, 0.125, -1.5+0.03, -1.5-0.02);
	dsp *G3=mkg(xloc, ploc, amp, saloc, NULL, 0, 1, 0.01, -0.02, 0);
	assert_int_equal(NY(G3), xloc->nloc);
	opd=mkopd(xloc, 0.25, 0.6, 0);
	g=spmulvec(G3, opd);
	check_slope(g, nsa, 0.25, 0.6, 1e-6, NULL);
	dfree(g); dfree(opd);
	dspfree(G3);
	/*same, without amplitude*/
	G3=mkg(xloc, ploc, NULL, saloc, NULL, 0, 1, 0, 0, 0);
	opd=mkopd(xloc, 0.25, 0.6, 0);
	g=spmulvec(G3, opd);
	check_slope(g, nsa, 0.25, 0.6, 1e-6, NULL);
	dfree(g); dfree(opd); dspfree(G3);
	locfree(xloc);
	dfree(amp); dfree(ampc);
	dspfree(G);
	locfree(saloc); locfree(ploc);
}
/**
   mkg_ngs(): NGS mode to gradient interaction matrix has unit response of
   tip/tilt and the analytic response for focus, astigmatism, and plate scale
   modes.
*/
static void wfs_mkg_ngs(void **state){
	(void)state;
	const long nsa=3;
	dmat *sacent=dnew(nsa, 2);
	for(long i=0; i<nsa; i++){
		P(sacent, i, 0)=0.1*(i+1);
		P(sacent, i, 1)=-0.2*i;
	}
	/*LTAO: tip, tilt, astig(2), focus*/
	dmat *pg=mkg_ngs(sacent, 0, 0, 0, 90e3, 0, 2, 4, 0);
	assert_int_equal(NX(pg), nsa*2);
	assert_int_equal(NY(pg), 5);
	for(long i=0; i<nsa; i++){
		real xm=P(sacent, i, 0), ym=P(sacent, i, 1);
		assert_float_equal(P(pg, i, 0), 1, 0);
		assert_float_equal(P(pg, i+nsa, 0), 0, 0);
		assert_float_equal(P(pg, i+nsa, 1), 1, 0);
		assert_float_equal(P(pg, i, 2), 2*xm, 1e-15);
		assert_float_equal(P(pg, i+nsa, 2), -2*ym, 1e-15);
		assert_float_equal(P(pg, i, 3), ym, 1e-15);
		assert_float_equal(P(pg, i+nsa, 3), xm, 1e-15);
		assert_float_equal(P(pg, i, 4), 2*xm, 1e-15);
		assert_float_equal(P(pg, i+nsa, 4), 2*ym, 1e-15);
	}
	dfree(pg);
	/*MCAO plate scale modes: with hc=0 they reduce to focus/astig like terms*/
	const real thetax=1e-5, thetay=-2e-5, hc=10e3, hs=INFINITY;
	pg=mkg_ngs(sacent, thetax, thetay, hc, hs, 2, 0, 5, 0);
	assert_int_equal(NY(pg), 6);
	for(long i=0; i<nsa; i++){
		real xm=P(sacent, i, 0), ym=P(sacent, i, 1);
		/*hs=inf: scale=1, scale1=0*/
		assert_float_equal(P(pg, i, 2), -2*thetax*hc, 1e-15);
		assert_float_equal(P(pg, i+nsa, 2), -2*thetay*hc, 1e-15);
		assert_float_equal(P(pg, i, 3), -2*thetax*hc, 1e-15);
		assert_float_equal(P(pg, i+nsa, 3), 2*thetay*hc, 1e-15);
		assert_float_equal(P(pg, i, 4), -thetay*hc, 1e-15);
		assert_float_equal(P(pg, i+nsa, 4), -thetax*hc, 1e-15);
		assert_float_equal(P(pg, i, 5), 2*xm, 1e-15);
		assert_float_equal(P(pg, i+nsa, 5), 2*ym, 1e-15);
	}
	dfree(pg);
	/*ahstfocus removes focus from the first plate scale mode*/
	pg=mkg_ngs(sacent, 0, 0, hc, 90e3, 2, 0, 5, 1);
	for(long i=0; i<nsa; i++){
		assert_float_equal(P(pg, i, 2), 0, 1e-15);
	}
	dfree(pg);
	dfree(sacent);
}
/**
   mkz(): the Z-tilt of a linear phase equals its slope (agreeing with the
   G-tilt of mkg()), constant amplitude does not change the fit, and Z-tilt of
   x^2 is twice the mean x of the points in the subaperture.
*/
static void wfs_mkz(void **state){
	(void)state;
	loc_t *ploc=mkploc();
	loc_t *saloc=mksaloc();
	const long nsa=saloc->nloc;
	dsp *Z=mkz(ploc, NULL, saloc, 1, 1, 0, 0);
	dsp *G=mkg(ploc, ploc, NULL, saloc, NULL, 0, 1, 0, 0, 0);
	assert_int_equal(NX(Z), nsa*2);
	assert_int_equal(NY(Z), ploc->nloc);
	dmat *opd=mkopd(ploc, 0.15, 0.35, 0);
	dmat *gz=spmulvec(Z, opd);
	dmat *gg=spmulvec(G, opd);
	check_slope(gz, nsa, 0.15, 0.35, TOL, NULL);
	assert_dmat_near(gz, gg, TOL);
	dfree(gz); dfree(gg); dfree(opd);
	/*x^2*/
	opd=mkopd(ploc, 0, 0, 1);
	gz=spmulvec(Z, opd);
	dmat *Zf=NULL;
	dspfull(&Zf, Z, 'n', 1);
	for(long isa=0; isa<nsa; isa++){
		real xs=0; long np=0;
		for(long i=0; i<ploc->nloc; i++){
			if(P(Zf, isa, i)!=0){
				xs+=ploc->locx[i]; np++;
			}
		}
		assert_true(np>0);
		assert_float_equal(P(gz, isa), 2*xs/np, 1e-8);
		assert_float_equal(P(gz, isa+nsa), 0, 1e-8);
	}
	dfree(Zf); dfree(gz); dfree(opd);
	/*constant amplitude weighting gives the same operator*/
	dmat *amp=dnew(ploc->nloc, 1); dset(amp, 2);
	dsp *Z2=mkz(ploc, P(amp), saloc, 1, 1, 0, 0);
	dsp *diff=dspadd2(Z, 1, Z2, -1);
	dmat *fd=NULL;
	dspfull(&fd, diff, 'n', 1);
	assert_float_equal(dmaxabs(fd), 0, 1e-8);
	dfree(fd); dspfree(diff); dspfree(Z2); dfree(amp);
	/*saloc as subaperture center (saorc=0): shifted saloc gives the same operator*/
	loc_t *saloc2=mksaloc();
	for(long i=0; i<nsa; i++){
		saloc2->locx[i]+=DSA*0.5;
		saloc2->locy[i]+=DSA*0.5;
	}
	Z2=mkz(ploc, NULL, saloc2, 0, 1, 0, 0);
	diff=dspadd2(Z, 1, Z2, -1);
	dspfull(&fd, diff, 'n', 1);
	assert_float_equal(dmaxabs(fd), 0, 1e-10);
	dfree(fd); dspfree(diff); dspfree(Z2); locfree(saloc2);
	/*transpose*/
	dsp *Zt=mkzt(ploc, NULL, saloc, 1, 1, 0, 0);
	assert_int_equal(NX(Zt), ploc->nloc);
	assert_int_equal(NY(Zt), nsa*2);
	dspfree(Zt);
	dspfree(G); dspfree(Z);
	locfree(saloc); locfree(ploc);
}
/**
   dtf_otf(): pixel transfer function is the product of sincs (separable),
   with value pdtheta^2 at zero frequency, gaussian pixel blur factor, 90 deg
   rotation transposes, treat=2 is an fftshift and treat=1 an fftshift with a
   checkerboard sign.
*/
static void wfs_dtf_otf(void **state){
	(void)state;
	const long nx=16, ny=12;
	const real pdx=2, pdy=1.5;
	dmat *nom=NULL;
	dtf_otf(&nom, nx, ny, pdx, pdy, 0, 0, 0);
	assert_int_equal(NX(nom), nx);
	assert_int_equal(NY(nom), ny);
	const long cx=nx/2, cy=ny/2;
	assert_float_equal(P(nom, cx, cy), pdx*pdy, 1e-14);
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			real ref=sinc((ix-cx)*pdx/nx)*sinc((iy-cy)*pdy/ny)*pdx*pdy;
			assert_float_equal(P(nom, ix, iy), ref, 1e-14);
			/*separable*/
			assert_float_equal(P(nom, ix, iy)*P(nom, cx, cy), P(nom, ix, cy)*P(nom, cx, iy), 1e-13);
		}
	}
	/*pixel blur multiplies by gaussian*/
	dmat *nomb=NULL;
	const real blur=0.3;
	dtf_otf(&nomb, nx, ny, pdx, pdy, 0, -blur, 0);/*sign of blur ignored*/
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			real jx=ix-cx, jy=iy-cy;
			real fac=exp(-2*pow(M_PI*blur, 2)*(pow(jx*pdx/nx, 2)+pow(jy*pdy/ny, 2)));
			assert_float_equal(P(nomb, ix, iy), P(nom, ix, iy)*fac, 1e-14);
		}
	}
	dfree(nomb);
	/*no binning, only blur: value 1 at center*/
	dtf_otf(&nomb, nx, ny, 0, 0, 0, 0, 0);
	assert_float_equal(dmin(nomb), 1, 0);
	assert_float_equal(dmax(nomb), 1, 0);
	dfree(nomb);
	/*90 degree rotation of a square grid transposes the pixel*/
	dmat *n0=NULL, *n90=NULL;
	dtf_otf(&n0, nx, nx, pdx, pdy, 0, blur, 0);
	dtf_otf(&n90, nx, nx, pdx, pdy, M_PI/2, blur, 0);
	for(long iy=1; iy<nx; iy++){
		for(long ix=1; ix<nx; ix++){
			/*rotated coordinate: ir=jy, ia=-jx*/
			assert_float_equal(P(n90, ix, iy), P(n0, iy, 2*cx-ix), 1e-12);
		}
	}
	/*treat=2 is fftshift*/
	dmat *ns=NULL;
	dtf_otf(&ns, nx, nx, pdx, pdy, 0, blur, 2);
	dmat *ref=ddup(n0);
	dfftshift(ref);
	assert_dmat_near(ns, ref, 1e-14);
	assert_float_equal(P(ns, 0), pdx*pdy, 1e-14);
	/*treat=1 is fftshift with checkerboard sign*/
	dmat *nt=NULL;
	dtf_otf(&nt, nx, nx, pdx, pdy, 0, blur, 1);
	for(long iy=0; iy<nx; iy++){
		for(long ix=0; ix<nx; ix++){
			assert_float_equal(P(nt, ix, iy), P(ref, ix, iy)*(((ix+iy)%2)?-1:1), 1e-12);
		}
	}
	dfree(nt); dfree(ns); dfree(ref); dfree(n0); dfree(n90); dfree(nom);
}
/**
   genotf(): diffraction limited OTF of an n*n square subaperture is the
   normalized overlap area (n-|dx|)(n-|dy|)/n^2; white noise phase reduces
   off-center values by exp(-k^2 sigma^2); a tilt OPD bias adds a linear phase
   and the OTF is Hermitian; Kolmogorov OTF decreases with separation and
   p/t/t removal increases it; multiple subapertures with area threshold.
*/
static void wfs_genotf(void **state){
	(void)state;
	const long n=4, npsf=2*n;
	const real dx=1./8;
	loc_t *loc=mksqloc_auto(n, n, dx, dx);
	const long nloc=loc->nloc;
	const real wvl=1e-6;
	dmat *cov=dnew(nloc, nloc);
	ccell *otf=NULL;
	genotf(&otf, loc, NULL, NULL, NULL, 0, wvl, cov, 0, 0, npsf, npsf, 1, 0);
	assert_non_null(otf);
	cmat *o=P(otf, 0);
	assert_int_equal(NX(o), npsf);
	for(long iy=0; iy<npsf; iy++){
		for(long ix=0; ix<npsf; ix++){
			long jx=labs(ix-npsf/2), jy=labs(iy-npsf/2);
			real ref=(jx<n&&jy<n)?(real)(n-jx)*(n-jy)/(n*n):0;
			assert_float_equal(creal(P(o, ix, iy)), ref, 1e-12);
			assert_float_equal(cimag(P(o, ix, iy)), 0, 1e-12);
		}
	}
	/*white noise phase with variance sigma2 (m^2)*/
	const real sigma2=pow(wvl/(2*M_PI), 2)*0.3;/*0.3 rad^2*/
	daddI(cov, sigma2);
	genotf(&otf, loc, NULL, NULL, NULL, 0, wvl, cov, 0, 0, npsf, npsf, 1, 0);
	o=P(otf, 0);
	assert_float_equal(creal(P(o, npsf/2, npsf/2)), 1, 1e-12);
	for(long iy=1; iy<npsf; iy++){
		for(long ix=1; ix<npsf; ix++){
			if(ix==npsf/2&&iy==npsf/2) continue;
			long jx=labs(ix-npsf/2), jy=labs(iy-npsf/2);
			real ref=(real)(n-jx)*(n-jy)/(n*n)*exp(-0.3);
			assert_float_equal(creal(P(o, ix, iy)), ref, 1e-12);
		}
	}
	/*p/t/t removal on white noise reduces the effective variance: larger OTF*/
	ccell *otf2=NULL;
	genotf(&otf2, loc, NULL, NULL, NULL, 0, wvl, cov, 0, 0, npsf, npsf, 1, 1);
	assert_true(creal(P(P(otf2, 0), npsf/2+1, npsf/2))>creal(P(o, npsf/2+1, npsf/2)));
	ccellfree(otf2);
	/*tilt bias: OTF picks up a linear phase; it is Hermitian*/
	dzero(cov);
	const real tilt=1e-6;
	dmat *bias=mkopd(loc, tilt, 0, 0);
	genotf(&otf, loc, NULL, bias, NULL, 0, wvl, cov, 0, 0, npsf, npsf, 1, 0);
	o=P(otf, 0);
	for(long iy=1; iy<npsf; iy++){
		for(long ix=1; ix<npsf; ix++){
			long jx=ix-npsf/2, jy=iy-npsf/2;
			real mag=(labs(jx)<n&&labs(jy)<n)?(real)(n-labs(jx))*(n-labs(jy))/(n*n):0;
			assert_float_equal(cabs(P(o, ix, iy)), mag, 1e-12);
			if(mag>0){
				real ph=-2*M_PI/wvl*tilt*jx*dx;
				assert_float_equal(creal(P(o, ix, iy)), mag*cos(ph), 1e-12);
				assert_float_equal(cimag(P(o, ix, iy)), mag*sin(ph), 1e-12);
			}
			/*Hermitian*/
			long kx=npsf-ix, ky=npsf-iy;
			assert_float_equal(creal(P(o, ix, iy)), creal(P(o, kx, ky)), 1e-12);
			assert_float_equal(cimag(P(o, ix, iy)), -cimag(P(o, kx, ky)), 1e-12);
		}
	}
	/*zero bias is treated as no bias*/
	dzero(bias);
	genotf(&otf, loc, NULL, bias, NULL, 0, wvl, cov, 0, 0, npsf, npsf, 1, 0);
	assert_float_equal(cimag(P(P(otf, 0), npsf/2+1, npsf/2)), 0, 0);
	dfree(bias);
	/*Kolmogorov from r0 (cov=NULL): decreasing with separation, more with pttr*/
	genotf(&otf, loc, NULL, NULL, NULL, 0, wvl, NULL, 0.2, 30, npsf, npsf, 1, 0);
	genotf(&otf2, loc, NULL, NULL, NULL, 0, wvl, NULL, 0.2, 30, npsf, npsf, 1, 1);
	o=P(otf, 0);
	real o0=creal(P(o, npsf/2, npsf/2));
	real o1=creal(P(o, npsf/2+1, npsf/2));
	real o2=creal(P(o, npsf/2+2, npsf/2));
	assert_float_equal(o0, 1, 1e-12);
	assert_true(o1<0.75&&o2<o1&&o2>0);/*DL values are 0.75, 0.5*/
	assert_true(creal(P(P(otf2, 0), npsf/2+2, npsf/2))>o2);
	ccellfree(otf2);
	/*multiple subapertures with amplitude and area: full ones copy the OTF,
	  partial ones are computed, tiny ones are skipped*/
	const long nsa=3;
	dmat *amp=dnew(nloc, nsa);
	dmat *area=dnew(nsa, 1);
	for(long i=0; i<nloc; i++){
		P(amp, i, 0)=1;
		P(amp, i, 1)=(loc->locx[i]<0)?1:0;/*half*/
		P(amp, i, 2)=0;
	}
	P(area, 0)=1; P(area, 1)=0.5; P(area, 2)=0;
	dzero(cov);
	ccellfree(otf);
	genotf(&otf, loc, amp, NULL, area, 0.9, wvl, cov, 0, 0, npsf, npsf, nsa, 0);
	assert_int_equal(NX(otf), nsa);
	assert_float_equal(creal(P(P(otf, 0), npsf/2+1, npsf/2)), 0.75, 1e-12);
	/*half subaperture (2 of 4 columns): normalized overlap*/
	assert_float_equal(creal(P(P(otf, 1), npsf/2, npsf/2)), 1, 1e-12);
	assert_float_equal(creal(P(P(otf, 1), npsf/2+1, npsf/2)), 0.5, 1e-12);
	assert_float_equal(creal(P(P(otf, 1), npsf/2+2, npsf/2)), 0, 1e-12);
	assert_float_equal(creal(P(P(otf, 1), npsf/2, npsf/2+1)), 0.75, 1e-12);
	assert_float_equal(cabs(P(P(otf, 2), npsf/2, npsf/2)), 0, 0);
	ccellfree(otf);
	/*invalid input is an error and returns without output*/
	{
		dmat *area2=dnew(nsa+1, 1);
		signal_caught=0;
		genotf(&otf, loc, amp, NULL, area2, 0.9, wvl, cov, 0, 0, npsf, npsf, nsa, 0);
		assert_int_equal(signal_caught, SIGUSR2);
		assert_null(otf);
		dmat *bias2=dnew(nloc, 1);
		signal_caught=0;
		genotf(&otf, loc, amp, bias2, area, 0.9, wvl, cov, 0, 0, npsf, npsf, nsa, 0);
		assert_int_equal(signal_caught, SIGUSR2);
		assert_null(otf);
		signal_caught=0;
		dfree(bias2);
		dfree(area2);
	}
	dfree(amp); dfree(area);
	dfree(cov);
	locfree(loc);
}
/**
   mk2dcov(): spatially averaged covariance of an identity (white) covariance
   is a delta, of a constant covariance is constant; norm=0 sums instead of
   averaging; amplitude threshold excludes points.
*/
static void wfs_mk2dcov(void **state){
	(void)state;
	const long n=4;
	loc_t *loc=mksqloc_auto(n, n, 0.5, 0.5);
	const long nloc=loc->nloc;
	dmat *cov=dnew(nloc, nloc);
	daddI(cov, 2);
	dmat *c2=mk2dcov(loc, NULL, 0, cov, 1);
	assert_int_equal(NX(c2), 2*(n-1));
	const long c=n-1;
	assert_float_equal(P(c2, c, c), 2, 1e-14);
	assert_float_equal(dsum(c2), 2, 1e-14);
	dfree(c2);
	/*sum instead of average*/
	c2=mk2dcov(loc, NULL, 0, cov, 0);
	assert_float_equal(P(c2, c, c), 2*nloc, 1e-12);
	dfree(c2);
	/*constant covariance*/
	dset(cov, 1);
	c2=mk2dcov(loc, NULL, 0, cov, 1);
	for(long iy=1; iy<NY(c2); iy++){
		for(long ix=1; ix<NX(c2); ix++){
			assert_float_equal(P(c2, ix, iy), 1, 1e-14);
		}
	}
	dfree(c2);
	/*amplitude threshold: only the left half counts. The count at zero
	  separation is the number of points included (norm=0).*/
	dmat *amp=dnew(nloc, 1);
	for(long i=0; i<nloc; i++) P(amp, i)=loc->locx[i]<0?1:0;
	c2=mk2dcov(loc, amp, 0.5, cov, 0);
	assert_float_equal(P(c2, c, c), nloc/2, 1e-12);
	/*no pair with separation of 2 columns along x in the left half*/
	assert_float_equal(P(c2, c+2, c), 0, 0);
	dfree(c2); dfree(amp);
	dfree(cov); locfree(loc);
}
/*Setup for the subaperture image tests.*/
#define NXSA 8   /*points across a subaperture*/
#define NOTF 16  /*OTF/PSF size: embedding factor 2*/
#define NPIX 8   /*pixels across subaperture image*/
static const real WVL=1e-6;
/**
   Diffraction limited short exposure PSF of a square subaperture (via genotf
   and gensepsf) for nsa identical subapertures.
*/
static dccell *mksepsf(long nsa, dmat *wvl){
	loc_t *loc=mksqloc_auto(NXSA, NXSA, DSA/NXSA, DSA/NXSA);
	dmat *cov=dnew(loc->nloc, loc->nloc);
	cccell *otfs=cccellnew(1, 1);
	genotf(&P(otfs, 0), loc, NULL, NULL, NULL, 0, WVL, cov, 0, 0, NOTF, NOTF, nsa, 1);
	dccell *sepsf=NULL;
	gensepsf(&sepsf, otfs, NULL, (const dmat *)NULL, wvl, 0, 0);
	cellfree(otfs);
	dfree(cov);
	locfree(loc);
	return sepsf;
}
/**
   mkdtf(): dimensions, PSF sampling, nominal and pixel sampling matrix, and
   frequency vectors. smooth(): resampling of sodium profile preserves the sum.
*/
static void wfs_mkdtf(void **state){
	(void)state;
	dmat *wvls=dnew(2, 1);
	P(wvls, 0)=WVL; P(wvls, 1)=2*WVL;
	const real dtheta0=WVL/(DSA*2);
	dtf_t *dtf=mkdtf(wvls, DSA, 2, NOTF, NOTF, NPIX, NPIX, dtheta0, dtheta0, NULL, NULL, 0.3, NULL);
	assert_int_equal(dtf[0].nwvl, 2);
	for(int iwvl=0; iwvl<2; iwvl++){
		const real dtheta=P(wvls, iwvl)/(DSA*2);
		assert_float_equal(dtf[iwvl].dtheta, dtheta, 1e-20);
		assert_int_equal(dtf[iwvl].radpix, 0);
		assert_int_equal(NX(dtf[iwvl].nominal), 1);
		dsp *si=P(dtf[iwvl].si, 0);
		assert_int_equal(NX(si), NPIX*NPIX);
		assert_int_equal(NY(si), NOTF*NOTF);
		/*nominal is treated with an fftshift: zero frequency at the corner*/
		dmat *nom=P(dtf[iwvl].nominal, 0);
		const real pd=dtheta0/dtheta;
		assert_float_equal(P(nom, 0, 0), pd*pd, 1e-12);
		/*frequency vector*/
		const real dux=1./(dtheta*NOTF);
		assert_float_equal(cimag(P(dtf[iwvl].Ux, 1)), -2*M_PI*dux, 1e-6*dux);
		assert_float_equal(cimag(P(dtf[iwvl].Ux, NOTF-1)), 2*M_PI*dux, 1e-6*dux);
		assert_float_equal(creal(P(dtf[iwvl].Uy, 3)), 0, 0);
		/*sampling a constant PSF gives a constant image (bilinear weights sum to 1)*/
		dmat *one=dnew(NOTF*NOTF, 1); dset(one, 1);
		dmat *img=spmulvec(si, one);
		for(long i=0; i<PN(img); i++){
			assert_float_equal(P(img, i), 1, 1e-12);
		}
		dfree(img); dfree(one);
	}
	dtf_free(dtf);
	assert_null(dtf);
	/*pixel offset and rotation create one dtf per subaperture*/
	dmat *pixoffx=dnew(2, 1); P(pixoffx, 1)=0.5;
	dmat *pixoffy=dnew(2, 1);
	dcell *pixrot=dcellnew(1, 1); P(pixrot, 0)=dnew(2, 1); P(P(pixrot, 0), 1)=M_PI/2;
	dtf=mkdtf(wvls, DSA, 2, NOTF, NOTF, NPIX, NPIX, dtheta0, dtheta0, pixoffx, pixoffy, 0.3, pixrot);
	assert_int_equal(dtf[0].radpix, 1);
	assert_int_equal(NX(dtf[0].si), 2);
	assert_int_equal(NX(dtf[0].nominal), 2);
	/*90 degree rotation of square pixels: same nominal up to round off*/
	dmat *n0=P(dtf[0].nominal, 0), *n1=P(dtf[0].nominal, 1);
	real d=0;
	for(long iy=0; iy<NOTF; iy++){
		for(long ix=0; ix<NOTF; ix++){
			d=MAX(d, fabs(P(n0, ix, iy)-P(n1, iy, ix)));
		}
	}
	assert_float_equal(d, 0, 1e-12);
	dtf_free(dtf);
	dfree(pixoffx); dfree(pixoffy); dcellfree(pixrot);
	dfree(wvls);
	/*smooth: coarser resampling preserves the sum*/
	const long nh=1001;
	/*Two profile columns. The second is twice the first.*/
	dmat *prof=dnew(nh, 3);
	for(long i=0; i<nh; i++){
		P(prof, i, 0)=80e3+i*10;
		P(prof, i, 1)=exp(-pow((P(prof, i, 0)-90e3)/1000, 2));
		P(prof, i, 2)=2*P(prof, i, 1);
	}
	dmat *ps=smooth(prof, 100);
	assert_int_equal(NX(ps), 100);
	assert_int_equal(NY(ps), 3);
	assert_float_equal(P(ps, 1, 0)-P(ps, 0, 0), 100, 1e-8);
	assert_float_equal(dvecsum(PCOL(ps, 1), NX(ps)), dvecsum(PCOL(prof, 1), nh), 1e-8);
	assert_float_equal(dvecsum(PCOL(ps, 2), NX(ps)), dvecsum(PCOL(prof, 2), nh), 1e-8);
	for(long i=0; i<NX(ps); i++){
		assert_float_equal(P(ps, i, 2), 2*P(ps, i, 1), 1e-12);
	}
	/*peak stays near 90 km*/
	long imax=0;
	for(long i=0; i<NX(ps); i++){
		if(P(ps, i, 1)>P(ps, imax, 1)) imax=i;
	}
	assert_true(fabs(P(ps, imax, 0)-90e3)<=100);
	dfree(ps);
	/*fine resampling returns a reference*/
	ps=smooth(prof, 15);
	assert_ptr_equal(P(ps), P(prof));
	dfree(ps);
	dcell *profs=dcellnew(2, 1);
	P(profs, 0)=dref(prof); P(profs, 1)=dref(prof);
	dcell *pss=smooth_cell(profs, 100);
	assert_int_equal(NX(P(pss, 1)), 100);
	dcellfree(pss); dcellfree(profs);
	dfree(prof);
}
/**
   mketf(): a sodium layer at the focus range or a zero LLT distance gives no
   elongation (ETF equals the DTF nominal); a single layer elsewhere is a pure
   phase; frequencies normal to the elongation are unaffected; interpolation
   and direct sum agree.
*/
static void wfs_mketf(void **state){
	(void)state;
	dmat *wvls=dnew(1, 1); P(wvls, 0)=WVL;
	const real dtheta=WVL/(DSA*2);
	dtf_t *dtf=mkdtf(wvls, DSA, 2, NOTF, NOTF, NPIX, NPIX, dtheta, dtheta, NULL, NULL, 0.3, NULL);
	const real hs=90e3;
	const long nsa=3;
	dcell *srot=dcellnew(1, 1); P(srot, 0)=dnew(nsa, 1);
	dcell *srsa=dcellnew(1, 1); P(srsa, 0)=dnew(nsa, 1);
	P(P(srsa, 0), 0)=0;   /*on axis with LLT*/
	P(P(srsa, 0), 1)=5;   /*elongated along x*/
	P(P(srsa, 0), 2)=5; P(P(srot, 0), 2)=M_PI/2;/*elongated along y*/
	dmat *nom=P(dtf[0].nominal, 0);
	/*single layer at hs: no elongation*/
	dmat *na1=dnew(1, 2); P(na1, 0, 0)=hs; P(na1, 0, 1)=1;
	etf_t *etf=mketf(dtf, na1, 0, srot, srsa, hs, 0, 0, 0);
	assert_int_equal(etf[0].nwvl, 1);
	for(long isa=0; isa<nsa; isa++){
		cmat *e=P(etf[0].etf, isa);
		for(long i=0; i<PN(e); i++){
			assert_float_equal(creal(P(e, i)), P(nom, i), 1e-10);
			assert_float_equal(cimag(P(e, i)), 0, 1e-10);
		}
	}
	etf_free(etf);
	/*single layer away from hs: pure phase*/
	P(na1, 0, 0)=hs+2000;
	etf=mketf(dtf, na1, 0, srot, srsa, hs, 0, 0, 0);
	cmat *e1=P(etf[0].etf, 1);
	for(long i=0; i<PN(e1); i++){
		assert_float_equal(cabs(P(e1, i)), fabs(P(nom, i)), 1e-10);
	}
	assert_true(fabs(cimag(P(e1, 1, 0)))>1e-3*fabs(P(nom, 1, 0)));
	etf_free(etf);
	dfree(na1);
	/*broad gaussian profile: interpolation vs direct sum*/
	const long nh=201;
	dmat *na=dnew(nh, 2);
	for(long i=0; i<nh; i++){
		P(na, i, 0)=hs-10e3+i*100;
		P(na, i, 1)=exp(-pow((P(na, i, 0)-hs)/1500, 2));
	}
	dvecnormalize_sumabs(PCOL(na, 1), nh, 1);
	etf_t *etfi=mketf(dtf, na, 0, srot, srsa, hs, 0, 0, 0);
	etf_t *etfd=mketf(dtf, na, 0, srot, srsa, hs, 0, 0, 1);
	for(long isa=0; isa<nsa; isa++){
		cmat *ei=P(etfi[0].etf, isa);
		cmat *ed=P(etfd[0].etf, isa);
		real maxd=0;
		for(long i=0; i<PN(ei); i++){
			maxd=MAX(maxd, cabs(P(ei, i)-P(ed, i)));
			/*elongation only reduces the transfer function*/
			assert_true(cabs(P(ed, i))<=fabs(P(nom, i))*(1+1e-10)+1e-12);
		}
		assert_true(maxd<0.02);
		/*frequencies normal to the elongation are not affected*/
		for(long j=0; j<NOTF; j++){
			cmat *e=P(etfd[0].etf, isa);
			if(isa==1){
				assert_float_equal(creal(P(e, 0, j)), P(nom, 0, j), 1e-10);
				assert_float_equal(creal(P(ei, 0, j)), P(nom, 0, j), 1e-8);
			} else if(isa==2){
				assert_float_equal(creal(P(e, j, 0)), P(nom, j, 0), 1e-10);
			}
		}
	}
	/*elongated subapertures are attenuated along the elongation*/
	assert_true(cabs(P(P(etfd[0].etf, 1), 2, 0))<0.9*fabs(P(nom, 2, 0)));
	assert_true(cabs(P(P(etfd[0].etf, 2), 0, 2))<0.9*fabs(P(nom, 0, 2)));
	/*rsa=0 in interpolation mode: no elongation*/
	for(long i=0; i<NOTF*NOTF; i++){
		assert_float_equal(creal(P(P(etfi[0].etf, 0), i)), P(nom, i), 1e-12);
	}
	etf_free(etfi); etf_free(etfd);
	/*profile in a cell with a second column selected by icol (wrapped)*/
	dcell *nac=dcellnew(1, 1); P(nac, 0)=ddup(na);
	etf_t *etfc=mketf(dtf, nac, 2, srot, srsa, hs, 0, 0, 1);
	assert_float_equal(creal(P(P(etfc[0].etf, 1), 0)), P(nom, 0), 1e-10);
	etf_free(etfc);
	dcellfree(nac);
	dfree(na);
	dcellfree(srot); dcellfree(srsa);
	dtf_free(dtf);
	dfree(wvls);
}
/** Sum of a column of dcell entries. */
static real flip_diff(const dmat *a, int sign, int alongx){
	real d=0;
	for(long iy=0; iy<NY(a); iy++){
		for(long ix=0; ix<NX(a); ix++){
			real b=alongx?P(a, NX(a)-1-ix, iy):P(a, ix, NY(a)-1-iy);
			d=MAX(d, fabs(P(a, ix, iy)-sign*b));
		}
	}
	return d;
}
/**
   gensepsf() and gensei(): the subaperture image of a centered diffraction
   limited PSF is symmetric, its derivatives are antisymmetric and agree with
   finite differences of images shifted by goff; i0scale normalizes to
   subaperture area; rotated derivatives; pixel offset; shift2center; the
   OTF output.
*/
static void wfs_gensei(void **state){
	(void)state;
	const long nsa=2;
	dmat *wvl=dnew(1, 1); P(wvl, 0)=WVL;
	dccell *sepsf=mksepsf(nsa, wvl);
	assert_int_equal(NX(sepsf), 1);
	dmat *psf=P(P(sepsf, 0), 0);
	assert_int_equal(NX(psf), NOTF);
	/*PSF sums to 1 and peaks in the center*/
	assert_float_equal(dsum(psf), 1, 1e-10);
	assert_float_equal(P(psf, NOTF/2, NOTF/2), dmax(psf), 0);
	/*padded output*/
	{
		dccell *sepsf2=NULL;
		loc_t *loc=mksqloc_auto(NXSA, NXSA, DSA/NXSA, DSA/NXSA);
		dmat *cov=dnew(loc->nloc, loc->nloc);
		cccell *otfs=cccellnew(1, 1);
		genotf(&P(otfs, 0), loc, NULL, NULL, NULL, 0, WVL, cov, 0, 0, NOTF, NOTF, 1, 1);
		dmat *saa=dnew(1, 1); P(saa, 0)=0.5;
		gensepsf(&sepsf2, otfs, otfs, saa, wvl, NOTF+4, NOTF+4);
		assert_int_equal(NX(P(P(sepsf2, 0), 0)), NOTF+4);
		/*uplink OTF identical: PSF is the autoconvolution; still sums to saa*/
		assert_float_equal(dsum(P(P(sepsf2, 0), 0)), 0.5, 1e-10);
		dfree(saa); cellfree(otfs); dfree(cov); locfree(loc); cellfree(sepsf2);
	}
	const real dtheta=WVL/(DSA*2);
	dtf_t *dtf=mkdtf(wvl, DSA, 2, NOTF, NOTF, NPIX, NPIX, dtheta, dtheta, NULL, NULL, 0.2, NULL);
	dmat *siglev=dnew(1, 1); P(siglev, 0)=100;
	dmat *wvlwts=dnew(1, 1); P(wvlwts, 0)=1;
	dcell *i0=NULL, *gx=NULL, *gy=NULL;
	cccell *fotf=NULL;
	gensei(&i0, &gx, &gy, &fotf, sepsf, dtf, NULL, NULL, NULL, siglev, wvlwts, NULL, 0, 0);
	assert_int_equal(NX(i0), nsa);
	assert_int_equal(NY(i0), 1);
	dmat *im=P(i0, 0);
	assert_int_equal(NX(im), NPIX);
	real isum=dsum(im);
	assert_true(isum>70&&isum<100.0001);
	assert_float_equal(flip_diff(im, 1, 1), 0, 1e-8*isum);
	assert_float_equal(flip_diff(im, 1, 0), 0, 1e-8*isum);
	assert_float_equal(flip_diff(P(gx, 0), -1, 1), 0, 1e-8*isum/dtheta);
	assert_float_equal(flip_diff(P(gx, 0), 1, 0), 0, 1e-8*isum/dtheta);
	assert_float_equal(flip_diff(P(gy, 0), -1, 0), 0, 1e-8*isum/dtheta);
	assert_true(dmin(im)>-1e-8*isum);
	/*fotf at zero frequency is sum of PSF (1) times nominal, with FFT normalization*/
	assert_non_null(fotf);
	assert_float_equal(creal(P(P(P(fotf, 0), 0, 0), 0)), P(P(dtf[0].nominal, 0), 0)/(NOTF*NOTF), 1e-12);
	/*finite difference with goff*/
	const real delta=0.02*dtheta;
	dcell *goff=dcellnew(1, 1); P(goff, 0)=dnew(nsa*2, 1);
	dcell *ip=NULL, *im2=NULL;
	P(P(goff, 0), 0)=delta;
	gensei(&ip, NULL, NULL, NULL, sepsf, dtf, NULL, NULL, NULL, siglev, wvlwts, goff, 0, 0);
	P(P(goff, 0), 0)=-delta;
	gensei(&im2, NULL, NULL, NULL, sepsf, dtf, NULL, NULL, NULL, siglev, wvlwts, goff, 0, 0);
	dmat *fd=NULL;
	dadd(&fd, 0, P(ip, 0), 0.5/delta);
	dadd(&fd, 1, P(im2, 0), -0.5/delta);
	real gmax=dmaxabs(P(gx, 0));
	assert_true(gmax>0);
	assert_dmat_near(fd, P(gx, 0), 1e-3*gmax);
	/*second subaperture not shifted*/
	assert_dmat_near(P(ip, 1), P(i0, 1), 1e-8*isum);
	dfree(fd);
	/*y shift*/
	dzero(P(goff, 0));
	P(P(goff, 0), nsa)=delta;
	gensei(&ip, NULL, NULL, NULL, sepsf, dtf, NULL, NULL, NULL, siglev, wvlwts, goff, 0, 0);
	P(P(goff, 0), nsa)=-delta;
	gensei(&im2, NULL, NULL, NULL, sepsf, dtf, NULL, NULL, NULL, siglev, wvlwts, goff, 0, 0);
	dadd(&fd, 0, P(ip, 0), 0.5/delta);
	dadd(&fd, 1, P(im2, 0), -0.5/delta);
	assert_dmat_near(fd, P(gy, 0), 1e-3*gmax);
	dfree(fd);
	dcellfree(ip); dcellfree(im2); dcellfree(goff);
	/*shift2center on a centered PSF does nothing; reuses output arrays*/
	dcell *i0b=dcelldup(i0);
	gensei(&i0b, NULL, NULL, NULL, sepsf, dtf, NULL, NULL, NULL, siglev, wvlwts, NULL, 0, 1);
	assert_dmat_near(P(i0b, 0), P(i0, 0), 1e-6*isum);
	dcellfree(i0b);
	/*rotated derivative: 90 degrees maps (gx,gy) to (gy,-gx)*/
	dcell *gxr=NULL, *gyr=NULL;
	dcell *rot=dcellnew(1, 1); P(rot, 0)=dnew(nsa, 1); P(P(rot, 0), 0)=M_PI/2;
	gensei(NULL, &gxr, &gyr, NULL, sepsf, dtf, NULL, NULL, rot, siglev, wvlwts, NULL, 0, 0);
	assert_dmat_near(P(gxr, 0), P(gy, 0), 1e-8*gmax);
	dmat *mgx=ddup(P(gx, 0)); dscale(mgx, -1);
	assert_dmat_near(P(gyr, 0), mgx, 1e-8*gmax);
	assert_dmat_near(P(gxr, 1), P(gx, 1), 1e-8*gmax);
	dfree(mgx);
	dcellfree(gxr); dcellfree(gyr); dcellfree(rot);
	/*i0scale: normalized to subaperture area, gradients scaled alike*/
	dcell *saa=dcellnew(1, 1); P(saa, 0)=dnew(nsa, 1);
	P(P(saa, 0), 0)=1; P(P(saa, 0), 1)=0.5;
	dcell *i0s=NULL, *gxs=NULL;
	gensei(&i0s, &gxs, NULL, NULL, sepsf, dtf, NULL, saa, NULL, siglev, wvlwts, NULL, 1, 0);
	assert_float_equal(dsum(P(i0s, 0)), 1, 1e-10);
	assert_float_equal(dsum(P(i0s, 1)), 0.5, 1e-10);
	dmat *gxsc=ddup(P(gx, 0)); dscale(gxsc, 1./isum);
	assert_dmat_near(P(gxs, 0), gxsc, 1e-8*gmax/isum);
	dfree(gxsc);
	dcellfree(i0s); dcellfree(gxs); dcellfree(saa);
	/*two signal levels give two sets of images scaled accordingly*/
	dmat *siglev2=dnew(2, 1); P(siglev2, 0)=100; P(siglev2, 1)=50;
	dcell *i02=NULL;
	gensei(&i02, NULL, NULL, NULL, sepsf, dtf, NULL, NULL, NULL, siglev2, wvlwts, NULL, 0, 0);
	assert_int_equal(NY(i02), 2);
	dmat *half=ddup(P(i0, 0)); dscale(half, 0.5);
	assert_dmat_near(P(i02, 0, 1), half, 1e-10*isum);
	dfree(half); dfree(siglev2); dcellfree(i02);
	/*half pixel offset: image centered on a pixel*/
	dmat *pixoffx=dnew(1, 1); P(pixoffx, 0)=0.5;
	dtf_t *dtf2=mkdtf(wvl, DSA, 2, NOTF, NOTF, NPIX, NPIX, dtheta, dtheta, pixoffx, NULL, 0.2, NULL);
	dcell *i0o=NULL;
	gensei(&i0o, NULL, NULL, NULL, sepsf, dtf2, NULL, NULL, NULL, siglev, wvlwts, NULL, 0, 0);
	dmat *io=P(i0o, 0);
	for(long iy=0; iy<NPIX; iy++){
		for(long k=1; k<NPIX/2; k++){
			assert_float_equal(P(io, NPIX/2+k, iy), P(io, NPIX/2-k, iy), 1e-8*isum);
		}
	}
	assert_true(P(io, NPIX/2, NPIX/2)>P(im, NPIX/2, NPIX/2));
	dcellfree(i0o); dtf_free(dtf2); dfree(pixoffx);
	/*with an ETF of no elongation the result is identical*/
	{
		dcell *srot=dcellnew(1, 1); P(srot, 0)=dnew(nsa, 1);
		dcell *srsa=dcellnew(1, 1); P(srsa, 0)=dnew(nsa, 1);
		dmat *na=dnew(1, 2); P(na, 0, 0)=90e3; P(na, 0, 1)=1;
		etf_t *etf=mketf(dtf, na, 0, srot, srsa, 90e3, 0, 0, 0);
		dcell *i0e=NULL;
		gensei(&i0e, NULL, NULL, NULL, sepsf, dtf, etf, NULL, NULL, siglev, wvlwts, NULL, 0, 0);
		assert_dmat_near(P(i0e, 0), im, 1e-8*isum);
		dcellfree(i0e);
		etf_free(etf); dfree(na); dcellfree(srot); dcellfree(srsa);
	}
	dfree(siglev); dfree(wvlwts);
	cellfree(fotf);
	dcellfree(i0); dcellfree(gx); dcellfree(gy);
	dtf_free(dtf);
	cellfree(sepsf);
	dfree(wvl);
}
/** Gaussian spot of nx*nx pixels centered at (cx,cy) offset (in pixels) from center, and its derivatives wrt shift. */
static void mkspot(dmat **i0, dmat **gx, dmat **gy, long nx, real sigma, real sx, real sy, real flux, real pixtheta){
	*i0=dnew(nx, nx);
	if(gx) *gx=dnew(nx, nx);
	if(gy) *gy=dnew(nx, nx);
	const real c=(nx-1)*0.5;
	for(long iy=0; iy<nx; iy++){
		for(long ix=0; ix<nx; ix++){
			real x=ix-c-sx, y=iy-c-sy;
			real v=flux/(2*M_PI*sigma*sigma)*exp(-(x*x+y*y)/(2*sigma*sigma));
			P(*i0, ix, iy)=v;
			/*d/d(shift) = -d/dx; per unit angle*/
			if(gx) P(*gx, ix, iy)=v*x/(sigma*sigma)/pixtheta;
			if(gy) P(*gy, ix, iy)=v*y/(sigma*sigma)/pixtheta;
		}
	}
}
static void mtch_apply(real g[2], const dmat *mtche, const dmat *i0){
	g[0]=g[1]=0;
	for(long i=0; i<PN(i0); i++){
		g[0]+=P(mtche, 0, i)*P(i0, i);
		g[1]+=P(mtche, 1, i)*P(i0, i);
	}
}
/**
   mtch(): the matched filter maps gx, gy to unit x/y responses and i0 (plus
   background) to zero; it estimates small shifts; NEA equals the inverse
   Fisher matrix block; constraint columns map shifted images to +/-1 pixel;
   qe, background, rotation, FFT derivative and negative pixels.
*/
static void wfs_mtch(void **state){
	(void)state;
	const long nx=8;
	const real pixtheta=0.5e-6;
	const real sigma=1.3, flux=500;
	dmat *i0=NULL, *gx=NULL, *gy=NULL;
	mkspot(&i0, &gx, &gy, nx, sigma, 0.2, -0.1, flux, pixtheta);
	const real rne=3, bkgrnd=10, bkgrndc=9;
	dmat *mt=NULL, *nea=NULL;
	mtch(&mt, &nea, i0, gx, gy, NULL, NULL, NULL, bkgrnd, bkgrndc, rne, pixtheta, pixtheta, 0, 0, 0);
	assert_int_equal(NX(mt), 2);
	assert_int_equal(NY(mt), nx*nx);
	real g[2];
	mtch_apply(g, mt, gx);
	assert_float_equal(g[0], 1, 1e-8); assert_float_equal(g[1], 0, 1e-8);
	mtch_apply(g, mt, gy);
	assert_float_equal(g[0], 0, 1e-8); assert_float_equal(g[1], 1, 1e-8);
	/*The measured image at zero shift, i0 plus the residual background, gives zero.*/
	dmat *ib=ddup(i0); dadds(ib, bkgrnd-bkgrndc);
	mtch_apply(g, mt, ib);
	assert_float_equal(g[0], 0, 1e-8); assert_float_equal(g[1], 0, 1e-8);
	dfree(ib);
	/*small shift is measured (no residual background)*/
	{
		dmat *mt0=NULL;
		mtch(&mt0, NULL, i0, gx, gy, NULL, NULL, NULL, bkgrnd, bkgrnd, rne, pixtheta, pixtheta, 0, 0, 0);
		mtch_apply(g, mt0, i0);
		assert_float_equal(g[0], 0, 1e-8*pixtheta); assert_float_equal(g[1], 0, 1e-8*pixtheta);
		const real s=0.05;
		dmat *i1=NULL;
		mkspot(&i1, NULL, NULL, nx, sigma, 0.2+s, -0.1-s, flux, pixtheta);
		mtch_apply(g, mt0, i1);
		assert_float_equal(g[0], s*pixtheta, 0.05*s*pixtheta);
		assert_float_equal(g[1], -s*pixtheta, 0.05*s*pixtheta);
		dfree(i1); dfree(mt0);
	}
	/*NEA is the 2x2 block of the inverse of the weighted normal matrix*/
	{
		long n=nx*nx;
		dmat *G=dnew(n, 3), *W=dnew(n, 1);
		for(long i=0; i<n; i++){
			P(G, i, 0)=P(gx, i);
			P(G, i, 1)=P(gy, i);
			P(G, i, 2)=(P(i0, i)+(bkgrnd-bkgrndc))/pixtheta;
			P(W, i)=1./(rne*rne+bkgrnd+P(i0, i));
		}
		dmat *GW=ddup(G);
		for(long j=0; j<3; j++) for(long i=0; i<n; i++) P(GW, i, j)*=P(W, i);
		dmat *F=NULL;
		dmm(&F, 0, G, GW, "tn", 1);
		dmat *Fi=dinv(F);
		assert_int_equal(NX(nea), 2);
		assert_int_equal(NY(nea), 2);
		for(int j=0; j<2; j++){
			for(int i=0; i<2; i++){
				assert_float_equal(P(nea, i, j), P(Fi, i, j), 1e-6*P(Fi, 0, 0));
			}
		}
		/*symmetric positive definite*/
		assert_true(P(nea, 0, 0)>0&&P(nea, 1, 1)>0);
		assert_float_equal(P(nea, 0, 1), P(nea, 1, 0), 1e-12*P(nea, 0, 0));
		dfree(G); dfree(W); dfree(GW); dfree(F); dfree(Fi);
	}
	/*higher read noise increases NEA*/
	{
		dmat *nea2=NULL;
		mtch(NULL, &nea2, i0, gx, gy, NULL, NULL, NULL, bkgrnd, bkgrndc, 3*rne, pixtheta, pixtheta, 0, 0, 0);
		assert_true(P(nea2, 0, 0)>P(nea, 0, 0));
		/*qe of 1 gives the same NEA*/
		dmat *qe=dnew(nx*nx, 1); dset(qe, 1);
		mtch(NULL, &nea2, i0, gx, gy, qe, NULL, NULL, bkgrnd, bkgrndc, rne, pixtheta, pixtheta, 0, 0, 0);
		assert_dmat_near(nea2, nea, 1e-6*P(nea, 0, 0));
		/*lower qe increases NEA*/
		dset(qe, 0.5);
		mtch(NULL, &nea2, i0, gx, gy, qe, NULL, NULL, bkgrnd, bkgrndc, rne, pixtheta, pixtheta, 0, 0, 0);
		assert_true(P(nea2, 0, 0)>P(nea, 0, 0));
		dfree(qe); dfree(nea2);
	}
	/*background images*/
	{
		dmat *b2=dnew(nx, nx), *b2c=dnew(nx, nx);
		for(long i=0; i<nx*nx; i++){P(b2, i)=1+0.1*i; P(b2c, i)=0.9+0.1*i;}
		dmat *mt2=NULL;
		mtch(&mt2, NULL, i0, gx, gy, NULL, b2, b2c, bkgrnd, bkgrndc, rne, pixtheta, pixtheta, 0, 0, 0);
		dmat *ib2=ddup(i0);
		dadds(ib2, bkgrnd-bkgrndc);
		dadd(&ib2, 1, b2, 1);
		dadd(&ib2, 1, b2c, -1);
		mtch_apply(g, mt2, ib2);
		assert_float_equal(g[0], 0, 1e-8); assert_float_equal(g[1], 0, 1e-8);
		mtch_apply(g, mt2, gx);
		assert_float_equal(g[0], 1, 1e-8);
		dfree(ib2); dfree(mt2); dfree(b2); dfree(b2c);
	}
	/*constraints: shifted images by one pixel map to +/- one pixel (mtch2 uses unit pixel size)*/
	dmat *i0p=NULL, *gxp=NULL, *gyp=NULL;
	mkspot(&i0p, &gxp, &gyp, nx, sigma, 0.2, -0.1, flux, 1);
	for(int cr=1; cr<=3; cr++){
		dmat *mt2=NULL;
		mtch2(&mt2, NULL, i0p, gxp, gyp, cr);
		mtch_apply(g, mt2, gxp);
		assert_float_equal(g[0], 1, 1e-6); assert_float_equal(g[1], 0, 1e-6);
		mtch_apply(g, mt2, i0p);
		assert_float_equal(g[0], 0, 1e-6); assert_float_equal(g[1], 0, 1e-6);
		if(cr!=3){
			dmat *sh=dnew(nx, nx);
			for(long iy=0; iy<nx; iy++) for(long ix=1; ix<nx; ix++) P(sh, ix, iy)=P(i0p, ix-1, iy);
			mtch_apply(g, mt2, sh);
			assert_float_equal(g[0], 1, 1e-6); assert_float_equal(g[1], 0, 1e-6);
			dfree(sh);
		}
		if(cr!=2){
			dmat *sh=dnew(nx, nx);
			for(long iy=0; iy<nx-1; iy++) for(long ix=0; ix<nx; ix++) P(sh, ix, iy)=P(i0p, ix, iy+1);
			mtch_apply(g, mt2, sh);
			assert_float_equal(g[0], 0, 1e-6); assert_float_equal(g[1], -1, 1e-6);
			dfree(sh);
		}
		dfree(mt2);
	}
	dfree(i0p); dfree(gxp); dfree(gyp);
	/*constraint with residual background and physical pixel size: the shifted
	  image plus the residual background maps to one pixel.*/
	{
		dmat *mt2=NULL;
		mtch(&mt2, NULL, i0, gx, gy, NULL, NULL, NULL, bkgrnd, bkgrndc, rne, pixtheta, pixtheta, 0, 0, 1);
		dmat *sh=dnew(nx, nx);
		for(long iy=0; iy<nx; iy++) for(long ix=1; ix<nx; ix++) P(sh, ix, iy)=P(i0, ix-1, iy);
		dadds(sh, bkgrnd-bkgrndc);
		mtch_apply(g, mt2, sh);
		assert_float_equal(g[0], pixtheta, 1e-6*pixtheta); assert_float_equal(g[1], 0, 1e-6*pixtheta);
		dfree(sh); dfree(mt2);
	}
	/*radial gradients with pixel rotation: output rotated to x/y*/
	{
		const real th=0.4;
		dmat *mt2=NULL, *nea2=NULL;
		mtch(&mt2, &nea2, i0, gx, gy, NULL, NULL, NULL, bkgrnd, bkgrndc, rne, pixtheta, pixtheta, th, 1, 0);
		mtch_apply(g, mt2, gx);
		assert_float_equal(g[0], cos(th), 1e-8); assert_float_equal(g[1], sin(th), 1e-8);
		mtch_apply(g, mt2, gy);
		assert_float_equal(g[0], -sin(th), 1e-8); assert_float_equal(g[1], cos(th), 1e-8);
		/*rotation preserves the trace of NEA*/
		assert_float_equal(P(nea2, 0, 0)+P(nea2, 1, 1), P(nea, 0, 0)+P(nea, 1, 1), 1e-8*P(nea, 0, 0));
		dfree(mt2); dfree(nea2);
	}
	/*negative pixels are weighted by read noise only; identities still hold*/
	{
		dmat *in=ddup(i0);
		P(in, 0)=-1;
		dmat *mt2=NULL;
		mtch2(&mt2, NULL, in, gx, gy, 0);
		mtch_apply(g, mt2, gy);
		assert_float_equal(g[1], 1, 1e-8);
		dfree(mt2); dfree(in);
	}
	dfree(mt); dfree(nea);
	dfree(i0); dfree(gx); dfree(gy);
}
/**
   derive_by_fft(): derivative with respect to a shift of a periodic image
   (-d/dx) is exact for a sinusoid; mtch() without gx/gy uses it.
*/
static void wfs_derive_by_fft(void **state){
	(void)state;
	const long nx=16, ny=8;
	dmat *im=dnew(nx, ny);
	const real kx=2*M_PI*2/nx, ky=2*M_PI*1/ny;
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			P(im, ix, iy)=3+cos(kx*ix)+0.5*sin(ky*iy);
		}
	}
	dmat *dx=derive_by_fft(im, 0);
	dmat *dy=derive_by_fft(im, M_PI/2);
	dmat *d45=derive_by_fft(im, M_PI/4);
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			assert_float_equal(P(dx, ix, iy), kx*sin(kx*ix), 1e-10);
			assert_float_equal(P(dy, ix, iy), -0.5*ky*cos(ky*iy), 1e-10);
			assert_float_equal(P(d45, ix, iy), (P(dx, ix, iy)+P(dy, ix, iy))*sqrt(0.5), 1e-10);
		}
	}
	dfree(dx); dfree(dy); dfree(d45);
	dfree(im);
	/*matched filter from FFT derivative maps the FFT derivatives correctly*/
	dmat *i0=NULL;
	mkspot(&i0, NULL, NULL, 16, 1.5, 0, 0, 1000, 1);
	dmat *mt=NULL;
	mtch2(&mt, NULL, i0, NULL, NULL, 0);
	dmat *gxf=derive_by_fft(i0, 0);
	real g[2];
	mtch_apply(g, mt, gxf);
	assert_float_equal(g[0], 1, 1e-8); assert_float_equal(g[1], 0, 1e-8);
	dfree(gxf); dfree(mt); dfree(i0);
}
/**
   mtch_cell(): per subaperture matched filter and NEA; dim subapertures are
   disabled; sigratio scales NEA; sums of intensity; automatic constraint.
*/
static void wfs_mtch_cell(void **state){
	(void)state;
	const long nx=8, nsa=3;
	const real pixtheta=1e-6;
	dcell *i0s=dcellnew(nsa, 1), *gxs=dcellnew(nsa, 1), *gys=dcellnew(nsa, 1);
	const real flux[3]={1000, 800, 10};
	const real sig[3]={1.2, 2.5, 1.2};
	for(long isa=0; isa<nsa; isa++){
		mkspot(&P(i0s, isa), &P(gxs, isa), &P(gys, isa), nx, sig[isa], 0, 0, flux[isa], pixtheta);
	}
	dcell *mtche=NULL, *sanea=NULL;
	dmat *i0sum=NULL, *i0sumsum=NULL;
	mtch_cell(&mtche, &sanea, &i0sum, &i0sumsum, i0s, gxs, gys, NULL, NULL, NULL, 0, 0, 1, pixtheta, pixtheta, NULL, 0, 0, 1);
	assert_int_equal(NX(mtche), nsa);
	assert_int_equal(NX(P(sanea, 0)), nsa);
	assert_int_equal(NY(P(sanea, 0)), 3);
	for(long isa=0; isa<nsa; isa++){
		assert_float_equal(P(i0sum, isa), dsum(P(i0s, isa)), 1e-10);
	}
	assert_float_equal(P(i0sumsum, 0), dsum(i0sum), 1e-8);
	/*bright subapertures: unit response; NEA from mtch*/
	for(long isa=0; isa<2; isa++){
		real g[2];
		mtch_apply(g, P(mtche, isa), P(gxs, isa));
		assert_float_equal(g[0], 1, 1e-8);
		dmat *nea=NULL;
		mtch(NULL, &nea, P(i0s, isa), P(gxs, isa), P(gys, isa), NULL, NULL, NULL, 0, 0, 1, pixtheta, pixtheta, 0, 0, 0);
		assert_float_equal(P(P(sanea, 0), isa, 0), P(nea, 0, 0), 1e-10*P(nea, 0, 0));
		assert_float_equal(P(P(sanea, 0), isa, 1), P(nea, 1, 1), 1e-10*P(nea, 0, 0));
		assert_float_equal(P(P(sanea, 0), isa, 2), P(nea, 1, 0), 1e-10*P(nea, 0, 0));
		dfree(nea);
	}
	/*broader spot has larger NEA*/
	assert_true(P(P(sanea, 0), 1, 0)>P(P(sanea, 0), 0, 0));
	/*dim subaperture is disabled*/
	assert_float_equal(P(P(sanea, 0), 2, 0), 0, 0);
	assert_float_equal(dmaxabs(P(mtche, 2)), 0, 0);
	/*sigratio scales NEA, auto constraint, rotation*/
	dcell *sanea2=NULL;
	dcell *pixrots=dcellnew(1, 1); P(pixrots, 0)=dnew(nsa, 1);
	mtch_cell(NULL, &sanea2, NULL, NULL, i0s, gxs, gys, NULL, NULL, NULL, 0, 0, 1, pixtheta, pixtheta, pixrots, 0, 0, 2);
	assert_float_equal(P(P(sanea2, 0), 0, 0), 0.5*P(P(sanea, 0), 0, 0), 1e-10*P(P(sanea, 0), 0, 0));
	mtch_cell(&mtche, &sanea2, NULL, NULL, i0s, gxs, gys, NULL, NULL, NULL, 0, 0, 1, pixtheta, pixtheta, pixrots, 0, -1, 1);
	assert_true(P(P(sanea2, 0), 0, 0)>0);
	real g[2];
	mtch_apply(g, P(mtche, 1), P(gxs, 1));
	assert_float_equal(g[0], 1, 1e-6);
	dcellfree(pixrots);
	dcellfree(sanea2);
	dcellfree(mtche); dcellfree(sanea); dfree(i0sum); dfree(i0sumsum);
	dcellfree(i0s); dcellfree(gxs); dcellfree(gys);
}
/** Field amplitude (sqrt of PSF) with peak in center for pupil amp and phase (radian), embedded by 2.*/
static dmat *mkamp2(const dmat *amp, const dmat *phi){
	long n=NX(amp);
	cmat *wvf=cnew(n*2, n*2);
	cembed_wvf(wvf, P(phi), P(amp), n, n, TWOPI, 0);
	cfft2(wvf, -1);
	cfftshift(wvf);
	dmat *amp2=dnew(n*2, n*2);
	real amax=0;
	for(long i=0; i<PN(wvf); i++){
		P(amp2, i)=cabs(P(wvf, i));
		if(P(amp2, i)>amax) amax=P(amp2, i);
	}
	dscale(amp2, 1./amax);
	cfree(wvf);
	return amp2;
}
/**
   gerchberg_saxton(): recovers a tilt (which has no twin-image ambiguity)
   from the pupil amplitude and focal plane amplitude, with and without modal
   constraint; starting from the true phase it stays there.
*/
static void phase_retrieval_case(long n, real cx, real cy){
	dmat *amp=dnew(n, n);
	real sx=0, sy=0, sa=0;
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			real x=ix-cx, y=iy-cy;
			real r2=(x*x+y*y)/(n*n*0.25);
			P(amp, ix, iy)=(r2<=1)?(1-0.2*r2):0;
			sx+=P(amp, ix, iy)*ix; sy+=P(amp, ix, iy)*iy; sa+=P(amp, ix, iy);
		}
	}
	sx/=sa; sy/=sa;/*centroid: tilt modes are orthogonal to piston*/
	dmat *mode=dnew(n*n, 2);
	dmat *phi=dnew(n, n);
	const real a=2*M_PI*0.25/n;/*a quarter wave across the pupil*/
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			real x=ix-sx, y=iy-sy;
			P(mode, ix+iy*n, 0)=x;
			P(mode, ix+iy*n, 1)=y;
			P(phi, ix, iy)=a*x;
		}
	}
	dmat *wt=dref_reshape(amp, n*n, 1);
	dmat *rmod=dpinv(mode, wt);
	dmat *amp2=mkamp2(amp, phi);
	/*modal constrained*/
	dmat *phi1=NULL, *mphi=NULL;
	gerchberg_saxton(&phi1, &mphi, amp, amp2, NULL, mode, rmod, 2, 20);
	assert_int_equal(PN(mphi), 2);
	assert_float_equal(P(mphi, 0), a, 1e-3*a);
	assert_float_equal(P(mphi, 1), 0, 1e-3*a);
	assert_int_equal(NX(phi1), n);
	/*unconstrained: fit tilt to the recovered phase inside the pupil*/
	dmat *phi2=NULL;
	gerchberg_saxton(&phi2, NULL, amp, amp2, NULL, NULL, NULL, 2, 40);
	dmat *phi2v=dref_reshape(phi2, n*n, 1);
	dmat *m2=NULL;
	dmm(&m2, 0, rmod, phi2v, "nn", 1);
	assert_float_equal(P(m2, 0), a, 0.05*a);
	assert_float_equal(P(m2, 1), 0, 0.05*a);
	/*start from the truth: modal update is ~0 and the phase is kept*/
	dmat *phi3=NULL, *m3=NULL;
	gerchberg_saxton(&phi3, &m3, amp, amp2, phi, mode, rmod, 2, 5);
	assert_float_equal(P(m3, 0), 0, 0.05*a);
	assert_float_equal(P(m3, 1), 0, 0.05*a);
	/*field amplitude with the same size as the embedded pupil; no outputs requested*/
	gerchberg_saxton(NULL, NULL, amp, amp2, NULL, mode, rmod, 2, 2);
	dfree(phi3); dfree(m3);
	dfree(m2); dfree(phi2v); dfree(phi2);
	dfree(phi1); dfree(mphi);
	dfree(amp2); dfree(rmod); dfree(wt);
	dfree(phi); dfree(mode); dfree(amp);
}
static void wfs_phase_retrieval(void **state){
	(void)state;
	//The pupil is slightly decentered and tapered.
	phase_retrieval_case(16, 7.3, 7.6);
	/*A symmetric pupil has exactly zero field at the Nyquist frequency, which
	  must not give NaN.*/
	phase_retrieval_case(16, 7.5, 7.5);
	//odd size: the pupil field and the amplitude constraint are embedded consistently.
	phase_retrieval_case(15, 6.8, 7.1);
}
/** Average gradient of a phase defined on (nx+1)*(ny+1) corners for each nx*ny subaperture. */
static void corner_grad(dmat **gx, dmat **gy, const dmat *phi, real dx){
	long nx=NX(phi)-1, ny=NY(phi)-1;
	*gx=dnew(nx, ny);
	*gy=dnew(nx, ny);
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			P(*gx, ix, iy)=(P(phi, ix+1, iy)-P(phi, ix, iy)+P(phi, ix+1, iy+1)-P(phi, ix, iy+1))/(2*dx);
			P(*gy, ix, iy)=(P(phi, ix, iy+1)-P(phi, ix, iy)+P(phi, ix+1, iy+1)-P(phi, ix+1, iy))/(2*dx);
		}
	}
}
/** max |a-b-mean(a-b)| over entries where mask (if not NULL) is nonzero. */
static real diff_nopiston(const dmat *a, const dmat *b, const lmat *mask){
	real s=0; long n=0;
	for(long i=0; i<PN(a); i++){
		if(mask&&!P(mask, i)) continue;
		s+=P(a, i)-P(b, i); n++;
	}
	s/=n;
	real m=0;
	for(long i=0; i<PN(a); i++){
		if(mask&&!P(mask, i)) continue;
		m=MAX(m, fabs(P(a, i)-P(b, i)-s));
	}
	return m;
}
/**
   cure(): the cumulative reconstructor recovers a linear phase exactly (up to
   piston) on a full square grid and on a circular pupil with invalid
   subapertures, and a smooth phase approximately. cure_loc() and cure_loc2()
   work on subaperture coordinates and gradient vectors.
*/
static void wfs_cure(void **state){
	(void)state;
	const long n=8;
	const real dx=0.5;
	dmat *phi=dnew(n+1, n+1);
	for(long iy=0; iy<=n; iy++){
		for(long ix=0; ix<=n; ix++){
			P(phi, ix, iy)=0.3*ix*dx-0.2*iy*dx;
		}
	}
	dmat *gx=NULL, *gy=NULL;
	corner_grad(&gx, &gy, phi, dx);
	dmat *rec=NULL;
	cure(&rec, gx, gy, dx);
	assert_int_equal(NX(rec), n+1);
	assert_int_equal(NY(rec), n+1);
	assert_float_equal(diff_nopiston(rec, phi, NULL), 0, 1e-10);
	/*1d version along x is also exact*/
	dmat *rec1=NULL;
	cure1d(&rec1, gx, gy, dx);
	assert_float_equal(diff_nopiston(rec1, phi, NULL), 0, 1e-10);
	/*too small preallocated output is resized*/
	dmat *rec2=dnew(3, 3);
	cure1d(&rec2, gx, gy, dx);
	assert_int_equal(NX(rec2), n+1);
	assert_float_equal(diff_nopiston(rec2, phi, NULL), 0, 1e-10);
	dfree(rec2); dfree(rec1);
	/*smooth quadratic phase: small relative error*/
	dmat *phiq=dnew(n+1, n+1);
	for(long iy=0; iy<=n; iy++){
		for(long ix=0; ix<=n; ix++){
			real x=(ix-n*0.5)*dx, y=(iy-n*0.5)*dx;
			P(phiq, ix, iy)=x*x+0.5*x*y-0.3*y*y;
		}
	}
	dmat *gxq=NULL, *gyq=NULL;
	corner_grad(&gxq, &gyq, phiq, dx);
	dmat *recq=NULL;
	cure(&recq, gxq, gyq, dx);
	real rq=dmaxabs(phiq);
	assert_true(diff_nopiston(recq, phiq, NULL)<0.1*rq);
	dfree(gxq); dfree(gyq); dfree(recq); dfree(phiq);
	/*Partial pupils: invalid subapertures are NaN. type 0 cuts the upper right
	  corner so that every row and column starts at index 0, and type 2 is a
	  disk with rows starting at index>0: the linear phase is recovered
	  exactly. type 1 is an annulus with two segments per row in the middle:
	  only approximate.*/
	lmat *valid=lnew(n, n);
	lmat *cmask=lnew(n+1, n+1);
	for(int type=2; type>=0; type--){
		lzero(valid); lzero(cmask);
		dfree(gx); dfree(gy);
		corner_grad(&gx, &gy, phi, dx);
		for(long iy=0; iy<n; iy++){
			for(long ix=0; ix<n; ix++){
				real x=ix+0.5-n*0.5, y=iy+0.5-n*0.5;
				real r2=x*x+y*y;
				if(type==0){
					P(valid, ix, iy)=(ix+iy<=n+2)?1:0;
				} else if(type==2){
					P(valid, ix, iy)=(r2<=n*n*0.25)?1:0;
				} else{
					P(valid, ix, iy)=(r2<=n*n*0.25&&r2>2)?1:0;
				}
				if(!P(valid, ix, iy)){
					P(gx, ix, iy)=NAN;
					P(gy, ix, iy)=NAN;
				} else{
					P(cmask, ix, iy)=P(cmask, ix+1, iy)=P(cmask, ix, iy+1)=P(cmask, ix+1, iy+1)=1;
				}
			}
		}
		cure(&rec, gx, gy, dx);
		if(type!=1){
			assert_float_equal(diff_nopiston(rec, phi, cmask), 0, 1e-10);
		} else{
			assert_true(diff_nopiston(rec, phi, cmask)<0.05*dmaxabs(phi));
		}
	}
	/*cure_loc: subaperture coordinates of valid subapertures*/
	long nsa=lsum(valid);
	loc_t *saloc=locnew(nsa, dx, dx);
	dmat *grad=dnew(nsa*2, 1);
	long count=0;
	for(long iy=0; iy<n; iy++){
		for(long ix=0; ix<n; ix++){
			if(P(valid, ix, iy)){
				saloc->locx[count]=ix*dx;
				saloc->locy[count]=iy*dx;
				P(grad, count)=P(gx, ix, iy);
				P(grad, count+nsa)=P(gy, ix, iy);
				count++;
			}
		}
	}
	dmat *rl=NULL;
	cure_loc(&rl, grad, saloc);
	assert_non_null(rl->keywords);
	assert_true(NX(rl)>=n+1);
	/*the 2xnsa layout is also accepted*/
	dmat *gradt=dtrans(grad);
	reshape(gradt, 2, nsa);
	dmat *gradt2=dnew(2, nsa);
	for(long i=0; i<nsa; i++){
		P(gradt2, 0, i)=P(grad, i);
		P(gradt2, 1, i)=P(grad, i+nsa);
	}
	dmat *rl2=NULL;
	dmat *g2=dnew(nsa, 2);
	for(long i=0; i<nsa; i++){P(g2, i, 0)=P(grad, i); P(g2, i, 1)=P(grad, i+nsa);}
	dmat *g2t=dtrans(g2);/*2 x nsa*/
	cure_loc(&rl2, g2t, saloc);
	assert_dmat_near(rl2, rl, 1e-12);
	/*invalid gradient dimension is an error and returns*/
	{
		dmat *gbad=dnew(nsa, 1);
		dmat *rbad=NULL;
		signal_caught=0;
		cure_loc(&rbad, gbad, saloc);
		assert_int_equal(signal_caught, SIGUSR2);
		signal_caught=0;
		assert_null(rbad);
		dfree(gbad);
	}
	dfree(rl2); dfree(g2); dfree(g2t); dfree(gradt); dfree(gradt2);
	/*cure_loc2: phase at subaperture lower-left corners*/
	dmat *ps=NULL;
	cure_loc2(&ps, grad, saloc);
	assert_int_equal(PN(ps), nsa);
	dmat *pref=dnew(nsa, 1);
	for(long i=0; i<nsa; i++){
		P(pref, i)=0.3*saloc->locx[i]-0.2*saloc->locy[i];
	}
	assert_float_equal(diff_nopiston(ps, pref, NULL), 0, 1e-10);
	/*preallocated output is overwritten*/
	cure_loc2(&ps, grad, saloc);
	assert_float_equal(diff_nopiston(ps, pref, NULL), 0, 1e-10);
	dfree(pref); dfree(ps); dfree(rl); dfree(grad); locfree(saloc);
	lfree(valid); lfree(cmask);
	dfree(rec); dfree(gx); dfree(gy); dfree(phi);
}
UNITTEST_MAIN(
	cmocka_unit_test(wfs_mkg),
	cmocka_unit_test(wfs_mkg_ngs),
	cmocka_unit_test(wfs_mkz),
	cmocka_unit_test(wfs_dtf_otf),
	cmocka_unit_test(wfs_genotf),
	cmocka_unit_test(wfs_mk2dcov),
	cmocka_unit_test(wfs_mkdtf),
	cmocka_unit_test(wfs_mketf),
	cmocka_unit_test(wfs_gensei),
	cmocka_unit_test(wfs_mtch),
	cmocka_unit_test(wfs_derive_by_fft),
	cmocka_unit_test(wfs_mtch_cell),
	cmocka_unit_test(wfs_phase_retrieval),
	cmocka_unit_test(wfs_cure)
)
