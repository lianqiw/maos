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
  Tests of loc_t utilities, bilinear weighting (mkw), Zernike modes, and
  related helpers (port of test_loc.c, whose tests were all disabled).
*/

/** Piston removed variance o'*W0*o-(W1'*o)^2 of opd o. */
static real wvar(const dsp *W0, const dmat *W1, const dmat *o){
	return dspwdot(o, W0, o)-pow(ddot(W1, o), 2);
}
/** OPD on loc: a + b*x + c*y + d*(x^2+y^2). */
static dmat *mkopd(const loc_t *loc, real a, real b, real c, real d){
	dmat *o=dnew(loc->nloc, 1);
	for(long i=0; i<loc->nloc; i++){
		real x=loc->locx[i], y=loc->locy[i];
		P(o, i)=a+b*x+c*y+d*(x*x+y*y);
	}
	return o;
}
/**
   Check the common properties of W0, W1: W0 symmetric, rows and columns of
   W0 sum to W1, sum(W1)=1, and piston has zero variance. Then check the
   variance of tip, tilt, and focus against the expected values.
*/
static void check_w(const dsp *W0, const dmat *W1, const loc_t *loc, real vartt, real varfocus, real tol, real tolfocus){
	long nloc=loc->nloc;
	assert_int_equal(NX(W0), nloc);
	assert_int_equal(NY(W0), nloc);
	assert_int_equal(PN(W1), nloc);
	assert_float_equal(dsum(W1), 1, 1e-10);
	assert_true(dmin(W1)>=0);
	dsp *W0t=dsptrans(W0);
	dsp *diff=dspadd2(W0, 1, W0t, -1);
	dmat *Fd=NULL;
	dspfull(&Fd, diff, 'n', 1);
	assert_float_equal(dmaxabs(Fd), 0, 1e-14);
	dfree(Fd); dspfree(diff); dspfree(W0t);
	dmat *s1=dspsum(W0, 1);
	dmat *s2=dspsum(W0, 2);
	for(long i=0; i<nloc; i++){
		assert_float_equal(P(s1, i), P(W1, i), 1e-14);
		assert_float_equal(P(s2, i), P(W1, i), 1e-14);
	}
	dfree(s1); dfree(s2);
	dmat *o=mkopd(loc, 3, 0, 0, 0);
	assert_float_equal(wvar(W0, W1, o), 0, 1e-12);
	dfree(o);
	o=mkopd(loc, 1, 1, 0, 0);
	assert_float_equal(wvar(W0, W1, o), vartt, tol);
	dfree(o);
	o=mkopd(loc, 1, 0, 1, 0);
	assert_float_equal(wvar(W0, W1, o), vartt, tol);
	dfree(o);
	o=mkopd(loc, 0, 0, 0, 1);
	assert_float_equal(wvar(W0, W1, o), varfocus, tolfocus);
	dfree(o);
}
/**
   mkw_amp (old test_w1 and test_wploc) with a gray circular amplitude, a
   binary (rounded) circular amplitude, and uniform (NULL) amplitude on a
   square grid. For a disk of radius R the tip variance is R^2/4 and focus
   (r^2) variance is R^4/12; for a square of side L it is L^2/12.
*/
static void loc_mkw_amp(void **state){
	(void)state;
	const long nn=48;
	const real dx=1./16.;
	loc_t *loc=mksqloc_auto(nn, nn, dx, dx);
	const real R=1;
	//gray pixel circular amplitude
	dmat *amp=dnew(loc->nloc, 1);
	loc_circle_add(amp, loc, 0, 0, R, 0, 1);
	assert_float_equal(dmax(amp), 1, 1e-12);
	assert_float_equal(dmin(amp), 0, 0);
	dnormalize_sum(amp, 1);
	assert_float_equal(dsum(amp), 1, 1e-12);
	dsp *W0=NULL;
	dmat *W1=NULL;
	mkw_amp(loc, P(amp), &W0, &W1);
	check_w(W0, W1, loc, R*R/4, pow(R, 4)/12, 0.005, 0.005);
	dspfree(W0); dfree(W1);
	//binary amplitude from rounding the gray amplitude (old test_wploc)
	real ampmax=dmax(amp);
	for(long i=0; i<PN(amp); i++){
		P(amp, i)=round(P(amp, i)/ampmax);
	}
	dnormalize_sum(amp, 1);
	mkw_amp(loc, P(amp), &W0, &W1);
	check_w(W0, W1, loc, R*R/4, pow(R, 4)/12, 0.01, 0.01);
	//W1 is zero far outside of the pupil
	for(long i=0; i<loc->nloc; i++){
		if(hypot(loc->locx[i], loc->locy[i])>R+2*dx){
			assert_float_equal(P(W1, i), 0, 0);
		}
	}
	dspfree(W0); dfree(W1);
	//uniform amplitude: bilinear interpolation of x is exact, so the tip
	//variance is that of the continuous square of side (nn-1)*dx. The
	//interpolation error of the quadratic focus is about 1%.
	mkw_amp(loc, NULL, &W0, &W1);
	real L=(nn-1)*dx;
	real L4=pow(L, 4);
	//var(x^2+y^2) on a square of side L: 2*var(x^2)=2*(L^4/80-L^4/144)
	check_w(W0, W1, loc, L*L/12, 2*(L4/80-L4/144), 1e-6, 0.01);
	dspfree(W0); dfree(W1);
	dfree(amp);
	locfree(loc);
}
/** mkw_circular (old test_wcir): exact circular aperture integration, centered and offset circles. */
static void loc_mkw_circular(void **state){
	(void)state;
	const real dx=1./16.;
	loc_t *loc=mksqloc_auto(48, 48, dx, dx);
	const real R=1;
	dsp *W0=NULL;
	dmat *W1=NULL;
	mkw_circular(loc, 0, 0, R, &W0, &W1);
	check_w(W0, W1, loc, R*R/4, pow(R, 4)/12, 2e-4, 2e-4);
	//Only points within one grid spacing of the circle have weight.
	for(long i=0; i<loc->nloc; i++){
		if(hypot(loc->locx[i], loc->locy[i])>R+1.5*dx){
			assert_float_equal(P(W1, i), 0, 0);
		}
	}
	dspfree(W0); dfree(W1);
	//An offset smaller circle: tilt variance about its own center is R^2/4.
	const real R2=0.6, cx=0.2, cy=-0.1;
	mkw_circular(loc, cx, cy, R2, &W0, &W1);
	dmat *ox=mkopd(loc, 0, 1, 0, 0);
	assert_float_equal(wvar(W0, W1, ox), R2*R2/4, 2e-4);
	//W1 gives the area weighted average: mean of x is cx.
	assert_float_equal(ddot(W1, ox), cx, 1e-4);
	dmat *oy=mkopd(loc, 0, 0, 1, 0);
	assert_float_equal(ddot(W1, oy), cy, 1e-4);
	dfree(ox); dfree(oy);
	dspfree(W0); dfree(W1);
	locfree(loc);
}
/**
   Area integration of a circle (old test_int_reg and test_int_rand): pixel
   counting on a regular grid, Monte Carlo with uniform random points, gray
   pixel loc_circle_add, and the number of points in mkannloc.
*/
static void loc_area(void **state){
	(void)state;
	const real cr=0.3;
	const real exact=M_PI*cr*cr/4;//quarter circle in the unit square
	//regular grid
	{
		const int n=100;
		const real dx=1./n;
		long area=0;
		for(int i=0; i<n; i++){
			real y=dx*(i+0.5);
			for(int j=0; j<n; j++){
				real x=dx*(j+0.5);
				if(x*x+y*y<cr*cr) area++;
			}
		}
		real ar=(real)area/(n*n);
		//error is bounded by the area of the pixels crossed by the arc
		assert_float_equal(ar, exact, M_PI*cr/2*dx*dx*2);
	}
	//Monte Carlo
	{
		const long n=100000;
		rand_t rstat;
		seed_rand(&rstat, 1);
		long area=0;
		for(long i=0; i<n; i++){
			real x=randu(&rstat);
			real y=randu(&rstat);
			if(x*x+y*y<=cr*cr) area++;
		}
		real ar=(real)area/n;
		assert_float_equal(ar, exact, 6*sqrt(exact*(1-exact)/n));
	}
	//gray pixel circle and annulus on a loc
	{
		const real dx=1./32.;
		loc_t *loc=mksqloc_auto(64, 64, dx, dx);
		dmat *amp=dnew(loc->nloc, 1);
		const real R=0.8, Rin=0.3;
		loc_circle_add(amp, loc, 0.05, -0.1, R, 0, 1);
		assert_float_equal(dsum(amp)*dx*dx, M_PI*R*R, 2*M_PI*R*dx*0.05);
		dzero(amp);
		loc_circle_add(amp, loc, 0, 0, R, Rin, 2);
		assert_float_equal(dsum(amp)*dx*dx, 2*M_PI*(R*R-Rin*Rin), 2*2*M_PI*(R+Rin)*dx*0.05);
		assert_float_equal(dmax(amp), 2, 1e-12);
		//center of the annulus is empty
		for(long i=0; i<loc->nloc; i++){
			real r=hypot(loc->locx[i], loc->locy[i]);
			if(r<Rin-dx) assert_float_equal(P(amp, i), 0, 1e-12);
			if(r>Rin+dx&&r<R-dx) assert_float_equal(P(amp, i), 2, 1e-12);
		}
		//loc_circle_mul masks
		dmat *ones=dnew(loc->nloc, 1);
		dset(ones, 3);
		loc_circle_mul(ones, loc, 0, 0, R, 0, 1);
		assert_float_equal(dsum(ones)*dx*dx, 3*M_PI*R*R, 3*2*M_PI*R*dx*0.05);
		dfree(ones);
		//dimension mismatch is rejected without change
		dmat *bad=dnew(3, 1);
		loc_circle_add(bad, loc, 0, 0, R, 0, 1);
		assert_float_equal(dsum(bad), 0, 0);
		loc_circle_add(NULL, loc, 0, 0, R, 0, 1);
		dfree(bad);
		dfree(amp);
		locfree(loc);
	}
	//mkannloc: points within the circle (rows are filled continuously, so
	//an annulus includes its central hole)
	{
		const real D=2, dx=1./32.;
		loc_t *loc=mkannloc(D, 0, dx, 0.5);
		real area=loc->nloc*dx*dx;
		assert_float_equal(area, M_PI*D*D/4, M_PI*D*dx);
		for(long i=0; i<loc->nloc; i++){
			assert_true(hypot(loc->locx[i], loc->locy[i])<D/2+dx);
		}
		loc_t *loc2=mkannloc(D, 0.5, dx, 0.5);
		assert_int_equal(loc2->nloc, loc->nloc);
		locfree(loc2);
		locfree(loc);
	}
}
/**
   Square grid creation and rotation (old test_sqlocrot; mksqlocrot no longer
   exists and is ported as mksqloc + locrot): geometry helpers, radii are
   preserved, loc_angle measures the rotation, and rotating back restores.
*/
static void loc_sqlocrot(void **state){
	(void)state;
	const long nn=64;
	const real dx=1./64.;
	loc_t *loc=mksqloc(nn, nn, dx, dx, -nn/2*dx, -nn/2*dx);
	assert_int_equal(loc->nloc, nn*nn);
	assert_float_equal(loc->locx[0], -nn/2*dx, 0);
	assert_float_equal(loc->locy[0], -nn/2*dx, 0);
	assert_float_equal(loc->locx[nn-1], (nn/2-1)*dx, 1e-15);
	assert_float_equal(loc->locy[nn*nn-1], (nn/2-1)*dx, 1e-15);
	long nx=0, ny=0;
	loc_nxny(&nx, &ny, loc);
	assert_int_equal(nx, nn);
	assert_int_equal(ny, nn);
	assert_float_equal(loc_diam(loc), (nn-1)*dx, 1e-12);
	int ic=loccenter(loc);
	assert_float_equal(loc->locx[ic], 0, 1e-15);
	assert_float_equal(loc->locy[ic], 0, 1e-15);
	real xm, ym;
	locmean(&xm, &ym, loc);
	assert_float_equal(xm, -dx/2, 1e-12);
	assert_float_equal(ym, -dx/2, 1e-12);
	//mksqloc_auto has the same origin convention
	loc_t *loca=mksqloc_auto(nn, nn, dx, dx);
	for(long i=0; i<loc->nloc; i++){
		assert_float_equal(loca->locx[i], loc->locx[i], 0);
		assert_float_equal(loca->locy[i], loc->locy[i], 0);
	}
	locfree(loca);

	const real theta=M_PI/10;
	loc_t *loc2=locdup(loc);
	locrot(loc2, 0);//no-op
	assert_float_equal(loc_angle(loc, loc2), 0, 1e-12);
	locrot(loc2, theta);
	for(long i=0; i<loc->nloc; i++){
		assert_float_equal(hypot(loc2->locx[i], loc2->locy[i]), hypot(loc->locx[i], loc->locy[i]), 1e-14);
	}
	//locrot rotates the coordinate CCW, i.e., the points CW by theta.
	assert_float_equal(loc_angle(loc, loc2), -theta, 1e-10);
	assert_float_equal(loc2->locx[0], loc->locx[0]*cos(theta)+loc->locy[0]*sin(theta), 1e-14);
	assert_float_equal(loc2->locy[0], -loc->locx[0]*sin(theta)+loc->locy[0]*cos(theta), 1e-14);
	//rotating back restores the original
	locrot(loc2, -theta);
	for(long i=0; i<loc->nloc; i++){
		assert_float_equal(loc2->locx[i], loc->locx[i], 1e-14);
		assert_float_equal(loc2->locy[i], loc->locy[i], 1e-14);
	}
	//shift
	locshift(loc2, 0.1, -0.2);
	locmean(&xm, &ym, loc2);
	assert_float_equal(xm, -dx/2+0.1, 1e-12);
	assert_float_equal(ym, -dx/2-0.2, 1e-12);
	assert_float_equal(loc_angle(NULL, loc2), 0, 0);
	locfree(loc2);
	locfree(loc);
}
/** A loc of the points of loc with radius in [rin, rout]. */
static loc_t *annulus(const loc_t *loc, real rin, real rout){
	loc_t *out=locnew(loc->nloc, loc->dx, loc->dy);
	long count=0;
	for(long i=0; i<loc->nloc; i++){
		real r=hypot(loc->locx[i], loc->locy[i]);
		if(r>=rin&&r<=rout){
			out->locx[count]=loc->locx[i];
			out->locy[count]=loc->locy[i];
			count++;
		}
	}
	locresize(out, count);
	return out;
}
/** Bilinear interpolation of opd=1+x-2y (exact for bilinear) from xloc must equal the same function on ploc. */
static void check_h(const dsp *H, const loc_t *xloc, const loc_t *ploc, int trans){
	dmat *ox=mkopd(xloc, 1, 1, -2, 0);
	dmat *op=dnew(ploc->nloc, 1);
	dspmv(op, H, ox, trans?'t':'n', 1);
	dmat *opref=mkopd(ploc, 1, 1, -2, 0);
	assert_dmat_near(op, opref, 1e-12);
	dfree(op); dfree(opref); dfree(ox);
}
/**
   loc_reduce_sp and loc_reduce_spcell (old test_loc_reduce_sp and
   test_loc_reduce_spcell) using a ray tracing operator H from a square xloc
   to an annular ploc made by mkh: uncoupled xloc points are removed with the
   matching rows/columns of H, and H still interpolates correctly. With
   cont=1 only points at the ends of each row are removed, so the uncoupled
   points in the annulus hole are kept.
*/
static void loc_reduce_sparse(void **state){
	(void)state;
	const real dx=1./16.;
	loc_t *grid=mksqloc_auto(32, 32, dx/2, dx/2);
	loc_t *ploc=annulus(grid, 0.3, 0.6);
	locfree(grid);
	for(int cont=0; cont<2; cont++){
		for(int dim=1; dim<=2; dim++){
			loc_t *xloc=mksqloc_auto(24, 24, dx, dx);
			long nloc0=xloc->nloc;
			dsp *H=mkh(ADDR((propdata_t){.locin=xloc, .locout=ploc}));
			assert_int_equal(NX(H), ploc->nloc);
			assert_int_equal(NY(H), nloc0);
			check_h(H, xloc, ploc, 0);
			//number of xloc points coupled to ploc
			dmat *cs=dspsumabs(H, 1);
			long ncoupled=0;
			for(long i=0; i<nloc0; i++) ncoupled+=P(cs, i)>EPS;
			dfree(cs);
			if(dim==1){//loc along rows
				dsp *Ht=dsptrans(H);
				loc_reduce_sp(xloc, Ht, 1, cont);
				assert_int_equal(NX(Ht), xloc->nloc);
				check_h(Ht, xloc, ploc, 1);
				dspfree(Ht);
			} else{//loc along columns
				loc_reduce_sp(xloc, H, 2, cont);
				assert_int_equal(NY(H), xloc->nloc);
				check_h(H, xloc, ploc, 0);
			}
			if(cont){
				assert_true(xloc->nloc>ncoupled);
				assert_true(xloc->nloc<nloc0);
			} else{
				assert_int_equal(xloc->nloc, ncoupled);
			}
			dspfree(H);
			locfree(xloc);
		}
	}
	//spcell version: two output locs share the same xloc.
	loc_t *ploc2=locdup(ploc);
	locshift(ploc2, 0.1, 0);
	for(int dim=1; dim<=2; dim++){
		loc_t *xloc=mksqloc_auto(24, 24, dx, dx);
		dsp *H1=mkh(ADDR((propdata_t){.locin=xloc, .locout=ploc}));
		dsp *H2=mkh(ADDR((propdata_t){.locin=xloc, .locout=ploc2}));
		//union of coupled points
		dmat *cs1=dspsumabs(H1, 1);
		dmat *cs2=dspsumabs(H2, 1);
		long ncoupled=0, n1=0;
		for(long i=0; i<xloc->nloc; i++){
			ncoupled+=(P(cs1, i)+P(cs2, i))>EPS;
			n1+=P(cs1, i)>EPS;
		}
		dfree(cs1); dfree(cs2);
		assert_true(ncoupled>n1);
		dspcell *Hc;
		if(dim==2){
			Hc=dspcellnew(2, 1);
			P(Hc, 0)=H1; P(Hc, 1)=H2;
		} else{
			Hc=dspcellnew(1, 2);
			P(Hc, 0)=dsptrans(H1); P(Hc, 1)=dsptrans(H2);
			dspfree(H1); dspfree(H2);
		}
		loc_reduce_spcell(xloc, Hc, dim, 0);
		assert_int_equal(xloc->nloc, ncoupled);
		check_h(P(Hc, 0), xloc, ploc, dim==1);
		check_h(P(Hc, 1), xloc, ploc2, dim==1);
		dspcellfree(Hc);
		locfree(xloc);
	}
	locfree(ploc2);
	//NULL input is ignored
	loc_reduce_sp(NULL, NULL, 1, 0);
	loc_reduce_spcell(NULL, NULL, 1, 0);
	locfree(ploc);
}
/**
   Zernike modes (old test_zernike) on a circular loc: piston, tip/tilt, and
   focus have the analytic form; the modes are nearly orthonormal on the grid
   and exactly so after Gram-Schmidt; selecting a single mode or radial modes
   gives the matching columns of the full set.
*/
static void loc_zernike(void **state){
	(void)state;
	const real D=1;
	loc_t *loc=mkannloc(D, 0, 1./64., 0.5);
	const long nloc=loc->nloc;
	const int rmax=10;
	const int nmod=(rmax+1)*(rmax+2)/2;
	//negative D: use D without clipping points outside of it.
	dmat *mod=zernike(loc, -D, 0, rmax, 0);
	assert_int_equal(NX(mod), nloc);
	assert_int_equal(NY(mod), nmod);
	const real R1=2./D;
	for(long i=0; i<nloc; i++){
		real x=loc->locx[i]*R1, y=loc->locy[i]*R1;
		assert_float_equal(P(mod, i, 0), 1, 1e-12);//piston
		assert_float_equal(P(mod, i, 1), 2*x, 1e-12);//Z2: x tilt
		assert_float_equal(P(mod, i, 2), 2*y, 1e-12);//Z3: y tilt
		assert_float_equal(P(mod, i, 3), sqrt(3.)*(2*(x*x+y*y)-1), 1e-12);//Z4: focus
	}
	//Nearly orthonormal: deviation (0.07 measured) is due to the discrete sampling.
	dmat *MM=NULL;
	dmm(&MM, 0, mod, mod, "tn", 1./nloc);
	dmat *Id=dnew(nmod, nmod);
	daddI(Id, 1);
	assert_dmat_near(MM, Id, 0.1);
	//Low order modes are better sampled (0.038 measured).
	dmat *MM15=dsub(MM, 0, 15, 0, 15);
	dmat *Id15=dsub(Id, 0, 15, 0, 15);
	assert_dmat_near(MM15, Id15, 0.05);
	dfree(MM15); dfree(Id15);
	//Single mode selection (flag<0) matches the column of Noll index -flag.
	for(int k=1; k<=nmod; k++){
		dmat *zk=zernike(loc, -D, 0, 0, -k);
		assert_int_equal(NY(zk), 1);
		dmat *col=dsub(mod, 0, nloc, k-1, 1);
		assert_dmat_near(zk, col, 0);
		dfree(col);
		dfree(zk);
	}
	//Radial only (flag>0): orders 2,4,...,rmax are the m=0 Noll modes 4, 11, 22, 37, 56.
	dmat *rad=zernike(loc, -D, 1, rmax, 1);
	assert_int_equal(NY(rad), rmax/2);
	const int noll_radial[]={4, 11, 22, 37, 56};
	for(int i=0; i<rmax/2; i++){
		dmat *col=dsub(mod, 0, nloc, noll_radial[i]-1, 1);
		dmat *rc=dsub(rad, 0, nloc, i, 1);
		assert_dmat_near(rc, col, 0);
		dfree(col); dfree(rc);
	}
	dfree(rad);
	//rmin skips lower orders
	dmat *mod2=zernike(loc, -D, 2, 3, 0);
	assert_int_equal(NY(mod2), 7);
	dmat *col=dsub(mod, 0, nloc, 3, 7);
	assert_dmat_near(mod2, col, 0);
	dfree(col);
	dfree(mod2);

	//Gram-Schmidt makes the modes exactly orthonormal without changing piston.
	dmat *gs=ddup(mod);
	dgramschmidt(gs, NULL);
	dfree(MM);
	dmm(&MM, 0, gs, gs, "tn", 1./nloc);
	assert_dmat_near(MM, Id, 1e-10);
	dmat *c0=dsub(gs, 0, nloc, 0, 1);
	dmat *m0=dsub(mod, 0, nloc, 0, 1);
	assert_dmat_near(c0, m0, 1e-12);
	dfree(c0); dfree(m0);
	//It is a small correction for the nearly orthonormal modes.
	dmat *gsd=ddup(gs);
	dadd(&gsd, 1, mod, -1);
	assert_true(sqrt(dsumsq(gsd)/nloc)<0.3*nmod);
	dfree(gsd);
	//Weighted Gram-Schmidt: <Mod|amp|Mod>/sum(amp)=I
	dmat *amp=dnew(nloc, 1);
	for(long i=0; i<nloc; i++) P(amp, i)=1+hypot(loc->locx[i], loc->locy[i]);
	dmat *gsw=ddup(mod);
	dgramschmidt(gsw, P(amp));
	dmat *wgs=ddup(gsw);
	for(long iy=0; iy<nmod; iy++){
		for(long ix=0; ix<nloc; ix++){
			P(wgs, ix, iy)*=P(amp, ix);
		}
	}
	dfree(MM);
	dmm(&MM, 0, gsw, wgs, "tn", 1./dsum(amp));
	assert_dmat_near(MM, Id, 1e-10);
	dfree(wgs); dfree(gsw); dfree(amp);
	dfree(gs);
	dfree(MM); dfree(Id);
	dfree(mod);
	//invalid input
	assert_null(zernike(NULL, D, 0, 2, 0));
	assert_null(zernike(loc, D, 3, 2, 0));
	locfree(loc);
}
/**
   dembed_locstat (old test_embed): embedding values on a circular loc into a
   2-D array places each value at the grid position of its point, centered in
   a larger preallocated array, with alpha/beta scaling, and the reverse
   operation extracts the values back.
*/
static void loc_embed_locstat(void **state){
	(void)state;
	const real dx=1./32.;
	loc_t *loc=mkannloc(2, 0, dx, 0);
	loc_create_stat(loc);
	locstat_t *st=loc->stat;
	long nx, ny;
	loc_nxny(&nx, &ny, loc);
	assert_int_equal(st->nx, nx);
	assert_int_equal(st->ny, ny);
	dmat *val=dnew(loc->nloc, 1);
	for(long i=0; i<loc->nloc; i++) P(val, i)=i+1;
	//allocated to the exact size
	dmat *opd=NULL;
	dembed_locstat(&opd, 0, loc, P(val), 1, 0);
	assert_int_equal(NX(opd), nx);
	assert_int_equal(NY(opd), ny);
	assert_float_equal(dsum(opd), dsum(val), 1e-6);
	for(long i=0; i<loc->nloc; i++){
		long ix=lround((loc->locx[i]-st->xmin)/dx);
		long iy=lround((loc->locy[i]-st->ymin)/dx);
		assert_float_equal(P(opd, ix, iy), P(val, i), 0);
	}
	//reverse: extract back
	dmat *back=dnew(loc->nloc, 1);
	dembed_locstat(&opd, 1, loc, P(back), 0, 1);
	assert_dmat_near(back, val, 0);
	//reverse with accumulation: back=back*beta+alpha*opd
	dembed_locstat(&opd, 2, loc, P(back), 1, 1);
	dmat *val3=ddup(val);
	dscale(val3, 3);
	assert_dmat_near(back, val3, 0);
	dfree(val3);
	//forward with alpha: opd=opd*alpha+val*beta
	dembed_locstat(&opd, 2, loc, P(val), -1, 0);
	assert_float_equal(dsum(opd), dsum(val), 1e-6);
	dfree(opd);
	//larger preallocated array: centered
	const long nbig=nx+5;
	opd=dnew(nbig, nbig);
	dembed_locstat(&opd, 0, loc, P(val), 1, 0);
	assert_int_equal(NX(opd), nbig);
	long off=(nbig-nx+1)/2;
	for(long i=0; i<loc->nloc; i++){
		long ix=lround((loc->locx[i]-st->xmin)/dx)+off;
		long iy=lround((loc->locy[i]-st->ymin)/dx)+off;
		assert_float_equal(P(opd, ix, iy), P(val, i), 0);
	}
	assert_float_equal(dsum(opd), dsum(val), 1e-6);
	dfree(opd);
	//NULL input is ignored
	opd=NULL;
	dembed_locstat(&opd, 0, loc, NULL, 1, 0);
	assert_null(opd);
	dfree(back);
	dfree(val);
	locfree(loc);
}
/** nextpow2 (old disabled loop in main): smallest power of 2 not less than n. */
static void loc_nextpow2(void **state){
	(void)state;
	for(long i=1; i<5000; i++){
		long p=nextpow2(i);
		assert_true(p>=i);
		assert_true(p<2*i);
		assert_int_equal(p&(p-1), 0);
	}
	assert_int_equal(nextpow2(1), 1);
	assert_int_equal(nextpow2(2), 2);
	assert_int_equal(nextpow2(3), 4);
	assert_int_equal(nextpow2(1024), 1024);
	assert_int_equal(nextpow2(1025), 2048);
	assert_int_equal(nextpow2(1L<<40), 1L<<40);
	assert_int_equal(nextpow2((1L<<40)+1), 1L<<41);
}
/** Assert that parse_poly(str) gives the terms (coefficient, power of x, power of y) in ref. */
static void check_poly(const char *str, int n, const real ref[][3]){
	dmat *c=parse_poly(str);
	assert_non_null(c);
	assert_int_equal(NX(c), 3);
	assert_int_equal(NY(c), n);
	for(int i=0; i<n; i++){
		for(int j=0; j<3; j++){
			assert_float_equal(P(c, j, i), ref[i][j], 1e-15);
		}
	}
	dfree(c);
}
/**
   parse_poly (old test_parse_poly, which parsed argv) on a set of
   polynomial strings, and loctransform using it. Invalid strings return NULL.
*/
static void loc_parse_poly(void **state){
	(void)state;
	check_poly("x", 1, (const real[][3]){{1, 1, 0}});
	check_poly("-2.5", 1, (const real[][3]){{-2.5, 0, 0}});
	check_poly("-x+y", 2, (const real[][3]){{-1, 1, 0}, {1, 0, 1}});
	check_poly("1.1*x-0.3*y+2", 3, (const real[][3]){{1.1, 1, 0}, {-0.3, 0, 1}, {2, 0, 0}});
	check_poly("x^2 + 2*y^3", 2, (const real[][3]){{1, 2, 0}, {2, 0, 3}});
	check_poly("0.2*x^2y", 1, (const real[][3]){{0.2, 2, 1}});
	check_poly("2*xy", 1, (const real[][3]){{2, 1, 1}});
	check_poly("x y", 1, (const real[][3]){{1, 1, 1}});
	check_poly("x*y", 1, (const real[][3]){{1, 1, 1}});
	check_poly("0.1+x+0.2*x^2*y", 3, (const real[][3]){{0.1, 0, 0}, {1, 1, 0}, {0.2, 2, 1}});//example in the documentation
	check_poly("3*x*2", 1, (const real[][3]){{6, 1, 0}});
	check_poly("1e-3*x", 1, (const real[][3]){{1e-3, 1, 0}});
	//more terms than the initial allocation of 5
	check_poly("1+x+y+x^2+y^2+x y+3", 7, (const real[][3]){{1, 0, 0}, {1, 1, 0}, {1, 0, 1}, {1, 2, 0}, {1, 0, 2}, {1, 1, 1}, {3, 0, 0}});
	assert_null(parse_poly(NULL));
	assert_null(parse_poly("x+z"));//invalid strings return NULL instead of looping
	assert_null(parse_poly("x^a"));
	assert_null(parse_poly("*x"));

	//loctransform with a polynomial string: shift, scale, and a quadratic term
	loc_t *loc=mksqloc_auto(8, 8, 0.1, 0.1);
	loc_t *loc2=loctransform(loc, "x+0.1;1.5*y-0.2+0.5*x^2y");
	assert_int_equal(loc2->nloc, loc->nloc);
	for(long i=0; i<loc->nloc; i++){
		real x=loc->locx[i], y=loc->locy[i];
		assert_float_equal(loc2->locx[i], x+0.1, 1e-14);
		assert_float_equal(loc2->locy[i], 1.5*y-0.2+0.5*x*x*y, 1e-14);
	}
	locfree(loc2);
	//loctransform2 with an n*4 coefficient matrix [px py cx cy]
	dmat *coeff=dnew(2, 4);
	P(coeff, 0, 0)=1; P(coeff, 0, 1)=0; P(coeff, 0, 2)=2; P(coeff, 0, 3)=0;//2x for new x
	P(coeff, 1, 0)=0; P(coeff, 1, 1)=1; P(coeff, 1, 2)=0; P(coeff, 1, 3)=-1;//-y for new y
	loc2=loctransform2(loc, coeff);
	for(long i=0; i<loc->nloc; i++){
		assert_float_equal(loc2->locx[i], 2*loc->locx[i], 1e-14);
		assert_float_equal(loc2->locy[i], -loc->locy[i], 1e-14);
	}
	locfree(loc2);
	dfree(coeff);
	assert_null(loctransform(NULL, "x;y"));
	assert_null(loctransform(loc, NULL));
	locfree(loc);
}

/** loc_convert with keywords, loc_keywords, loc_dxdy, locref, lochash, mk1dloc, mk1dloc_vec, loc2mat. */
static void loc_convert_keywords(void **state){
	(void)state;
	loc_t *loc=mksqloc(4, 3, 0.5, 0.25, -1, 2);
	loc->ht=1000; loc->iac=0.3; loc->dratio=1.2; loc->aoi=0.05;
	loc->dmat->make_keywords((cell *)loc);
	assert_non_null(loc->keywords);
	dmat *A=ddup(loc->dmat);
	free(A->keywords);
	A->keywords=strdup(loc->keywords);
	loc_t *loc2=loc_convert(A);
	assert_int_equal(loc2->id, M_LOC);
	assert_int_equal(loc2->nloc, 12);
	assert_float_equal(loc2->dx, 0.5, 1e-15);
	assert_float_equal(loc2->dy, 0.25, 1e-15);
	assert_float_equal(loc2->ht, 1000, 0);
	assert_float_equal(loc2->iac, 0.3, 1e-15);
	assert_float_equal(loc2->dratio, 1.2, 1e-15);
	assert_float_equal(loc2->aoi, 0.05, 1e-15);
	assert_float_equal(loc2->locy[11], 2.5, 1e-15);
	locfree(loc2);
	//without keywords dx, dy are derived from the data
	A=ddup(loc->dmat);
	free(A->keywords); A->keywords=NULL;
	loc2=loc_convert(A);
	assert_float_equal(loc2->dx, 0.5, 1e-12);
	assert_float_equal(loc2->dy, 0.25, 1e-12);
	//a wrong dx is corrected with a warning
	loc2->dx=0.3; loc2->dy=0;
	loc_dxdy(loc2);
	assert_float_equal(loc2->dx, 0.5, 1e-12);
	assert_float_equal(loc2->dy, 0.25, 1e-12);
	locfree(loc2);
	assert_null(loc_convert(NULL));
	//locref shares coordinates and copies the metadata
	loc_t *lr=locref(loc);
	assert_ptr_equal(lr->locx, loc->locx);
	assert_ptr_equal(lr->locy, loc->locy);
	assert_float_equal(lr->ht, 1000, 0);
	assert_float_equal(lr->aoi, 0.05, 0);
	assert_int_equal(lochash(lr, 3), lochash(loc, 3));
	assert_int_equal(lochash(NULL, 3), 3);
	assert_int_equal(lochash(loc, 3), dhash(loc->dmat, 3));
	locfree(lr);
	assert_null(locref(NULL));
	//loc2mat with and without piston
	dmat *m=loc2mat(loc, 1);
	assert_int_equal(NY(m), 3);
	for(long i=0; i<12; i++){
		assert_float_equal(P(m, i, 0), 1, 0);
		assert_float_equal(P(m, i, 1), loc->locx[i], 0);
		assert_float_equal(P(m, i, 2), loc->locy[i], 0);
	}
	dfree(m);
	m=loc2mat(loc, 0);
	assert_int_equal(NY(m), 2);
	assert_float_equal(P(m, 5, 1), loc->locy[5], 0);
	dfree(m);
	assert_null(loc2mat(NULL, 0));
	//1-d loc
	loc_t *l1=mk1dloc(-1, 0.25, 9);
	assert_int_equal(l1->nloc, 9);
	assert_float_equal(l1->locx[8], 1, 1e-15);
	assert_float_equal(l1->locy[8], 0, 0);
	real xv[4]={0, 0.2, 0.4, 0.6};
	loc_t *l2=mk1dloc_vec(xv, 4);
	assert_float_equal(l2->dx, 0.2, 1e-15);
	assert_float_equal(l2->locx[3], 0.6, 0);
	assert_null(mk1dloc_vec(NULL, 0));
	locfree(l1); locfree(l2);
	locfree(loc);
}
/** loc_mcc_ptt, loc_calc_ptt, loc_calc_mod, loc_add_ptt, loc_sub_ptt and loc_remove_ptt. */
static void loc_ptt(void **state){
	(void)state;
	loc_t *loc=mkannloc(2, 0.4, 1./16, 0.5);
	const long nloc=loc->nloc;
	dmat *amp=dnew(nloc, 1);
	for(long i=0; i<nloc; i++) P(amp, i)=1+0.5*loc->locx[i]*loc->locx[i];
	dnormalize_sum(amp, 1);
	//mcc=[1 x y]'*diag(amp)*[1 x y]
	dmat *mcc=loc_mcc_ptt(loc, P(amp));
	real ref[3][3]={{0}};
	for(long i=0; i<nloc; i++){
		real v[3]={1, loc->locx[i], loc->locy[i]};
		for(int a=0; a<3; a++) for(int b=0; b<3; b++) ref[a][b]+=P(amp, i)*v[a]*v[b];
	}
	for(int a=0; a<3; a++) for(int b=0; b<3; b++) assert_float_equal(P(mcc, a, b), ref[a][b], 1e-12);
	assert_null(loc_mcc_ptt(NULL, NULL));
	//pure PTT opd: coefficients are recovered, PTT removed residual is zero
	dmat *opd=dnew(nloc, 1);
	const real ptt[3]={0.3, -2, 1.5};
	loc_add_ptt(opd, ptt, loc);
	for(long i=0; i<nloc; i++) assert_float_equal(P(opd, i), 0.3-2*loc->locx[i]+1.5*loc->locy[i], 1e-14);
	real rms[3], coeff[3];
	loc_calc_ptt(rms, coeff, loc, 0, NULL, P(amp), P(opd));
	for(int i=0; i<3; i++) assert_float_equal(coeff[i], ptt[i], 1e-10);
	assert_float_equal(rms[2], 0, 1e-12);//PTT removed
	//PR=TT for pure PTT
	assert_float_equal(rms[0], rms[1], 1e-12);
	real wvar=0, wmean=0;
	for(long i=0; i<nloc; i++){ wmean+=P(amp, i)*P(opd, i); wvar+=P(amp, i)*P(opd, i)*P(opd, i); }
	assert_float_equal(rms[0], wvar-wmean*wmean, 1e-12);
	//with precomputed imcc
	dmat *imcc=dinvspd(mcc);
	real rms2[3], coeff2[3];
	loc_calc_ptt(rms2, coeff2, loc, 1./P(mcc, 0, 0), imcc, P(amp), P(opd));
	for(int i=0; i<3; i++){
		assert_float_equal(coeff2[i], coeff[i], 1e-10);
		assert_float_equal(rms2[i], rms[i], 1e-10);
	}
	//without amplitude: uniform weighting
	dmat *ones=dnew(nloc, 1);
	dset(ones, 1);
	loc_calc_ptt(NULL, coeff2, loc, 0, NULL, NULL, P(opd));
	for(int i=0; i<3; i++) assert_float_equal(coeff2[i], ptt[i], 1e-10);
	//add focus: PTT removed residual is the weighted focus variance
	dmat *opd2=ddup(opd);
	loc_add_focus(opd2, loc, 0.7);
	real piston=0;
	for(long i=0; i<nloc; i++) piston+=P(opd2, i)-P(opd, i);
	assert_float_equal(piston, 0, 1e-10);//focus is piston free
	loc_calc_ptt(rms2, coeff2, loc, 0, NULL, P(amp), P(opd2));
	assert_true(rms2[2]>0);
	//sub_ptt removes the PTT part
	loc_sub_ptt(opd2, coeff2, loc);
	loc_calc_ptt(NULL, coeff, loc, 0, NULL, P(amp), P(opd2));
	for(int i=0; i<3; i++) assert_float_equal(coeff[i], 0, 1e-10);
	//focus with offset center
	dmat *opd3=dnew(nloc, 1);
	loc_add_focus_offset(opd3, loc, 1, 0.2, -0.1);
	real mean3=dsum(opd3)/nloc;
	assert_float_equal(mean3, 0, 1e-12);
	assert_float_equal(P(opd3, 0)-P(opd3, 1), (pow(loc->locx[0]-0.2, 2)+pow(loc->locy[0]+0.1, 2))-(pow(loc->locx[1]-0.2, 2)+pow(loc->locy[1]+0.1, 2)), 1e-12);
	//loc_remove_ptt on several columns and on a covariance matrix
	dmat *opds=dnew(nloc, 2);
	for(long i=0; i<nloc; i++){
		P(opds, i, 0)=P(opd, i)+loc->locx[i]*loc->locy[i];
		P(opds, i, 1)=2*P(opd, i);
	}
	loc_remove_ptt(opds, loc, P(amp), NULL, 0);
	for(int ic=0; ic<2; ic++){
		loc_calc_ptt(NULL, coeff, loc, 0, NULL, P(amp), PCOL(opds, ic));
		for(int i=0; i<3; i++) assert_float_equal(coeff[i], 0, 1e-10);
	}
	//covariance: C=o*o'; removing PTT on both sides gives zero
	loc_t *small=mksqloc(4, 4, 0.5, 0.5, -0.75, -0.75);
	dmat *o=dnew(16, 1);
	loc_add_ptt(o, ptt, small);
	dmat *cov=NULL;
	dmm(&cov, 0, o, o, "nt", 1);
	loc_remove_ptt(cov, small, NULL, NULL, 1);
	assert_float_equal(dmaxabs(cov), 0, 1e-10);
	dfree(cov); dfree(o); locfree(small);
	//loc_calc_mod with orthonormal modes
	dmat *mod=dnew(nloc, 2);
	for(long i=0; i<nloc; i++){ P(mod, i, 0)=1; P(mod, i, 1)=loc->locx[i]; }
	dgramschmidt(mod, P(amp));
	dscale(mod, 1);//<mod|amp|mod>=sum(amp)=1
	real rmsm[2], coeffm[2];
	loc_calc_mod(rmsm, coeffm, mod, P(amp), P(opd));
	real tot=0, c0=0, c1=0;
	for(long i=0; i<nloc; i++){
		tot+=P(amp, i)*P(opd, i)*P(opd, i);
		c0+=P(amp, i)*P(opd, i)*P(mod, i, 0);
		c1+=P(amp, i)*P(opd, i)*P(mod, i, 1);
	}
	assert_float_equal(coeffm[0], c0, 1e-12);
	assert_float_equal(coeffm[1], c1, 1e-12);
	assert_float_equal(rmsm[0], tot-c0*c0, 1e-12);
	assert_float_equal(rmsm[1], tot-c0*c0-c1*c1, 1e-12);
	dfree(mod);
	dfree(opds); dfree(opd3); dfree(opd2); dfree(ones); dfree(imcc);
	dfree(opd); dfree(mcc); dfree(amp);
	locfree(loc);
}
/** pts_mcc_ptt and pts_ztilt recover the tilt of each subaperture, loc_remove_focus_grad. */
static void loc_pts_ztilt(void **state){
	(void)state;
	const int nsa=3, nx=4;
	const real dsa=0.5, dx=dsa/nx;
	pts_t *pts=ptsnew(nsa, dsa, dsa, nx, nx, dx, dx);
	for(int isa=0; isa<nsa; isa++){
		pts->origx[isa]=isa*dsa-0.75;
		pts->origy[isa]=-0.25+0.1*isa;
	}
	const long np=nsa*nx*nx;
	dmat *amp=dnew(np, 1);
	dset(amp, 1);
	P(amp, 0)=0.5;//gray pixel
	dmat *opd=dnew(np, 1);
	//global tilt plus a different piston in each subaperture
	const real tx=0.3, ty=-0.7;
	for(int isa=0; isa<nsa; isa++){
		for(int iy=0; iy<nx; iy++){
			for(int ix=0; ix<nx; ix++){
				real x=pts->origx[isa]+ix*dx, y=pts->origy[isa]+iy*dx;
				P(opd, ix+iy*nx+isa*nx*nx)=isa+tx*x+ty*y;
			}
		}
	}
	dcell *mcc=pts_mcc_ptt(pts, P(amp));
	assert_int_equal(NX(mcc), nsa);
	//check the first subaperture against a direct sum
	real a00=0, a11=0, a12=0;
	for(int iy=0; iy<nx; iy++) for(int ix=0; ix<nx; ix++){
		real x=pts->origx[0]+ix*dx, y=pts->origy[0]+iy*dx;
		real a=P(amp, ix+iy*nx);
		a00+=a; a11+=a*x*x; a12+=a*x*y;
	}
	assert_float_equal(P(P(mcc, 0), 0, 0), a00, 1e-12);
	assert_float_equal(P(P(mcc, 0), 1, 1), a11, 1e-12);
	assert_float_equal(P(P(mcc, 0), 2, 1), a12, 1e-12);
	for(int isa=0; isa<nsa; isa++) dinvspd_inplace(P(mcc, isa));
	dmat *out=NULL;
	pts_ztilt(&out, pts, mcc, P(amp), P(opd));
	assert_int_equal(NX(out), 2*nsa);
	for(int isa=0; isa<nsa; isa++){
		assert_float_equal(P(out, isa), tx, 1e-10);
		assert_float_equal(P(out, isa+nsa), ty, 1e-10);
	}
	pts_ztilt(&out, pts, mcc, P(amp), P(opd));//accumulates
	assert_float_equal(P(out, 0), 2*tx, 1e-10);
	assert_null(pts_mcc_ptt(NULL, NULL));
	dfree(out);
	dcellfree(mcc);
	dfree(opd); dfree(amp);
	ptsfree(pts);
	//focus mode in gradients is proportional to the subaperture location
	loc_t *saloc=mksqloc(3, 3, 1, 1, -1, -1);
	dmat *grad=dnew(18, 1);
	for(int i=0; i<9; i++){
		P(grad, i)=0.4*saloc->locx[i]+(i==4?0.1:0);
		P(grad, i+9)=0.4*saloc->locy[i];
	}
	real focus=loc_remove_focus_grad(grad, saloc, 1);
	assert_float_equal(focus, 0.2, 1e-12);
	for(int i=0; i<9; i++){
		assert_float_equal(P(grad, i), i==4?0.1:0, 1e-12);
		assert_float_equal(P(grad, i+9), 0, 1e-12);
	}
	P(grad, 0)=0.4*saloc->locx[0];
	real f2=loc_remove_focus_grad(grad, saloc, 0);//only estimate
	assert_true(f2>0);
	assert_float_equal(P(grad, 0), 0.4*saloc->locx[0], 0);
	dfree(grad);
	locfree(saloc);
}
/** loc_ellipse_add area, locstretch and loctransform with non-integer powers or a dmat. */
static void loc_ellipse_stretch(void **state){
	(void)state;
	const real dx=1./32;
	loc_t *loc=mksqloc_auto(96, 96, dx, dx);
	dmat *phi=dnew(loc->nloc, 1);
	loc_ellipse_add(phi, loc, 0.1, -0.2, 0.8, 0.4, 2);
	assert_float_equal(dsum(phi)*dx*dx, 2*M_PI*0.8*0.4, 0.002);
	assert_float_equal(dmax(phi), 2, 1e-12);
	//a circle gives the same as loc_circle_add
	dmat *phi2=dnew(loc->nloc, 1);
	dzero(phi);
	loc_ellipse_add(phi, loc, 0, 0, 0.5, 0.5, 1);
	loc_circle_add(phi2, loc, 0, 0, 0.5, 0, 1);
	assert_float_equal(dsum(phi), dsum(phi2), 0.01*dsum(phi2));
	dfree(phi2); dfree(phi);
	//locstretch: rotate, scale x by frac, rotate back
	loc_t *l=mksqloc(3, 2, 0.5, 0.5, -0.5, 0);
	loc_t *ls=locdup(l);
	locstretch(ls, M_PI/2, 2);//stretch along y
	for(long i=0; i<l->nloc; i++){
		assert_float_equal(ls->locx[i], l->locx[i], 1e-12);
		assert_float_equal(ls->locy[i], 2*l->locy[i], 1e-12);
	}
	locfree(ls);
	//non-integer and negative powers use pow()
	dmat *coeff=dnew(2, 4);
	P(coeff, 0, 0)=0.5; P(coeff, 0, 1)=0; P(coeff, 0, 2)=3; P(coeff, 0, 3)=0;//new x=3*x^0.5
	P(coeff, 1, 0)=0; P(coeff, 1, 1)=1; P(coeff, 1, 2)=0; P(coeff, 1, 3)=1;//new y=y
	loc_t *lp=mksqloc(3, 2, 0.5, 0.5, 0.25, 1);
	loc_t *lt=loctransform2(lp, coeff);
	for(long i=0; i<lp->nloc; i++){
		assert_float_equal(lt->locx[i], 3*sqrt(lp->locx[i]), 1e-12);
		assert_float_equal(lt->locy[i], lp->locy[i], 1e-12);
	}
	locfree(lt);
	//dmat passed in place of the string
	lt=loctransform(lp, (const char *)coeff);
	assert_float_equal(lt->locx[2], 3*sqrt(lp->locx[2]), 1e-12);
	locfree(lt);
	dfree(coeff);
	locfree(lp);
	locfree(l);
	locfree(loc);
}
/** loc_reduce with continuity requirements. */
static void loc_reduce_cont(void **state){
	(void)state;
	//two columns of 12 points
	const real a0[12]={1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1};
	const real a1[12]={0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0};
	for(int cont=0; cont<3; cont++){
		loc_t *loc=mksqloc(12, 2, 1, 1, 0, 0);
		dmat *amp=dnew(24, 1);
		for(int i=0; i<12; i++){ P(amp, i)=a0[i]; P(amp, i+12)=a1[i]; }
		int *skip=NULL;
		loc_reduce(loc, amp, 0.5, cont, &skip);
		int nkeep=0;
		for(int i=0; i<24; i++) if(!skip[i]) nkeep++;
		assert_int_equal(loc->nloc, nkeep);
		assert_int_equal(PN(amp), nkeep);
		switch(cont){
		case 0://any point below threshold
			assert_int_equal(nkeep, 8+8);
			break;
		case 1://only at the ends of each column
			assert_int_equal(nkeep, 12+9);
			assert_int_equal(skip[4], 0);
			assert_int_equal(skip[12], 1);
			assert_int_equal(skip[19], 0);
			assert_int_equal(skip[23], 1);
			break;
		case 2://switch after 4 consecutive points
			assert_int_equal(skip[4], 1);
			assert_int_equal(skip[7], 1);
			assert_int_equal(skip[8], 0);
			assert_int_equal(skip[12], 1);
			assert_int_equal(skip[13], 1);
			assert_int_equal(skip[14], 0);
			assert_int_equal(skip[19], 0);//isolated gap is kept
			assert_int_equal(nkeep, 8+10);
			break;
		}
		free(skip);
		dfree(amp);
		locfree(loc);
	}
}

UNITTEST_MAIN(
	cmocka_unit_test(loc_mkw_amp),
	cmocka_unit_test(loc_mkw_circular),
	cmocka_unit_test(loc_area),
	cmocka_unit_test(loc_sqlocrot),
	cmocka_unit_test(loc_reduce_sparse),
	cmocka_unit_test(loc_zernike),
	cmocka_unit_test(loc_embed_locstat),
	cmocka_unit_test(loc_nextpow2),
	cmocka_unit_test(loc_parse_poly),
	cmocka_unit_test(loc_convert_keywords),
	cmocka_unit_test(loc_ptt),
	cmocka_unit_test(loc_pts_ztilt),
	cmocka_unit_test(loc_ellipse_stretch),
	cmocka_unit_test(loc_reduce_cont),
)
