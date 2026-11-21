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
  Tests of lib/libmisc.c (noise, CoG, polynomial fitting, wrapping and
  interpolation helpers, distortion calibration, WFE vs field), lib/hyst.c
  (DM hysteresis), lib/comb.c (star combinations), lib/cn2est.c (SLODAR Cn2
  estimation from synthesized gradients) and the remaining parts of
  lib/turbulence.c.
*/
#include "unittest.h"

/*mean and variance of an array*/
static void meanvar(const dmat *A, real *mean, real *var){
	real s=0, s2=0;
	for(long i=0; i<PN(A); i++){
		s+=P(A, i);
		s2+=P(A, i)*P(A, i);
	}
	*mean=s/PN(A);
	*var=s2/PN(A)-(*mean)*(*mean);
}

/** addnoise: mean and variance of photon, read out, background and excess noise, and QE flat fielding. */
static void misc_addnoise(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 1);
	const long n=40000;
	dmat *A=dnew(n, 1);
	real mean, var;
	//photon + background (calibrated) + read out noise
	dset(A, 100);
	addnoise(A, &rstat, 10, 10, NULL, NULL, NULL, 3, 1);
	meanvar(A, &mean, &var);
	assert_float_equal(mean, 100, 0.5);
	assert_float_equal(var, 110+9, 0.05*119);
	//excess noise factor doubles the photon noise variance. Values below 1 are treated as 1.
	dset(A, 100);
	addnoise(A, &rstat, 0, 0, NULL, NULL, NULL, 0, 2);
	meanvar(A, &mean, &var);
	assert_float_equal(mean, 100, 0.5);
	assert_float_equal(var, 200, 0.05*200);
	dset(A, 100);
	addnoise(A, &rstat, 0, 0, NULL, NULL, NULL, 0, 0.5);
	meanvar(A, &mean, &var);
	assert_float_equal(var, 100, 0.05*100);
	/*pixel dependent background and its calibration. Unit QE.*/
	dmat *bk2=dnew(n, 1);
	dmat *bk2c=dnew(n, 1);
	dmat *qe=dnew(n, 1);
	dset(bk2, 20);
	dset(bk2c, 15);
	dset(qe, 1);
	dset(A, 100);
	addnoise(A, &rstat, 0, 0, bk2, bk2c, qe, 1, 1);
	meanvar(A, &mean, &var);
	assert_float_equal(mean, 105, 0.5);
	assert_float_equal(var, 120+1, 0.05*121);
	/*QE<1 reduces the detected photons and is flat fielded: the mean is
	  preserved and the variance becomes (signal*qe+rne^2)/qe^2.*/
	dset(qe, 0.5);
	dset(A, 100);
	addnoise(A, &rstat, 0, 0, bk2, bk2c, qe, 1, 1);
	meanvar(A, &mean, &var);
	assert_float_equal(mean, 105, 0.5);
	assert_float_equal(var, (120*0.5+1)/0.25, 0.05*244);
	//noise free with zero signal
	dzero(A);
	addnoise(A, &rstat, 0, 0, NULL, NULL, NULL, 0, 1);
	assert_float_equal(dsumabs(A), 0, 0);
	dfree(qe);
	dfree(bk2c);
	dfree(bk2);
	dfree(A);
}

/** addnoise_grad: the noise covariance is L*L' with L the lower triangular factor in neal. */
static void misc_addnoise_grad(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 2);
	const long nsa=40000;
	dmat *neal=dnew(nsa, 3);
	const real a=2, b=1, c=0.5;
	for(long i=0; i<nsa; i++){
		P(neal, i, 0)=a;
		P(neal, i, 1)=b;
		P(neal, i, 2)=c;
	}
	dmat *grad=dnew(nsa*2, 1);
	addnoise_grad(grad, neal, &rstat);
	real sxx=0, syy=0, sxy=0;
	for(long i=0; i<nsa; i++){
		real gx=P(grad, i), gy=P(grad, i+nsa);
		sxx+=gx*gx; syy+=gy*gy; sxy+=gx*gy;
	}
	assert_float_equal(sxx/nsa, a*a, 0.05*a*a);
	assert_float_equal(syy/nsa, b*b+c*c, 0.05*(b*b+c*c));
	assert_float_equal(sxy/nsa, a*c, 0.05*a*a);
	dfree(grad);
	dfree(neal);
}

/*Gaussian spot of sigma at (cx,cy) added to im*/
static void add_spot(dmat *im, real cx, real cy, real sigma, real peak){
	for(long iy=0; iy<NY(im); iy++){
		for(long ix=0; ix<NX(im); ix++){
			real r2=pow(ix-cx, 2)+pow(iy-cy, 2);
			P(im, ix, iy)+=peak*exp(-r2/(2*sigma*sigma));
		}
	}
}

/** cog_multi: CoG of symmetric spots near nominal positions. */
static void misc_cog_multi(void **state){
	(void)state;
	dmat *im=dnew(64, 64);
	const real pos[3][2]={{10, 12}, {30.5, 20}, {45, 50.5}};
	for(int i=0; i<3; i++){
		add_spot(im, pos[i][0], pos[i][1], 1.2, 100);
	}
	dmat *loc=dnew(3, 2);
	for(int i=0; i<3; i++){//nominal positions are off by a couple of pixels
		P(loc, i, 0)=pos[i][0]+2;
		P(loc, i, 1)=pos[i][1]-1;
	}
	dmat *cg=NULL;
	assert_int_equal(cog_multi(&cg, im, loc, 9), 0);
	assert_int_equal(NX(cg), 3);
	assert_int_equal(NY(cg), 3);
	for(int i=0; i<3; i++){
		assert_float_equal(P(cg, i, 0), pos[i][0], 1e-10);
		assert_float_equal(P(cg, i, 1), pos[i][1], 1e-10);
		assert_true(P(cg, i, 2)>100);
	}
	dfree(cg);
	dfree(loc);
	dfree(im);
}

/** cog_nea: Monte Carlo NEA of a photon noise limited Gaussian spot is sigma^2/N. */
static void misc_cog_nea(void **state){
	(void)state;
	const real sigma=1.5, nph=1e4;
	dmat *ints=dnew(20, 20);
	add_spot(ints, 9.5, 9.5, sigma, 1);
	dscale(ints, nph/dsum(ints));
	rand_t rstat;
	real nea[4];
	cog_nea(nea, ints, NULL, 0, 0, 2000, &rstat, 0, 0, NULL, NULL, 0);
	real ref=sigma*sigma/nph;
	assert_float_equal(nea[0], ref, 0.15*ref);
	assert_float_equal(nea[3], ref, 0.15*ref);
	assert_float_equal(nea[1], nea[2], 0);
	assert_true(fabs(nea[1])<0.15*ref);
	//read out noise increases NEA. Result is reproducible as the seed is reset.
	real nea2[4];
	cog_nea(nea2, ints, NULL, 0, 0, 2000, &rstat, 1, 1, NULL, NULL, 3);
	assert_true(nea2[0]>nea[0]*1.2);
	real nea3[4];
	cog_nea(nea3, ints, NULL, 0, 0, 2000, &rstat, 1, 1, NULL, NULL, 3);
	assert_float_equal(nea3[0], nea2[0], 0);
	dfree(ints);
}

/** poly2fit recovers a polynomial coordinate transform exactly; polyfit/polyval round trip. */
static void misc_polyfit(void **state){
	(void)state;
	rand_t rstat;
	seed_rand(&rstat, 3);
	const long np=50;
	dmat *in=dnew(np, 2);
	dmat *out=dnew(np, 2);
	for(long i=0; i<np; i++){
		real x=randu(&rstat)*2-1, y=randu(&rstat)*2-1;
		P(in, i, 0)=x; P(in, i, 1)=y;
		P(out, i, 0)=0.1+1.01*x+0.02*y+0.003*x*x-0.004*x*y;
		P(out, i, 1)=-0.2+0.99*y-0.01*x+0.005*y*y;
	}
	dmat *coeff=poly2fit(in, out, 2);
	assert_int_equal(NX(coeff), 6);
	assert_int_equal(NY(coeff), 4);
	for(long i=0; i<np; i++){
		real fx=0, fy=0;
		for(long ic=0; ic<NX(coeff); ic++){
			real t=pow(P(in, i, 0), P(coeff, ic, 0))*pow(P(in, i, 1), P(coeff, ic, 1));
			fx+=P(coeff, ic, 2)*t;
			fy+=P(coeff, ic, 3)*t;
		}
		assert_float_equal(fx, P(out, i, 0), 1e-10);
		assert_float_equal(fy, P(out, i, 1), 1e-10);
	}
	//constant term
	for(long ic=0; ic<NX(coeff); ic++){
		if(P(coeff, ic, 0)==0&&P(coeff, ic, 1)==0){
			assert_float_equal(P(coeff, ic, 2), 0.1, 1e-10);
			assert_float_equal(P(coeff, ic, 3), -0.2, 1e-10);
		}
	}
	dfree(coeff);
	//negative order limits the power of each: bilinear terms
	for(long i=0; i<np; i++){
		P(out, i, 0)=1+P(in, i, 0)*P(in, i, 1);
		P(out, i, 1)=P(in, i, 1);
	}
	coeff=poly2fit(in, out, -1);
	assert_int_equal(NX(coeff), 4);
	for(long ic=0; ic<4; ic++){
		real xo=P(coeff, ic, 0), yo=P(coeff, ic, 1);
		assert_float_equal(P(coeff, ic, 2), (xo==yo)?1:0, 1e-10);
		assert_float_equal(P(coeff, ic, 3), (xo==0&&yo==1)?1:0, 1e-10);
	}
	dfree(coeff);
	//mismatch
	dmat *bad=dnew(np, 3);
	signal_caught=0;
	assert_null(poly2fit(in, bad, 1));
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	dfree(bad);

	//1-d polynomial
	dmat *x=dlinspace(-2, 0.1, 41);
	dmat *c0=dnew(3, 1);
	P(c0, 0)=2; P(c0, 1)=-3; P(c0, 2)=1;
	dmat *y=polyval(x, c0, 0);
	for(long i=0; i<PN(x); i++){
		real xi=P(x, i);
		assert_float_equal(P(y, i), 2*xi*xi-3*xi+1, 1e-12);
	}
	dmat *ys=polyval(x, c0, 1);//separate terms
	assert_int_equal(NY(ys), 3);
	for(long i=0; i<PN(x); i++){
		assert_float_equal(P(ys, i, 0)+P(ys, i, 1)+P(ys, i, 2), P(y, i), 1e-12);
		assert_float_equal(P(ys, i, 2), 1, 0);
	}
	dmat *c1=polyfit(x, y, 2);
	assert_dmat_near(c1, c0, 1e-9);
	dfree(c1);
	//row vector input
	reshape(x, 1, PN(x));
	c1=polyfit(x, y, 2);
	assert_dmat_near(c1, c0, 1e-9);
	dfree(c1);
	//not vectors
	dmat *m=dnew(2, 2);
	signal_caught=0;
	assert_null(polyfit(m, m, 1));
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	dfree(m);
	dfree(ys);
	dfree(y);
	dfree(c0);
	dfree(x);
	dfree(out);
	dfree(in);
}

/** wrap_seq, wrap2range, pchip_wt and atm_interp. */
static void misc_wrap_interp(void **state){
	(void)state;
	const int seq[]={0, 1, 2, 3, 2, 1, 0, 1, 2, 3};
	for(int i=0; i<10; i++){
		assert_int_equal(wrap_seq(i, 4), seq[i]);
	}
	assert_int_equal(wrap_seq(-1, 4), 1);
	assert_int_equal(wrap_seq(-6, 4), 0);
	assert_int_equal(wrap_seq(5, 1), 0);
	assert_int_equal(wrap_seq(5, 0), 0);
	assert_float_equal(wrap2range(370, -180, 180), 10, 1e-12);
	assert_float_equal(wrap2range(-190, -180, 180), 170, 1e-12);
	assert_float_equal(wrap2range(0.5, 0, TWOPI), 0.5, 1e-12);
	assert_float_equal(wrap2range(7, 0, TWOPI), 7-TWOPI, 1e-12);
	assert_float_equal(wrap2range(3, 2, 2), 2, 0);
	//pchip weights sum to 1, interpolate linear functions, and are 1 at the nodes.
	for(int iu=0; iu<=10; iu++){
		float u=iu*0.1f;
		real s=0, s1=0;
		for(int i=-1; i<=2; i++){
			s+=pchip_wt(i, u);
			s1+=pchip_wt(i, u)*i;
		}
		assert_float_equal(s, 1, 1e-6);
		assert_float_equal(s1, u, 1e-6);
	}
	assert_float_equal(pchip_wt(0, 0), 1, 0);
	assert_float_equal(pchip_wt(1, 0), 0, 0);
	assert_float_equal(pchip_wt(1, 1), 1, 0);
	signal_caught=0;
	pchip_wt(3, 0.5);
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	//atm_interp
	real wt=-1;
	assert_int_equal(atm_interp(&wt, 0, 9, 4, 5, 0), 2);
	assert_float_equal(wt, 1, 0);//no interpolation
	assert_int_equal(atm_interp(&wt, 1, 9, 4, 5, 0), 3);
	assert_float_equal(wt, 0, 0);
	for(int interp=1; interp<=2; interp++){
		real w0, w1;
		int i0=atm_interp(&w0, 0, 5, 4, 5, interp);
		int i1=atm_interp(&w1, 1, 5, 4, 5, interp);
		assert_int_equal(i0, 1);
		assert_int_equal(i1, 2);
		assert_float_equal(w0+w1, 1, 1e-12);
		if(interp==1){
			assert_float_equal(w1, 0.25, 1e-12);
		}else{
			assert_float_equal(w1, pow(sin(0.25*M_PI*0.5), 2), 1e-12);
		}
	}
	real wsum=0;
	for(int ips=0; ips<4; ips++){
		real w;
		int iframe=atm_interp(&w, ips, 6, 4, 10, 3);
		assert_int_equal(iframe, ips);//frame 1 plus ips-1
		wsum+=w;
	}
	assert_float_equal(wsum, 1, 1e-6);
}

/** wfe_fov_prep sorts directions by radius; wfe_fov_fill computes RMS WFE and averages equal radii. */
static void misc_wfe_fov(void **state){
	(void)state;
	const int nevl=5;
	dmat *thetax=dnew(nevl, 1);
	dmat *thetay=dnew(nevl, 1);
	//radius (arcsec): 20, 0, 10, 10, 30. Radii within 1 unit are grouped.
	const real tx[5]={20, 0, 10, 0, 0}, ty[5]={0, 0, 0, -10, 30};
	for(int i=0; i<nevl; i++){
		P(thetax, i)=tx[i]*AS2RAD;
		P(thetay, i)=ty[i]*AS2RAD;
	}
	dmat *one=dnew(1, 1);
	assert_null(wfe_fov_prep(one, one, RAD2AS, 1, 0));
	dfree(one);
	dccell *res=wfe_fov_prep(thetax, thetay, RAD2AS, 2, 1);
	assert_int_equal(NX(res), 4);
	const int order[5]={1, 2, 3, 0, 4};
	for(int i=0; i<nevl; i++){
		assert_int_equal((int)P(P(P(res, 0), 0), i), order[i]);
		real r=RSS(tx[order[i]], ty[order[i]]);
		assert_float_equal(P(P(P(res, 1), 1), i, 0), r, 1e-9);
		assert_float_equal(P(P(P(res, 3), 0, 1), i, 0), r, 1e-9);
	}
	//clep of each direction: 3 x nsim with constant WFE variance
	const int nsim=100;
	dcell *clep=dcellnew(nevl, 1);
	const real var[5]={4, 1, 9, 16, 25};//nm^2*1e-18
	for(int i=0; i<nevl; i++){
		P(clep, i)=dnew(3, nsim);
		for(int is=0; is<nsim; is++){
			P(P(clep, i), 0, is)=var[i]*1e-18;
			P(P(clep, i), 2, is)=var[i]*1e-18/4;
		}
	}
	wfe_fov_fill(res, clep, 0, 0);
	dmat *r0=P(P(res, 1), 0);
	dmat *rho=P(P(res, 2), 0);
	assert_float_equal(P(r0, 0, 1), 1, 1e-9);//direction 1
	//directions 2 and 3 have the same radius: averaged
	real avg=sqrt((9+16)*0.5);
	assert_float_equal(P(r0, 1, 1), avg, 1e-9);
	assert_float_equal(P(r0, 2, 1), avg, 1e-9);
	assert_float_equal(P(rho, 1, 1), avg/2, 1e-9);
	assert_float_equal(P(r0, 3, 1), 2, 1e-9);
	assert_float_equal(P(r0, 4, 1), 5, 1e-9);
	//NULL inputs are ignored. Using fewer steps.
	wfe_fov_fill(NULL, clep, 0, 0);
	wfe_fov_fill(res, clep, 1, 10);
	assert_float_equal(P(P(P(res, 1), 1), 4, 1), 5, 1e-9);
	dcellfree(clep);
	cellfree(res);
	dfree(thetax);
	dfree(thetay);
}

/** loc_calib determines the shift between actuator and subaperture grids from the interaction matrix. */
static void misc_loc_calib(void **state){
	(void)state;
	const int nsa=40;
	const real dsa=0.5;
	loc_t *saloc=mksqloc(nsa, nsa, dsa, dsa, -nsa*dsa/2, -nsa*dsa/2);
	loc_t *aloc=mksqloc(nsa+1, nsa+1, dsa, dsa, -nsa*dsa/2, -nsa*dsa/2);
	const real sx[3]={0, 0.1, 0.1}, sy[3]={0, 0, -0.2};
	for(int ic=0; ic<3; ic++){
		dsp *GA=mkg(aloc, aloc, NULL, saloc, NULL, 0, 1, sx[ic], sy[ic], 1);
		propdata_t pd={0};
		pd.locin=aloc;
		pd.locout=saloc;
		dmat *c=loc_calib(GA, &pd, 1);
		assert_non_null(c);
		assert_int_equal(NX(c), 3);
		for(int i=0; i<NX(c); i++){
			real xo=P(c, i, 0), yo=P(c, i, 1);
			if(xo==0&&yo==0){//shift
				assert_float_equal(P(c, i, 2), sx[ic], 0.01);
				assert_float_equal(P(c, i, 3), sy[ic], 0.01);
			}else{//identity
				assert_float_equal(P(c, i, 2), xo, 1e-6);
				assert_float_equal(P(c, i, 3), yo, 1e-6);
			}
		}
		dfree(c);
		dspfree(GA);
	}
	//Layer above the guide star returns NULL
	dsp *GA=mkg(aloc, aloc, NULL, saloc, NULL, 0, 1, 0, 0, 1);
	aloc->ht=1000;
	propdata_t pd={0};
	pd.locin=aloc;
	pd.locout=saloc;
	pd.hs=500;
	assert_null(loc_calib(GA, &pd, 1));
	dspfree(GA);
	locfree(aloc);
	locfree(saloc);
}

/** Hysteresis: identity without hysteresis, closed loop with the specified width, reset, dcell and NULL model. */
static void misc_hyst(void **state){
	(void)state;
	const real stroke=1e-6;
	const int nt=400;
	dmat *cmd=dnew(2, nt);//sinusoid of amplitude stroke
	for(int it=0; it<nt; it++){
		P(cmd, 0, it)=stroke*sin(TWOPI*it/100.);
		P(cmd, 1, it)=-P(cmd, 0, it);
	}
	dmat *real0=hyst_test(0, 2, stroke, cmd);
	assert_dmat_near(real0, cmd, 1e-20);
	dfree(real0);
	const real hysteresis=0.1;
	dmat *res=hyst_test(hysteresis, -1, stroke, cmd);
	//in the last cycle, compare ascending and descending branches at x=0
	real yasc=P(res, 0, 300), ydesc=P(res, 0, 350);
	assert_float_equal(P(cmd, 0, 300), 0, 1e-20);
	assert_true(ydesc>yasc);
	assert_float_equal((ydesc-yasc)/(2*stroke), hysteresis, 0.2*hysteresis);
	//odd symmetry: the second actuator mirrors the first
	for(int it=0; it<nt; it++){
		assert_float_equal(P(res, 1, it), -P(res, 0, it), 1e-20);
	}
	//the output is periodic after the first cycle
	for(int it=200; it<300; it++){
		assert_float_equal(P(res, 0, it), P(res, 0, it+100), 1e-3*stroke);
	}
	//The tips of the loop are near the stroke
	assert_float_equal(P(res, 0, 325), stroke, 0.2*stroke);
	dfree(res);

	//hyst_dmat in place with state continuing across calls; reset restores
	hyst_t *h=hyst_new(hysteresis, 2, INFINITY, 2);
	dmat *out=dnew(2, nt);
	hyst_dmat(h, out, cmd);
	dmat *out2=dnew(2, nt);
	hyst_reset(h);
	hyst_dmat(h, out2, cmd);
	assert_dmat_near(out, out2, 0);
	//dcell with one model per DM
	hyst_reset(h);
	dcell *cmdc=dcellnew(1, 1);
	P(cmdc, 0)=dref(cmd);
	dcell *outc=dcellnew_same(1, 1, 2, nt);
	hyst_dcell(&h, outc, cmdc);
	assert_dmat_near(P(outc, 0), out, 0);
	//NULL model copies
	hyst_dcell(NULL, outc, cmdc);
	assert_dmat_near(P(outc, 0), cmd, 0);
	hyst_dmat(NULL, out2, cmd);
	assert_dmat_near(out2, cmd, 0);
	hyst_free(h);
	hyst_free(NULL);
	hyst_reset(NULL);
	//invalid hysteresis
	signal_caught=0;
	h=hyst_new(1.5, 2, stroke, 1);
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	hyst_free(h);
	cellfree(outc);
	cellfree(cmdc);
	dfree(out2);
	dfree(out);
	dfree(cmd);
}

/** comb_stars enumerates star combinations for one or two WFS types. */
static void misc_comb(void **state){
	(void)state;
	lmat *nwfs=lnew(1, 1);
	P(nwfs, 0)=2;
	lmat *valid=lnew(4, 1);
	//all combinations of 2 from 4 in lexicographic order
	lmat *res=comb_stars(nwfs, valid, 1);
	assert_int_equal(NX(res), 2);
	assert_int_equal(NY(res), 6);
	const int ref[6][2]={{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
	for(int i=0; i<6; i++){
		assert_int_equal(P(res, 0, i), ref[i][0]);
		assert_int_equal(P(res, 1, i), ref[i][1]);
	}
	lfree(res);
	//flag 0: only the first
	res=comb_stars(nwfs, valid, 0);
	assert_int_equal(NY(res), 1);
	assert_int_equal(P(res, 1, 0), 1);
	lfree(res);
	//invalid star is skipped
	P(valid, 1)=-1;
	res=comb_stars(nwfs, valid, 1);
	assert_int_equal(NY(res), 3);
	for(long i=0; i<PN(res); i++){
		assert_true(P(res, i)!=1);
	}
	lfree(res);
	//not enough valid stars: reduce the number of WFS
	lmat *valid2=lnew(2, 1);
	P(valid2, 1)=-1;
	res=comb_stars(nwfs, valid2, 1);
	assert_int_equal(NX(res), 1);
	assert_int_equal(NY(res), 1);
	assert_int_equal(P(res, 0), 0);
	lfree(res);
	//no valid star
	P(valid2, 0)=-1;
	res=comb_stars(nwfs, valid2, 1);
	assert_null(res);
	lfree(valid2);
	//more WFS than stars
	P(nwfs, 0)=5;
	P(valid, 1)=0;
	res=comb_stars(nwfs, valid, 1);
	assert_int_equal(NX(res), 4);
	assert_int_equal(NY(res), 1);
	lfree(res);
	lfree(valid);
	lfree(nwfs);

	//two WFS types: 1 for the first, 2 for the second, from 4 stars
	nwfs=lnew(2, 1);
	P(nwfs, 0)=1; P(nwfs, 1)=2;
	valid=lnew(4, 2);
	res=comb_stars(nwfs, valid, 1);
	assert_int_equal(NX(res), 3);
	assert_int_equal(NY(res), 12);
	for(int ic=0; ic<12; ic++){
		//all distinct, the last two sorted
		assert_true(P(res, 0, ic)!=P(res, 1, ic));
		assert_true(P(res, 0, ic)!=P(res, 2, ic));
		assert_true(P(res, 1, ic)<P(res, 2, ic));
	}
	lfree(res);
	//flag 2: the first type uses the brightest star of each set
	res=comb_stars(nwfs, valid, 2);
	assert_int_equal(NY(res), 4);
	for(int ic=0; ic<4; ic++){
		assert_true(P(res, 0, ic)<P(res, 1, ic));
		assert_true(P(res, 1, ic)<P(res, 2, ic));
	}
	lfree(res);
	//a star invalid for the first type only
	P(valid, 0, 0)=-1;
	res=comb_stars(nwfs, valid, 1);
	for(int ic=0; ic<NY(res); ic++){
		assert_true(P(res, 0, ic)!=0);
	}
	assert_int_equal(NY(res), 9);
	lfree(res);
	//errors
	lmat *valid1=lnew(4, 1);
	signal_caught=0;
	assert_null(comb_stars(nwfs, valid1, 1));
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	lmat *nwfs3=lnew(3, 1);
	lmat *valid3=lnew(4, 3);
	assert_null(comb_stars(nwfs3, valid3, 1));
	assert_int_equal(signal_caught, SIGUSR2);
	signal_caught=0;
	lfree(valid3);
	lfree(nwfs3);
	lfree(valid1);
	lfree(valid);
	lfree(nwfs);
}

/*Average gradient of each subaperture of size ovs pixels (nsa x nsa) from a
  screen, using the corners, starting at pixel offset (ox, oy).*/
static void grad_screen(real *gx, real *gy, const map_t *s, int nsa, int ovs, int ox, int oy, int add){
	const real dsa=s->dx*ovs;
	const dmat *m=s->dmat;
	for(int iy=0; iy<nsa; iy++){
		for(int ix=0; ix<nsa; ix++){
			int x0=ox+ix*ovs, y0=oy+iy*ovs, x1=x0+ovs, y1=y0+ovs;
			real vx=(P(m, x1, y0)+P(m, x1, y1)-P(m, x0, y0)-P(m, x0, y1))/(2*dsa);
			real vy=(P(m, x0, y1)+P(m, x1, y1)-P(m, x0, y0)-P(m, x1, y0))/(2*dsa);
			if(add){
				gx[ix+iy*nsa]+=vx;
				gy[ix+iy*nsa]+=vy;
			}else{
				gx[ix+iy*nsa]=vx;
				gy[ix+iy*nsa]=vy;
			}
		}
	}
}

/** cn2est recovers r0 and the layer weights of a synthesized two layer
	atmosphere from WFS gradients of two guide stars. */
static void misc_cn2est(void **state){
	(void)state;
	const int nsa=20, ovs=4, nstep=80;
	const real dsa=0.5, dx=dsa/ovs, r0=0.2, wt1=0.3;
	loc_t *saloc=mksqloc(nsa, nsa, dsa, dsa, -nsa*dsa/2, -nsa*dsa/2);
	const real dtheta=dsa/5000.;//one subaperture shift at 5 km
	dmat *wfspair=dnew(2, 1);
	P(wfspair, 0)=0; P(wfspair, 1)=1;
	dmat *theta=dnew(2, 2);
	P(theta, 1, 0)=dtheta;
	dmat *hs=dnew(1, 1);
	P(hs, 0)=INFINITY;
	dmat *htrecon=dnew(3, 1);
	P(htrecon, 1)=5000; P(htrecon, 2)=10000;
	assert_null(cn2est_new(NULL, theta, saloc, NULL, 0, hs, htrecon, 1, 30));
	cn2est_t *c=cn2est_new(wfspair, theta, saloc, NULL, 0, hs, htrecon, 1, 30);
	assert_non_null(c);
	assert_int_equal(c->nwfspair, 1);
	//native slodar heights are multiples of 5 km
	dmat *ht=P(c->ht, 0);
	for(long i=0; i<NX(ht); i++){
		assert_float_equal(P(ht, i), 5000*(c->pair[0].iht0+i), 1e-6);
	}
	//gradients: wfs1 sees the upper layer shifted by one subaperture
	const long ng=nsa*nsa*2;
	dcell *grad=dcellnew_same(2, 1, ng, 1);
	dcell *gradall=dcellnew_same(2, 1, ng, nstep);//all steps
	const real r0a=pow(1-wt1, -3./5.)*r0, r0b=pow(wt1, -3./5.)*r0;
	const long nsa2=nsa*nsa;
	for(int is=0; is<nstep; is++){
		map_t *s0=genatm_simple(r0a, 30, -11./3., dx, 128, is*2+1);
		map_t *s1=genatm_simple(r0b, 30, -11./3., dx, 128, is*2+2);
		real *g0=P(P(grad, 0)), *g1=P(P(grad, 1));
		grad_screen(g0, g0+nsa2, s0, nsa, ovs, 8, 8, 0);
		grad_screen(g1, g1+nsa2, s0, nsa, ovs, 8, 8, 0);
		grad_screen(g0, g0+nsa2, s1, nsa, ovs, 8, 8, 1);
		grad_screen(g1, g1+nsa2, s1, nsa, ovs, 8+ovs, 8, 1);
		cn2est_push(c, grad);
		memcpy(PCOL(P(gradall, 0), is), g0, ng*sizeof(real));
		memcpy(PCOL(P(gradall, 1), is), g1, ng*sizeof(real));
		mapfree(s0);
		mapfree(s1);
	}
	assert_int_equal(c->count, nstep);
	cn2est_est(c, 1);
	assert_float_equal(P(c->r0, 0), r0, 0.1*r0);
	assert_float_equal(c->r0m, P(c->r0, 0), 1e-12);
	dmat *wt=P(c->wt, 0);
	assert_float_equal(dsum(wt), 1, 1e-10);
	//layer at 0, 5 km and above
	long i0=-c->pair[0].iht0;
	assert_float_equal(P(wt, i0), 1-wt1, 0.1);
	assert_float_equal(P(wt, i0+1), wt1, 0.1);
	dmat *wtr=P(c->wtrecon, 0);
	assert_float_equal(dsum(wtr), 1, 1e-10);
	assert_float_equal(P(wtr, 0), 1-wt1, 0.1);
	assert_float_equal(P(wtr, 1), wt1, 0.1);

	//cn2est_all with all steps at once (wfstheta in arcsec) gives the same answer.
	dmat *theta_as=ddup(theta);
	dscale(theta_as, RAD2AS);
	dmat *theta_as0=ddup(theta_as);
	cn2est_t *c2=cn2est_all(wfspair, theta_as, saloc, NULL, 0, hs, htrecon, 1, 30, gradall);
	assert_dmat_near(theta_as, theta_as0, 0);//input is not modified
	assert_int_equal(c2->count, nstep);
	assert_float_equal(P(c2->r0, 0), P(c->r0, 0), 1e-6*P(c->r0, 0));
	assert_dmat_near(P(c2->wt, 0), wt, 1e-6);
	cn2est_free(c2);
	//reset
	cn2est_reset(c);
	assert_int_equal(c->count, 0);
	cn2est_free(c);
	cn2est_free(NULL);

	//keepht=2 estimates directly on htrecon. Self pair only estimates r0.
	dmat *wfspair2=dnew(4, 1);
	P(wfspair2, 0)=0; P(wfspair2, 1)=1; P(wfspair2, 2)=0; P(wfspair2, 3)=0;
	c=cn2est_new(wfspair2, theta, saloc, NULL, 0, hs, htrecon, 2, 30);
	assert_int_equal(c->pair[0].nht, 3);
	assert_int_equal(c->pair[1].nht, 1);
	cn2est_push(c, gradall);
	cn2est_est(c, 0);
	assert_float_equal(P(c->r0, 1), r0, 0.15*r0);
	assert_float_equal(dsum(P(c->wtrecon, 0)), 1, 1e-10);
	cn2est_free(c);

	//saa thresholding removes all subapertures: not possible.
	dmat *saa=dnew(nsa2, 1);
	P(saa, 0)=1;
	assert_null(cn2est_new(wfspair, theta, saloc, saa, 0.5, hs, htrecon, 1, 30));
	dfree(saa);
	dfree(wfspair2);
	dfree(theta_as0);
	dfree(theta_as);
	dcellfree(gradall);
	dcellfree(grad);
	dfree(htrecon);
	dfree(hs);
	dfree(theta);
	dfree(wfspair);
	locfree(saloc);
}

/** Turbulence statistics helpers: calc_aniso, calc_aniso2, calc_greenwood, turbcov, turbpsd/spatial_psd. */
static void turb_calc(void **state){
	(void)state;
	real ht[2]={10000, 0}, wt[2]={1, 0.001};
	//single layer (weights below 0.01 are ignored)
	assert_float_equal(calc_aniso(0.2, 2, ht, wt), 0.3144*0.2/10000, 1e-12);
	real ws[1]={10}, wt1[1]={1};
	assert_float_equal(calc_greenwood(0.2, 1, ws, wt1), 0.426/0.2*10, 1e-10);
	//two layers: greenwood frequency scales with weighted v^(5/3)
	real ws2[2]={10, 20}, wt2[2]={0.5, 0.5};
	real fg=0.426/0.2*pow(0.5*pow(10, 5./3.)+0.5*pow(20, 5./3.), 3./5.);
	assert_float_equal(calc_greenwood(0.2, 2, ws2, wt2), fg, 1e-10);
	//a layer conjugated to one of the two DMs is corrected perfectly
	real ht1[1]={8000};
	assert_true(isinf(calc_aniso2(0.2, 1, ht1, wt1, 0, 8000)));
	//better than a single DM.
	real ht2[1]={5000};
	real t2a=calc_aniso2(0.2, 1, ht2, wt1, 0, 8000);
	assert_true(t2a>calc_aniso(0.2, 1, ht2, wt1));
	//symmetric in the order of the conjugates
	assert_float_equal(calc_aniso2(0.2, 1, ht2, wt1, 8000, 0), t2a, 1e-12*t2a);
	//two layers at the DM conjugates
	real ht3[2]={0, 8000}, wt3[2]={0.6, 0.4};
	assert_true(isinf(calc_aniso2(0.2, 2, ht3, wt3, 0, 8000)));

	//turbcov
	dmat *r=dlinspace(0, 0.5, 20);
	const real r0=0.2, L0=30;
	dmat *cov=turbcov(r, 10, r0, L0);
	real k2=pow(0.5e-6/TWOPI, 2);
	real var=0.0863*pow(L0/r0, 5./3.)*k2;//von Karman phase variance
	assert_float_equal(P(cov, 0), var, 0.02*var);
	for(long i=1; i<PN(r); i++){
		assert_true(P(cov, i)<P(cov, i-1));
	}
	//continuous at r=0
	dmat *rs=dnew(2, 1);
	P(rs, 0)=1e-5;
	P(rs, 1)=1e-3;
	dmat *covs=turbcov(rs, 10, r0, L0);
	assert_float_equal(P(covs, 0), P(cov, 0), 1e-3*P(cov, 0));
	/*structure function D(r)=2(C(0)-C(r)) approaches 6.88(r/r0)^(5/3) for
	  r<<L0. The von Karman correction is about 1-1.485*(r/L0)^(1/3).*/
	real d=2*(P(cov, 0)-P(covs, 1));
	real dref=6.88*pow(P(rs, 1)/r0, 5./3.)*k2*(1-1.485*pow(P(rs, 1)/L0, 1./3.));
	assert_float_equal(d, dref, 0.02*dref);
	dfree(covs);
	dfree(cov);
	//kolmogorov
	cov=turbcov(r, 10, r0, INFINITY);
	real coeff=0.5*6.88*k2*pow(r0, -5./3.);
	assert_float_equal(P(cov, 0), coeff*pow(10, 5./3.), 1e-10*P(cov, 0));
	for(long i=0; i<PN(r); i++){
		assert_float_equal(P(cov, 0)-P(cov, i), coeff*pow(P(r, i), 5./3.), 1e-10*P(cov, 0));
	}
	dfree(cov);
	dfree(rs);
	dfree(r);

	//turbpsd follows the von Karman PSD with zero piston
	const long nx=32;
	const real dx=0.25;
	dmat *psd=turbpsd(nx, nx, dx, r0, L0, -11./3., 1);
	real strength=0.0229*pow(r0, -5./3.)*k2;
	real df=1./(nx*dx);
	assert_float_equal(P(psd, 0), 0, 0);
	real f2=pow(3*df, 2)+pow(2*df, 2);
	real ref=strength*pow(f2+pow(L0, -2), -11./6.)*df*df;
	assert_float_equal(P(psd, 3, 2), ref, 1e-10*ref);
	assert_float_equal(P(psd, nx-3, nx-2), ref, 1e-10*ref);//negative frequency
	//power 0.5 is the square root
	dmat *psd2=turbpsd(nx, nx, dx, r0, L0, 0, 0.5);
	assert_float_equal(P(psd2, 3, 2), sqrt(ref), 1e-10*sqrt(ref));
	dfree(psd2);
	//reuse the output when the size matches.
	dmat *psdold=psd;
	spatial_psd(&psd, nx, nx, dx, 1, 0, 0, 0, -2, 1);
	assert_ptr_equal(psd, psdold);
	assert_float_equal(P(psd, 0), 0, 0);
	assert_float_equal(P(psd, 3, 0), pow(3*df, -2)*df*df, 1e-10);
	//frequency cut offs. The reused output has no stale values outside of the cut offs.
	spatial_psd(&psd, nx, nx, dx, 1, 0, 2*df, 5*df, -2, 1);
	assert_ptr_equal(psd, psdold);
	assert_float_equal(P(psd, 1, 0), 0, 0);//below minfreq
	assert_float_equal(P(psd, 3, 0), pow(3*df, -2)*df*df, 1e-10);
	assert_float_equal(P(psd, 6, 0), 0, 0);//above maxfreq
	//different size reallocates
	spatial_psd(&psd, 8, 4, dx, 1, 0, 0, 0, -2, 1);
	assert_int_equal(NX(psd), 8);
	assert_int_equal(NY(psd), 4);
	dfree(psd);
}

/** genscreen with spatially varying r0, genatm_loc, and genscreen_str with rms, zernike and petal modes. */
static void turb_genscreen(void **state){
	(void)state;
	const long nx=64;
	real wt=1, L0=30;
	rand_t rstat;
	seed_rand(&rstat, 10);
	genatm_t cfg={&rstat, &wt, 0.2, &L0, 1./16, 0, 0, -11./3., nx, nx, 1, 0, 0, 0, 0};
	mapcell *s0=genscreen(&cfg);
	//zero strength of log(r0) variation gives the same screen
	dmat *r0log=dnew(2, 1);
	P(r0log, 0)=-2; P(r0log, 1)=0;
	seed_rand(&rstat, 10);
	cfg.r0logpsds=r0log;
	mapcell *s1=genscreen(&cfg);
	assert_dmat_near(P(s1, 0)->dmat, P(s0, 0)->dmat, 1e-20);
	cellfree(s1);
	//nonzero strength scales the screen by a positive factor
	dresize(r0log, 4, 1);
	P(r0log, 1)=1e-2; P(r0log, 2)=0.1; P(r0log, 3)=2;
	seed_rand(&rstat, 10);
	s1=genscreen(&cfg);
	real maxdiff=0;
	for(long i=0; i<PN(P(s0, 0)); i++){
		real a=P(P(s0, 0), i), b=P(P(s1, 0), i);
		if(fabs(a)>1e-12){
			assert_true(b/a>0);
		}
		maxdiff=MAX(maxdiff, fabs(a-b));
	}
	assert_true(maxdiff>0);
	cellfree(s1);
	cellfree(s0);
	dfree(r0log);

	//genatm_loc
	loc_t *loc=mksqloc_auto(32, 32, 1./16, 1./16);
	dmat *opd=genatm_loc(loc, 0.2, 30, 0, 1);
	assert_int_equal(NX(opd), loc->nloc);
	real rms=sqrt(dsumsq(opd)/PN(opd));
	assert_true(rms>1e-9&&rms<1e-5);
	dmat *opd2=genatm_loc(loc, 0.2, 30, 0, 1);
	assert_dmat_near(opd, opd2, 0);
	dfree(opd2);
	dfree(opd);
	locfree(loc);

	//rms specification: same as specifying the corresponding r0
	mapcell *ss=genscreen_str("rms=100;L0=4;nx=64;dx=0.0625;seed=3");
	assert_int_equal(NX(P(ss, 0)), 64);
	assert_true(dsumsq(P(ss, 0)->dmat)>0);
	cellfree(ss);
	//zernike mode
	ss=genscreen_str("rms=100;mode=4;nx=32;dx=0.0625;ht=500;vx=3");
	assert_float_equal(P(ss, 0)->ht, 500, 0);
	assert_float_equal(P(ss, 0)->vx, 3, 0);
	loc=mksqloc_auto(32, 32, 0.0625, 0.0625);
	dmat *z=zernike(loc, 32*0.0625, 0, 0, -4);
	dscale(z, 100e-9);
	for(long i=0; i<PN(z); i++){
		assert_float_equal(P(P(ss, 0), i), P(z, i), 1e-15);
	}
	dfree(z);
	locfree(loc);
	cellfree(ss);
	//a single petal
	ss=genscreen_str("rms=50;petal=-1;npetal=6;nx=32;dx=0.0625;r0=0.2");
	dmat *mod=dnew(6, 1);
	P(mod, 1)=1;
	dadds(mod, -1./6.);
	dscale(mod, 50e-9);
	dmat *popd=petal_mkopd(32, 32, 16, 16, 6, 0, mod);
	assert_dmat_near(P(ss, 0)->dmat, popd, 1e-15);
	dfree(popd);
	dfree(mod);
	cellfree(ss);
	//random petals: piston free
	ss=genscreen_str("r0=0.2;petal=1;npetal=4;nx=32;dx=0.0625;seed=2;L0=10");
	assert_true(dsumsq(P(ss, 0)->dmat)>0);
	cellfree(ss);
}

UNITTEST_MAIN(
	cmocka_unit_test(misc_addnoise),
	cmocka_unit_test(misc_addnoise_grad),
	cmocka_unit_test(misc_cog_multi),
	cmocka_unit_test(misc_cog_nea),
	cmocka_unit_test(misc_polyfit),
	cmocka_unit_test(misc_wrap_interp),
	cmocka_unit_test(misc_wfe_fov),
	cmocka_unit_test(misc_loc_calib),
	cmocka_unit_test(misc_hyst),
	cmocka_unit_test(misc_comb),
	cmocka_unit_test(misc_cn2est),
	cmocka_unit_test(turb_calc),
	cmocka_unit_test(turb_genscreen)
)
