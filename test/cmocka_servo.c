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
  Tests of lib/servo.c (servo gain optimization, transfer functions, time
  domain filtering, SHO) and the remaining parts of lib/kalman.c (SDE fitting
  variants, multi-rate kalman_test2).
*/
#include "unittest.h"

#if CPU_SINGLE
#define TOL 1e-4
#else
#define TOL 1e-9
#endif

/*Open loop transfer function of an integrator with gain g used by
  servo_calc_init for dtrat=1: Hol=g*z^-1*exp(-s*dt*(1+al))/(1-z^-1).*/
static dcomplex hol_int(real f, real dt, real al, real g){
	dcomplex s=COMPLEX(0, TWOPI*f);
	dcomplex zinv=cexp(-s*dt);
	return g*zinv*cexp(-s*dt*(1+al))/(1.-zinv);
}
/*A power law PSD (two columns of frequency and PSD) on a linear grid.*/
static dmat *psd_powerlaw(real fmin, real fmax, long nf, real beta, real alpha){
	dmat *psd=dnew(nf, 2);
	for(long i=0; i<nf; i++){
		real f=fmin+(fmax-fmin)*i/(nf-1);
		P(psd, i, 0)=f;
		P(psd, i, 1)=beta*pow(f, alpha);
	}
	return psd;
}
/*1x1 dmat*/
static dmat *dscalar(real v){
	dmat *A=dnew(1, 1);
	P(A, 0)=v;
	return A;
}

/** servo_optim with a PSD input: the optimal gain decreases and the noise
	residual increases with noise; servo_residual reproduces the residuals; the
	total error is minimal at the returned gain. */
static void servo_optim_psd(void **state){
	(void)state;
	const real dt=1./1000;
	dmat *psd=psd_powerlaw(0.1, 500, 5000, 1e-14, -8./3.);
	dmat *sigma2n=dnew(3, 1);
	P(sigma2n, 0)=0;
	P(sigma2n, 1)=1e-15;
	P(sigma2n, 2)=1e-13;
	dcell *gm=servo_optim(dt, 1, 1, 0, 0, 0, 1, psd, sigma2n);
	assert_int_equal(NX(gm), 3);
	assert_int_equal(NY(gm), 1);
	real glast=INFINITY, rnlast=-1;
	for(int i=0; i<3; i++){
		dmat *gi=P(gm, i);
		assert_int_equal(PN(gi), 3);//[g, res_sig, res_n]
		const real g=P(gi, 0);
		assert_true(g>0&&g<1);
		assert_true(g<=glast*(1+1e-6));
		glast=g;
		assert_true(P(gi, 2)>=rnlast);//noise residual grows with noise
		rnlast=P(gi, 2);
		assert_true(P(gi, 1)>0);
		//consistent with servo_residual()
		real namp=0;
		dmat *gain=dscalar(g);
		real res=servo_residual(&namp, psd, dt, 1, 1, gain, 1);
		assert_float_equal(res, P(gi, 1), 1e-6*res);
		assert_float_equal(namp*P(sigma2n, i), P(gi, 2), 1e-6*P(gi, 2)+1e-30);
		//optimum: a lower gain is worse. A higher gain is either worse or violates the phase margin.
		const real tot=P(gi, 1)+P(gi, 2);
		const real gmax=servo_optim_margin(dt, 1, 1, M_PI/4, 0, 0);
		for(int j=-1; j<=1; j+=2){
			P(gain, 0)=g*(1+0.2*j);
			real res2=servo_residual(&namp, psd, dt, 1, 1, gain, 1);
			if(P(gain, 0)<gmax){
				assert_true(res2+namp*P(sigma2n, i)>=tot*(1-1e-6));
			}
		}
		dfree(gain);
	}
	//without noise, the gain is limited by the phase margin only.
	real gmargin=servo_optim_margin(dt, 1, 1, M_PI/4, 0, 0);
	assert_true(P(P(gm, 0), 0)<=gmargin*1.01);
	assert_true(P(P(gm, 0), 0)>gmargin*0.3);
	//The residual is smaller than the input variance
	assert_true(P(P(gm, 0), 1)<psd_inte2(psd));
	dcellfree(gm);

	//NULL sigma2n gives one result
	gm=servo_optim(dt, 1, 1, 0, 0, 0, 1, psd, NULL);
	assert_int_equal(PN(gm), 1);
	assert_float_equal(P(P(gm, 0), 2), 0, 0);
	dcellfree(gm);
	//A SHO in the loop reduces the gain
	dcell *gm2=servo_optim(dt, 1, 1, 0, 50, 0.5, 1, psd, NULL);
	gm=servo_optim(dt, 1, 1, 0, 0, 0, 1, psd, NULL);
	assert_true(P(P(gm2, 0), 0)<P(P(gm, 0), 0));
	dcellfree(gm2);
	dcellfree(gm);
	//dtrat>1
	gm=servo_optim(dt, 4, 0, 0, 0, 0, 1, psd, NULL);
	gm2=servo_optim(dt, 1, 0, 0, 0, 0, 1, psd, NULL);
	assert_true(P(P(gm, 0), 1)>P(P(gm2, 0), 1));//slower loop has more residual
	dcellfree(gm2);
	dcellfree(gm);
	dfree(sigma2n);
	dfree(psd);
}

/** Type II servo optimization returns [g, a, T, res_sig, res_n], with a lead
	filter (0<a<=1) whose signal residual is reproduced by servo_residual. */
static void servo_optim_type2(void **state){
	(void)state;
	const real dt=1./1000;
	dmat *psd=psd_powerlaw(0.1, 500, 5000, 1e-14, -8./3.);
	dmat *sigma2n=dscalar(1e-16);
	dcell *gm=servo_optim(dt, 1, 0, M_PI/4, 0, 0, 2, psd, sigma2n);
	dmat *g=P(gm, 0);
	assert_int_equal(PN(g), 5);
	assert_true(P(g, 0)>0);
	assert_true(P(g, 1)>0&&P(g, 1)<=1);
	assert_true(P(g, 2)>=0);
	real namp=0;
	dmat *gain=dsub(g, 0, 3, 0, 1);
	real res=servo_residual(&namp, psd, dt, 1, 0, gain, 2);
	assert_float_equal(res, P(g, 3), 1e-6*res);
	/*The lead filter gives exactly the requested phase margin. The loop must
	  be accepted as stable: no instability penalty in res_n and a reasonable gain.*/
	assert_true(P(g, 0)>0.01);
	assert_true(P(g, 4)<psd_inte2(psd));
	assert_true(P(g, 3)<psd_inte2(psd));
	dfree(gain);
	dcellfree(gm);
	dfree(sigma2n);
	dfree(psd);
}

/** servo_cl2ol inverts the rejection transfer function applied analytically. */
static void servo_cl2ol_test(void **state){
	(void)state;
	const real dt=1./800, al=0.5, g=0.4, sigma2n=1e-16;
	const long nf=200;
	dmat *psdcl=dnew(nf, 2);
	dmat *psdol=dnew(nf, 2);
	for(long i=0; i<nf; i++){
		real f=1+i*1.5;
		dcomplex hol=hol_int(f, dt, al, g);
		real hrej2=pow(cabs(1./(1.+hol)), 2);
		P(psdol, i, 0)=f;
		P(psdol, i, 1)=1e-12*pow(f, -2);
		P(psdcl, i, 0)=f;
		//servo_cl2ol implements psdol=psdcl/|Hrej|^2-sigma2n*2*dt
		P(psdcl, i, 1)=(P(psdol, i, 1)+sigma2n*2*dt)*hrej2;
		if(i==nf-1){//less than the noise floor: clipped to zero
			P(psdol, i, 1)=0;
			P(psdcl, i, 1)=sigma2n*dt*hrej2;
		}
	}
	dmat *out=servo_cl2ol(psdcl, dt, 1, al, g, sigma2n);
	assert_int_equal(NX(out), nf);
	assert_int_equal(NY(out), 2);
	for(long i=0; i<nf; i++){
		assert_float_equal(P(out, i, 0), P(psdol, i, 0), 0);
		assert_float_equal(P(out, i, 1), P(psdol, i, 1), 1e-6*P(psdol, i, 1));
	}
	dfree(out);
	//zero gain is a copy
	out=servo_cl2ol(psdcl, dt, 1, al, 0, sigma2n);
	assert_dmat_near(out, psdcl, 0);
	dfree(out);
	dfree(psdcl);
	dfree(psdol);
}

/*run a servo with constant error signal for nstep and return the output history of mode 0.*/
static void run_const(servo_t *st, dmat *merr, int nstep, real *hist){
	dmat *out=NULL;
	for(int i=0; i<nstep; i++){
		servo_filter(st, merr);
		servo_output(st, &out);
		hist[i]=P(out, 0);
	}
	dfree(out);
}

/** Time domain integrator: step response, latency, fractional latency, leaky
	integrator, per mode gain, PID, external gain, servo_add and servo_reset. */
static void servo_time_domain(void **state){
	(void)state;
	const real dt=1e-3;
	const int nstep=10;
	real hist[10];
	dmat *merr=dnew(2, 1);
	P(merr, 0)=1; P(merr, 1)=-2;
	//pure integrator, no latency: out_k=g*(k+1)
	dmat *ep=dscalar(0.3);
	servo_t *st=servo_new(merr, NULL, 0, dt, ep);
	assert_int_equal(st->initialized, 1);
	run_const(st, merr, nstep, hist);
	for(int i=0; i<nstep; i++){
		assert_float_equal(hist[i], 0.3*(i+1), TOL);
	}
	assert_float_equal(P(P(st->mintd, 0), 1), -2*0.3*nstep, TOL);
	//servo_add adds to the integrator
	dmat *madj=dnew(2, 1); P(madj, 0)=1; P(madj, 1)=1;
	servo_add(st, madj, 0.5);
	assert_float_equal(P(P(st->mintd, 0), 0), 0.3*nstep+0.5, TOL);
	//servo_reset zeros the state
	servo_reset(st);
	assert_float_equal(dsumabs(P(st->mintd, 0)), 0, 0);
	servo_free(st);

	//latency of 1.5 steps: integer part delays, fractional part interpolates.
	st=servo_new(merr, NULL, 1.5, dt, ep);
	assert_int_equal(st->alint, 1);
	assert_float_equal(st->alfrac, 0.5, 0);
	run_const(st, merr, nstep, hist);
	//step k (0 based) integrates k errors; output averages last two values
	for(int i=0; i<nstep; i++){
		real cur=0.3*i;
		real prev=i>0?0.3*(i-1):0;
		assert_float_equal(hist[i], 0.5*cur+0.5*prev, TOL);
	}
	//No error signal yet with latency returns 0.
	servo_free(st);
	st=servo_new(NULL, NULL, 1, dt, ep);
	assert_int_equal(st->initialized, 0);
	assert_int_equal(servo_filter(st, merr), 0);//delayed
	assert_int_equal(servo_filter(st, merr), 1);
	assert_int_equal(st->initialized, 1);
	servo_free(st);

	//leaky integrator: x=0.9*x+0.5*e converges to 5
	st=servo_new_scalar(merr, 0.9, 0, dt, 0.5);
	real hist2[200];
	run_const(st, merr, 200, hist2);
	assert_float_equal(hist2[0], 0.5, TOL);
	assert_float_equal(hist2[1], 0.5*0.9+0.5, TOL);
	assert_float_equal(hist2[199], 5, 1e-7);
	servo_free(st);
	//two term modified integrator x_k=a0*x_{k-1}+a1*x_{k-2}+g*e
	dmat *ap=dnew(2, 1); P(ap, 0)=0.5; P(ap, 1)=0.25;
	st=servo_new(merr, ap, 0, dt, ep);
	run_const(st, merr, 4, hist);
	real x1=0, x2=0;
	for(int i=0; i<4; i++){
		real x=0.5*x1+0.25*x2+0.3;
		x2=x1; x1=x;
		assert_float_equal(hist[i], x, TOL);
	}
	servo_free(st);
	dfree(ap);

	//per mode gain. ep 4x1 is reshaped to 1x4.
	dmat *merr4=dnew(4, 1);
	dmat *ep4=dnew(4, 1);
	for(int i=0; i<4; i++){
		P(merr4, i)=1;
		P(ep4, i)=0.1*(i+1);
	}
	st=servo_new(merr4, NULL, 0, dt, ep4);
	dmat *out=NULL;
	for(int i=0; i<3; i++){
		servo_filter(st, merr4);
	}
	servo_output(st, &out);
	for(int i=0; i<4; i++){
		assert_float_equal(P(out, i), 3*0.1*(i+1), TOL);
	}
	servo_free(st);
	dfree(ep4);

	//PID: ep=[a,b]: each step adds a*e_k+b*(e_k-e_{k-1})
	dmat *ep2=dnew(2, 1); P(ep2, 0)=0.2; P(ep2, 1)=0.05;
	st=servo_new(merr, NULL, 0, dt, ep2);
	run_const(st, merr, nstep, hist);
	for(int i=0; i<nstep; i++){
		assert_float_equal(hist[i], 0.2*(i+1)+0.05, TOL);
	}
	servo_free(st);
	dfree(ep2);

	//ep is NULL: unit gain
	st=servo_new(merr4, NULL, 0, dt, NULL);
	servo_filter(st, merr4);
	servo_filter(st, merr4);
	servo_output(st, &out);
	for(int i=0; i<4; i++){
		assert_float_equal(P(out, i), 2, TOL);
	}
	servo_free(st);

	//SHO output converges to the integrator value
	dfree(out);
	st=servo_new_sho(merr, NULL, 0.5, dt, ep, 2000, 0.7);
	assert_non_null(st->sho);
	for(int i=0; i<300; i++){
		servo_filter(st, i<5?merr:NULL);
		servo_output(st, &out);
	}
	/*With fractional latency, the output keeps averaging the last two
	  integrator values when no new error signal arrives: 0.5*(1.2+1.5)*/
	assert_float_equal(P(out, 0), 1.35, 1e-3);
	assert_float_equal(P(out, 1), -2.7, 2e-3);
	servo_free(st);
	//no SHO for infinite frequency
	st=servo_new_sho(merr, NULL, 0, dt, ep, INFINITY, 0.7);
	assert_null(st->sho);
	servo_free(st);
	servo_free(NULL);

	dfree(out);
	dfree(merr4);
	dfree(madj);
	dfree(ep);
	dfree(merr);
}

/** Type II time domain filter follows the lead filter recursion, for dmat and dcell inputs. */
static void servo_typeII(void **state){
	(void)state;
	const real dt=1e-3, g=0.4, a=0.3, T=0.002;
	dmat *ep=dnew(3, 1);
	P(ep, 0)=g; P(ep, 1)=a; P(ep, 2)=T;
	dcell *merr=dcellnew_same(2, 1, 3, 1);
	servo_t *st=servo_new(merr, NULL, 0, dt, ep);
	//converted gains
	const real e1a=exp(-dt/(a*T)), e1=exp(-dt/T);
	assert_float_equal(P(st->ep, 1), e1a, TOL);
	assert_float_equal(P(st->ep, 2), e1, TOL);
	real lead=0, pre=0, mint=0, elast=0;
	dcell *out=NULL;
	rand_t rstat;
	seed_rand(&rstat, 5);
	for(int istep=0; istep<20; istep++){
		real e=randu(&rstat)-0.5;
		for(int ic=0; ic<2; ic++){
			for(int i=0; i<3; i++){
				P(P(merr, ic), i)=e*(ic+1)*(i+1);
			}
		}
		servo_filter(st, merr);
		lead=e1a*lead+g*(1-e1a)/(1-e1)*(e-e1*elast);
		elast=e;
		pre+=lead;
		mint+=pre;
		servo_output(st, &out);
		for(int ic=0; ic<2; ic++){
			for(int i=0; i<3; i++){
				assert_float_equal(P(P(out, ic), i), mint*(ic+1)*(i+1), 1e-8*(1+fabs(mint)));
			}
		}
	}
	servo_free(st);
	//Per mode type II gain: 3xnmod
	dmat *ep3=dnew(3, 3);
	for(int i=0; i<3; i++){
		P(ep3, 0, i)=g*(i+1); P(ep3, 1, i)=a; P(ep3, 2, i)=T;
	}
	dmat *merr1=dnew(3, 1);
	dset(merr1, 1);
	st=servo_new(merr1, NULL, 0, dt, ep3);
	servo_filter(st, merr1);
	dmat *out1=NULL;
	servo_output(st, &out1);
	real lead1=g*(1-e1a)/(1-e1);
	for(int i=0; i<3; i++){
		assert_float_equal(P(out1, i), lead1*(i+1), TOL);
	}
	//servo_reset clears all states
	servo_reset(st);
	servo_output(st, &out1);
	assert_float_equal(dsumabs(out1), 0, 0);
	servo_free(st);
	dfree(out1);
	dfree(merr1);
	dfree(ep3);
	dcellfree(out);
	dcellfree(merr);
	dfree(ep);
}

/** servo_test (closed loop with two step delay) rejects a constant, attenuates
	a sinusoid as |1/(1+Hol)| and propagates noise as predicted by servo_residual. */
static void servo_test_closed_loop(void **state){
	(void)state;
	const real dt=1e-3, g=0.4;
	const int nstep=4000;
	dmat *gain=dscalar(g);
	dmat *sigma2n=dnew(1, 1);//zero: no noise
	//constant input is rejected by the integrator
	dmat *in=dnew(nstep, 1);
	dset(in, 1);
	dmat *res=servo_test(in, dt, 1, sigma2n, gain);
	assert_int_equal(NX(res), 1);
	assert_int_equal(NY(res), nstep);
	assert_float_equal(P(res, 0), 1, 0);
	assert_float_equal(P(res, nstep-1), 0, 1e-10);
	dfree(res);
	//sinusoid
	const real f=20;
	for(int i=0; i<nstep; i++){
		P(in, 0, i)=sin(TWOPI*f*i*dt);
	}
	res=servo_test(in, dt, 1, sigma2n, gain);
	real s1=0, s2=0;
	for(int i=nstep/2; i<nstep; i++){
		s1+=pow(P(res, i), 2);
		s2+=pow(P(in, 0, i), 2);
	}
	real ratio=sqrt(s1/s2);
	real hrej=cabs(1./(1.+hol_int(f, dt, 0, g)));
	assert_float_equal(ratio, hrej, 0.01*hrej);
	dfree(res);
	//noise propagation compared with servo_residual
	dzero(in);
	P(sigma2n, 0)=1e-2;
	res=servo_test(in, dt, 1, sigma2n, gain);
	real var=dsumsq(res)/PN(res);
	dmat *psd=psd_powerlaw(0.1, 500, 1000, 1e-14, -2);
	real namp=0;
	servo_residual(&namp, psd, dt, 1, 0, gain, 1);
	assert_true(namp<1);
	assert_float_equal(var, namp*P(sigma2n, 0), 0.1*namp*P(sigma2n, 0));
	dfree(res);
	//multiple modes with 2x2 noise covariance and dtrat=2
	dmat *in2=dnew(2, 2000);
	for(int i=0; i<2000; i++){
		P(in2, 0, i)=1; P(in2, 1, i)=-1;
	}
	dmat *s2n=dnew(2, 2);
	daddI(s2n, 1e-8);
	res=servo_test(in2, dt, 2, s2n, gain);
	assert_int_equal(NX(res), 2);
	assert_float_equal(P(res, 0, 1999), 0, 1e-3);
	assert_float_equal(P(res, 1, 1999), 0, 1e-3);
	assert_float_equal(P(res, 0, 0), 1, 0);
	dfree(res);
	dfree(s2n);
	dfree(in2);
	dfree(psd);
	dfree(in);
	dfree(sigma2n);
	dfree(gain);
}

/** SHO: unit DC gain, overshoot of an underdamped step response, NULL SHO copies input. */
static void servo_sho(void **state){
	(void)state;
	assert_null(sho_new(0, 0.5));
	assert_null(sho_new(INFINITY, 0.5));
	const real f0=50, zeta=0.3, dt=1e-4;
	const int n=2000;
	dmat *x=dnew(1, n);//row vector
	dset(x, 1);
	dmat *y=sho_filter(x, dt, f0, zeta, 0);
	assert_int_equal(NX(y), n);
	real ymax=dmax(y);
	real overshoot=exp(-M_PI*zeta/sqrt(1-zeta*zeta));
	//forward Euler integration overestimates the overshoot slightly.
	assert_float_equal(ymax-1, overshoot, 0.1*overshoot);
	//peak time is half of the damped period
	long imax=0;
	for(long i=0; i<n; i++){
		if(P(y, i)==ymax){ imax=i; break; }
	}
	real tpeak=0.5/(f0*sqrt(1-zeta*zeta));
	assert_float_equal(imax*dt, tpeak, 0.05*tpeak);
	assert_float_equal(P(y, n-1), 1, 1e-3);
	//averaged output with larger steps also settles to 1
	dmat *x2=dnew(400, 2);
	dset(x2, 2);
	dmat *y2=sho_filter(x2, 1e-3, f0, 0.7, 1);
	assert_float_equal(P(y2, 399, 0), 2, 1e-3);
	assert_float_equal(P(y2, 399, 1), 2, 1e-3);
	//With average, the first output is smaller than the instantaneous one
	dmat *y3=sho_filter(x2, 1e-3, f0, 0.7, 0);
	assert_true(P(y2, 0)<P(y3, 0));
	//NULL sho copies
	dcell *xi=dcellnew_same(1, 1, 3, 1);
	dset(P(xi, 0), 1.5);
	dcell *xo=NULL;
	sho_step(&xo, NULL, xi, dt, 1);
	assert_float_equal(P(P(xo, 0), 2), 1.5, 0);
	//sho_reset
	sho_t *sho=sho_new(f0, zeta);
	sho_step(&xo, sho, xi, 1e-3, 0);
	assert_true(fabs(P(P(xo, 0), 0))>0);
	sho_reset(sho);
	assert_float_equal(dsumabs(dmat_cast(P(sho->x, 0))), 0, 0);
	sho_free(sho);
	sho_free(NULL);
	dcellfree(xi);
	dcellfree(xo);
	dfree(y3);
	dfree(y2);
	dfree(x2);
	dfree(y);
	dfree(x);
}

/*SDE coefficient of a damped oscillator.*/
static dmat *sde_coeff(real f0, real zeta, real sigma){
	dmat *coeff=dnew(3, 1);
	real omega0=TWOPI*f0;
	P(coeff, 0)=2*zeta*omega0;
	P(coeff, 1)=omega0*omega0;
	P(coeff, 2)=sigma;
	return coeff;
}
static dmat *sde_psd_table(const dmat *coeff, real df, real fmax){
	long nf=(long)round(fmax/df);
	dmat *f=dnew(nf, 1);
	for(long i=0; i<nf; i++){
		P(f, i)=df*(i+1);
	}
	dmat *psd=sde_psd2(f, coeff);
	dmat *out=dnew(nf, 2);
	for(long i=0; i<nf; i++){
		P(out, i, 0)=P(f, i);
		P(out, i, 1)=P(psd, i);
	}
	dfree(f);
	dfree(psd);
	return out;
}

/** sde_psd: second order specialization agrees with the general order path,
	per mode output, and 3rd order model. */
static void kalman_sde_psd(void **state){
	(void)state;
	dmat *f=dlinspace(0.5, 0.5, 100);
	dmat *coeff=dnew(3, 2);
	P(coeff, 0, 0)=10; P(coeff, 1, 0)=4000; P(coeff, 2, 0)=1;
	P(coeff, 0, 1)=20; P(coeff, 1, 1)=9000; P(coeff, 2, 1)=2;
	dmat *psd=sde_psd2(f, coeff);
	dmat *psd2=dnew(100, 2);//separate modes
	sde_psd(&psd2, f, P(coeff), 3, 2);
	for(long i=0; i<100; i++){
		assert_float_equal(P(psd, i), P(psd2, i, 0)+P(psd2, i, 1), 1e-12*P(psd, i));
		real w2=pow(TWOPI*P(f, i), 2);
		real ref=1./(pow(w2-4000, 2)+100*w2);
		assert_float_equal(P(psd2, i, 0), ref, 1e-12*ref);
	}
	//order 3 path: coeff [c1 c2 c3 sigma] => sigma^2/|s^3+c1 s^2+c2 s+c3|^2
	dmat *c3=dnew(4, 1);
	P(c3, 0)=30; P(c3, 1)=500; P(c3, 2)=2000; P(c3, 3)=3;
	dmat *psd3=sde_psd2(f, c3);
	for(long i=0; i<100; i++){
		dcomplex s=COMPLEX(0, TWOPI*P(f, i));
		real ref=9./pow(cabs(s*s*s+30.*s*s+500.*s+2000.), 2);
		assert_float_equal(P(psd3, i), ref, 1e-10*ref);
	}
	//order 3 with separate columns
	dmat *psd4=dnew(100, 1);
	sde_psd(&psd4, f, P(c3), 4, 1);
	assert_dmat_near(psd4, psd3, 1e-12*dmax(psd3));
	dfree(psd4);
	dfree(psd3);
	dfree(c3);
	dfree(psd2);
	dfree(psd);
	dfree(coeff);
	dfree(f);
}

/** sde_fit with covariance fitting and sde_fit_auto recover the oscillator parameters. */
static void kalman_sde_fit_cov(void **state){
	(void)state;
	const real f0=5, zeta=0.5;
	dmat *coeff0=sde_coeff(f0, zeta, 1);
	dmat *psd=sde_psd_table(coeff0, 0.05, 20);
	//covariance fitting
	dmat *coeff=sde_coeff(f0*1.3, zeta*0.8, 1);
	real diff=sde_fit(&coeff, psd, 0.5, 0);
	assert_true(diff<1e-2);
	real f1=sqrt(P(coeff, 1))/TWOPI;
	real z1=P(coeff, 0)/(2*sqrt(P(coeff, 1)));
	assert_float_equal(f1, f0, 0.1*f0);
	assert_float_equal(z1, zeta, 0.2*zeta);
	dfree(coeff);
	//auto fitting from scratch on a single PSD
	coeff=NULL;
	diff=sde_fit_auto(&coeff, psd, 0);
	assert_non_null(coeff);
	f1=sqrt(P(coeff, 1))/TWOPI;
	assert_float_equal(f1, f0, 0.2*f0);
	dfree(coeff);
	//multiple PSDs as columns
	dmat *psdm=dnew(NX(psd), 3);
	for(long i=0; i<NX(psd); i++){
		P(psdm, i, 0)=P(psd, i, 0);
		P(psdm, i, 1)=P(psd, i, 1);
		P(psdm, i, 2)=P(psd, i, 1)*4;
	}
	diff=sde_fit_auto(&coeff, psdm, 0);
	assert_int_equal(NX(coeff), 3);
	assert_int_equal(NY(coeff), 2);
	//second PSD is 4 times stronger: sigma doubles
	assert_float_equal(P(coeff, 2, 1), 2*P(coeff, 2, 0), 0.05*P(coeff, 2, 1));
	assert_float_equal(P(coeff, 1, 1), P(coeff, 1, 0), 1e-3*P(coeff, 1, 0));
	dfree(coeff);
	//cell of PSDs
	dcell *psdc=dcellnew(1, 1);
	P(psdc, 0)=dref(psd);
	diff=sde_fit_auto(&coeff, psdc, 0);
	assert_int_equal(NY(coeff), 1);
	assert_float_equal(sqrt(P(coeff, 1))/TWOPI, f0, 0.2*f0);
	(void)diff;
	dfree(coeff);
	dcellfree(psdc);
	dfree(psdm);
	dfree(psd);
	dfree(coeff0);
}

/** sde_fit with vibration identification fits the peak and the background separately. */
static void kalman_sde_fit_vib(void **state){
	(void)state;
	//A narrow vibration at 30 Hz on a flat floor.
	dmat *cv=sde_coeff(30, 0.02, 2e-4);
	dmat *psd=sde_psd_table(cv, 0.01, 60);
	for(long i=0; i<NX(psd); i++){
		P(psd, i, 1)+=1e-20;
	}
	//psd_vibid finds the single peak
	dmat *vibs=psd_vibid(psd);
	assert_int_equal(NY(vibs), 1);
	assert_float_equal(P(vibs, 0, 0), 30, 1);
	dfree(vibs);
	dmat *coeff=NULL;
	real diff=sde_fit(&coeff, psd, 0, 1);
	assert_true(isfinite(diff));
	assert_int_equal(NX(coeff), 3);
	assert_true(NY(coeff)>=2);
	//one column is the vibration, the last is the background
	//The vibration is fitted accurately.
	int found=0;
	for(int iy=0; iy<NY(coeff)-1; iy++){
		real fv=sqrt(P(coeff, 1, iy))/TWOPI;
		if(fabs(fv-30)<0.5){
			found++;
			real zeta=P(coeff, 0, iy)/(2*sqrt(P(coeff, 1, iy)));
			assert_float_equal(zeta, 0.02, 0.001);
			assert_float_equal(P(coeff, 2, iy), 2e-4, 2e-5);
		}
	}
	assert_int_equal(found, 1);
	for(long i=0; i<PN(coeff); i++){
		assert_true(isfinite(P(coeff, i)));
	}
	dfree(coeff);
	dfree(psd);
	dfree(cv);
}

/** kalman_test2 with two WFS at different rates, mdirect, projection and goff.
	All three servo flags reduce a slowly varying input. */
static void kalman_multirate(void **state){
	(void)state;
	const real dthi=1./800;
	dmat *coeff=dnew(3, 2);
	for(int im=0; im<2; im++){
		dmat *c=sde_coeff(1+im, 0.7, 1);
		memcpy(PCOL(coeff, im), P(c), 3*sizeof(real));
		dfree(c);
	}
	lmat *dtrat=lnew(2, 1);
	P(dtrat, 0)=1; P(dtrat, 1)=4;
	lmat *mdirect=lnew(2, 1);
	P(mdirect, 1)=1;
	dcell *Gwfs=dcellnew(2, 1);
	P(Gwfs, 0)=dnew(2, 2); daddI(P(Gwfs, 0), 1);
	P(Gwfs, 1)=dnew(2, 2); daddI(P(Gwfs, 1), 1);
	dcell *Cnn=dcellnew(2, 2);
	for(int i=0; i<4; i++){
		P(Cnn, i)=dnew(2, 2);
		daddI(P(Cnn, i), 1e-4);
	}
	dmat *Proj=dnew(2, 2);
	daddI(Proj, 1);
	dcell *goff=dcellnew(2, 1);
	P(goff, 0)=dnew(2, 1);
	P(goff, 1)=dnew(2, 1);
	//input: smoothly varying
	const int nstep=1000;
	dmat *input=dnew(2, nstep);
	for(int i=0; i<nstep; i++){
		P(input, 0, i)=sin(TWOPI*1*i*dthi);
		P(input, 1, i)=cos(TWOPI*0.5*i*dthi);
	}
	real var_in=dsumsq(input);
	for(int flag=0; flag<3; flag++){
		dmat *res=kalman_test2(coeff, dthi, dtrat, mdirect, Gwfs, Cnn, Proj, goff, input, flag);
		assert_int_equal(NX(res), 2);
		assert_int_equal(NY(res), nstep);
		real var_res=0;
		for(int i=nstep/2; i<nstep; i++){
			var_res+=pow(P(res, 0, i), 2)+pow(P(res, 1, i), 2);
		}
		assert_true(var_res<0.1*var_in);
		dfree(res);
	}
	//kalman_update with a dcell that is not contiguous (meas->m not set)
	kalman_t *k=sde_kalman(coeff, dthi, dtrat, mdirect, Gwfs, Cnn, NULL);
	assert_non_null(k);
	kalman_init(k);
	dcell *meas=dcellnew(2, 1);
	P(meas, 0)=dnew(2, 1); P(meas, 1)=dnew(2, 1);
	P(P(meas, 0), 0)=1;
	kalman_update(k, meas, 0);
	assert_true(dsumabs(P(k->xhat, 0))>0);
	dmat *out=NULL;
	kalman_output(k, &out, 0, 1, 0);//out is allocated
	assert_int_equal(NX(out), 2);
	assert_true(dsumabs(out)>0);
	dfree(out);
	dcellfree(meas);
	kalman_free(k);

	dfree(input);
	dcellfree(goff);
	dfree(Proj);
	dcellfree(Cnn);
	dcellfree(Gwfs);
	lfree(mdirect);
	lfree(dtrat);
	dfree(coeff);
}

UNITTEST_MAIN(
	cmocka_unit_test(servo_optim_psd),
	cmocka_unit_test(servo_optim_type2),
	cmocka_unit_test(servo_cl2ol_test),
	cmocka_unit_test(servo_time_domain),
	cmocka_unit_test(servo_typeII),
	cmocka_unit_test(servo_test_closed_loop),
	cmocka_unit_test(servo_sho),
	cmocka_unit_test(kalman_sde_psd),
	cmocka_unit_test(kalman_sde_fit_cov),
	cmocka_unit_test(kalman_sde_fit_vib),
	cmocka_unit_test(kalman_multirate)
)
