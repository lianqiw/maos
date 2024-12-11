/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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



#include "servo.h"
#include "psd.h"
typedef struct SERVO_CALC_T{
	cmat* s;
	cmat* Hol;      /**<End-to-End open loop transfer function*/
	cmat* Hint;     /**<Descrete Integrator transfer function*/
	cmat* Hsys;     /**<System transfer function, including Hwfs, Hlag, Hdac, Hmir, Hint. without gain*/
	cmat* Hwfs;
	dmat* nu;       /**<[out] The frequency*/
	dmat* psd;      /**<[out] The PSD defined on nu*/
	real sigma2n;  /**<[out] Noise variance*/
	real pmargin; /*phase margin. default: M_PI/4*/
	real fny;     /**Nyqust frequency for dtrat*/
	real var_sig; /**Integration of PSD*/
	/*Output. */
	real res_sig; /**<[out] Residual signal*/
	real res_n;   /**<[out] Residual noise propagation*/
	real gain_n;  /**<[out] noise amplification*/
	real g;
	real a;
	real T;
	int type;
	int keepnu;
}SERVO_CALC_T;

/*Compute phase between -2pi, and 0*/
static inline real phase(comp val){
	real ang=atan2(cimag(val), creal(val));
	if(ang>0){
		ang=ang-TWOPI;
	}
	return ang;
}
/**
   With given open loop transfer function, compute its cross over frequency and
   phase.

   Written: 2010-06-11

   Tested OK: 2010-06-11
*/
static real phase_at_gain(real* fcross, /**<[out] Cross over frequency*/
	const dmat* nu, /**<[in] frequency grid*/
	const cmat* Hol, /**<[in] open loop transfer function defined on nu*/
	real gain/**<[in] compute phase at this gain*/
){
	/*if(cabs(P(Hol,0))<gain){
	error("Hol is less than 1 from the beginning\n");
	}*/

	int found=0;
	real phi=0;
	for(long i=1; i<NX(nu); i++){
		real val=cabs(P(Hol,i));
		if(val<gain){
			real valpre=cabs(P(Hol,i-1));
			real logvalpre=log10(valpre);/*positive. log10(val) is negative. */
			real rat=(logvalpre-log10(gain))/(logvalpre-log10(val));
			*fcross=pow(10, log10(P(nu,i-1))*(1-rat)+log10(P(nu,i))*rat);
			real phi0=phase(P(Hol,i-1));
			real phi1=phase(P(Hol,i));
			//Unwrap the differeance
			real diff=(phi1-phi0)/(2*M_PI);/*positive. */
			diff=diff-round(diff);
			phi1=phi0+diff*2*M_PI;
			phi=phi0*(1-rat)+phi1*rat;
			//info("valpre=%g, val=%g, phi0=%g, phi1=%g, phi=%g\n", valpre, val, phi0, phi1, phi);
			found=1;
			break;
		}
	}
	if(!found){/*warning("Hol does not decrease to 1.\n");*/
		phi=phase(P(Hol,nu->nx-1));
		*fcross=P(nu,nu->nx-1);
	}
	real nphi=phi/(2*M_PI);
	nphi=nphi-floor(nphi)-1;
	phi=nphi*2*M_PI;/*from -2*M_PI to 0. */
	return phi;/*not phase margin*/
}

/**
 * Determine the gain in dB between Hol and 1 when phase of Hol is angle
 * */
static real gain_at_phase(real* fcross, /**<[out] Cross over frequency*/
	const dmat* nu, /**<[in] frequency grid*/
	const cmat* Hol, /**<[in] open loop transfer function defined on nu*/
	real angle     /**<[in] compute gain at this phase*/
){
	if(angle<-TWOPI||angle>0){
		error("angle=%g is invalid\n", angle);
	}
	long i;
	//Skip leading terms if they are belowangle
	for(i=0; i<NX(nu); i++){
		if(phase(P(Hol,i))>angle){
			break;
		}
	}
	int found=0;
	real gain=0;
	for(; i<NX(nu); i++){
		real phi1=phase(P(Hol,i));
		if(phi1<angle){
			real phi0=phase(P(Hol,i-1));/*how much above -pi*/
			real rat=(phi0-angle)/(phi0-phi1);
			real val=cabs(P(Hol,i));
			real valpre=cabs(P(Hol,i-1));
			*fcross=P(nu,i-1)+(P(nu,i)-P(nu,i-1))*rat;
			gain=20*log10(valpre+(val-valpre)*rat);
			found=1;
			break;
		}
	}
	if(!found){/*Loop is unstable*/
		*fcross=0;
		gain=NAN;
	}
	/*if(0){
		info("gain=%g, angle=%g\n", gain, angle);
		static int saved=0;
		if(!saved){
			saved=1;
			writebin(Hol, "Hol");
			writebin(nu, "nu");
		}
	}*/
	return gain;
}

/**
   Make basic arrays for servo analysis.
   al: additional latency
*/
static void servo_calc_init(SERVO_CALC_T* st, const dmat* psdin, real dt, long dtrat, real al){
	real Ts=dt*dtrat;
	st->fny=0.5/Ts;
	st->type=1;//integrator
	if(st->keepnu && psdin){
		st->nu=drefcols(psdin, 0, 1);
		st->psd=drefcols(psdin, 1, 1);
	}else{
		real nu0=-3;
		if(psdin){
			nu0=log10(P(psdin, 0, 0)==0?P(psdin, 1, 0):P(psdin, 0, 0));
		}
		st->nu=dlogspace(nu0, log10(0.5/dt), 1000);//must interpolate to log space. servo_calc_do depends on it.
		if(psdin){
			if(NY(psdin)!=2){
				error("psdin should have two columns\n");
			}
			st->psd=dinterp1(psdin, 0, st->nu, NAN);
		}
	}
	if(psdin){
		st->var_sig=psd_inte2(psdin);
	}
	comp pi2i=COMPLEX(0, TWOPI);
	if(st->Hsys||st->Hwfs||st->Hint||st->s){
		error("Already initialized\n");
	}
	st->Hsys=cnew(NX(st->nu), 1);
	st->Hwfs=cnew(NX(st->nu), 1);
	st->Hint=cnew(NX(st->nu), 1);
	st->s=cnew(NX(st->nu), 1);
	for(long i=0; i<NX(st->nu); i++){
		comp s=P(st->s, i)=pi2i*P(st->nu, i);
		comp zInv=cexp(-s*Ts);
		comp Hint=P(st->Hint,i)=1./(1-zInv);
		comp Hwfs, Hdac;
		if(dtrat==1){//we have a pure delay
			Hwfs=P(st->Hwfs,i)=zInv;
			Hdac=1;
		} else{
			Hwfs=P(st->Hwfs,i)=(1-zInv)/(Ts*s);
			Hdac=Hwfs;
		}
		comp Hlag=cexp(-s*dt*(1+al));/*lag due to readout/computation*/
		comp Hmir=1;/*DM */
		P(st->Hsys,i)=Hwfs*Hlag*Hdac*Hmir*Hint;
	}
	if(!P(st->nu,0)){
		P(st->Hsys,0)=P(st->Hsys,1);
		P(st->s,0)=P(st->s,1);
	}
}
/**
 * Multiply SHO transfer function to Hsys.
 * */
static void
servo_calc_sho(SERVO_CALC_T *st, real f0, real zeta){
	if(f0<=0||isinf(f0)) return;
	real omega0=f0*2*M_PI;
	real omega02=omega0*omega0;
	real zeta2o=2.*zeta*omega0;
	for(long i=0; i<NX(st->nu); i++){
		comp s=P(st->s, i);
		comp Hfsm=omega02/(s*s+omega02+zeta2o*s);
		P(st->Hsys, i)*=Hfsm;
	}
}
static void servo_calc_free(SERVO_CALC_T* st){
	cfree(st->s);
	cfree(st->Hol);
	cfree(st->Hint);
	cfree(st->Hsys);
	cfree(st->Hwfs);
	dfree(st->nu);
	dfree(st->psd);
}
static int servo_isstable(const dmat* nu, const cmat* Hol, real pmargin_thres){
	real fp, fg;
	real pmargin=phase_at_gain(&fp, nu, Hol, 1)+M_PI;
	real gmargin=0-gain_at_phase(&fg, nu, Hol, -M_PI);
	int isstable=pmargin>pmargin_thres && fp<fg && gmargin>0;
	/*if(!isstable){
	   info("Unstable: phase margin: %4.0f deg at %3.0f Hz. Gain margin: %2.0fdB at %3.0f Hz.\n",
		 pmargin*180/M_PI, fp, gmargin, fg);
		 }*/
	return isstable;
}
/*
static real servo_calc_max_gain(SERVO_CALC_T *st, real pmargin){
	if(!st->Hol){
	st->Hol=cnew(NX(st->nu),1);
	}
	cadd(&st->Hol, 0, st->Hsys, 1);
	real angle=pmargin-M_PI;
	if(st->type==2){
	ccwm(st->Hol, st->Hint);
	angle-=M_PI/2.4;//allow more due to lead filter
	}
	real fcross;
	real dB=gain_at_phase(&fcross, st->nu, st->Hol, angle);
	real gain=!isinf(dB)?pow(10, dB/20):0;
	return gain;
	}*/
/**
   Calculate total error. If g0 is 0, use st->g, otherwith use g0 to figure out g, a T.
*/
static real servo_calc_do(SERVO_CALC_T* st, real g0){
	const dmat* nu=st->nu;
	if(!st->Hol){
		st->Hol=cnew(NX(nu), 1);
	}
	/*Compute Hol with the first integrator and gain.*/
	if(g0>EPS){
		st->g=g0;
	}
	cadd(&st->Hol, 0, st->Hsys, st->g);
	if(st->type==2){//type II controller
		real g2=1;/*additional g to multiply. !=1 if g0 is nonzero and type is 2.*/
		ccwm(st->Hol, st->Hint);/*multiply the second integrator*/
		if(fabs(g0)>EPS){/*figure out a, T from new g0*/
			real margin, fcross;
			margin=phase_at_gain(&fcross, nu, st->Hol, 1)+M_PI;/*see how much phase lead is needed*/
			real phineed=st->pmargin-margin;
			real a, T;
			if(phineed*2.2>M_PI){/*lead filter is not suitable*/
				a=1;
				T=0;
			} else{
				a=(1-sin(phineed))/(1+sin(phineed));
				real f0=fcross*sqrt(a);
				T=1./(2.*M_PI*f0);
			}
			/* Hlead is multiplied by sqrt(a) so it has unit gain at cross over frequency*/
			g2=sqrt(a);
			st->g=g0*g2;
			st->a=a;
			st->T=T;
			//dbg("g0=%g, g2=%g, phineed=%.1f\n", g0, g2, phineed*180/M_PI);
		}
		real a=st->a;
		real T=st->T;
		for(int i=0; i<NX(nu); i++){
			comp Hlead=(1+T*P(st->s,i))/(1+a*T*P(st->s,i))*g2;
			P(st->Hol,i)*=Hlead;
		}
	}
	real res_sig=0;
	real sum_n=0, sum_1=0;
	const dmat* psd=st->psd;
	for(int i=0; i<NX(nu); i++){
		comp Hol=P(st->Hol,i);
		comp Hrej=1./(1.+Hol);
		//The gain never reach below -50dB
		res_sig+=P(psd,i)*creal(Hrej*conj(Hrej)+1e-5)*P(nu,i);
		if(P(st->nu,i)<st->fny){
			//compute noise prop only within nyqust frequency
			comp Hcl=Hol*Hrej;
			comp Hwfs=P(st->Hwfs,i);
			comp Hn=Hcl/Hwfs;
			sum_n+=creal(Hn*conj(Hn))*P(nu,i);
			sum_1+=P(nu,i);
		}
	}
	real dlognu=(log(P(nu,nu->nx-1))-log(P(nu,0)))/(nu->nx-1);
	res_sig*=dlognu;
	st->res_sig=res_sig;
	st->gain_n=sum_n/sum_1;
	if(!servo_isstable(nu, st->Hol, st->pmargin)){
		st->gain_n=100;
		/*put a high penalty to drive down the gain*/
		st->res_n=10*(1+g0)*(st->var_sig+st->sigma2n);
		/*warning("Unstable: g0=%g, res_sig=%g, res_n=%g, tot=%g, gain_n=%g sigma2n=%g\n",
			 g0, res_sig, st->res_n, st->res_n+res_sig, st->gain_n, st->sigma2n);*/
	} else{
		if(st->gain_n>1){
			st->gain_n=pow(st->gain_n, 3);/*a fudge factor to increase the penalty*/
		}
		st->res_n=st->sigma2n*st->gain_n;
		/*info("  Stable: g0=%g, res_sig=%g, res_n=%g, tot=%g, gain_n=%g sigma2n=%g\n",
	  	g0, res_sig, st->res_n, st->res_n+res_sig, st->gain_n, st->sigma2n);*/
	}
	return res_sig+st->res_n;
}

/**
   Optimize the servo gains by balancing errors due to noise and signal propagation.

   Written: 2010-06-11

   Tested OK: 2010-06-11

   2011-01-17: The optimization process is quite slow, but the result is only
   dependent on the sigma2n and fs. psdin does not change during the
   simulation. Built a lookup table in skyc using various sigma2n and interpolate
   to get ress, resn, and gain.

   2012-03-28: Cleaned up this routine. groupped similar calculations together.
   Change g2=a from g=sqrt(a)

   2013-06-27: The maximum gain should not be limited to 0.5 beause it is later scaled by sqrt(a);

   sigma2n is a dmat array of all wanted sigma2n.
   Returns a zfarray of a dmat of [g0, a, T, res_sig, res_n]

   The upper end of frequency must be nyquist freq so that noise transfer can be
   computed. But we need to capture the turbulence PSD beyond nyquist freq,
   which are uncorrectable.

   2022-03-04: added SHO response parameter.

*/
dcell* servo_optim(real dt, long dtrat, real al, real pmargin, real f0, real zeta,
	int servo_type, const dmat *psdin, const dmat *sigma2n){
	if(!pmargin) pmargin=M_PI/4;
	SERVO_CALC_T st={0}; //memset(&st, 0, sizeof(SERVO_CALC_T));
	servo_calc_init(&st, psdin, dt, dtrat, al);
	st.type=servo_type;
	st.pmargin=pmargin;
	if(f0>0){
		servo_calc_sho(&st, f0, zeta);
	}
	int ng=1;
	switch(servo_type){
	case 1:
		ng=1;break;
	case 2:
		ng=3; break;
	default:
		error("Invalid servo_type=%d\n", servo_type);
	}
	if(!psdin) sigma2n=NULL;
	dcell* gm=dcellnew(sigma2n?NX(sigma2n):1, sigma2n?NY(sigma2n):1);
	if(psdin){//minimize total error.
		real g0_step=1e-6;
		real g0_min=1e-6;/*the minimum gain allowed.*/
		real g0_max=2.0;
		for(long ins=0; ins<NX(gm)*NY(gm); ins++){
			st.sigma2n=sigma2n?P(sigma2n,ins):0;
			real g0=golden_section_search((golden_section_fun)servo_calc_do, &st, g0_min, g0_max, g0_step);
			servo_calc_do(&st, g0);
			P(gm,ins)=dnew(ng+2, 1);
			P(P(gm,ins),0)=st.g;
			if(servo_type==2){
				P(P(gm,ins),1)=st.a;
				P(P(gm,ins),2)=st.T;
			}
			P(P(gm,ins),ng)=st.res_sig;
			P(P(gm,ins),ng+1)=st.res_n;
			/*info("g0=%.1g, g2=%.1g, res_sig=%.1g, res_n=%.1g, tot=%.1g, gain_n=%.1g sigma2n=%.1g\n",
			g0, st.g, st.res_sig, st.res_n, st.res_n+st.res_sig, st.gain_n, st.sigma2n);*/
		}/*for ins. */
	}else{//only optimize for phase margin
		real fcross;
		real gdB=gain_at_phase(&fcross, st.nu, st.Hsys, -M_PI+pmargin);
		P(gm,0)=dnew(1,1);
		P(P(gm, 0), 0)=1./pow(10, gdB*0.05);
	}
	servo_calc_free(&st);
	if(isnan(P(P(gm,0),0))){
		error("g=%g\n", P(P(gm, 0), 0));
	}
	return gm;
}
/**
 * Optimize integrator control with optional SHO response (output) to have correct phase margin (pi/4) when abs(Hol)==1.
 * */
real servo_optim_margin(real dt, long dtrat, real al, real pmargin, real f0, real zeta){
	dcell *g=servo_optim(dt, dtrat, al, pmargin, f0, zeta, 1, NULL, 0);
	real gain=P(P(g,0),0);
	dcellfree(g);
	return gain;
	/*SERVO_CALC_T st={0};
	servo_calc_init(&st, NULL, dt, dtrat, al);
	st.type=1;
	st.pmargin=pmargin;
	if(f0>0){
		servo_calc_sho(&st, f0, zeta);
	}
	real fcross;
	
	real gdB=gain_at_phase(&fcross, st.nu, st.Hsys, -M_PI+pmargin);
	servo_calc_free(&st);
	if(0){
		writebin(st.nu, "nu");
		writebin(st.Hsys, "Hsys");
		info("gain=%gdB. fcross=%g\n", gdB, fcross);
	}
	return 1./pow(10, gdB*0.05);*/
}
/**
   Convert Closed loop residual PSD back to OL psd using rejection transfer function:
   PSD_OL=PSD_CL/Hrej-sigma2n/F_nyquist;
   Only implemented for integrator.
 */
dmat* servo_cl2ol(const dmat* psdcl, real dt, long dtrat, real al, real gain, real sigma2n){
	if(fabs(gain)<1e-15){
		return ddup(psdcl);
	}
	SERVO_CALC_T st=(SERVO_CALC_T){.keepnu=1}; //avoid interpolation
	servo_calc_init(&st, psdcl, dt, dtrat, al);
	const dmat* nu=st.nu;
	const dmat* psd=st.psd;
	dmat* psdol=dnew(NX(psd), psd->ny+1);
	real psdn=sigma2n*(dt*dtrat*2);
	cadd(&st.Hol, 0, st.Hsys, gain);
	for(int i=0; i<NX(nu); i++){
		comp Hol=P(st.Hol,i);
		comp Hrej=1./(1.+Hol);
		real normHrej=creal(Hrej*conj(Hrej));
		//comp Hcl=Hol*Hrej;
		//comp Hwfs=P(st.Hwfs,i);
		//comp Hn=Hcl/Hwfs;
		//real normHn=Hn*conj(Hn);
		P(psdol, i, 0)=P(nu,i);//frequency.
		for(int icol=0; icol<NY(psd); icol++){
			P(psdol, i, icol+1)=P(psd, i, icol)/(normHrej)-psdn;
			if(P(psdol, i, icol+1)<0){
				P(psdol, i, icol+1)=0;
			}
		}
	}
	servo_calc_free(&st);
	return psdol;
}
/**
   Compute the residual error after servo rejection of a type II servo given the gain.

   Written: 2010-06-11

   Tested OK: 2010-06-11
*/
real servo_residual(real* noise_amp, const dmat* psdin, real dt, long dtrat, real al, const dmat* gain, int servo_type){
	SERVO_CALC_T st={0}; //memset(&st, 0, sizeof(st));
	servo_calc_init(&st, psdin, dt, dtrat, al);
	st.type=servo_type;
	switch(servo_type){
	case 1:
		st.g=P(gain,0);
		break;
	case 2:
		st.g=P(gain,0);
		st.a=P(gain,1);
		st.T=P(gain,2);
		break;
	default:
		error("Invalid type\n");
	}
	servo_calc_do(&st, 0);
	*noise_amp=st.gain_n;
	servo_calc_free(&st);
	return st.res_sig;
}
/**
   Apply type II servo filter on measurement error and output integrator.  gain
   must be 3x1 or 3x5.  */
static inline void
servo_typeII_filter(const dmat *gain, cell *mpreintc, cell *mleadc, cell *merrlastc, const cell *merrc){
	if(!merrc) return;
	int indmul=0;
	if(NX(gain)!=3){
		error("Wrong format in gain\n");
	}
	real gg, e1a, e1;
	if(iscell(merrc)){
		for(int ic=0; ic<PN(merrc); ic++){
			servo_typeII_filter(gain, P(mpreintc, ic), P(mleadc, ic), P(merrlastc, ic), P(merrc, ic));
		}
	}else{
		dmat* merr=dmat_cast(merrc);
		dmat* mlead=dmat_cast(mleadc);
		dmat* merrlast=dmat_cast(merrlastc);
		int nmod=0;/*number of modes*/
		nmod=PN(merr);
		if(NY(merr)==1){
			nmod=NX(merr);
		} else{
			if(NX(merr)!=1){
				error("Don't handle this case\n");
			}
			nmod=NY(merr);
		}
		if(NY(gain)==1){
			indmul=0;
		} else if(NY(gain)==nmod){
			indmul=1;
		} else{
			error("Wrong format in gain\n");
		}
		/**2013-12-03: fix lead filter implementation. gg is relocated
		   2013-12-05: Implemented a more accurate and robust lead filter
		 */
		for(int imod=0; imod<nmod; imod++){
			int indm=imod*indmul;
			gg=P(gain, 0, indm);
			e1a=P(gain, 1, indm);
			e1=P(gain, 2, indm);
			P(mlead,imod)=e1a*P(mlead,imod)+gg*(1-e1a)/(1-e1)*(P(merr,imod)-e1*P(merrlast,imod));
		}
		dcelladd(&merrlastc, 0, merrc, 1);//save for next time.
		dcelladd(&mpreintc, 1, mleadc, 1);//integrator
	}
}
static void servo_init(servo_t* st, const anyarray merr_){
	cell* merr=merr_.c;
	if(!merr) return;
	if(NX(st->ep)>1){
		dcelladd(&st->mpreint, 0, merr, EPS);
	}
	if(NX(st->ep)==3){
		dcelladd(&st->mlead, 0, merr, EPS);
		dcelladd(&st->merrlast, 0, merr, EPS);
	}
	for(int i=0; i<NX(st->mint); i++){
		dcelladd(&P(st->mint, i), 0, merr, EPS);
	}
	st->initialized=1;
}
/**
   Update servo parameters
*/
static void servo_update_ep(servo_t* st, const dmat* ep){
	dfree(st->ep);
	if(NX(ep)!=3){//type I
		st->ep=ddup(ep);
		if(NY(st->ep)==1&&NX(st->ep)>3){
			reshape(st->ep, 1, NX(st->ep));
		}
	} else{//type II. convert data format
		st->ep=dnew(NX(ep), NY(ep));
		for(int i=0; i<NY(ep); i++){
			P(st->ep,0,i)=P(ep,0,i);
			real a=P(ep, 1, i);
			real T=P(ep, 2, i);
			P(st->ep,1,i)=exp(-st->dt/(a*T));
			P(st->ep,2,i)=exp(-st->dt/T);
		}
	}
}
/**
   Initialize servo_t. al is additional latency
*/
servo_t* servo_new(anyarray merr, const dmat* ap, real al, real dt, const dmat* ep){
	servo_t* st=mycalloc(1, servo_t);
	if(ap){
		st->ap=ddup(ap);
	} else{
		st->ap=dnew(2, 1);
		P(st->ap,0)=1;
	}
	if(NX(st->ap)<2){
		dresize(st->ap, 2, 1);//2 element to ensure we keep integrator history.
	}
	st->mint=cellnew(NX(st->ap), 1);
	st->dt=dt;//this is simulation time step. not loop update time.
	st->alint=(int)floor(al);
	st->alfrac=al-floor(al);
	st->merrhist=cellnew(st->alint+1, 1);
	servo_update_ep(st, ep);
	servo_init(st, merr);
	return st;
}
/**
 * Initialize servo_t with scalar ap and ep
 * */
servo_t *servo_new_scalar(anyarray merr, real ap, real al, real dt, real ep){
	dmat *ap2=dnew(1,1); P(ap2,0)=ap;
	dmat *ep2=dnew(1,1); P(ep2,0)=ep;
	servo_t *st=servo_new(merr, ap2, al, dt, ep2);
	dfree(ap2); dfree(ep2); 
	return st;
}
/**
 * Initialize servo_t with an sho
 * */
servo_t *servo_new_sho(anyarray merr, const dmat *ap, real al, real dt, const dmat *ep, real f0, real zeta){
	servo_t *st=servo_new(merr, ap, al, dt, ep);
	if(f0>0 && !isinf(f0)){
		st->sho=sho_new(f0, zeta);
	}
	return st;
}
/**
   prepare the modified integrator by shifting commands. similar to laos.
   P(inte,0)=P(inte,0)*ap[0]+P(inte,1)*ap[1]+...
*/
static void servo_shift_ap(servo_t* st){
	const dmat* ap=st->ap;
	if(!ap) return; //no need to shift.
	if(NX(st->mint)<NX(ap)){
		cellresize(st->mint, NX(ap), 1);
	}
	if(!st->initialized) return;
	cell *recycle=P(st->mint,ap->nx-1);
	dcellscale(recycle, P(ap, ap->nx-1));
	for(int iap=ap->nx-2; iap>=0; iap--){
		dcelladd(&recycle, 1, P(st->mint,iap), P(ap, iap));
	}
	for(int iap=ap->nx-1; iap>0; iap--){
		P(st->mint, iap)=P(st->mint,iap-1);/*shifting */
	}
	P(st->mint,0)=recycle;/*new command. */
}
/*A FIFO queue to add delay*/
static const cell* servo_shift_al(servo_t* st, const cell* merr){
	if(!st->alint){
		return merr;
	} else{
		long nhist=NX(st->merrhist);
		cell* recycle=P(st->merrhist,0);
		for(int i=0; i<nhist-1; i++){
			P(st->merrhist,i)=P(st->merrhist,i+1);
		}
		P(st->merrhist,nhist-1)=recycle;
		if(!merr){
			cellfree(P(st->merrhist,nhist-1));
		} else{
			dcelladd(&P(st->merrhist,nhist-1), 0, merr, 1);
		}
		return P(st->merrhist,0);
	}
}
/*
	Scale each row by each value of ep
*/
static void scale_row_ep(cell *A, dmat *ep){
	if(!A) return;
	if(iscell(A)){
		for(int ic=0; ic<PN(A); ic++){
			scale_row_ep(P(A, ic), ep);
		}
	}else{
		dmat *Ad=dmat_cast(A);
		if(NX(Ad)!=PN(ep)){
			error("Mismatch: A is %ldx%ld, ep is %ldx%ld\n", NX(A), NY(A), NX(ep), NY(ep));
			return;
		}
		for(int iy=0; iy<NY(Ad); iy++){
			for(int ix=0; ix<NX(Ad); ix++){
				P(Ad, ix, iy)*=P(ep, ix);
			}
		}
	}
}
/**
   Applies type I or type II filter based on number of entries in gain.
*/
int servo_filter(servo_t* st, const anyarray _merr){
	if(!st->initialized&&_merr.c){
		servo_init(st, _merr);
	}
	const cell* merr=servo_shift_al(st, _merr.c);
	if(!merr) return 0;
	servo_shift_ap(st);
	if(!st->mint){
		error("servo_t must be created using servo_new()\n");
	}
	switch(NX(st->ep)){
	case 1://type I
		if(NY(st->ep)==1){//single gain
			dcelladd(&st->mpreint, 0, merr, P(st->ep,0));//just record what is added.
		} else{//gain per mode (row)
			dcelladd(&st->mpreint, 0, merr, 1);
			scale_row_ep(st->mpreint, st->ep);
		}
		break;
	case 2:{//PID controller
		if(NY(st->ep)!=1) error("not supported\n");
		real g1=P(st->ep,0)+P(st->ep,1);
		real g2=-P(st->ep,1);
		dcelladd(&st->mpreint, 0, merr, g1);
		dcelladd(&st->mpreint, 1, st->merrlast, g2);
		dcelladd(&st->merrlast, 0, merr, 1);
	}
		  break;
	case 3://type II
		servo_typeII_filter(st->ep, st->mpreint, st->mlead, st->merrlast, merr);
		break;
	default:
		error("Invalid: NX(st->ep)=%ld\n", NX(st->ep));
	}
	dcelladd(&P(st->mint, 0), 1, st->mpreint, 1);
	return 1;
}
/**
   Adjust integrator content without shift.
*/
void servo_add(servo_t* st, const anyarray madj, real alpha){
	dcelladd(&P(st->mint, 0), 1, madj.c, alpha);
}
/**
   Create servo output. It handles st->alfrac. and outputs the averaged position over the integration period.
 */
void servo_output(const servo_t* st, panyarray out_){
	cell **out=out_.c;
	dcellscale(*out, 0);
	if(st->sho){//filter output using SHO
		if(st->alfrac){//from previous output
			sho_step(&st->sho->ytmp, st->sho, P(st->mint, 1), st->alfrac*st->dt, 1);
			dcelladd(out, 1, st->sho->ytmp, st->alfrac);
		}
		sho_step(&st->sho->ytmp, st->sho, P(st->mint, 0), (1.-st->alfrac)*st->dt, 1);
		dcelladd(out, 1, st->sho->ytmp, 1.-st->alfrac);
	}else{
		if(st->alfrac){//previous output
			dcelladd(out, 1, P(st->mint, 1), st->alfrac);
		}
		dcelladd(out, 1, P(st->mint, 0), 1.-st->alfrac);//current output
	}
}

/**
   test type I/II filter with ideal measurement to make sure it is implemented correctly.
*/
dmat* servo_test(dmat* input, real dt, int dtrat, dmat* sigma2n, dmat* gain){
	if(NY(input)==1){/*single mode. each column is for a mode.*/
		reshape(input, 1, NX(input));
	}
	int nmod=NX(input);
	dmat* pinput=input;
	dmat* merr=dnew(nmod, 1);
	dcell* mreal=dcellnew(1, 1);
	dmat* mres=dnew(nmod, NY(input));
	dmat* sigman=NULL;
	if(dnorm(sigma2n)>0){
		sigman=dchol(sigma2n);
	}
	dcell* meas=dcellnew(1, 1);
	dmat* noise=dnew(nmod, 1);
	servo_t* st2t=servo_new(NULL, NULL, 0, dt*dtrat, gain);
	rand_t rstat;
	seed_rand(&rstat, 1);
	dmat* pmres=mres;
	/*two step delay is ensured with the order of using, copy, acc*/
	for(int istep=0; istep<NY(input); istep++){
		memcpy(P(merr), PCOL(pinput, istep), nmod*sizeof(real));
		dadd(&merr, 1, P(mreal,0), -1);
		memcpy(PCOL(pmres, istep), P(merr), sizeof(real)*nmod);
		if(istep%dtrat==0){
			dzero(P(meas,0));
		}
		dadd(&P(meas,0), 1, merr, 1);/*average the error. */
		dcellcp(&mreal, P(st2t->mintc,0));
		if((istep+1)%dtrat==0){
			if(dtrat!=1) dscale(P(meas,0), 1./dtrat);
			if(sigman){
				drandn(noise, 1, &rstat);
				if(NX(sigman)>0){
					dmm(&P(meas,0), 1, sigman, noise, "nn", 1);
				} else{
					dadd(&P(meas,0), 1, noise, P(sigman,0));
				}
			}
			servo_filter(st2t, meas);
		}
	}
	dfree(sigman);
	dfree(merr);
	dcellfree(mreal);
	dcellfree(meas);
	servo_free(st2t);
	return mres;
}
/**
   reset the data to 0.
*/
void servo_reset(servo_t* st){
	dcellscale(st->mlead, 0);
	dcellscale(st->merrlast, 0);
	dcellscale(st->mpreint, 0);
	dcellscale(st->merrhist, 0);
	dcellscale(st->mint, 0);
}
/**
   Free servo_t struct
*/
void servo_free(servo_t* st){
	if(!st) return;
	cellfree(st->mlead);
	cellfree(st->merrlast);
	cellfree(st->mpreint);
	cellfree(st->merrhist);
	cellfree(st->mint);
	dfree(st->ap);
	dfree(st->ep);
	sho_free(st->sho);
	free(st);
}
/**
   Second harmonic oscillator. Initialization.
   Equation is is d^x/dt^2+c1*dx/dt+c2*x=c2*xi
   State space model is 
   |x'|'   | -c1 -c2 | |x'|      |xi|
   |x |  = |   1   0 | | x| + c2 |0 |
   where x' is dx/dt and written as dx in struct.
   we use c2*xi so that it is in the same unit and scale as x.
 */
sho_t* sho_new(real f0,   /**<Resonance frequency*/
	real zeta  /**<Damping*/
	){
	if(!f0||isinf(f0)) return NULL;
	sho_t* out=mycalloc(1, sho_t);
	const real omega0=2*M_PI*f0;
	out->dt=0.01/f0;
	out->c1=2*zeta*omega0;
	out->c2=omega0*omega0;
	return out;
}
/**
 * */
void sho_free(sho_t *sho){
	if(!sho) return;
	cellfree(sho->dx);
	cellfree(sho->x);
	cellfree(sho->ddx);
	cellfree(sho->ytmp);
	free(sho);
}
/**
   Second harmonic oscillator state space update.
   xi is driving position. 
 */

/**
   Second harmonic oscillator state space update.
   xi is driving position. Update the position after time dt.
   Also compute the average position between 0 and dt (for WFS integration) if avg is set.
	
	//sho->x +=dt*sho->dx;//update position
	//ddx=sho->dx*(-sho->c1)+sho->x*(-sho->c2)+xi;//update speed
	//sho->dx+=dt*ddx;	 //update speed
 */
void sho_step(panyarray xout_, sho_t *sho, anyarray xi_, real dt, int avg){
	cell* xi=xi_.c;
	cell** xout=xout_.c;
	if(!sho){//null filter. just copy input to output.
		dcelladd(xout, 0, xi, 1);
		return;
	}
	//divide dt to multiple time to do proper integration.
	long nover=(long)ceil(dt/sho->dt);
	if(nover==1) avg=0;//no need for averaging.
	real scale=1./nover;
	dt*=scale;
	if(avg)	dcellscale(*xout, 0);
	for(long i=0; i<nover; i++){
		dcelladd(&sho->ddx, 0, sho->dx, -sho->c1);
		dcelladd(&sho->ddx, 1, sho->x,  -sho->c2);
		dcelladd(&sho->ddx, 1, xi,       sho->c2);
		dcelladd(&sho->x,   1, sho->dx,  dt);//update position
		dcelladd(&sho->dx,  1, sho->ddx, dt);//update position
		if(avg) dcelladd(xout, 1, sho->x, scale);
	}
	if(!avg){
		dcelladd(xout, 0, sho->x, 1);
	}
}
/**
   Second harmonic oscillator. Reset.
*/
void sho_reset(sho_t* sho){
	dcellscale(sho->dx, 0);
	dcellscale(sho->x, 0);
}
/**
   Second harmonic oscillator. Filter a time series for testing.
*/
dmat* sho_filter(const dmat* x,/**<Input time series*/
	real dt,     /**<Input time series sampling*/
	real f0,    /**<Resonance frequency*/
	real zeta,  /**<Damping*/
	int avg /**<Average or not*/
	){
	sho_t* sho=sho_new(f0, zeta);
	dmat* xi=dref(x);
	if(NX(xi)==1){
		reshape(xi, NY(xi), 1);
	}
	dmat* yo=dnew(NX(xi), NY(xi));
	dcell *xtmp=dcellnew_same(1,1,1,1);
	dcell *ytmp=dcellnew_same(1, 1, 1, 1);
	for(long iy=0; iy<NY(xi); iy++){
		sho_reset(sho);
		for(long ix=0; ix<NX(xi); ix++){
			P(xtmp,0,0,0,0)=P(xi, ix, iy);
			sho_step(&ytmp, sho, xtmp, dt, avg);
			P(yo, ix, iy)=P(ytmp,0,0,0,0);
		}
	}
	dcellfree(xtmp);
	dcellfree(ytmp);
	dfree(xi);
	sho_free(sho);
	return yo;
}
