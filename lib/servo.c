/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

typedef struct SERVO_CALC_T{
    cmat *s;
    cmat *Hol;      /**<End-to-End open loop transfer function*/
    cmat *Hint;     /**<Descrete Integrator transfer function*/
    cmat *Hsys;     /**<System transfer function, including Hwfs, Hlag, Hdac, Hmir, Hint. without gain*/
    cmat *Hwfs;
    dmat *nu;       /**<[out] The frequency*/
    dmat *psd;      /**<[out] The PSD defined on nu*/
    double sigma2n;  /**<[out] Noise variance*/
    double pmargin; /*phase margin. default: M_PI/4*/
    double fny;     /**Nyqust frequency for dtrat*/
    double var_sig; /**Integration of PSD*/
    /*Output. */
    double res_sig; /**<[out] Residual signal*/
    double res_n;   /**<[out] Residual noise propagation*/
    double gain_n;  /**<[out] noise amplification*/
    double g;
    double a;
    double T;
    int type;
}SERVO_CALC_T;
#define TWOPI 6.283185307179586
/*Compute phase between -2pi, and 0*/
INLINE double phase(dcomplex val){
    double ang=atan2(cimag(val), creal(val));
    if(ang>0){
	ang=ang-TWOPI;
    }
    return ang;
}
/**
   With given open loop transfer function, compute its cross over frequency and
   phase margin.

   Written: 2010-06-11

   Tested OK: 2010-06-11
*/
/*Determine the phase difference between Hol and -180 when abs(Hol)==1*/
static double phase_at_gain(double *fcross, /**<[out] Cross over frequency*/
			    const dmat *nu, /**<[in] frequency grid*/
			    const cmat *Hol, /**<[in] open loop transfer function defined on nu*/
			    double gain/**<[in] compute phase at this gain*/
    ){
    /*if(cabs(Hol->p[0])<gain){
	error("Hol is less than 1 from the beginning\n");
	}*/

    int found=0;
    double phi=0;
    for(long i=1; i<nu->nx; i++){
	double val=cabs(Hol->p[i]);
	if(val<gain){
	    double valpre=cabs(Hol->p[i-1]);
	    double logvalpre=log10(valpre);/*positive. log10(val) is negative. */
	    double rat=(logvalpre-log10(gain))/(logvalpre-log10(val));
	    *fcross=pow(10,log10(nu->p[i-1])*(1-rat)+log10(nu->p[i])*rat);
	    double phi0=phase(Hol->p[i-1]);
	    double phi1=phase(Hol->p[i]);
	    //Unwrap the differeance
	    double diff=(phi1-phi0)/(2*M_PI);/*positive. */
	    diff=diff-round(diff);
	    phi1=phi0+diff*2*M_PI;
	    phi=phi0*(1-rat)+phi1*rat;
	    //info2("valpre=%g, val=%g, phi0=%g, phi1=%g, phi=%g\n", valpre, val, phi0, phi1, phi);
	    found=1;
	    break;
	}
    }
    if(!found){/*warning("Hol does not decrease to 1.\n");*/
	phi=phase(Hol->p[nu->nx-1]);
	*fcross=nu->p[nu->nx-1];
    }
    double nphi=phi/(2*M_PI);
    nphi=nphi-floor(nphi)-1;
    phi=nphi*2*M_PI;/*from -2*M_PI to 0. */
    return phi+M_PI;/**/
}

/*Determine the ratio in dB between Hol and 1 when phase of Hol is -180*/
static double gain_at_phase(double *fcross, /**<[out] Cross over frequency*/
			    const dmat *nu, /**<[in] frequency grid*/
			    const cmat *Hol, /**<[in] open loop transfer function defined on nu*/
			    double angle     /**<[in] compute gain at this phase*/
			       ){
    if(angle<-TWOPI || angle>0){
	error("angle=%g is invalid\n", angle);
    }
    long i;
    //Skip leading terms if they are belowangle
    for(i=0; i<nu->nx; i++){
	if(phase(Hol->p[i])>angle){
	    break;
	}
    }
    int found=0;
    double gain=0;
    for(; i<nu->nx; i++){
	double phi1=phase(Hol->p[i]);
	if(phi1<angle){
	    double phi0=phase(Hol->p[i-1]);/*how much above -pi*/
	    double rat=(phi0+M_PI)/(phi0-phi1);
	    double val=cabs(Hol->p[i]);
	    double valpre=cabs(Hol->p[i-1]);
	    *fcross=nu->p[i-1]+(nu->p[i]-nu->p[i-1])*rat;
	    gain=-20*log10(valpre+(val-valpre)*rat);
	    found=1;
	    break;
	}
    }
    if(!found){/*Loop is unstable*/
	*fcross=0;
	gain=NAN;
    }
    /*{
	info2("gain=%g, angle=%g\n", gain, angle);
	static int saved=0;
	if(!saved){
	    saved=1;
	    cwrite(Hol, "Hol");
	    dwrite(nu, "nu");
	}
	}*/
    return gain;/**/
}

/**
   Make basic arrays for servo analysis.
*/
static void servo_calc_init(SERVO_CALC_T *st, const dmat *psdin, double dt, long dtrat){
    if(psdin->ny!=2){
	error("psdin should have two columns\n");
    }
    double Ts=dt*dtrat;
    st->fny=0.5/Ts;
    dmat *nu=st->nu=dlogspace(-3,log10(0.5/dt),1000);
    st->psd=dinterp1(psdin, 0, nu);
    st->var_sig=psd_inte2(psdin);
    dcomplex pi2i=TWOPI*I;
    if(st->Hsys || st->Hwfs || st->Hint || st->s){
	error("Already initialized\n");
    }
    st->Hsys=cnew(nu->nx, 1);
    st->Hwfs=cnew(nu->nx, 1);
    st->Hint=cnew(nu->nx, 1);
    st->s=cnew(nu->nx,1);
    for(long i=0; i<nu->nx; i++){
	dcomplex s=st->s->p[i]=pi2i*nu->p[i];
	dcomplex zInv=cexp(-s*Ts);
	dcomplex Hint=st->Hint->p[i]=1./(1-zInv);
	dcomplex Hwfs, Hdac;
	if(dtrat==1){//we have a pure delay
	    Hwfs=st->Hwfs->p[i]=zInv; 
	    Hdac=1;
	}else{
	    Hwfs=st->Hwfs->p[i]=(1-zInv)/(Ts*s);
	    Hdac=Hwfs;
	}
	dcomplex Hlag=cexp(-s*dt);/*lag due to readout/computation*/
	dcomplex Hmir=1;/*DM */
	st->Hsys->p[i]=Hwfs*Hlag*Hdac*Hmir*Hint;
    }
}
static void servo_calc_free(SERVO_CALC_T *st){
    cfree(st->s);
    cfree(st->Hol);
    cfree(st->Hint);
    cfree(st->Hsys);
    cfree(st->Hwfs);
    dfree(st->nu);
    dfree(st->psd);
}
static int servo_isstable(dmat *nu, cmat *Hol){
    double fp, fg;
    double pmargin=phase_at_gain(&fp, nu, Hol, 1);
    double gmargin=gain_at_phase(&fg, nu, Hol, -M_PI);
    int isstable=pmargin>M_PI/4.1 && fp<fg && gmargin>0;
    /*if(!isstable){
       info2("Unstable: phase margin: %4.0f deg at %3.0f Hz. Gain margin: %2.0fdB at %3.0f Hz.\n",
	     pmargin*180/M_PI, fp, gmargin, fg);
	     }*/
    return isstable;
}
/*
static double servo_calc_max_gain(SERVO_CALC_T *st, double pmargin){
    if(!st->Hol){
	st->Hol=cnew(st->nu->nx,1);
    }
    cadd(&st->Hol, 0, st->Hsys, 1);
    double angle=pmargin-M_PI;
    if(st->type==2){
	ccwm(st->Hol, st->Hint);
	angle-=M_PI/2.4;//allow more due to lead filter
    }
    double fcross;
    double dB=gain_at_phase(&fcross, st->nu, st->Hol, angle);
    double gain=isfinite(dB)?pow(10, dB/20):0;
    return gain;
    }*/
/**
   Calculate total error. If g0 is 0, use st->g, otherwith use g0 to figure out g, a T.
*/
static double servo_calc_do(SERVO_CALC_T *st, double g0){
    dmat *nu=st->nu;
    if(!st->Hol){
	st->Hol=cnew(nu->nx,1);
    }
    /*Compute Hol with the first integrator and gain.*/
    if(g0>EPS){
	st->g=g0;
    }
    cadd(&st->Hol, 0, st->Hsys, st->g);
    double g2=1;/*additional g to multiply. !=1 if g0 is nonzero and type is 2.*/
    if(st->type==2){
	ccwm(st->Hol, st->Hint);/*multiply the second integrator*/
        if(fabs(g0)>EPS){/*figure out a, T from new g0*/
	    double margin, fcross;
	    margin=phase_at_gain(&fcross, nu, st->Hol, 1);/*see how much phase lead is needed*/
	    double phineed=st->pmargin-margin;
	    double a,T;
	    if(phineed*2.2>M_PI){/*lead filter is not suitable*/
		a=1;
		T=0;
	    }else{
		 a=(1-sin(phineed))/(1+sin(phineed));
		 double f0=fcross*sqrt(a);
		 T=1./(2.*M_PI*f0);
	    }
	    /* Hlead is multiplied by sqrt(a) so it has unit gain at cross over frequency*/
	    g2=sqrt(a);
	    st->g=g0*g2;
	    st->a=a;
	    st->T=T;
	    //info("g0=%g, g2=%g, phineed=%.1f\n", g0, g2, phineed*180/M_PI);
	}
	double a=st->a;
	double T=st->T;
	for(int i=0; i<nu->nx; i++){
	    dcomplex Hlead=(1+T*st->s->p[i])/(1+a*T*st->s->p[i])*g2;
	    st->Hol->p[i]*=Hlead;
	}
    }
    double res_sig=0;
    double sum_n=0, sum_1=0;
    dmat *psd=st->psd;
    for(int i=0; i<nu->nx; i++){
	dcomplex Hol=st->Hol->p[i];
	dcomplex Hrej=1./(1.+Hol);
	//The gain never reach below -50dB
	res_sig+=psd->p[i]*creal(Hrej*conj(Hrej)+1e-5)*nu->p[i];
	if(st->nu->p[i]<st->fny){
	    //compute noise prop only within nyqust frequency
	    dcomplex Hcl=Hol*Hrej;
	    dcomplex Hwfs=st->Hwfs->p[i];
	    dcomplex Hn=Hcl/Hwfs;
	    sum_n+=creal(Hn*conj(Hn))*nu->p[i];
	    sum_1+=nu->p[i];
	}
    }
    double dlognu=(log(nu->p[nu->nx-1])-log(nu->p[0]))/(nu->nx-1);
    res_sig*=dlognu;
    st->res_sig=res_sig;
    st->gain_n=sum_n/sum_1;
    if(!servo_isstable(nu, st->Hol)){
	st->gain_n=100;
	/*put a high penalty to drive down the gain*/
	st->res_n=10*(1+g0)*(st->var_sig+st->sigma2n);
	/*warning2("Unstable: g0=%g, g2=%g, res_sig=%g, res_n=%g, tot=%g, gain_n=%g sigma2n=%g\n",
		 g0, g2, res_sig, st->res_n, st->res_n+res_sig, st->gain_n, st->sigma2n);*/
    }else{
	if(st->gain_n>1){
	    st->gain_n=pow(st->gain_n,3);/*a fudge factor to increase the penalty*/
	}
        st->res_n=st->sigma2n*st->gain_n;
        /*info2("  Stable: g0=%g, g2=%g, res_sig=%g, res_n=%g, tot=%g, gain_n=%g sigma2n=%g\n",
	  g0, g2, res_sig, st->res_n, st->res_n+res_sig, st->gain_n, st->sigma2n);*/
    }
    return res_sig+st->res_n;
}

/**
   Optimize the type II servo gains by balancing errors due to noise and
   signal propagation.

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
   Returns a cellarray of a dmat of [g0, a, T, res_sig, res_n]
*/
dcell* servo_optim(const dmat *psdin,  double dt, long dtrat, double pmargin,
		   const dmat* sigma2n, int servo_type){
    /*The upper end must be nyquist freq so that noise transfer can be
      computed. But we need to capture the turbulence PSD beyond nyquist freq,
      which are uncorrectable.
    */
    SERVO_CALC_T st={0};
    servo_calc_init(&st, psdin, dt, dtrat);
    st.type=servo_type;
    st.pmargin=pmargin;
    int ng=1;
    switch(servo_type){
    case 1:
	ng=1;break;
    case 2:
	ng=3; break;
    default:
	error("Invalid servo_type=%d\n", servo_type);
    }
    dcell *gm=dcellnew(sigma2n->nx, sigma2n->ny);
    double g0_step=1e-6;
    double g0_min=1e-6;/*the minimum gain allowed.*/
    double g0_max=2.0;
    for(long ins=0; ins<sigma2n->nx*sigma2n->ny; ins++){
	st.sigma2n=sigma2n->p[ins];
	double g0=golden_section_search((golden_section_fun)servo_calc_do, &st, g0_min, g0_max, g0_step);
	servo_calc_do(&st, g0);
	gm->p[ins]=dnew(ng+2,1);
	gm->p[ins]->p[0]=st.g;
	if(servo_type==2){
	    gm->p[ins]->p[1]=st.a;
	    gm->p[ins]->p[2]=st.T;
	}
	gm->p[ins]->p[ng]=st.res_sig;
	gm->p[ins]->p[ng+1]=st.res_n;
	/*info2("g0=%.1g, g2=%.1g, res_sig=%.1g, res_n=%.1g, tot=%.1g, gain_n=%.1g sigma2n=%.1g\n",
	  g0, st.g, st.res_sig, st.res_n, st.res_n+st.res_sig, st.gain_n, st.sigma2n);*/
    }/*for ins. */
    servo_calc_free(&st);
    return gm;
}
/**
   Compute the residual error after servo rejection of a type II servo given the gain.
   
   Written: 2010-06-11
   
   Tested OK: 2010-06-11
*/
 double servo_residual(double *noise_amp, const dmat *psdin, double dt, long dtrat, const dmat *gain, int servo_type){
     SERVO_CALC_T st={0};
    servo_calc_init(&st, psdin, dt, dtrat);
    st.type=servo_type;
    switch(servo_type){
    case 1:
	st.g=gain->p[0];
	break;
    case 2:
	st.g=gain->p[0];
	st.a=gain->p[1];
	st.T=gain->p[2];
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
servo_typeII_filter(SERVO_T *st, dcell *merrc){
    if(!merrc) return;
    const dmat *gain=st->ep;
    PDMAT(gain,pgain);
    int indmul=0;
    if(gain->nx!=3){
	error("Wrong format in gain\n");
    }
    double gg,e1a, e1;
    for(int ic=0; ic<merrc->nx*merrc->ny; ic++){
	dmat *merr=merrc->p[ic];
	dmat *mlead=st->mlead->p[ic];
	dmat *merrlast=st->merrlast->p[ic];
	int nmod=0;/*error. */
	if(merr->ny==1){
	    nmod=merr->nx;
	}else{
	    if(merr->nx!=1){
		error("Don't handle this case\n");
	    }
	    nmod=merr->ny;
	}
	if(gain->ny==1){
	    indmul=0;
	}else if(gain->ny==nmod){
	    indmul=1;
	}else{
	    error("Wrong format in gain\n");
	}
	/**2013-12-03: fix lead filter implementation. gg is relocated
	   2013-12-05: A more accurate and robust lead filter
	 */
	for(int imod=0; imod<nmod; imod++){
	    int indm=imod * indmul;
	    gg=pgain[indm][0];
	    e1a=pgain[indm][1];
	    e1=pgain[indm][2];
	    mlead->p[imod] = e1a*mlead->p[imod]+gg*(1-e1a)/(1-e1)*(merr->p[imod]-e1*merrlast->p[imod]);
	}
    }
    dcellcp(&st->merrlast, merrc);
    dcelladd(&st->mpreint,1, st->mlead, 1);
}
static void servo_init(SERVO_T *st, dcell *merr){
    if(!merr || st->initialized){
	error("merr must be valid and SERVO_T must be not yet initialized\n");
    }
    if(st->ep->nx>1){
	st->mpreint=dcellnew2(merr);
    }
    if(st->ep->nx==3){
	st->mlead=dcellnew2(merr);
	st->merrlast=dcellnew2(merr);
    }
    for(int i=0; i<st->nmint; i++){
	st->mint[i]=dcellnew2(merr); 
    }
    st->initialized=1;
}
/**
   Initialize. al is additional latency
*/
SERVO_T *servo_new(dcell *merr, const dmat *ap, int al, double dt, const dmat *ep){
    SERVO_T *st=calloc(1, sizeof(SERVO_T));
    st->nmint=ap?MAX(2,ap->nx):1;
    st->mint=calloc(st->nmint, sizeof(dcell));
    st->ap=dref(ap);
    if(ep->nx!=3){//type I
	st->ep=dref(ep);
    }else{//type II. preprocess
	st->ep=dnew(ep->nx, ep->ny);
	for(int i=0; i<ep->ny; i++){
	    st->ep->p[i*3]=ep->p[i*3];
	    double a=ep->p[1+i*3];
	    double T=ep->p[2+i*3];
	    st->ep->p[1+i*3]=exp(-dt/(a*T));
	    st->ep->p[2+i*3]=exp(-dt/T);
	}
    }
    st->dt=dt;
    st->al=al;
    st->merrhist=calloc(st->al, sizeof(dcell*));
    if(merr && merr->nx!=0 && merr->ny!=0 && merr->p[0]){
	servo_init(st, merr);
    }
    return st;
}
/**
   prepare the integrator by shifting commands. similar to laos.
   inte->p[0]=inte->p[0]*ap[0]+inte->p[1]*ap[1]+...
*/
static void servo_shift_ap(SERVO_T *st){
    const dmat *ap=st->ap;
    if(!ap) return; //no need to shift.
    if(st->nmint<ap->nx){
	st->nmint=ap->nx;
	st->mint=realloc(st->mint, sizeof(dcell*));
    }
    if(!st->initialized) return;
    dcell **inte=st->mint;
    dcell *tmp=NULL;
    dcell *keepjunk=inte[ap->nx-1];
    for(int iap=ap->nx-1; iap>=0; iap--){
	dcelladd(&tmp,1,inte[iap],ap->p[iap]);
	if(iap>0){
	    inte[iap]=inte[iap-1];/*shifting */
	}else{
	    inte[iap]=tmp;/*new command. */
	}
    }
    dcellfree(keepjunk);
}
/*A FIFO queue to add delay*/
static dcell*servo_shift_al(SERVO_T *st, dcell *merr){
    if(!st->al){
	return dcellref(merr);
    }else{
	dcell *out=st->merrhist[0];
	for(int i=0; i<st->al-1; i++){
	    st->merrhist[i]=st->merrhist[i+1];
	}
	st->merrhist[st->al-1]=dcelldup(merr);
	return out;
    }
}
/**
   Applies type I or type II filter based on number of entries in gain.
*/
int servo_filter(SERVO_T *st, dcell *_merr){
    if(!st->initialized && _merr){
	servo_init(st, _merr);
    }
    dcell *merr=servo_shift_al(st, _merr);
    if(!merr) return 0;
    servo_shift_ap(st);
    if(!st->mint){
	error("SERVO_T must be created using servo_new()\n");
    }
    switch(st->ep->nx){
    case 1://type I
	if(st->ep->ny==1){
	    dcelladd(&st->mpreint, 0, merr, st->ep->p[0]);//just record what is added.
	}else{
	    if(!st->mpreint){
		st->mpreint=dcellnew2(merr);
	    }
	    for(int ic=0; ic<merr->nx; ic++){
		assert(merr->p[ic]->nx==st->ep->ny);
		for(long i=0; i<merr->p[ic]->nx; i++){
		    st->mpreint->p[ic]->p[i]=st->ep->p[i]*merr->p[ic]->p[i];
		}
	    }
	}
	break;
    case 2:{//PID controller
	if(st->ep->ny!=1) error("not supported\n");
	double g1=st->ep->p[0]+st->ep->p[1];
	double g2=-st->ep->p[1];
	dcelladd(&st->mpreint, 0, merr, g1);
	dcelladd(&st->mpreint, 1, st->merrlast, g2);
	dcellcp(&st->merrlast, merr);
    }
	break;
    case 3://type II
	servo_typeII_filter(st, merr);
	break;
    default:
	error("Invalid: st->ep->nx=%ld", st->ep->nx);
    }
    dcelladd(st->mint, 1, st->mpreint, 1);
    dcellfree(merr);
    return 1;
}


/**
   test type I/II filter with ideal measurement to make sure it is implemented correctly.
*/
dmat* servo_test(dmat *input, double dt, int dtrat, dmat *sigma2n, dmat *gain){
    if(input->ny==1){/*single mode. each column is for a mode.*/
	input->ny=input->nx;
	input->nx=1;
    }
    int nmod=input->nx;
    PDMAT(input,pinput);
    dmat *merr=dnew(nmod,1);
    dcell *mreal=dcellnew(1,1);
    dmat *mres=dnew(nmod,input->ny);
    dmat *sigman=NULL;
    if(dnorm2(sigma2n)>0){
	sigman=dchol(sigma2n);
    }
    dcell *meas=dcellnew(1,1);
    dmat *noise=dnew(nmod, 1);
    SERVO_T *st2t=servo_new(NULL, NULL, 0, dt*dtrat, gain);
    rand_t rstat;
    seed_rand(&rstat, 1);
    PDMAT(mres,pmres);
    /*two step delay is ensured with the order of using, copy, acc*/
    for(int istep=0; istep<input->ny; istep++){
	memcpy(merr->p, pinput[istep], nmod*sizeof(double));
	dadd(&merr, 1, mreal->p[0], -1);
	memcpy(pmres[istep],merr->p,sizeof(double)*nmod);
	if(istep % dtrat == 0){
	    dzero(meas->p[0]);
	}
	dadd(&meas->p[0], 1, merr, 1);/*average the error. */
	dcellcp(&mreal, st2t->mint[0]);
	if((istep+1) % dtrat == 0){
	    if(dtrat!=1) dscale(meas->p[0], 1./dtrat);
	    if(sigman){
		drandn(noise, 1, &rstat);
		if(sigman->nx>0){
		    dmm(&meas->p[0], 1, sigman, noise, "nn", 1);
		}else{
		    dadd(&meas->p[0], 1, noise, sigman->p[0]);
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
void servo_reset(SERVO_T *st){
    dcellzero(st->mlead);
    dcellzero(st->merrlast);
    dcellzero(st->mpreint);
    if(st->merrhist){
	for(int i=0; i<st->al; i++){
	    dcellzero(st->merrhist[i]);
	}
    }
    if(st->mint){
	for(int i=0; i<st->nmint; i++){
	    dcellzero(st->mint[i]);
	}
    }
}
/**
   Free SERVO_T struct
*/
void servo_free(SERVO_T *st){
    if(!st) return;
    dcellfree(st->mlead);
    dcellfree(st->merrlast);
    dcellfree(st->mpreint);
    dcellfreearr(st->merrhist, st->al);
    dcellfreearr(st->mint, st->nmint);
    dfree(st->ap);
    dfree(st->ep);
    free(st);
}
/**
   Integrated a PSF that defines on linear or logrithmically spaced grid nu.
*/
double psd_inte(const double *nu, const double *psd, long n){
    double dnu=(nu[n-1]-nu[0])/(n-1);
    double dlognu=(log(nu[n-1])-log(nu[0]))/(n-1);
    double res_sig=0;
    if(fabs(nu[1]-nu[0]-dnu)<dnu*1.e-4){
	for(long i=0; i<n; i++){
	    res_sig+=psd[i];
	}
	res_sig*=dnu;
    }else if((log(nu[1])-log(nu[0])-dlognu)<dlognu*1.e-4){
        for(long i=0; i<n; i++){
	    res_sig+=psd[i]*nu[i];
	}	
	res_sig*=dlognu;
    }
    return res_sig;
}
/**
   wraps psd_inte
*/
double psd_inte2(const dmat *psdin){
    if(psdin->ny!=2){
	error("psdin  should have two columns\n");
    }
    double *nu=psdin->p;
    long n=psdin->nx;
    double *psd=nu+n;
    return psd_inte(nu, psd, n);
}


/**
   Convert PSD into time series.*/
dmat* psd2time(const dmat *psdin, rand_t *rstat, double dt, int nstepin){
    if(psdin->ny!=2){
	error("psdin should have two columns\n");
    }
    long nstep=nextpow2(nstepin);
    double df=1./(dt*nstep);
    dmat *fs=dlinspace(0, df, nstep);
    dmat *psd=NULL;
    double var=psd_inte2(psdin);
    psd=dinterp1(psdin, 0, fs);
    psd->p[0]=0;/*disable pistion. */
    cmat *wshat=cnew(nstep, 1);
    cfft2plan(wshat, -1);
    for(long i=0; i<nstep; i++){
	wshat->p[i]=sqrt(psd->p[i]*df)*(randn(rstat)+I*randn(rstat));
    }
    cfft2(wshat, -1);
    dmat *out=NULL;
    creal2d(&out, 0, wshat, 1);
    cfree(wshat);
    dfree(psd);
    dfree(fs);
    dresize(out, nstepin, 1);
    double var2=dinn(out,out)/out->nx;
    //info2("Input psd has variance of %g.\nTime series has variance of %g\n",var,var2);
    dscale(out, sqrt(var/var2));
    return out;
}

/**
   Add two PSDs that doesn't have the same frequency. the first column of each
   dmat is the frequency nu, and the second column is PSD. Bug discovered on
   2013-03-24:only psd2 was added to to psd.*/
static dmat *add_psd_nomatch(const dmat *psd1,const dmat *psd2){
    dmat *nu1=dsub(psd1,0,psd1->nx,0,1);
    dmat *p2ynew=dinterp1(psd2, 0, nu1);
    dmat *psd=dnew(nu1->nx,2);
    double *py=psd->p+psd->nx;
    const double *p1y=psd1->p+psd1->nx;
    for(long i=0; i<psd->nx; i++){
	psd->p[i]=nu1->p[i];
	py[i]=p1y[i]+p2ynew->p[i];
    }
    dfree(nu1);
    dfree(p2ynew);
    return psd;
}
/**
   Add two PSDs. the first column of each dmat is the frequency nu, and the
   second column is PSD*/
dmat *add_psd(const dmat *psd1, const dmat *psd2){
    if(psd1->nx!=psd2->nx){
	//warning("The two PSDs have different length\n");
	return add_psd_nomatch(psd1, psd2);
    }
    dmat *psd=dnew(psd1->nx,2);
    double *restrict pp=psd->p+psd->nx;
    const double *p1=psd1->p+psd1->nx;
    const double *p2=psd2->p+psd2->nx;
    for(long i=0; i<psd->nx; i++){
	if(fabs(psd1->p[i]-psd2->p[i])>1.e-2){
	    warning("The two PSDs doesn't have the same freq.");
	    dfree(psd);
	    return add_psd_nomatch(psd1,psd2);
	    /*todo: apply interp1 to interpolate the second PSD. */
	}
	psd->p[i]=psd1->p[i];
	pp[i]=p1[i]+p2[i];
    }
    return psd;
}
/*
  Add a PSD to another.
*/
void add_psd2(dmat **out, const dmat *in){
    if(!*out){
	*out=ddup(in);
    }else{
	dmat *tmp=*out;
	*out=add_psd(tmp, in);
	dfree(tmp);
    }
}

/**
   contains data related to DM hysterisis modeling for all the common DMs (not
MOAO). let input be x, and output of each mode be y.  */
struct HYST_T{
    dmat *coeff;      /**<contains data from parms->dm.hyst*/
    dmat *xlast;      /**<Record x from last time step*/
    dmat *ylast;      /**<Record y from last time step*/
    dmat *dxlast;     /**<Change of x*/
    dmat *x0;         /**<Initial x (before dx changes sign)*/
    dmat *y0;         /**<Initial y*/
    double stroke;    /**<Surface stroke of DM. coeff->alpha is scaled by this*/
    double power;     /**<Slope of changing alpha*/
    double xscale;    /**<Scaling of x to fix overal slope error.*/
};
/*
  Hysteresis modeling. Create the hysteresis model
*/
HYST_T *hyst_new(dmat *coeff, int naloc){
    HYST_T *hyst=calloc(1, sizeof(HYST_T));
    int nhmod=coeff->ny;
    if(coeff->nx!=3 || nhmod<1){
	error("DM hystereis has wrong format. Expect 3 rows, but has %ldx%ld\n", coeff->nx, coeff->ny);
    }
    hyst->coeff=dref(coeff);
    hyst->xlast=dnew(naloc,1);
    hyst->ylast=dnew(nhmod,naloc);
    hyst->dxlast=dnew(naloc,1);
    hyst->x0=dnew(naloc,1);
    hyst->y0=dnew(nhmod,naloc);
    hyst->xscale=1;
    if(coeff->header){
	double stroke=search_header_num(coeff->header, "stroke");
	if(isfinite(stroke)){
	    warning("stroke=%g\n", stroke);
	    for(int i=0; i<coeff->ny; i++){
		coeff->p[i*coeff->nx+1]*=stroke;
	    }
	    hyst->stroke=stroke;
	}
	double alpha_power=search_header_num(coeff->header, "alpha_power");
	if(isfinite(alpha_power)){
	    warning("alpha_power=%g\n", alpha_power);
	    hyst->power=alpha_power;
	}
    }
    return hyst;
}
void hyst_reset(HYST_T *hyst){
    dzero(hyst->xlast);
    dzero(hyst->ylast);
    dzero(hyst->dxlast);
    dzero(hyst->x0);
    dzero(hyst->y0);
}
void hyst_free(HYST_T *hyst){
    dfree(hyst->coeff);
    dfree(hyst->xlast);
    dfree(hyst->ylast);
    dfree(hyst->dxlast);
    dfree(hyst->x0);
    dfree(hyst->y0);
    free(hyst);
}
void hyst_dmat(HYST_T *hyst, dmat *dmreal, const dmat *dmcmd){
    double *restrict x=dmcmd->p;
    double *restrict xout=dmreal->p;
    double *restrict xlast=hyst->xlast->p;
    double *restrict dxlast=hyst->dxlast->p;
    double *restrict x0=hyst->x0->p;
    PDMAT(hyst->ylast, ylast);
    PDMAT(hyst->y0, py0);
    PDMAT(hyst->coeff, coeff);
    int nmod=hyst->coeff->ny;
    int naloc=dmcmd->nx;
    for(int ia=0; ia<naloc; ia++){
	double xia=x[ia]*hyst->xscale;
	double dx=xia-xlast[ia];
	if(fabs(dx)>1e-14){/*There is change in command */
	    if(dx*dxlast[ia]<0){
		/*Changes in moving direction, change the initial condition */
		x0[ia]=xlast[ia];
		for(int imod=0; imod<nmod; imod++){
		    py0[ia][imod]=ylast[ia][imod];
		}
	    }
	    double alphasc=dx>0?1:-1;/*To revert the sign of alpha when dx<0 */
	    if(hyst->power){
		alphasc*=pow((hyst->stroke-x[ia])/(hyst->stroke*2.), hyst->power);
	    }
	    for(int imod=0; imod<nmod; imod++){
		const double alpha=alphasc*coeff[imod][1];
		const double alphabeta=alpha*coeff[imod][2];
		ylast[ia][imod]=xia-alphabeta+(py0[ia][imod]-x0[ia]+alphabeta)*exp(-(xia-x0[ia])/alpha);
	    }
	    xlast[ia]=xia;
	    dxlast[ia]=dx;
	}/*else: no change in voltage, no change in output. */
	/*update output. */
	double y=0;
	for(int imod=0; imod<nmod; imod++){
	    y+=ylast[ia][imod]*coeff[imod][0];
	}
	xout[ia]=y;
    }
}
void hyst_dcell(HYST_T **hyst, dcell *dmreal, const dcell *dmcmd){
    if(!hyst) return;
    for(int idm=0; idm<dmcmd->nx*dmcmd->ny; idm++){
	hyst_dmat(hyst[idm], dmreal->p[idm], dmcmd->p[idm]);
    }
}
/*Calibrate the slope of hysteresis curve.*/
void hyst_calib(HYST_T *hyst, int ih){
    int N=1000;
    dmat *cmd=dnew(N,1);
    dmat *real=dnew(N,1);
    dmat *cmd0=dnew(1,1);
    dmat *real0=dnew(1,1);
    double stroke=hyst->stroke;
    if(!hyst->stroke){
	stroke=10e-6;
    }
    for(int i=0; i<N; i++){
	cmd0->p[0]=cmd->p[i]=sin(i*M_PI*0.01)*stroke;
	hyst_dmat(hyst, real0, cmd0);
	real->p[i]=real0->p[0];
    }
    hyst_reset(hyst);
    hyst->xscale*=(dmax(cmd)-dmin(cmd))/(dmax(real)-dmin(real));
    warning("x would be multiplied by %g\n", hyst->xscale);
    dwrite(real, "hyst%d_real", ih);
    dwrite(cmd,  "hyst%d_cmd", ih);
    for(int i=0; i<N; i++){
	cmd0->p[0]=cmd->p[i];
	hyst_dmat(hyst, real0, cmd0);
	real->p[i]=real0->p[0];
    }
    hyst_reset(hyst);
    dwrite(real, "hyst%d_realcalib", ih);
    dfree(cmd);
    dfree(real);
    dfree(cmd0);
    dfree(real0);
}
