/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "mathmisc.h"
#define CONTINUOUS 1 /*1: use real continuous/discrete AO servo model. 0: Assume 2 cycle delay + discrete integrator*/
/**
   \file servo.c
   Routines for servo optimization, filtering, etc.
*/


typedef struct SERVO_CALC_T{
    cmat *s;
    cmat *Hol;      /**<End-to-End open loop transfer function*/
    cmat *Hint;     /**<Descrete Integrator transfer function*/
    cmat *Hsys;     /**<System transfer function, including Hwfs, Hlag, Hdac, Hmir, Hint. without gain*/
    cmat *Hwfs;
    dmat *nu;       /**<[out] The frequency*/
    dmat *psd;      /**<[out] The PSD defined on nu*/
    double sigman;  /**<[out] Noise vairance*/
    double pmargin; /*phase margin. default: M_PI/4*/
    /*Output. */
    double rms2_sig;/**<[out] Uncorrectable signal*/
    double rms_sig; /**<[out] signal pass through*/
    double rms_n;
    double gain_n;  /**<[out] noise amplification*/
    double g;
    double a;
    double T;
    int type;
}SERVO_CALC_T;

/**
   With given open loop transfer function, compute its cross over frequency and
   phase margin.

   Written: 2010-06-11

   Tested OK: 2010-06-11
*/

static double servo_phi_margin(double *fcross, /**<[out] Cross over frequency*/
			       const dmat *nu, /**<[out] frequency grid*/
			       const cmat *Hol /**<[in] open loop transfer function defined on nu*/
			       ){
    if(cabs(Hol->p[0])<1){
	error("Hol is less than 1 from the beginning\n");
    }

    int found=0;
    double fc=0;
    double phi=0;
    for(long i=0; i<nu->nx; i++){
	double val=cabs(Hol->p[i]);
	if(val<1){
	    double valpre=cabs(Hol->p[i-1]);
	    double logvalpre=log10(valpre);/*positive. log10(val) is negative. */
	    double rat=logvalpre/(logvalpre-log10(val));
	    fc=pow(10,log10(nu->p[i-1])*(1-rat)+log10(nu->p[i])*rat);
	    double phi0=atan2(cimag(Hol->p[i-1]), creal(Hol->p[i-1]));
	    double phi1=atan2(cimag(Hol->p[i]), creal(Hol->p[i]));
	    /*if(fabs(phi0-phi1)>2*pi) */
	    double diff=(phi1-phi0)/(2*M_PI);/*positive. */
	    diff=diff-round(diff);
	    phi1=phi0+diff*2*M_PI;
	    phi=pow(10,log10(phi0)*(1-rat)+log10(phi1)*rat);
	    double nphi=phi/(2*M_PI);
	    nphi=nphi-floor(nphi)-1;
	    phi=nphi*2*M_PI;/*from -2*M_PI to 0. */
	    found=1;
	    break;
	}
    }
    if(!found){
	cwrite(Hol, "Hol");
	error("Hol does not decrease to 1.\n");
    }
    *fcross=fc;
    return phi+M_PI;/**/
}
/**
   Make basic arrays for servo analysis.
*/
static void servo_calc_init(SERVO_CALC_T *st, const dmat *psdin, double dt, long dtrat){
    double Ts=dt*dtrat;
    double fs=1./Ts;
    dmat *psdf=dnew_ref(psdin->nx,1,psdin->p);
    dmat *psdval=dnew_ref(psdin->nx,1,psdin->p+psdin->nx);  
    
    /*compute error in uncorrectable region.*/
    dmat *nu2=dlogspace(log10(fs/2),3,1000);/*Frequencies that no correction can be made. */
    dmat *psd2=dinterp1log(psdf,psdval,nu2);
    st->rms2_sig=psd_intelog(nu2->p, psd2->p, nu2->nx);
    dfree(nu2); dfree(psd2);

    /*Should not go beyond Nyquist freq. Hrej=1 for nu>fs/2. */
    dmat *nu=st->nu=dlogspace(-3,log10(fs/2),1000);
    st->psd=dinterp1log(psdf,psdval,nu);
    dfree(psdf);
    dfree(psdval);
    dcomplex pi2i=2*M_PI*I;
    if(st->Hsys || st->Hwfs || st->Hint || st->s){
	error("Already initialized\n");
    }
    st->Hsys=cnew(nu->nx, 1);
    st->Hwfs=cnew(nu->nx, 1);
    st->Hint=cnew(nu->nx, 1);
    st->s=cnew(nu->nx,1);
    for(long i=0; i<nu->nx; i++){
	dcomplex s=st->s->p[i]=pi2i*nu->p[i];
	dcomplex expsTs=1-cexp(-s*Ts);
#if CONTINUOUS == 1
	dcomplex Hint=st->Hint->p[i]=1./expsTs;
	dcomplex Hwfs=st->Hwfs->p[i]=expsTs/(Ts*s);
	dcomplex Hdac=Hwfs;
	dcomplex Hlag=cexp(-s*dt);/*lag due to readout/computation*/
	dcomplex Hmir=1;/*DM */
	st->Hsys->p[i]=Hwfs*Hlag*Hdac*Hmir*Hint;
#else
	st->Hwfs->p[i]=1;
	st->Hsys->p[i]=cexp(-s*(dt+dt*dtrat));
#endif
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
/**
   Calculate Hol. If g0 is 0, use st->g, otherwith use g0 to figure out g, a T.
*/
static double servo_calc_do(SERVO_CALC_T *st, double g0){
    dmat *nu=st->nu;
    if(!st->Hol){
	st->Hol=cnew(nu->nx,1);
    }
    /*Compute Hol with the first integrator and gain.*/
    if(fabs(g0)>EPS){
	st->g=g0;
    }
    cadd(&st->Hol, 0, st->Hsys, st->g);
    double g2=1;/*additional g to multiply. !=1 if g0 is nonzero and type is 2.*/
    if(st->type==2){
	ccwm(st->Hol, st->Hint);/*multiply the second integrator*/
        if(fabs(g0)>EPS){/*use supplied g0, figure out a, T*/
	    double margin, fcross;
	    margin=servo_phi_margin(&fcross, nu, st->Hol);/*see how much phase lead is needed*/
	    double phineed=st->pmargin-margin;
	    double a=(1-sin(phineed))/(1+sin(phineed));
	    double f0=fcross*sqrt(a);
	    double T=1./(2.*M_PI*f0);
	    /*
	      According to page 22 of
	      http://wwwhome.math.utwente.nl/~meinsmag/dmcs/docs/DMCSn2.pdf A
	      lead filter should have the form: C(s)=k*(1+sT)/(1+sTa).  T is
	      determined by the cross-over frequency. T=1/(2*pi*fcross*sqrt(a));
	      And a is determined by the necessary phase lead. k may take the
	      form of 1, sqrt(a) or a. k is fused to g for
	      backward-compatibility.
	     */
	    g2=sqrt(a);
	    st->g=g0*g2;
	    st->a=a;
	    st->T=T;
	}
	double a=st->a;
	double T=st->T;
	for(int i=0; i<nu->nx; i++){
	    dcomplex Hlead=(1+T*st->s->p[i])/(1+a*T*st->s->p[i])*g2;
	    st->Hol->p[i]*=Hlead;
	}
    }
    double rms_sig=0;
    double sum_n=0, sum_1=0;
    dmat *psd=st->psd;
    for(int i=0; i<nu->nx; i++){
	dcomplex Hol=st->Hol->p[i];
	dcomplex Hwfs=st->Hwfs->p[i];
	dcomplex Hrej=1./(1.+Hol);
	dcomplex Hcl=Hol*Hrej;
	dcomplex Hn=Hcl/Hwfs;
	rms_sig+=psd->p[i]*creal(Hrej*conj(Hrej))*nu->p[i];
	sum_n+=creal(Hn*conj(Hn))*nu->p[i];
	sum_1+=nu->p[i];
    }
    double dlognu=(log(nu->p[nu->nx-1])-log(nu->p[0]))/(nu->nx-1);
    rms_sig*=dlognu;
    st->rms_sig=rms_sig;
    st->gain_n=sum_n/sum_1;
    st->rms_n=st->sigman*st->gain_n;
    /*info2("g0=%g, g2=%g, rms_sig=%g, rms_n=%g, tot=%g, gain_n=%g\n",
      g0, g2, rms_sig, st->rms_n, st->rms_n+rms_sig, st->gain_n);*/
    return rms_sig+st->rms_n+st->rms2_sig;
}

/**
   Optimize the type II servo gains by balancing errors due to noise and
   signal propagation.

   Written: 2010-06-11

   Tested OK: 2010-06-11

   2011-01-17: The optimization process is quite slow, but the result is only
   dependent on the sigman and fs. psdin does not change during the
   simulation. Built a lookup table in skyc using various sigman and interpolate
   to get ress, resn, and gain.
   
   2012-03-28: Cleaned up this routine. groupped similar calculations together. 
   Change g2=a from g=sqrt(a)
   Limit maximum gain to 0.5. 
   
   sigman is a dmat array of all wanted sigman.
   Returns a cellarray of a dmat of [g0, a, T, res_n, res_sig]
*/
dcell* servo_optim(const dmat *psdin,  double dt, long dtrat, double pmargin,
		   const dmat* sigman, int servo_type){
    /*The upper end must be nyquist freq so that noise transfer can be
      computed. But we need to capture the turbulence PSD beyond nyquist freq,
      which are uncorrectable.
    */
    SERVO_CALC_T st={0};
    servo_calc_init(&st, psdin, dt, dtrat);
    st.type=servo_type;
    st.pmargin=pmargin;

    dcell *gm=dcellnew(sigman->nx, sigman->ny);
    double g0_min=0.001;/*the minimum gain allowed. was 0.001. changed on 2012-03-27*/
    double g0_max=0.5;
    /*{
	st.sigman=sigman->p[0];
	dmat *res=dnew(100,3);PDMAT(res,pres);
	for(int i=0; i<100; i++){
	    double g0=i*0.01+0.001;
	    servo_calc_do(&st, g0);
	    pres[0][i]=g0;
	    pres[1][i]=st.rms_sig;
	    pres[2][i]=st.rms_n;
	}
	dwrite(res, "res");
	}*/

    for(long ins=0; ins<sigman->nx*sigman->ny; ins++){
	st.sigman=sigman->p[ins];
	double g0=golden_section_search((golden_section_fun)servo_calc_do, &st, g0_min, g0_max, 1e-3);
	servo_calc_do(&st, g0);
	gm->p[ins]=dnew(5,1);
	gm->p[ins]->p[0]=st.g;
	gm->p[ins]->p[1]=st.a;
	gm->p[ins]->p[2]=st.T;
	gm->p[ins]->p[3]=st.rms_sig+st.rms2_sig;
	gm->p[ins]->p[4]=st.rms_n;
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
    return st.rms_sig;
}

/**
   Apply type II servo filter on measurement error and output integrator.  gain
   must be 3x1 or 3x5.  */
static inline void 
servo_typeII_filter(SERVO_T *st, dcell *merrc, double dt, const dmat *gain){
    if(!merrc) return;
    PDMAT(gain,pgain);
    int indmul=0;
    if(gain->nx!=3){
	error("Wrong format in gain\n");
    }
    double dt1=1./dt;
    double gg,ga,gs;
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
	    error("Wrong format\n");
	}
	if(indmul>0){
	    for(int imod=0; imod<nmod; imod++){
		int indm=imod * indmul;
		gg=pgain[indm][0];
		ga=pgain[indm][1];
		gs=pgain[indm][2]*dt1;
		mlead->p[imod] = (gg/(2*ga*gs+1))*(mlead->p[imod]*(2*ga*gs-1)
						   +merr->p[imod]*(2*gs+1)
						   -merrlast->p[imod]*(2*gs-1));
	    }
	}else{
	    gg=pgain[0][0];
	    ga=pgain[0][1];
	    gs=pgain[0][2]*dt1;
	    for(int imod=0; imod<nmod; imod++){
		mlead->p[imod] = (gg/(2*ga*gs+1))*(mlead->p[imod]*(2*ga*gs-1)
						   +merr->p[imod]*(2*gs+1)
						   -merrlast->p[imod]*(2*gs-1));
	    }
	}
    }
    dcellcp(&st->merrlast, merrc);
    dcelladd(&st->mpreint,1, st->mlead,1);
}
static void servo_init(SERVO_T *st, dcell *merr, const dmat *gain){
    if(!merr || st->initialized){
	error("merr must be valid and SERVO_T must be not yet initialized\n");
    }
    st->initialized=1;
    if(gain->nx>1){
	st->mpreint=dcellnew2(merr);
    }
    if(gain->nx>2){
	st->mlead=dcellnew2(merr);
	st->merrlast=dcellnew2(merr);
    }
    st->mint[0]=dcellnew2(merr); 
}
/**
   Initialize.
*/
SERVO_T *servo_new(dcell *merr, const dmat *gain){
    SERVO_T *st=calloc(1, sizeof(SERVO_T));
    st->mint=calloc(2, sizeof(dcell));
    st->nmint=2;
    if(merr && merr->nx!=0 && merr->ny!=0 && merr->p[0]){
	st->initialized=1;
	if(gain->nx>1){
	    st->mpreint=dcellnew2(merr);
	}
	if(gain->nx>2){
	    st->mlead=dcellnew2(merr);
	    st->merrlast=dcellnew2(merr);
	}
	st->mint[0]=dcellnew2(merr);
    }
    return st;
}
/**
   Applies type I or type II filter based on number of entries in gain.
*/
void servo_filter(SERVO_T *st, dcell *merr, double dt, const dmat *gain){
    if(!merr) return;
    if(!st->mint){
	error("SERVO_T must be created using servo_new()\n");
    }
    if(!st->initialized){
	servo_init(st, merr, gain);
    }
    switch(gain->nx){
    case 1://type I
	if(gain->ny!=1) error("not supported\n");
	dcelladd(&st->mpreint, 0, merr, gain->p[0]);//just record what is added.
	break;
    case 2:{//PID controller
	if(gain->ny!=1) error("not supported\n");
	double g1=gain->p[0]+gain->p[1];
	double g2=-gain->p[1];
	dcelladd(&st->mpreint, 0, merr, g1);
	dcelladd(&st->mpreint, 1, st->merrlast, g2);
	dcellcp(&st->merrlast, merr);
    }
	break;
    case 3://type II
	servo_typeII_filter(st, merr, dt, gain);
	break;
    default:
	error("Invalid");
    }
    dcelladd(st->mint, 1, st->mpreint, 1);
}
/**
   prepare the integrator by shifting commands. similar to laos.
   inte->p[0]=inte->p[0]*ap[0]+inte->p[1]*ap[1]+...
*/
void servo_shift(SERVO_T *st, dmat *ap){
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

/**
   test type I/II filter with ideal measurement to make sure it is implemented correctly.
*/
dmat* servo_test(dmat *input, double dt, int dtrat, double sigma2n, dmat *gain){
    if(input->ny==1){/*single mode. each column is for a mode.*/
	input->ny=input->nx;
	input->nx=1;
    }
    int nmod=input->nx;
    PDMAT(input,pinput);
    dmat *merr=dnew(nmod,1);
    dcell *mreal=dcellnew(1,1);
    dmat *mres=dnew(nmod,input->ny);
    double sigma=sqrt(sigma2n);
    dcell *meas=dcellnew(1,1);
    dmat *noise=dnew(nmod, 1);
    SERVO_T *st2t=servo_new(NULL, gain);
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
	    drandn(noise, sigma, &rstat);
	    dadd(&meas->p[0], 1, noise, 1);
	    servo_filter(st2t, meas, dt*dtrat, gain);
	}
    }
    dfree(merr);
    dcellfree(mreal);
    dcellfree(meas);
    servo_free(st2t);
    return mres;
}

/**
   Generate random time series using temporal PSD.
*/
dmat *psd2temp(dmat *psdin, double dt, double N, rand_t* rstat){
    double df=1./(N*dt);
    dmat *f=dlinspace(df,df,N);
    dmat *psdf=dnew_ref(psdin->nx,1,psdin->p);
    dmat *psdval=dnew_ref(psdin->nx,1,psdin->p+psdin->nx);
    dmat *psd2=dinterp1log(psdf,psdval,f);
    dfree(f);
    double sqdf=sqrt(df);
    dfree(psdf);
    dfree(psdval);
    cmat *psdc=cnew(N,1);
    cfft2plan(psdc, -1);
    for(long i=0; i<psd2->nx; i++){
	psdc->p[i]=sqrt(psd2->p[i])*sqdf*(randn(rstat)+I*randn(rstat));
    }
    cfft2(psdc,-1);
    creal2d(&psd2,0,psdc,1);
    cfree(psdc);
    /*transpose. */
    psd2->ny=psd2->nx;
    psd2->nx=1;
    return psd2;
}
/**
   Free SERVO_T struct
*/
void servo_free(SERVO_T *st){
    if(!st) return;
    dcellfree(st->mlead);
    dcellfree(st->merrlast);
    dcellfree(st->mpreint);
    dcellfreearr(st->mint, st->nmint);
    free(st);
}
/**
   Integrated a PSF that defines on logrithmically spaced grid nu.
*/
double psd_intelog(double *nu, double *psd, long n){
    double dlognu=(log(nu[n-1])-log(nu[0]))/(n-1);
    if((log(nu[1])-log(nu[0])-dlognu)>1.e-4){
	error("Nu is not logirthmic spaced\n");
    }
    double rms_sig=0;
    for(long i=0; i<n; i++){
	rms_sig+=psd[i]*nu[i];
    }	
    rms_sig*=dlognu;
    return rms_sig;
}
/**
   wraps psd_intelog
*/
double psd_intelog2(dmat *psdin){
    if(psdin->ny!=2){
	error("Wrong format\n");
    }
    double *nu=psdin->p;
    long n=psdin->nx;
    double *psd=nu+n;
    return psd_intelog(nu, psd, n);
}


/**
   Convert PSD into time series.*/
dmat* psd2time(dmat *psdin, rand_t *rstat, double dt, int nstepin){
    if(psdin->ny!=2){
	error("psd should have two columns\n");
    }
    dmat *psdx=dnew_ref(psdin->nx,1,psdin->p);
    dmat *psdy=dnew_ref(psdin->nx,1,psdin->p+psdin->nx);
    long nstep=nextpow2(nstepin);
    double df=1./(dt*nstep);
    dmat *fs=dlinspace(0, df, nstep);
    dmat *psd=NULL;
    double var=psd_intelog2(psdin);
    info2("Input psd has variance of %g m^2\n",var*1e18);
    if(fabs(psdx->p[0]*psdx->p[2]-pow(psdx->p[1],2))<EPS){/*log spaced */
	psd=dinterp1log(psdx, psdy, fs);
    }else if(fabs(psdx->p[0]+psdx->p[2]-psdx->p[1]*2.)<EPS){/*linear spaced */
	psd=dinterp1(psdx, psdy, fs);
    }else{
	error("Frequency in PSD must be log or linearly spaced\n");
    }
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
    dfree(psdx);
    dfree(psdy);
    dfree(psd);
    dfree(fs);
    dresize(out, nstepin, 1);
    double var2=dinn(out,out)/out->nx;
    info2("Time series has variance of %g m^2\n",var2*1e18);
    dscale(out, sqrt(var/var2));
    return out;
}
