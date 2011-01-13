/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include "skyc.h"
#include "types.h"
#include "servo.h"
/**
   \file servo.c
   Routines for servo optimization, filtering, etc.
*/
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
	    double logvalpre=log10(valpre);//positive. log10(val) is negative.
	    double rat=logvalpre/(logvalpre-log10(val));
	    fc=pow(10,log10(nu->p[i-1])*(1-rat)+log10(nu->p[i])*rat);
	    double phi0=atan2(cimag(Hol->p[i-1]), creal(Hol->p[i-1]));
	    double phi1=atan2(cimag(Hol->p[i]), creal(Hol->p[i]));
	    //if(fabs(phi0-phi1)>2*pi)
	    double diff=(phi1-phi0)/(2*M_PI);//positive.
	    diff=diff-round(diff);
	    phi1=phi0+diff*2*M_PI;
	    phi=pow(10,log10(phi0)*(1-rat)+log10(phi1)*rat);
	    double nphi=phi/(2*M_PI);
	    nphi=nphi-floor(nphi)-1;
	    phi=nphi*2*M_PI;//from -2*M_PI to 0.
	    found=1;
	    break;
	}
    }
    if(!found){
	error("Hol doesn't fall to below 1.\n");
    }
    *fcross=fc;
    return phi+M_PI;
}

/**
   Optimize the type II servo gains by balancing errors due to noise and
   signal propagation.

   Written: 2010-06-11

   Tested OK: 2010-06-11
*/
dmat* servo_typeII_optim(double *ress, double *resn, const dmat *psdin,
			 double fs, double lgsdt, double sigman){
    /*The upper end must be nyquist freq so that noise transfer can be
      computed. But we need to captuer the turbulence PSD beyond nyquist freq,
      which are uncorrectable.
     */
    dmat *psdx=dnew_ref(psdin->nx,1,psdin->p);
    dmat *psdy=dnew_ref(psdin->nx,1,psdin->p+psdin->nx);
     
    //Compute error in un-corretable part of the PSD
    dmat *nu2=dlogspace(log10(fs/2),3,1000);//Frequencies that no correction can be made.
    dmat *psd2=dinterp1log(psdx,psdy,nu2);
  
    double rms2_sig=psd_intelog(nu2->p, psd2->p, nu2->nx);
    dfree(nu2); dfree(psd2);

    dmat *nu=dlogspace(-3,log10(fs/2),1000);
    dmat *psd=dinterp1log(psdx,psdy,nu);
    cmat *s=cnew(nu->nx, 1);
    dcomplex pi2i=2*M_PI*I;
    double Ts=1./fs;
    cmat *Hsys=cnew(nu->nx,1);
    cmat *Hwfs=cnew(nu->nx,1);
    dcomplex expsTs,Hdac,Hmir,Hlag,Hint;
    Hmir=1;//DM
    for(long i=0; i<s->nx; i++){
	s->p[i]=pi2i*nu->p[i];
	expsTs=1-cexp(-s->p[i]*Ts);
	Hwfs->p[i]=expsTs/(Ts*s->p[i]);
	Hdac=Hwfs->p[i];
	Hlag=cexp(-s->p[i]*lgsdt);//lag
	Hint=1./expsTs;
	Hsys->p[i]=Hwfs->p[i]*Hlag*Hint*Hdac*Hmir*Hint;
    }
    double g0_min;
    if(fs<100){
	g0_min=0.01;
    }else{
	g0_min=0.01;
    }
    double g0=g0_min;
    cmat *Hol=cnew(nu->nx,1);
    double margin, fcross;
    double rms_tot_save=INFINITY, rms_sig_save=0, rms_n_save=0;
    double g0_save=0, a_save=0, T_save=0;

    while(1){
	//Open loop transfer function without phase lead
	cadd(&Hol, 0, Hsys, g0);
	margin=servo_phi_margin(&fcross, nu, Hol);
	double phineed=M_PI/4-margin;//want to ensure PI/4 margin
	double a=(1-sin(phineed))/(1+sin(phineed));
	double f0=fcross*sqrt(a);
	double T=1./(2.*M_PI*f0);
	double g=sqrt(a);

	double rms_sig=0, sum_n=0, sum_1=0;
	for(long i=0; i<s->nx; i++){
	    dcomplex Hlead=(1+T*s->p[i])/(1+a*T*s->p[i]);
	    dcomplex Holt=Hol->p[i]*Hlead*g;
	    dcomplex Hrej=1./(1.+Holt);
	    dcomplex Hcl=Holt*Hrej;
	    dcomplex Hn=Hcl/Hwfs->p[i];
	    //fixme: cabs(Hrej)*cabs(Hrej) or cabs(Href*Href)?
	    rms_sig+=psd->p[i]*creal(Hrej*conj(Hrej))*nu->p[i];
	    sum_n+=pow(cabs(Hn),2)*nu->p[i];
	    sum_1+=nu->p[i];
	}
	
	double dlognu=(log(nu->p[nu->nx-1])-log(nu->p[0]))/(nu->nx-1);
	rms_sig*=dlognu;
	sum_n*=dlognu;
	sum_1*=dlognu;
	double gain_n=(sum_n/sum_1);
	double rms_n=gain_n*sigman;
	double rms_tot=rms_sig+rms_n;
	if(rms_tot<rms_tot_save && gain_n<=1){
	    g0_save=g0*g;//g is an additional factor.
	    a_save=a;
	    T_save=T;
	    g0=g0+0.01;
	    rms_tot_save=rms_tot;
	    rms_sig_save=rms_sig;
	    rms_n_save=rms_n;
	}else{
	    break;
	}
    }

    dfree(nu);
    dfree(psdx);
    dfree(psdy);
    dfree(psd);
    cfree(Hsys);
    cfree(Hwfs);
    cfree(Hol);
    cfree(s);
    dmat *gm=dnew(3,1);
    gm->p[0]=g0_save;
    gm->p[1]=a_save;
    gm->p[2]=T_save;
    *resn=rms_n_save;
    *ress=rms_sig_save+rms2_sig;
    return gm;
}
/**
   Compute adaptive optics open loop transfer function of the type II servo with lead filter.
 */
cmat *servo_typeII_Hol(const dmat *gain, double fs, double lgsdt){
    dmat *nu=dlogspace(-3,3,1000);

    dcomplex pi2i=2*M_PI*I;
    double Ts=1./fs;
    double g0=gain->p[0];
    double a=gain->p[1];
    double T=gain->p[2];
    dcomplex s,expsTs,Hdac,Hmir,Hlag,Hint,Hsys,Hwfs,Hlead;
    Hmir=1;//DM
    cmat *Hol=cnew(nu->nx, 3);
    PCMAT(Hol, pHol);
    for(long i=0; i<nu->nx; i++){
	s=pi2i*nu->p[i];
	expsTs=1-cexp(-s*Ts);
	Hwfs=expsTs/(Ts*s);
	Hdac=Hwfs;
	Hlag=cexp(-s*lgsdt);//lag
	Hint=1./expsTs;
	Hsys=Hwfs*Hlag*Hint*Hdac*Hmir*Hint;
	Hlead=(1+T*s)/(1+a*T*s);
	pHol[0][i]=nu->p[i];
	pHol[1][i]=Hsys*g0*Hlead;
	pHol[2][i]=Hwfs;
	//Hol=Hsys*g0*Hlead;
	//Hrej=1./(1+Hol);
    }
    dfree(nu);
    return Hol;
}

/**
   Compute the residual error after servo rejection of a type II servo given the gain.
   
   Written: 2010-06-11
   
   Tested OK: 2010-06-11
*/
double servo_typeII_residual(const dmat *gain, const dmat *psdin, double fs, double lgsdt){
    dmat *nu=dlogspace(-3,3,1000);//Should go beyond Nyquist freq. Hrej=1 for nu>fs/2.
    dmat *psdx=dnew_ref(psdin->nx,1,psdin->p);
    dmat *psdy=dnew_ref(psdin->nx,1,psdin->p+psdin->nx);  
    dmat *psd=dinterp1log(psdx,psdy,nu);
    dcomplex pi2i=2*M_PI*I;
    double Ts=1./fs;

    dcomplex s,expsTs,Hdac,Hmir,Hlag,Hint,Hsys,Hwfs,Hol,Hrej,Hlead;
    Hmir=1;//DM
    double g0=gain->p[0];
    double a=gain->p[1];
    double T=gain->p[2];
    double rms_sig=0;
    for(long i=0; i<nu->nx; i++){
	s=pi2i*nu->p[i];
	expsTs=1-cexp(-s*Ts);
	Hwfs=expsTs/(Ts*s);
	Hdac=Hwfs;
	Hlag=cexp(-s*lgsdt);//lag
	Hint=1./expsTs;
	Hsys=Hwfs*Hlag*Hint*Hdac*Hmir*Hint;
	Hlead=(1+T*s)/(1+a*T*s);
	Hol=Hsys*g0*Hlead;
	Hrej=1./(1+Hol);
	rms_sig+=psd->p[i]*pow(cabs(Hrej),2)*nu->p[i];//we integrate f(nu)nu d(log(nu))
    }
    double dlognu=(log(nu->p[nu->nx-1])-log(nu->p[0]))/(nu->nx-1);
    rms_sig*=dlognu;
    dfree(psdx);
    dfree(psdy);
    dfree(nu);
    dfree(psd);
    return rms_sig;
}
/**
   Simple integrator filter.
 */
void servo_typeI_filter(SERVO_S *st, dmat *merr, double gain){
    if(!st->initialized){
	st->initialized=1;
	st->mlead=dnew(merr->nx,1);
	st->merrlast=dnew(merr->nx,1);
	st->mintfirst=dnew(merr->nx,1);
	st->mint=dnew(merr->nx,1);
    }
    dadd(&st->mint, 1, merr, gain);
}
/**
   Apply type II servo filter on measurement error and output integrator.  gain
   must be 3x1 or 3x5.  */
void servo_typeII_filter(SERVO_S *st, dmat *merr, double dtngs, const dmat *gain){

    if(!merr) return;
    PDMAT(gain,pgain);
    if(!st->initialized){
	st->initialized=1;
	st->mlead=dnew(merr->nx,1);
	st->merrlast=dnew(merr->nx,1);
	st->mintfirst=dnew(merr->nx,1);
	st->mint=dnew(merr->nx,1);
    }

    int indmul=0;
    if(gain->nx!=3){
	error("Wrong format in gain\n");
    }
    int nmod=0;//error.
    if(merr->ny==1){
	nmod=merr->nx;
    }else{
	if(merr->nx!=1){
	    error("Don't handle this case\n");
	}
	nmod=merr->ny;
    }

    if(gain->ny==nmod){
	indmul=1;
    }else if(gain->ny==1){
	indmul=0;
    }else{
	error("Wrong format\n");
    }

    double gg,ga,gs;
    for(int imod=0; imod<nmod; imod++){
	int indm=imod * indmul;
	gg=pgain[indm][0];
	ga=pgain[indm][1];
	gs=pgain[indm][2]/dtngs;
	
	st->mlead->p[imod] = (gg/(2*ga*gs+1))*(st->mlead->p[imod]*(2*ga*gs-1)
					       +merr->p[imod]*(2*gs+1)
					       -st->merrlast->p[imod]*(2*gs-1));
    }
    dcp(&st->merrlast, merr);
    dadd(&st->mintfirst,1, st->mlead,1);
    dadd(&st->mint, 1,st->mintfirst,1);
}
/**
   test type II filter with ideal measurement to make sure it is implemented correctly.
 */
dmat* servo_typeII_test(dmat *mideal, dmat *gain, double dtlgs, int dtrat){
    int nmod=mideal->nx;
    PDMAT(mideal,pmideal);
    dmat *merr=dnew(nmod,1);
    dmat *mreal=NULL;
    dmat *mres=dnew(nmod,mideal->ny);
    dmat *meas=NULL;
    SERVO_S *st2t=calloc(1, sizeof(SERVO_S));
    PDMAT(mres,pmres);
    for(int istep=0; istep<mideal->ny; istep++){
	memcpy(merr->p, pmideal[istep], nmod*sizeof(double));
	dadd(&merr, 1, mreal, -1);
	memcpy(pmres[istep],merr->p,sizeof(double)*nmod);
	if(istep % dtrat == 0){
	    dzero(meas);
	}
	dadd(&meas, 1, merr, 1);//average the error.
	dcp(&mreal, st2t->mint);
	if((istep+1) % dtrat == 0){
	    dscale(meas, 1./dtrat);
	    servo_typeII_filter(st2t, meas, dtlgs*dtrat, gain);
	    //servo_typeI_filter(st2t, meas, .5);
	}
    }
    dfree(merr);
    dfree(mreal);
    dfree(meas);
    servo_free(st2t);
    return mres;
}

/**
   Generate random time series using temporal PSD.
*/
dmat *psd2temp(dmat *psdin, double dt, double N, rand_t* rstat){
    double df=1./(N*dt);
    dmat *f=dlinspace(df,df,N);
    dmat *psdx=dnew_ref(psdin->nx,1,psdin->p);
    dmat *psdy=dnew_ref(psdin->nx,1,psdin->p+psdin->nx);
    dmat *psd2=dinterp1log(psdx,psdy,f);
    dfree(f);
    double sqdf=sqrt(df);
    dfree(psdx);
    dfree(psdy);
    cmat *psdc=cnew(N,1);
    cfft2plan(psdc, -1);
    for(long i=0; i<psd2->nx; i++){
	psdc->p[i]=sqrt(psd2->p[i])*sqdf*(randn(rstat)+I*randn(rstat));
    }
    cfft2(psdc,-1);
    creal2d(&psd2,0,psdc,1);
    cfree(psdc);
    //transpose.
    psd2->ny=psd2->nx;
    psd2->nx=1;
    return psd2;
}
/**
   Free SERVO_S struct
 */
void servo_free(SERVO_S *st){
    dfree(st->mlead);
    dfree(st->merrlast);
    dfree(st->mintfirst);
    dfree(st->mint);
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
   standalone routine that does servo filtering.
 */
#ifdef TEST
int main(int argc, char **argv){
    if(argc<6){
	info2("Usage: %s psd.bin fs dt sigma gainout.bin\n", argv[0]);
	exit(1);
    }
    dmat *psd=dread("%s",argv[1]);
    double fs=strtod(argv[2],NULL);
    double dt=strtod(argv[3],NULL);
    double sigma=strtod(argv[4],NULL);//m^2
    double res,resn;
    dmat *gain=servo_typeII_optim(&res,&resn,psd,fs,dt,sigma);
    dwrite(gain,"%s",argv[5]);
}
#endif
