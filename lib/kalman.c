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
#include "kalman.h"
#include "psd.h"
#include "nr.h"
typedef struct{
    dmat *f1;
    dmat *f2;
    dmat *psd_in;
    dmat *cov_in;
    dmat *psd_sde;
    dmat *cov_sde;
    double min;
    double max;
    double ratio;
    int ncoeff;
    int nmod;
    int ncov;
}sde_fit_t;
static void do_fft(dmat *psd){
    dfft1plan_r2hc(psd, -1);
    fft2(psd->fft, -1);
}
static void sde_psd(dmat **psd, dmat *f, double *coeff, int ncoeff, int nmod){
    if(!psd){
	*psd=dnew(f->nx,1);
    }else{
	dzero(*psd);
    }
    dmat *ppsd=*psd;
    int order=ncoeff-1;
    dcomplex twopii=2*M_PI*I;
    for(int is=1; is<f->nx; is++){
	dcomplex s=f->p[is]*twopii;
	for(int im=0; im<nmod; im++){
	    dcomplex denom=cpow(s, order);
	    for(int ic=0; ic<order-1; ic++){
		denom+=coeff[ic+im*ncoeff]*cpow(s, order-ic-1);
	    }
	    denom+=coeff[order-1+im*ncoeff];
	    double val=pow(coeff[order+im*ncoeff],2)/cabs2(denom);
	    if(ppsd->ny==nmod){//seperate
		ppsd->p[is+im*f->nx]=val;
	    }else{
		ppsd->p[is]+=val;
	    }
	}	
    }
}
static int coeff_isgood(double *coeff, int ncoeff, double min, double max){
    for(int icoeff=0; icoeff<ncoeff; icoeff++){
	if(coeff[icoeff]<min || coeff[icoeff]>max){
	    //warning("coeff[%d]=%g, beyond [%g %g]\n", icoeff, coeff[icoeff], min, max);
	    return coeff[icoeff];
	}
    }
    return 1;
}
/**
   receives the difference in psd
*/
static double sde_diff_psd(double *coeff, void *pdata){
    sde_fit_t *data=pdata;
    if(!coeff_isgood(coeff, data->ncoeff*data->nmod, data->min, data->max)){
	return INFINITY;
    }
    sde_psd(&data->psd_sde, data->f1, coeff, data->ncoeff, data->nmod);
    dadd(&data->psd_sde, 1, data->psd_in, -1);
    double diff=dnorm2(data->psd_sde);
    return diff;
}
/**
   receives the difference in cov
*/
static double sde_diff_cov(double *coeff, void *pdata){
    sde_fit_t *data=pdata;
    double check;
    if(!(check=coeff_isgood(coeff, data->ncoeff*data->nmod, data->min, data->max))){
	return (1+fabs(check))*1e50;
    }
    sde_psd(&data->cov_sde, data->f2, coeff, data->ncoeff, data->nmod);
    do_fft(data->cov_sde);
    dmat *cov1=dnew_ref(data->ncov, 1, data->cov_sde->p);
    dmat *cov2=dnew_ref(data->ncov, 1, data->cov_in->p);
    data->ratio=cov2->p[0]/cov1->p[0];
    /*scale to same max value*/
    dadd(&cov1, 1./cov1->p[0], cov2, -1./cov2->p[0]);
    double diff=dnorm2(cov1);
    dfree(cov1); dfree(cov2);
    return diff;
}
/**
   Fit a PSD with SDE model. Use initial guess from coeff0 and return final
   answer in the same format.
*/
dmat* sde_fit(const dmat *psdin, const dmat *coeff0, double tmax_fit, double min, double max){
    if(psdin->ny!=2){
	error("psd must contain nu and psd\n");
    }
    int nf=1024*4; /**array length for FFT*/
    info("sdefit: nf=%d\n", nf);
    double maxf=psdin->p[psdin->nx-1];
    double df=maxf/(nf/2-1);
    double dt=1./maxf;
    dmat *f2=dnew(nf, 1);
    for(long i=0; i<nf/2; i++){
	f2->p[i]=df*i;
	f2->p[i+nf/2]=df*(nf/2-i);
    }
    dmat *cov=psdinterp1(psdin, f2);
    dmat *f1=dnew_ref(nf/2, 1, f2->p);
    dmat *psd=psdinterp1(psdin, f1);
    int ncoeff=coeff0->nx;
    int nmod=coeff0->ny;
    dmat *coeff=dnew(ncoeff, nmod);dcp(&coeff, coeff0);
    dmat *scale=dnew(ncoeff, nmod);dadd(&scale, 0, coeff0, 0.01);
    dmat *psd_sde=dnew(nf/2, 1);
    dmat *cov_sde=dnew(nf, 1);
    do_fft(cov);
    int ncov=round(tmax_fit/dt);
    if(ncov>nf/2){
	ncov=nf/2;
    }
    sde_fit_t data={f1, f2, psd, cov, psd_sde, cov_sde, min, max, 0, ncoeff, nmod, ncov};
    if(data.ncov>0){//covariance fitting
	double tol=1e-15;
	dminsearch(coeff->p, scale->p, ncoeff*nmod, tol, sde_diff_cov, &data);
	for(int im=0; im<nmod; im++){
	    coeff->p[(1+im)*ncoeff-1]*=sqrt(data.ratio);
	}
    }else{//PSD fitting
	double tol=dnorm2(psd)*1e-15;
	dminsearch(coeff->p, scale->p, ncoeff*nmod, tol, sde_diff_psd, &data);
    }
    dfree(f1);
    dfree(f2);
    dfree(cov);
    dfree(psd);
    dfree(psd_sde);
    dfree(cov_sde);
    dfree(scale);
    return coeff;
}
static void reccati(dmat **Mout, dmat **Pout, 
		    dmat *A, dmat *Qn, dmat *C, dmat *Rn){
    double diff=1, diff2=1, lastdiff=INFINITY;
    dmat *P=dnew(A->nx, A->ny); daddI(P, 1);//Initialize P to identity
    dmat *AP=0, *CP=0, *P2=0, *CPAt=0, *CPCt=0, *tmp=0;
    while(diff>1e-13 && diff2>1e-12){
	dmm(&AP, 0, A, P, "nn", 1);
	dmm(&CP, 0, C, P, "nn", 1);
	dmm(&CPAt, 0, CP, A, "nt", 1);
	dmm(&CPCt, 0, CP, C, "nt", 1);
	dadd(&CPCt, 1, Rn, 1);
	dsvd_pow(CPCt, -1, 0);
	dcp(&P2, Qn);
	dmm(&P2, 1, AP, A, "nt", 1);
	dmm(&tmp, 0, CPCt, CPAt, "nn", 1);
	dmm(&P2, 1, CPAt, tmp, "tn", -1);
	diff=dnorm2(P);
	dadd(&P, 1, P2, -1);
	diff=sqrt(dnorm2(P)/diff);
	diff2=fabs(diff-lastdiff);
	lastdiff=diff;
	dcp(&P, P2);
    }
    dmm(Mout, 0, CP, CPCt, "tn", 1);
    if(Pout) dcp(Pout, P);
    dfree(AP);dfree(CP); dfree(P); dfree(P2); 
    dfree(CPAt); dfree(CPCt); dfree(tmp);
}
/*Kalman filter based on SDE model*/
kalman_t* sde_kalman(dmat *coeff, /**<SDE coefficients*/
		     double dthi, /**<Loop frequency*/
		     int dtrat,   /**<WFS frequency as a fraction of loop*/
		     dmat *Gwfs,  /**<WFS measurement from modes. Can be identity*/
		     dmat *Rwfs,  /**<WFS measurement noise covariance*/
		     dmat *Proj   /**<Project modes in statespace to DM/correction space*/){
    int nblock=coeff->ny;
    int order=coeff->nx-1;
    /*Ac is block diagonal matrix for continuous domain state evolution. 
      For third order, each block is 
      |-c1 -c2 -c3 |
      |1    0  0  |
      |0    1  0  |
    */

    int nmod=nblock*order;/*number of modes in state*/
    dmat *Ac=dnew(nmod, nmod);
    dmat *Sigma_ep=dnew(nmod, nmod);
    dmat *Pd0=dnew(nblock, nmod);
    for(int iblock=0; iblock<nblock; iblock++){
	double *pcoeff=coeff->p+(order+1)*iblock;
	double *Aci=Ac->p+iblock*order*(nmod+1);
	double *si=Sigma_ep->p+iblock*order*(nmod+1);
	for(int i=0; i<order; i++){
	    Aci[i*nmod]=-pcoeff[i]; /*evolution*/
	    if(i+1<order){
		Aci[i*(nmod+1)+1]=1;/*keep*/
	    }
	}
	si[0]=pow(pcoeff[order],2);
	Pd0->p[((iblock+1)*order-1)*nblock+iblock]=1;
    }
    dmat *Pd=0; 
    if(Proj){
	dmm(&Pd, 0, Proj, Pd0, "nn", 1); 
    }else{
	Pd=dref(Pd0);
    }
    dfree(Pd0);
    dmat* Sigma_varep=0;/*state noise*/
    dmat* Sigma_zeta=0;/*state measurement noise*/
    dmat *AcI=ddup(Ac); 
    dsvd_pow(AcI, -1, 0);
    double dT=dthi*dtrat;
    {
	/*Evaluate noise covariance matrix of discrete state, and measurement of the
	 * state*/
	int nsec=dtrat*10;
	double dT2=dT/nsec;
	dmat *expAj=0;
	dmat *expAjn=0;
	dmat *tmp1=0, *tmp2=0;
	for(int i=0; i<nsec; i++){
	    double ti=(i+0.5)*dT2;
	    dexpm(&expAj, 0, Ac, ti);
	    dexpm(&expAjn, 0, Ac, -ti);
	 
	    dmm(&tmp1, 0, expAj, Sigma_ep, "nn", 1);
	    dmm(&Sigma_varep, 1, tmp1, expAj, "nt", 1);
	    
	    dscale(expAjn, -1); daddI(expAjn, 1); //1-exp(-A*tj)
	    dmm(&tmp1, 0, expAjn, AcI, "nn", 1);
	    dmm(&tmp2, 0, tmp1, Sigma_ep, "nn", 1);
	    dmm(&Sigma_zeta, 1, tmp2, tmp1, "nt", 1);
	}
	dscale(Sigma_varep, dT2);
	dscale(Sigma_zeta, dT2/(dT*dT));
	dfree(tmp1); dfree(tmp2);
	dfree(expAj); dfree(expAjn);
    }
    dmat *Ad=0; dexpm(&Ad, 0, Ac, dT);  /*discrete state propagation at dT*/
    dmat *AdM=0; dexpm(&AdM, 0, Ac, dthi);/*discrete state propagation at dthi*/
    dmat *Radd=0;/*measurement error*/
    {
	dmat *tmp=0;
	dmm(&tmp, 0, Sigma_zeta, Pd, "nt", 1);
	dmm(&Radd, 0, Pd, tmp, "nn", 1);
	dfree(tmp);
    }
    dmat *Cd=0; /*From discrete state to WFS measurement*/
    dmat *Fd=0; /*From discrete state to averaged mode dT*/
    dmat *FdM=0;/*From discrete state to averaged mode for dthi*/
    {
	dmat *Xi=0;
	dmat *XiM=0;
	dmat *tmp=0;
	dexpm(&tmp, 0, Ac, -dT);
	dscale(tmp, -1); 
	daddI(tmp, 1);
	dmm(&Xi, 0, tmp, AcI, "nn", 1./dT);
	
	dexpm(&tmp, 0, Ac, -dthi);
	dscale(tmp, -1); 
	daddI(tmp, 1);
	dmm(&XiM, 0, tmp, AcI, "nn", 1./dthi);
	
	dmm(&Fd, 0, Pd, Xi, "nn", 1);
	dmm(&FdM, 0, Pd, XiM, "nn", 1);
	dmm(&Cd, 0, Gwfs, Fd, "nn", 1);
	dfree(Xi);
	dfree(XiM);
	dfree(tmp);
    }

    dmat *Rn=ddup(Rwfs);
    {
	/*compute Rn=Rn+Gwfs*Radd*Gwfs' */
	dmat *tmp=0;
	dmm(&tmp, 0, Radd, Gwfs, "nt", 1);
	dmm(&Rn, 1, Gwfs, Radd, "nn", 1);
	dfree(tmp);
    }
    dmat *P=0, *M=0;
    reccati(&M, &P, Ad, Sigma_varep, Cd, Rn);
    {
	/*computed estimation error in modes */
	dmat *tmp1=0;
	dmm(&tmp1, 0, Pd, P, "nn", 1);
	dfree(P);
	dmm(&P, 0, tmp1, Pd, "nt", 1);
	dfree(tmp1);
    }
    kalman_t *res=calloc(1, sizeof(kalman_t));
    res->Ad=Ad;
    res->Fd=Fd;
    res->Cd=Cd;
    res->AdM=AdM;
    res->FdM=FdM;
    res->M=M;
    res->P=P;
    res->dthi=dthi;
    res->dtrat=dtrat;
    res->Gwfs=ddup(Gwfs);
    res->Rwfs=ddup(Rwfs);
    dfree(Ac);
    dfree(Sigma_ep);
    dfree(Pd);
    dfree(Sigma_varep);
    dfree(Sigma_zeta);
    dfree(AcI);
    dfree(Radd);
    dfree(Rn);
    return res;
}
void kalman_free(kalman_t *kalman){
    if(!kalman) return;
    dfree(kalman->Ad);
    dfree(kalman->Fd);
    dfree(kalman->Cd);
    dfree(kalman->AdM);
    dfree(kalman->FdM);
    dfree(kalman->M);
    dfree(kalman->P);
    dfree(kalman->Gwfs);
    dfree(kalman->Rwfs);
    dfree(kalman->xhat);
    dfree(kalman->xhat2);
    dfree(kalman->xhat3);
}
/*Initialize kalman filter state*/
void kalman_init(kalman_t *kalman){
    if(!kalman->xhat){
	kalman->xhat=dnew(kalman->Ad->nx, 1);
	kalman->xhat2=dnew(kalman->Ad->nx, 1);
	kalman->xhat3=dnew(kalman->Ad->nx, 1);
    }else{
	dzero(kalman->xhat);
	dzero(kalman->xhat2);
	dzero(kalman->xhat3);
    }
}
/*Update state vector when there is a measurement. It modifies meas*/
void kalman_update(kalman_t *kalman, dmat *meas){
    /*differene between real and predicted measurement*/
    dmm(&meas, 1, kalman->Cd, kalman->xhat, "nn", -1);
    /*updated estimate*/
    dmm(&kalman->xhat, 1, kalman->M, meas, "nn", 1);
    dmm(&kalman->xhat2,0, kalman->AdM, kalman->xhat, "nn", 1);/*correction for this time step*/
    dcp(&kalman->xhat3, kalman->xhat);
    dmm(&kalman->xhat, 0, kalman->Ad, kalman->xhat3, "nn", 1);/*predicted state*/
}
/*Output correction*/
void kalman_output(kalman_t *kalman, dmat **out, double alpha, double beta){
    dcp(&kalman->xhat3, kalman->xhat2);
    dmm(&kalman->xhat2, 0, kalman->AdM, kalman->xhat3, "nn", 1);
    dmm(out, alpha, kalman->FdM, kalman->xhat2, "nn", beta);
}
/**
   Test the performance of kalman filter (LQG controller). Derived from servo_test()
*/
dmat *kalman_test(kalman_t *kalman, dmat *input){
    if(input->ny==1){/*single mode. each column is for a mode.*/
	input->ny=input->nx;
	input->nx=1;
    }
    int dtrat=kalman->dtrat;
    const dmat *Gwfs=kalman->Gwfs;
    dmat *rmsn=0;
    if(dnorm2(kalman->Rwfs)>0){
	rmsn=dchol(kalman->Rwfs);
    }
    double dtrat1=1./dtrat;
    int nmod=input->nx;
    PDMAT(input,pinput);
    dmat *mres=ddup(input);
    PDMAT(mres, pmres);
    dmat *meas=0;
    dmat *noise=dnew(Gwfs->nx, 1);
    kalman_init(kalman);
    rand_t rstat;
    seed_rand(&rstat, 1);
    dmat *ini=dnew_ref(Gwfs->nx,1,(double*)1);//wrapper
    dmat *outi=dnew_ref(nmod, 1, (double*)1);//wraper
    for(int istep=0; istep<input->ny; istep++){
	ini->p=(double*)(pinput+istep);
	outi->p=(double*)(pmres+istep);
	kalman_output(kalman, &outi, 1, -1);
	if((istep) % dtrat == 0){/*There is measurement from last step*/
	    if(rmsn){
		drandn(noise, 1, &rstat);
		if(rmsn->nx>1){
		    dmm(&meas, 1, rmsn, noise, "nn", 1);/*add measurement error*/
		}else{
		    dadd(&meas, 1, noise, rmsn->p[0]);
		}
	    }
	    kalman_update(kalman, meas);
	    dzero(meas);
	}
	dmm(&meas, 1, Gwfs, ini, "nn", dtrat1);/*OL measurement*/
    }
    dfree(rmsn);
    dfree(meas);
    dfree(noise);
    dfree(ini);
    dfree(outi);
    return mres;
}
