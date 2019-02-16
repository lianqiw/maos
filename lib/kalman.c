/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
int KALMAN_IN_SKYC=0;
typedef struct{
    double df;
    dmat *freq;
    dmat *psdcov_in;
    dmat *psdcov_sde;
    double norm_in;
    int count;
    int ncoeff;
    int nmod;
    int ncov;
}sde_fit_t;
/**
   Do FFT to convert PSD to Covariance
 */
static void psd2cov(dmat *psd, double df){
    dscale(psd, df);
    psd->p[0]=0;//zero dc. critical
    dfft1plan_r2hc(psd, -1);
    dfft1(psd, -1);
}
/**
   Computes PSD of SDE with coefficients (coeff) and frequency (f)
 */
void sde_psd(dmat **psd, const dmat *f, const double *coeff, int ncoeff, int nmod){
    int order=ncoeff-1;
    if(!*psd){
	*psd=dnew(f->nx,1);
    }else{
	dzero(*psd);
    }
    dmat *ppsd=*psd;
    if(ppsd->nx != f->nx){
	warning("f and psd does not match\n");
    }
    if(order==2){//special version for order 2. faster.
	double twopi=2*M_PI;
	for(int im=0; im<nmod; im++){
	    double c1sq=pow(coeff[0+im*3],2);
	    double c2=coeff[1+im*3];
	    double sigma2=pow(coeff[2+im*3],2);
	    double *pp;
	    if(ppsd->ny==nmod){//seperate
		pp=ppsd->p+im*ppsd->nx;
	    }else{
		pp=ppsd->p;
	    }
	    for(int is=0; is<f->nx; is++){
		double omega2=pow(f->p[is]*twopi, 2);
		pp[is]+=sigma2/(pow(omega2-c2, 2)+c1sq*omega2);//2017-07-21: was =
	    }
	}
    }else{
	dcomplex twopii=COMPLEX(0, 2*M_PI);
	for(int is=0; is<f->nx; is++){
	    dcomplex s=f->p[is]*twopii;
	    for(int im=0; im<nmod; im++){
		dcomplex denom=cpow(s, order);
		for(int ic=0; ic<order-1; ic++){
		    denom+=coeff[ic+im*ncoeff]*cpow(s, order-ic-1);
		}
		denom+=coeff[order-1+im*ncoeff];
		double val=pow(coeff[order+im*ncoeff],2)/cabs2(denom);
		if(ppsd->ny==nmod){//seperate
		    P(ppsd, is, im)+=val;
		}else{
		    ppsd->p[is]+=val;
		}
	    }	
	}
    }
}
static double coeff_isbad(const double *coeff, int ncoeff){
    for(int icoeff=0; icoeff<ncoeff; icoeff++){
	if(coeff[icoeff]<=0){
	    return 1-coeff[icoeff];
	}
    }
    return 0;
}
/**

   Compute the difference between the input PSD and SDE PSD, or the covariance of ncov>0.

   2007-07-20: Updated the routine for robustness
   
   We no longer require the total variance to be preserved to be more flexible during fitting.
   
   It is better to fit covariance than PSD because what the LQG
   cares is how the variable involves with time within the next time
   steps. Covariance describes this better than the PSD.
*/

static double sde_diff(const double *coeff, void *pdata){
    sde_fit_t *data=(sde_fit_t*)pdata;
    double tmp;
    double diff;
    data->count++;
    if((tmp=coeff_isbad(coeff, data->ncoeff*data->nmod))){
	diff=(tmp)*10+10;
    }else{
	sde_psd(&data->psdcov_sde, data->freq, coeff, data->ncoeff, data->nmod);
	if(!data->ncov){//PSD fitting
	    dadd(&data->psdcov_sde, 1, data->psdcov_in, -1);
	    /*Must use dnorm here. dsumabs does not work well for peaks when
	     * initial guess of c1 is not close.*/
	    diff=dnorm(data->psdcov_sde)/data->norm_in;
	}else{//Covariance fitting.
	    /*if(data->count<20 || data->count%100==0){
		writebin(data->psdcov_sde, "sde_psd_%04d", data->count);
		}*/
	    psd2cov(data->psdcov_sde, data->df);
	    /*if(data->count<20 || data->count%100==0){
		writebin(data->psdcov_sde, "sde_cov_%04d", data->count);
		}*/
	    diff=0;
#define METHOD 1
#if METHOD==1
	    for(long i=0; i<data->ncov; i++){
		double val1=data->psdcov_in->p[i];
		double val2=data->psdcov_sde->p[i];
		diff+=pow(val1-val2,2);
	    }
	    diff=diff/((abs2(data->psdcov_in->p[0]+data->psdcov_sde->p[0]))*data->ncov);
#elif METHOD==2
	    double scale=data->psdcov_in->p[0]/data->psdcov_sde->p[0];
	     for(long i=0; i<data->ncov; i++){
		double val1=data->psdcov_in->p[i];
		double val2=data->psdcov_sde->p[i];
		diff+=pow(val1-val2*scale,2);
	    }
	     diff+=abs2(data->psdcov_in->p[0]-data->psdcov_sde->p[0])*data->ncov;
	     diff=diff/((abs2(data->psdcov_in->p[0]))*data->ncov);
#endif
	}
    }
    if(!isfinite(diff)){
	error("Diff is not finite\n");
    }
    /*info("coeff=");
    for(int imod=0; imod<data->nmod; imod++){
	info("[%g %.2f %g] ", coeff[3*imod], coeff[3*imod+1], coeff[3*imod+2]);
    }
    info("ncov=%d. cov0=[%g, %g] diff=%g. count=%d\n", data->ncov, data->psdcov_in->p[0], data->psdcov_sde->p[0], diff, data->count);
    */
    return diff;
}
/**
   Scale the third coefficient (sigma) to have the same energy
*/
static void sde_scale_coeff(dmat *coeff, double var_in, dmat *psdcov_sde, const dmat *freq, int ncov){
    for(int iy=0; iy<coeff->ny; iy++){
	if(coeff->p[2+iy*3]<=0){
	    coeff->p[2+iy*3]=1;
	}
    }
    double ratio;
    sde_psd(&psdcov_sde, freq, coeff->p, coeff->nx, coeff->ny);
    double var_sde;
    if(ncov){
	//writebin(psdcov_sde, "sde_psd");
	psd2cov(psdcov_sde, freq->p[2]-freq->p[1]);
	//writebin(psdcov_sde, "sde_cov");
	var_sde=psdcov_sde->p[0];
    }else{
	var_sde=dtrapz(freq, psdcov_sde);
    }
    ratio=sqrt(var_in/var_sde);
    for(int iy=0; iy<coeff->ny; iy++){
	coeff->p[2+iy*3]*=ratio;
    }
}
/**
   Fit a PSD with SDE model. Use initial guess from coeff0 and return final
   answer in the same format.
*/
static dmat* sde_fit_do(const dmat *psdin, const dmat *coeff0, double tmax_fit){
    if(psdin->ny!=2){
	error("psd must contain nu and psd\n");
    }
    double df=0.01;//ensure we can sample narrow peaks when interpolating psd
    dmat *freq=0, *psdcov_in=0, *psdcov_sde=0;
    int ncov=0;
    double var_in=0;//strength of input
    /* 2014-07-31: Always upsample the PSD to 0.01 Hz sampling. Otherwise SDE
       PSD may have unresolved peaks, causing fitting to misbehave.
    */
    double maxf=psdin->p[psdin->nx-1];
    double norm_in;
    if(tmax_fit>0){//Use covariance fit
	/* Create Frequency vector with proper fft indexing [0:1:nf/2-1 * nf/2:-1:1]. 
	   Negative index is mapped to positive. */
	long nf=round(maxf/df)*2;
	ncov=MIN(round(tmax_fit*maxf), nf/2);
	freq=dnew(nf, 1);
	for(long i=0; i<nf/2; i++){
	    freq->p[i]=df*i;
	    freq->p[i+nf/2]=df*(nf/2-i);
	}
	psdcov_in=psdinterp1(psdin, freq, 0);
	
	//writebin(psdcov_in, "in_psd");
	psd2cov(psdcov_in, df);
	var_in=psdcov_in->p[0];
	if(!isfinite(psdcov_in->p[0])){
	    warning("covariance is not finite at 0\n");
	    writebin(psdcov_in, "bad_cov");
	}
	//writebin(psdcov_in, "in_cov");
	psdcov_sde=dnew(nf, 1);
	norm_in=ncov*abs2(psdcov_in->p[0]);
    }else{//Use PSD fit
	if((psdin->p[1]-psdin->p[0])<df*2){
	    freq=drefcols(psdin, 0, 1);
	    psdcov_in=drefcols(psdin, 1, 1);
	}else{//Resampling.
	    double minf=psdin->p[0];
	    long nf=round((maxf-minf)/df);
	    freq=dnew(nf, 1);
	    for(long i=0; i<nf; i++){
		freq->p[i]=minf+df*i;
	    }
	    psdcov_in=psdinterp1(psdin, freq, 0);
	}
	psdcov_sde=dnew(freq->nx, 1);
	var_in=dtrapz(freq, psdcov_in);
	norm_in=dnorm(psdcov_in);
    }
    int ncoeff=coeff0->nx;
    int nmod=coeff0->ny;
    dmat *coeff=dnew(ncoeff, nmod);dcp(&coeff, coeff0);
    //Scale to make sure total energy is preserved.
    sde_scale_coeff(coeff, var_in, psdcov_sde, freq, ncov);

    sde_fit_t data={df, freq, psdcov_in, psdcov_sde, norm_in, 0, ncoeff, nmod, ncov};
    double diff0=sde_diff(coeff->p, &data);
    double tol=1e-10;
    int nmax=2000;
    dminsearch(coeff->p, ncoeff*nmod, tol, nmax, sde_diff, &data);
    //Do not scale coeff after the solution.
    double diff1=sde_diff(coeff->p, &data);
    info("sde_fit: %d interations: %g->%g.\n", data.count, diff0, diff1);
    //Scale to make sure total energy is preserved.
    /*
      if(diff1>0.2 && diff1>diff0*0.75){
	static int count=0;
	writebin(psdin, "sde_fit_psdin_%d_%g", count, tmax_fit);
	writebin(coeff0,"sde_fit_coeff_%d_%g", count, tmax_fit);
	count++;
	if(tmax_fit>0){
	    info("Redo with PSD fitting.\n");
	    dfree(coeff);
	    coeff=sde_fit_do(psdin, coeff0, 0);
	}else{
	    warning("Failed to converge.\n");
	}
    }else{
	info("\n");
	}*/
    dfree(freq);
    dfree(psdcov_in);
    dfree(psdcov_sde);
    return coeff;
}

/**
   Estiamte the total PSD power for vibration peaks using FWHM*peak
 */
static double sde_vib_est(double c1, double c2){
    double sqrt4ac=sqrt(c1*c1+c2*4);
    double omega1=sqrt(0.5*((c2*2-c1*c1)-sqrt4ac));
    double omega2=sqrt(0.5*((c2*2-c1*c1)+sqrt4ac));
    double fwhm=(omega2-omega1)/(2*M_PI);
    double peak=1./(c1*c1*c2);
    return fwhm*peak;
}
/**
   If coeff0 is not null, use it immediately, otherwise, do vibration identification
 */
dmat* sde_fit(const dmat *psdin, const dmat *coeff0, double tmax_fit, int vibid){
    if (coeff0){
	return sde_fit_do(psdin, coeff0, tmax_fit);
    }else{
	//Do vibration identification
	dmat *vibs=vibid?psd_vibid(psdin):NULL;
	dcell *coeffs=dcellnew(1, vibs?(1+vibs->ny):1);
	dmat *coeffi=dnew(3,1);
	dmat *psd2=ddup(psdin);
	if(vibs && vibs->ny>0){
	    dbg("\nnvib=%ld\n", vibs->ny);
	    for(int ivib=0; ivib<vibs->ny; ivib++){
		double fi=vibs->p[ivib*4+0];
		int i1=vibs->p[ivib*4+2];
		int i2=vibs->p[ivib*4+3]+1;
		if(i2-i1<3){
		    i2++;
		    i1--;
		}
		dmat *psdi=dsub(psdin, i1, i2-i1, 0, 0);//extract peak
		coeffi->p[0]=1;
		coeffi->p[1]=pow(2*M_PI*fi,2);
		coeffi->p[2]=0;//sde_fit_do will figure it out.
		coeffs->p[ivib]=sde_fit_do(psdi, coeffi, tmax_fit);
		if(fabs(sqrt(coeffi->p[1])-sqrt(coeffs->p[ivib]->p[1]))>2*M_PI){
		    warning("Fitting failed for %d, Freq=%g, %g\n", ivib, 
			    sqrt(coeffi->p[1])/(2*M_PI), sqrt(coeffs->p[ivib]->p[1])/(2*M_PI));
		    writebin(psdi, "psdi_%d", ivib);
		    writebin(coeffi, "coeffi_%d", ivib);
		    writebin(coeffs->p[ivib], "coeffo_%d", ivib);
		    if(fabs(sqrt(coeffi->p[1])-sqrt(coeffs->p[ivib]->p[1]))>2*M_PI*5){
			dcp(&coeffs->p[ivib], coeffi);
			warning("Use input\n");
		    }
		}
		double* psd2p=psd2->p+psd2->nx;
		for(int i=i1; i<i2; i++){//replace the peak by a linear interpolation
		    psd2p[i]=(psd2p[i1-1]*(i2+1-i)+psd2p[i2+1]*(i-(i1-1)))/(i2-i1+2);
		}
	    }
	    dfree(vibs);
	}
	coeffi->p[0]=1;
	coeffi->p[1]=1;
	coeffi->p[2]=0;//sde_fit_do will figure it out.
	coeffs->p[coeffs->ny-1]=sde_fit_do(psd2, coeffi, tmax_fit);
	dfree(coeffi);
	dfree(psd2);
	dmat *coeff=dcell2m(coeffs);
	dcellfree(coeffs);
	return coeff;
    }
}
/**
   Compute the reccati equation.
 */
/*Notice that P may not be symmetric due to round off errors.*/
#define RECCATI_CALC				\
    dcp(&P, P2);				\
    dmm(&AP, 0, A, P, "nn", 1);			\
    dmm(&CP, 0, C, P, "nn", 1);			\
    dmm(&CPAt, 0, CP, A, "nt", 1);		\
    dmm(&APCt, 0, AP, C, "nt", 1);		\
    dmm(&CPCt, 0, CP, C, "nt", 1);		\
    dadd(&CPCt, 1, Rn, 1);			\
    dsvd_pow(CPCt, -1, 1e-12);			\
    dmm(&P2, 0, AP, A, "nt", 1);		\
    dmm(&tmp, 0, CPCt, CPAt, "nn", 1);		\
    dmm(&P2, 1, APCt, tmp, "nn", -1);		\
    dadd(&P2, 1, Qn, 1)

dmat* reccati(dmat **Pout, const dmat *A, const dmat *Qn, const dmat *C, const dmat *Rn){
    double diff=1, diff2=1, lastdiff=INFINITY;
    dmat *P2=dnew(A->nx, A->ny); daddI(P2, Qn->p[0]);//Initialize P to identity
    dmat *AP=0, *CP=0, *P=0, *CPAt=0, *CPCt=0, *APCt=0, *tmp=0;
    int count=0;
    double thres=1e-14;
    const int maxcount=10000;
    while(diff>thres && diff2>thres && count++<maxcount){
	RECCATI_CALC;//P2 has the new result.
	diff=dsumsq(P);
	dadd(&P, 1, P2, -1);
	diff=sqrt(dsumsq(P)/diff);
	diff2=fabs(diff-lastdiff);
	lastdiff=diff;
    }
    if(count>=maxcount){
	warning_once("recatti: count=%d, diff=%g, diff2=%g, thres=%g\n", count, diff, diff2, thres);
    }
    dmat *Mout=0;
    dmm(&Mout, 0, CP, CPCt, "tn", 1);
    if(Pout) dcp(Pout, P2);
    dfree(AP);dfree(CP); dfree(P); dfree(P2); 
    dfree(CPAt); dfree(CPCt); dfree(APCt); dfree(tmp);
    return Mout;
}
#if 1
dcell* reccati_cell(dmat **Pout, const dmat *A, const dmat *Qn, const dcell *Cs, const dcell *Rns){
    int nk=Cs->nx;
    dcell *Mout=dcellnew(nk, 1);
    if(nk==1){
	Mout->p[0]=reccati(Pout&&!*Pout?Pout:0, A, Qn, Cs->p[0], Rns->p[0]);
	return Mout;
    }
    double diff=1, diff2=1, lastdiff=INFINITY;
    dmat *P2=dnew(A->nx, A->ny); daddI(P2, Qn->p[0]);//Initialize P to identity
    dmat *AP=0, *CP=0, *P=0, *CPAt=0, *CPCt=0, *APCt=0, *tmp=0;
    int count=0;
    int ik=-1;
    const double thres=1e-14;
    const int maxcount=10000;
    while((diff>thres && diff2>thres) && count++<maxcount){
	do{
	    ik=(ik+1)%nk;
	}while(!Cs->p[ik]);
	dmat *C=Cs->p[ik];
	dmat *Rn=Rns->p[ik];
	RECCATI_CALC;
	diff=dsumsq(P);
	dadd(&P, 1, P2, -1);
	diff=sqrt(dsumsq(P)/diff);
	diff2=fabs(diff-lastdiff);
	lastdiff=diff;
    }
    if(count>=maxcount){
	warning_once("reccati: count=%d, diff=%g, diff2=%g, thres=%g\n", count, diff, diff2, thres);
    }
   
    for(ik=0; ik<nk; ik++){
	dmat *C=Cs->p[ik];
	dmat *Rn=Rns->p[ik];
	if(!C) continue;
	RECCATI_CALC;
	dmm(&Mout->p[ik], 0, CP, CPCt, "tn", 1);
    }
    if(Pout) dcp(Pout, P);
    dfree(AP);dfree(CP); dfree(P); dfree(P2); 
    dfree(CPAt); dfree(CPCt); dfree(APCt); dfree(tmp);
    return Mout;
}
#else //the above one is slightly better.
dcell* reccati_cell(dmat **Pout, const dmat *A, const dmat *Qn, const dcell *Cs, const dcell *Rns){
    dcell *Mout=cellnew(Cs->nx, 1);
    for(int ik=0; ik<Cs->nx; ik++){
	if(Cs->p[ik]){
	    Mout->p[ik]=reccati(Pout&&!*Pout?Pout:0, A, Qn, Cs->p[ik], Rns->p[ik]);
	}
    }
    return Mout;
}
#endif
/**
   Kalman filter based on SDE model
*/
kalman_t* sde_kalman(const dmat *coeff, /**<SDE coefficients*/
		     const double dthi, /**<Loop frequency*/
		     const lmat *dtrat_wfs,   /**<WFS frequency as a fraction of loop*/
		     const dcell *Gwfs,  /**<WFS measurement from modes. Can be identity*/
		     const dcell *Rwfs,  /**<WFS measurement noise covariance*/
		     const dmat *Proj   /**<Project modes in statespace to DM/correction space*/){
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
    dmat *Pd0=dnew(nblock, nmod);//Project from state vector to state.
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
    dmat *AcI=ddup(Ac); 
    dsvd_pow(AcI, -1, 0);
    const int nwfs=Gwfs->nx;
    if(dtrat_wfs->nx!=nwfs || dtrat_wfs->ny!=1){
	error("dtrat_wfs should have size %ldx1\n", dtrat_wfs->nx);
    }
    int ndtrat=0;
    lmat *dtrats=lnew(nwfs, 1);//unique dtrats
    int dtrat_prod=1;
    ndtrat=0;
    for(int iwfs=0; iwfs<dtrat_wfs->nx; iwfs++){
	int found=0;
	for(int jdtrat=0; jdtrat<ndtrat; jdtrat++){
	    if(dtrats->p[jdtrat]==dtrat_wfs->p[iwfs]){
		found=1; break;
	    }
	}
	if(!found){
	    dtrats->p[ndtrat++]=dtrat_wfs->p[iwfs];
	    dtrat_prod*=dtrat_wfs->p[iwfs];
	}
    }
    lresize(dtrats, ndtrat, 1);
    lsort(dtrats, 1);
    const int dtrat_min=dtrats->p[0];
    dmat*  Sigma_varep=0;/*state noise in discrete space*/
    dcell* Sigma_zeta =dcellnew(ndtrat,1);/*state measurement noise*/
    dcell *Radd       =dcellnew(ndtrat,1);
    dcell *Raddchol   =dcellnew(ndtrat,1);
    dcell *Xi         =dcellnew(ndtrat,1);
    for(int idtrat=0; idtrat<ndtrat; idtrat++){
	/*Evaluate noise covariance matrix of discrete state, and measurement of the
	 * state. This block takes most of the preparation step*/
	int dtrat=dtrats->p[idtrat];
	double dT=dthi*dtrat;
	int nsec=dtrat*10;
	double dT2=dT/nsec;
	dmat *expAj=0;
	dmat *expAjn=0;
	dmat *tmp1=0, *tmp2=0;
	for(int i=0; i<nsec; i++){
	    double ti=(i+0.5)*dT2;
	    if(idtrat==0){
		dexpm(&expAj, 0, Ac, ti);
		dmm(&tmp1, 0, expAj, Sigma_ep, "nn", 1);
		dmm(&Sigma_varep, 1, tmp1, expAj, "nt", 1);
	    }
	    dexpm(&expAjn, 0, Ac, -ti);
	    dscale(expAjn, -1); daddI(expAjn, 1); //1-exp(-A*tj)
	    dmm(&tmp1, 0, expAjn, AcI, "nn", 1);
	    dmm(&tmp2, 0, tmp1, Sigma_ep, "nn", 1);
	    dmm(&Sigma_zeta->p[idtrat], 1, tmp2, tmp1, "nt", 1);
	}
	if(idtrat==0){
	    dscale(Sigma_varep, dT2);
	}
	dscale(Sigma_zeta->p[idtrat], dT2/(dT*dT));
	{
	    dmat *tmp=0;
	    dmm(&tmp, 0, Sigma_zeta->p[idtrat], Pd, "nt", 1);
	    dmm(&Radd->p[idtrat], 0, Pd, tmp, "nn", 1);
	    Raddchol->p[idtrat]=dchol(Radd->p[idtrat]);
	}
	{
	    dmat *tmp=0;
	    dexpm(&tmp, 0, Ac, -dT);
	    dscale(tmp, -1); 
	    daddI(tmp, 1);
	    dmm(&Xi->p[idtrat], 0, tmp, AcI, "nn", 1./dT);
	    dfree(tmp);
	}   
	dfree(tmp1); dfree(tmp2);
	dfree(expAj); dfree(expAjn);
    }
    dmat *Ad=0; dexpm(&Ad, 0, Ac, dthi*dtrat_min);  /*discrete state propagation at dT*/
    dmat *AdM=0; dexpm(&AdM, 0, Ac, dthi);/*discrete state propagation at dthi*/
    dmat *FdM=0;/*From discrete state to averaged mode for dthi*/
    {//Compute FdM=Pd*[exp(-Ac*dthi)/dthi]*Ac^-1.
	dmat *XiM=0;
	dmat *tmp=0;
	dexpm(&tmp, 0, Ac, -dthi);
	dscale(tmp, -1); 
	daddI(tmp, 1);
	dmm(&XiM, 0, tmp, AcI, "nn", 1./dthi);
	dmm(&FdM, 0, Pd, XiM, "nn", 1);
	dfree(XiM);
	dfree(tmp);
    }
    int nkalman=(1<<(dtrat_wfs->nx))-1;
    kalman_t *res=mycalloc(1,kalman_t);
    res->Ad=Ad;
    res->AdM=AdM;
    res->FdM=FdM;
    res->dthi=dthi;
    res->dtrat=ldup(dtrat_wfs);
    res->Gwfs=dcelldup(Gwfs);//Used for simulation. Do not mess with it.
    res->Rwfs=dcelldup(Rwfs);
    res->Rn=dcellnew(nkalman, 1);
    res->Qn=dref(Sigma_varep);
    res->Cd=dcellnew(nkalman, 1);
    
    /*Loop over first set of steps to find needed kalman filter.*/
    for(int istep=0; istep<dtrat_prod; istep++){
	int indk=0;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int dtrat=dtrat_wfs->p[iwfs];
	    if((istep+1) % dtrat == 0){
		indk|=1<<iwfs;/*this is how we compute the index into kalman*/
	    }
	}
	if(indk && !res->Cd->p[indk-1]){
	    dcell *Rnadd=dcellnew(nwfs, nwfs);
	    dcell *Cd=dcellnew(nwfs, 1);
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		dmat *Gwfsi=0;//Use for reconstruction. 
		if(indk==1 && nmod==12 && KALMAN_IN_SKYC){
	            /*Our sky coverage case with single WFS, make 3rd column zero*/
	            Gwfsi=ddup(Gwfs->p[iwfs]);
		    memset(Gwfsi->p+2*Gwfsi->nx, 0, Gwfsi->nx*sizeof(double));
		}else{
		    Gwfsi=dref(Gwfs->p[iwfs]);
		}
		int dtrat=dtrat_wfs->p[iwfs];
		Cd->p[iwfs]=dnew(Gwfsi->nx, Ac->ny);
		if((istep+1) % dtrat == 0){
		    //WFS active. Locate the index into Xi, Radd.
		    int jdtrat;
		    for(jdtrat=0; jdtrat<ndtrat; jdtrat++){
			if(dtrats->p[jdtrat]==dtrat){
			    break;
			}
		    }
		    if(jdtrat<ndtrat){
			dmat *Fd=0;
			dmm(&Fd, 0, Pd, Xi->p[jdtrat], "nn", 1);
			dmm(&Cd->p[iwfs], 1, Gwfsi, Fd, "nn", 1);
			dfree(Fd);
			//Here Radd is chol of covariance
			dmat *tmp=0;
			dmm(&tmp, 0, Gwfsi, Radd->p[jdtrat], "nn", 1);
			dmm(&Rnadd->p[iwfs+iwfs*nwfs], 1, tmp, Gwfsi, "nt", 1);
		    }else{
			error("not found\n");
		    }
		}
		dfree(Gwfsi);
	    }
	    dcell *Rn=0;
	    /*compute Rn=Rwfs+Gwfs*Radd*Gwfs' */
	    dcelladd(&Rn, 1, Rwfs, 1);
	    dcelladd(&Rn, 1, Rnadd, 1);
	   
	    res->Rn->p[indk-1]=dcell2m(Rn);
	    res->Cd->p[indk-1]=dcell2m(Cd);
	    dcellfree(Cd);
	    dcellfree(Rnadd);
	    dcellfree(Rn);
	}
    }
    res->M=reccati_cell(&res->P, Ad, res->Qn, res->Cd, res->Rn);
    {
	/*convert estimation error in modes */
	dmat *tmp1=0;
	dmm(&tmp1, 0, Pd, res->P, "nn", 1);
	dfree(res->P);
	dmm(&res->P, 0, tmp1, Pd, "nt", 1);
	dfree(tmp1);
    }
    lfree(dtrats);
    dcellfree(Xi);
    dfree(Ac);
    dfree(Sigma_ep);
    dfree(Pd);
    dfree(Sigma_varep);
    dcellfree(Sigma_zeta);
    dfree(AcI);
    dcellfree(Radd);
    dcellfree(Raddchol);
    return res;
}
void kalman_free(kalman_t *kalman){
    if(!kalman) return;
    dfree(kalman->Ad);
    dcellfree(kalman->Cd);
    dfree(kalman->AdM);
    dfree(kalman->FdM);
    dfree(kalman->Qn);
    dcellfree(kalman->M);
    dfree(kalman->P);
    lfree(kalman->dtrat);
    dcellfree(kalman->Gwfs);
    dcellfree(kalman->Rwfs);
    dcellfree(kalman->Rn);
    dfree(kalman->xhat);
    dfree(kalman->xhat2);
    dfree(kalman->xhat3);
    free(kalman);
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
/*Update state vector when there is a measurement. It modifies meas.
    xhat3=xhat+M*(meas-Cd*xhat)
    xhat2=AdM*xhat3;
    xhat =Ad*xhat3;
 */
void kalman_update(kalman_t *kalman, dmat *meas, int ik){
    /*difference between real and predicted measurement*/
    dmm(&meas, 1, kalman->Cd->p[ik], kalman->xhat, "nn", -1);
    /*updated estimate*/
    dcp(&kalman->xhat3, kalman->xhat);
    dmm(&kalman->xhat3, 1, kalman->M->p[ik], meas, "nn", 1);
    dmm(&kalman->xhat, 0, kalman->Ad, kalman->xhat3, "nn", 1);/*predicted state*/
    dmm(&kalman->xhat2,0, kalman->AdM, kalman->xhat3, "nn", 1);/*correction for this time step*/
}
/*Output correction: out=out*alpha+Fdm*(AdM*xhat2)*beta */
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
    dcell *Gwfs=kalman->Gwfs;
    const int nwfs=Gwfs->nx;
    dcell *rmsn=dcellnew(nwfs,1);
    dcell *noise=dcellnew(nwfs,1);
    long ngs[nwfs];
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ng=Gwfs->p[iwfs]->nx;
	if(dsumsq(kalman->Rwfs->p[iwfs+nwfs*iwfs])>0){
	    rmsn->p[iwfs]=dchol(kalman->Rwfs->p[iwfs+nwfs*iwfs]);
	    noise->p[iwfs]=dnew(ng, 1);
	}
	ngs[iwfs]=ng;
    }
    dcell *acc=dcellnew3(nwfs, 1, ngs, 0);
    dcell *meas=dcellnew3(nwfs, 1, ngs, 0);
    int nmod=input->nx;
    dmat *mres=ddup(input);
    kalman_init(kalman);
    rand_t rstat;
    seed_rand(&rstat, 1);
    dmat *outi=dnew_ref(nmod, 1, (double*)1);//wraper
    dcell *inic=dcellnew(1,1); inic->p[0]=dnew_ref(nmod,1,(double*)1);//wrapper
    dmat *ini=inic->p[0];
    for(int istep=0; istep<input->ny; istep++){
	ini->p=PCOL(input, istep);
	outi->p=PCOL(mres,istep);
	kalman_output(kalman, &outi, 1, -1);
	int indk=0;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int dtrat=(int)kalman->dtrat->p[iwfs];
	    if(istep && (istep) % dtrat == 0){/*There is measurement from last step*/
		indk|=1<<iwfs;
		//Average the measurement.
		dadd(&meas->p[iwfs], 0, acc->p[iwfs], 1./dtrat); dzero(acc->p[iwfs]);
		//Add noise
		if(rmsn->p[iwfs]){
		    drandn(noise->p[iwfs], 1, &rstat);
		    if(rmsn->p[iwfs]->nx>1){
			dmm(&meas->p[iwfs], 1, rmsn->p[iwfs], noise->p[iwfs], "nn", 1);
		    }else{
			dadd(&meas->p[iwfs], 1, noise->p[iwfs], rmsn->p[iwfs]->p[0]);
		    }
		}
	    }
	}
	if(indk){
	    kalman_update(kalman, meas->m, indk-1);
	}
	dcellmm(&acc, Gwfs, inic, "nn", 1);/*OL measurement*/
    }
    dcellfree(rmsn);
    dcellfree(meas);
    dcellfree(acc);
    dcellfree(noise);
    dfree(outi);
    dcellfree(inic);
    return mres;
}
/**
   Save kalman_t to file
*/
void kalman_write(kalman_t *kalman, const char *format, ...){
    format2fn;
    file_t *fp=zfopen(fn, "wb");
    if(kalman){
	header_t header={MCC_ANY, 12, 1, (char*)"type=struct"};
	write_header(&header, fp);
	char *tmp;
#define WRITE_KEY(fp, str, key)			\
	tmp=str->key->header;			\
	str->key->header=(char*)#key;		\
	writebindata(fp, str->key);		\
	str->key->header=tmp;

	WRITE_KEY(fp, kalman, Ad);
	WRITE_KEY(fp, kalman, Cd);
	WRITE_KEY(fp, kalman, AdM);
	WRITE_KEY(fp, kalman, FdM);
	
	WRITE_KEY(fp, kalman, Qn);
	WRITE_KEY(fp, kalman, M);
	WRITE_KEY(fp, kalman, P);
	WRITE_KEY(fp, kalman, dtrat);
	
	WRITE_KEY(fp, kalman, Gwfs);
	WRITE_KEY(fp, kalman, Rwfs);
	WRITE_KEY(fp, kalman, Rn);
	writearr(fp, 0, sizeof(double), M_DBL, "dthi", &kalman->dthi, 1, 1);
    }
    zfclose(fp);
}
