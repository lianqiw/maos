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
    double df;
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
static void psd2cov(dmat *psd, double df){
    dfft1plan_r2hc(psd, -1);
    dscale(psd, df);
    psd->p[0]=0;//zero dc. critical
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
	    return 0;
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
    int check;
    double diff;
    if(!(check=coeff_isgood(coeff, data->ncoeff*data->nmod, data->min, data->max))){
	diff=(1+fabs(check))*1e50;
    }else{
	sde_psd(&data->cov_sde, data->f2, coeff, data->ncoeff, data->nmod);
	psd2cov(data->cov_sde, data->df);
	dmat *cov1=dnew_ref(data->ncov, 1, data->cov_sde->p);
	dmat *cov2=dnew_ref(data->ncov, 1, data->cov_in->p);
	data->ratio=cov2->p[0]/cov1->p[0];
	/*scale to same max value*/
	dadd(&cov1, 1./cov1->p[0], cov2, -1./cov2->p[0]);
	diff=dnorm2(cov1)/data->ncov;
	/*dwrite(cov1, "diff");*/
	/*{
	    static int count=0; count++;
	    info("n=%d, diffrms=%g, ratio=%g\n", data->ncov, diff, data->ratio);
	}*/
	dfree(cov1); dfree(cov2);
    }
    return diff;
}
/**
   Fit a PSD with SDE model. Use initial guess from coeff0 and return final
   answer in the same format.
*/
dmat* sde_fit(const dmat *psdin, const dmat *coeff0, double tmax_fit, double min, double max, double df){
    if(psdin->ny!=2){
	error("psd must contain nu and psd\n");
    }
    int nf;
    double maxf=psdin->p[psdin->nx-1];
    if(df<EPS){
	nf=1024*4; /**array length for FFT*/
	df=maxf/(nf/2-1);
    }else{
	nf=round(maxf/df+1)*2;
	info2("nf=%d\n", nf);
    }
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
    dmat *scale=dnew(ncoeff, nmod);dadd(&scale, 0, coeff0, 1);
    dmat *psd_sde=dnew(nf/2, 1);
    dmat *cov_sde=dnew(nf, 1);
    psd2cov(cov, df);
    if(!isfinite(cov->p[0])){
	dmat *psd0=psdinterp1(psdin, f2);
	dwrite(cov, "bad_cov");
	dwrite(psd0, "bad_psd");
    }
    int ncov=round(tmax_fit/dt);
    if(ncov>nf/2){
	ncov=nf/2;
    }
    sde_fit_t data={df, f1, f2, psd, cov, psd_sde, cov_sde, min, max, 0, ncoeff, nmod, ncov};
    if(data.ncov>0){//covariance fitting
	double tol=1e-15;
	dminsearch(coeff->p, scale->p, ncoeff*nmod, tol, sde_diff_cov, &data);
	sde_diff_cov(coeff->p, &data);//needed to set data.ratio with latest result.
	for(int im=0; im<nmod; im++){
	    coeff->p[(1+im)*ncoeff-1]*=sqrt(data.ratio);
	}
	/*info("fitting done\n");
	sde_diff_cov(coeff->p, &data);
	dwrite(data.cov_in, "cov_psd");
	sde_psd(&data.cov_sde, data.f2, coeff->p, data.ncoeff, data.nmod);
	dwrite(data.cov_sde, "psd_sde");
	psd2cov(data.cov_sde, df);
	dwrite(data.cov_sde, "cov_sde");*/
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
void reccati(dmat **Mout, dmat **Pout, 
	     const dmat *A, const dmat *Qn, const dmat *C, const dmat *Rn){
    double diff=1, diff2=1, lastdiff=INFINITY;
    dmat *P2=dnew(A->nx, A->ny); daddI(P2, Qn->p[0]);//Initialize P to identity
    dmat *AP=0, *CP=0, *P=0, *CPAt=0, *CPCt=0, *APCt=0, *tmp=0;
    int count=0;
    while(diff>1e-13 && diff2>1e-12 && count++<1000){
	/*Notice that P may not be symmetric due to round off errors.*/
	dcp(&P, P2);
	dmm(&AP, 0, A, P, "nn", 1);
	dmm(&CP, 0, C, P, "nn", 1);
	dmm(&CPAt, 0, CP, A, "nt", 1);
	dmm(&APCt, 0, AP, C, "nt", 1);
	dmm(&CPCt, 0, CP, C, "nt", 1);
	dadd(&CPCt, 1, Rn, 1);
	dsvd_pow(CPCt, -1, 1e-12);
	dmm(&P2, 0, AP, A, "nt", 1);
	dmm(&tmp, 0, CPCt, CPAt, "nn", 1);
	dmm(&P2, 1, APCt, tmp, "nn", -1);
	dadd(&P2, 1, Qn, 1);//new P
	diff=dnorm2(P);
	dadd(&P, 1, P2, -1);
	diff=sqrt(dnorm2(P)/diff);
	diff2=fabs(diff-lastdiff);
	lastdiff=diff;
	//info2("diff=%g, diff2=%g\n", diff, diff2);
    }
    dmm(Mout, 0, CP, CPCt, "tn", 1);
    if(Pout) dcp(Pout, P);
    dfree(AP);dfree(CP); dfree(P); dfree(P2); 
    dfree(CPAt); dfree(CPCt); dfree(APCt); dfree(tmp);
}
/*Kalman filter based on SDE model*/
kalman_t* sde_kalman(dmat *coeff, /**<SDE coefficients*/
		     double dthi, /**<Loop frequency*/
		     dmat *dtrat_wfs,   /**<WFS frequency as a fraction of loop*/
		     dcell *Gwfs,  /**<WFS measurement from modes. Can be identity*/
		     dcell *Rwfs,  /**<WFS measurement noise covariance*/
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
    dmat *AcI=ddup(Ac); 
    dsvd_pow(AcI, -1, 0);
    //make sure dtrat_wfs is column vector
    dtrat_wfs->nx=dtrat_wfs->nx*dtrat_wfs->ny; 
    dtrat_wfs->ny=1;
    const int nwfs=Gwfs->nx;
    int ndtrat=0;
    dmat *dtrats=dnew(nwfs, 1);
    if(dtrat_wfs->nx==1 || dtrat_wfs->nx==nwfs){
	ndtrat=0;
	for(int iwfs=0; iwfs<dtrat_wfs->nx; iwfs++){
	    int found=0;
	    for(int jdtrat=0; jdtrat<ndtrat; jdtrat++){
		if((int)dtrats->p[jdtrat]==(int)dtrat_wfs->p[iwfs]){
		    found=1; break;
		}
	    }
	    if(!found){
		dtrats->p[ndtrat++]=dtrat_wfs->p[iwfs];
	    }
	}
	dresize(dtrats, ndtrat, 1);
	dsort(dtrats, 1);
	for(int idtrat=1; idtrat<ndtrat; idtrat++){
	    if((int)dtrats->p[idtrat]%(int)dtrats->p[idtrat-1]!=0){
		error("dtrat_wfs is in the wrong format\n");
	    }
	}
    }else{
	error("dtrat_wfs is in wrong format. Should have length of either 1 or %d\n", nwfs);
    }
    const int dtrat_min=dtrats->p[0];
    dcell* Sigma_varep=dcellnew(ndtrat,1);/*state noise in discrete space*/
    dcell* Sigma_zeta =dcellnew(ndtrat,1);/*state measurement noise*/
    dcell *Radd       =dcellnew(ndtrat,1);
    dcell *Xi=dcellnew(ndtrat,1);
    for(int idtrat=0; idtrat<ndtrat; idtrat++){
	/*Evaluate noise covariance matrix of discrete state, and measurement of the
	 * state*/
	int dtrat=(int)dtrats->p[idtrat];
	double dT=dthi*dtrat;
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
	    dmm(&Sigma_varep->p[idtrat], 1, tmp1, expAj, "nt", 1);
	    
	    dscale(expAjn, -1); daddI(expAjn, 1); //1-exp(-A*tj)
	    dmm(&tmp1, 0, expAjn, AcI, "nn", 1);
	    dmm(&tmp2, 0, tmp1, Sigma_ep, "nn", 1);
	    dmm(&Sigma_zeta->p[idtrat], 1, tmp2, tmp1, "nt", 1);
	}
	dscale(Sigma_varep->p[idtrat], dT2);
	dscale(Sigma_zeta->p[idtrat], dT2/(dT*dT));
	dfree(tmp1); dfree(tmp2);
	dfree(expAj); dfree(expAjn);
	{
	    dmat *tmp=0;
	    dmm(&tmp, 0, Sigma_zeta->p[idtrat], Pd, "nt", 1);
	    dmm(&Radd->p[idtrat], 0, Pd, tmp, "nn", 1);
	    dfree(tmp);
	    tmp=dchol(Radd->p[idtrat]);
	    dfree(Radd->p[idtrat]);
	    Radd->p[idtrat]=tmp; tmp=0;
	}
	{
	    dmat *tmp=0;
	    dexpm(&tmp, 0, Ac, -dT);
	    dscale(tmp, -1); 
	    daddI(tmp, 1);
	    dmm(&Xi->p[idtrat], 0, tmp, AcI, "nn", 1./dT);
	    dfree(tmp);
	}   
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
    const int nkalman=nwfs;
    kalman_t *res=calloc(1, sizeof(kalman_t));
    res->Ad=Ad;
    res->AdM=AdM;
    res->FdM=FdM;
    res->M=dcellnew(nkalman, 1);
    res->P=dcellnew(nkalman, 1);
    res->dthi=dthi;
    res->dtrat=ddup(dtrat_wfs);
    res->Gwfs=dcelldup(Gwfs);
    res->Rwfs=dcelldup(Rwfs);
    res->Rn=dcellnew(nkalman, 1);
    res->Qn=dref(Sigma_varep->p[0]);
    res->Cd=dcellnew(nkalman, 1);
    for(int idtrat=0; idtrat<ndtrat; idtrat++){
	int dtrat=(int)dtrats->p[idtrat];
	dcell *Rn=dcelldup(Rwfs);
	dcell *Rnadd=dcellnew(nwfs, 1);
	dcell *Cd=dcellnew(nwfs, 1);
	int active_wfs=0;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    Cd->p[iwfs]=dnew(Gwfs->p[iwfs]->nx, Ac->ny);
	    if(dtrat % (int)dtrat_wfs->p[iwfs]==0){
		active_wfs++;
		int jdtrat;
		for(jdtrat=0; jdtrat<ndtrat; jdtrat++){
		    if((int)dtrats->p[jdtrat]==(int)dtrat_wfs->p[iwfs]){
			break;
		    }
		}
		if(jdtrat<ndtrat){
		    dmat *Fd=0;
		    dmm(&Fd, 0, Pd, Xi->p[jdtrat], "nn", 1);
		    dmm(&Cd->p[iwfs], 1, Gwfs->p[iwfs], Fd, "nn", 1);
		    dfree(Fd);
		    dmm(&Rnadd->p[iwfs], 1, Gwfs->p[iwfs], Radd->p[jdtrat], "nn", 1);
		}else{
		    error("not found\n");
		}
	    }
	    //Here Radd is chol of covariance
	}
	/*compute Rn=Rwfs+Gwfs*Radd*Gwfs' */
	dcellmm(&Rn, Rnadd, Rnadd, "nt", 1);
	int ik=active_wfs-1;
	dmat *Rnm=res->Rn->p[ik]=dcell2m(Rn);
	dmat *Cdm=dcell2m(Cd);
	dmat *M=0, *P=0;
	reccati(&M, &P, Ad, res->Qn, Cdm, Rnm);
	{
	    /*convert estimation error in modes */
	    dmat *tmp1=0;
	    dmm(&tmp1, 0, Pd, P, "nn", 1);
	    dfree(P);
	    dmm(&P, 0, tmp1, Pd, "nt", 1);
	    dfree(tmp1);
	}
	res->M->p[ik]=M;M=0;
	res->P->p[ik]=P;P=0;
	res->Cd->p[ik]=dref(Cdm);
	dcellfree(Cd);
	dcellfree(Rnadd);
	dcellfree(Rn);
	dfree(Cdm);
    }
    dfree(dtrats);
    dcellfree(Xi);
    dfree(Ac);
    dfree(Sigma_ep);
    dfree(Pd);
    dcellfree(Sigma_varep);
    dcellfree(Sigma_zeta);
    dfree(AcI);
    dcellfree(Radd);
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
    dcellfree(kalman->P);
    dfree(kalman->dtrat);
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
    dcell *rmsn=dcellnew(nwfs,nwfs);
    dcell *noise=dcellnew(nwfs,1);
    long ngs[nwfs];
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ng=Gwfs->p[iwfs]->nx;
	if(dnorm2(kalman->Rwfs->p[iwfs+nwfs*iwfs])>0){
	    dmat *tmp=dcell2m(kalman->Rwfs);
	    rmsn->p[iwfs+nwfs*iwfs]=dchol(kalman->Rwfs->p[iwfs+nwfs*iwfs]);
	    dfree(tmp);
	    noise->p[iwfs]=dnew(ng, 1);
	}
	ngs[iwfs]=ng;
    }
    dcell *acc=dcellnew3(nwfs, 1, ngs, 0);
    dcell *meas=dcellnew3(nwfs, 1, ngs, 0);
    int nmod=input->nx;
    PDMAT(input,pinput);
    dmat *mres=ddup(input);
    PDMAT(mres, pmres);
    kalman_init(kalman);
    rand_t rstat;
    seed_rand(&rstat, 1);
    dmat *outi=dnew_ref(nmod, 1, (double*)1);//wraper
    dcell *inic=dcellnew(1,1); inic->p[0]=dnew_ref(nmod,1,(double*)1);//wrapper
    dmat *ini=inic->p[0];
    for(int istep=0; istep<input->ny; istep++){
	ini->p=(double*)(pinput+istep);
	outi->p=(double*)(pmres+istep);
	kalman_output(kalman, &outi, 1, -1);
	int active_wfs=0;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int dtrat=(int)kalman->dtrat->p[iwfs];
	    if(istep && (istep) % dtrat == 0){/*There is measurement from last step*/
		dadd(&meas->p[iwfs], 0, acc->p[iwfs], 1./dtrat); dzero(acc->p[iwfs]);
		active_wfs++;
		if(rmsn->p[iwfs+nwfs*iwfs]){
		    drandn(noise->p[iwfs], 1, &rstat);
		    if(rmsn->p[iwfs+nwfs*iwfs]->nx>1){
			dmm(&meas->p[iwfs], 1, rmsn->p[iwfs+nwfs*iwfs], noise->p[iwfs], "nn", 1);
		    }else{
			dadd(&meas->p[iwfs], 1, noise->p[iwfs], rmsn->p[iwfs+nwfs*iwfs]->p[0]);
		    }
		}
	    }
	}
	if(active_wfs){
	    kalman_update(kalman, meas->m, active_wfs-1);
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
