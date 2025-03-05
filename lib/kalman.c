/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
	real df;
	dmat* freq;
	dmat* psdcov_in;
	dmat* psdcov_sde;
	real norm_in;
	int count;
	int ncoeff;
	int nmod;
	int ncov;
}sde_fit_t;
/**
   Do FFT to convert PSD to Covariance
 */
static void psd2cov(dmat* psd, real df){
	dscale(psd, df);
	P(psd,0)=0;//zero dc. critical
	dfft1plan_r2hc(psd, -1);
	dfft1(psd, -1);
}
/**
   Computes PSD of SDE with coefficients (coeff) and frequency (f)
 */
void sde_psd(dmat** psd, const dmat* f, const real* coeff, int ncoeff, int nmod){
	int order=ncoeff-1;
	if(!*psd){
		*psd=dnew(NX(f), 1);
	} else{
		dzero(*psd);
	}
	dmat* ppsd=*psd;
	if(NX(ppsd)!=NX(f)){
		warning("f and psd does not match\n");
	}
	if(order==2){//special version for order 2. faster.
		real twopi=2*M_PI;
		for(int im=0; im<nmod; im++){
			real c1sq=pow(coeff[0+im*3], 2);
			real c2=coeff[1+im*3];
			real sigma2=pow(coeff[2+im*3], 2);
			real* pp;
			if(NY(ppsd)==nmod){//seperate
				pp=PCOL(ppsd, im);
			} else{
				pp=P(ppsd);
			}
			for(int is=0; is<NX(f); is++){
				real omega2=pow(P(f,is)*twopi, 2);
				pp[is]+=sigma2/(pow(omega2-c2, 2)+c1sq*omega2);//2017-07-21: was =
			}
		}
	} else{
		comp twopii=COMPLEX(0, 2*M_PI);
		for(int is=0; is<NX(f); is++){
			comp s=P(f,is)*twopii;
			for(int im=0; im<nmod; im++){
				comp denom=cpow(s, order);
				for(int ic=0; ic<order-1; ic++){
					denom+=coeff[ic+im*ncoeff]*cpow(s, order-ic-1);
				}
				denom+=coeff[order-1+im*ncoeff];
				real val=pow(coeff[order+im*ncoeff], 2)/cabs2(denom);
				if(NY(ppsd)==nmod){//seperate
					P(ppsd, is, im)+=val;
				} else{
					P(ppsd,is)+=val;
				}
			}
		}
	}
}
/**
   A convenient wrapper for sde_psd()
*/
dmat* sde_psd2(const dmat* ff, const dmat* coeff){
	dmat* psd=0;
	sde_psd(&psd, ff, P(coeff), NX(coeff), NY(coeff));
	return psd;
}
static real coeff_isbad(const real* coeff, int ncoeff){
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

static real sde_diff(const real* coeff, void* pdata){
	sde_fit_t* data=(sde_fit_t*)pdata;
	real tmp;
	real diff;
	data->count++;
	if((tmp=coeff_isbad(coeff, data->ncoeff*data->nmod))){
		diff=(tmp)*10+10;
	} else{
		sde_psd(&data->psdcov_sde, data->freq, coeff, data->ncoeff, data->nmod);
		if(!data->ncov){//PSD fitting
			dadd(&data->psdcov_sde, 1, data->psdcov_in, -1);
			/*Must use dnorm here. dsumabs does not work well for peaks when
			 * initial guess of c1 is not close.*/
			diff=dnorm(data->psdcov_sde)/data->norm_in;
		} else{//Covariance fitting.
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
				real val1=P(data->psdcov_in,i);
				real val2=P(data->psdcov_sde,i);
				diff+=pow(val1-val2, 2);
			}
			diff=diff/((abs2(P(data->psdcov_in,0)+P(data->psdcov_sde,0)))*data->ncov);
#elif METHOD==2
			real scale=P(data->psdcov_in,0)/P(data->psdcov_sde,0);
			for(long i=0; i<data->ncov; i++){
				real val1=P(data->psdcov_in,i);
				real val2=P(data->psdcov_sde,i);
				diff+=pow(val1-val2*scale, 2);
			}
			diff+=abs2(P(data->psdcov_in,0)-P(data->psdcov_sde,0))*data->ncov;
			diff=diff/((abs2(P(data->psdcov_in,0)))*data->ncov);
#endif
		}
	}
	if(isinf(diff)){
		error("Diff is not finite\n");
	}
	/*info2("coeff=");
	for(int imod=0; imod<data->nmod; imod++){
	info2("[%g %.2f %g] ", coeff[3*imod], coeff[3*imod+1], coeff[3*imod+2]);
	}
	info2("ncov=%d. cov0=[%g, %g] diff=%g. count=%d\n", data->ncov, P(data->psdcov_in,0), P(data->psdcov_sde,0), diff, data->count);
	*/
	return diff;
}
/**
   Scale the third coefficient (sigma) to have the same energy
*/
static void sde_scale_coeff(dmat* coeff, real var_in, dmat* psdcov_sde, const dmat* freq, int ncov){
	for(int iy=0; iy<NY(coeff); iy++){
		if(P(coeff,2,iy)<=0){
			P(coeff,2,iy)=1;
		}
	}
	real ratio;
	sde_psd(&psdcov_sde, freq, P(coeff), NX(coeff), NY(coeff));
	real var_sde;
	if(ncov){
	//writebin(psdcov_sde, "sde_psd");
		psd2cov(psdcov_sde, P(freq,2)-P(freq,1));
		//writebin(psdcov_sde, "sde_cov");
		var_sde=P(psdcov_sde,0);
	} else{
		var_sde=dtrapz(freq, psdcov_sde);
	}
	ratio=sqrt(var_in/var_sde);
	for(int iy=0; iy<NY(coeff); iy++){
		P(coeff,2,iy)*=ratio;
	}
}
/**
   Fit a PSD with SDE model. Use initial guess from coeff0 and return final
   answer in the same format.
*/
static dmat* sde_fit_do(const dmat* psdin, const dmat* coeff0, real tmax_fit){
	if(NY(psdin)!=2){
		error("psd must contain nu and psd\n");
	}
	real df=0.01;//ensure we can sample narrow peaks when interpolating psd
	dmat* freq=0, * psdcov_in=0, * psdcov_sde=0;
	int ncov=0;
	real var_in=0;//strength of input
	/* 2014-07-31: Always upsample the PSD to 0.01 Hz sampling. Otherwise SDE
	   PSD may have unresolved peaks, causing fitting to misbehave.
	*/
	real maxf=P(psdin, psdin->nx-1);
	real norm_in;
	if(tmax_fit>0){//Use covariance fit
	/* Create Frequency vector with proper fft indexing [0:1:nf/2-1 * nf/2:-1:1].
	   Negative index is mapped to positive. */
		long nf=round(maxf/df)*2;
		ncov=MIN(round(tmax_fit*maxf), nf/2);
		freq=dnew(nf, 1);
		for(long i=0; i<nf/2; i++){
			P(freq,i)=df*i;
			P(freq, i+nf/2)=df*(nf/2-i);
		}
		psdcov_in=psd_interp1(psdin, freq, 0);

		//writebin(psdcov_in, "in_psd");
		psd2cov(psdcov_in, df);
		var_in=P(psdcov_in,0);
		if(isinf(P(psdcov_in,0))){
			warning("covariance is not finite at 0\n");
			writebin(psdcov_in, "bad_cov");
		}
		//writebin(psdcov_in, "in_cov");
		psdcov_sde=dnew(nf, 1);
		norm_in=ncov*abs2(P(psdcov_in,0));
	} else{//Use PSD fit
		if((P(psdin,1)-P(psdin,0))<df*2){
			freq=drefcols(psdin, 0, 1);
			psdcov_in=drefcols(psdin, 1, 1);
		} else{//Resampling.
			real minf=P(psdin,0);
			long nf=round((maxf-minf)/df);
			freq=dnew(nf, 1);
			for(long i=0; i<nf; i++){
				P(freq,i)=minf+df*i;
			}
			psdcov_in=psd_interp1(psdin, freq, 0);
		}
		psdcov_sde=dnew(NX(freq), 1);
		var_in=dtrapz(freq, psdcov_in);
		norm_in=dnorm(psdcov_in);
	}
	int ncoeff=NX(coeff0);
	int nmod=NY(coeff0);
	dmat* coeff=dnew(ncoeff, nmod);dcp(&coeff, coeff0);
	//Scale to make sure total energy is preserved.
	sde_scale_coeff(coeff, var_in, psdcov_sde, freq, ncov);

	sde_fit_t data={df, freq, psdcov_in, psdcov_sde, norm_in, 0, ncoeff, nmod, ncov};
	real diff0=sde_diff(P(coeff), &data);
	real tol=1e-10;
	int nmax=2000;
	dminsearch(P(coeff), ncoeff*nmod, tol, nmax, sde_diff, &data);
	//Do not scale coeff after the solution.
	real diff1=sde_diff(P(coeff), &data);
	info2("sde_fit: %d interations: %g->%g.\n", data.count, diff0, diff1);
	//Scale to make sure total energy is preserved.
	/*
	  if(diff1>0.2 && diff1>diff0*0.75){
	static int count=0;
	writebin(psdin, "sde_fit_psdin_%d_%g", count, tmax_fit);
	writebin(coeff0,"sde_fit_coeff_%d_%g", count, tmax_fit);
	count++;
	if(tmax_fit>0){
		info2("Redo with PSD fitting.\n");
		dfree(coeff);
		coeff=sde_fit_do(psdin, coeff0, 0);
	}else{
		warning("Failed to converge.\n");
	}
	}else{
	info2("\n");
	}*/
	dfree(freq);
	dfree(psdcov_in);
	dfree(psdcov_sde);
	return coeff;
}

/**
   Estiamte the total PSD power for vibration peaks using FWHM*peak
 */
/*static real sde_vib_est(real c1, real c2){
	real sqrt4ac=sqrt(c1*c1+c2*4);
	real omega1=sqrt(0.5*((c2*2-c1*c1)-sqrt4ac));
	real omega2=sqrt(0.5*((c2*2-c1*c1)+sqrt4ac));
	real fwhm=(omega2-omega1)/(2*M_PI);
	real peak=1./(c1*c1*c2);
	return fwhm*peak;
	}*/
/**
   If coeff0 is not null, use it immediately, otherwise, do vibration identification
 */
dmat* sde_fit(const dmat* psdin, const dmat* coeff0, real tmax_fit, int vibid){
	if(coeff0){
		return sde_fit_do(psdin, coeff0, tmax_fit);
	} else{
	//Do vibration identification
		dmat* vibs=vibid?psd_vibid(psdin):NULL;
		dcell* coeffs=dcellnew(1, vibs?(1+NY(vibs)):1);
		dmat* coeffi=dnew(3, 1);
		dmat* psd2=ddup(psdin);
		if(vibs&&NY(vibs)>0){
			dbg("\nnvib=%ld\n", NY(vibs));
			for(int ivib=0; ivib<NY(vibs); ivib++){
				real fi=P(vibs,ivib*4+0);
				int i1=P(vibs,ivib*4+2);
				int i2=P(vibs,ivib*4+3)+1;
				if(i2-i1<3){
					i2++;
					i1--;
				}
				dmat* psdi=dsub(psdin, i1, i2-i1, 0, 0);//extract peak
				P(coeffi,0)=1;
				P(coeffi,1)=pow(2*M_PI*fi, 2);
				P(coeffi,2)=0;//sde_fit_do will figure it out.
				P(coeffs,ivib)=sde_fit_do(psdi, coeffi, tmax_fit);
				if(fabs(sqrt(P(coeffi,1))-sqrt(P(P(coeffs,ivib),1)))>2*M_PI){
					warning("Fitting failed for %d, Freq=%g, %g\n", ivib,
						sqrt(P(coeffi,1))/(2*M_PI), sqrt(P(P(coeffs,ivib),1))/(2*M_PI));
					writebin(psdi, "psdi_%d", ivib);
					writebin(coeffi, "coeffi_%d", ivib);
					writebin(P(coeffs,ivib), "coeffo_%d", ivib);
					if(fabs(sqrt(P(coeffi,1))-sqrt(P(P(coeffs,ivib),1)))>2*M_PI*5){
						dcp(&P(coeffs,ivib), coeffi);
						warning("Use input\n");
					}
				}
				real* psd2p=PCOL(psd2, 1);
				for(int i=i1; i<i2; i++){//replace the peak by a linear interpolation
					psd2p[i]=(psd2p[i1-1]*(i2+1-i)+psd2p[i2+1]*(i-(i1-1)))/(i2-i1+2);
				}
			}
			dfree(vibs);
		}
		P(coeffi,0)=1;
		P(coeffi,1)=1;
		P(coeffi,2)=0;//sde_fit_do will figure it out.
		P(coeffs,coeffs->ny-1)=sde_fit_do(psd2, coeffi, tmax_fit);
		dfree(coeffi);
		dfree(psd2);
		dmat* coeff=dcell2m(coeffs);
		dcellfree(coeffs);
		return coeff;
	}
}

/**
   Compute the reccati equation.
 */
dmat* reccati(dmat** Pout, const dmat* A, const dmat* Qn, const dmat* C, const dmat* Rn){
	real diff=1, diff2=1, lastdiff=INFINITY;
	dmat* P2=dnew(NX(A), NY(A)); daddI(P2, P(Qn,0));//Initialize P to identity
	dmat* AP=0, * CP=0, * P=0, * CPAt=0, * CCinv=0, * APCt=0, * tmp=0;
	int count=0;
	real thres=1e-14;
	const int maxcount=10000;
	while(diff>thres&&diff2>thres&&count++<maxcount){
		//Compute Sigma_infty aka P
		//P=(A*P*A')-(A*P*C')*(C*P*C'+Rn)^-1*(C*P*A')+Qn
		//Notice that P may not be symmetric due to round off errors.
		dcp(&P, P2);//same old result
		dmm(&AP, 0, A, P, "nn", 1);
		dmm(&CP, 0, C, P, "nn", 1);
		dmm(&CPAt, 0, CP, A, "nt", 1);
		dmm(&APCt, 0, AP, C, "nt", 1);
		dmm(&CCinv, 0, CP, C, "nt", 1);	
		dadd(&CCinv, 1, Rn, 1);	 
		dsvd_pow(CCinv, -1);
		dmm(&P2, 0, AP, A, "nt", 1);
		dmm(&tmp, 0, CCinv, CPAt, "nn", 1);
		dmm(&P2, 1, APCt, tmp, "nn", -1);
		dadd(&P2, 1, Qn, 1);
		diff=dsumsq(P);
		dadd(&P, 1, P2, -1);
		diff=sqrt(dsumsq(P)/diff);
		diff2=fabs(diff-lastdiff);
		lastdiff=diff;
	}
	if(count>=maxcount){
		warning_once("recatti: count=%d, diff=%g, diff2=%g, thres=%g\n", count, diff, diff2, thres);
	}
	dmat* Kinf=0;//Kalman gain K_\infty
	dmm(&Kinf, 0, CP, CCinv, "tn", 1);
	if(Pout) dcp(Pout, P2);
	dfree(AP);dfree(CP); dfree(P); dfree(P2);
	dfree(CPAt); dfree(CCinv); dfree(APCt); dfree(tmp);
	return Kinf;
}

/**
 * @brief Compute Q=(1/dT)*(1-exp(-Ac*dT))*Ac^-1 : averages the state over a sampling period.
 * 
 * @param pQ 	Return the Result
 * @param Ac 	The continuous domain state evolution matrix
 * @param AcI 	Ac^-1
 * @param dT 	The sampling period
 */
static void calc_Q(dmat **pQ, dmat *Ac, dmat *AcI, real dT){
	dmat* tmp=0;
	dexpm(&tmp, 0, Ac, -dT);
	dscale(tmp, -1);
	daddI(tmp, 1);
	dmm(pQ, 0, tmp, AcI, "nn", 1./dT);
	dfree(tmp);
}

/**
 * @brief Compute 3 matrix multipliation: out=out*alpha+beta*op(A)*op(B)*op(C)
 */
static void dmm3(dmat **pout, real alpha, const dmat *A, const dmat *B, const dmat *C, const char trans[4], real beta){
	dmat *tmp=0;
	dmm(&tmp, 0, A, B, trans, 1);
	char trans2[3]; trans2[0]='n'; trans2[1]=trans[2]; trans2[2]=0;
	dmm(pout, alpha, tmp, C, trans2, beta);
	dfree(tmp);
}
/**
 * @brief Compute the discrete state noise covariance from continuous state noise covariance
 * Sigma_kappa: Qn=\int_ti exp(Ac*ti)*Sigma_ep*exp(Ac'*ti) 
 * 
 * @param Ac 		The continuous domain state evolution matrix
 * @param Sigma_ep 	The continuous state noise covariance
 * @param dT 		The sampling period
 * @return dmat* 
 */
static dmat* calc_Qn(dmat *Ac, const dmat *Sigma_ep, real dT, int nsec){
	dmat *expAj=0;
	dmat *Qn=0;//Sigma_kappa
	real dT2=dT/nsec;
	for(int i=0; i<nsec; i++){
		real ti=(i+0.5)*dT2;
		dexpm(&expAj, 0, Ac, ti);//exp(Ac*ti)
		dmm3(&Qn, 1, expAj, Sigma_ep, expAj, "nnt", dT2);
	}
	dfree(expAj);
	return Qn;
}
/**
 * @brief Fine unique dtrats and sort from small to large
 * 
 * @param dtrat_wfs 	dtrat for WFS
 * @return lmat* 
 */
static lmat *unique_dtrats(const lmat *dtrat_wfs){
	int nwfs=NX(dtrat_wfs);
	lmat *dtrats=lnew(nwfs, 1);//unique dtrats
	int ndtrat=0;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int found=0;
		if(!P(dtrat_wfs,iwfs)){
			error("dtrat_wfs[%d]=0\n", iwfs);
		}
		for(int jdtrat=0; jdtrat<ndtrat; jdtrat++){
			if(P(dtrats,jdtrat)==P(dtrat_wfs,iwfs)){
				found=1; break;
			}
		}
		if(!found){
			P(dtrats,ndtrat)=P(dtrat_wfs,iwfs);
			ndtrat++;
		}
	}
	lresize(dtrats, ndtrat, 1);
	lsort(dtrats, 1);//from fast to slow
	return dtrats;
}
/**
 * @brief Kalman filter based on SDE model
 * 
 * @param coeff 	SDE coefficients
 * @param dthi 		High order loop frequency
 * @param dtrat_wfs 	WFS dtrat over high order loop. Each wfs may have a different rate. But there should be only 2 rates
 * @param Gwfs 		WFS measurement from modes. Can be identity if reconstructor is already applied.
 * @param sanea 		WFS measurement noise covariance for each dtrat.
 * @param Proj 		Project modes in statespace to DM/correction space. Can be identity or NULL.
 * @return kalman_t* 
 */
kalman_t* sde_kalman(const dmat *coeff, const real dthi, const lmat* dtrat_wfs, 
	const dcell *Gwfs, const dcell *sanea, const dmat *Proj){
	int nblock=NY(coeff);
	int order=coeff->nx-1;
	/*Ac is block diagonal matrix for continuous domain state evolution.
	  For third order, each block is
	  |-c1 -c2 -c3 |
	  |1    0  0  |
	  |0    1  0  |
	*/

	int nmod=nblock*order;/*number of modes in state*/
	dmat* Ac=dnew(nmod, nmod);//Continous time state propagation vector
	dmat* Sigma_ep=dnew(nmod, nmod);//Continous time state noise 
	dmat* Pd0=dnew(nblock, nmod);//Project from state vector to state.
	for(int iblock=0; iblock<nblock; iblock++){
		real* pcoeff=P(coeff)+(order+1)*iblock;
		real* Aci=P(Ac)+iblock*order*(nmod+1);
		real* si=P(Sigma_ep)+iblock*order*(nmod+1);
		for(int i=0; i<order; i++){
			Aci[i*nmod]=-pcoeff[i]; /*evolution*/
			if(i+1<order){
				Aci[i*(nmod+1)+1]=1;/*keep*/
			}
		}
		si[0]=pow(pcoeff[order], 2);
		P(Pd0, iblock, (iblock+1)*order-1)=1;
	}
	dmat* Pd=0;
	if(Proj){
		dmm(&Pd, 0, Proj, Pd0, "nn", 1);
	} else{
		Pd=dref(Pd0);
	}
	dfree(Pd0);
	dmat* AcI=ddup(Ac);
	dsvd_pow(AcI, -1);//Inverse of Ac
	const int nwfs=NX(Gwfs);
	if(NX(dtrat_wfs)!=nwfs||NY(dtrat_wfs)!=1){
		error("dtrat_wfs should have size %ldx1\n", NX(dtrat_wfs));
	}
	kalman_t* res=mycalloc(1, kalman_t);
	res->dtrats=unique_dtrats(dtrat_wfs);
	int ndtrat=PN(res->dtrats); //(1<<(NX(dtrat_wfs)))-1;
	if(ndtrat>2){
		error("More than two rates are not yet implemented.\n");
	}else if(P(res->dtrats, ndtrat-1)%P(res->dtrats,0)!=0){
		error("The two dtrats must be multiple of each other\n");
	}
	const int dtrat_min=P(res->dtrats,0);
	//const int dtrat_max=P(res->dtrats, ndtrat-1);
	
	res->dthi=dthi;
	res->dtrat_wfs=ldup(dtrat_wfs);
	res->Gwfs=dcelldup(Gwfs);//Used for simulation. Do not change.
	res->sanea=dcelldup(sanea);//WFS measurement noise due to photon and read out noise
	res->Qn=calc_Qn(Ac, Sigma_ep, dthi*dtrat_min, dtrat_min*10);//for fast rate
	
	dexpm(&res->AdM, 0, Ac, dthi); /*discrete state propagation at dthi*/
	//Compute BM=Pd*[1-exp(-Ac*dthi)]*Ac^-1/dthi.
	dmat* QM=NULL;
	calc_Q(&QM, Ac, AcI, dthi);
	dmm(&res->BM, 0, Pd, QM, "nn", 1);//BM=Pd*[1-exp(-Ac*dthi)]*Ac^-1/dthi
	dfree(QM);
	
	res->Rn=dcellnew(ndtrat, 1);//Total WFS noise covariance: state noise plus measurement noise. Sigma_eta
	res->Ad=dcellnew(ndtrat, 1);//State to WFS measurement discrete.
	res->Cd=dcellnew(ndtrat, 1);//State to WFS measurement discrete.
	res->M=dcellnew(ndtrat, 1);
	res->P=dcellnew(ndtrat, 1);//State residual
	res->Rlsq=dccellnew(ndtrat, 1); //for least square reconstruction (testing)
	/*Loop over first set of steps to find needed kalman filter.*/
	for(int idtrat=0; idtrat<ndtrat; idtrat++){
		int dtrat=P(res->dtrats, idtrat);
		int indk=0;
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			if(dtrat%P(dtrat_wfs,iwfs)==0){
				indk|=1<<iwfs;/*this is how we compute the index into kalman*/
			}
		}
		if(P(res->Cd,idtrat)) continue;//case already done
		dmat* Sigma_zeta=0;//covariance of state noise integrated over WFS exposure time
		dmat* Radd=0;
		dmat* Qwfs=0;//State averaging vector (integration over the WFS expoure time): Sigma_kappa

		real dT=dthi*dtrat;
		int nsec=dtrat*10;
		real dT2=dT/nsec;
		dmat* expAjn=0;
		dmat* TQj=0;
		//Compute discrete time state and measurement noise covariance matrix using integral. 
		//This block takes most of the preparation step
		for(int i=0; i<nsec; i++){
			real ti=(i+0.5)*dT2;
			//state noise in WFS: Sigma_zeta=\int_ti (1/T^2)*(1-exp(-Ac*ti))*Ac^-1*Sigma_ep*Ac^-T*(1-exp(-Ac'*ti))
			dexpm(&expAjn, 0, Ac, -ti);	dscale(expAjn, -1); daddI(expAjn, 1); //1-exp(-Ac*ti)
			dmm(&TQj, 0, expAjn, AcI, "nn", 1); //[1-exp(-Ac*ti)]*Ac^-1
			dmm3(&Sigma_zeta, 1, TQj, Sigma_ep, TQj, "nnt", dT2/(dT*dT));
		}
		dfree(TQj);dfree(expAjn);
		dmm3(&Radd, 1, Pd, Sigma_zeta, Pd, "nnt", 1);//Radd = Pd * Sigma_zeta * Pd': state noise in WFS
		calc_Q(&Qwfs, Ac, AcI, dT);//Qwfs = (1-exp(-Ac*dT))*Ac^-1 / dT: state averaging vector. 

		/*compute Rn=Gwfs*Radd*Gwfs'+sanea */
		dcell* Rn=dcellnew(nwfs, nwfs);//WFS measurement noise (photon/rne + state noise)
		dcell* Cd=dcellnew(nwfs, 1);//State to WFS measurement discrete.
		dexpm(&P(res->Ad, idtrat), 0, Ac, dthi*dtrat);  /*discrete state propagation at dT*/
		dcell *GwfsU=dcellnew(nwfs, 1);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			dmat* Gwfsi=0;//Use for reconstruction. 
			if(indk==1&&nmod==12&&KALMAN_IN_SKYC){
				/*Special case for sky coverage with single WFS, make 3rd column zero*/
				Gwfsi=ddup(P(Gwfs,iwfs));
				memset(PCOL(Gwfsi, 2), 0, NX(Gwfsi)*sizeof(real));
			} else{
				Gwfsi=dref(P(Gwfs,iwfs));
			}
			P(Cd,iwfs)=dnew(NX(Gwfsi), NY(Ac));
			P(Rn, iwfs, iwfs)=dnew(NX(Gwfsi), NX(Gwfsi));
			if(dtrat%P(dtrat_wfs,iwfs)==0){
				dmm3(&P(Cd,iwfs), 0, Gwfsi, Pd, Qwfs, "nnn", 1);//Cd=Gwfs*Pd*Qwfs
				dmm3(&P(Rn, iwfs, iwfs), 0, Gwfsi, Radd, Gwfsi, "nnt", 1);//Rn=Gwfsi*Radd*Gwfsi'
				dadd(&P(Rn, iwfs, iwfs), 1, P(sanea, iwfs, iwfs), sqrt((real)P(dtrat_wfs, iwfs)/dtrat));
				P(GwfsU, iwfs)=Gwfsi; 
			}
			//dfree(Gwfsi);
		}
		warning_once("Todo: make sanea per dtrat\n");
		P(res->Rn,idtrat)=dcell2m(Rn); dcellfree(Rn);
		P(res->Cd,idtrat)=dcell2m(Cd); dcellfree(Cd);
		P(res->Rlsq, idtrat)=dcellpinv(GwfsU, dcellsum(res->sanea)>0?res->sanea:0);
		dcellfree(GwfsU);
		dmat *Sigma_infty=0;//Sigma_infty, estimation error covariance. convert to error in modes.
		P(res->M, idtrat)=reccati(&Sigma_infty, P(res->Ad,idtrat), res->Qn, P(res->Cd, idtrat), P(res->Rn, idtrat));
		dmm3(&P(res->P, idtrat), 0, Pd, Sigma_infty, Pd, "nnt", 1);
		dfree(Sigma_infty);
		dfree(Qwfs);
		dfree(Sigma_zeta);
		dfree(Radd);
	}

	dfree(Ac);
	dfree(AcI);
	dfree(Pd);
	dfree(Sigma_ep);
	
	return res;
}
/**free the struct*/
void kalman_free(kalman_t* kalman){
	if(!kalman) return;
	dfree(kalman->AdM);
	dfree(kalman->BM);
	dfree(kalman->Qn);
	lfree(kalman->dtrat_wfs);
	lfree(kalman->dtrats);

	cellfree(kalman->Ad);
	cellfree(kalman->Cd);
	cellfree(kalman->Gwfs);
	cellfree(kalman->sanea);
	cellfree(kalman->Rn);
	cellfree(kalman->M);
	cellfree(kalman->P);
	cellfree(kalman->Rlsq);

	cellfree(kalman->xhat);
	cellfree(kalman->xhat2);
	cellfree(kalman->xhat3);
	cellfree(kalman->psol);
	free(kalman);
}
/**Initialize kalman filter state*/
void kalman_init(kalman_t* kalman){
	if(!kalman->xhat){
		long nx=NX(kalman->AdM);
		kalman->xhat=dcellnew_same(NX(kalman->M),1, nx, 1);
		kalman->xhat2=dcellnew_same(NX(kalman->M),1, nx, 1);
		kalman->xhat3=dnew(nx, 1);
	} else{
		dcellzero(kalman->xhat);
		dcellzero(kalman->xhat2);
		dzero(kalman->xhat3);
		dcellfree(kalman->psol);
	}
}
/**
   Update state vector when there is a measurement. It modifies meas.
	
	xhat3=xhat+M*(meas-Cd*xhat) //xhat is x_k|k-1. xhat3 is x_k|k
	xhat =Ad*xhat3; //xhat is now x_k+1|k
	xhat2=AdM*xhat3; //xhat2 is x_k+1/M|k
	
 */
void kalman_update(kalman_t* kalman, dcell* meas, int idtrat){
	//First, compute PSOL commands
	if(kalman->psol) dcellmm(&meas, kalman->Gwfs, kalman->psol, "nn", 1);
	dmat *measv=meas->m;
	if(!measv){
		warning_once("meas->m is not set\n");
		measv=dcell2m(meas);
	}
	/*difference between real and predicted measurement*/
	dmm(&measv, 1, P(kalman->Cd,idtrat), P(kalman->xhat,idtrat), "nn", -1);
	/*updated estimate*/
	dcp(&kalman->xhat3, P(kalman->xhat,idtrat));
	dmm(&kalman->xhat3, 1, P(kalman->M,idtrat), measv, "nn", 1);
	dmm(&P(kalman->xhat,idtrat), 0, P(kalman->Ad,idtrat), kalman->xhat3, "nn", 1);//xhat is x_k+1|k
	dmm(&P(kalman->xhat2,idtrat), 0, kalman->AdM, kalman->xhat3, "nn", 1);//xhat2 is x_k+1/M|k
	if(!meas->m){
		dfree(measv);
	}
}
/**
   Output correction
   xhat2=AdM*xhat2 //evolve correction forward
   out=out*alpha+(BM*xhat2-psol)*beta
*/
void kalman_output(kalman_t* kalman, dmat** out, real alpha, real beta, int idtrat){
	dcp(&kalman->xhat3, P(kalman->xhat2, idtrat));
	dmm(&P(kalman->xhat2, idtrat), 0, kalman->AdM, kalman->xhat3, "nn", 1);//xhat2 is x_k+2/M for u_k+1/m
	if(idtrat==0&&(alpha==1||beta==1)&&*out){//if we are running integrator, save PSOL.
		if(!kalman->psol) kalman->psol=dcellnew(1,1);
		dcp(&P(kalman->psol, 0), *out);
	}
	dmm(out, alpha, kalman->BM, P(kalman->xhat2, idtrat), "nn", beta);
	if(alpha==1){
		if(kalman->psol) dadd(out, 1, P(kalman->psol, 0), -beta);//somehow this is better than use previous psol.
	}
}
/**
   Test the performance of kalman filter (LQG controller). Derived from servo_test()
   Gwfs2 is a different gradient interaction matrix to model the aliasing in system
*/
dmat* kalman_test(kalman_t* kalman, const dcell *goff, const dmat* input, int flag){
	if(NX(input)>1 && NY(input)==1){/*single mode. each column is for a mode.*/
		reshape((dmat*)input, 1, NX(input));
	}
	const dcell* Gwfs=kalman->Gwfs;
	const int nwfs=NX(Gwfs);
	dcell* rmsn=dcellnew(nwfs, 1);
	dcell* noise=dcellnew(nwfs, 1);
	long ngs[nwfs];
	int dtrat_slow=P(kalman->dtrat_wfs, 0);
	int dtrat_fast=P(kalman->dtrat_wfs, 0);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ng=P(Gwfs,iwfs)->nx;
		if(dsumsq(P(kalman->sanea, iwfs, iwfs))>0){
			P(rmsn,iwfs)=dchol(P(kalman->sanea, iwfs, iwfs));
			P(noise,iwfs)=dnew(ng, 1);
		}
		ngs[iwfs]=ng;
		dtrat_slow=MAX(dtrat_slow, P(kalman->dtrat_wfs, iwfs));
		dtrat_fast=MIN(dtrat_fast, P(kalman->dtrat_wfs, iwfs));
	}
	dcell* acc=dcellnew3(nwfs, 1, ngs, 0);//measurement accumulation
	dcell* meas_fast=dcellnew3(nwfs, 1, ngs, 0);//measurement at fast rate
	dcell* meas_slow=dcellnew3(nwfs, 1, ngs, 0);//measurement at slow rate
	//int nmod=NX(input);
	dmat* mres=ddup(input);
	kalman_init(kalman);
	rand_t rstat;
	seed_rand(&rstat, 1);
	dcell *outc=dcellnew(1,1); 
	dmat *outi=P(outc,0)=dnew_do(NX(mres), 1, (void*)1, 0);//used to wrap column of output
	dcell *corr=dcellnew_same(1,1,NX(input), 1);//corrector status
	dcell *merr_fast=dcellnew(1,1);//fast rate error signal
	dcell *corr_slow=dcellnew(1,1);//slow rate correction status
	real gain=0.5;
	int ndtrat=PN(kalman->M);
	int has_output=0;
	int int_fast=flag==1;//fast loop is integrator
	int int_slow=flag>0;//slow loop is integrator
	//notice the 2-step delay between WFS and DM.
	for(int istep=0; istep<NY(input); istep++){
		outi->p=PCOL(mres,istep);
		dadd(&outi, 1, P(corr,0), -1);//residual
		dcellmm(&acc, Gwfs, outc, "nn", 1);/*CL measurement*/
		if(goff){
			dcelladd(&acc, 1, goff, 1);
		}
		if(!int_fast){
			//handles PSOL saving and usage automatically in LQG integrator mode.
			kalman_output(kalman, &P(corr, 0), 1, gain, 0);//always output fast loop. PSOL estimate
		}else if(has_output){
			dadd(&P(corr, 0), 1, P(merr_fast, 0), gain); 
			has_output=0; //for next step
		}

		int indk=0;

		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int dtrat_wfs=P(kalman->dtrat_wfs,iwfs);
			if((1+istep)%dtrat_wfs==0){//There is measurement output available
				indk|=1<<iwfs;//mark available WFS
				//Average the measurement and zer out the averager.
				dadd(&P(meas_fast,iwfs), 0, P(acc,iwfs), 1./dtrat_wfs); dzero(P(acc,iwfs));
				if(P(rmsn,iwfs)){//Add noise
					drandn(P(noise,iwfs), 1, &rstat);
					if(P(rmsn,iwfs)->nx>1){
						dmm(&P(meas_fast,iwfs), 1, P(rmsn,iwfs), P(noise,iwfs), "nn", 1);
					} else{
						dadd(&P(meas_fast,iwfs), 1, P(noise,iwfs), P(P(rmsn,iwfs),0));
					}
				}
				if(ndtrat>1){
					dadd(&P(meas_slow, iwfs), 1, P(meas_fast, iwfs), (float)dtrat_wfs/dtrat_slow);//accumulate to the slow rate
				}
			}
		}

		if(indk){
			int idtrat=(indk+1)==(1<<nwfs)?(ndtrat-1):0;//0: fast rate, 1: slow rate
			{//fast rate loop always have output
				if(ndtrat>1 && P(corr_slow,0)){//use slow as offset to fast loop
					for(int iwfs=0; iwfs<nwfs; iwfs++){
						if(dtrat_fast==P(kalman->dtrat_wfs,iwfs)){
							dcellmm(&P(meas_fast, iwfs), P(Gwfs, iwfs), P(corr_slow, 0), "nn", 1);
						}
					}
				}
				if(int_fast){
					dcellzero(merr_fast);
					dcellmm(&merr_fast, P(kalman->Rlsq, 0), meas_fast, "nn", 1);
					has_output=1;
				}else{
					kalman_update(kalman, meas_fast, 0);
				}
			}
			if(idtrat!=0){//slow loop has output
				if(int_slow){
					dcellmm(&corr_slow, P(kalman->Rlsq, idtrat), meas_slow, "nn", 0.5);//accumulate
				}else{
					kalman_update(kalman, meas_slow, idtrat);//slow loop has update
					kalman_output(kalman, &P(corr_slow, 0), 1, 0.5, idtrat);
				}
				dzero(meas_slow->m);
			}
		}
	}
	dcellfree(rmsn);
	dcellfree(meas_fast);
	dcellfree(acc);
	dcellfree(noise);
	dfree(outi);
	cellfree(merr_fast);
	cellfree(corr_slow);
	return mres;
}
dmat* kalman_test2(const dmat* coeff, /**<SDE coefficients*/
	const real dthi, /**<Loop frequency*/
	const lmat* dtrat_wfs,   /**<WFS frequency as a fraction of loop*/
	const dcell* Gwfs,  /**<WFS measurement from modes. Can be identity*/
	const dcell* sanea,  /**<WFS measurement noise covariance*/
	const dmat* Proj,
	const dcell *Gwfs2, 
	const dmat* input,
	int flag
){
	kalman_t *kalman=sde_kalman(coeff, dthi, dtrat_wfs, Gwfs, sanea, Proj);
	dmat* res=kalman_test(kalman, Gwfs2, input, flag);
	kalman_free(kalman);
	return res;
}
/**
   Save kalman_t to file
*/
void kalman_write(kalman_t* kalman, const char* format, ...){
	format2fn;
	file_t* fp=zfopen(fn, "wb");
	if(kalman){
		header_t header={MCC_ANY, 12, 1, (char*)"type=struct"};
		write_header(&header, fp);
		char* tmp;
#define WRITE_KEY(fp, str, key)			\
			tmp=str->key->keywords;			\
			str->key->keywords=(char*)#key;		\
			writedata(fp, (cell*)str->key, 0);	\
			str->key->keywords=tmp;

		WRITE_KEY(fp, kalman, AdM);
		WRITE_KEY(fp, kalman, BM);
		WRITE_KEY(fp, kalman, Qn);
		WRITE_KEY(fp, kalman, dtrat_wfs);
		WRITE_KEY(fp, kalman, dtrats);
		WRITE_KEY(fp, kalman, Ad);
		WRITE_KEY(fp, kalman, Cd);
		WRITE_KEY(fp, kalman, Gwfs);
		WRITE_KEY(fp, kalman, sanea);
		WRITE_KEY(fp, kalman, Rn);
		WRITE_KEY(fp, kalman, M);
		WRITE_KEY(fp, kalman, P);
		writearr(fp, 0, sizeof(real), M_REAL, "dthi", &kalman->dthi, 1, 1);
	}
	zfclose(fp);
}
