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
#include "servo.h"
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
			real c1sq=pow(coeff[0+im*3], 2);//c1 is 2*zeta*omega0: damping
			real c2=coeff[1+im*3]; 			//c2 is omega0^2: resonance
			real sigma2=pow(coeff[2+im*3], 2);//strength
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
			for(long i=1; i<data->ncov; i++){
				real val1=P(data->psdcov_in,i);
				real val2=P(data->psdcov_sde,i);
				diff+=pow(val1-val2, 2);
			}
			diff=diff/((abs2(P(data->psdcov_in,0)+P(data->psdcov_sde,0)))*data->ncov);
#elif METHOD==2
			real scale=P(data->psdcov_in,0)/P(data->psdcov_sde,0);
			for(long i=1; i<data->ncov; i++){
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
   Fit a PSD with SDE model. Use initial guess from coeff and return final
   answer in the same format. Functions returns difference in PSD to indicate residual
*/
static real sde_fit_do(dmat **pcoeff, const dmat* psdin, real tmax_fit, int print){
	if(!pcoeff) error("pcoeff must not be NULL\n");
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
	if(!*pcoeff){
		*pcoeff=dnew(3,1);
	}
	if(!P(*pcoeff, 0)||!P(*pcoeff, 1)){
		P(*pcoeff, 0)=1;
		P(*pcoeff, 1)=1;
	}
	dmat *coeff=*pcoeff;
	int ncoeff=NX(coeff);
	int nmod=NY(coeff);
	//Scale to make sure total energy is preserved.
	sde_scale_coeff(coeff, var_in, psdcov_sde, freq, ncov);

	sde_fit_t data={df, freq, psdcov_in, psdcov_sde, norm_in, 0, ncoeff, nmod, ncov};
	//real diff0=sde_diff(P(coeff), &data);
	real tol=1e-8;
	int nmax=2000;
	dminsearch(P(coeff), ncoeff*nmod, tol, nmax, sde_diff, &data);
	//Do not scale coeff after the solution.
	real diff1=sde_diff(P(coeff), &data);
	if(print){
		real zeta=P(coeff,0)/sqrt(P(coeff,1))/2;
		real f0=P(coeff,0)/(2*zeta*2*M_PI);
		info("sde_fit_auto: tmax=%4.2f, f0=%4.1f, zeta=%4.2f, diff=%.2e.\n", tmax_fit, f0, zeta, diff1);
	}
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
	return diff1;
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
real sde_fit(dmat** pcoeff, const dmat* psdin, real tmax_fit, int vibid){
	if(!pcoeff) error("pcoeff must not be NULL\n");
	if(!vibid){
		return sde_fit_do(pcoeff, psdin, tmax_fit, 1);
	} else{
		//Do vibration identification
		dmat* vibs=vibid?psd_vibid(psdin):NULL;
		dcell* coeffs=dcellnew(1, vibs?(1+NY(vibs)):1);
		
		dmat* psd2=ddup(psdin);
		if(vibs&&NY(vibs)>0){
			info("\nnvib=%ld\n", NY(vibs));
			for(int ivib=0; ivib<NY(vibs); ivib++){
				real omegai=2*M_PI*P(vibs,ivib*4+0);
				int i1=P(vibs,ivib*4+2);
				int i2=P(vibs,ivib*4+3)+1;
				if(i2-i1<3){
					i2++;
					i1--;
				}
				dmat* psdi=dsub(psdin, i1, i2-i1, 0, 0);//extract peak
				dmat* coeffi=P(coeffs,ivib)=dnew(3, 1);

				P(coeffi,0)=2*omegai*0.1;
				P(coeffi,1)=omegai*omegai;
				real diff=sde_fit_do(&coeffi, psdi, tmax_fit, 1);
				if(diff>1e-4 || fabs(sqrt(P(coeffi,1))-sqrt(P(P(coeffs,ivib),1)))>2*M_PI){
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
				dfree(psdi);
				real* psd2p=PCOL(psd2, 1);
				for(int i=i1; i<i2; i++){//replace the peak by a linear interpolation
					psd2p[i]=(psd2p[i1-1]*(i2+1-i)+psd2p[i2+1]*(i-(i1-1)))/(i2-i1+2);
				}
			}
			dfree(vibs);
		}
		dmat* coeffi=P(coeffs, coeffs->ny-1)=dnew(3, 1);
		real diff=sde_fit_do(&coeffi, psd2, tmax_fit, 1);
		
		dfree(psd2);
		*pcoeff=dcell2m(coeffs);
		dcellfree(coeffs);
		return diff;
	}
}
/**
 * @brief Fit PSD to SHO by maximum the damping ratio while keeping it less than
 * 5. Does not support vibration identification which have different damping
 * characteristics.
 *
 * @param psdin     The input PSD
 * @return dmat*    The SHO coefficients
 */
real sde_fit_auto(dmat **pcoeff, const_anyarray psdin_, real tfit){
	if(!iscell(psdin_.dc) && NY(psdin_.dm)==2){//single PSD
		dmat *psdin=dmat_cast(psdin_.dm);
		real fmax=P(psdin, NX(psdin)-1, 0);
		if(!tfit) tfit=0.2;
		real zeta_best=0;//record the best zeta
		real diff_best=0;
		real tfit_best=0;
		dmat *coeff=dnew(3,1);
		//maximize the damping while keeping it below 5.
		for(int i=0; i<10; i++){
			P(coeff, 0)=1;
			P(coeff, 1)=1;
			real diff=sde_fit_do(&coeff, psdin, tfit, 0);
			real zeta=P(coeff,0)/sqrt(P(coeff,1))/2;
			real f0=sqrt(P(coeff, 1))/2/M_PI;
			if(zeta>zeta_best && zeta<2 && f0<fmax){
				zeta_best=zeta;
				diff_best=diff;
				tfit_best=tfit;
				dcp(pcoeff, coeff);
			}
			//info("sde_fit_auto: tmax=%4.2f, f0=%4.2f, zeta=%4.2f, diff=%g (%d).\n", tfit, f0, zeta, diff, tfit==tfit_best);
			tfit*=2;
		}
		dfree(coeff);
		if(!zeta_best){
			error("Fitting failed\n");
		}
		real f0=sqrt(P(*pcoeff, 1))/2/M_PI;
		info("sde_fit_auto: tmax=%4.2f, f0=%4.2f, zeta=%4.2f, diff=%.2e.\n", tfit_best, f0, zeta_best, diff_best);
		return diff_best;
	}else{
		int npsd=0;
		if(iscell(psdin_.dc)){//cell array of PSDs
			npsd=PN(psdin_.dc);
		}else if(NY(psdin_.dm)>2){//each column beyond first one is a PSD.
			npsd=NY(psdin_.dm)-1;
		}else{
			error("Invalid input for.\n");
		}
		dcell* coeff=dcellnew(1, npsd);
		real diff_worst=0;
		dmat* psd=NULL;
		for(int ipsd=0; ipsd<npsd; ipsd++){
			if(iscell(psdin_.dc)){
				psd=P(psdin_.dc, ipsd);
			}else{
				if(!psd){
					psd=dnew(NX(psdin_.dm), 2);
					memcpy(PCOL(psd,0), PCOL(psdin_.dm, 0), NX(psdin_.dm)*sizeof(real));
				}
				memcpy(PCOL(psd,1), PCOL(psdin_.dm, 1+ipsd), NX(psdin_.dm)*sizeof(real));
			}
			real diff=sde_fit_auto(&P(coeff, ipsd), psd, tfit);
			if(diff>diff_worst) diff_worst=diff;
		}
		if(!iscell(psdin_.dc)){
			dfree(psd);
		}
		*pcoeff=dcell2m(coeff);
		dcellfree(coeff);
		return diff_worst;
	}
}

/**
 * @brief Compute the reccati equation. When noise (Rn) is too large, the iteration will fail to converge.
 * @param Kinf  [Out] Asymptotic Kalman gain K_\infty
 * @param Pout  [Out] The Sigma_infty: Estimation error covariance matrix
 * @param A 	Block diagonal stage evolution matrix
 * @param Qn 	Discrete state noise covariance
 * @param C 	Measurement interaction matrix
 * @param Rn 	Measurement error covariance matrix
 * @return diff	Converegence indicator
 */
real reccati(dmat**Kinf, dmat** Pout, const dmat* A, const dmat* Qn, const dmat* C, const dmat* Rn){
	real diff=100;
	dmat* P2=dnew(NX(A), NY(A)); daddI(P2, P(Qn,0));//Initialize P to identity
	dmat* AP=0, * CP=0, * P=0, * CPAt=0, * CCinv=0, * APCt=0, * tmp=0;
	int count=0;
	real thres=1e-14;//threshold for the diff: ||P-P2||/||P||
	const int maxcount=10000;
	while(diff>thres&&count++<maxcount){
		//Compute Sigma_infty aka P
		//P=(A*P*A')-(A*P*C')*(C*P*C'+Rn)^-1*(C*P*A')+Qn
		//Notice that P may not be symmetric due to round off errors.
		dcp(&P, P2);//save old result
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
		if(diff>0.1 && count>100){
			dbg("count=%d, diff=%g. not converging, break.\n", count, diff);
			break;
		}
	}
	dmm(Kinf, 0, CP, CCinv, "tn", 1);
	if(Pout) dcp(Pout, P2);
	dfree(AP);dfree(CP); dfree(P); dfree(P2);
	dfree(CPAt); dfree(CCinv); dfree(APCt); dfree(tmp);
	return diff;
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
 * @brief Kalman filter based on SDE model
 * 
 * @param coeff 	SDE coefficients
 * @param dthi 		High order loop frequency
 * @param dtrat_wfs WFS dtrat over high order loop. Each wfs may have a different rate. Only 1 or 2 rates are supported
 * @param mdirect	Mode indices that are directly controlled by the slower loop.
 * @param Gwfs 		WFS measurement from modes. Can be identity if reconstructor is already applied.
 * @param Cnn 		WFS measurement noise covariance in radian^2. nwfs*ndtrat.
 * @param Proj 		Project modes in statespace to DM/correction space. Can be identity or NULL.
 
 * @return kalman_t* 
 */
kalman_t* sde_kalman(const dmat *coeff, const real dthi, const lmat* dtrat_wfs,const lmat *mdirect,
	const dcell *Gwfs, const dcell *Cnn, const dmat *Proj){
	const long nmod=NY(coeff);//number of modes
	const long order=coeff->nx-1;
	/*Ac is block diagonal matrix for continuous domain state evolution.
	  For third order, each block is
	  |-c1 -c2 -c3 |
	  |1    0  0  |
	  |0    1  0  |
	*/

	const long nstate=nmod*order;/*number of states*/
	dmat* Ac=dnew(nstate, nstate);//Continous time state propagation vector
	dmat* Sigma_ep=dnew(nstate, nstate);//Continous time state noise 
	dmat* Pd0=dnew(nmod, nstate);//Project from state vector to state.
	for(int iblock=0; iblock<nmod; iblock++){
		real* pcoeff=P(coeff)+(order+1)*iblock;
		real* Aci=P(Ac)+iblock*order*(nstate+1);
		real* si=P(Sigma_ep)+iblock*order*(nstate+1);
		for(int i=0; i<order; i++){
			Aci[i*nstate]=-pcoeff[i]; /*evolution*/
			if(i+1<order){
				Aci[i*(nstate+1)+1]=1;/*keep*/
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
	const long nwfs=NX(Gwfs);
	if(NX(dtrat_wfs)!=nwfs||NY(dtrat_wfs)!=1){
		error("dtrat_wfs should have size %ldx1\n", NX(dtrat_wfs));
	}
	kalman_t* res=mycalloc(1, kalman_t);
	res->dtrats=lunique(dtrat_wfs, 1);
	const long ndtrat=PN(res->dtrats); //(1<<(NX(dtrat_wfs)))-1;
	if(ndtrat>2){
		error("More than two rates are not yet implemented.\n");
	}else if(P(res->dtrats, ndtrat-1)%P(res->dtrats,0)!=0){
		error("The two dtrats must be multiple of each other\n");
	}
	
	if(Cnn && NY(Cnn)<ndtrat){
		error("Cnn should have at least %ld columns\n", ndtrat);
	}
	if(mdirect && PN(mdirect)!=nmod){
		error("mdirect should have %ld elements.\n", nmod);
	}
	if(PN(Gwfs)!=nwfs){
		error("Gwfs shall have %ld elements\n", nwfs);
	}
	if(Proj && (NX(Proj)!=nmod || NY(Proj)!=nmod)){
		error("Proj shall have dimenion of %ldx%ld\n", nmod, nmod);
	}
	//const int dtrat_min=P(res->dtrats,0);
	//const int dtrat_max=P(res->dtrats, ndtrat-1);
	
	res->dthi=dthi;
	res->dtrat_wfs=ldup(dtrat_wfs);
	res->Gwfs=dcelldup(Gwfs);//Used for simulation. Do not change.
	res->Cnn=dcelldup(Cnn);//WFS measurement noise due to photon and read out noise
	res->mdirect=mdirect?lref(mdirect):NULL;
	dexpm(&res->AdM, 0, Ac, dthi); /*discrete state propagation at dthi*/
	//Compute BM=Pd*[1-exp(-Ac*dthi)]*Ac^-1/dthi.
	dmat* QM=NULL;
	calc_Q(&QM, Ac, AcI, dthi);
	dmm(&res->BM, 0, Pd, QM, "nn", 1);//BM=Pd*[1-exp(-Ac*dthi)]*Ac^-1/dthi
	dfree(QM);
	res->Rn=dcellnew(ndtrat, 1);//Total WFS noise covariance: state noise plus measurement noise. Sigma_eta
	res->Ad=dcellnew(ndtrat, 1);//State to WFS measurement discrete.
	res->Cd=dcellnew(ndtrat, 1);//State to WFS measurement discrete.
	res->M=dcellnew(ndtrat, 1);//Kalman gain
	res->P=dcellnew(ndtrat, 1);//residual
	res->Rlsq=dccellnew(ndtrat, 1); //for least square reconstruction (testing)
	int failed=0;
	/*Loop over first set of steps to find needed kalman filter.*/
	for(int idtrat=0; idtrat<ndtrat; idtrat++){
		int dtrat=P(res->dtrats, idtrat);
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
		dcell *Mproj=NULL;
		if(mdirect){
			Mproj=dcellnew_same(1, 1, nmod, nmod);
			if(idtrat==0){//fast loop
				for(int i=0; i<nmod; i++){
					P(P(Mproj,0),i,i)=P(mdirect,i)?0:1;
				}
			}else{
				dcellmm(&Mproj, P(res->Rlsq, 0), Gwfs, "nn", -1);
				daddI(P(Mproj, 0), 1);
				//writebin(P(res->Rlsq, 0), "Rlsq0");
				//writebin(Gwfs, "Gwfs0");
			}
			//dshow(P(Mproj, 0), "Mproj_%d", idtrat);
		}
		dmm3(&Radd, 1, Pd, Sigma_zeta, Pd, "nnt", 1);//Radd = Pd * Sigma_zeta * Pd': state noise in WFS
		calc_Q(&Qwfs, Ac, AcI, dT);//Qwfs = (1-exp(-Ac*dT))*Ac^-1 / dT: state averaging vector. 
		dmat *Qn=calc_Qn(Ac, Sigma_ep, dthi*dtrat, dtrat*10);
		/*compute Rn=Gwfs*Radd*Gwfs'+Cnn */
		dcell* Rn=dcellnew(nwfs, nwfs);//WFS measurement noise (photon/rne + state noise)
		dcell* Cd=dcellnew(nwfs, 1);//State to WFS measurement discrete.
		dexpm(&P(res->Ad, idtrat), 0, Ac, dthi*dtrat);  /*discrete state propagation at dT*/
		dcell *GwfsU=dcellnew(nwfs, 1);
		dcell *neai=dcellnew(nwfs, nwfs);
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int ng=NX(P(Gwfs, iwfs));
			P(Cd,iwfs)=dnew(ng, nstate);
			P(Rn, iwfs, iwfs)=dnew(ng, ng);
			if(dtrat%P(dtrat_wfs,iwfs)==0){
				dmat* Gwfsi=0;//Use for reconstruction. 
				if(Mproj){
					dmm(&Gwfsi, 0, P(Gwfs, iwfs), P(Mproj,0), "nn", 1);
				} else{
					Gwfsi=dref(P(Gwfs,iwfs));
				}
				P(GwfsU, iwfs)=Gwfsi; 
				dmm3(&P(Cd,iwfs), 0, Gwfsi, Pd, Qwfs, "nnn", 1);//Cd=Gwfs*Pd*Qwfs
				dmm3(&P(Rn, iwfs, iwfs), 0, Gwfsi, Radd, Gwfsi, "nnt", 1);//Rn=Gwfsi*Radd*Gwfsi'
				if(P(Cnn, iwfs, idtrat)){
					dadd(&P(Rn, iwfs, iwfs), 1, P(Cnn, iwfs, idtrat), sqrt((real)P(dtrat_wfs, iwfs)/dtrat));
					P(neai, iwfs, iwfs)=dnew(ng, ng);
					for(int ig=0; ig<ng; ig++){
						P(P(neai, iwfs, iwfs), ig, ig)=1./P(P(Cnn, iwfs, idtrat), ig, ig); //in radian^-2 for weighting
					}
					//dshow(P(neai, iwfs, iwfs), "neai");
				}
			}
		}
		P(res->Rn,idtrat)=dcell2m(Rn); dcellfree(Rn);
		P(res->Cd,idtrat)=dcell2m(Cd); dcellfree(Cd);
		P(res->Rlsq, idtrat)=dcellpinv(GwfsU, neai);
		if(Mproj){
			dcell *tmp=P(res->Rlsq, idtrat);
			P(res->Rlsq,idtrat)=NULL;
			dcellmm(&P(res->Rlsq,idtrat), Mproj, tmp, "nn", 1);
			dcellfree(tmp);
		}
		dcellfree(GwfsU);
		dcellfree(neai);
		dmat *Sigma_infty=0;//Sigma_infty, estimation error covariance. convert to error in modes.
		real diff=reccati(&P(res->M, idtrat), &Sigma_infty, P(res->Ad,idtrat), Qn, P(res->Cd, idtrat), P(res->Rn, idtrat));
		dmm3(&P(res->P, idtrat), 0, Pd, Sigma_infty, Pd, "nnt", 1);
		if(diff>0.1){//not converging. mark failed.
			failed=1;
		}
		//info("trace(P)=%g\n", dtrace(P(res->P, idtrat)));
		dfree(Sigma_infty);
		dfree(Qwfs);
		dfree(Sigma_zeta);
		dfree(Radd);
		dfree(Qn);
		cellfree(Mproj);
	}//for idtrat

	dfree(Ac);
	dfree(AcI);
	dfree(Pd);
	dfree(Sigma_ep);
	if(failed){
		kalman_free(res);
		res=NULL;
	}
	return res;
}
/**free the struct*/
void kalman_free(kalman_t* kalman){
	if(!kalman) return;
	dfree(kalman->AdM);
	dfree(kalman->BM);
	lfree(kalman->dtrat_wfs);
	lfree(kalman->dtrats);
	lfree(kalman->mdirect);

	cellfree(kalman->Ad);
	cellfree(kalman->Cd);
	cellfree(kalman->Gwfs);
	cellfree(kalman->Cnn);
	cellfree(kalman->Rn);
	cellfree(kalman->M);
	cellfree(kalman->P);
	cellfree(kalman->Rlsq);

	cellfree(kalman->xhat);
	cellfree(kalman->xhat2);
	cellfree(kalman->xhat3);
	cellfree(kalman->xout);
	cellfree(kalman->psol);
	free(kalman);
}
/**Initialize kalman filter state*/
void kalman_init(kalman_t* kalman){
	if(!kalman->xhat){
		long nstat=NX(kalman->AdM);
		long ndtrat=NX(kalman->M);
		long nmod=NX(kalman->BM);
		kalman->xhat=dcellnew_same(ndtrat, 1, nstat, 1);
		kalman->xhat2=dcellnew_same(ndtrat, 1, nstat, 1);
		kalman->xhat3=dnew(nstat, 1);
		kalman->xout=dcellnew(2,1);
		kalman->psol=dcellnew_same(ndtrat, 1, nmod, 1);
	} else{
		dcellzero(kalman->xhat);
		dcellzero(kalman->xhat2);
		dzero(kalman->xhat3);
		dcellzero(kalman->xout);
		dcellzero(kalman->psol);
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
	for(int iwfs=0; iwfs<NX(kalman->Gwfs); iwfs++){
		dmm(&P(meas, iwfs), 1, P(kalman->Gwfs, iwfs), P(kalman->psol, idtrat), "nn", 1);
	}
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
	if((alpha==1||beta==1)&&*out){//save PSOL.
		dcp(&P(kalman->psol, idtrat), *out);
	}
	dcp(&kalman->xhat3, P(kalman->xhat2, idtrat));
	dmm(&P(kalman->xhat2, idtrat), 0, kalman->AdM, kalman->xhat3, "nn", 1);//xhat2 is x_k+2/M for u_k+1/m
	dmm(&P(kalman->xout, idtrat), 0, kalman->BM, P(kalman->xhat2, idtrat), "nn", 1);
	if(alpha==1){
		dadd(&P(kalman->xout, idtrat), 1, P(kalman->psol, idtrat), -1);//form closed loop signal
	}
	if(!*out){
		*out=dnew(NX(kalman->xout, idtrat), 1);
	}
	dadd(out, alpha, P(kalman->xout, idtrat), beta);
}
/**
   Test the performance of kalman filter (LQG controller). Derived from servo_test()
   Gwfs2 is a different gradient interaction matrix to model the aliasing in system
   flag controlls the servo type:
   0: both LQG
   1: both integrator
   2: LQG (fast) + integrator (slow)
*/
dmat* kalman_test(kalman_t* kalman, const dcell *goff, const dmat* input, int flag){
	if(NX(input)>1 && NY(input)==1){/*single mode. each column is for a mode.*/
		reshape((dmat*)input, 1, NX(input));
	}
	const dcell* Gwfs=kalman->Gwfs;
	const int nwfs=NX(Gwfs);
	dcell* rmsn=dcellnew(nwfs, 1);//RMS noise
	long ngs[nwfs];
	int dtrat_fast=P(kalman->dtrats, 0);
	int dtrat_slow=P(kalman->dtrats, PN(kalman->dtrats)>1?1:0);
	for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ng=P(Gwfs,iwfs)->nx;
		int idtrat=(P(kalman->dtrat_wfs, iwfs)==dtrat_fast)?0:1;
		if(dsumsq(P(kalman->Cnn, iwfs, idtrat))>0){
			P(rmsn, iwfs)=dnew(ng, 1);
			for(int ig=0; ig<ng; ig++){
				P(P(rmsn,iwfs), ig)=sqrt(P(P(kalman->Cnn, iwfs, idtrat),ig,ig));
			}
		}
		ngs[iwfs]=ng;
	}
	dcell* acc=dcellnew3(nwfs, 1, ngs, 0);//measurement accumulation
	dcell* meas_fast=dcellnew3(nwfs, 1, ngs, 0);//measurement at fast rate
	dcell* meas_slow=dcellnew3(nwfs, 1, ngs, 0);//measurement at slow rate
	const int nmod=NX(input);
	dmat* mres=ddup(input);	//store correction residual
	kalman_init(kalman);
	rand_t rstat;
	seed_rand(&rstat, 1);
	dcell *outc=dcellnew(1,1); 
	dmat *outi=P(outc,0)=dnew_do(NX(mres), 1, (void*)1, 0);//used to wrap column of output
	dcell *merr_fast=dcellnew_same(1,1, nmod, 1);//fast rate error signal
	dcell *merr_slow=dcellnew_same(1,1, nmod, 1);//slow rate error signal
	dmat *mreal=dnew(nmod, 1);//corrector status
	dmat *mreal_offset=dnew(nmod, 1);//offset to fast rate
	dmat *mreal_slow=dnew(nmod, 1);//slow loop output
	real gain=0.5;
	int ndtrat=PN(kalman->M);
	int int_fast=flag==1;//fast loop is integrator
	int int_slow=flag>0;//slow loop is integrator
	dmat *ep=dnew(1,1); P(ep, 0)=gain;
	servo_t *st_fast=int_fast?servo_new(mreal, NULL, 0, kalman->dthi*dtrat_fast, ep):NULL;
	P(ep, 0)=0.1;
	servo_t *st_slow=int_slow?servo_new(mreal_slow, NULL, 0, kalman->dthi*dtrat_slow, ep):NULL;
	dfree(ep);
	//notice the 2-step delay between WFS and DM.
	for(int istep=0; istep<NY(input); istep++){
		outi->p=PCOL(mres,istep);
		dadd(&outi, 1, mreal, -1);//residual
		dcellmm(&acc, Gwfs, outc, "nn", 1);//WFS measurement
		if(goff){
			dcelladd(&acc, 1, goff, 1);//WFS measurement bias
		}
		//Process servo update (2 cycle delay)
		if(st_fast){
			servo_output(st_fast, &mreal);
		}else{
			kalman_output(kalman, &mreal, 1, gain, 0);//output fast loop at every step.
		}
		if(dtrat_slow!=dtrat_fast){
			if(st_slow){
				servo_output(st_slow, &mreal_slow);
			}else{
				kalman_output(kalman, &mreal_slow, 1, 0.5, 1);
			}
			for(int i=0; i<nmod; i++){
				if(kalman->mdirect && P(kalman->mdirect, i)){
					P(mreal, i)=P(mreal_slow, i);
				}else{
					P(mreal_offset, i)=P(mreal_slow, i);
				}
			}
		}
		//Process WFS measurement output.
		int indk=0;
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			int dtrat_wfs=P(kalman->dtrat_wfs,iwfs);
			if((1+istep)%dtrat_wfs==0){//WFS has output
				indk|=1<<iwfs;//mark WFS that has output
				//Average the measurement and zero out the averager.
				dadd(&P(meas_fast,iwfs), 0, P(acc,iwfs), 1./dtrat_wfs); 
				dzero(P(acc,iwfs));
				for(int ig=0; ig<NX(P(meas_fast,iwfs)); ig++){//add measurement noise
					P(P(meas_fast,iwfs), ig)+=randn(&rstat)*PR(P(rmsn, iwfs), ig);
				}
				if(ndtrat>1){//accumulate to the slow rate output
					dadd(&P(meas_slow, iwfs), 1, P(meas_fast, iwfs), (float)dtrat_wfs/dtrat_slow);
				}
			}
		}

		if(indk){
			int idtrat=(indk+1)==(1<<nwfs)?(ndtrat-1):0;//0: fast rate, 1: slow rate
			//fast rate loop always have output
			if(dtrat_slow!=dtrat_fast){//use slow loop as offset to fast loop
				for(int iwfs=0; iwfs<nwfs; iwfs++){
					if(dtrat_fast==P(kalman->dtrat_wfs,iwfs)){//fast WFS
						dcellmm(&P(meas_fast, iwfs), P(Gwfs, iwfs), mreal_offset, "nn", 1);
					}
				}
			}
			if(st_fast){
				dcellzero(merr_fast);
				dcellmm(&merr_fast, P(kalman->Rlsq, 0), meas_fast, "nn", 1);
			}else{
				kalman_update(kalman, meas_fast, 0);
			}
			
			if(idtrat!=0){//slow loop has output
				if(st_slow){
					dcellzero(merr_slow);
					dcellmm(&merr_slow, P(kalman->Rlsq, idtrat), meas_slow, "nn", 1);
					servo_filter(st_slow, P(merr_slow,0));
				}else{
					kalman_update(kalman, meas_slow, idtrat);//slow loop has update
				}
				dzero(meas_slow->m);
			}
		}//if indk
		if(st_fast){
			servo_filter(st_fast, indk?P(merr_fast,0):NULL);
		}		
	}
	dcellfree(rmsn);
	dcellfree(meas_fast);
	dcellfree(acc);
	dfree(outi);
	cellfree(merr_fast);
	cellfree(mreal_slow);
	cellfree(mreal_offset);
	return mres;
}
dmat* kalman_test2(const dmat* coeff, /**<SDE coefficients*/
	const real dthi, /**<Loop frequency*/
	const lmat* dtrat_wfs,   /**<WFS frequency as a fraction of loop*/
	const lmat* mdirect,
	const dcell* Gwfs,  /**<WFS measurement from modes. Can be identity*/
	const dcell* Cnn,  /**<WFS measurement noise covariance in radian^2*/
	const dmat* Proj,
	const dcell *goff, 
	const dmat* input,
	int flag
){
	kalman_t *kalman=sde_kalman(coeff, dthi, dtrat_wfs, mdirect, Gwfs, Cnn, Proj);
	dmat* res=kalman_test(kalman, goff, input, flag);
	kalman_free(kalman);
	return res;
}
dmat* kalman_test3(const dmat* input, int flag, const char* fn){
	kalman_t *kalman=kalman_read("%s", fn);
	dmat* res=kalman_test(kalman, NULL, input, flag);
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
		header_t header={MCC_ANY, 14, 1, (char*)"type=struct"};
		write_header(&header, fp);
		char* tmp;
#define WRITE_KEY(fp, str, key)	if(str->key){\
			tmp=str->key->keywords;		\
			str->key->keywords=(char*)#key;\
			writedata(fp, (cell*)str->key, 0);	\
			str->key->keywords=tmp;\
		}else{\
			header_t header2={MCC_ANY, 0, 0, (char*)#key};\
			write_header(&header2, fp);\
		}

		WRITE_KEY(fp, kalman, AdM);
		WRITE_KEY(fp, kalman, BM);
		
		WRITE_KEY(fp, kalman, dtrat_wfs);
		WRITE_KEY(fp, kalman, dtrats);
		WRITE_KEY(fp, kalman, mdirect);

		WRITE_KEY(fp, kalman, Ad);
		WRITE_KEY(fp, kalman, Cd);
		WRITE_KEY(fp, kalman, Gwfs);
		WRITE_KEY(fp, kalman, Cnn);
		WRITE_KEY(fp, kalman, Rn);
		WRITE_KEY(fp, kalman, M);
		WRITE_KEY(fp, kalman, P);
		WRITE_KEY(fp, kalman, Rlsq);
		writearr(fp, 0, sizeof(real), M_REAL, "dthi", &kalman->dthi, 1, 1);
	}
	zfclose(fp);
}

kalman_t* kalman_read(const char* format, ...){
	format2fn;
	cell *tmp=readbin("%s", fn);
	kalman_t *kalman=mycalloc(1, kalman_t);
	int ic=0;
	kalman->AdM=dmat_cast(tmp->p[ic++]);
	kalman->BM=dmat_cast(tmp->p[ic++]);
	
	kalman->dtrat_wfs=lmat_cast(tmp->p[ic++]);
	kalman->dtrats=lmat_cast(tmp->p[ic++]);
	kalman->mdirect=lmat_cast(tmp->p[ic++]);

	kalman->Ad=dcell_cast(tmp->p[ic++]);
	kalman->Cd=dcell_cast(tmp->p[ic++]);
	kalman->Gwfs=dcell_cast(tmp->p[ic++]);
	kalman->Cnn=dcell_cast(tmp->p[ic++]);
	kalman->Rn=dcell_cast(tmp->p[ic++]);
	kalman->M=dcell_cast(tmp->p[ic++]);
	kalman->P=dcell_cast(tmp->p[ic++]);
	kalman->Rlsq=(dccell*)tmp->p[ic++];
	kalman->dthi=dmat_cast(tmp->p[ic++])->p[0];
	if(ic!=PN(tmp)){
		error("kalman_t: read and write Mismatch\n");
	}
	kalman_init(kalman);
	return kalman;
}
