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
#include <unistd.h>
#include "zernike.h"
#include "accphi.h"
#include "turbulence.h"
#include "mkh.h"
/**
   Generating Zernike Rnm for radial order ir, azimuthal or im.  used by zernike().
   */
static dmat* zernike_Rnm(const dmat* locr, int ir, int im){
	if(ir<0||im < 0||im>ir||(ir-im)%2!=0)
		error("Invalid ir, im (%d, %d)\n", ir, im);
	const long nloc=NX(locr)*NY(locr);
	dmat* Rnm=dnew(NX(locr), NY(locr));
	const int ns=(ir-im)/2+1;
	real coeff[ns];
	int power[ns];
	for(int s=0; s<ns; s++){
		coeff[s]=factorial((ir+im)/2-s+1, ir-s)/factorial(1, s)/factorial(1, (ir-im)/2-s)*pow(-1, s);
		power[s]=ir-2*s;
	}
OMP_FOR(NTHREAD)
	for(long iloc=0; iloc<nloc; iloc++){
		real tmp=0;
		real r=P(locr,iloc);
		for(int s=0; s<ns; s++){
			tmp+=coeff[s]*pow(r, power[s]);
		}
		P(Rnm,iloc)=tmp;
	}
	return Rnm;
}

/**
   Create Zernike modes on loc with diameter D and radial order from rmin to rmax
   r=0 is piston
   r=1 is tip/tilt
   r=2 is quadratic modes
   if nopiston is set, skip piston mode.
   if flag is >0, only radial mode is used.
   if flag is  0, all modes between rmin and rmax.
   if flag is <0, only generate mode -flag. rmin and rmax is irrelevant
   set D to -|D| to avoid clipping points outside of D.
*/
dmat* zernike(const loc_t* loc, real D, int rmin, int rmax, int flag){
	if(!loc) return NULL;
	if(flag<0){
		rmin=ceil((sqrt(8.*(-flag)+1)-3)*0.5);
		rmax=rmin;
	}
	int nr3=(int)floor((sqrt(loc->nloc*8+1)-3)*0.5);
	if(rmax>nr3){
		warning("Reduce rmax=%d to %d\n", rmax, nr3);
		rmax=nr3;
	}
	if(rmin>rmax){
		warning("Invalid rmin=%d, rmax=%d\n", rmin, rmax);
		return NULL;
	}
	const long nloc=loc->nloc;
	int nmod=0;
	if(flag>0){//radial only
		nmod=rmax/2-(rmin+1)/2+1;
	}else if(flag<0){//a specific mode
		nmod=1;
	}else{
		nmod=(rmax+1)*(rmax+2)/2-(rmin)*(rmin+1)/2;
	}
	real D2=loc_diam(loc);
	int limitr=1;//make opd 0 outside of D
	if(D<0){
		D=-D;
		limitr=0;
	}else if(!D||D>D2){
		D=D2;
	}
	if(fabs(D-D2)>(D2)*0.5){
		if(limitr){
			warning("Specified diameter D=%g is smaller than half of loc D=%g\n", D, D2);
		}
	}

	dmat* restrict opd=dnew(nloc, nmod);
	dmat* restrict locr=dnew(nloc, 1);
	dmat* restrict locs=dnew(nloc, 1);
	const real* restrict locx=loc->locx;
	const real* restrict locy=loc->locy;
	const real R1=2./D;
	long nover=0; real rover=1;
	for(long iloc=0; iloc<nloc; iloc++){
		P(locr,iloc)=sqrt(pow(locx[iloc], 2)+pow(locy[iloc], 2))*R1;
		if(P(locr,iloc)>1 && limitr){
			nover++;
			if(P(locr,iloc)>rover){
				rover=P(locr,iloc);
			}
			P(locr,iloc)=1;//2020-03-06: prevent r from above 1.
		}
		P(locs,iloc)=atan2(locy[iloc], locx[iloc]);
	}
	if(rover>1.5){
	//if(nover > (M_PI*D/loc->dx)*1.5 ){
		info("%ld (of %ld) points are outside D=%g circle with maximum radius %gm\n",
			nover, nloc, D, rover*D*0.5);
	}
	int cmod=0;//index into opd
	int imod=0;//Noll's count of modes
	for(int ir=rmin; ir<=rmax; ir++){
		imod=(ir)*(ir+1)/2+1;
		for(int im=0; im<=ir; im++){
			if((ir-im)%2!=0) continue;//invalid combo
			if(flag>0&&im!=0) continue;//we want radial only
			
			dmat* Rnm=zernike_Rnm(locr, ir, im);
			if(im==0){/*Radial*/
				if(flag>=0 || (flag+imod)==0){
					real coeff=sqrt(ir+1.);
					real* restrict pmod=PCOL(opd, cmod); 
					OMP_FOR(NTHREAD)
					for(long iloc=0; iloc<nloc; iloc++){
						pmod[iloc]=P(Rnm,iloc)*coeff;
					}
					cmod++;
				}
				imod++;
			} else if(flag<=0){
				real coeff=sqrt(2*(ir+1.));
				int off1=(imod)%2==1?0:1;
				if(!flag||(imod+off1+flag)==0){
					real* restrict pmods=PCOL(opd, flag?0:(cmod+off1));//odd imod for sin
					OMP_FOR(NTHREAD)
					for(long iloc=0; iloc<nloc; iloc++){
						pmods[iloc]=P(Rnm,iloc)*coeff*sin(im*P(locs,iloc));
					}
				}
				if(!flag||(imod+1-off1+flag)==0){
					real* restrict pmodc=PCOL(opd, flag?0:(cmod+1-off1));//even imod for cos
					OMP_FOR(NTHREAD)
					for(long iloc=0; iloc<nloc; iloc++){
						pmodc[iloc]=P(Rnm,iloc)*coeff*cos(im*P(locs,iloc));
					}
					cmod+=2;
				}
				imod+=2;
			}
			dfree(Rnm);
		}
	}
	dfree(locr);
	dfree(locs);
	return opd;
}
/**
   return zernike index of radial mode ir, and azimuthal mode im. Maximum radial
   order is nr (inclusive)
 */
static lmat* zernike_index(int nr){
	if(nr<0) return 0;
	lmat* out=lnew(nr+1, nr+1);
	for(long ix=0; ix<(nr+1)*(nr+1); ix++){
		P(out,ix)=-1;
	}
	int count=0;
	for(int ir=0; ir<=nr; ir++){
		for(int im=0; im<=ir; im++){
			if((ir-im)%2!=0) continue;
			P(out,im,ir)=count;
			if(im==0){//single, pure radial mode
				count++;
			} else{//real mode as a pair
				count+=2;
			}
		}
	}
	return out;
}
/**
   Covariance of zernike modes in Kolmogorov Turbulence. Only modes with the
   same m have non-zero covariance.

   Based on Eq 3.14 in Adaptive Optics in Astronomy (Roddier 1999).
   Verified against values in J.Y.Wang 1978, table II(a,b).

   Notice that (D/r0)^(5/3) has been factored out.

   2015-06-05: This is not longer being used. It is inaccurate for annular
   pupil. Analytic calculation using the power spectrum is used instead.
*/
dmat* zernike_cov_kolmogorov(int nr){
	int nmod=(nr+1)*(nr+2)/2;
	dmat* res=dnew(nmod, nmod);
	dmat* pres=res;
	lmat* zind=zernike_index(nr);
	for(int ir=0; ir<=nr; ir++){
		for(int im=0; im<=ir; im++){
			if((ir-im)%2!=0) continue;
			long ict0=P(zind,im,ir);
			long icts=0, ictc=0;
			if(ict0%2==1){
				ictc=ict0;//cos term
				icts=ict0+1;
			} else{
				icts=ict0;
				ictc=ict0+1;
			}
			for(int jr=0; jr<=ir; jr++){
				if(ir==0||jr==0) continue;//ignore piston
				long jct0=P(zind,im,jr);
				long jcts=0, jctc=0;
				if((jct0)%2==1){
					jctc=jct0;
					jcts=jct0+1;
				} else{
					jcts=jct0;
					jctc=jct0+1;
				}

				if(jct0==-1){//doesn't have the same azimuthal mode
					continue;
				}
				real tmp=7.2e-3*sqrt((ir+1.)*(jr+1.))*pow(-1, (ir+jr-2*im)*0.5)*pow(M_PI, 8./3.)
					*(tgamma(14./3.)*tgamma((ir+jr-5./3.)*0.5))
					/(tgamma((ir-jr+17./3.)*0.5)*tgamma((jr-ir+17./3.)*0.5)*tgamma((ir+jr+23./3.)*0.5));
				if(im==0){
					P(pres, ict0, jct0)=P(pres, jct0, ict0)=tmp;
				} else{
					P(pres, jctc, ictc)=P(pres, ictc, jctc)=tmp;
					P(pres, jcts, icts)=P(pres, icts, jcts)=tmp;
				}
			}
		}
	}
	//P(pres,0,0)=0; //piston covariance
	lfree(zind);
	return res;
}

/**
   Compute covariance matrix of modes in von Karman turbulence.
   @param loc   The location grid
   @param modz  The modal matrix (zernike or zonal)
   @param L0    The outer scale
   @returns     The mode coefficient covariance matrix

   The turbulence OPD y can be expressed as \f$y=M*x\f$ where \f$x\f$ is the
   modal coefficient and \f$M\f$ is the modal matrix (parameter modz). The OPD
   covariance \f$<yy'>\f$ can be expressed as followes: \f$<yy'>=M*<xx'>*M'\f$
   where \f$<xx'>\f$ is the modal coefficient covariance matrix. We have
   \f$M'<yy'>M=(M'M)<xx'>(M'M)=CC<xx'>CC\f$ where \f$CC=(M'M)\f$. Therefore
   \f$<xx'>=CC^{-1} M'<yy'>M CC^{-1}\f$

   The function first computes \f$DD=M'<yy'>M\f$ with \f$DD[i,j]=\sum_{k,m}
   M[k,i]*yy[k,m]*M[m,j]\f$ which can be expressed in Fourier domain as
   \f$DD=\sum_{k',m'} Mhat[k',i] Mhat[m',j] \Phi[k',m']\f$ where \f$\Phi\f$ is
   the turbulence power spectrum. 

   For modes \f$M\f$ defined on the DM, its influence function should be
   includes. \f$<yy'>=H*M*<xx'>*M'H'\f$. 
*/
dmat* cov_vonkarman(const loc_t* loc, const dmat* modz,	real L0){
	if(!loc|| !modz) return NULL;
	dmat* CC=0;
	dmm(&CC, 1, modz, modz, "tn", 1);
	long nembed=0;
	lmat* embed=loc_create_embed(&nembed, loc, 2, 0);
	int nmod=NY(modz);
	ccell* spect=ccellnew(nmod, 1);
	OMP_FOR(NTHREAD)
	for(long ic=0; ic<nmod; ic++){
		P(spect,ic)=cnew(nembed, nembed);
		for(long ix=0; ix<loc->nloc; ix++){
			P(P(spect,ic),P(embed,ix))=P(modz, ix, ic);
		}
		cfftshift(P(spect,ic));
		cfft2(P(spect,ic), -1);
	}
	dmat* turbspec=turbpsd(nembed, nembed, loc->dx, 0.2, L0, -11./3., 1);
	P(turbspec,0)=0;//remove piston.
	dmat* DD=dnew(nmod, nmod);
	OMP_FOR(NTHREAD)
	for(long ic=0; ic<nmod; ic++){
		for(long id=0; id<=ic; id++){
			real tmp=0;
			for(long ip=0; ip<nembed*nembed; ip++){
				tmp+=creal(P(P(spect,ic),ip)*conj(P(P(spect,id),ip)))*P(turbspec,ip);
			}
			P(DD, ic, id)=P(DD, id, ic)=tmp;
		}
	}
	dmat* CCi=dpinv(CC, 0);
	dmat* tmp=0;
	dmm(&tmp, 0, CCi, DD, "nn", 1);
	dmm(&DD, 0, tmp, CCi, "nt", 1);
	dfree(tmp);
	lfree(embed);
	dfree(CCi);
	dfree(turbspec);
	ccellfree(spect);
	return DD;
}

/**
	Diagnolize the covariance matrix and transform the mode.
	Random vector is expressed as y=mod*x. Its covariance is 
	<yy'>=mod*<x*x'>*mod=mod*cov*mod'; 
	The covariance of the coefficient x can be decomposed with svd: cov=U*S*V'. So we have
	<yy'>=mod*<x*x'>*mod=mod*U*S*v'*mod'; Refine the bases with KL=mod*U, we have
	<yy'>=KL*S*KL' where the covariance of coefficient is diagonal S.
	returns the diagnolized mode and the covariance of each modal coefficient.
*/
dcell* cov_diagnolize(const dmat* mod, /**<Input mode*/
	const dmat* cov  /**<Covariance of modes*/
){
	if(!mod || !cov) return NULL;
	dmat* U=0, * Vt=0, * S=0;
	dsvd(&U, &S, &Vt, cov);
	dmat* kl=0;
	dmm(&kl, 0, mod, U, "nn", 1);
	
	//Drop modes after a sudden drop in strength of over 100x.
	long i;
	for(i=0; i<S->nx-1; i++){
		if(P(S,i+1)<P(S,i)*1e-2){
			break;
		}
	}
	i++;
	if(i<S->nx){
		dbg("Drop %ld last columns.\n", S->nx-i);
		dresize(kl, NX(kl), i);
	}
	dfree(U);
	dfree(Vt);
	dcell *res=dcellnew(2,1);
	P(res,0)=kl;
	P(res,1)=S;
	return res;
}

/**
   see KL_vonkarman()
 */
static dcell* KL_vonkarman_do(const loc_t* loc, int ttr, real iac, real L0){
	if(!loc) return NULL;
	dcell *kl=NULL;
	dmat *modz=NULL;
	dmat *cov=NULL;
	modz=dnew(loc->nloc, loc->nloc);
	real val=sqrt(loc->nloc); //this ensures rms wfe is 1, or orthonormal.
	daddI(modz, val);
	dadds(modz, -val/loc->nloc);//this ensure every column sum to 0 (no piston)
	if(ttr){
		loc_remove_ptt(modz, loc, NULL, NULL, 0);
		/*const dmat *tt=DMAT(loc);
		dmat *ptt=dpinv(tt, NULL, 1e-7);
		dmat *tmp=NULL;
		dmm(&tmp, 0, ptt, modz, "nn", 1);
		dmm(&modz, 1, tt, tmp, "nn", -1);
		dfree(ptt);
		dfree(tmp);*/
	}
	if(iac<=0){//not cubic
		cov=cov_vonkarman(loc, modz, L0);
	}else{//handle cubic influence function with over sampling
		dbg("KL_vonkarman_do: bicubic influence function\n");
		real xmin, xmax, ymin, ymax;
		dvecmaxmin(loc->locx, loc->nloc, &xmax, &xmin);
		dvecmaxmin(loc->locy, loc->nloc, &ymax, &ymin);
		/*real dx=loc->dx*0.5;
		real dy=loc->dy*0.5;
		loc_t *ploc=mksqloc((xmax-xmin+2*loc->dx)/dx+1, (ymax-ymin+loc->dy*2)/dy+1, dx, dy, xmin-loc->dx, ymin-loc->dy);*/
		dsp *h=mkh_cubic(loc, loc, 0, 0, 1, 0, iac);
		dmat *modu=NULL;
		dcellmm(&modu, h, modz, "nn", 1);//convert to dense modal matrix.
		dspfree(h);
		cov=cov_vonkarman(loc, modu, L0);
		//locfree(ploc);
	}
	kl=cov_diagnolize(modz, cov);
	dfree(modz);
	dfree(cov);
	return kl;
}
/**
   Create Karhunen-Loeve modes for which each mode coefficient is statistically
   independent in von Karman spectrum.

   The original method is to compute covariance of zernike modes in von Karman
   spectrum and diagnolize the covariance matrix. But this methods suffers for
   higher order systems because high order zernike modes are hard to generate on
   account of limited resolution of real precision floating point numbers. To
   make matters worse, to generate m KL modes, around 2*m Zernike modes are
   needed because each KL mode is a linear combination of all Zernike modes with
   same azimuthal (m) but higher radial (r) order. A work around of the
   precisiosn issue is to limit the absolute value of the zernike modes to
   within 100 to avoid excessive norm due to round off errors.

   The new method is to use zonal modes (identity matrix). The advantage is that
   the zonal modes spans the whole range of possible modes and are
   straightforward to generate.

   In theory, the KL modes computes should be independent on the starting mode,
   as long as the modes are independent and span the whole vector space for the
   coordinate.
*/
dcell *KL_vonkarman_full(const loc_t *loc, int ttr, real iac, real L0){
	if(!loc) return 0;
	dcell *res=0;
	if(loc->nloc<500){
		res=KL_vonkarman_do(loc, ttr, iac, L0);
	}else{
		uint32_t key=lochash(loc, 0);
		char fn[PATH_MAX];
		snprintf(fn, sizeof(fn), "KL/KL_vonkarman_%g_%ld_%d_%g_%u", L0, loc->nloc, ttr, iac, key);
		CACHE_FILE(res, fn, dcellread, ({res=KL_vonkarman_do(loc, ttr, iac, L0);}), writebin);
	}
	return res;
}
/*A convenient wrapper for simplified useage.*/
dmat *KL_vonkarman(const loc_t *loc, int ttr, real L0){
	if(!loc) return NULL;
	//Note 2025-01-06: KL modal matrix with cubic influenction function needs further
	//optimization. In current implementation, it does not perform as well as linear
	//version and needs more mode truncation in PWFS 2xDM case.
	dcell *res=KL_vonkarman_full(loc, ttr, 0, L0);
	dmat *kl=P(res, 0);P(res, 0)=NULL;
	dcellfree(res);
	return kl;
}
/**
   Generate FFT mode with period.
 */
dmat* fft_mode(const loc_t* loc, real D, real px, real py){
	if(!loc) return NULL;
	dmat* opd=dnew(loc->nloc, 1);
	real tx=2*M_PI*px/D;
	real ty=2*M_PI*py/D;
	for(long iloc=0; iloc<loc->nloc; iloc++){
		P(opd, iloc)=cos(loc->locx[iloc]*tx)*cos(loc->locy[iloc]*ty);
	}
	return opd;
}
