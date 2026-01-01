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
#include "common.h"
#include "sim.h"
#define TIMING 0
#if TIMING
#define TIM0 static real tk1=0,tk2=0,tk3=0,tk4=0;static int tkct=0;real tk=0,tk0=myclockd();tkct++;
#define TIM(A) tk=myclockd(); tk##A+=tk-tk0;tk0=tk;
#define TIM1 info2("wfsints timing: copy %.3f fft %.3f cwm %.3f spmul %.3f tot %.3f\n", tk1/tkct, tk2/tkct, tk3/tkct, tk4/tkct, (tk1+tk2+tk3+tk4)/tkct)
#else
#define TIM0
#define TIM(A)
#define TIM1
#endif

/*
   Contains wfsints() that computes physical optics WFS
   subaperture images from OPDs.  */
/**
 * out=out*A*(B*wtb+C*wtc)
 * A, B and C are optional except that when C is present, B must be.
*/
static void apply_dtf_etf(cmat *out, dmat *A, cmat *B, real wtb, cmat *C, real wtc){
	if(C){
		if(A){
			for(long i=0; i<PN(out); i++){
				P(out, i)*=P(A, i)*(P(B, i)*wtb+P(C, i)*wtc);
			}
		}else{
			for(long i=0; i<PN(out); i++){
				P(out, i)*=(P(B, i)*wtb+P(C, i)*wtc);
			}
		}
	}else if(B){
		if(A){
			for(long i=0; i<PN(out); i++){
				P(out, i)*=P(A, i)*P(B, i);
			}
		}else{
			for(long i=0; i<PN(out); i++){
				P(out, i)*=P(B, i);
			}
		}
	}else if(A){
		for(long i=0; i<PN(out); i++){
			P(out, i)*=P(A, i);
		}
	}
}
/**
   compute physical optics images and add to ints.  be careful that fftw inverse
   fft doesn't have 1/n^2 scaling.  this function lives inside the threading
   routine wfsgradx

   Notice that the uplink have to have the same spatial sampling as the subaperture (dx).
   uplink psf has sampling lambda/(nlwvf*dx)
   subaperture has sampling lambda/(nwvf*dx)
   uplink otf has the sampling dx/lambda sampling as subaperture. 

   WHen notf is bigger than nwvf*embfac, the PSF is embeded which preserves the sampling.
   The uplink OTF need to be upsampled to dx/lambda *(nwvf/notf)
   
*/
void* wfsints(thread_t* thread_data){
	/* first, unwrap the data */
	wfsints_t* data=(wfsints_t*)thread_data->data;
	const parms_t* parms=global->parms;
	const powfs_t* powfs=global->powfs;
	const int iwfs=data->iwfs;
	const int ipowfs=parms->wfs[iwfs].powfs;
	const dmat* opd=data->opd;
	const dmat* lltopd=data->lltopd;
	const dmat* wvlwts=parms->wfs[iwfs].wvlwts;
	const int isa_start=thread_data->start;
	const int isa_end=thread_data->end;
	const int wfsind=P(parms->powfs[ipowfs].wfsind,iwfs);
	const int hasllt=(parms->powfs[ipowfs].llt!=NULL);
	const int illt=hasllt?P(parms->powfs[ipowfs].llt->i,wfsind):0;
	const int nsa=powfs[ipowfs].saloc->nloc;
	const int notfx=powfs[ipowfs].notfx;/*necessary size to build detector image. */
	const int notfy=powfs[ipowfs].notfy;
	const int nopd=powfs[ipowfs].pts->nxsa;
	const int saoffset=nopd*nopd;
	const int nwvf=nopd*parms->powfs[ipowfs].embfac;
	const int use1d=(nwvf>2*notfy?1:0)&&(!hasllt)&&(!lltopd);
	const int nwvl=parms->powfs[ipowfs].nwvl;
	dcell* ints=data->ints;
	dcell* pistatout=data->pistatout;
	cmat* wvf=NULL;
	cmat* psf=NULL;
	cmat* psftmp=NULL;
	cmat* lotfc=NULL;
	cmat* lwvf=NULL;
	//this coefficient normalize the complex psf so its abs2 sum to 1.
	//nopd scales the amplitude map, nwvf canceling the scaling of FFT.
	real norm_psf=sqrt(powfs[ipowfs].areascale)/(real)(nopd*nwvf);
	/*normalized pistat. notf is due to a pair of FFT on psf. */
	real norm_pistat=norm_psf*norm_psf/((real)notfx*notfy);
	/*this notfx*notfy is due to cfft2 after cwm and detector transfer function. */
	real norm_ints=parms->wfs[iwfs].sigsim*norm_psf*norm_psf/((real)notfx*notfy);
	/*wvf first contains the complex wavefront, embed to get nyquist sampling.*/
	wvf=cnew(nwvf, nwvf);
	/* psf contains the psf/otf necessary to cover the focal plane. square */
	if(nwvf!=notfx || notfx!=notfy){
		psf=cnew(notfx, notfy);
	} else{
		psf=wvf;
	}

	cmat* fftpsfout=NULL;
	/* hold the complex psf to save to file. */
	ccell* psfout=data->psfout;
	/* need to output psf time history */
	if(psfout){
		fftpsfout=cnew(NX(psf), NY(psf));
	}
	real* gx=NULL; real* gy=NULL;
	/* need to output pixel itnensity averages */
	if(pistatout){
		assert(NX(pistatout)==nsa&&NY(pistatout)==nwvl);
		psftmp=cnew(NX(psf), NY(psf));
		/* the gradient reference for pistatout*/
		if(data->gradref){
			gx=P(data->gradref);
			gy=gx+nsa;
		}
	}
	real* amp=P(PR(powfs[ipowfs].amp,wfsind));
	TIM0;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
		const real wvl=P(parms->powfs[ipowfs].wvl,iwvl);
		const real dtheta1=(nwvf*powfs[ipowfs].pts->dx)/wvl;//: 1/dtheta
		/* uplink llt opd*/
		if(lltopd){
			const int nlx=powfs[ipowfs].llt->pts->nxsa;
			const int nlwvf=nlx*parms->powfs[ipowfs].embfac;
			/*embed opd to compute complex pupil function*/
			if(!lwvf) {
				lwvf=cnew(nlwvf, nlwvf); 
			}else{
				czero(lwvf);
			}
			cembed_wvf(lwvf, P(lltopd), P(powfs[ipowfs].llt->amp), nlx, nlx, wvl, 0);
			/*turn to complex psf*/
			cfft2(lwvf, -1);
			/*turn to psf*/
			cabs2toreal(lwvf, 1./((real)nlwvf*nlwvf));//sum to 1. scaling is independent of notf
			/*if(isa_start==0){//print uplink spot size.
				cfftshift(lwvf);
				real fwhm=1.1774*cgauss_width(lwvf, 0.01)/dtheta1*RAD2AS;
				info("LGS %d on sky fwhm is %g\n", iwfs, fwhm);
				cfftshift(lwvf);
			}*/
			int islotf=0;//whether lotfc contains OTF
			cmat *lotfc2=NULL;
			if(nlwvf!=nwvf){//need to embed/crop uplink OTF
				/* uplink has different aperture size than LGS subaperture, but
				the same spatial and OTF sampling. Crop or pad OTF toget to get
				the same PSF/OTF sampling.*/
				islotf=1;cfft2(lwvf, -1);//turn to OTF
				lotfc2=cnew(nwvf, nwvf);//make uplink PSF(OTF) same size and sampling as downlink
				ccpcorner(lotfc2, lwvf, C_FULL);//crop PSF
			}else{
				lotfc2=lwvf; 
			}
			if(nwvf!=notfx||notfx!=notfy){//need to embed/crop uplink PSF
				if(islotf){
					cfft2(lotfc2, 1);//turn back to PSF
					cscale(lotfc2, 1./((real)nwvf*nwvf));
					islotf=0;
				}
				if(!lotfc){
					lotfc=cnew(notfx, notfy);
				}else{
					czero(lotfc);
				}
				ccpcorner(lotfc, lotfc2, C_FULL);//clip or expand PSF.
				if(lotfc2!=lwvf) cfree(lotfc2);
			}else{
				lotfc=lotfc2; lotfc2=NULL;
			}
			if(!islotf) cfft2(lotfc, -1);/*turn to otf. */
			/*lotfc has peak in lower left corner. */
		}
		int multi_nominal=(NX(powfs[ipowfs].dtf[iwvl].si)==nsa);
		/* nominal and si are used to sampled PSF onto detector pixels */
		dmat* nominal=NULL;
		dsp* si=NULL;
		if(!multi_nominal){
			/*true only if 1) no elongation, 2) no radpix */
			if(!hasllt){
				nominal=P(powfs[ipowfs].dtf[iwvl].nominal,0);
			}
			si=P(powfs[ipowfs].dtf[iwvl].si,0);
		}
		/* elongation due to uplink projection */
		ccell* petf1=NULL, * petf2=NULL;
		real etf1wt=1;
		real etf2wt=0;
		if(hasllt){
			petf1=powfs[ipowfs].etfsim[iwvl].etf;
			if(powfs[ipowfs].etfsim2){
				petf2=powfs[ipowfs].etfsim2[iwvl].etf;
				const int dtrat=parms->powfs[ipowfs].llt->coldtrat?parms->powfs[ipowfs].llt->coldtrat:parms->powfs[ipowfs].zoomdtrat;
				etf2wt=(real)(data->isim%dtrat)/(real)dtrat;
				etf1wt=1.-etf2wt;
			}
		}
		/* now big loop over subapertures */
		for(int isa=isa_start; isa<isa_end; isa++){
			if(multi_nominal){
				if(!hasllt){
					nominal=P(powfs[ipowfs].dtf[iwvl].nominal,isa,illt);
				}
				si=P(powfs[ipowfs].dtf[iwvl].si,isa,illt);
			}
			int ioffset=isa*saoffset;
			/*embed amp/opd to complex wvf with a embedding factor of 2. */
			cembed_wvf(wvf, P(opd)+ioffset,amp+ioffset, nopd, nopd, P(parms->powfs[ipowfs].wvl,iwvl), 0);
			TIM(1);//1 is copy
			if(use1d){ /*use 1d fft */
				cfft2partial(wvf, notfy, -1);
			} else{
				cfft2(wvf, -1); /*use 2d fft to form PSF. */
			}
			TIM(2);//2 is fft
			if(psf!=wvf){/*copy the peaks (at corner) from wvf to psf */
				ccpcorner(psf, wvf, C_FULL);
			}

			/*output complex pupil function to use in skyc*/
			if(psfout){
				ccp(&fftpsfout, psf);
				/*notf * notf to cancel out the effect of fft pair (one here, one later in skyc)*/
				cscale(fftpsfout, norm_psf/((real)nwvf*nwvf));
				/*peak in corner. become WVF in center.*/
				cfft2(fftpsfout, 1);
				/*output center of the complex pupil function.*/
				cembed(P(psfout, isa, iwvl), fftpsfout, 0);
			}
			/* form PSF with peak in corner*/
			cabs2toreal(psf, 1);
			TIM(1);
			/* need to turn to otf to add llt contribution or output pixel intensities.*/
			cfft2(psf, -1);   /*turn to otf. peak in corner */
			TIM(2);
			if(pistatout){  /*The pistat does not include uplink effect*/
				/*copy to temporary array. peak in in corner*/
				ccp(&psftmp, psf);
				if(gx){      /*remove tip/tilt from OTF using tilt reference. */
					ctilt(psftmp, -gx[isa]*dtheta1, -gy[isa]*dtheta1, 0);
				}
				cfft2(psftmp, 1);   /*back to psf. */
				cfftshift(psftmp); /*peak in center. */
				cabstoreal(psftmp);/*take abs.*/
				if(!gx){           /*no gradient reference. compute using cog and shift.*/
					cshift2center(psftmp, 0.5, 0.5);
				}
				/*accumulate the real part of psf*/
				creal2d(&P(pistatout, isa, iwvl), 1, psftmp, norm_pistat);
			}
			if(lltopd){            /*add uplink otf */
				ccwmc(psf, lotfc);   /*normalization done in gen of lotfc. */
			}
			/* we have otf here in psf*/
			if(ints){
				apply_dtf_etf(psf, nominal, petf1?P(petf1, isa, illt):0, etf1wt, petf2?P(petf2, isa, illt):0, etf2wt); TIM(3);
				/*max(otf) is 1 after multiply with norm. peak in corner  */
				cfft2(psf, 1); TIM(2);
				/*Now peak in center because nominal is pre-treated.  */
				dspmulcreal(P(P(ints,isa)), si, P(psf), P(wvlwts,iwvl)*norm_ints); TIM(4);
			}
		}/*isa */
	}/*iwvl */
	if(psf!=wvf) cfree(psf);
	cfree(wvf);
	cfree(psftmp);
	cfree(fftpsfout);
	if(lltopd){
		if(lwvf!=lotfc){
			cfree(lwvf);
		}
		cfree(lotfc);
	}
	TIM1;
	return NULL;
}
