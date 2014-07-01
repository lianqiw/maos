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

#include <math.h>
#include <time.h>
#include "maos.h"
#include "sim.h"
#define TIMING 1
/**
   \file wfsints.c Contains wfsints() that computes physical optics WFS
   subaperture images from OPDs.  */

/**
   compute physical optics images and add to ints.  be careful that fftw inverse
   fft doesn't have 1/n^2 scaling.  this function lives inside the threading
   routine wfsgradx  
   
   Notice that the uplink psf have to have the same sampling as the subaperture psf.
   uplink psf has sampling lambda/(nlwvf*ldx)
   subaperture has sampling lambda/(nwvf*dx)
   When the uplink pupil is larger than subaperture, we need to scale ldx accordingly.
*/
void wfsints(thread_t *thread_data){
    /* first, unwrap the data */
    WFSINTS_T *data=thread_data->data;
    const PARMS_T *parms=data->parms;
    const POWFS_T *powfs=data->powfs;
    const int iwfs=data->iwfs;
    const dmat *opd=data->opd;
    const dmat *lltopd=data->lltopd;
    const int isa_start=thread_data->start;
    const int isa_end=thread_data->end;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int wfsind=parms->powfs[ipowfs].wfsind->p[iwfs];
    const int hasllt=(parms->powfs[ipowfs].llt!=NULL);
    const int illt=hasllt?parms->powfs[ipowfs].llt->i->p[wfsind]:0;
    const double *srot=(hasllt && parms->powfs[ipowfs].radrot)?powfs[ipowfs].srot->p[illt]->p:NULL;
    const int nsa=powfs[ipowfs].pts->nsa;
    const int ncompx=powfs[ipowfs].ncompx;/*necessary size to build detector image. */
    const int ncompy=powfs[ipowfs].ncompy;
    const int notf=MAX(ncompx,ncompy);
    const int nopd=powfs[ipowfs].pts->nx;
    const int nxsa=nopd*nopd;
    const int nwvf=nopd*parms->powfs[ipowfs].embfac;
    const int use1d=(nwvf>2*ncompy?1:0)&&(!hasllt)&&(!lltopd);
    const int nwvl=parms->powfs[ipowfs].nwvl;
    dcell *ints=data->ints;
    dcell *pistatout=data->pistatout; 
    const int isotf=(lltopd || pistatout);
    cmat *wvf=NULL;
    cmat *psf=NULL;
    cmat *psftmp=NULL;
    cmat *otf=NULL;
    cmat *lotfc=NULL;
    cmat *lwvf=NULL;
    /*this coefficient normalize the complex psf so its abs2 sum to 1.*/
    double norm_psf=sqrt(powfs[ipowfs].areascale)/(double)(powfs[ipowfs].pts->nx*nwvf);
    /*normalized pistat. notf is due to a pair of FFT on psf. */
    double norm_pistat=norm_psf*norm_psf/((double)notf*notf);
    /*this ncompx*ncompy is due to cfft2 after cwm and detector transfer function. */
    double norm_ints=parms->wfs[iwfs].siglevsim*norm_psf*norm_psf/((double)ncompx*ncompy);
    /*wvf first contains the complex wavefront, embed to get nyquist sampling.*/
    wvf=cnew(nwvf,nwvf);
    if(use1d){/* use 1d fft for NGS if nwvl >> notf*/
	cfft2partialplan(wvf,notf, -1);
    }else{
	cfft2plan(wvf, -1);
    }
    /* psf contains the psf/otf necessary to cover the focal plane. square */
    if(nwvf!=notf){
	psf=cnew(notf,notf);
    }else{
	psf=wvf;
    }
    cfft2plan(psf,-1);
    cfft2plan(psf,1);
    /* otf contains the psf/otf used to generate detecter image. maybe rectangular*/
    if(notf!=ncompx || notf!=ncompy || srot){
	otf=cnew(ncompx,ncompy);
	cfft2plan(otf,-1);
	cfft2plan(otf,1);
	if(isotf){/*there is an additional pair of FFT*/
	    norm_ints/=(double)(notf*notf);
	}
    }else{
	otf=psf;
    }
    /* there is uplink beam */
    if(lltopd){
	const int nlx=powfs[ipowfs].llt->pts->nx;
	const int nlwvf=nlx*parms->powfs[ipowfs].embfac;
	lwvf=cnew(nlwvf,nlwvf);
	cfft2plan(lwvf,-1);
	cfft2plan(lwvf,1);
	if(nlwvf != notf){
	    lotfc=cnew(notf,notf);
	    cfft2plan(lotfc,-1);
	    cfft2plan(lotfc,1);
	}else{
	    lotfc=lwvf;
	}
    }
    cmat *fftpsfout=NULL;
    /* hold the complex psf to save to file. */
    cmat *(*ppsfout)[nsa]=NULL;
    /* need to output psf time history */
    if(data->psfout){
	ppsfout=(void*)data->psfout->p;
	fftpsfout=cnew(psf->nx, psf->ny);
	cfft2plan(fftpsfout,1);
    }
    dmat *(*ppistatout)[nsa]=NULL;
    double *gx=NULL; double *gy=NULL;
    /* need to output pixel itnensity averages */
    if(pistatout){
	assert(pistatout->nx==nsa && pistatout->ny==nwvl);
	ppistatout=(void*)pistatout->p;
	psftmp=cnew(psf->nx,psf->ny);
	cfft2plan(psftmp,1);
	/* the gradient reference for pistatout*/
	if(data->gradref){
	    gx=data->gradref->p;
	    gy=gx+nsa;
	}
    }
    double *realamp=powfs[ipowfs].realamp->p[wfsind]->p;
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	const double wvl=parms->powfs[ipowfs].wvl->p[iwvl];
	const double dtheta1=(nwvf*powfs[ipowfs].pts->dx)/wvl;
	/* uplink llt opd*/
	if(lltopd){
	    const int nlx=powfs[ipowfs].llt->pts->nx;
	    const int nlwvf=nlx*parms->powfs[ipowfs].embfac;
	    /*embed opd to compute complex pupil function*/
	    cembed_wvf(lwvf,lltopd->p,powfs[ipowfs].llt->amp->p,nlx,nlx,wvl,0);
	    /*turn to complex psf*/
	    cfft2(lwvf,-1);
	    /*turn to psf*/
	    cabs2toreal(lwvf);
	    /*turn to otf. */
	    cfft2(lwvf, 1);
	    if(lwvf != lotfc){
		/* uplink llt has different aperture size than LGS subaperture*/
		ccpcorner(lotfc, lwvf, C_FULL);
	    }
	    /*lotfc has peak nlwvf*nlwvf in corner. */
	    cscale(lotfc,1./(double)((long)nlwvf*nlwvf));/*max of 1 */
	}
	int multi_nominal=(powfs[ipowfs].dtf[iwvl].si->nx==nsa);
	/* nominal and si are used to sampled PSF onto detector pixels */
	cmat *nominal=NULL;
	dsp *si=NULL;
	if(!multi_nominal){
	    /*true only if 1) no elongation, 2) no radpix or 3) radpix && radrot is true*/
	    if(!powfs[ipowfs].dtf[iwvl].fused){
		nominal=powfs[ipowfs].dtf[iwvl].nominal->p[0];
	    }
	    si=powfs[ipowfs].dtf[iwvl].si->p[0];
	}
	/* elongation due to uplink projection */
	cmat *(*petf)[nsa]=NULL;
	void (*pccwm)(cmat*,const cmat*,const cmat*)=NULL;
	if(hasllt){
	    if(powfs[ipowfs].etfsim[iwvl].p1){
		petf=(void*)powfs[ipowfs].etfsim[iwvl].p1->p;
		pccwm=ccwm3col;
	    }else{
		petf=(void*)powfs[ipowfs].etfsim[iwvl].p2->p;
		pccwm=ccwm3;
	    }
	}
	/* now big loop over subapertures */
	for(int isa=isa_start; isa<isa_end; isa++){
	    if(multi_nominal){
		if(!powfs[ipowfs].dtf[iwvl].fused){
		    nominal=powfs[ipowfs].dtf[iwvl].nominal->p[isa+nsa*illt];
		}
		si=powfs[ipowfs].dtf[iwvl].si->p[isa+nsa*illt];
	    }
	    int ioffset=isa*nxsa;
	    /*embed amp/opd to complex wvf with a embedding factor of 2. */
	    cembed_wvf(wvf,opd->p+ioffset, 
		       realamp+ioffset,nopd,nopd,
		       parms->powfs[ipowfs].wvl->p[iwvl],0);
	    if(use1d){ /*use 1d fft */
		cfft2partial(wvf,notf, -1);
	    }else{
		cfft2(wvf,-1); /*use 2d fft to form PSF. */
	    }
	    if(psf!=wvf){
		/*copy the peaks (at corner) from wvf to psf */
		ccpcorner(psf, wvf, C_FULL);
	    }

	    /*output complex pupil function to use in skyc*/
	    if(ppsfout){
		ccp(&fftpsfout, psf);
		/*notf * notf to cancel out the effect of fft pair (one here, one later in skyc)*/
		cscale(fftpsfout, norm_psf/((double)notf*notf));
		/*peak in corner. become WVF in center.*/
		cfft2(fftpsfout,1);
		/*output center of the complex pupil function.*/
		cembedc(ppsfout[iwvl][isa], fftpsfout, 0, C_FULL);
	    }
	    /* form PSF with peak in corner*/
	    cabs2toreal(psf);
	    /* need to turn to otf to add llt contribution or output pixel intensities.*/
	    if(isotf){
		cfft2(psf,-1);   /*turn to otf. peak in corner */
		if(ppistatout){  /*The pistat does not include uplink effect*/
		                 /*copy to temporary array. peak in in corner*/
		    ccp(&psftmp,psf);
		    if(gx){      /*remove tip/tilt from OTF using tilt reference. */
			ctilt(psftmp,-gx[isa]*dtheta1,-gy[isa]*dtheta1,0);
		    }
		    cfft2(psftmp,1);   /*back to psf. */
		    cfftshift(psftmp); /*peak in center. */
		    cabstoreal(psftmp);/*take abs.*/
		    if(!gx){           /*no gradient reference. compute using cog and shift.*/
			cshift2center(psftmp,0.5,0.5);
		    }
		                       /*accumulate the real part of psf*/
		    creal2d(&ppistatout[iwvl][isa],1,psftmp,norm_pistat);
		}
		if(lltopd){            /*add uplink otf */
		    ccwm(psf,lotfc);   /*normalization done in gen of lotfc. */
		}
		/* we have otf here in psf*/
	    }/* else: we have psf here in psf*/
	    if(ints){
		if(!isotf /* Need to turn PSF to OTF*/ || otf!=psf /* Need to embed*/ ){
		    if(isotf){ /*turn back to PSF for embedding. peak in corner.*/
			cfft2(psf,1);
		    }
		    /* now we have PSF with peak in corner */
		    if(srot){
			cfftshift(psf);/*peak in center */
			cembedc(otf,psf,-srot[isa],C_REAL);/*notice otf and psf may have different size */
			cfftshift(otf);/*peak in corner */
		    }else if(otf!=psf){/*copy the corner (peak)*/
			ccpcorner(otf, psf, C_FULL);
		    }
		    cfft2(otf,-1);/*turn to OTF. peak in corner */
		}
		if(hasllt){/*has llt, multiply with DTF and ETF.*/
		    (*pccwm)(otf,nominal,petf[illt][isa]);
		}else{/*no uplink, multiply with DTF only.*/
		    ccwm(otf,nominal);
		}
		/*max(otf) is 1 after multiply with norm. peak in corner  */
		cfft2(otf,1);
		/*Now peak in center because nominal is pre-treated.  */
		spmulcreal(ints->p[isa]->p,si, otf->p, parms->wfs[iwfs].wvlwts->p[iwvl]*norm_ints);
	    }
	}/*isa */
    }/*iwvl */
    if(otf!=psf) cfree(otf);
    if(psf!=wvf) cfree(psf);
    cfree(wvf);
    cfree(psftmp);
    cfree(fftpsfout);
    if(lltopd){
	if(lotfc!=lwvf) cfree(lotfc);
	cfree(lwvf);
    }
}
