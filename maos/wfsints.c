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

#include <math.h>
#include <time.h>
#include "maos.h"
#include "sim.h"
#ifndef ROT_OTF
#define ROT_OTF 0
#endif
#define TIMING 1
/**
   \file wfsints.c Contains wfsints() that computes physical optics WFS
   subaperture images from OPDs.  */

/**
   compute physical optics images and add to ints.  be careful that fftw inverse
   fft doesn't have 1/n^2 scaling.  this function lives inside the threading
   routine wfsgradx  */
void wfsints(thread_t *thread_data){
    /*void wfsints(dcell *ints, ccell *psfout, dcell *pistatout,
	     const dmat *gradref,const PARMS_T *parms,
	     const POWFS_T *powfs,int iwfs,
	     const dmat *opd, const dmat *lltopd){
    */
    /**
       first, unwrap the data
     */
    WFSINTS_T *data=thread_data->data;
    dcell *ints=data->ints;
    ccell *psfout=data->psfout;
    dcell *pistatout=data->pistatout;
    const dmat *gradref=data->gradref;
    const PARMS_T *parms=data->parms;
    const POWFS_T *powfs=data->powfs;
    const int iwfs=data->iwfs;
    const dmat *opd=data->opd;
    const dmat *lltopd=data->lltopd;
    const int isa_start=thread_data->start;
    const int isa_end=thread_data->end;
    const int ipowfs=parms->wfs[iwfs].powfs;
    const int wfsind=parms->powfs[ipowfs].wfsind[iwfs];
    const int hasllt=(parms->powfs[ipowfs].llt!=NULL);
    const int nsa=powfs[ipowfs].pts->nsa;
    const int nx=powfs[ipowfs].pts->nx;
    const int ncompx=powfs[ipowfs].ncompx;//necessary size to build detector image.
    const int ncompy=powfs[ipowfs].ncompy;
    const int npsf=nx*parms->powfs[ipowfs].embfac;

    const int nwvl=parms->powfs[ipowfs].nwvl;

    cmat *psflarge=NULL;
    cmat *psf=NULL;
    cmat *psftmp=NULL;
    cmat *otf=cnew(ncompx,ncompy);
    cfft2plan(otf,-1);
    cfft2plan(otf,1);
    //make sure we use 2d for LGS, need to multiply to lotfc.
    const int use1d=(npsf>ncompy?1:0)&&(!hasllt)&&(!lltopd);
    const int ncompm=ncompx>ncompy?ncompx:ncompy;
    const int nxsa=nx*nx;
    if(use1d){
	psflarge=cnew(npsf,npsf);
	cfft2partialplan(psflarge,ncompm, -1);
	//there should be no inverse fft2 here.
	psf=cnew(ncompx,ncompy);
    }else{
	psf=cnew(npsf,npsf);
    }
    cfft2plan(psf,-1);
    cfft2plan(psf,1);
  
    int illt=0;
    if(hasllt){
	//has llt
	const int indwfs=parms->powfs[ipowfs].wfsind[iwfs];
	illt=parms->powfs[ipowfs].llt->i[indwfs];
    }
    cmat *fftpsfout=NULL;
    cmat *(*ppsfout)[nsa]=NULL;
    if(psfout){
	assert(psfout->nx==nsa && psfout->ny==nwvl);
	ppsfout=(void*)psfout->p;
	fftpsfout=cnew(psf->nx, psf->ny);
	cfft2plan(fftpsfout,1);
    }
    dmat *(*ppistatout)[nsa]=NULL;
    double *gx=NULL; double *gy=NULL;
    if(pistatout){
	assert(pistatout->nx==nsa && pistatout->ny==nwvl);
	ppistatout=(void*)pistatout->p;
	psftmp=cnew(psf->nx,psf->ny);
	cfft2plan(psftmp,1);
	if(gradref){
	    gx=gradref->p;
	    gy=gradref->p+nsa;
	}
    }
    double *realamp=powfs[ipowfs].realamp[wfsind];
    
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	
	double wvl=parms->powfs[ipowfs].wvl[iwvl];
	const double dtheta1=(powfs[ipowfs].pts->nx*parms->powfs[ipowfs].embfac
			      *powfs[ipowfs].pts->dx)/wvl;
	cmat *lotfc=NULL;
	//const double lltxc=parms->powfs[ipowfs].llt[iwfs].ox;
	if(lltopd){
	    const int nlx=powfs[ipowfs].llt->pts->nx;
	    const int nlpsf=nlx*parms->powfs[ipowfs].embfac;

	    lotfc=cnew(nlpsf,nlpsf);
	    cfft2plan(lotfc,-1);
	    cfft2plan(lotfc,1);
	    //build otf. use same config as subaperture
	    cembed_wvf(lotfc,lltopd->p,powfs[ipowfs].llt->amp->p,nlx,nlx,wvl,0);
	    cfft2(lotfc,-1);
	    if(npsf != nlpsf){//moved to before the second fft on 2011-08-07
		cmat *tmp=cnew(npsf, npsf);
		cfft2plan(tmp, 1);
		ccpcorner(tmp, lotfc, C_FULL);
		cfree(lotfc);
		lotfc=tmp;
	    }
	    cabs2toreal(lotfc);
	    cfft2(lotfc,1);//lotfc has peak nlpsf*nlpsf in corner.
	    //lotfc need to be scaled by 1/(npsf*npsf).
	    cscale(lotfc,1./(double)((long)nlpsf*nlpsf));//max of 1

	}
#if ROT_OTF == 1
	double xscale=(double)npsf/(double)ncompx;
	double yscale=(double)npsf/(double)ncompy;
#endif
	/*

	  Originally norm=(powfs[ipowfs].pts->area[isa]/maxarea)
	  /(powfs[ipowfs].pts->sumamp2[isa]*npsf*npsf*ncompx*ncompy);
	  with maxarea=max(powfs[ipowfs].pts->area[isa]);
	  where 
	  (ncompx*ncompy) is to cancel out the scaling in DTF nominal.
	  *npsf*npsf is to cancel out psf normalization

	  The area is normalized by square subaperture area.
	  We have powfs[ipowfs].pts->area[isa]=powfs[ipowfs].pts->sumamp2[isa]*dx*dx/(dsa*dsa);
	  therefore powfs[ipowfs].pts->area[isa]/powfs[ipowfs].pts->sumamp2[isa])
	  =dsa*dsa/(dx*dx)=pts->nx*pts->nx;
	  therefore norm=(1/maxarea)/(pts->nx*pts->nx*npsfx*npsfy*ncompx*ncompy);
	  we set 1/maxarea to powfs[ipowfs].areascale .
	  so the new norm becomes
	  norm=norm_otf/(ncompx*nocmpy);
	  with norm_otf=areascale/(pts->nx*pts->nx*npsf*npsf));

	  
	  Notice that fftw ifft doesn't normalize. a forward
	  and backward fft doesn't bring data back to
	  original value, instead Nx*Ny times larger.

	  Notice: the non-rotated case is handled properly
	  in ROT_OTF or ROT_PSF mode without executing extra
	  FFT.
	 */

	//normalize complex psf
	double norm_psf=sqrt(powfs[ipowfs].areascale)/((double)powfs[ipowfs].pts->nx*npsf);
	double norm_otf;
#if ROT_OTF == 1
	//norm_otf noramlized otf before mul ccwm to 1.
	norm_otf=norm_psf*norm_psf;
#else
	if(lltopd || ppistatout){
	    norm_otf=norm_psf*norm_psf/((double)psf->nx*psf->ny);
	}else{
	    norm_otf=norm_psf*norm_psf;
	}
#endif
	//this ncompx*ncompy is due to cfft2 after 
	//cwm and detector transfer function.
	double norm=norm_otf/((double)ncompx*ncompy);
	//normalized pistat
	double norm_pistat=norm_psf*norm_psf/((double)psf->nx*psf->ny);
	int multi_nominal=(powfs[ipowfs].dtf[iwvl].si->nx==nsa);
	cmat *nominal=NULL;
	dsp *si=NULL;
	if(!multi_nominal){
	    /*true only if 1) no elongation, 2) no radpix or 3) radpix && radrot is true*/
	    if(!powfs[ipowfs].dtf[iwvl].fused){
		nominal=powfs[ipowfs].dtf[iwvl].nominal->p[0];
	    }
	    si=powfs[ipowfs].dtf[iwvl].si->p[0];
	}
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
	int rotpsfotf=(hasllt && parms->powfs[ipowfs].radrot);
	double angle=0;
	for(int isa=isa_start; isa<isa_end; isa++){
	    if(multi_nominal){
		if(!powfs[ipowfs].dtf[iwvl].fused){
		    nominal=powfs[ipowfs].dtf[iwvl].nominal->p[isa+nsa*illt];
		}
		si=powfs[ipowfs].dtf[iwvl].si->p[isa+nsa*illt];
	    }
	    if(rotpsfotf){
		angle=powfs[ipowfs].srot->p[illt]->p[isa];
	    }
	    int ioffset=isa*nxsa;
	    //embed amp/opd to complex wvf with a embedding factor of 2.
	    if(use1d){
		cembed_wvf(psflarge,opd->p+ioffset, 
			   realamp+ioffset,nx,nx,
			   parms->powfs[ipowfs].wvl[iwvl],0);
		//fft
		cfft2partial(psflarge,ncompm, -1);
		//move from psflarge to psf
		ccpcorner(psf,psflarge,C_FULL);
	    }else{
		cembed_wvf(psf,opd->p+ioffset, 
			   realamp+ioffset,nx,nx,
			   parms->powfs[ipowfs].wvl[iwvl],0);
		cfft2(psf,-1);
		//form PSF.
	    }
	    //sum(psf) is sumamp2*npsf*npsf
	    if(ppsfout){
		ccp(&fftpsfout, psf);
		cscale(fftpsfout, norm_psf);
		cifft2(fftpsfout,1);//peak in corner. WVF. cifft2 is normalized.
		cembed(ppsfout[iwvl][isa], fftpsfout, 0, C_FULL);
		/*
		  ccp(&ppsfout[iwvl][isa], psf);
		  cscale(ppsfout[iwvl][isa], norm_psf);
		  cfftshift(ppsfout[iwvl][isa]);//peak in center.
		*/
	    }
	    cabs2toreal(psf);//peak in corner
#if ROT_OTF == 1//deprecated
	    cfft2(psf,-1);
#elif ROT_OTF == 0
	    if(lltopd || ppistatout){
		cfft2(psf,-1);//turn to otf. peak in corner
	    }
#endif
	    if(lltopd){//lotfc
		ccwm(psf,lotfc);//normalization done in gen of lotfc.
	    }
	    if(ppistatout){
		//remove tip/tilt from OTF using tilt reference.
		ccp(&psftmp,psf);//peak is in the corner
		if(gradref){
		    ctilt(psftmp,-gx[isa]*dtheta1,-gy[isa]*dtheta1,0);
		}
		cfft2(psftmp,1);//back to psf.
		cfftshift(psftmp);
		cabstoreal(psftmp);//take abs. peak at center.
		if(!gradref){
		    cshift2center(psftmp,0.5,0.5);
		}
		creal2d(&ppistatout[iwvl][isa],1,psftmp,norm_pistat);
	    }
	    if(ints){
#if ROT_OTF == 1 //rotate OTF preserves sum(PSF)
		cfftshift(psf);//OTF with peak in center
		//rotate from xy to ra, so negative angle.
		//Rotate and scale otf (in psf var) into otf.
		cembedscaleout(otf,psf,xscale,yscale,-angle,C_FULL);
		cfftshift(otf);//put otf peak in corner
#elif ROT_OTF == 0 //rotate PSF. preserves maximum PSF (strehl)
		if(lltopd || ppistatout){
		    //turn back to PSF. peak in corner. should be real
		    cfft2(psf,1);
		}
		cfftshift(psf);//peak in center
		cembed(otf,psf,-angle,C_REAL);//notice otf and psf may have different size
		cfftshift(otf);//peak in corner
		cfft2(otf,-1);//turn to OTF. peak in corner
#endif
		//fft to otf. peak in corner
		if(hasllt){//has llt, ETF
		    (*pccwm)(otf,nominal,petf[illt][isa]);
		}else{
		    ccwm(otf,nominal);
		}
		//max(otf) is 1. peak in corner 
		//Now peak in center because nominal has excessive freq. 
		cfft2(otf,1);
		spmulcreal(ints->p[isa]->p,si,
			   otf->p, parms->wfs[iwfs].wvlwts[iwvl]*norm);
	    }
	}//isa
	if(lotfc) cfree(lotfc);
    }//iwvl
    cfree(psf);
    cfree(psflarge);
    cfree(psftmp);
    cfree(otf);
    if(psfout){
	cfree(fftpsfout);
    }
}
