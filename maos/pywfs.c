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

/**
   \file pywfs.c
   Routines related to pyramid WFS
*/
#include "common.h"
#include "setup_powfs.h"
void pywfs_fft(dcell **pupim, PYWFS_T *pywfs, dmat *opd);
/**
   Complex wavefront created by pyramid 
*/

void pywfs_setup(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs){
    pywfs_free(powfs[ipowfs].pywfs);
    PYWFS_T *pywfs=powfs[ipowfs].pywfs=calloc(1, sizeof(PYWFS_T));
    map_t *map=0;
    double dx=parms->powfs[ipowfs].dx; 
    create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, 0, dx, dx, 0, 0, 0, 0, 0, 0);
    pywfs->loc=map2loc(map);
    mapfree(map);
    pywfs->amp=mkamp(pywfs->loc, aper->ampground, 
		     parms->misreg.pupil->p[0],parms->misreg.pupil->p[1], 
		     parms->aper.d, parms->aper.din);
    pywfs->locfft=locfft_init(pywfs->loc, pywfs->amp, 0, parms->powfs[ipowfs].wvl, 0);
    pywfs->wvlwts=ddup(parms->powfs[ipowfs].wvlwts);
    pywfs->modulate=parms->powfs[ipowfs].modulate;
    pywfs->order=parms->powfs[ipowfs].order;
    long nembed=pywfs->locfft->nembed->p[0];
    int nwvl=pywfs->locfft->wvl->nx;
    double wvlmin, wvlmax;
    dmaxmin(parms->powfs[ipowfs].wvl->p, nwvl, &wvlmax, &wvlmin);
    double dtheta_min=wvlmin/(dx*nembed);
    //size of the part of the PSF captured by pyramid
    long ncomp=nembed;
    if(parms->powfs[ipowfs].fov){
	ncomp=ceil(parms->powfs[ipowfs].fov/dtheta_min*0.5)*2;
	if(ncomp>nembed){
	    error("nembed=%ld is smaller than ncomp=%ld\n", nembed, ncomp);
	}
    }
    long ncomp2=ncomp/2;
    pywfs->pyramid=ccellnew(nwvl, 1);
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	pywfs->pyramid->p[iwvl]=cnew(ncomp, ncomp);
	PCMAT(pywfs->pyramid->p[iwvl], pp);
	dcomplex coeff=M_PI*0.5*I;
	long skip=0;
	if(parms->powfs[ipowfs].fov){//Limit fov per wvl
	    double dtheta=parms->powfs[ipowfs].wvl->p[iwvl]/(dx*nembed);
	    int nstop=ceil(parms->powfs[ipowfs].fov/dtheta*0.5)*2;
	    skip=(ncomp-nstop)/2;
	}
	for(long iy=skip; iy<ncomp-skip; iy++){
	    for(long ix=skip; ix<ncomp-skip; ix++){
		pp[iy][ix]=cexp((abs(iy-ncomp2)+abs(ix-ncomp2))*coeff);
	    }
	}
    }

    cmat *nominal=pywfs->nominal=cnew(ncomp, ncomp);
    cfft2plan(nominal, -1);
    cfft2plan(nominal, 1);
    PCMAT(nominal, pn);
    long npix=pywfs->order;
    double pixmeter=parms->aper.d/npix;//size of detector pixel in meter
    double dx2=dx*nembed/ncomp;//sampling of pupil after inverse fft
    double du=1./(dx2*ncomp);
    double dup=pixmeter*du;
    double pdmeter=pow(pixmeter/dx2, 2);
    double pixblur=parms->powfs[ipowfs].pixblur;
    double e0b=-2*pow(M_PI*pixblur*pixmeter*du, 2);
    for(int iy=0; iy<ncomp; iy++){
	int jy=iy-ncomp2;
	for(int ix=0; ix<ncomp; ix++){
	    int jx=ix-ncomp2; 
	    pn[iy][ix]=sinc(jy*dup)*sinc(jx*dup)*pdmeter;
	    if(pixblur){
		pn[iy][ix]*=exp(e0b*(jx*jx+jy*jy));
	    }
	}
    }
    cfftshift(nominal);
    pywfs->si=cellnew(4,1);//for each quadrant.
    //Make loc_t symmetric to ensure proper sampling onto detector
    powfs[ipowfs].saloc=mksqloc(npix, npix, pixmeter, pixmeter, 
				(-npix*0.5+0.5)*pixmeter, (-npix*0.5+0.5)*pixmeter);
    loc_t *loc_fft=mksqloc(ncomp, ncomp, dx2, dx2, (-ncomp2+0.5)*dx2, (-ncomp2+0.5)*dx2);
    for(int iy=0; iy<2; iy++){
	for(int ix=0; ix<2; ix++){
	    pywfs->si->p[ix+iy*2]=mkh(loc_fft, powfs[ipowfs].saloc, NULL, 
				      ((ix-0.5)*dx2*ncomp2), ((iy-0.5)*dx2*ncomp2),
				      1, 0, 0);	    
	}
    }
    locfree(loc_fft);

    {
	dcell *pupim=0;
	int nn=3;
	dmat *opds=zernike(pywfs->locfft->loc, parms->aper.d, 3);
	cellarr *pupsave=cellarr_init(nn,opds->ny,"pupim");
	for(int im=0; im<opds->ny; im++){
	    for(int j=0; j<nn; j++){
		info2("j=%d\n", j);
		dmat *opd=dsub(opds, 0, 0, im, 1);
		dscale(opd, (j+1)*1.e-7);
		dcellzero(pupim);
		pywfs_fft(&pupim, powfs[ipowfs].pywfs, opd);
		cellarr_push(pupsave, j+im*nn, pupim);
		dfree(opd);
	    }
	}
	cellarr_close(pupsave);
	cellfree(pupim);
	dfree(opds);
	exit(0);
    }
    if(parms->save.setup){
	cellwrite(pywfs->loc, "pywfs_loc");
	cellwrite(pywfs->amp, "pywfs_amp");
	cellwrite(pywfs->locfft->embed, "pywfs_embed");
	cellwrite(pywfs->pyramid, "pyramid");
	cellwrite(nominal, "nominal");
	cellwrite(pywfs->si, "si");
    }
}
/**
   Perform FFT on each quadrant of the PSF separately.
*/
void pywfs_fft(dcell **pupim, PYWFS_T *pywfs, dmat *opd){
    locfft_t *locfft=pywfs->locfft;
    ccell *psfs=0;
    locfft_psf(&psfs, locfft, opd, NULL, 1);//PSF sum to 1.
    int nwvl=locfft->wvl->nx;
    double dx=locfft->loc->dx;
    long nembed=locfft->nembed->p[0];
    long nembed2=nembed/2;
    dmat *wvlwts=pywfs->wvlwts;
    //position of pyramid for modulation
    int pos_n=1;
    double pos_r=pywfs->modulate*0.5;
    if(pos_r){
	pos_n=8;
    }
    long ncomp=pywfs->nominal->nx;
    long ncomp2=ncomp/2;
    cmat *otf=cnew(ncomp, ncomp);
    cfft2plan(otf, -1);
    cfft2plan(otf, 1);
    dmat *pupraw=dnew(ncomp, ncomp);
    for(int ipos=0; ipos<pos_n; ipos++){
	double theta=2*M_PI*ipos/pos_n;
	double posx=cos(theta)*pos_r;
	double posy=sin(theta)*pos_r;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    double dtheta=locfft->wvl->p[iwvl]/(dx*nembed);
	    long offy=(long)round(posy/dtheta);
	    long offy2=nembed2+offy-ncomp2;
	    long iy0=MAX(-offy2, 0);
	    long ny2=MIN(ncomp+offy2, nembed);

	    long offx=(long)round(posx/dtheta);
	    long offx2=nembed/2+offx-ncomp2;
	    long ix0=MAX(-offx2, 0);
	    long nx2=MIN(ncomp+offx2, nembed);

	    czero(otf);
	    dcomplex *pyramid=pywfs->pyramid->p[iwvl]->p;
	    for(long iy=iy0; iy<ny2-offy2; iy++){
		for(long ix=ix0; ix<nx2-offx2; ix++){
		    long indin=ix+offx2+(iy+offy2)*nembed;
		    long indout=ix+iy*ncomp;
		    otf->p[indout]=psfs->p[iwvl]->p[indin]*pyramid[indout];
		}
	    }
	    cfft2(otf, 1);
	    cabs22d(&pupraw, 1., otf, wvlwts->p[iwvl]/(ncomp*ncomp*pos_n));
	}//for iwvl
    }//for ipos
    if(!(*pupim)){
	(*pupim)=dcellnew(2, 2);
    }

    ccpd(&otf, pupraw);
    cfft2(otf, -1);
    ccwm(otf, pywfs->nominal);
    cfft2(otf, 1);
    for(int i=0; i<4; i++){
	if(!(*pupim)->p[i]){
	    (*pupim)->p[i]=dnew(pywfs->order, pywfs->order);
	}
	dspmulcreal((*pupim)->p[i]->p, pywfs->si->p[i], otf->p, 1./(ncomp*ncomp));
    }
    cellfree(psfs);
    dfree(pupraw);
    cfree(otf);
}
void pywfs_free(PYWFS_T *pywfs){
    if(!pywfs) return;
    locfree(pywfs->loc);
    dfree(pywfs->amp);
    locfft_free(pywfs->locfft);
    cellfree(pywfs->pyramid);
    cfree(pywfs->nominal);
    dspcellfree(pywfs->si);

}
