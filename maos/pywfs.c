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
void setup_pywfs(POWFS_T *powfs, const PARMS_T *parms, APER_T *aper, int ipowfs){
    pywfs_free(powfs[ipowfs].pywfs);
    PYWFS_T *pywfs=powfs[ipowfs].pywfs=calloc(1, sizeof(PYWFS_T));
    map_t *map=0;
    double dx=parms->powfs[ipowfs].dx; 
    create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, 0, dx, dx, 0.5, dx*0.5, 0, 0, 0, 0);
    pywfs->loc=map2loc(map);
    mapfree(map);
    pywfs->amp=mkamp(pywfs->loc, aper->ampground, 
		     parms->misreg.pupil->p[0],parms->misreg.pupil->p[1], 
		     parms->aper.d, parms->aper.din);
    pywfs->locfft=locfft_init(pywfs->loc, pywfs->amp, 0, 
			      parms->powfs[ipowfs].wvl, parms->powfs[ipowfs].fov);
    long nembed=pywfs->locfft->nembed->p[0];
    int nwvl=pywfs->locfft->wvl->nx;
    double wvlmin, wvlmax;
    dmaxmin(parms->powfs[ipowfs].wvl->p, nwvl, &wvlmax, &wvlmin);
    double dtheta_min=wvlmin/(dx*nembed);
    //size of the part of the PSF captured by quadrant of pyramid
    long nquad;
    if(parms->powfs[ipowfs].fieldstop){
	nquad=ceil(parms->powfs[ipowfs].fieldstop/dtheta_min*0.25)*2;
	if(nquad*2<nembed){
	    error("nembed=%ld is too small. nquad=%ld\n", nembed, nquad);
	}
    }else{
	if(nembed&1){
	    error("nembed=%ld should be even instead\n", nembed);
	}
	nquad=nembed/2;
    }
    pywfs->psfquad=cnew(nquad, nquad);
    cfft2plan(pywfs->psfquad, 1);
    cfft2plan(pywfs->psfquad, -1);
    cmat *nominal=pywfs->nominal=cnew(nquad, nquad);
    cfft2plan(nominal, -1);
    cfft2plan(nominal, 1);
    PCMAT(nominal, pn);
    long nquad2=nquad/2;
    long npix=parms->powfs[ipowfs].order;
    double pixmeter=parms->aper.d/npix;//size of detector pixel in meter
    double dx2=dx*nembed/nquad;//sampling of pupil after inverse fft
    double du=1./(dx2*nquad);
    double dup=pixmeter*du;
    double pdmeter=pow(pixmeter/dx2, 2);
    double pixblur=parms->powfs[ipowfs].pixblur;
    double e0b=-2*pow(M_PI*pixblur*pixmeter*du, 2);
    for(int iy=0; iy<nquad; iy++){
	int jy=iy-nquad2;
	for(int ix=0; ix<nquad; ix++){
	    int jx=ix-nquad2; 
	    pn[iy][ix]=sinc(jy*dup)*sinc(jx*dup)*pdmeter;
	    if(pixblur){
		pn[iy][ix]*=exp(e0b*(jx*jx+jy*jy));
	    }
	}
    }
    cfftshift(nominal);
    loc_t *loc_fft=mksqloc(nquad, nquad, dx2, dx2, -nquad2*dx2, -nquad2*dx2);
    loc_t *loc_ccd=mksqloc(npix, npix, pixmeter, pixmeter, (-npix+1)*0.5*pixmeter, (-npix+1)*0.5*pixmeter);
    cwrite(nominal, "nominal");
    spwrite(pywfs->si, "si");
    locwrite(loc_fft, "loc_fft");
    locwrite(loc_ccd, "loc_ccd");
    pywfs->si=mkh(loc_fft, loc_ccd, NULL, 0, 0, 1, 0, 0);
    locfree(loc_fft);
    locfree(loc_ccd);

    {
	dmat *opd=dnew(pywfs->locfft->loc->nloc, 1);
	dcell *pupim=dcellnew(2,2);
	pywfs_fft(&pupim, powfs, parms, 0, opd);
	dcellwrite(pupim, "pupim_zero");
	dcellzero(pupim);
	for(int i=0; i<opd->nx; i++){
	    opd->p[i]=pywfs->locfft->loc->locx[i]*1e-6;
	}
	pywfs_fft(&pupim, powfs, parms, 0, opd);
	dcellwrite(pupim, "pupim_ttx");
	dcellzero(pupim);
	for(int i=0; i<opd->nx; i++){
	    opd->p[i]=pywfs->locfft->loc->locy[i]*1e-6;
	}
	pywfs_fft(&pupim, powfs, parms, 0, opd);
	dcellwrite(pupim, "pupim_tty");
	dcellzero(pupim);
	dfree(opd); exit(0);
    }
}
void pywfs_fft(dcell **pupim, POWFS_T *powfs, const PARMS_T *parms, int ipowfs, dmat *opd){
    PYWFS_T *pywfs=powfs[ipowfs].pywfs;
    locfft_t *locfft=pywfs->locfft;
    ccell *psfs=locfft_psf(locfft, opd, NULL, 1);//PSF sum to 1.
    dcell *pupraw=dcellnew(4,1);
    int nwvl=locfft->wvl->nx;
    double dx=locfft->loc->dx;
    int nembed=locfft->nembed->p[0];
    dmat *wvlwts=parms->powfs[ipowfs].wvlwts;
    //position of pyramid
    int pos_n=1;
    double pos_r=parms->powfs[ipowfs].modulate*0.5;
    if(pos_r){
	pos_n=8;
    }
    long nquad=pywfs->psfquad->nx;
    for(int ipos=0; ipos<pos_n; ipos++){
	double theta=2*M_PI*ipos/pos_n;
	double posx=cos(theta)*pos_r;
	double posy=sin(theta)*pos_r;
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    double dtheta=locfft->wvl->p[iwvl]/(dx*nembed);
	    long offx=(long)round(posx/dtheta);
	    long offy=(long)round(posy/dtheta);
	    //Loop over quadrant
	    for(int iqy=0; iqy<2; iqy++){
		long offy2=nembed/2*iqy+offy;
		long ny2=pywfs->psfquad->ny+offy2;
		if(ny2>nembed){
		    ny2=nembed;
		}
		for(int iqx=0; iqx<2; iqx++){
		    long offx2=nembed/2*iqx+offx;
		    long nx2=pywfs->psfquad->nx+offx2;
		    if(nx2>nembed){
			nx2=nembed;
		    }
		    czero(pywfs->psfquad);
		    for(long iy=0; iy<ny2-offy2; iy++){
			for(long ix=0; ix<nx2-offx2; ix++){
			    pywfs->psfquad->p[ix+iy*pywfs->psfquad->nx]
				=psfs->p[iwvl]->p[ix+offx2+(iy+offy2)*nembed];
			}
		    }
		    cfft2(pywfs->psfquad, 1);
		    cabs22d(&pupraw->p[iqx+iqy*2], 1., pywfs->psfquad, wvlwts->p[iwvl]/(nquad*nquad));
		}
	    }
	}//for iwvl
    }//for ipos
    for(int i=0; i<4; i++){
	if(!(*pupim)->p[i]){
	    (*pupim)->p[i]=dnew(pywfs->si->nx, 1);
	}
	ccpd(&pywfs->psfquad, pupraw->p[i]);
	cfft2(pywfs->psfquad, -1);
	ccwm(pywfs->psfquad, pywfs->nominal);
	cfft2(pywfs->psfquad, 1);
	spmulcreal((*pupim)->p[i]->p, pywfs->si, pywfs->psfquad->p, 1./(nquad*nquad*pos_n));
    }
    cellfree(psfs);
    cellfree(pupraw);

}
void pywfs_free(PYWFS_T *pywfs){
    if(!pywfs) return;
    locfree(pywfs->loc);
    dfree(pywfs->amp);
    locfft_free(pywfs->locfft);
    cfree(pywfs->psfquad);
    cfree(pywfs->nominal);
    spfree(pywfs->si);
}
