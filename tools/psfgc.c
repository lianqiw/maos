/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
 * Standalone code to sample IRIS PSFs onto detector for the Galactic center simulations.
 */
#include <unistd.h>
#include <getopt.h>

#include "../lib/aos.h"

static void usage(){
    info2("Usage: psfgc output.fits seed exposure x y npix pixsize pixoffx pixoffy blur\n"
	  "output.fits: stores the output\n"
	  "seed:        the atmospheric seed. valid: integers from 1 to 10\n"
	  "exposure:    in seconds. valid: 0, 4, 8, 12, 16, 20. 0 is for all length\n"
	  "x y:         the focal plane coordinate in arcsec. valid: from -14 to 14 at step of 2.\n"
	  "npix:        number of pixels in the detector. 0 will output PSF.\n"
	  "pixsize:     size of detector pixel in arcsec.\n"
	  "pixoffx:     offset of the center of image from center of detector along x axis in unit of detector pixel. 1 will offset by 1 pixel.\n"
	  "pixoffy:     offset along y axis.\n"
	  "Notice that the implementation error of 110 nm is accounted in\n"
	  );
}

typedef struct psfiris_t{
    int notf;
    int nwvl;
    double dx;
    double sumpsf;
    int npix;
    double pixsize;
    double pixoffx;
    double pixoffy;
    double blur;
    double imperr;
    double *wvls;
    dcell *psf_lgs;
    dcell *output;
    char *msg;
}psfiris_t;

static void psfiris_do(thread_t *info){
    psfiris_t *data=info->data;
    int ipsf=info->start;
    int nwvl=data->nwvl;
    int iwvl=ipsf%nwvl;
    int notf=data->notf;
    int npix=data->npix;
    double sumpsf=data->sumpsf;
    double dx=data->dx;
    double pixsize=data->pixsize;
    double pixoffx=data->pixoffx;
    double pixoffy=data->pixoffy;
    double blur=data->blur;
    double imperr=data->imperr;
    double *wvls=data->wvls;
    dcell *psf_lgs=data->psf_lgs;
    dcell *output=data->output;
    char *msg=data->msg;

    cmat *otf=cnew(notf, notf);
    cfft2plan(otf,1);
    cfft2plan(otf,-1);
    info2("%d ",ipsf);
    /*first create OTF of tt/ps modes on coarse sampling.*/
    double wvl=wvls[iwvl]*1e-6;
    double dtheta=wvl/(notf*dx);
    ccpd(&otf, psf_lgs->p[ipsf]);
    dfree(psf_lgs->p[ipsf]);
    cfftshift(otf);
    cfft2(otf, 1);
    double impst=exp(-pow(2*M_PI/wvl*imperr*1e-9,2))/(notf*notf);
    if(npix>0){
	ccell *dtf=NULL;
	spcell *si=NULL;
	mkdtf(&dtf, &si, notf, notf, dtheta, npix, npix, pixsize, pixsize, pixoffx, pixoffy, blur*pixsize, blur*pixsize, wvl, NULL);
	ccwm(otf, dtf->p[0]);
	cfft2(otf,-1);
	output->p[ipsf]=dnew(npix, npix);
	spmulcreal(output->p[ipsf]->p, si->p[0], otf->p, impst/sumpsf);
	ccellfree(dtf);
	spcellfree(si);
    }else{
	cfft2(otf, -1);
	cfftshift(otf);
	creal2d(&output->p[ipsf], 0, otf, impst);
    }
    char header[500];
    snprintf(header, 500, "%s"
	     "Wavelength: %g\n"
	     "PSF Sampling: %g\"\n"
	     ,msg, wvl, dtheta*206265);
    output->p[ipsf]->header=strdup(header);

    cfree(otf);

}


int main(int argc, char *argv[]){
    /*parameters defined by MAOS simulation */

    const int nwvl=5;
    double wvls[15]={1.908,2.067,2.12,2.173,2.332};
    const double imperr=122.22;

    enum{
	P_OUTFILE=1,
	P_SEED,
	P_EXP,
	P_X,
	P_Y,
	P_NPIX,
	P_PIXSIZE,
	P_PIXOFFX,
	P_PIXOFFY,
	P_BLUR,
	P_TOT
    };

    if(argc!=P_TOT){
	usage();
	warning2("Invalid input\n");
	_exit(1);
    }
    char *outfile=NULL;
    if(!check_suffix(argv[P_OUTFILE],".fits")){
	outfile=stradd(argv[P_OUTFILE],".fits",NULL);
    }else{
	outfile=strdup(argv[P_OUTFILE]);
    }
    int seed=strtol(argv[P_SEED], NULL, 10);
    int exposure=strtol(argv[P_EXP], NULL, 10);
    double thetax=strtod(argv[P_X], NULL);
    double thetay=strtod(argv[P_Y], NULL);
    int npix=strtol(argv[P_NPIX], NULL, 10);
    double pixsize=strtod(argv[P_PIXSIZE], NULL)/206265.;
    double pixoffx=strtod(argv[P_PIXOFFX], NULL);
    double pixoffy=strtod(argv[P_PIXOFFY], NULL);
    double blur=strtod(argv[P_BLUR], NULL);
 
    char msg[400];
    snprintf(msg,400,
	     "%s\n"
	     "Focal plane location: (%g, %g)\"\n"
	     "Pixel size: %g\"\n"
	     "Image offset from center: (%g, %g) of a pixel\n"
	     "Charge diffusion (blur): %g of a pixel\n"
	     , npix?"Detector image":"Raw psf", 
	     thetax, thetay,
	     pixsize*206265, pixoffx, pixoffy, blur
	     );
    info2("%s",msg);
    dcell *psf_lgs=dcellread("PSF%d/evlpsfcl_%d_x%g_y%g.fits", seed, seed, thetax, thetay);
    int nexp=0; /*number of exposure*/
    if(exposure%4!=0 || exposure <0 || exposure>20){
	error("exposure time must be multiple of 4. 0 for all length, 4,8,12,16,20");
    }
    if(exposure>0){
	int i1, i10=(exposure/4-1)*nwvl;
	nexp=1;
	for(i1=0; i1<nexp*nwvl; i1++){
	    psf_lgs->p[i1]=psf_lgs->p[i1+i10];
	}
	psf_lgs->nx=nexp*nwvl;
	psf_lgs->ny=1;
    }
    int npsf=nexp*nwvl;
    double image_size=pixsize*npix;
    double dx, sumpsf;
    int notf;
    if(image_size<wvls[0]/64*512){
	dx=1./16;
	notf=1024;/*size of otf/psf saved by maos*/
	sumpsf=6.34772385106148;
    }else{
	int i1=0;
	dx=1./64;
	notf=4096;/*size of otf/psf saved by maos*/
	sumpsf=6.34772385106148; /*use sumpsf value as 1/16 sampled one because we emphasize the center more.*/
	info("Enlarging PSF\n");
	dcell *psf_large=dcellread("evlpsfcl_4096.fits");
	for(i1=0; i1<npsf; i1++){
	    dmat *dtmp=ddup(psf_large->p[i1%nwvl]);
	    dmat *dtmp2=dsub(psf_lgs->p[i1], 256, 512, 256, 512);
	    dblend(dtmp, dtmp2, 10);
	    dfree(psf_lgs->p[i1]);
	    dfree(dtmp2);
	    psf_lgs->p[i1]=dtmp;
	}
	/*don't change sumpsf after blending.*/
	dcellfree(psf_large);
    }
    dcell *output=dcellnew(nwvl,nexp);
    info2("%d: ", nwvl);
    psfiris_t data={notf, nwvl, dx, sumpsf, npix, pixsize, pixoffx, pixoffy, blur, imperr, wvls, psf_lgs, output, msg};
    thread_t *info=calloc(npsf, sizeof(thread_t));
    thread_prep(info, 0, npsf, npsf, psfiris_do, &data);
    THREAD_POOL_INIT(NCPU);
    CALL_THREAD(info, npsf, 0);
    info2(" done\n");
    dcellwrite(output, "%s", outfile);

    dcellfree(psf_lgs);
    dcellfree(output);
    free(outfile);
    exit_success=1;
}
