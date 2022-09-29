/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file psfgc.c
   Standalone code to sample IRIS PSFs onto detector for the Galactic center simulations.
 */

#include <getopt.h>
#include <unistd.h>
#include "../lib/aos.h"

static void usage(){
	info("Usage: psfgc output.fits seed exposure x y npix pixsize pixoffx pixoffy blur\n"
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
	double* wvls;
	dcell* psf_lgs;
	dcell* output;
	char* msg;
}psfiris_t;

static void* psfiris_do(thread_t* info){
	psfiris_t* data=(psfiris_t*)info->data;
	int ipsf=info->start;
	int nwvl=data->nwvl;
	int iwvl=ipsf%nwvl;
	int notf=data->notf;
	int npix=data->npix;
	double sumpsf=data->sumpsf;
	double dx=data->dx;
	double pixsize=data->pixsize;
	dmat* pixoffx=dnew(1, 1); P(pixoffx,0)=data->pixoffx;
	dmat* pixoffy=dnew(1, 1); P(pixoffy,0)=data->pixoffy;
	double blur=data->blur;
	double imperr=data->imperr;
	double* wvls=data->wvls;
	dcell* psf_lgs=data->psf_lgs;
	dcell* output=data->output;
	char* msg=data->msg;
	int ncomp=P(psf_lgs,0)->nx;
	cmat* otf=cnew(ncomp, ncomp);
	//cfft2plan(otf,1);
	//cfft2plan(otf,-1);
	info("%d ", ipsf);
	/*first create OTF of tt/ps modes on coarse sampling.*/
	double wvl=wvls[iwvl];
	double dtheta=wvl/(notf*dx);
	cembedd(otf, P(psf_lgs,ipsf), 0);
	dfree(P(psf_lgs,ipsf));
	cfftshift(otf);
	cfft2(otf, 1);
	double impst=exp(-pow(2*M_PI/wvl*imperr*1e-9, 2))/(ncomp*ncomp);
	if(npix>0){
		dmat* wvlmat=dnew(1, 1);
		P(wvlmat,0)=wvl;
		double dxsa=30;//30 meter
		double embfac=wvl/dtheta/dxsa;
		dbg("embfac=%g\n", embfac);
		dtf_t* dtf=mkdtf(wvlmat, dxsa, embfac, ncomp, ncomp, npix, npix, pixsize, pixsize, pixoffx, pixoffy, blur, NULL);
		ccwm(otf, P(dtf->nominal,0));
		cfft2(otf, -1);
		P(output,ipsf)=dnew(npix, npix);
		dspmulcreal(P(P(output,ipsf)), P(dtf->si,0), P(otf), impst/sumpsf);
		dtf_free(dtf);
	} else{
		cfft2(otf, -1);
		cfftshift(otf);
		creal2d(&P(output,ipsf), 0, otf, impst);
	}
	char header[500];
	snprintf(header, 500, "%s"
		"Wavelength: %g\n"
		"PSF Sampling: %g\"\n"
		, msg, wvl, dtheta*RAD2AS);
	P(output,ipsf)->header=strdup(header);

	cfree(otf);
	dfree(pixoffx);
	dfree(pixoffy);
	return NULL;
}


int main(int argc, char* argv[]){
	/*parameters defined by MAOS simulation */

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
		warning("Invalid input\n");
		_exit(1);
	}
	char* outfile=NULL;
	if(!check_suffix(argv[P_OUTFILE], ".fits")){
		outfile=stradd(argv[P_OUTFILE], ".fits", NULL);
	} else{
		outfile=strdup(argv[P_OUTFILE]);
	}
	int seed=strtol(argv[P_SEED], NULL, 10);
	double exposure=strtod(argv[P_EXP], NULL);
	double thetax=strtod(argv[P_X], NULL);
	double thetay=strtod(argv[P_Y], NULL);
	int npix=strtol(argv[P_NPIX], NULL, 10);
	double pixsize=strtod(argv[P_PIXSIZE], NULL)*AS2RAD;
	double pixoffx=strtod(argv[P_PIXOFFX], NULL);
	double pixoffy=strtod(argv[P_PIXOFFY], NULL);
	double blur=strtod(argv[P_BLUR], NULL);

	char msg[400];
	snprintf(msg, 400,
		"%s\n"
		"Focal plane location: (%g, %g)\"\n"
		"Pixel size: %g\"\n"
		"Image offset from center: (%g, %g) of a pixel\n"
		"Charge diffusion (blur): %g of a pixel\n"
		, npix?"Detector image":"Raw psf",
		thetax, thetay,
		pixsize*RAD2AS, pixoffx, pixoffy, blur
	);
	info("%s", msg);
	dcell* psf_lgs=NULL;
	if(zfexist("evlpsfcl_%d_x%g_y%g.fits", seed, thetax, thetay)){
		psf_lgs=dcellread("evlpsfcl_%d_x%g_y%g.fits", seed, thetax, thetay);
	} else if(zfexist("PSF%d/evlpsfcl_%d_x%g_y%g.fits", seed, seed, thetax, thetay)){
		psf_lgs=dcellread("PSF%d/evlpsfcl_%d_x%g_y%g.fits", seed, seed, thetax, thetay);
	}
	if(psf_lgs->nx*psf_lgs->nx==0){
		warning("Array is empty, nothing to do\n"); exit(0);
	}
	int nwvl=0; /*number of wavelength*/
	int nexp=0; /*number of exposure*/
	double* wvls=NULL;
	if(psf_lgs->ny!=1){
		nwvl=psf_lgs->nx;
		nexp=psf_lgs->ny;
	} else{
		for(int i=0; i<psf_lgs->nx*psf_lgs->ny; i++){
			double wvl=search_header_num_valid(P(psf_lgs,i)->header, "Wavelength");
			int iwvl;
			for(iwvl=0; iwvl<nwvl; iwvl++){
				if(fabs(wvls[iwvl]-wvl)<1e-14){
					break;
				}
			}
			if(iwvl==nwvl){//not found.
				nwvl++;
				wvls=myrealloc(wvls, nwvl, double);
				wvls[nwvl-1]=wvl;
			}
		}
		nexp=psf_lgs->nx/nwvl;
		if(psf_lgs->nx!=nexp*nwvl){
			error("nwvl=%d, nexp=%d, not matching ntot=%ld\n", nwvl, nexp, psf_lgs->nx);
		}
	}
	double texp=search_header_num_valid(P(psf_lgs,0)->header, "Exposure");
	double tmp;
	if(fabs(modf(exposure/texp, &tmp))>1e-5||exposure <0||exposure>texp*nexp){
		error("exposure time must be multiple of %g and less than or equal to %g. 0 for all length\n",
			texp, texp*nexp);
	}
	if(exposure>0){
		int i1;
		int i10=(int)round(exposure/texp-1)*nwvl;
		nexp=1;
		for(i1=0; i1<nexp*nwvl; i1++){
			P(psf_lgs,i1)=P(psf_lgs,i1+i10);
		}
		psf_lgs->nx=nexp*nwvl;
		psf_lgs->ny=1;
	}
	int npsf=nexp*nwvl;
	int psfsize1=P(psf_lgs,0)->nx;
	double image_size=pixsize*npix;
	int notf=(int)search_header_num_valid(P(psf_lgs,0)->header, "FFT Grid");
	double sumpsf=search_header_num_valid(P(psf_lgs,0)->header, "PSF Sum to");
	double dx=search_header_num_valid(P(psf_lgs,0)->header, "OPD Sampling");
	double dtheta1=search_header_num_valid(P(psf_lgs,0)->header, "PSF Sampling")*AS2RAD;
	int psfsizevalid=MIN(psfsize1, notf/2);/*valid psf range*/
	if(image_size>dtheta1*psfsizevalid){/*need blending*/
		dbg("Enlarging PSF\n");
		dcell* psf_large=dcellread("evlpsfcl_4096.fits");
		dx=search_header_num_valid(P(psf_large,0)->header, "OPD Sampling");
		notf=(int)search_header_num_valid(P(psf_large,0)->header, "FFT Grid");
		if(psf_large->nx!=nwvl){
			error("PSF large has incorrect dimension\n");
		}
		for(int i1=0; i1<npsf; i1++){
			double wvl2=search_header_num_valid(P(psf_large,i1%nwvl)->header, "Wavelength");
			if(fabs(wvl2-wvls[i1%nwvl])>1e-10){
				error("Wavelenght mismatch.\n");
			}
			dmat* dtmp=ddup(P(psf_large,i1%nwvl));
			dmat* dtmp2=dsub(P(psf_lgs,i1),
				(psfsize1-psfsizevalid)/2, psfsizevalid,
				(psfsize1-psfsizevalid)/2, psfsizevalid);
			dblend(dtmp, dtmp2, 10);
			dfree(P(psf_lgs,i1));
			dfree(dtmp2);
			P(psf_lgs,i1)=dtmp;
		}
		/*don't change sumpsf after blending.*/
		dcellfree(psf_large);
	}
	dcell* output=dcellnew(nwvl, nexp);
	info("%d: ", nwvl);
	psfiris_t data={notf, nwvl, dx, sumpsf, npix, pixsize, pixoffx, pixoffy, blur, imperr, wvls, psf_lgs, output, msg};
	thread_t* tdata=thread_prep(0, npsf, npsf, psfiris_do, &data);
	THREAD_POOL_INIT(NCPU);
	CALL_THREAD(tdata, 0);
	info(" done\n");
	writebin(output, "%s", outfile);
	free(tdata);
	dcellfree(psf_lgs);
	dcellfree(output);
	free(outfile);
	free(wvls);
}
