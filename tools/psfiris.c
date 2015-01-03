/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
   \file psfiris.c
   
   Standalone code to sample IRIS PSFs onto detector.
 */
#include <unistd.h>
#include <getopt.h>

#include "../lib/aos.h"

static void usage(){
    info2("Usage: psfiris output.fits prof za ngscount skyp idir npix pixsize("") pixoffx pixoffy blur\n");
}

typedef struct psfiris_t{
    int npix;
    int notf1;
    int notf2;
    double dx1;
    double dx2;
    double pixsize;
    double pixoffx;
    double pixoffy;
    double blur;
    loc_t *ploc;
    dmat *pamp;
    dmat *cc_opd;
    dmat *cc_zero;
    double imperr;
    double *wvls;
    dcell *psf_lgs;
    dcell *output;
    char *msg;
}psfiris_t;

static void psfiris_do(thread_t *info){
    psfiris_t *data=info->data;
    int iwvl=info->start;
    int npix=data->npix;
    int notf1=data->notf1;
    int notf2=data->notf2;
    double dx1=data->dx1;
    double dx2=data->dx2;
    double pixsize=data->pixsize;
    double pixoffx=data->pixoffx;
    double pixoffy=data->pixoffy;
    double blur=data->blur;
    loc_t *ploc=data->ploc;
    dmat *pamp=data->pamp;
    dmat *cc_opd=data->cc_opd;
    dmat *cc_zero=data->cc_zero;
    double imperr=data->imperr;
    double *wvls=data->wvls;
    dcell *psf_lgs=data->psf_lgs;
    dcell *output=data->output;
    char *msg=data->msg;
    cmat *otf=cnew(notf2, notf2);
    //cfft2plan(otf,1);
    //cfft2plan(otf,-1);
    info2("%d ",iwvl);
    /*first create OTF of tt/ps modes on coarse sampling.*/
    cmat *otf0=NULL;
    cmat *otf2=NULL;
    double wvl=wvls[iwvl]*1e-6;
    double dtheta=wvl/(notf1*dx1);/*nyquist sampling*/
    double dtheta2=wvl/(notf2*dx2);
    genotf(&otf2, ploc, pamp, NULL, 0, 0, wvl, dtheta, cc_opd, 0, 0, notf1, notf1, 1, 0);
    genotf(&otf0, ploc, pamp, NULL, 0, 0, wvl, dtheta, cc_zero, 0, 0, notf1, notf1, 1, 0);
    ccwdiv(otf2, otf0, 0);
    dmat *otf3=NULL;
    creal2d(&otf3, 0, otf2, 1);/*otf2 should be real. confirmed.*/
    /*used to up-sample OTF.*/
    map_t *otf_coarse=mapnew(notf1, notf1, dx1/wvl,dx1/wvl, otf3->p); otf3->p=NULL;
    map_t *otf_fine=mapnew(notf2, notf2, dx2/wvl,dx2/wvl, NULL);
    prop_grid_map(otf_coarse, otf_fine, 1, 0, 0, 1, 0, 0, 0);
    mapfree(otf_coarse);
    cfree(otf0);
    cfree(otf2);
    dfree(otf3);
	
    dfftshift((dmat*)otf_fine);/*peak in corner*/
    ccpd(&otf, psf_lgs->p[iwvl]);
    dfree(psf_lgs->p[iwvl]);
    cfftshift(otf);
    cfft2(otf, 1);
    ccwmd(otf, (dmat*)otf_fine, 1); 
    mapfree(otf_fine);
    double sumpsf=otf->p[0];
    double impst=exp(-pow(2*M_PI/wvl*imperr*1e-9,2))/(notf2*notf2);
    if(npix>0){
	dmat *wvlmat=dnew(1,1);
	wvlmat->p[0]=wvl;
	double dxsa=30;//30 meter
	double embfac=wvl/dtheta2/dxsa;
	DTF_T *dtf=mkdtf(wvlmat, dxsa, embfac, notf2, notf2, npix, npix, pixsize, pixsize, pixoffx, pixoffy, blur, NULL, 0, 0);
	ccwm(otf, dtf->nominal->p[0]);
	cfft2(otf,-1);
	output->p[iwvl]=dnew(npix, npix);
	dspmulcreal(output->p[iwvl]->p, dtf->si->p[0], otf->p, impst/sumpsf);
	dtf_free(dtf, 1);
    }else{
	cfft2(otf, -1);
	cfftshift(otf);
	creal2d(&output->p[iwvl], 0, otf, impst);
    }
    char header[500];
    snprintf(header, 500, "%s"
	     "Wavelength: %g\n"
	     "PSF Sampling: %g\"\n"
	     ,msg, wvl, dtheta2*206265);
    output->p[iwvl]->header=strdup(header);

    cfree(otf);

}


int main(int argc, char *argv[]){
    enum{
	P_OUTFILE=1,
	P_PROF,
	P_ZA,
	P_NGSCOUNT,
	P_SKYP,
	P_IDIR,
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
    char *outfile;
    if(!check_suffix(argv[P_OUTFILE],".fits")){
	outfile=stradd(argv[P_OUTFILE],".fits",NULL);
    }else{
	outfile=strdup(argv[P_OUTFILE]);
    }
    
    int prof=strtol(argv[P_PROF], NULL, 10);
    int za=strtol(argv[P_ZA], NULL, 10);
    int ngscount=strtol(argv[P_NGSCOUNT], NULL, 10);
    int skyp=strtol(argv[P_SKYP], NULL, 10);
    int idir=strtol(argv[P_IDIR], NULL, 10);
    int npix=strtol(argv[P_NPIX], NULL, 10);
    double pixsize=strtod(argv[P_PIXSIZE], NULL)/206265.;
    double pixoffx=strtod(argv[P_PIXOFFX], NULL);
    double pixoffy=strtod(argv[P_PIXOFFY], NULL);
    double blur=strtod(argv[P_BLUR], NULL);
    const int nsky=500;
    int iza=0;
    double imperr=0;
    switch(za){
    case 0:
	iza=0;
	imperr=107.56;
	break;
    case 30:
	iza=1;
	imperr=116.63;
	break;
    case 45:
	iza=2;
	imperr=122.22;
	break;
    default:
	error("Invalid za\n");
    }
    int ingscount=0;
    switch(ngscount){
    case 2300: ingscount=0; break;
    case 2800: ingscount=1; break;
    case 3700: ingscount=2; break;
    case 4500: ingscount=3; break;
    case 5500: ingscount=4; break;
    case 6500: ingscount=5; break;
    case 7500: ingscount=6; break;
    default: error("Invalid ngscount\n");
    }
    double thetax[]={0, 8.5, 8.5};
    double thetay[]={0, 0, 8.5};
    char msg[400];
    snprintf(msg,400,
	     "%s\n"
	     "Cn2 Profile:  %d%%\n"
	     "Zenith Angle: %ddeg \n"
	     "J<=19 star count: %d /deg^2\n"
	     "Sky coverage: %d%%\n"
	     "Focal plane location: (%g, %g)\"\n"
	     "Pixel size: %g\"\n"
	     "Image offset from center: (%g, %g) of a pixel\n"
	     "Charge diffusion (blur): %g of a pixel\n"
	     , npix?"Detector image":"Raw psf", 
	     prof, za, ngscount, skyp, thetax[idir], thetay[idir],
	     pixsize*206265, pixoffx, pixoffy, blur
	     );
    info2("%s",msg);
    double tt=0,ps=0;
    if(skyp>0){
	dcell *skycres=dcellread("resfull_%dp.bin", prof);
	if(skycres->nx!=3 || skycres->ny!=7){
	    error("format is wrong\n");
	}
	PDCELL(skycres, pskycres);
	dmat *skyresi=pskycres[ingscount][iza];
	PDMAT(skyresi, pskyresi);
	int ind=(int)round(skyp*(nsky/100.))-1;
	tt=pskyresi[1][ind]+pskyresi[2][ind];
	ps=pskyresi[0][ind]-tt;
	if(ps<0) ps=0;
	dcellfree(skycres);
    }
  
    info2("NGS mode wavefront error: \nTip/tilt:%g nm PS: %g nm\n", sqrt(tt)*1e9, sqrt(ps)*1e9);
    loccell *aloc=loccellread("setup/setup/aloc");
    int naloc=aloc->nx;
    dcell *mode_aloc=dcellread("setup/setup/ahst_Modes");
    int nmod=mode_aloc->p[0]->ny;
    
    const double D=32;
    const double dx1=1.;
    const double dx2=1./16.;
    const int notf1=64;
    const int notf2=1024;/*size of otf/psf saved by maos*/
    
    loc_t *ploc=mksqloc_auto((int)ceil(D/dx1), (int)ceil(D/dx1), dx1, dx1);
    dmat  *pamp=dnew(ploc->nloc, 1);
    locannular(pamp->p, ploc, 0, 0, 15, 1.8, 1);
    dmat *pwt=ddup(pamp); normalize_sum(pwt->p, pwt->nx, 1);
    dmat *mode_ploc=dnew(ploc->nloc, nmod);
    for(int imod=0; imod<nmod; imod++){
	for(int ialoc=0; ialoc<naloc; ialoc++){
	    prop_nongrid_cubic(aloc->p[ialoc], mode_aloc->p[ialoc]->p+imod*mode_aloc->p[ialoc]->nx, 
			       ploc, pamp->p, mode_ploc->p+ploc->nloc*imod, 
			       1, thetax[idir]/206265., thetay[idir]/206265., 1, 0.3, 0, 0);
	}
	double inp=dotdbl(mode_ploc->p+ploc->nloc*imod, mode_ploc->p+ploc->nloc*imod, pwt->p, ploc->nloc);
	dmat *dtmp=dnew_ref(ploc->nloc, 1, mode_ploc->p+ploc->nloc*imod);
	dscale(dtmp, sqrt(1./inp));
	dfree(dtmp);
    }
    dfree(pwt);
    cellfree(aloc);
    dcellfree(mode_aloc);
    dmat *cc_mode=dnew(5,5);
    PDMAT(cc_mode, pcc_mode);
    pcc_mode[0][0]=pcc_mode[1][1]=tt/2.;
    pcc_mode[2][2]=pcc_mode[3][3]=pcc_mode[4][4]=ps/3.;
    int nploc=ploc->nloc;
    dmat *tmp=NULL;
    dmat *cc_opd=NULL;
    dmm(&tmp,0, mode_ploc, cc_mode, "nn", 1);
    dmm(&cc_opd, 0,tmp, mode_ploc, "nt", 1);
    dfree(cc_mode);
    dfree(mode_ploc);
    dfree(tmp);
    const int nwvl=15;
    double wvls[15]={0.928, 0.840, 1.026, 1.092, 0.988, 1.206, 1.270, 1.149, 1.403, 1.629, 1.474, 1.800, 2.182, 1.975, 2.412};
    dmat *cc_zero=dnew(nploc, nploc);    
    
    dcell *psf_lgs=dcellread("za%d_%dp/evlpsfcl_ngsr_1_x%g_y%g.fits", za, prof, thetax[idir], thetay[idir]);
    
    dcell *output=cellnew(nwvl,1);
    info2("%d: ", nwvl);
    psfiris_t data={npix, notf1, notf2, dx1, dx2, pixsize, pixoffx, pixoffy, blur, ploc, pamp, cc_opd, cc_zero, imperr, wvls, psf_lgs, output, msg};
    thread_t *info=calloc(nwvl, sizeof(thread_t));
    thread_prep(info, 0, nwvl, nwvl, psfiris_do, &data);
    THREAD_POOL_INIT(NCPU);
    CALL_THREAD(info, 0);
    writebin(output, "%s", outfile);

    locfree(ploc);
    dfree(pamp);
    dfree(cc_opd);
    dfree(cc_zero);
    dcellfree(psf_lgs);
    dcellfree(output);
    free(outfile);
}
