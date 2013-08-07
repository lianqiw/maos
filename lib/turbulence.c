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

#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include "../sys/sys.h"
#include "turbulence.h"
#include "dmat.h"
#include "cmat.h"
#include "loc.h"
#include "fft.h"
#include "fractal.h"
#include "mathmisc.h"
#include "nr.h"
#include "locbin.h"
#include "cellarr.h"

/**
   Contains routines to generate atmospheric turbulence screens
 */
enum{
    T_VONKARMAN=0,
    T_FRACTAL,
    T_BIHARMONIC
};

/**
 * hash the data to get a unique file name
 */
static char *get_fnatm(GENSCREEN_T *data){
    uint32_t key;
    key=hashlittle(data->rstat, sizeof(rand_t), 0);/*contains seed */
    key=hashlittle(data->wt, sizeof(double)*data->nlayer, key);
    key=hashlittle(&data->dx, sizeof(double), key);
    key=hashlittle(&data->r0, sizeof(double), key);
    key=hashlittle(&data->l0, sizeof(double), key);
    key=hashlittle(&data->nx, sizeof(long), key);
    key=hashlittle(&data->ny, sizeof(long), key);
    key=hashlittle(&data->nlayer, sizeof(long), key);
    key=hashlittle(&data->ninit, sizeof(long), key);

    char diratm[PATH_MAX];
    snprintf(diratm,PATH_MAX,"%s/.aos/atm", HOME);
    if(!exist(diratm)) mymkdir("%s", diratm);
    char fnatm[PATH_MAX];
    char *types[]={"vonkarman","fractal","biharmonic"};
    snprintf(fnatm,PATH_MAX,"%s/maos_%s_%ld_%ldx%ld_%g_%ud.bin",
	     diratm,types[data->method],data->nlayer,data->nx,data->ny,data->dx,key);
    if(zfexist(fnatm)) zftouch(fnatm);
    remove_file_older(diratm, 30*24*3600);
    long avail=available_space(diratm);
    long need=data->nx*data->ny*data->nlayer*sizeof(double)+500000000;
    if(avail>need){
	return strdup(fnatm);
    }else{
	return NULL;
    }
}

/**
 * Geneate the screens sequentially and appends to file if fc is not
 * NULL. Handles large screens well without using the full storage.
 */
static void spect_screen_save(cellarr *fc, GENSCREEN_T *data){
    if(!data->spect){
	info2("Generating spect..."); TIC; tic;
	switch(data->method){
	case T_VONKARMAN:
	    data->spect=turbpsd(data->nx,data->ny,data->dx,data->r0,data->l0,0.5);
	    break;
	case T_BIHARMONIC:
	    data->spect=turbpsd_full(data->nx,data->ny,data->dx,data->r0,data->l0,-2,0.5); 
	    break;
	case T_FRACTAL:
	    break;
	}
	toc2("done");
    }
    rand_t *rstat = data->rstat;
    dmat *spect   = data->spect;
    double* wt    = data->wt;
    int nlayer    = data->nlayer;
    dcell *dc     = dcellnew(2,1);
    long nx = data->nx;
    long ny = data->ny;
    double dx=data->dx;
    dc->p[0] = dnew(nx, ny);
    dc->p[1] = dnew(nx, ny);
    fft_t *fft=dcell_fft2plan(dc, -1, data->nthread);
    double *restrict p1=dc->p[0]->p;
    double *restrict p2=dc->p[1]->p;
    char header[1024];
    double ox=-nx/2*dx;
    double oy=-ny/2*dx;
    snprintf(header, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
	     ox, oy, dx, 0., 0., 0.);
    dc->p[0]->header=strdup(header);
    dc->p[1]->header=strdup(header);
    for(int ilayer=0; ilayer<nlayer; ilayer+=2){
	double tk1=myclockd();
	for(long i=0; i<nx*ny; i++){
	    p1[i]=randn(rstat)*spect->p[i];/*real */
	    p2[i]=randn(rstat)*spect->p[i];/*imag */
	}
	double tk2=myclockd();
	fft2(fft, -1);
	dscale(dc->p[0], sqrt(wt[ilayer]));
	if(ilayer+1<nlayer){
	    dscale(dc->p[1], sqrt(wt[ilayer+1]));
	}
	double tk3=myclockd();
	if(fc){/*save to file. */
	    cellarr_dmat(fc, ilayer, dc->p[0]);
	    if(ilayer+1<nlayer){
		cellarr_dmat(fc, ilayer+1, dc->p[1]);
	    }
	}else{
	    dcp((dmat**)&data->screen[ilayer], dc->p[0]);
	    if(ilayer+1<nlayer){
		dcp((dmat**)&data->screen[ilayer+1], dc->p[1]);
	    }
	}
	double tk4=myclockd();
	info2("Layer %d: Randn: %.2f FFT: %.2f %s: %.2f seconds.\n", 
	      ilayer, tk2-tk1, tk3-tk2, fc?"Save":"Copy",tk4-tk3);
    }
    dcellfree(dc);
    fft_free_plan(fft);
}
/**
 *   Generate turbulence screens all in memory
 */
static void spect_screen_do(GENSCREEN_T *data){
    spect_screen_save(NULL, data);
}
/**
 * Generates multiple screens from spectrum. Note that if data->share=1, the
 * atmosphere will be different from data->share=0 due to different algorithms
 * used.
 */
static map_t** create_screen(GENSCREEN_T *data, 
			     void (*funsave)(cellarr *fc, GENSCREEN_T *data),
			     void (*funmem)(GENSCREEN_T *data)){
    map_t **screen;
    long nlayer=data->nlayer;
    char *fnatm=NULL;
    if(data->share){/*shared with file */
	fnatm=get_fnatm(data);
    }
    if(fnatm){
	dcell *in=NULL;
	while(!in){
	    if(exist(fnatm)){
		info2("Reading %s\n", fnatm);
		in=dcellread_mmap(fnatm);
	    }else{
		char fnlock[PATH_MAX];
		snprintf(fnlock, PATH_MAX, "%s.lock", fnatm);
		/*non blocking exclusive lock. */
		int fd=lock_file(fnlock, 0, 0);
		if(fd>=0){/*succeed to lock file. */
		    char fntmp[PATH_MAX];
		    snprintf(fntmp, PATH_MAX, "%s.partial.bin", fnatm);
		    int disable_save_save=disable_save;
		    disable_save=0;//temporarily disable this feature.
		    cellarr *fc = cellarr_init(nlayer, 1, "%s", fntmp); 
		    funsave(fc, data);
		    disable_save=disable_save_save;
		    cellarr_close(fc);
		    if(rename(fntmp, fnatm)){
			error("Unable to rename %s\n", fnlock);
		    }
		    close(fd); remove(fnlock);
		}else{/*wait for the previous lock to release.*/
		    warning("Waiting for previous lock to release ...");
		    fd=lock_file(fnlock, 1, 0);
		    warning2("OK\n");
		    close(fd); remove(fnlock);
		}
	    }
	}
	int nlayer2;
	screen=dcell2map(&nlayer2, in);
	assert(nlayer==nlayer2);
	dcellfree(in);
	free(fnatm);
    }else{
  	screen=calloc(nlayer,sizeof(map_t*));
	long nx = data->nx;
	long ny = data->ny;
	double dx = data->dx;
	for(int ilayer=0; ilayer<nlayer; ilayer++){
	    screen[ilayer]=mapnew(nx, ny, dx, dx, NULL);
	}
	data->screen=screen;
	funmem(data);
    }
    return screen;
}
/**
 *   Generate vonkarman screens from turbulence statistics.
 */
map_t** vonkarman_screen(GENSCREEN_T *data){
    data->method=T_VONKARMAN;
    map_t **screen=create_screen(data, spect_screen_save, spect_screen_do);
    return(screen);
}

/**
 *  Generate screens from PSD with power of 12/3 instead of 11/3.
 */
map_t** biharmonic_screen(GENSCREEN_T *data){
    data->method=T_BIHARMONIC;
    map_t **screen=create_screen(data, spect_screen_save, spect_screen_do);
    return(screen);
}
/**
 * Generate one screen at a time and save to file
 */
static void fractal_screen_save(cellarr *fc, GENSCREEN_T *data){
    long nx=data->nx;
    long ny=data->ny;
    dmat *dm = dnew(data->nx, data->ny);
    for(int ilayer=0; ilayer<data->nlayer; ilayer++){
	drandn(dm, 1, data->rstat);
	double r0i=data->r0*pow(data->wt[ilayer], -3./5.);
	fractal(dm->p, nx, ny, data->dx, r0i, data->l0, data->ninit);
	remove_piston(dm->p, nx*ny);
	cellarr_dmat(fc, ilayer, dm);
    }
    dfree(dm);
}
static void fractal_screen_thread(GENSCREEN_T *data){
    rand_t *rstat=data->rstat;
    map_t** screen=data->screen;
    const double *wt=data->wt;
    long nx=screen[0]->nx;
    long ny=screen[0]->ny;
 repeat:
    LOCK(data->mutex_ilayer);
    int ilayer=data->ilayer;
    data->ilayer++;
    if(ilayer>=data->nlayer){
	UNLOCK(data->mutex_ilayer);
	return;
    }
    drandn((dmat*)screen[ilayer], 1, rstat);
    UNLOCK(data->mutex_ilayer);
    double r0i=data->r0*pow(wt[ilayer], -3./5.);
    /*info("r0i=%g\n", r0i); */
    fractal(screen[ilayer]->p, nx, ny, screen[0]->dx, r0i, data->l0, data->ninit);
    remove_piston(screen[ilayer]->p, nx*ny);
    goto repeat;
}
static void fractal_screen_do(GENSCREEN_T *data){
    PINIT(data->mutex_ilayer);
    CALL(fractal_screen_thread, data, data->nthread,1);
}

/**
 * Generate Fractal screens. Not good statistics.
 */

map_t **fractal_screen(GENSCREEN_T *data){
    data->method=T_FRACTAL;
    return create_screen(data, fractal_screen_save, fractal_screen_do);
}

/**
 *  Compute the covariance for separation of r, and put the values in cov. In
 *  kolmogorov spectrum, the variance are defined as half of the structure
 *  function between two points separated by rmax.
 */

dmat* turbcov(dmat *r, double rmax, double r0, double L0){
    double tg1=tgamma(11./6) * pow(24./5 * tgamma(6./5.), 5./6.)
	* pow(2 * M_PI/0.5e-6, -2) / pow(M_PI, 8./3.);
    double vkcoeff  = tg1 / pow(2, 5./6.);    
    double vkcoeff0 = tg1 * tgamma(5./6.) / 2 ;/*for variance */
    dmat *cov=dnew(r->nx, r->ny);
    long n=r->nx*r->ny;
    if(!isfinite(L0)){/*kolmogorov. */
	const double power=5./3.;
	double coeff=6.88*pow(2*M_PI/0.5e-6, -2) * pow(r0, -power);
	double sigma2=0.5*coeff*pow(rmax, power);
	for(long i=0; i<n; i++){
	    cov->p[i]=sigma2-0.5*coeff*pow(r->p[i], power);
	}
    }else{/*von karman. */
	const double f0=1./L0;
	const double r0f0p=pow(r0*f0, -5./3.);
	double ri, rk, rip, rkp;
	double r2pif0;	
	for(long i=0; i<n; i++){
	    if(fabs(r->p[i])<EPS){
		cov->p[i]=vkcoeff0*r0f0p;
	    }else{
		r2pif0=r->p[i]*2*M_PI*f0;
		bessik(r2pif0, 5./6., &ri, &rk, &rip, &rkp);
		cov->p[i]=vkcoeff*r0f0p*pow(r2pif0, 5./6.)*rk;
	    }
	}
    }
    return cov;
}

/**
 * Compute the turbulence spectrum at size nx*ny, with spacing dx. Notice that
 * the zero frequency component is in the corner psd->p[0].
 */

dmat *turbpsd_full(long nx,      /**<The size*/
		   long ny,      /**<The size*/
		   double dx,    /**<The sampling of spatial coordinate.*/
		   double r0,    /**<The Fried parameter*/
		   double L0,    /**<The outer scale*/
		   double slope, /**<should be -11/6 for von karman or kolmogorov
				    screens, or -2 for biharmonic screen (just
				    testing only).*/
		   double power  /**< optionally do a power of psd.*/
		   ){
    if(nx & 1 || ny & 1){
	warning("Screen is odd size.");
    }
    slope*=power;
    const double dfx=1./(nx*dx);
    const double dfy=1./(ny*dx);
    const double dfx2=dfx*dfx;
    const double L02=pow(L0,-2);
    const double scrnstr=pow(0.0229*pow(r0,-5./3.)*pow((0.5e-6)/(2.*M_PI),2)*(dfx*dfy),power);
    const int nx2=nx/2;
    const int ny2=ny/2;
    dmat *psd=dnew(nx,ny);
    for(int i=0;i<ny;i++){
	double r2y=pow((i<ny2?i:i-ny)*dfy,2);/* to avoid fft shifting. */
	double *psd1=psd->p+i*nx;
	for(int j=0;j<nx2;j++){
	    double r2x=j*j*dfx2;
	    psd1[j] = pow(r2x+r2y+L02,slope)*scrnstr;
	}
	for(int j=nx2;j<nx;j++){
	    double r2x=(j-nx)*(j-nx)*dfx2;
	    psd1[j] = pow(r2x+r2y+L02,slope)*scrnstr;
	}
    }
    if(!isfinite(L0)) psd->p[0]=0;  //remove infinite piston mode if l0 is infinity.
    return psd;
}

