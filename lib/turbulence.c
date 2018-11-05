/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include <dirent.h>
#include "../math/mathdef.h"
#include "turbulence.h"
#include "fractal.h"
#include "accphi.h"
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
static char *create_fnatm(GENATM_T *data){
    uint32_t key;
    key=hashlittle(data->rstat, sizeof(rand_t), 0);/*contains seed */
    key=hashlittle(data->wt, sizeof(double)*data->nlayer, key);
    key=hashlittle(&data->r0, sizeof(double), key);
    key=hashlittle(data->L0, sizeof(double)*data->nlayer, key);
    key=hashlittle(&data->dx, sizeof(double), key);
    key=hashlittle(&data->fmin, sizeof(double), key);
    key=hashlittle(&data->fmax, sizeof(double), key);
    key=hashlittle(&data->nx, sizeof(long), key);
    key=hashlittle(&data->ny, sizeof(long), key);
    key=hashlittle(&data->nlayer, sizeof(long), key);
    key=hashlittle(&data->ninit, sizeof(long), key);
    if(data->r0logpsds){
	key=hashlittle(data->r0logpsds->p, sizeof(double)*data->r0logpsds->nx, key);
    }
    char diratm[PATH_MAX];
    snprintf(diratm,PATH_MAX,"%s/atm", CACHE);
    if(!exist(diratm)) mymkdir("%s", diratm);
    char fnatm[PATH_MAX];
    const char *types[]={"vonkarman","fractal","biharmonic"};
    snprintf(fnatm,PATH_MAX,"%s/%s_%ld_%ldx%ld_%g_%ud.bin",
	     diratm,types[data->method],data->nlayer,data->nx,data->ny,data->dx,key);
    long avail=available_space(diratm);
    long need=data->nx*data->ny*data->nlayer*sizeof(double)+500000000;
    if(avail>need){
	return strdup(fnatm);
    }else{
	return NULL;
    }
}

/**
 * Geneate the screens two layer at a time and appends to file if fc is not
 * NULL. Handles large screens well without using the full storage.
 */
static void spect_screen_do(zfarr *fc, GENATM_T *data){
    double slope=0;
    switch(data->method){
    case T_VONKARMAN:
	slope=-11./3.;
	break;
    case T_BIHARMONIC:
	slope=-4.;
	break;
    default:
	error("Invalid call\n");
    }
   
    rand_t *rstat = data->rstat;
    double* wt    = data->wt;
    int nlayer    = data->nlayer;
    dcell *dc     = dcellnew(2,1);
    long nx = data->nx;
    long ny = data->ny;
    double dx=data->dx;
    dc->p[0] = dnew(nx, ny);
    dc->p[1] = dnew(nx, ny);
    double *restrict p1=dc->p[0]->p;
    double *restrict p2=dc->p[1]->p;
    char header[1024];
    double ox=-nx/2*dx;
    double oy=-ny/2*dx;
    snprintf(header, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
	     ox, oy, dx, 0., 0., 0.);
    dc->p[0]->header=strdup(header);
    dc->p[1]->header=strdup(header);
    //For create spatially varying r0
    dmat *spect2=0;
    dcell *dc2=0;
    rand_t rstat2;//don't consume rstat
    if(data->r0logpsds){//Scale r0 across the screen.
	double slope2=data->r0logpsds->p[0];
	double strength=data->r0logpsds->p[1]*(5./3.);//strength of log(wt)
	double minfreq=0, maxfreq=0;
	if(data->r0logpsds->nx>2){//low frequency end is converted to outscale
	    minfreq=data->r0logpsds->p[2];
	}
	if(data->r0logpsds->nx>3){
	    maxfreq=data->r0logpsds->p[3];
	}
	spatial_psd(&spect2, nx, ny, dx, strength, INFINITY, minfreq, maxfreq, slope2, 0.5);
	dc2=dcellnew(2,1);
	dc2->p[0] = dnew(nx, ny);
	dc2->p[1] = dnew(nx, ny);
	seed_rand(&rstat2, rstat->statevec[0]);
	//writebin(spect2, "spect_r0log");
    }
    dmat *prev_screen=0;
    dmat *this_screen=0;
    dmat *prev_scale=0;
    dmat *this_scale=0;
    dmat *spect=0;
    double L0=0;
    for(int ilayer=0; ilayer<nlayer; ilayer++){
	double tk1=myclockd();
	double tk2=tk1;
	if(!prev_screen){
	    double strength=0.0229*pow(data->r0,-5./3.)*pow((0.5e-6)/(2.*M_PI),2);;
	    if(!spect || fabs(data->L0[ilayer]-L0)>EPS){
		L0=data->L0[ilayer];
		spatial_psd(&spect, data->nx,data->ny,data->dx,strength,L0,data->fmin,data->fmax,slope, 0.5);
	    }
	    for(long i=0; i<nx*ny; i++){//don't parallelize this one
		p1[i]=randn(rstat)*spect->p[i];/*real */
		p2[i]=randn(rstat)*spect->p[i];/*imag */
	    }
	    tk2=myclockd();
	    dcell_fft2(dc, -1);
	    this_screen=dc->p[0];
	    if(ilayer+1<nlayer && fabs(data->L0[ilayer+1]-L0)<EPS){//matched L0.
		prev_screen=dc->p[1];
	    }
	}else{
	    this_screen=prev_screen; prev_screen=0;
	}
	dscale(this_screen, sqrt(wt[ilayer]));
	
	if(dc2){
	    if(!prev_scale){
		for(long i=0; i<nx*ny; i++){
		    dc2->p[0]->p[i]=randn(&rstat2)*spect2->p[i];/*real */
		    dc2->p[1]->p[i]=randn(&rstat2)*spect2->p[i];/*imag */
		}
		dcell_fft2(dc2, -1);
		dcwexp(dc2->p[0], 1);
		dcwexp(dc2->p[1], 1);
		this_scale=dref(dc2->p[0]);
		prev_scale=dref(dc2->p[1]);
	    }else{
		this_scale=prev_scale; prev_scale=0;
	    }
	    dcwm(this_screen, this_scale);
	}
	double tk3=myclockd();
	if(fc){/*save to file. */
	    zfarr_dmat(fc, ilayer, this_screen);
	}else{
	    dcp((dmat**)&data->screen->p[ilayer], this_screen);
	}
	double tk4=myclockd();
	info("Layer %d: Randn: %.2f FFT: %.2f %s: %.2f seconds.\n", 
	     ilayer, tk2-tk1, tk3-tk2, fc?"Save":"Copy",tk4-tk3);
    }
    dcellfree(dc);
    dcellfree(dc2);
    dfree(spect2);
    dfree(spect);
}

/**
 * Generate one screen at a time and save to file
 */
static void fractal_screen_do(zfarr *fc, GENATM_T *data){
    const long nx=data->nx;
    const long ny=data->ny;
    if(fc){
	dmat *screen = dnew(data->nx, data->ny);
	for(int ilayer=0; ilayer<data->nlayer; ilayer++){
	    drandn(screen, 1, data->rstat);
	    double r0i=data->r0*pow(data->wt[ilayer], -3./5.);
	    fractal_do(screen, data->dx, r0i, data->L0[ilayer], data->ninit);
	    remove_piston(screen->p, nx*ny);
	    zfarr_dmat(fc, ilayer, screen);
	}
	dfree(screen);
    }else{
	map_t **screen=(map_t**)data->screen->p;
	for(int ilayer=0; ilayer<data->nlayer; ilayer++){
	    drandn((dmat*)screen[ilayer], 1, data->rstat);
	}
	OMPTASK_FOR(ilayer, 0, data->nlayer){
	    double r0i=data->r0*pow(data->wt[ilayer], -3./5.);
	    fractal_do((dmat*)screen[ilayer], screen[0]->dx, r0i, data->L0[ilayer], data->ninit);
	    remove_piston(screen[ilayer]->p, nx*ny);
	}
	OMPTASK_END;
    }
}

/**
 * Generates multiple screens from spectrum. Note that if data->share=1, the
 * atmosphere will be different from data->share=0 due to different algorithms
 * used.
 */
static mapcell* create_screen(GENATM_T *data, void (*atmfun)(zfarr *fc, GENATM_T *data)){
    mapcell* screen;
    long nlayer=data->nlayer;
    char *fnatm=NULL;
    if(data->share){/*shared with file */
	fnatm=create_fnatm(data);
    }
    if(fnatm){
	char fnlock[PATH_MAX];
	snprintf(fnlock, PATH_MAX, "%s.lock", fnatm);
	dcell *in=NULL;
	while(!in){
	    if(exist(fnatm)){
		info("Usinging %s\n", fnatm);
		in=dcellread_mmap("%s",fnatm);
	    }else{
		/*non blocking exclusive lock. */
		int fd=lock_file(fnlock, 0, 0);
		if(fd>=0){/*succeed to lock file. */
		    char fntmp[PATH_MAX];
		    snprintf(fntmp, PATH_MAX, "%s.partial.bin", fnatm);
		    zfarr *fc = zfarr_init(nlayer, 1, "%s", fntmp); 
		    atmfun(fc, data);
		    zfarr_close(fc);
		    if(rename(fntmp, fnatm)){
			error("Unable to rename %s to %s\n", fntmp, fnatm);
		    }
		}else{/*wait for the previous lock to release.*/
		    warning("Waiting for previous lock to release ...");
		    fd=lock_file(fnlock, 1, 0);
		}
		close(fd); remove(fnlock);
	    }
	}
	screen=dcell2map(in);
	dcellfree(in);
	free(fnatm);
    }else{
  	screen=(mapcell*)cellnew(nlayer, 1);
	long nx = data->nx;
	long ny = data->ny;
	double dx = data->dx;
	for(int ilayer=0; ilayer<nlayer; ilayer++){
	    screen->p[ilayer]=mapnew(nx, ny, dx, dx, NULL);
	}
	data->screen=screen;
	atmfun(0, data);
    }
    return screen;
}
/**
 *  Generate vonkarman screens from turbulence statistics.
 */
mapcell* vonkarman_screen(GENATM_T *data){
    data->method=T_VONKARMAN;
    mapcell *screen=create_screen(data, spect_screen_do);
    return screen;
}

/**
 * Generate screens from PSD with power of 12/3 instead of 11/3.
 */
mapcell* biharmonic_screen(GENATM_T *data){
    data->method=T_BIHARMONIC;
    mapcell *screen=create_screen(data, spect_screen_do);
    return screen;
}

/**
 * Generate Fractal screens. Not good statistics.
 */
mapcell *fractal_screen(GENATM_T *data){
    data->method=T_FRACTAL;
    return create_screen(data, fractal_screen_do);
}

map_t *genatm_simple(double r0, double L0, double dx, long nx){
    rand_t rstat;
    seed_rand(&rstat, 1);
    double wt[1]; wt[0]=1;
    GENATM_T cfg={&rstat, wt, r0, &L0, dx, nx, nx, 1, 0, 0, 0, 0, 0, 0};
    mapcell *screens=vonkarman_screen(&cfg);
    map_t *out=mapref(screens->p[0]);
    cellfree(screens);
    return out;
}
/**
   Generate atmosphere and map onto loc.
*/
dmat *genatm_loc(loc_t *loc, double r0, double dsa){
    double D=loc_diam(loc);
    dmat *opd=dnew(loc->nloc, 1);
    map_t *atm=genatm_simple(r0, dsa, loc->dx, ceil(D/loc->dx)*2);
    prop_grid(atm, loc, opd->p, 1, 0,0,1,1,0,0);
    mapfree(atm);
    return opd;
}

/**
 * Compute the covariance for separation of r, and put the values in cov. In
 * kolmogorov spectrum, the variance are defined as half of the structure
 * function between two points separated by rmax.
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
		dbessik(r2pif0, 5./6., &ri, &rk, &rip, &rkp);
		cov->p[i]=vkcoeff*r0f0p*pow(r2pif0, 5./6.)*rk;
	    }
	}
    }
    return cov;
}

/**
   Creates 2-d PSD at size nx*ny: psd=(strength*(f^2+L0^-2)^(slope/2))^power.
   Zero frequency component is in the corner psd->p[0].
 */
void spatial_psd(dmat **pout,  /**<Output*/
		 long nx,      /**<The size*/
		 long ny,      /**<The size*/
		 double dx,    /**<The sampling of spatial coordinate.*/
		 double strength, /**<Strength coefficient*/
		 double outerscale, /**<Outerscale */
		 double minfreq,  /**<Low end frequency cut off*/
		 double maxfreq,  /**<High end frequency cut off*/
		 double slope, /**<should be -11/3 for von karman or kolmogorov
				  screens, or -4 for biharmonic screen (just
				  testing only).*/
		 double power  /**< optionally do a power of psd.*/
    ){
    if(slope==0) slope=-11./3.;//Kolmogorov
    slope*=power/2.;
    if(maxfreq==0) maxfreq=INFINITY;
    const double dfx=1./(nx*dx);
    const double dfy=1./(ny*dx);
    const double zerofreq2=outerscale==0?0:pow(outerscale, -2);
    const double minfreq2=minfreq*minfreq;
    const double maxfreq2=(maxfreq==0?INFINITY:(maxfreq*maxfreq));
    const double scrnstr=pow(strength*(dfx*dfy),power);
    const int nx2=nx/2;
    const int ny2=ny/2;
    if(*pout && ((*pout)->nx!=nx || (*pout)->ny!=ny)){
	dfree(*pout);
    }
    if(!*pout){
	*pout=dnew(nx, ny);
    }
    dmat *psd=*pout;
#pragma omp parallel for
    for(int iy=0;iy<ny;iy++){
	double r2y=pow((iy<ny2?iy:iy-ny)*dfy,2);/* to avoid fft shifting. */
	for(int ix=0;ix<nx;ix++){
	    double r2=pow((ix<nx2?ix:ix-nx)*dfx,2)+r2y;
	    if(r2<=maxfreq2 && r2>=minfreq2){
		psd->p[ix+iy*nx] = pow(r2+zerofreq2,slope)*scrnstr;
	    }
	}
    }
    if(minfreq==0) psd->p[0]=0;  //remove infinite piston mode if minfreq is zero (L0 is inf).
}
dmat* turbpsd(long nx, long ny, double dx, double r0, double L0, double slope, double power){
    double strength=0.0229*pow(r0,-5./3.)*pow((0.5e-6)/(2.*M_PI),2);
    dmat *out=0;
    spatial_psd(&out, nx, ny, dx, strength, L0, 0, INFINITY, slope, power);
    return out;
}
/**
   Estimate anisoplanatic angle theta0 from Fried parameter r0, layer height and
   weights.  */
double calc_aniso(double r0, int nps, double *ht, double *wt){
    double wh=0;
    for(int ips=0; ips<nps; ips++){
	if(wt[ips]>0.01){//only account for positive and significant layers, slodar may give negative results.
	    wh+=pow(fabs(ht[ips]),5./3.)*wt[ips];
	}
    }
    return 0.3144*r0*pow(wh,-3./5.);
}
/**
   Estimate Green wood frequency
*/
double calc_greenwood(double r0, int nps, double *ws, double *wt){
    double wv=0;
    for(int ips=0; ips<nps; ips++){
	wv+=pow(fabs(ws[ips]), 5./3.)*wt[ips];
    }
    return 0.426/r0*pow(wv, 3./5.);
}
/**
   Estimate generalized aniso angle theta2 from Fried parameter r0, and layer
   height and weights, and deformable mirror conjugation heights hc1 hc2 of the
   ground and altitude DMs. */
double calc_aniso2(double r0, int nps, double *ht, double *wt, double hc1, double hc2){
    double wh=0;
    double hh=pow(hc2-hc1,5./3.);
    for(int ips=0; ips<nps; ips++){
	double t1=0.5*pow(fabs(ht[ips]-hc1),5./3.)+0.5*pow(fabs(ht[ips]-hc2),5./3.);
	double t2=-0.25*hh-0.25/hh*pow(pow(fabs(ht[ips]-hc1),5./3.)-pow(fabs(ht[ips]-hc2),5./3.),2);
	wh+=wt[ips]*(t1+t2);
    }
    return 0.3144*r0*pow(wh,-3./5.);
}

