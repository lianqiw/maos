/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>

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


#include <unistd.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <dirent.h>
#include "../math/mathdef.h"
#include "turbulence.h"
#include "fractal.h"
#include "accphi.h"
#include "zernike.h"
#include "petal.h"

/**
 * hash the data to get a unique file name
 */
static char* create_fnatm(genatm_t* data){
	uint32_t key;
	key=hashlittle(data->rstat, sizeof(rand_t), 0);/*contains seed */
	key=hashlittle(data->wt, sizeof(real)*data->nlayer, key);
	key=hashlittle(&data->r0, sizeof(real), key);
	key=hashlittle(data->L0, sizeof(real)*data->nlayer, key);
	key=hashlittle(&data->dx, sizeof(real), key);
	key=hashlittle(&data->fmin, sizeof(real), key);
	key=hashlittle(&data->fmax, sizeof(real), key);
	key=hashlittle(&data->nx, sizeof(long), key);
	key=hashlittle(&data->ny, sizeof(long), key);
	key=hashlittle(&data->nlayer, sizeof(long), key);
	key=hashlittle(&data->ninit, sizeof(long), key);
	if(data->r0logpsds){
		key=hashlittle(P(data->r0logpsds), sizeof(real)*NX(data->r0logpsds), key);
	}
	char fnatm[PATH_MAX];
	snprintf(fnatm, PATH_MAX, "%s/atm", CACHE);
	if(!exist(fnatm)) mymkdir("%s", fnatm);
	long avail=available_space(fnatm);
	long need=NX(data)*NY(data)*data->nlayer*sizeof(real)+500000000;
	if(avail<need){
		return NULL;
	}else{
		const char* prefix=NULL;
		if(fabs(data->slope+11./3.)<EPS){
			prefix="vonkarman";
		} else if(fabs(data->slope+4.)<EPS){
			prefix="biharmonic";
		} else if(data->slope==0){
			prefix="fractal";
		} else{
			prefix="spect";
		}
		snprintf(fnatm, sizeof(fnatm), "%s/atm/%s_%ld_%ldx%ld_%g_%ud.bin",
			CACHE, prefix, data->nlayer, NX(data), NY(data), data->dx, key);
		return strdup(fnatm);
	}
}

/**
 * Geneate the screens two layer at a time and appends to file if fc is not
 * NULL. Handles large screens well without using the full storage.
 */
static void spect_screen_do(zfarr* fc, genatm_t* data){
	real slope=data->slope;
	if(!slope) slope=data->slope=-11./3.;
	rand_t* rstat=data->rstat;
	real* wt=data->wt;
	int nlayer=data->nlayer;
	long nx=data->nx;
	long ny=data->ny;
	real dx=data->dx;
	dcell* dc=dcellnew_same(2, 1, nx, ny);
	real* restrict p1=P(P(dc, 0));
	real* restrict p2=P(P(dc, 1));
	char keywords[1024];
	real ox=-nx/2*dx;
	real oy=-ny/2*dx;
	snprintf(keywords, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
		ox, oy, dx, 0., 0., 0.);
	P(dc, 0)->keywords=strdup(keywords);
	P(dc, 1)->keywords=strdup(keywords);
	//For create spatially varying r0
	dmat* spect2=0;
	dcell* dc2=0;//for scaling.
	rand_t rstat2;//don't consume rstat
	if(data->r0logpsds){//Scale r0 across the screen.
		real slope2=P(data->r0logpsds, 0);
		real strength=P(data->r0logpsds, 1)*(5./3.);//strength of log(wt)
		real minfreq=0, maxfreq=0;
		if(NX(data->r0logpsds)>2){//low frequency end is converted to outscale
			minfreq=P(data->r0logpsds, 2);
		}
		if(NX(data->r0logpsds)>3){
			maxfreq=P(data->r0logpsds, 3);
		}
		spatial_psd(&spect2, nx, ny, dx, strength, INFINITY, minfreq, maxfreq, slope2, 0.5);
		dc2=dcellnew(2, 1);
		P(dc2, 0)=dnew(nx, ny);
		P(dc2, 1)=dnew(nx, ny);
		seed_rand(&rstat2, rstat->statevec[0]);
		//writebin(spect2, "spect_r0log");
	}
	dmat* prev_screen=0;
	dmat* this_screen=0;
	dmat* prev_scale=0;
	dmat* this_scale=0;
	dmat* spect=0;
	if(!fc){//output to screen
		data->screen=dcellnew(nlayer,1);
	}
	real L0=0;
	for(int ilayer=0; ilayer<nlayer; ilayer++){
		real tk1=myclockd();
		real tk2=tk1;
		if(!prev_screen){
			real strength=0.0229*pow(data->r0, -5./3.)*pow((0.5e-6)/(2.*M_PI), 2);;
			if(!spect||fabs(data->L0[ilayer]-L0)>EPS){
				L0=data->L0[ilayer];
				spatial_psd(&spect, NX(data), NY(data), data->dx, strength, L0, data->fmin, data->fmax, slope, 0.5);
			}
			for(long i=0; i<nx*ny; i++){//don't parallelize this one
				p1[i]=randn(rstat)*P(spect, i);/*real */
				p2[i]=randn(rstat)*P(spect, i);/*imag */
			}
			tk2=myclockd();
			dcell_fft2(dc, -1);
			this_screen=P(dc, 0);
			if(ilayer+1<nlayer&&fabs(data->L0[ilayer+1]-L0)<EPS){//matched L0.
				prev_screen=P(dc, 1);
			}
		} else{
			this_screen=prev_screen; prev_screen=0;
		}
		dscale(this_screen, sqrt(wt[ilayer]));

		if(dc2){
			if(!prev_scale){
				for(long i=0; i<nx*ny; i++){
					P(P(dc2, 0), i)=randn(&rstat2)*P(spect2, i);/*real */
					P(P(dc2, 1), i)=randn(&rstat2)*P(spect2, i);/*imag */
				}
				dcell_fft2(dc2, -1);
				dcwexp(P(dc2, 0), 1);
				dcwexp(P(dc2, 1), 1);
				this_scale=dref(P(dc2, 0));
				prev_scale=dref(P(dc2, 1));
			} else{
				this_scale=prev_scale; prev_scale=0;
			}
			dcwm(this_screen, this_scale);
		}
		real tk3=myclockd();
		if(fc){/*save to file. */
			zfarr_push(fc, ilayer, this_screen);
		} else{
			dcp(&P(data->screen, ilayer), this_screen);
		}
		real tk4=myclockd();
		info("Layer %d: Randn: %.2f FFT: %.2f %s: %.2f seconds.\n",
			ilayer, tk2-tk1, tk3-tk2, fc?"Save":"Copy", tk4-tk3);
	}
	dcellfree(dc);
	dcellfree(dc2);
	dfree(spect2);
	dfree(spect);
}

/**
 * Generate one screen at a time and save to file
 */
static void fractal_screen_do(zfarr* fc, genatm_t* data){
	const long nx=NX(data);
	const long ny=NY(data);
	char keywords[1024];
	snprintf(keywords, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
		-data->nx/2*data->dx, -data->ny/2*data->dx, data->dx, 0., 0., 0.);
	if(fc){
		dmat* screen=dnew(NX(data), NY(data));
		screen->keywords=strdup(keywords);
		for(int ilayer=0; ilayer<data->nlayer; ilayer++){
			drandn(screen, 1, data->rstat);
			real r0i=data->r0*pow(data->wt[ilayer], -3./5.);
			fractal_do(screen, data->dx, r0i, data->L0[ilayer], data->ninit);
			//remove_piston(P(screen), nx*ny);
			dadds(screen, -dsum(screen)/(nx*ny));
			zfarr_push(fc, ilayer, screen);
		}
		dfree(screen);
	} else{
		data->screen=dcellnew(data->nlayer, 1);
		dmat** screen=P(data->screen);
		for(int ilayer=0; ilayer<data->nlayer; ilayer++){
			drandn(screen[ilayer], 1, data->rstat);
			screen[ilayer]->keywords=strdup(keywords);
		}
OMP_TASK_FOR(4)
		for(long ilayer=0; ilayer<data->nlayer; ilayer++){
			real r0i=data->r0*pow(data->wt[ilayer], -3./5.);
			fractal_do(screen[ilayer], data->dx, r0i, data->L0[ilayer], data->ninit);
			//remove_piston(P(screen[ilayer]), nx*ny);
			dadds(screen[ilayer], -dsum(screen[ilayer])/(nx*ny));
		}
	}
}
static void genscreen_do(zfarr* fc, genatm_t* data){
	if(fabs(data->slope+1)<EPS){
		fractal_screen_do(fc, data);
	} else{
		spect_screen_do(fc, data);
	}
}
/**
 * Generates multiple screens from spectrum. Note that if data->share=1, the
 * atmosphere will be different from data->share=0 due to different algorithms
 * used.
 */
mapcell* genscreen(genatm_t* data){
	dcell *in=NULL;
	mapcell *screen=NULL;
	long nlayer=data->nlayer;
	char* fnatm=NULL;
	if(data->share){/*shared with file */
		fnatm=create_fnatm(data);
	}
	if(fnatm){
		CACHE_FILE(in, fnatm, ({in=dcellread_mmap("%s", fnatm);}),
					({zfarr*fc=zfarr_init(nlayer, 1, "%s", fnatm);
						genscreen_do(fc, data);
						zfarr_close(fc);
					}),
					{});
	}else{
		genscreen_do(NULL, data);
		in=data->screen; data->screen=0;
	}
	screen=dcell2map(in);
	cellfree(in);
	free(fnatm);fnatm=0;
	return screen;
}
/**
 * Calculate RMS WFE in m from r0 and L0
*/
static real calc_rms(real r0, real L0, real slope){
	real kn2=pow(2*M_PI/0.5e-6, -2);
	if(slope==0) slope=-11./3.;
	real coeff=(slope==-4)?0.0229:0.02748;
	return sqrt(coeff*M_PI*kn2*pow(r0, -5./3.)*pow(L0, -slope-2));
}
/**
 * Calculate r0 from RMS WFE in m and L0
*/
static real calc_r0(real rms, real L0, real slope){
	real kn2=pow(2*M_PI/0.5e-6, -2);
	if(slope==0) slope=-11./3.;
	real coeff=(slope==-4)?0.0229:0.02748;
	return pow(rms*rms/(coeff*M_PI*kn2*pow(L0, -slope-2)), -3./5.);
}
/**
 * Generate screen according to a header. It has multiple possible options.
 * 1) If is a filename. load from file
 * 2) if mode and petal is not set. generate screen from PSD with the following keys
	* L0 (Outer scale)
	* slope (optional, default is -11/3)
	* seed (optional)
 * 3) If mode is set. generate from a zernike mode with the following keys
  	* mode (zernike mode)
 * 4) If petal is set. generate randomized petal modes each with rms value.
 	* petal sets total number of petals
	* theta0 sets petal gap orientation angular offset in radian. 0 means a gap is along y axis.
	* seed
	* piston: generates piston mode if set
	* tip: generate tip mode if set
	* tilt: generates tilt mode if set
 The following common parameters are used:
	* r0 (Fried parameter in m) or rms (in nm)
	* dx (sampling)
	* nx (number of points)
 */
mapcell* genscreen_str(const char* keywords){
	mapcell* surfs=NULL;
	if(zfexist("%s",keywords)){
		info2("Loading surface OPD from %s\n", keywords);
		surfs=mapcellread("%s", keywords);
	} else{
		info2("Generating surface OPD from %s\n", keywords);
		real r0=search_keyword_num(keywords, "r0");
		real rmsnm=search_keyword_num(keywords, "rms");//in nm.
		real mode=search_keyword_num(keywords, "mode");
		real dx=search_keyword_num_default(keywords, "dx", 1./64.);
		real nx=search_keyword_num_default(keywords, "nx", 2048);
		real ny=search_keyword_num_default(keywords, "nx", nx);
		real ht=search_keyword_num_default(keywords, "ht", 0);
		real vx=search_keyword_num_default(keywords, "vx", 0);
		real vy=search_keyword_num_default(keywords, "vy", 0);
		real L0=search_keyword_num_default(keywords, "L0", nx*dx);
		real slope=search_keyword_num_default(keywords, "slope", -11./3.);

		static real seed=0;//avoid using the same seed
		seed=search_keyword_num_default(keywords, "seed", seed+1);
		real petal=search_keyword_num(keywords, "petal");

		if(isnan(r0) && isnan(rmsnm)){
			error("Either r0 or rms should be specified\n");
		}else if(isnan(mode)&&isnan(petal)){//mode is not specified
			if(isnan(r0)){
				r0=calc_r0(rmsnm*1e-9, L0, slope);
				info("r0=%g is computed from rms=%g nm\n", r0, rmsnm);
			}
			info("Generating screen with r0=%g, L0=%g, dx=%g, slope=%g, nx=%g, seed=%g\n",
				r0, L0, dx, slope, nx, seed);
			rand_t rstat;
			seed_rand(&rstat, (int)seed);
			real wt=1;
			genatm_t cfg={&rstat, &wt, r0, &L0, dx, 0, 0, slope, (long)nx, (long)ny, 1, 0, 0, 0, 0};
			surfs=genscreen(&cfg);
		}else{
			if(isnan(rmsnm)){
				rmsnm=calc_rms(r0, L0, slope)*1e9;
				info("rms=%g nm is computed from r0=%g m\n", rmsnm, r0);
			}
			surfs=mapcellnew(1, 1);
			P(surfs, 0)=mapnew(nx, ny, dx, dx);
			if(!isnan(mode)){
				info("Generating screen with zernike mode %d for %g nm RMS.\n", (int)mode, rmsnm);
				loc_t *loc=mksqloc_auto(nx, ny, dx, dx);
				dmat *opd=zernike(loc, nx*dx, 0, 0, -(int)mode);
				dscale(opd, rmsnm*1e-9);
				reshape(opd, nx, ny);
				dcp((dmat**)&P(surfs,0), opd);
				dfree(opd);
				locfree(loc);
			}else if(!isnan(petal)){
				int npetal=(int)search_keyword_num_default(keywords, "npetal", 6);
				info("Generating petal modes for %g nm RMS.\n", rmsnm);
				
				//real dtheta=TWOPI/nseg;
				real theta0=search_keyword_num_default(keywords, "rotdeg", 0)*M_PI/180;
				real cx=search_keyword_num_default(keywords, "cx", nx/2);
				real cy=search_keyword_num_default(keywords, "cy", ny/2);
				dmat *mod=dnew(npetal, 1);
				if(petal<0 && petal+npetal>0){
					P(mod,-(int)petal,0)=1;//a single petal
				}else{//randomized
					rand_t rstat;
					seed_rand(&rstat, (int)seed);
					drandn(mod, 1, &rstat);
				}
				dadds(mod, -dmean(mod));//remove average piston
				dshow(mod,"petal mode");
				dscale(mod, rmsnm*1e-9);
				petal_opd(P(surfs, 0), cx, cy, npetal, theta0, mod);
				dfree(mod);
			}
		}
		for(int i=0; i<PN(surfs); i++){
			if(P(surfs,i)){
				char* old=P(surfs, i)->keywords;
				P(surfs, i)->keywords=stradd(keywords, old, NULL);
				P(surfs,i)->h=ht;
				P(surfs, i)->vx=vx;
				P(surfs, i)->vx=vy;
				if(old) free(old);
				real rmsnmi=sqrt(dsumsq(DMAT(P(surfs,i)))/PN(surfs,i))*1e9;
				if(isnan(rmsnm)){
					rmsnm=calc_rms(r0, L0, slope)*1e9;
				}
				info("Layer %d has %g nm rms wfe, expected %g nm.\n", i, rmsnmi, rmsnm);
			}
		}
	}
	return surfs;
}
/**
   A simpler interface to gerate a single screen.
 */
map_t* genatm_simple(real r0, real L0, real slope, real dx, long nx, int seed){
	rand_t rstat;
	seed_rand(&rstat, seed);
	real wt=1.;
	if(slope==0) slope=-11./3.;
	genatm_t cfg={&rstat, &wt, r0, &L0, dx, 0, 0, slope, nx, nx, 1, 0, 0, 0, 0};
	mapcell* screens=genscreen(&cfg);
	map_t* out=mapref(P(screens, 0));
	cellfree(screens);
	return out;
}
/**
   Generate atmosphere and map onto loc.
*/
dmat* genatm_loc(loc_t* loc, real r0, real dsa, real slope, int seed){
	real D=loc_diam(loc);
	dmat* opd=dnew(loc->nloc, 1);
	map_t* atm=genatm_simple(r0, dsa, slope, loc->dx, ceil(D/loc->dx)*2, seed);
	prop_grid(atm, loc, P(opd), 1, 0, 0, 1, 1, 0, 0);
	mapfree(atm);
	return opd;
}

/**
 * Compute the covariance for separation of r, and put the values in cov. In
 * kolmogorov spectrum, the variance are defined as half of the structure
 * function between two points separated by rmax.
 */

dmat* turbcov(dmat* r, real rmax, real r0, real L0){
	real tg1=tgamma(11./6)*pow(24./5*tgamma(6./5.), 5./6.)
		*pow(2*M_PI/0.5e-6, -2)/pow(M_PI, 8./3.);
	real vkcoeff=tg1/pow(2, 5./6.);
	real vkcoeff0=tg1*tgamma(5./6.)/2;/*for variance */
	dmat* cov=dnew(NX(r), NY(r));
	long n=NX(r)*NY(r);
	if(isinf(L0)){/*kolmogorov. */
		const real power=5./3.;
		real coeff=0.5*6.88*pow(2*M_PI/0.5e-6, -2)*pow(r0, -power);
		real sigma2=coeff*pow(rmax, power);
		for(long i=0; i<n; i++){
			P(cov, i)=sigma2-coeff*pow(P(r, i), power);
		}
	} else{/*von karman. */
		const real f0=1./L0;
		const real r0f0p=pow(r0*f0, -5./3.);
		real ri, rk, rip, rkp;
		real r2pif0;
		for(long i=0; i<n; i++){
			if(fabs(P(r, i))<EPS){
				P(cov, i)=vkcoeff0*r0f0p;
			} else{
				r2pif0=P(r, i)*2*M_PI*f0;
				dbessik(r2pif0, 5./6., &ri, &rk, &rip, &rkp);
				P(cov, i)=vkcoeff*r0f0p*pow(r2pif0, 5./6.)*rk;
			}
		}
	}
	return cov;
}

/**
   Creates 2-d PSD at size nx*ny: psd=(strength*(f^2+L0^-2)^(slope/2))^power.
   Zero frequency component is in the corner P(psd,0).
 */
void spatial_psd(dmat** pout,  /**<Output*/
	long nx,      /**<The size*/
	long ny,      /**<The size*/
	real dx,    /**<The sampling of spatial coordinate.*/
	real strength, /**<Strength coefficient*/
	real outerscale, /**<Outerscale */
	real minfreq,  /**<Low end frequency cut off*/
	real maxfreq,  /**<High end frequency cut off*/
	real slope, /**<should be -11/3 for von karman or kolmogorov
			 screens, or -4 for biharmonic screen (just
			 testing only).*/
	real power  /**< optionally do a power of psd.*/
){
	if(slope==0) slope=-11./3.;//Kolmogorov
	slope*=power/2.;
	if(maxfreq==0) maxfreq=INFINITY;
	const real dfx=1./(nx*dx);
	const real dfy=1./(ny*dx);
	const real zerofreq2=outerscale==0?0:pow(outerscale, -2);
	const real minfreq2=minfreq*minfreq;
	const real maxfreq2=(maxfreq==0?INFINITY:(maxfreq*maxfreq));
	const real scrnstr=pow(strength*(dfx*dfy), power);
	const int nx2=nx/2;
	const int ny2=ny/2;
	if(*pout&&((*pout)->nx!=nx||(*pout)->ny!=ny)){
		dfree(*pout);
	}
	if(!*pout){
		*pout=dnew(nx, ny);
	}
	dmat* psd=*pout;
	OMP_TASK_FOR(4)
	for(int iy=0;iy<ny;iy++){
		real r2y=pow((iy<ny2?iy:iy-ny)*dfy, 2);/* to avoid fft shifting. */
		for(int ix=0;ix<nx;ix++){
			real r2=pow((ix<nx2?ix:ix-nx)*dfx, 2)+r2y;
			if(r2<=maxfreq2&&r2>=minfreq2){
				P(psd, ix, iy)=pow(r2+zerofreq2, slope)*scrnstr;
			}
		}
	}
	if(minfreq==0) P(psd, 0)=0;  //remove infinite piston mode if minfreq is zero (L0 is inf).
}
/**
   Compute spatial PSD of turbulence spectrum.
 */
dmat* turbpsd(long nx, long ny, real dx, real r0, real L0, real slope, real power){
	real strength=0.0229*pow(r0, -5./3.)*pow((0.5e-6)/(2.*M_PI), 2);
	dmat* out=0;
	spatial_psd(&out, nx, ny, dx, strength, L0, 0, INFINITY, slope, power);
	return out;
}
/**
   Estimate anisoplanatic angle theta0 from Fried parameter r0, layer height and
   weights.  */
real calc_aniso(real r0, int nps, real* ht, real* wt){
	real wh=0;
	for(int ips=0; ips<nps; ips++){
		if(wt[ips]>0.01){//only account for positive and significant layers, slodar may give negative results.
			wh+=pow(fabs(ht[ips]), 5./3.)*wt[ips];
		}
	}
	return 0.3144*r0*pow(wh, -3./5.);
}
/**
   Estimate Green wood frequency
*/
real calc_greenwood(real r0, int nps, real* ws, real* wt){
	real wv=0;
	for(int ips=0; ips<nps; ips++){
		wv+=pow(fabs(ws[ips]), 5./3.)*wt[ips];
	}
	return 0.426/r0*pow(wv, 3./5.);
}
/**
   Estimate generalized aniso angle theta2 from Fried parameter r0, and layer
   height and weights, and deformable mirror conjugation heights hc1 hc2 of the
   ground and altitude DMs. */
real calc_aniso2(real r0, int nps, real* ht, real* wt, real hc1, real hc2){
	real wh=0;
	real hh=pow(hc2-hc1, 5./3.);
	for(int ips=0; ips<nps; ips++){
		real t1=0.5*pow(fabs(ht[ips]-hc1), 5./3.)+0.5*pow(fabs(ht[ips]-hc2), 5./3.);
		real t2=-0.25*hh-0.25/hh*pow(pow(fabs(ht[ips]-hc1), 5./3.)-pow(fabs(ht[ips]-hc2), 5./3.), 2);
		wh+=wt[ips]*(t1+t2);
	}
	return 0.3144*r0*pow(wh, -3./5.);
}

