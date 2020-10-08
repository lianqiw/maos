/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

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

#include "mtch.h"

/**
   shift i0 along x or ywithout wraping into i0s1 (+sx/sy) and i0s2 (-sx/sy)
*/
static void mki0sh(real* i0x1, real* i0x2, const dmat* i0, real scale, long sx, long sy){
	int nx=i0->nx;
	typedef real pcol[nx];
	pcol* i0x1p=(pcol*)i0x1;
	pcol* i0x2p=(pcol*)i0x2;
	for(int iy=0; iy<i0->ny-sy; iy++){
		for(int ix=0; ix<i0->nx-sx; ix++){
			i0x1p[iy+sy][ix+sx]=P(i0, ix, iy)*scale;
			i0x2p[iy][ix]=P(i0, ix+sx, iy+sy)*scale;
		}
	}
}

/**
   Compute the derivative using FFT method along direction specified by theta (radian).
 */
dmat* derive_by_fft(const dmat* i0, real theta){
	cmat* otf=0;
	ccpd(&otf, i0);
	cfft2(otf, -1);
	real sx=cos(theta)*2.*M_PI/i0->nx;
	real sy=sin(theta)*2.*M_PI/i0->ny;
	long ny2=i0->ny/2;
	long nx2=i0->nx/2;
	for(long iy=0; iy<i0->ny; iy++){
		for(long ix=0; ix<i0->nx; ix++){
			P(otf, ix, iy)*=-I*((ix<nx2?ix:(ix-i0->nx))*sx+(iy<ny2?iy:(iy-i0->ny))*sy);
		}
	}
	cfft2(otf, 1);
	dmat* gx=0;
	creal2d(&gx, 0, otf, 1./(i0->nx*i0->ny));
	cfree(otf);
	return gx;
}

/**
   add two vectors: out=out*alpha+in*beta+theta
*/
void addvec(real* restrict out, real alpha,
	const real* restrict in, int N, real beta, real theta){

	if(fabs(alpha)<EPS){
		memset(out, 0, sizeof(real)*N);
		if(in){
			for(int i=0; i<N; i++){
				out[i]=in[i]*beta+theta;
			}
		} else{
			for(int i=0; i<N; i++){
				out[i]=theta;
			}
		}
	} else{
		if(in){
			for(int i=0; i<N; i++){
				out[i]=out[i]*alpha+in[i]*beta+theta;
			}
		} else{
			for(int i=0; i<N; i++){
				out[i]=out[i]*alpha+theta;
			}
		}
	}
}


/**
   Generating matched filter from averaged short exposure images.
*/
void mtch(dmat** mtche,   /**<[out] the matched filter*/
	dmat** neaout,  /**<[out] the subaperture noise equivalent angle*/
	const dmat* i0, /**<[in] Averaged subaperture image*/
	const dmat* gx, /**<[in] derivative of i0 along x (r)*/
	const dmat* gy, /**<[in] derivative of i0 along y (a)*/
	const dmat* qe, /**<[in] non uniform quantum efficiency (optional)*/
	const dmat* dbkgrnd2,  /**<[in] background*/
	const dmat* dbkgrnd2c, /**<[in] background calibration*/
	real bkgrnd,    /**<[in] global background*/
	real bkgrndc,   /**<[in] global background calibration*/
	real rne,       /**<[in] Detector read noise*/
	real pixthetax, /**<[in] Size of pixel along x*/
	real pixthetay, /**<[in] Size of pixel along y*/
	real pixrot,    /**<[in] Rotation (CCW, radian) of pixel island 0 for cartesian*/
	int radgx,      /**<[in] 1: gx/gy is along r/a coord.*/
	int cr          /**<Constraint flag 0: disable, 1: both axis, 2: x only, 3: y only*/
){
	const real* bkgrnd2=dbkgrnd2?dbkgrnd2->p:0;
	const real* bkgrnd2c=dbkgrnd2c?dbkgrnd2c->p:0;
	const real bkgrnd_res=bkgrnd-bkgrndc;
	const real kpx=1./pixthetax;
	const real kpy=1./pixthetay;
	int nmod=3;
	int mtchcrx=0;
	int mtchcry=0;

	if(cr==1||cr==2){
		mtchcrx=nmod;
		nmod+=2;
	}
	if(cr==1||cr==3){
		mtchcry=nmod;
		nmod+=2;
	}
	const int i0n=i0->nx*i0->ny;
	dmat* i0m=dnew(2, nmod);
	dmat* i0g=dnew(i0n, nmod);
	dmat* wt=dnew(i0n, 1);
	/*Derivative is along r/a or x/y*/
	P(i0m, 0, 0)=1;
	P(i0m, 1, 1)=1;
	real theta=radgx?0:pixrot;
	if(mtchcrx){/*constrained x(radial) */
		P(i0m, 0, mtchcrx)=cos(theta);
		P(i0m, 1, mtchcrx)=sin(theta);
		P(i0m, 0, mtchcrx+1)=-P(i0m, 0, mtchcrx);
		P(i0m, 1, mtchcrx+1)=-P(i0m, 1, mtchcrx);
	}
	if(mtchcry){/*constrained y(azimuthal). */
		P(i0m, 0, mtchcry)=-sin(theta);
		P(i0m, 1, mtchcry)=cos(theta);
		P(i0m, 0, mtchcry+1)=-P(i0m, 0, mtchcry);
		P(i0m, 1, mtchcry+1)=-P(i0m, 1, mtchcry);
	}
	dmat* gx2=0, * gy2=0;
	if(!gx){
		warning_once("Compute derivative using FFT\n");
		gx=gx2=derive_by_fft(i0, theta); dscale(gx2, kpx);
		gy=gy2=derive_by_fft(i0, theta+M_PI/2); dscale(gy2, kpy);
	}
	addvec(PCOL(i0g, 0), 1, gx->p, i0n, 1, 0);
	addvec(PCOL(i0g, 1), 1, gy->p, i0n, 1, 0);
	addvec(PCOL(i0g, 2), 1, i0->p, i0n, kpx, bkgrnd_res);
	addvec(PCOL(i0g, 2), 1, bkgrnd2, i0n, 1, bkgrnd_res);
	addvec(PCOL(i0g, 2), 1, bkgrnd2c, i0n, -1, 0);/*subtract calibration */
	if(mtchcrx){
		mki0sh(PCOL(i0g, mtchcrx), PCOL(i0g, mtchcrx+1), i0, kpx, 1, 0);
		addvec(PCOL(i0g, mtchcrx), 1, bkgrnd2, i0n, 1, bkgrnd_res);
		addvec(PCOL(i0g, mtchcrx), 1, bkgrnd2c, i0n, -1, 0);
		addvec(PCOL(i0g, mtchcrx+1), 1, bkgrnd2, i0n, 1, bkgrnd_res);
		addvec(PCOL(i0g, mtchcrx+1), 1, bkgrnd2c, i0n, -1, 0);
	}
	if(mtchcry){
		mki0sh(PCOL(i0g, mtchcry), PCOL(i0g, mtchcry+1), i0, kpy, 0, 1);
		addvec(PCOL(i0g, mtchcry), 1, bkgrnd2, i0n, 1, bkgrnd_res);
		addvec(PCOL(i0g, mtchcry), 1, bkgrnd2c, i0n, -1, 0);
		addvec(PCOL(i0g, mtchcry+1), 1, bkgrnd2, i0n, 1, bkgrnd_res);
		addvec(PCOL(i0g, mtchcry+1), 1, bkgrnd2c, i0n, -1, 0);
	}

	/*adding rayleigh backscatter poisson noise. */
	real rne2=rne*rne;
	for(int i=0; i<i0n; i++){/*noise weighting. */
		if(i0->p[i]<0){//ignore negative pixels.
			wt->p[i]=1./rne2;
		} else if(qe){
			wt->p[i]=qe->p[i]/(rne2/(qe->p[i])+bkgrnd+i0->p[i]+(bkgrnd2?bkgrnd2[i]:0));
		} else{
			wt->p[i]=1./(rne2+bkgrnd+i0->p[i]+(bkgrnd2?bkgrnd2[i]:0));
		}
	}

	dmat* tmp=dpinv(i0g, wt);
	dmat* mtche0=0;
	if(!mtche) mtche=&mtche0;
	dmm(mtche, 0, i0m, tmp, "nn", 1);
	dfree(tmp);

	for(int i=0; i<i0n; i++){/*noise weighting. */
		wt->p[i]=1./wt->p[i];
	}
	dmat* nea2=dtmcc(*mtche, wt);

	if(radgx&&pixrot){
	//Rotate mtched filter to x/y
		drotvect(*mtche, pixrot);
		//Rotate NEA to (x/y)
		if(neaout) drotvecnn(neaout, nea2, pixrot);
	} else{//Already in x/y
		if(neaout) dcp(neaout, nea2);
	}
	dfree(nea2);
	dfree(i0m);
	dfree(i0g);
	dfree(wt);
	dfree(gx2);
	dfree(gy2);
	if(mtche0) dfree(mtche0);
}
/**
   A simplified wrapper for mtch
*/
void mtch2(dmat** mtche, dmat** nea, const dmat* i0, const dmat* gx, const dmat* gy, int cr){
	mtch(mtche, nea, i0, gx, gy, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, cr);
}
