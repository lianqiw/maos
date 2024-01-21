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



#include "mtch.h"

/**
   shift i0 along x or ywithout wraping into i0s1 (+sx/sy) and i0s2 (-sx/sy)
*/
static void mki0sh(real* i0x1, real* i0x2, const dmat* i0, real scale, long sx, long sy){
	int nx=NX(i0);
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
	real sx=cos(theta)*2.*M_PI/NX(i0);
	real sy=sin(theta)*2.*M_PI/NY(i0);
	long ny2=NY(i0)/2;
	long nx2=NX(i0)/2;
	for(long iy=0; iy<NY(i0); iy++){
		for(long ix=0; ix<NX(i0); ix++){
			P(otf, ix, iy)*=-I*((ix<nx2?ix:(ix-NX(i0)))*sx+(iy<ny2?iy:(iy-NY(i0)))*sy);
		}
	}
	cfft2(otf, 1);
	dmat* gx=0;
	creal2d(&gx, 0, otf, 1./(PN(i0)));
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
	real pixrot,    /**<[in] Rotation (CCW, radian) of pixel island. 0 for cartesian*/
	int radgx,      /**<[in] 1: gx/gy is along r/a coord. 0 for cartesian.*/
	int cr          /**<Constraint flag 0: disable, 1: both axis, 2: x only, 3: y only*/
){
	const real* bkgrnd2=dbkgrnd2?P(dbkgrnd2):0;
	const real* bkgrnd2c=dbkgrnd2c?P(dbkgrnd2c):0;
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
	const int i0n=NX(i0)*NY(i0);
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
	addvec(PCOL(i0g, 0), 1, P(gx), i0n, 1, 0);
	addvec(PCOL(i0g, 1), 1, P(gy), i0n, 1, 0);
	addvec(PCOL(i0g, 2), 1, P(i0), i0n, kpx, bkgrnd_res);
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
		if(P(i0,i)<0){//ignore negative pixels.
			P(wt,i)=1./rne2;
		} else if(qe){
			P(wt,i)=P(qe,i)/(rne2/(P(qe,i))+bkgrnd+P(i0,i)+(bkgrnd2?bkgrnd2[i]:0));
		} else{
			P(wt,i)=1./(rne2+bkgrnd+P(i0,i)+(bkgrnd2?bkgrnd2[i]:0));
		}
	}

	dmat* tmp=dpinv(i0g, wt);
	dmat* mtche0=0;
	if(!mtche) mtche=&mtche0;
	dmm(mtche, 0, i0m, tmp, "nn", 1);
	dfree(tmp);

	for(int i=0; i<i0n; i++){/*noise weighting. */
		P(wt,i)=1./P(wt,i);
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

/**
   The routine used to generate matched filter from WFS mean short exposure
   pixel intensities.
 */
void mtch_cell(
	dcell** pmtche,  /**<[out] matched filter. */
	dcell** psanea,  /**<[out] subaperture noise equivalent angle*/
	dmat** pi0sum,   /**<[out] sum of subaperture intensity*/
	dmat** pi0sumsum,/**<[out] sum of all subaperture intensity*/
	const dcell* i0s, /**<Subaperture intensities*/
	const dcell* gxs, /**<Subaperture intensity gradients along x (radial)*/
	const dcell* gys, /**<Subaperture intensity gradients along y (azimuthal)*/
	const dmat* qe, /**<Quantum efficiency of each pixel*/
	const dcell* bkgrnd2, /**<bkgrnd image*/
	const dcell* bkgrnd2c,/**<bkgrnd correction image*/
	real bkgrnd,  /**<bkgrnd per pixel*/
	real bkgrndc, /**<bkgrnd correction per pixel*/
	real rne,     /**<Read out noise per pixel*/
	real pixthetax, /**<Pixel size along x (radial)*/
	real pixthetay, /**<Pixel size along y (azimuthal)*/
	dcell* pixrots,   /**<subaperture pixel rotation.*/
	int radgx,     /**<Leave gradients at radial direction */
	int mtchcr,    /**<constraint. -1: auto*/
	real sigratio /**<scale signal level to increase NEA. default: 1*/
){
	int ni0=NY(i0s);
	const int nsa=NX(i0s);
	if(psanea){
		dcellfree(*psanea);
		*psanea=dcellnew_same(ni0, 1, nsa, 3);
	}
	if(!sigratio) sigratio=1;
	dmat *i0sum=dcellsum_each(i0s);
	real i0max=dmax(i0sum);
	real i0thres=MAX(i0max*0.1, 10*rne);
	if(pi0sum){
		dfree(*pi0sum);
		*pi0sum=dref(i0sum);
	}

	if(pi0sumsum){
		dfree(*pi0sumsum);
		*pi0sumsum=dnew(ni0, 1);
	}

	const int npix=PN(i0s, 0, 0);
	if(pmtche){
		dcellfree(*pmtche);
		*pmtche=dcellnew_same(nsa, ni0, 2, npix);
	}
	real sigratior=1./sigratio;
	for(int ii0=0; ii0<ni0; ii0++){
		real nea2thres=pixthetax*pixthetay*100;
		//P(sanea,ii0)=dnew(nsa, 3);
		real i0sumsum=0;
		int ncrdisable=0;
		dmat* nea2=0;
		for(int isa=0; isa<nsa; isa++){
			real pixrot=pixrots?P(PR(pixrots, ii0, 0), isa):0;//pixel rotation
			int cr=mtchcr;
			if(cr==-1){
				long fwhm=dfwhm_gauss(P(i0s, isa, ii0));
				if(fwhm>4){
					cr=1;
				} else{
					cr=0;
					ncrdisable++;
				}
			}
			dmat* bkgrnd2i=bkgrnd2?PR(bkgrnd2, isa, ii0):NULL;
			dmat* bkgrnd2ci=bkgrnd2c?PR(bkgrnd2c, isa, ii0):NULL;

			mtch(&P(*pmtche, isa, ii0), &nea2, P(i0s, isa, ii0),
				gxs?P(gxs, isa, ii0):0, gys?P(gys, isa, ii0):0, qe,
				bkgrnd2i, bkgrnd2ci, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
				pixrot, radgx, cr);
			if(fabs(sigratio-1)>1e-5){
				if(1){//new simplified scaling
					dcellscale(nea2, 1./sigratio);
				}else{
					dscale(P(i0s, isa, ii0), sigratio);
					if(gxs) dscale(P(gxs, isa, ii0), sigratio);
					if(gys) dscale(P(gys, isa, ii0), sigratio);
					mtch(NULL, &nea2, P(i0s, isa, ii0),
						gxs?P(gxs, isa, ii0):0, gys?P(gys, isa, ii0):0,
						qe, bkgrnd2i, bkgrnd2ci, bkgrnd, bkgrndc, rne, pixthetax, pixthetay,
						pixrot, radgx, cr);
					dscale(P(i0s, isa, ii0), sigratior);
					if(gxs) dscale(P(gxs, isa, ii0), sigratior);
					if(gys) dscale(P(gys, isa, ii0), sigratior);
				}
			}

			i0sumsum+=P(i0sum, isa, ii0);
			if(P(i0sum, isa, ii0)<i0thres||P(nea2, 0)>nea2thres||P(nea2, 3)>nea2thres){
			//Signal level too low or error to high.
				P(nea2, 0)=P(nea2, 3)=nea2thres;
				P(nea2, 1)=P(nea2, 2)=0;
				dset(P(*pmtche, isa, ii0), 0);
			}
			if(psanea){
				P(P(*psanea, ii0), isa, 0)=P(nea2, 0);
				P(P(*psanea, ii0), isa, 1)=P(nea2, 3);
				P(P(*psanea, ii0), isa, 2)=P(nea2, 1);
			}
		}/*isa  */
		dfree(nea2);

		if(mtchcr==-1){
			info("Mtched filter contraint are disabled for %d subaps out of %d.\n", ncrdisable, nsa);
		}
		if(pi0sumsum){
			P(*pi0sumsum, ii0)=i0sumsum;
		}
	}/*ii0 */
	dfree(i0sum);
}