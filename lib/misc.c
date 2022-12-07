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
   Miscellaneous routines.
*/
#include "misc.h"
#include "cure.h"
/**
   add photon and read out noise.  pcalib part of bkgrnd is calibrated
   out. pcalib2 part of bkgrnd2 is calibrated out.  */
void addnoise(dmat* A,              /**<[in/out]The pixel intensity array*/
	rand_t* rstat,        /**<[in]The random stream*/
	const real bkgrnd,  /**<[in]Real background in PDEs per pixel per frame*/
	const real bkgrndc, /**<[in]Removed background in PDEs per pixel per frame*/
	const dmat* bkgrnd2,  /**<[in]Real background in PDEs of each pixel per frame.*/
	const dmat* bkgrnd2c, /**<[in]Removed background in PDEs of each pixel per frame.*/
	const dmat* qe,       /**<[in]Pixel dependent Quantum Efficiency*/
	const real rne,     /**<[in]Read out noise per pixel per read*/
	real excess   /**<[in]Excess noise factor*/
){
	long np=NX(A)*NY(A);
	assert(!bkgrnd2||NX(bkgrnd2)*NY(bkgrnd2)==np);
	assert(!bkgrnd2c||NX(bkgrnd2c)*NY(bkgrnd2c)==np);
	if(excess<1) excess=1;
	for(int ix=0; ix<np; ix++){
		real tot=P(A,ix)+bkgrnd+(bkgrnd2?P(bkgrnd2,ix):0);
		real corr=bkgrndc+(bkgrnd2c?P(bkgrnd2c,ix):0);
		real scale=1;
		if(qe){//the second qe factor is flat-field correction.
			tot*=P(qe,ix);
			scale=1./P(qe,ix);
		}
		P(A,ix)=(randp(rstat, tot*excess)+tot*(1.-excess)+rne*randn(rstat))/scale-corr;
	}
}
/**
   Add noise to gradients according to neal, which is LL' decomposition of the sanea
 */
void addnoise_grad(dmat* grad, const dmat* neal, rand_t* srand){
	const long nsa=NX(neal);
	assert(nsa*2==NX(grad)*NY(grad));
	const real* neax=PCOL(neal, 0);
	const real* neay=PCOL(neal, 1);
	const real* neaxy=PCOL(neal, 2);
	real* restrict ggx=P(grad);
	real* restrict ggy=ggx+nsa;
	for(long isa=0; isa<nsa; isa++){
	/*Preserve the random sequence. */
		real n1=randn(srand);
		real n2=randn(srand);
		real errx=neax[isa]*n1;
		real erry=neay[isa]*n2+neaxy[isa]*n1;/*cross term. */
		ggx[isa]+=errx;
		ggy[isa]+=erry;
	}
}
/**
   Determine the CoG on multiple locations near the nominal position.
*/
int cog_multi(
	dmat** cg,      /**<[out]n*3 The determined position and total intensity.*/
	const dmat* im, /**<[in] 2-d array contains the image*/
	const dmat* loc,/**<[in] n*2 array. contains the nominal position of each peak*/
	int wsize       /**<[in] window size*/
){
	if(!*cg){
		*cg=dnew(NX(loc), 3);
	}
	dmat* pcg=*cg;
	int wsize2=(wsize-1)/2;
	for(int iloc=0; iloc<NX(loc); iloc++){
		int cx=P(loc, iloc, 0);
		int cy=P(loc, iloc, 1);
		int cx2=cx, cy2=cy;
		real  cmax=-INFINITY, cmin=INFINITY;
		for(int jy=cy-wsize2; jy<=cy+wsize2; jy++){
			for(int jx=cx-wsize2; jx<=cx+wsize2; jx++){
				if(P(im, jx, jy)>cmax){
					cmax=P(im, jx, jy);
					cx2=jx;
					cy2=jy;
				}
				if(P(im, jx, jy)<cmin){
					cmin=P(im, jx, jy);
				}
			}
		}
		//Invalid subimage
		if(abs(cx2-cx)>wsize2||abs(cy2-cy)>wsize2){
			P(pcg, iloc, 0)=NAN;
			P(pcg, iloc, 1)=NAN;
			continue;
		}
		real thres=cmax*0.2+cmin*0.8;
		real sum=0, sumx=0, sumy=0;
		for(int jy=cy2-wsize2; jy<=cy2+wsize2; jy++){
			for(int jx=cx2-wsize2; jx<=cx2+wsize2; jx++){
				real val=P(im, jx, jy);
				if(val>thres){
					sum+=val;
					sumx+=val*jx;
					sumy+=val*jy;
				}
			}
		}
		P(pcg, iloc, 2)=sum;
		if(sum>thres){
			P(pcg, iloc, 0)=sumx/sum;
			P(pcg, iloc, 1)=sumy/sum;
		}
	}
	return 0;
}
/**
   Determine the polynomial coefficients that transforms in to out.
*/
dmat* poly2fit(const dmat* in,  /**<[in] input grid*/
	const dmat* out, /**<[in] distorted grid*/
	int maxorder     /**<[in] Maximum order*/
){
	if(!(NX(in)==NX(out)&&NY(in)==2&&NY(out)==2)){
		error("poly2fit: in and out mismatch\n");
		return 0;
	}
	const int limitmax=maxorder>0?1:0;
	maxorder=abs(maxorder);
	const int nmod=(maxorder+1)*(maxorder*(limitmax?1:2)+2)/2;
	dmat* coeff=dnew(nmod, 4);
	int ic=0;
	for(int order=0; order<=maxorder; order++){
		int xmax=limitmax?order:maxorder;
		for(int xo=0; xo<=xmax; xo++){
			int yo=order+(limitmax?(-xo):0);
			P(coeff, ic, 0)=xo;//order of x
			P(coeff, ic, 1)=yo;//order of y.
			ic++;
		}
	}
	if(ic!=nmod){
		error("incorrect calculation %d vs %d\n", ic, nmod);
	}
	dmat* M=dnew(NX(in), nmod);
	for(long ip=0; ip<NX(in); ip++){
		const real x=P(in, ip, 0);
		const real y=P(in, ip, 1);
		for(ic=0; ic<nmod; ic++){
			const int xo=P(coeff, ic, 0);
			const int yo=P(coeff, ic, 1);
			P(M, ip, ic)=pow(x, xo)*pow(y, yo);
		}
	}

	dmat* MI=dpinv(M, 0);
	dmulvec(PCOL(coeff, 2), MI, PCOL(out, 0), 1);//x
	dmulvec(PCOL(coeff, 3), MI, PCOL(out, 1), 1);//y

	dfree(M);
	dfree(MI);
	return coeff;
}
/**
   Calibrate the distortion as measured using interaction matrix.
*/
dmat* loc_calib(const dsp* GA,     /**<[in] Measured interaction matrix*/
	const loc_t* aloc, /**<[in] Actuator grid*/
	const loc_t* saloc,/**<[in] Subaperture grid*/
	real dispx,      /**<[in] Beam displacement along x*/
	real dispy,      /**<[in] Beam displacement along y*/
	real scale,      /**<[in] Beam cone effect*/
	int maxorder       /**<[in] Maximum power of x/y. Negative to limit total power*/
){
	//static int count=-1; count++;
	if(NX(GA)!=saloc->nloc*2||NY(GA)!=aloc->nloc){
		error("GA, aloc, and saloc does not match\n");
	}

	const int period0=10;//period of poking on saloc
	const real stroke=1e-6; //stroke of poke

	const real sigma=period0/6.*saloc->dx;//width of poke (on saloc)
	const real denom=2*pow(sigma*scale, 2);
	const int period=round(period0*scale);
	const real periodx=aloc->dx*period;
	const real periody=aloc->dy*period;

	//Determine actuators offset from origin.
	const real xoff=fmod(aloc->locx[0], aloc->dx);
	const real yoff=fmod(aloc->locy[0], aloc->dy);

	const real diam=loc_diam(saloc); //size of aperture

	const real thres=pow((diam-MAX(periodx, periody))*0.5, 2);
	const real thres2=pow(aloc->dx*0.0001, 2);
	real saminx=saloc->nloc?saloc->locx[0]:0;
	real saminy=saloc->nloc?saloc->locy[0]:0; //origin of saloc
	for(int iloc=0; iloc<saloc->nloc; iloc++){
		if(saminx>saloc->locx[iloc]){
			saminx=saloc->locx[iloc];
		}
		if(saminy>saloc->locy[iloc]){
			saminy=saloc->locy[iloc];
		}
	}
	//First, create a poke pattern and record valid poking locations
	dmat* opd=dnew(aloc->nloc, 1);
	dmat* gloc=dnew(diam*diam*scale*scale/(periodx*periody), 2);
	int ng=0;
	//From saloc to aloc: scaloc*scale+displace
	for(int iloc=0; iloc<aloc->nloc; iloc++){
	//Location of closest poke. Make sure it is on actuator.
		real locx0=round((aloc->locx[iloc]-xoff)/periodx)*periodx+xoff;
		real locy0=round((aloc->locy[iloc]-yoff)/periody)*periody+yoff;
		//Gaussian peak
		real R1=pow(aloc->locx[iloc]-locx0, 2)+pow(aloc->locy[iloc]-locy0, 2);
		P(opd, iloc)=stroke*exp(-R1/denom);

		if(R1<thres2){//This is the actuator being poked, compute poke index on saloc+cure
			real sax=(locx0-dispx)/scale;
			real say=(locy0-dispy)/scale;
			real RR=sax*sax+say*say;
			if(RR<thres){//Only use if fully inside aperture
				P(gloc, ng, 0)=(sax-saminx)/saloc->dx;
				P(gloc, ng, 1)=(say-saminy)/saloc->dy;
				ng++;
				if(ng==NX(gloc)){
					dresize(gloc, NX(gloc)*2, NY(gloc));
				}
			}
		}
	}
	dresize(gloc, ng, 2);

	//Calculate gradients for this poke pattern.
	dmat* gg=NULL;
	dspmm(&gg, GA, opd, "nn", 1);
	dfree(opd);
	dmat* opd2=NULL;
	//Use cure to reconstruct OPD of poke pattern
	cure_loc(&opd2, gg, saloc);
	dfree(gg);
	dmat* cg=NULL;
	//Determine actual poke position using CoG.
	cog_multi(&cg, opd2, gloc, period-2);

	dfree(opd2);
	//convert to metrix coordinate. Remove invalid apertures

	real imax;//max of sub-image intensity
	dmaxmin(PCOL(cg, 2), NX(cg), &imax, 0);
	imax*=0.4;//threshold to keep
	int jloc=0;
	for(int iloc=0; iloc<NX(gloc); iloc++){
		if(P(cg, iloc, 2)>imax){//only keep valid apertures
			P(gloc, jloc, 0)=P(gloc, iloc, 0)*saloc->dx+saminx;
			P(gloc, jloc, 1)=P(gloc, iloc, 1)*saloc->dy+saminy;
			P(cg, jloc, 0)=P(cg, iloc, 0)*saloc->dx+saminx;
			P(cg, jloc, 1)=P(cg, iloc, 1)*saloc->dy+saminy;
			jloc++;
		}
	}
	dresize(gloc, jloc, 2);
	dresize(cg, jloc, 2);

	//Determine distortion parameter.
	dmat* coeff=poly2fit(cg, gloc, maxorder);

	dfree(gloc);
	dfree(cg);

	return coeff;
}
