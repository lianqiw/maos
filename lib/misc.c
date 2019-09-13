/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file misc.c
   
   Miscellaneous routines.
*/
#include "misc.h"
#include "cure.h"
/**
   add photon and read out noise.  pcalib part of bkgrnd is calibrated
   out. pcalib2 part of bkgrnd2 is calibrated out.  */
void addnoise(dmat *A,              /**<[in/out]The pixel intensity array*/
	      rand_t* rstat,        /**<[in]The random stream*/
	      const double bkgrnd,  /**<[in]Real background in PDEs per pixel per frame*/
	      const double bkgrndc, /**<[in]Removed background in PDEs per pixel per frame*/
	      const dmat *bkgrnd2,  /**<[in]Real background in PDEs of each pixel per frame.*/
	      const dmat *bkgrnd2c, /**<[in]Removed background in PDEs of each pixel per frame.*/
	      const dmat *qe,       /**<[in]Pixel dependent Quantum Efficiency*/
	      const double rne,     /**<[in]Read out noise per pixel per read*/
	      double excess   /**<[in]Excess noise factor*/
    ){
    long np=A->nx*A->ny;
    assert(!bkgrnd2 || bkgrnd2->nx*bkgrnd2->ny==np);
    assert(!bkgrnd2c || bkgrnd2c->nx*bkgrnd2c->ny==np);
    if(excess<1) excess=1;
    for(int ix=0; ix<np; ix++){
	double tot=A->p[ix]+bkgrnd+(bkgrnd2?bkgrnd2->p[ix]:0);
	double corr=bkgrndc+(bkgrnd2c?bkgrnd2c->p[ix]:0);
	double scale=1;
	if(qe){//the second qe factor is flat-field correction.
	    tot*=qe->p[ix];
	    scale=1./qe->p[ix];
	}
	A->p[ix]=(randp(rstat, tot*excess)+tot*(1.-excess)+rne*randn(rstat))/scale-corr;
    }
}
/**
   Determine the CoG on multiple locations near the nominal position.
*/
int cog_multi(
    dmat **cg,      /**<[out]n*3 The determined position and total intensity.*/
    const dmat *im, /**<[in] 2-d array contains the image*/
    const dmat *loc,/**<[in] n*2 array. contains the nominal position of each peak*/
    int wsize       /**<[in] window size*/
    ){
    if(!*cg){
	*cg=dnew(loc->nx, 3);
    }
    dmat *pcg=*cg;
    int wsize2=(wsize-1)/2;
    for(int iloc=0; iloc<loc->nx; iloc++){
	int cx=P(loc, iloc, 0);
	int cy=P(loc, iloc, 1);
	int cx2=cx, cy2=cy;
	double  cmax=-INFINITY, cmin=INFINITY;
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
	if(abs(cx2-cx)>wsize2 || abs(cy2-cy)>wsize2){
	    P(pcg,iloc,0)=NAN;
	    P(pcg,iloc,1)=NAN;
	    continue;
	}
	double thres=cmax*0.2+cmin*0.8;
	double sum=0, sumx=0, sumy=0;
	for(int jy=cy2-wsize2; jy<=cy2+wsize2; jy++){
	    for(int jx=cx2-wsize2; jx<=cx2+wsize2; jx++){
		double val=P(im, jx, jy);
		if(val>thres){
		    sum+=val;
		    sumx+=val*jx;
		    sumy+=val*jy;
		}
	    }
	}
	P(pcg,iloc,2)=sum;
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
dmat *poly2fit(const dmat *in,  /**<[in] input grid*/
	       const dmat *out, /**<[in] distorted grid*/
	       int maxorder     /**<[in] Maximum order*/
    ){
    if(!(in->nx==out->nx && in->ny==2 && out->ny==2)){
	error("polyfit: in and out mismatch\n");
	return 0;
    }
    const int limitmax=maxorder>0?1:0;
    maxorder=abs(maxorder);
    const int nmod=(maxorder+1)*(maxorder*(limitmax?1:2)+2)/2;
    dmat *coeff=dnew(nmod,4);
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
    dmat *M=dnew(in->nx, nmod);
    for(long ip=0; ip<in->nx; ip++){
	const double x=P(in, ip, 0);
	const double y=P(in, ip, 1);
	for(ic=0; ic<nmod; ic++){
	    const int xo=P(coeff, ic, 0);
	    const int yo=P(coeff, ic, 1);
	    P(M, ip, ic)=pow(x, xo)*pow(y, yo);
	}
    }

    dmat *MI=dpinv(M, 0); 
    dmulvec(PCOL(coeff, 2), MI, PCOL(out, 0), 1);//x
    dmulvec(PCOL(coeff, 3), MI, PCOL(out, 1), 1);//y

    dfree(M);
    dfree(MI);
    return coeff;
}
/**
   Calibrate the distortion as measured using interaction matrix.
*/
dmat *loc_calib(const dsp *GA,     /**<[in] Measured interaction matrix*/
		const loc_t *aloc, /**<[in] Actuator grid*/
		const loc_t *saloc,/**<[in] Subaperture grid*/
		double dispx,      /**<[in] Beam displacement along x*/
		double dispy,      /**<[in] Beam displacement along y*/
		double scale,      /**<[in] Beam cone effect*/
		int maxorder       /**<[in] Maximum power of x/y. Negative to limit total power*/
    ){
    static int count=-1; count++;
    if(GA->nx!=saloc->nloc*2 || GA->ny!=aloc->nloc){
	error("GA, aloc, and saloc does not match\n");
    }

    const int period0=10;//period of poking on saloc
    const double stroke=1e-6; //stroke of poke

    const double sigma=period0/6.*saloc->dx;//width of poke (on saloc)
    const double denom=2*pow(sigma*scale,2);
    const int period=round(period0*scale);
    const double periodx=aloc->dx*period;
    const double periody=aloc->dy*period;

    //Determine actuators offset from origin.
    const double xoff=fmod(aloc->locx[0], aloc->dx);
    const double yoff=fmod(aloc->locy[0], aloc->dy);

    const double diam=loc_diam(saloc); //size of aperture

    const double thres=pow((diam-MAX(periodx,periody))*0.5,2);
    const double thres2=pow(aloc->dx*0.0001,2);
    double saminx=INT_MAX, saminy=INT_MAX; //origin of saloc
    for(int iloc=0; iloc<saloc->nloc; iloc++){
	if(saminx>saloc->locx[iloc]){
	    saminx=saloc->locx[iloc];
	}
	if(saminy>saloc->locy[iloc]){
	    saminy=saloc->locy[iloc];
	}
    }
    //First, create a poke pattern and record valid poking locations
    dmat *opd=dnew(aloc->nloc,1);
    dmat *gloc=dnew(diam*diam*scale*scale/(periodx*periody), 2);
    int ng=0;
    //From saloc to aloc: scaloc*scale+displace
    for(int iloc=0; iloc<aloc->nloc; iloc++){
	//Location of closest poke. Make sure it is on actuator.
	double locx0=round((aloc->locx[iloc]-xoff)/periodx)*periodx+xoff;
	double locy0=round((aloc->locy[iloc]-yoff)/periody)*periody+yoff;
	//Gaussian peak
	double R1=pow(aloc->locx[iloc]-locx0,2)+pow(aloc->locy[iloc]-locy0,2);
	P(opd, iloc)=stroke*exp(-R1/denom);

	if(R1<thres2){//This is the actuator being poked, compute poke index on saloc+cure
	    double sax=(locx0-dispx)/scale;
	    double say=(locy0-dispy)/scale;
	    double RR=sax*sax+say*say;
	    if(RR<thres){//Only use if fully inside aperture
		P(gloc, ng, 0)=(sax-saminx)/saloc->dx;
		P(gloc, ng, 1)=(say-saminy)/saloc->dy;
		ng++;
		if(ng==gloc->nx){
		    dresize(gloc, gloc->nx*2, gloc->ny);
		}
	    }
	}
    }
    dresize(gloc, ng, 2);

    //Calculate gradients for this poke pattern.
    dmat *gg=0;
    dspmm(&gg, GA, opd, "nn", 1);
    dfree(opd);
    dmat *opd2=0;
    //Use cure to reconstruct OPD of poke pattern
    cure_loc(&opd2, gg, saloc);
    dfree(gg);
    dmat *cg=0;
    //Determine actual poke position using CoG.
    cog_multi(&cg, opd2, gloc, period-2);

    dfree(opd2);
    //convert to metrix coordinate. Remove invalid apertures

    double imax;//max of sub-image intensity
    dmaxmin(PCOL(cg,2), cg->nx, &imax, 0);
    imax*=0.4;//threshold to keep
    int jloc=0;
    for(int iloc=0; iloc<gloc->nx; iloc++){
	if(P(cg, iloc, 2)>imax){//only keep valid apertures
	    P(gloc, jloc, 0)= P(gloc, iloc, 0)*saloc->dx+saminx;
	    P(gloc, jloc, 1)= P(gloc, iloc, 1)*saloc->dy+saminy;
	    P(cg, jloc, 0)= P(cg, iloc, 0)*saloc->dx+saminx;
	    P(cg, jloc, 1)= P(cg, iloc, 1)*saloc->dy+saminy;
	    jloc++;
	}
    }
    dresize(gloc, jloc, 2);
    dresize(cg, jloc, 2);

    //Determine distortion parameter.
    dmat*coeff=poly2fit(cg, gloc, maxorder);

    dfree(gloc);
    dfree(cg);
   
    return coeff;
}
