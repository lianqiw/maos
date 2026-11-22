/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "libmisc.h"
#include "cure.h"
#include "accphi.h"
/**arcsec in unit of radian*/
const real AS2RAD=4.848136811095360e-06;
/**arcsec in unit of radian*/
const real MAS2RAD=4.848136811095360e-09;
/**radian in unit of arcsec*/
const real RAD2AS=206264.8062470964;
/**radian to milli-arcsecond*/
const real RAD2MAS=206264806.2470964;

/**
   add photon and read out noise.  pcalib part of bkgrnd is calibrated
   out. pcalib2 part of bkgrnd2 is calibrated out.

   @param[in,out] A The pixel intensity array
   @param[in,out] rstat The random stream
   @param[in] bkgrnd Real background in PDEs per pixel per frame
   @param[in] bkgrndc Removed background in PDEs per pixel per frame
   @param[in] bkgrnd2 Real background in PDEs of each pixel per frame. May be NULL.
   @param[in] bkgrnd2c Removed background in PDEs of each pixel per frame. May be NULL.
   @param[in] qe Pixel dependent Quantum Efficiency. May be NULL.
   @param[in] rne Read out noise per pixel per read
   @param[in] excess Excess noise factor. Values below 1 are treated as 1.
*/
void addnoise(dmat* A, rand_t* rstat, const real bkgrnd, const real bkgrndc, const dmat* bkgrnd2, const dmat* bkgrnd2c,
	const dmat* qe, const real rne, real excess){
	long np=NX(A)*NY(A);
	assert(!bkgrnd2||NX(bkgrnd2)*NY(bkgrnd2)==np);
	assert(!bkgrnd2c||NX(bkgrnd2c)*NY(bkgrnd2c)==np);
	assert(!qe||NX(qe)*NY(qe)==np);
	if(excess<1) excess=1;
	for(int ix=0; ix<np; ix++){
		real tot=P(A,ix)+bkgrnd+(bkgrnd2?P(bkgrnd2,ix):0);
		real corr=bkgrndc+(bkgrnd2c?P(bkgrnd2c,ix):0);
		real scale=1;
		if(qe){//detected photons are reduced by qe. scale is the flat-field correction.
			tot*=P(qe,ix);
			scale=1./P(qe,ix);
		}
		P(A,ix)=(randp(rstat, tot*excess)+tot*(1.-excess)+rne*randn(rstat))*scale-corr;
	}
}
/**
   Add noise to gradients according to neal, which is LL' decomposition of the sanea

   @param[in,out] grad Gradients of 2*nsa with all x gradients first. Noise is added to it.
   @param[in] neal nsa*3: the xx, yy, and xy elements of the lower triangular cholesky factor of NEA (radian).
   @param[in,out] srand The random stream.
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
		//if(isfinite(neax[isa])&&isfinite(neay[isa])){
			ggx[isa]+=errx;
			ggy[isa]+=erry;
		//}else{//zero out gradients if below threshold.
		//	ggx[isa]=0;
		//	ggy[isa]=0;
		//}
	}
}
/**
   Determine the CoG on multiple locations near the nominal position.

   @param[in,out] cg n*3 The determined position (in pixel) and total intensity. Allocated if NULL.
   Position is NaN if the peak is outside of the window.
   @param[in] im 2-d array contains the image
   @param[in] loc n*2 array. contains the nominal position (in pixel) of each peak
   @param[in] wsize window size
   @return 0.
*/
int cog_multi(dmat** cg, const dmat* im, const dmat* loc, int wsize){
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
   compute cog NEA using Monte Carlo realizations of noise

   @param[out] nea The NEA covariance (xx, xy, yx, yy) in pixel^2.
   @param[in] ints Noise free subaperture image.
   @param[in] cogmask Mask for CoG computation. May be NULL.
   @param[in] cogthres Threshold for CoG.
   @param[in] cogoff Offset (subtracted) for CoG.
   @param[in] ntry Number of Monte Carlo realizations.
   @param[in,out] rstat The random stream. Reseeded with 1 at start.
   @param[in] bkgrnd Real background in PDEs per pixel per frame.
   @param[in] bkgrndc Removed background in PDEs per pixel per frame.
   @param[in] bkgrnd2i Real background in PDEs of each pixel per frame. May be NULL.
   @param[in] bkgrnd2ic Removed background in PDEs of each pixel per frame. May be NULL.
   @param[in] rne Read out noise per pixel per read.
*/
void cog_nea(real* nea, const dmat* ints, const dmat* cogmask, real cogthres, real cogoff, int ntry,
	rand_t* rstat, real bkgrnd, real bkgrndc, const dmat* bkgrnd2i, const dmat* bkgrnd2ic, real rne
){
	dmat* ints2=dnew(NX(ints), NY(ints));
	real gnf[2]={0,0};
	real gny[2]={0,0};
	dcog(gnf, ints, 0, 0, cogthres, cogoff, 0, cogmask);
	seed_rand(rstat, 1);/*reset the seed each time.*/
	nea[0]=0; nea[1]=0; nea[2]=0; nea[3]=0;
	for(int i=0; i<ntry; i++){
		dcp(&ints2, ints);
		addnoise(ints2, rstat, bkgrnd, bkgrndc, bkgrnd2i, bkgrnd2ic, 0, rne, 1);
		dcog(gny, ints2, 0, 0, cogthres, cogoff, 0, cogmask);
		real errx=gny[0]-gnf[0];
		real erry=gny[1]-gnf[1];
		nea[0]+=errx*errx;
		nea[1]+=errx*erry;
		nea[3]+=erry*erry;
	}
	dfree(ints2);
	real stry=1./ntry;
	nea[0]=nea[0]*stry;
	nea[3]=nea[3]*stry;
	nea[1]=nea[1]*stry;
	nea[2]=nea[1];
}

/**
   Determine the polynomial coefficients that transforms 2-d coordinate (of a grid) in to out.

   @param[in] in input grid. n*2
   @param[in] out distorted grid, n*2
   @param[in] maxorder Maximum order. Positive to limit the total power of x and y; negative to limit the power of each.
   @return Newly allocated nmod*4 array with the power of x, power of y, coefficient for x, and coefficient for y
   of each term. The caller frees it.
*/
dmat* poly2fit(const dmat* in, const dmat* out, int maxorder){
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

   @param[in] GA Measured interaction matrix
   @param[in,out] propdata See accphi.h. locin is the actuator grid and locout the subaperture grid.
   Defaults are filled by prop_prep().
   @param[in] maxorder Maximum order. See poly2fit().
   @return Newly allocated distortion coefficients as from poly2fit(), or NULL if the layer is above the guide star.
   The caller frees it.
*/
dmat* loc_calib(const dsp* GA, propdata_t *propdata, int maxorder){
	const loc_t *aloc=propdata->locin;
	const loc_t *saloc=propdata->locout;

	//static int count=-1; count++;
	if(NX(GA)!=saloc->nloc*2||NY(GA)!=aloc->nloc){
		error("GA, aloc, and saloc does not match\n");
		return NULL;
	}
	real dispx, dispy, alpha, scale;
	if(!prop_prep(&dispx, &dispy, &alpha, &scale, propdata)){
		dbg("Layer is above guide star\n");
		return NULL;
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
	dvecmaxmin(PCOL(cg, 2), NX(cg), &imax, 0);
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
/**
   Evaluate the polynomial with coefficients coeff (descending order, as from
   polyfit) at x.

   @param[in] x input vector
   @param[in] coeff the coefficient from polyfit
   @param[in] separate if set, separately for each order
   @return Newly allocated column vector, or one column per term if separate is set. The caller frees it.
*/
dmat* polyval(const dmat *x, const dmat *coeff, const int separate){
	int mul=separate?1:0;//column multiplier
	long nmod=separate?PN(coeff):1;
	dmat *out=dnew(PN(x), nmod);
	for(long iy=0; iy<PN(coeff); iy++){
		long icol=iy*mul;
		real pc=P(coeff,iy); //coefficient
		real po=PN(coeff)-1-iy;//descending order
		for(long ix=0; ix<PN(x); ix++){
			//not optimized for speed.
			P(out, ix, icol)+=pow(P(x,ix),po)*pc;
		}
	}
	return out;
}

/**
   Determine the polynomial coefficients that transforms vector in to out.
   Both column and row vectors are allowed.

   @param[in] x input vector
   @param[in] y output
   @param[in] maxorder maximum order to fit
   @return Newly allocated coefficients in descending order (for polyval). The caller frees it.
*/
dmat *polyfit(const dmat *x, const dmat *y, const int maxorder){
	if(!(PN(x)==PN(y) && (NX(x)==1 || NY(x)==1) && (NX(y)==1 || NY(y)==1))){
		error("polyfit: in and out must be vector\n");
		return NULL;
	}
	dmat *coeff=dnew(maxorder+1,1);
	dset(coeff, 1);
	dmat *mod=polyval(x, coeff, 1);
	dmat *my=0;//M^T * y
	dmat *mm=0;//M^T * M
	dmm(&my, 1, mod, y, "tn", 1);
	dmm(&mm, 1, mod, mod, "tn", 1);
	dinvspd_inplace(mm);
	dzero(coeff);
	dmm(&coeff, 1, mm, my, "nn", 1);
	dfree(mod);
	dfree(my);
	dfree(mm);
	return coeff;
}



/**
   Wrap the index for dataset with total of n frames for continuity. For example, if n is 4, the result is 
   0 1 2 3 2 1 0 1 2 3 2 1 0 ...

   @param[in] index Index to wrap. May be negative.
   @param[in] n Total number of frames.
   @return Wrapped index between 0 and n-1. 0 if n<2.
*/
long wrap_seq(long index, long n){
	if(n<2) return 0;
	long m=n*2-2;
	index=(index%m+m)%m;
	if(index>=n) index=m-index;
	return index;
}

/**
   wrap val to between low and high
   Notice the difference between remainder() and fmod().
   remainder(a,b) wraps a to [-b/2, b/2)
   fmod(a,b) wraps a to [-b,b)

   @param[in] val Value to wrap.
   @param[in] low Lower bound.
   @param[in] high Upper bound.
   @return val wrapped to [low, high] (centered at (low+high)/2) where range=|high-low|; low if low==high.
*/
real wrap2range(real val, real low, real high){
	if(low==high){
		return low;
	}
	real med=0.5*(low+high);
	real range=fabs(high-low);
	return med+remainder(val-med, range);
}
/**
   compute p-chip interpolation weighting
   See https://en.wikipedia.org/wiki/Cubic_Hermite_spline

   @param[in] i index into p[], p[i] with i from -1 to 2
   @param[in] u fractional distance from i. between [0,1)
   @return the weighting function
*/
real pchip_wt(int i, float u){
	if(i<-1||i>2||u<0||u>1){
		error("Invalid i=%d, u=%g\n", i, u);
		return 0;
	}
	real ans=0;
	switch(i){
	case -1:
		ans=0.5*(u*u*(2-u)-u); break;
	case 0:
		ans=0.5*(u*u*(3*u-5)+2); break;
	case 1:
		ans=0.5*(u*u*(4-3*u)+u); break;
	case 2:
		ans=0.5*(u*u*(u-1)); break;
	}
	return ans;
}
/**
   Compute dome seeing turbulence interpolation parameters

   @param[out] wt Returned weight
   @param[in] ips Current layer index
   @param[in] isim Current simulation frame index
   @param[in] dtrat Turbulence dtrat
   @param[in] natm Turbulence frame count
   @param[in] interp Interpolation method: 0: no interpolation. 1: linear. 2: sin^2. 3: p-chip.
   @return Turbulence frame to use
*/
int atm_interp(real *wt, int ips, int isim, int dtrat, int natm, int interp){
	if(interp==3) ips=ips-1;
	int iframe=(int)wrap_seq(isim/dtrat+ips, natm);
	if(interp<=0){
		*wt=0;
	} else{
		*wt=(real)(isim%dtrat)/dtrat;
		if(interp==2){
			*wt=pow(sin(*wt*M_PI*0.5), 2);
		} else if(interp==3){
			*wt=pchip_wt(ips, *wt);
		}
	}
	if(ips==0&&interp<3) *wt=1.-*wt;
	return iframe;
}
/**
   Compare the first element of two double pairs for qsort in ascending order.

   @param[in] a Pointer to the first pair.
   @param[in] b Pointer to the second pair.
   @return -1, 0, or 1 if the first element of a is less than, equal to, or greater than that of b.
*/
static int compare_dbl2_ascend(const void *a, const void *b){
	const real va=((const real *)a)[0];
	const real vb=((const real *)b)[0];
	return (va>vb)-(va<vb);
}
/**
   Prepare data structure to plot WFE vs field angle.
   Use together with wfe_fov_fill()

   @param[in] thetax Field angle along x of each evaluation direction.
   @param[in] thetay Field angle along y of each evaluation direction.
   @param[in] unit Scaling to convert field angle to the plotting unit.
   @param[in] npath Number of simulation directories
   @param[in] nwvl Number of wavelengths for Strehl plotting. 0 to skip.
   @return Newly allocated cell with sorted direction index, and radius vs WFE (total and high order) and Strehl
   for each path, or NULL if there is only one direction. The caller frees it.
*/
dccell *wfe_fov_prep(const dmat *thetax, const dmat *thetay, real unit, int npath, int nwvl){
	const int nevl=(int)PN(thetax);
	if(nevl<=1) return NULL;
	dccell *reseach=dccellnew(3+(nwvl?1:0), 1);
	P(reseach, 0)=dcellnew_same(1, 1, nevl, 1);//CL vs dir indexing
	P(reseach, 1)=dcellnew_same(npath, 1, nevl, 2);//CL vs dir WEF
	P(reseach, 2)=dcellnew_same(npath, 1, nevl, 2);//CL vs dir WEF
	if(nwvl) P(reseach, 3)=dcellnew_same(nwvl, npath, nevl, 2);//CL vs dir Strehl
	// sort evaluation directions
	dmat *tmp=dnew(2, nevl);
	for(int ievl=0; ievl<nevl; ievl++){
		real r=RSS(P(thetax, ievl), P(thetay, ievl))*unit;
		P(tmp, 0, ievl)=r;
		P(tmp, 1, ievl)=ievl;//mark direction
	}
	qsort(P(tmp), NY(tmp), 2*sizeof(real), compare_dbl2_ascend);
	for(int ievl=0; ievl<nevl; ievl++){
		P(P(P(reseach, 0), 0), ievl, 0)=P(tmp, 1, ievl);//index
		for(int ipath=0; ipath<npath; ipath++){
			P(P(P(reseach, 1), ipath), ievl, 0)=P(tmp, 0, ievl);	
			P(P(P(reseach, 2), ipath), ievl, 0)=P(tmp, 0, ievl);	
			for(int iwvl=0; iwvl<nwvl; iwvl++){//for Strehl plotting
				P(P(P(reseach, 3), iwvl, ipath), ievl, 0)=P(tmp, 0, ievl);
			}
		}
	}
	dfree(tmp);
	
	return reseach;
}
/**
   Compute RMS WFE for each direction.
   Use together with wfe_fov_prep()

   @param[in,out] reseach Data structure from wfe_fov_prep(). The WFE (nm) of path ipath is filled.
   @param[in] clep Closed loop wavefront error time history of each direction.
   @param[in] ipath Index of the simulation directory.
   @param[in] isim Number of time steps to use. 0 or negative to use all.
*/
void wfe_fov_fill(dccell *reseach, const dcell *clep, int ipath, int isim){
	if(!reseach || !clep) return;
	if(isim<=0) isim=(int)NY(P(clep,0));
	dmat *res=P(P(reseach, 1), ipath);
	dmat *res_ho=P(P(reseach, 2), ipath);
	const int nevl=(int)NX(res);
	int istart=MAX((int)(0.2*isim), 20);
	if(istart>=isim) istart=0;
	for(int ievl=0; ievl<nevl; ievl++){
		int jevl=(int)P(P(P(reseach, 0), 0), ievl);
		P(res, ievl, 1)=0;
		P(res_ho, ievl, 1)=0;
		for(int i=istart; i<isim; i++){
			P(res, ievl, 1)+=P(P(clep, jevl),0,i);
			P(res_ho, ievl, 1)+=P(P(clep, jevl),2,i);
		}
		P(res, ievl, 1)=sqrt(P(res, ievl, 1)/(isim-istart))*1e9;
		P(res_ho, ievl, 1)=sqrt(P(res_ho, ievl, 1)/(isim-istart))*1e9;
	}
	//check directions of the same radius and average them.
	for(int ievl=0; ievl<nevl-1; ievl++){
		int jevl=ievl+1;
		for(; jevl<nevl; jevl++){
			if(fabs(P(res, ievl, 0)-P(res, jevl, 0))>1){
				break;
			}
		}
		if(jevl>ievl+1){//directions of the same radius
			real sum2=0;
			real sum2_ho=0;
			for(int i=ievl; i<jevl; i++){
				sum2+=pow(P(res, i, 1),2);
				sum2_ho+=pow(P(res_ho, i, 1),2);
			}
			sum2=sqrt(sum2/(jevl-ievl));
			sum2_ho=sqrt(sum2_ho/(jevl-ievl));
			for(int i=ievl; i<jevl; i++){
				P(res, i, 1)=sum2;
				P(res_ho, i, 1)=sum2_ho;
			}
		}
	}
	
}
/**
   Create WFS amplitude map from coordinate, masked with annular defined by (D,Din).

   @param[in] loc Coordinate of the amplitude map.
   @param[in] aper Aperture map. If not NULL, it is ray traced to loc instead of using (D,Din).
   @param[in] misregx Misregistration along x (m).
   @param[in] misregy Misregistration along y (m).
   @param[in] D Outer diameter (m).
   @param[in] Din Inner diameter (m).
   @return Newly allocated amplitude; the caller frees it.
*/
dmat* mkamp(const loc_t* loc, const map_t* aper, real misregx, real misregy, real D, real Din){
	dmat* amp=dnew(loc->nloc, 1);
	if(aper){
		prop(ADDR((propdata_t){.mapin=aper, .locout=loc, .phiout=P(amp), .misregx=misregx, .misregy=misregy}));
	}else{
		loc_circle_add(amp, loc, -misregx, -misregy, D*0.5, Din*0.5, 1);
	}
	/*if(misregx || misregy){//apply centered circular mask
		loc_circle_mul(amp, loc, 0, 0, D*0.5, 0, 1);
	}*/
	return amp;
}
/**
   Fit 3 points around the peak of a 1-d profile.

   @param[in] corr The 1-d profile.
   @return Peak location (in pixels) relative to the center.
*/
real parabolic_peak_1d(dmat* corr){
	if(corr->nx<3) return 0;
	real valmax=-INFINITY;
	int jx=1;
	for(int ix=1; ix<corr->nx-1; ix++){
		if(P(corr, ix)>valmax){
			jx=ix;
			valmax=P(corr, ix);
		}
	}
	//Parabolic fit.
	real px[2];
	px[0]=(P(corr, jx+1)+P(corr, jx-1))*0.5-P(corr, jx);
	px[1]=(P(corr, jx+1)-P(corr, jx-1))*0.5;
	return px[0]==0?0:(-px[1]/(2*px[0])+jx-(corr->nx-1)*0.5);
}
/**
   First sum along 1 dimension, then fit 3 points around the peak. More robust than the old method.

   @param[out] grad Peak location along x and y (in pixels) relative to the center.
   @param[in] corr The 2-d correlation.
   @param[in] nbox Size of the central box to use. Use the full array if not positive or too large.
*/
void parabolic_peak_sum(real* grad, dmat* corr, int nbox){
	const long nx=NX(corr);
	const long ny=NY(corr);
	if(nbox<=0||nbox>nx) nbox=(int)nx;
	if(nbox<=0||nbox>ny) nbox=(int)ny;
	dmat* corrx=dnew(nbox, 1);
	dmat* corry=dnew(nbox, 1);
	const long offx=(nx-nbox)/2;
	const long offy=(ny-nbox)/2;
	for(int iy=0; iy<nbox; iy++){
		for(int ix=0; ix<nbox; ix++){
			P(corrx, ix)+=P(corr, ix+offx, iy+offy);
			P(corry, iy)+=P(corr, ix+offx, iy+offy);
		}
	}
	grad[0]=parabolic_peak_1d(corrx);
	grad[1]=parabolic_peak_1d(corry);
	if(0){
		info2("grad=%g, %g\n", grad[0], grad[1]);
		writebin(corr, "corr");
		exit(0);
	}
}
