/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "genotf.h"
#include "turbulence.h"
/**
private data struct to mark valid pairs of points.  */
typedef struct T_VALID{
	long n;
	long(*loc)[2];
}T_VALID;
/**
 Wrap the data to genotf to have multi-thread capability.*/
typedef struct GENOTF_T{
	ccell* otf;
	loc_t* loc;     /**<the common aperture grid*/
	const dmat* amp;    /**<The amplitude map of all the (sub)apertures*/
	const dmat* opdbias;/**<The static OPD bias. */
	const dmat* area;   /**<area of a (sub)aperture*/
	real thres;/**<threshold to consider an (sub)aperture as full*/
	real wvl;  /**<The wavelength. only needef if opdbias is not null*/
	long npsfx;   /**<Size of PSF/OTF*/
	long npsfy;   /**<Size of PSF/OTF*/
	long nsa;    /**<Number of (sub)apertures*/
	long pttr;   /**<Remove piston/tip/tilt*/
	const dmat* B;
	const T_VALID* pval;
	long isafull;
	const cmat* otffull;
}GENOTF_T;
/**
   Remove tip/tilt from the covariance matrix.
*/
static dmat* pttr_B(const dmat* B,   /**<The B matrix. */
	loc_t* loc,       /**<The aperture grid*/
	const real* amp /**<The amplitude map*/
){
	if(!amp) error("amplitude map has to be not empty to remove pistion/tip/tilt\n");
	real* locx=loc->locx;
	real* locy=loc->locy;
	int nloc=loc->nloc;

	dmat* B2=dnew(nloc, nloc);
	dmat* BP=B2;

	real* mod[3];
	dmat* mcc=dnew(3, 3);/*modal cross coupling matrix. */

	mod[0]=NULL;
	mod[1]=locx;
	mod[2]=locy;
	for(int im=0; im<3;im++){
		for(int jm=im;jm<3;jm++){
			P(mcc, im, jm)=P(mcc, jm, im)=dvecdot(mod[im], mod[jm], amp, nloc);
		}
	}
	dinvspd_inplace(mcc);
	dmat* M=dnew(nloc, 3);/*The tip/tilt modal matrix */
	dmat* MW=dnew(nloc, 3);/*M*W */
	dmat* MCC=dnew(3, nloc);/*M*inv(M'*W*M) */
	dmat* Mtmp=dnew(3, nloc);/*B'*MW; */

	for(long iloc=0; iloc<nloc; iloc++){
		P(M,iloc,0)=1;
	}
	memcpy(PCOL(M, 1), locx, nloc*sizeof(real));
	memcpy(PCOL(M, 2), locy, nloc*sizeof(real));
	for(long iloc=0; iloc<nloc; iloc++){
		P(MW,iloc,0)=amp[iloc];
		P(MW,iloc,1)=amp[iloc]*locx[iloc];
		P(MW,iloc,2)=amp[iloc]*locy[iloc];
	}
	/* MCC = - cci' *M' */
	dmm(&MCC, 0, mcc, M, "tt", -1);
	dmat* pMCC=MCC;
	/* Mtmp =  MW' * B  */
	dmm(&Mtmp, 0, MW, B, "tn", 1);
	/*Remove tip/tilt from left side*/
	dmat* pMtmp=Mtmp;
	for(long iloc=0; iloc<nloc; iloc++){
		real tmp1=P(pMtmp, 0, iloc);
		real tmp2=P(pMtmp, 1, iloc);
		real tmp3=P(pMtmp, 2, iloc);
		for(long jloc=0; jloc<nloc; jloc++){
			P(BP, jloc, iloc)=P(B, jloc, iloc)+
				(P(pMCC, 0, jloc)*tmp1
					+P(pMCC, 1, jloc)*tmp2
					+P(pMCC, 2, jloc)*tmp3);
		}
	}
	/* Mtmp = MW' * BP' */
	dmm(&Mtmp, 0, MW, B2, "tt", 1);
	/*Remove tip/tilt from right side*/
	for(long iloc=0; iloc<nloc; iloc++){
		real tmp1=P(pMCC, 0, iloc);
		real tmp2=P(pMCC, 1, iloc);
		real tmp3=P(pMCC, 2, iloc);
		for(long jloc=0; jloc<nloc; jloc++){
			P(BP, jloc, iloc)+=
				tmp1*P(pMtmp, 0, jloc)
				+tmp2*P(pMtmp, 1, jloc)
				+tmp3*P(pMtmp, 2, jloc);
		}
	}
	dfree(mcc);
	dfree(M);
	dfree(MW);
	dfree(MCC);
	dfree(Mtmp);
	return B2;
}
/**
   Generate OTF from the B or tip/tilted removed B matrix. Notice that tip/tilt
   in opdbias is NOT removed.
*/
static void genotf_do(cmat** otf, long pttr, long npsfx, long npsfy,
	loc_t* loc, const real* amp, const real* opdbias, real wvl,
	const dmat* B, const T_VALID* pval){
		{
			real ampsum=dvecsum(amp, loc->nloc);
			real ampmax;
			dmaxmin(amp, loc->nloc, &ampmax, 0);
			if(ampsum<=ampmax*0.1*sqrt((real)loc->nloc)){
				warning("genotf_do: amplitude may be too sparse\n");
			}
		}
		long nloc=loc->nloc;
		dmat* BP;
		if(pttr){/*remove p/t/t from the B matrix */
			BP=pttr_B(B, loc, amp);
		} else{
			BP=ddup(B);/*duplicate since we need to modify it. */
		}

		if(!*otf){
			*otf=cnew(npsfx, npsfy);
		}
		/*Do the exponential.*/
		real k2=pow(2*M_PI/wvl, 2);
		dmat* BPD=dnew(nloc, 1);
		for(long iloc=0; iloc<nloc; iloc++){
			for(long jloc=0; jloc<nloc; jloc++){
				P(BP, jloc, iloc)=exp(k2*P(BP, jloc, iloc));
			}
			P(BPD, iloc)=pow(P(BP, iloc, iloc), -0.5);
		}
		real otfnorm=0;
		if(amp){
			for(long iloc=0; iloc<nloc; iloc++){
				otfnorm+=amp[iloc]*amp[iloc];
			}
		} else{
			otfnorm=nloc;
		}
		otfnorm=1./otfnorm;

		struct T_VALID(*qval)[npsfx]=(struct T_VALID(*)[npsfx])pval;

		comp wvk=COMPLEX(0, 2.*M_PI/wvl);
		for(long jm=0; jm<npsfy; jm++){
			for(long im=0; im<npsfx; im++){
				long(*jloc)[2]=qval[jm][im].loc;
				comp tmpc=0.;
				real tmpr=0.;
				for(long iloc=0; iloc<qval[jm][im].n; iloc++){
					long iloc1=jloc[iloc][0];/*iloc1 is continuous. */
					long iloc2=jloc[iloc][1];/*iloc2 is not continuous. */
					real tmp12=P(BPD, iloc1)*P(BPD, iloc2)*P(BP, iloc2, iloc1);
					if(amp){
						tmp12*=amp[iloc1]*amp[iloc2];
					}
					if(opdbias){
						comp tmp3=cexp(wvk*(opdbias[iloc1]-opdbias[iloc2]));
						tmpc+=tmp12*tmp3;
					} else{
						tmpr+=tmp12;
					}
				}
				P(*otf, im, jm)=(tmpc+tmpr)*otfnorm;
			}
		}
		dfree(BPD);
		dfree(BP);
}
/**
   A wrapper to execute pttr parallel in pthreads
 */
static void genotf_wrap(thread_t* info){
	GENOTF_T* data=(GENOTF_T*)info->data;
	const int nsa=data->nsa;
	ccell* otf=data->otf;
	loc_t* loc=data->loc;
	const long nxsa=loc->nloc;
	const real wvl=data->wvl;
	const dmat* area=data->area;
	const real thres=data->thres;
	const cmat* otffull=data->otffull;
	const dmat* amp=data->amp;
	const dmat* opdbias=data->opdbias;
	const long pttr=data->pttr;
	const dmat* B=data->B;
	const T_VALID* pval=data->pval;
	if(!(check(!area||NX(area)*NY(area)==nsa) 
		|| check(!amp||NX(amp)*NY(amp)==nxsa*nsa) 
		||	check(!opdbias||NX(opdbias)*NY(opdbias)==nxsa*nsa))){
		error("Invalid input\n");
	}
	for(int isa=info->start; isa<info->end; isa++){
		if(!detached&&nsa>10&&info->ithread==0){
			info_console("%6ld of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", isa*(nsa/info->end), nsa);
		}
		const real* opdbiasi=NULL;
		if(opdbias){
			opdbiasi=P(opdbias)+isa*nxsa;
		} else{
			opdbiasi=NULL;
		}
		if(otffull&&(!area||P(area,isa)>thres)){
			ccp(&P(otf,isa), otffull);/*just copy the full array */
		} else if(!area||P(area,isa)>0.01){
			genotf_do(&P(otf,isa), pttr, data->npsfx, data->npsfy, loc, amp?P(amp)+isa*nxsa:NULL, opdbiasi, wvl, B, pval);
		}
	}
	//if(!detached && nsa>10) info2("Thread %ld done\n", info->ithread);
}
/**
   Generate pairs of overlapping points for structure function.

   2010-11-08: removed amp. It caused wrong otf because it uses the amp of the
   first subaperture to build pval, but this one is not fully illuminated.
 */
static T_VALID* gen_pval(long npsfx, long npsfy, loc_t* loc){
	const long nloc=loc->nloc;
	const real* locx=loc->locx;
	const real* locy=loc->locy;
	const long pvaltot=npsfx*npsfy*nloc*2;
	typedef long long2[2];
	long2* pval0=mymalloc(pvaltot, long2);
	if(!pval0){
		error("malloc for %ld failed\n", pvaltot);
	}
	T_VALID* pval=mymalloc(npsfx*npsfy, T_VALID);
	T_VALID(*restrict qval)[npsfx]=(T_VALID(*)[npsfx])(pval);
	long count=0, count2;
	loc_create_map(loc);
	map_t* map=loc->map;
	long npsfx2=npsfx/2;
	long npsfy2=npsfy/2;
	real dx1=1./loc->dx;
	real dy1=1./loc->dy;
	lmat* mloc=lnew(2, loc->nloc);
	for(long iloc=0; iloc<loc->nloc; iloc++){
		P(mloc, 0, iloc)=(long)round((locx[iloc]-map->ox)*dx1);
		P(mloc, 1, iloc)=(long)round((locy[iloc]-map->oy)*dy1);
	}
	for(long jm=0; jm<npsfy; jm++){
		long jm2=(jm-npsfy2);/*peak in the center */
		for(long im=0; im<npsfx; im++){
			long im2=(im-npsfx2);
			count2=count;
			for(long iloc=0; iloc<loc->nloc; iloc++){
				long iloc2=(long)loc_map_get(map, P(mloc, 0, iloc)+im2, P(mloc, 1, iloc)+jm2);
				if(iloc2>0){
					pval0[count][0]=iloc;
					pval0[count][1]=iloc2-1;
					count++;
				}
			}
			qval[jm][im].loc=pval0+count2;
			qval[jm][im].n=count-count2;
		}
	}
	if(count>pvaltot){
		error("count=%ld > pvaltot=%ld\n", count, pvaltot);
	}
	/*loc_free_map(loc);*//*do not free map. dangerous in multi-threaded envorionment. where other threads may be visiting loc->map.*/
	/*pval0=myrealloc(pval0,count*2,int); //do not realloc. will change position. */
	return pval;
}
/**
   Compute the separation vector of each point in loc.
*/
static dmat* loc_sep(const loc_t* loc){
	long nloc=loc->nloc;
	real* locx=loc->locx;
	real* locy=loc->locy;
	dmat* sep=dnew(nloc, nloc);
	for(long i=0; i<nloc; i++){
		for(long j=i; j<nloc; j++){
			P(sep, i, j)=P(sep, j, i)=sqrt(pow(locx[i]-locx[j], 2)+pow(locy[i]-locy[j], 2));
		}
	}
	return sep;
}
/**
   Generate OTFs for an aperture or multiple subapertures. ALl these apertures
   must share the same geometry, but may come with different amplitude map and/or
   OPD biasas. if pttr is 1, the OTF will have tip/tilt removed. make r0 to
   infinity to build diffraction limited OTF. make r0 to infinity and opdbias to
   none null to build OTF for a static map.

   2020-01-21: Compute OTF using nyquist sampling and then upsample with FFT.
*/

void genotf(ccell** potf,    /**<The otf array for output*/
	loc_t* loc,    /**<the aperture grid (same for all apertures)*/
	const dmat* amp,    /**<The amplitude map of all the (sub)apertures*/
	const dmat* opdbias,/**<The static OPD bias (complex part of amp). */
	const dmat* area,   /**<normalized area of the (sub)apertures*/
	real thres,  /**<The threshold to consider a (sub)aperture as full*/
	real wvl,    /**<The wavelength. only needef if opdbias is not null*/
	const dmat* cov,/**<The covariance. If not supplied use r0 for kolmogorov spectrum.*/
	real r0,     /**<Fried parameter*/
	real l0,     /**<Outer scale*/
	long npsfx,   /**<Size of PSF*/
	long npsfy,   /**<Size of PSF*/
	long nsa,      /**<Number of (sub)apertures*/
	long pttr      /**<Remove piston/tip/tilt*/
){
	if(amp&&loc->nloc*nsa!=NX(amp)*NY(amp)){
		error("loc and amp mismatch. loc->nloc=%ld, amp is %ldx%ld, nsa=%ld\n", loc->nloc, NX(amp), NY(amp), nsa);
	} else if(cov&&(NX(amp)!=NX(cov)||NX(cov)!=NY(cov))){
		error("loc and cov mismatch\n");
	} else if(nsa<1||npsfx<1||npsfy<1){
		error("nsa, notfx, notfy has to be at least 1\n");
	}
	/*creating pairs of points that both exist with given separation*/
	T_VALID* pval=gen_pval(npsfx, npsfy, loc);/*returns T_VALID array. */
	/* Generate the B matrix. */
	dmat* B;
	if(cov){
		B=(dmat*)cov;
	} else{
		dmat* sep=loc_sep(loc);
		B=turbcov(sep, 0, r0, l0);
		dfree(sep);
	}
	cmat* otffull=NULL;
	const long nloc=loc->nloc;
	long isafull=-1;
	if(opdbias&&dsumsq(opdbias)==0){
		opdbias=0;
	}
	if(!opdbias&&nsa>1){
		real maxarea=0;
		for(long isa=0; isa<nsa; isa++){
			if(P(area,isa)>maxarea){
				maxarea=P(area,isa);
				isafull=isa;
			}
		}
		if(isafull>0){
			genotf_do(&otffull, pttr, npsfx, npsfy, loc, amp?P(amp)+isafull*nloc:NULL, NULL, wvl, B, pval);
		}
	}
	if(!*potf){
		*potf=ccellnew_same(nsa,1,npsfx,npsfy);
	}
	GENOTF_T data={*potf, loc, amp, opdbias, area, thres, wvl, npsfx, npsfy, nsa, pttr, B, pval, isafull, otffull};

	thread_t info[NCPU];
	thread_prep(info, 0, nsa, NCPU, genotf_wrap, &data);
	CALL_THREAD(info, 1);
	cfree(otffull);
	if(!cov) dfree(B);
	free(pval[0].loc);
	free(pval);
}

/**
   Average spatially the 4-d covariance function to create a 2-d covariance
   function. For OPD f defined on points x (2-d coordinate), the 4-d covariance
   is simply <f'f> where f is vector form of the OPD and the average is over
   time. The 2-d covariance is additionally averaged over all the points so that
   B(r)=<f(x)'f(x+r)>_x,t To compute B, we first figure out the number of
   overlapping pairs of points for each r and then compute the averaging. When
   the amplitude is less than the threshold, the point does not count.*/

dmat* mk2dcov(loc_t* loc, const dmat* amp, real ampthres, const dmat* cov, int norm){
	if(loc->nloc!=NX(cov)||loc->nloc!=NY(cov)){
		error("loc and cov does not match. loc->nloc=%ld, cov is %ldx%ld\n", loc->nloc, NX(cov), NY(cov));
	}
	real xmin, xmax, ymin, ymax;
	long nloc=loc->nloc;
	real* locx=loc->locx;
	real* locy=loc->locy;
	dmaxmin(locx, nloc, &xmax, &xmin);
	dmaxmin(locy, nloc, &ymax, &ymin);
	real dx1=1./loc->dx;
	real dy1=1./loc->dy;
	long ncovx=(long)round((xmax-xmin)*dx1)*2;
	long ncovy=(long)round((ymax-ymin)*dy1)*2;
	dmat* cov2d=dnew(ncovx, ncovy);
	/*the following is adapted from gen_pval*/
	loc_create_map(loc);
	map_t* map=loc->map;
	long ncovx2=ncovx/2;
	long ncovy2=ncovy/2;
	long* map_x=mymalloc(loc->nloc, long);
	long* map_y=mymalloc(loc->nloc, long);
	for(long iloc=0; iloc<loc->nloc; iloc++){
		map_x[iloc]=(long)round((locx[iloc]-map->ox)*dx1);
		map_y[iloc]=(long)round((locy[iloc]-map->oy)*dy1);
	}
	for(long jm=0; jm<ncovy; jm++){
		long jm2=(jm-ncovy2);//peak in the center 
		/*long jm2=jm<ncovy2?jm:jm-ncovy;//peak in the corner */
		for(long im=0; im<ncovx; im++){
			long im2=(im-ncovx2);//peak in the center 
			/*long im2=im<ncovx2?im:im-ncovx; //peak in the corner */
			long count=0;
			real acc=0;
			for(long iloc=0; iloc<loc->nloc; iloc++){
				if(amp&&P(amp,iloc)<ampthres) continue;
				long ix=map_x[iloc]+im2;
				long iy=map_y[iloc]+jm2;
				long iloc2=(long)loc_map_get(map, ix, iy);
				if(iloc2>0&&(!amp||P(amp,iloc2)>=ampthres)){
					acc+=P(cov, iloc2-1, iloc);
					count++;
				}
			}
			if(count>0){
				if(norm){/*compute the covariance*/
					P(cov2d, im, jm)=acc/count;
				} else{/*compute approximate PSD.*/
					P(cov2d, im, jm)=acc;
				}
			}
		}
	}
	free(map_x);
	free(map_y);
	return cov2d;
}

