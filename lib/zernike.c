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

#include "zernike.h"
/**
   Generating Zernike Rnm for radial order ir, azimuthal or im.  used by loc_zernike.
   */
static dmat *genRnm(dmat *locr, int ir, int im){
    if(ir<0 || im < 0|| im>ir || (ir-im)%2!=0)
	error("Invalid ir, im (%d, %d)\n", ir, im);
    
    const long nloc=locr->nx;
    dmat *Rnm=dnew(nloc,1);
    for(int s=0; s<=(ir-im)/2; s++){
	double coeff=pow(-1,s)*factorial(ir-s)/factorial(s)/factorial((ir+im)/2-s)/factorial((ir-im)/2-s);
	int power=ir-2*s;
	if(power==0){
	    for(long iloc=0; iloc<nloc; iloc++){
		Rnm->p[iloc]+=coeff;
	    }
	}else if(power==1){
	    for(long iloc=0; iloc<nloc; iloc++){
		Rnm->p[iloc]+=coeff*locr->p[iloc];
	    }
	}else{
	    for(long iloc=0; iloc<nloc; iloc++){
		Rnm->p[iloc]+=coeff*pow(locr->p[iloc],power);
	    }
	}
    }
    return Rnm;
}

/**
   Create Zernike modes on loc with radius R and radial order upto nr
   nr=0 is piston
   nr=1 is tip/tilt
   nr=2 is quadratic modes
   
*/
dmat* zernike(loc_t *loc, double R, int nr){
    if(nr<0) error("Invalid nr\n");
    int nmod=(nr+1)*(nr+2)/2;
    const long nloc=loc->nloc;
    dmat *restrict MOD=dnew(nloc,nmod);
    
    dmat *restrict locr=dnew(nloc,1);
    dmat *restrict locs=dnew(nloc,1);
    const double *restrict locx=loc->locx;
    const double *restrict locy=loc->locy;
    const double R1=1./R;
    for(long iloc=0; iloc<nloc; iloc++){
	locr->p[iloc]=sqrt(pow(locx[iloc],2)+pow(locy[iloc],2))*R1;
	locs->p[iloc]=atan2(locy[iloc], locx[iloc]);
    }
    int cmod=0;
    for(int ir=0; ir<=nr; ir++){
	for(int im=0; im<=ir; im++){
	    if((ir-im)%2!=0) continue;
	    dmat *Rnm=genRnm(locr, ir, im);
	    /*dwrite(Rnm,"Rnm_%d_%d",ir,im); */
	    if(im==0){
		double coeff=sqrt(ir+1.);
		double *restrict pmod=MOD->p+nloc*cmod;
		for(long iloc=0; iloc<nloc; iloc++){
		    pmod[iloc]=Rnm->p[iloc]*coeff;
		}
		cmod++;
	    }else{
		double coeff=sqrt(2*(ir+1.));
		double *restrict pmodc;
		double *restrict pmods;
		if((cmod+1) % 2 == 1){
		    pmods=MOD->p+nloc*cmod;
		    pmodc=MOD->p+nloc*(cmod+1);
		}else{
		    pmodc=MOD->p+nloc*cmod;
		    pmods=MOD->p+nloc*(cmod+1);
		}
		for(long iloc=0; iloc<nloc; iloc++){
		    pmods[iloc]=Rnm->p[iloc]*coeff*sin(im*locs->p[iloc]);
		    pmodc[iloc]=Rnm->p[iloc]*coeff*cos(im*locs->p[iloc]);
		}
		cmod+=2;
	    }
	    dfree(Rnm);
	}
    }
    dfree(locr);
    dfree(locs);
    return MOD;
}
/**
   return zernike index of radial mode ir, and azimuthal mode im. Maximum radial
   order is nr (inclusive)
 */
imat *zernike_index(int nr){
    if(nr<0) return 0;
    imat *out=inew(nr+1, nr+1);
    for(long ix=0; ix<(nr+1)*(nr+1); ix++){
	out->p[ix]=-1;
    }
    int count=0;
    for(int ir=0; ir<=nr; ir++){
	for(int im=0; im<=ir; im++){
	    if((ir-im)%2!=0) continue;
	    out->p[im+ir*(nr+1)]=count;
	    if(im==0){//single, pure radial mode
		count++;
	    }else{//double mode as a pair
		count+=2;
	    }
	}
    }
    return out;
}
/**
   Covariance of zernike modes in Kolmogorov Turbulence. Only modes with the
   same m have non-zero covariance.

   Based on Eq 3.14 in Adaptive Optics in Astronomy (Roddier 1999).
   Verified against values in J.Y.Wang 1978, table II(a,b).
*/
dmat *zernike_turb_cov(int nr){
    int nmod=(nr+1)*(nr+2)/2;
    dmat *res=dnew(nmod, nmod);
    PDMAT(res, pres);
    imat *zind=zernike_index(nr);
    for(int ir=0; ir<=nr; ir++){
	for(int im=0; im<=ir; im++){
	    if((ir-im)%2!=0) continue;
	    long ict0=zind->p[im+ir*(nr+1)];
	    long icts=0, ictc=0;
	    if(ict0%2==1){
		ictc=ict0;//cos term
		icts=ict0+1;
	    }else{
		icts=ict0;
		ictc=ict0+1;
	    }
	    for(int jr=0; jr<=ir; jr++){
		long jct0=zind->p[im+jr*(nr+1)];
		long jcts=0, jctc=0;
		if((jct0)%2==1){
		    jctc=jct0;
		    jcts=jct0+1;
		}else{
		    jcts=jct0;
		    jctc=jct0+1;
		}

		if(jct0==-1){//doesn't have the same azimuthal mode
		    continue;
		}
		double tmp=7.2e-3*sqrt((ir+1.)*(jr+1.))*pow(-1, (ir+jr-2*im)*0.5)*pow(M_PI, 8./3.)
		    *(tgamma(14./3.)*tgamma((ir+jr-5./3.)*0.5))
		    /(tgamma((ir-jr+17./3.)*0.5)*tgamma((jr-ir+17./3.)*0.5)*tgamma((ir+jr+23./3.)*0.5));
		if(im==0){
		    pres[jct0][ict0]=pres[ict0][jct0]=tmp;
		}else{
		    pres[ictc][jctc]=pres[jctc][ictc]=tmp;
		    pres[icts][jcts]=pres[jcts][icts]=tmp;
		}
	    }
	}
    }
    ifree(zind);
    return res;
}
