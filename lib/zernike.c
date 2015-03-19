/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "accphi.h"
/**
   Generating Zernike Rnm for radial order ir, azimuthal or im.  used by loc_zernike.
   */
static dmat *genRnm(const dmat *locr, int ir, int im){
    if(ir<0 || im < 0|| im>ir || (ir-im)%2!=0)
	error("Invalid ir, im (%d, %d)\n", ir, im);
    
    const long nloc=locr->nx;
    dmat *Rnm=dnew(nloc,1);
    for(int s=0; s<=(ir-im)/2; s++){
	double coeff=pow(-1,s)*factorial((ir+im)/2-s+1, ir-s)/factorial(1,s)/factorial(1,(ir-im)/2-s);
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
   Create Zernike modes on loc with diameter D and radial order upto nr
   nr=0 is piston
   nr=1 is tip/tilt
   nr=2 is quadratic modes
   if nopiston is set, skip piston mode.
   if flag is >0, only radial mode is used.
   
   if flag is is negative, only generate mode -flag. rmin and rmax is irrelevant
*/
dmat* zernike(const loc_t *loc, double D, int rmin, int rmax, int flag){
    if(flag<0){
	rmin=ceil((sqrt(8*(-flag)+1)-3)/2);
	rmax=rmin;
    }
    int nr3=(int)floor((sqrt(loc->nloc*8+1)-3)*0.5);
    if(rmax>nr3){
	warning("Reduce rmax=%d to %d\n", rmax, nr3);	
	rmax=nr3;
    }
    if(rmin>rmax) error("Invalid rmin=%d, rmax=%d\n", rmin, rmax);
    const long nloc=loc->nloc;
    int nmod=0;
    if(flag>0){//radial only
	nmod=rmax-rmin;
    }else{
	nmod=(rmax+1)*(rmax+2)/2-(rmin)*(rmin+1)/2;
    }
    double D2=loc_diam(loc);
    if(D<=0){
	D=D2;
    }else if(fabs(D-D2)>D*0.5){
	writebin(loc, "loc_wrongD");
	warning("specified diameter is incorrect. D=%g, loc D=%g\n", D, D2);
    }
 
    dmat *restrict opd=dnew(nloc,nmod);
    dmat *restrict locr=dnew(nloc,1);
    dmat *restrict locs=dnew(nloc,1);
    const double *restrict locx=loc->locx;
    const double *restrict locy=loc->locy;
    const double R1=2./D;
    for(long iloc=0; iloc<nloc; iloc++){
	locr->p[iloc]=sqrt(pow(locx[iloc],2)+pow(locy[iloc],2))*R1;
	locs->p[iloc]=atan2(locy[iloc], locx[iloc]);
    }
    int cmod=0;
    for(int ir=rmin; ir<=rmax; ir++){
	for(int im=0; im<=ir; im++){
	    if((ir-im)%2!=0) continue;
	    dmat *Rnm=genRnm(locr, ir, im);
	    if(im==0){/*Radial*/
		double coeff=sqrt(ir+1.);
		double *restrict pmod=opd->p+nloc*cmod;
		for(long iloc=0; iloc<nloc; iloc++){
		    pmod[iloc]=Rnm->p[iloc]*coeff;
		}
		cmod++;
	    }else if(flag<=0){
		double coeff=sqrt(2*(ir+1.));
		double *restrict pmodc;
		double *restrict pmods;
		if((cmod+1) % 2 == 1){
		    pmods=opd->p+nloc*cmod;
		    pmodc=opd->p+nloc*(cmod+1);
		}else{
		    pmodc=opd->p+nloc*cmod;
		    pmods=opd->p+nloc*(cmod+1);
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
    if(flag<0){//select only one mode
	dmat *opd0=opd;
	opd=dsub(opd0, 0, 0, -flag-(rmin)*(rmin+1)/2-1, 1);
	dfree(opd0);
    }else if(nmod>cmod){
	dresize(opd, nloc, cmod);
    }else if(nmod<cmod){
	error("over flow\n");
    }
    dfree(locr);
    dfree(locs);
    return opd;
}
/**
   return zernike index of radial mode ir, and azimuthal mode im. Maximum radial
   order is nr (inclusive)
 */
static lmat *zernike_index(int nr){
    if(nr<0) return 0;
    lmat *out=lnew(nr+1, nr+1);
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

   Notice that (D/r0)^(5/3) has been factored out.
*/
dmat *zernike_cov_kolmogorov(int nr){
    int nmod=(nr+1)*(nr+2)/2;
    dmat *res=dnew(nmod, nmod);
    PDMAT(res, pres);
    lmat *zind=zernike_index(nr);
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
		if(ir==0 || jr==0) continue;//ignore piston
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
    //pres[0][0]=0; //piston covariance
    lfree(zind);
    return res;
}

/**
   Diagnolize the covariance matrix and transform the mode.
*/
dmat *diag_mod_cov(const dmat *mz, /**<Modes in zernike space*/
		   const dmat *cov, /**<Covariance of zernike modes in kolmogorov (or any other)*/
		   int nmod       /**<keep nmod leading modes*/
    ){
    dmat *U=0, *Vt=0, *S=0;
    dsvd(&U, &S, &Vt, cov);
    //writebin(U,"KL_U");
    //writebin(S,"KL_S");
    dmat *kl=0;
    dmat *U2=0;
    if(nmod && nmod<U->ny){
	U2=dsub(U,0,0,0,nmod);
    }else{
	U2=dref(U);
    }
    dmm(&kl, 0, mz, U2, "nn", 1);
    dfree(U2);
    dfree(U);
    dfree(Vt);
    dfree(S);
    return kl;
}
/**
   Return cashed KL modes up to maxr radial order
 */
cell *KL_kolmogorov_cached(int maxr){
    int nmod=(maxr+1)*(maxr+2)/2;
    char fn[PATH_MAX];
    snprintf(fn, PATH_MAX, "%s/.aos/KL_kolmogorov.bin", HOME);
    char fnlock[PATH_MAX];
    snprintf(fnlock, PATH_MAX, "%s.lock", fn);
    cell *kl=0;
  redo:
    if(!exist(fnlock) && zfexist(fn)){
	kl=(cell*)readbin(fn);
    }
    if(!kl || kl->p[1]->ny<nmod){
	cellfree(kl);
	int fd=lock_file(fnlock, 0, 0);
	if(fd>=0){//succeed
	    kl=cellnew(2,1);
	    int KL_OVER=10;
	    READ_ENV_INT(KL_OVER, 0, 100);
	    int maxr2=maxr+10;
	    int nmod2=(maxr2+1)*(maxr2+2)/2;
	    int nx=ceil(sqrt(nmod2*KL_OVER/M_PI))*2+1;
	    info("KL_OVER=%d, nx=%d\n", KL_OVER,nx);
	    kl->p[0]=(cell*)mkcirloc(2, 2./(nx-1));
	    //writebin(kl->p[0], "KL_loc");
	    dmat *cov=zernike_cov_kolmogorov(maxr2);
	    dmat *modz=zernike((loc_t*)kl->p[0], 2., 0, maxr2, 0);
	    //writebin(modz, "KL_modz");
	    kl->p[1]=(cell*)diag_mod_cov(modz, cov, nmod);
	    //writebin(kl->p[1], "KL_kl");
	    dfree(modz);
	    dfree(cov);
	    writebin(kl, fn);
	    close(fd); remove(fnlock);
	}else{
	    info("waiting to lock %s:", fnlock);
	    fd=lock_file(fnlock, 1, 0);
	    info2("locked\n");
	    close (fd); remove(fnlock);
	    goto redo;
	}
    }
    return kl;
}
/**
   Create Karhunen-Loeve modes for which each mode is statistically independent
   in Kolmogorov spectrum. It first compute the covariance of zernike modes in
   kolmogorov spectrum, and then diagnolize the covariance of these modes using
   SVD and use the Unitery matrix to transform the zernike modes. Notice the
   ordering of the master zernike modes may be reordering in the SVD process.

   Since each KL mode is linear combination of zenike modes with same m, but
   equal or HIGHER radial order, we have to compute more zernike modes to have
   better accuracy. nr2 controls this overshoot.

   If the required zernike order is close to or more than the sampling of the
   grid, the process will not work or give poor results. Consequently, the
   result has to be computed in a finer grid and interpolated back to the
   coarser grid.

   We cache KL modes for a certain maximum order and interpolate to return results.
 */
dmat *KL_kolmogorov(const loc_t *loc, double D, int maxr){
    cell *klcache=KL_kolmogorov_cached(maxr);
    int nmod=(maxr+1)*(maxr+2)/2;
    dmat *modkl=dnew(loc->nloc, nmod);
    loc_t *locin=(loc_t*)klcache->p[0];
    dmat *opdin=(dmat*)klcache->p[1];
    for(int imod=0; imod<nmod; imod++){
	prop_nongrid((loc_t*)locin, opdin->p+locin->nloc*imod,
		     loc, NULL, modkl->p+loc->nloc*imod, 
		     1, 0, 0, 2./D, 0, 0);
    }
    cellfree(klcache);
    return modkl;
}
