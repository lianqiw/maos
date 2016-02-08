/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "turbulence.h"
/**
   Generating Zernike Rnm for radial order ir, azimuthal or im.  used by loc_zernike.
   */
dmat *zernike_Rnm(const dmat *locr, int ir, int im){
    if(ir<0 || im < 0|| im>ir || (ir-im)%2!=0)
	error("Invalid ir, im (%d, %d)\n", ir, im);
    const long nloc=locr->nx*locr->ny;
    dmat *Rnm=dnew(locr->nx, locr->ny);
    const int ns=(ir-im)/2+1;
    double coeff[ns];
    int power[ns];
    for(int s=0; s<ns; s++){
	coeff[s]=factoriall((ir+im)/2-s+1, ir-s)/factoriall(1,s)/factoriall(1,(ir-im)/2-s)*pow(-1,s);
	power[s]=ir-2*s;
    }
    for(long iloc=0; iloc<nloc; iloc++){
	double tmp=0;
	double r=locr->p[iloc];
	for(int s=0; s<ns; s++){
	    tmp+=coeff[s]*pow(r,power[s]);
	}
	Rnm->p[iloc]=tmp;
    }
    return Rnm;
}

/**
   Create Zernike modes on loc with diameter D and radial order from rmin to rmax
   r=0 is piston
   r=1 is tip/tilt
   r=2 is quadratic modes
   if nopiston is set, skip piston mode.
   if flag is >0, only radial mode is used.
   if flag is  0, all modes between rmin and rmax.
   if flag is <0, only generate mode -flag. rmin and rmax is irrelevant
*/
dmat* zernike(const loc_t *loc, double D, int rmin, int rmax, int flag){
    if(flag<0){
	rmin=ceil((sqrt(8.*(-flag)+1)-3)*0.5);
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
    }else{//a specific mode
	nmod=(rmax+1)*(rmax+2)/2-(rmin)*(rmin+1)/2;
    }
    double D2=loc_diam(loc);
    if(D<=0){
	D=D2;
    }else if(fabs(D-D2)>D*0.5){
	warning("specified diameter is incorrect. D=%g, loc D=%g\n", D, D2);
    }
 
    dmat *restrict opd=dnew(nloc,nmod);
    dmat *restrict locr=dnew(nloc,1);
    dmat *restrict locs=dnew(nloc,1);
    const double *restrict locx=loc->locx;
    const double *restrict locy=loc->locy;
    const double R1=2./D;
    long nover=0; double rover=1;
    for(long iloc=0; iloc<nloc; iloc++){
	locr->p[iloc]=sqrt(pow(locx[iloc],2)+pow(locy[iloc],2))*R1;
	if(locr->p[iloc]>1){
	    nover++;
	    if(locr->p[iloc]>rover){
		rover=locr->p[iloc];
	    }
	}
	locs->p[iloc]=atan2(locy[iloc], locx[iloc]);
    }
    if(nover > (M_PI*D/loc->dx)){
	warning("%ld/%ld points outside unit circle with maximum radius %g\n",
		nover, nloc, rover);
    }
    int cmod=0;
    for(int ir=rmin; ir<=rmax; ir++){
	for(int im=0; im<=ir; im++){
	    if((ir-im)%2!=0) continue;
	    if(flag>0 && im!=0) continue;//we want radial only
	    dmat *Rnm=zernike_Rnm(locr, ir, im);
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

   2015-06-05: This is not longer being used. It is inaccurate for annular
   pupil. Analytic calculation using the power spectrum is used instead.
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
   Compute covariance matrix of zernike modes in von Karman turbulence
*/
dmat *cov_vonkarman(const loc_t *loc, /**<The location grid*/
		    const dmat *modz, /**<Zernike modes*/
		    double L0 /**<Outer scale*/){
    dmat *CC=0;
    dmm(&CC, 1, modz, modz, "tn", 1);//the covariance of the modes
    long nembed=0;
    lmat *embed=loc_create_embed(&nembed, loc, 2, 0);
    int nmod=modz->ny;
    ccell *spect=ccellnew(nmod, 1);
#pragma omp parallel for
    for(long ic=0; ic<nmod; ic++){
	spect->p[ic]=cnew(nembed, nembed);
	for(long ix=0; ix<loc->nloc; ix++){
	    spect->p[ic]->p[embed->p[ix]]=IND(modz, ix, ic);
	}
	cfftshift(spect->p[ic]);
	cfft2(spect->p[ic], -1);
    }
    dmat *turbspec=turbpsd(nembed, nembed, loc->dx, 0.2, L0, -11./3., 1);
    turbspec->p[0]=0;//remove piston.
    dmat *DD=dnew(nmod, nmod);
#pragma omp parallel for
    for(long ic=0; ic<nmod; ic++){
	for(long id=0; id<=ic; id++){
	    double tmp=0;
	    for(long ip=0; ip<nembed*nembed; ip++){
		tmp+=spect->p[ic]->p[ip]*conj(spect->p[id]->p[ip])*turbspec->p[ip];
	    }
	    IND(DD, ic, id)=IND(DD,id,ic)=tmp;//*scale;
	}
    }
    dmat *CCi=dpinv(CC, 0);
    dmat *tmp=0;
    dmm(&tmp, 0, CCi, DD, "nn", 1);
    dmm(&DD, 0, tmp, CCi, "nt", 1);
    dfree(tmp);
    lfree(embed);
    dfree(CCi);
    dfree(turbspec);
    ccellfree(spect);
    return DD;
}

/**
   Diagnolize the covariance matrix and transform the mode.
*/
dmat *cov_diagnolize(const dmat *mod, /**<Input mode*/
		     const dmat *cov  /**<Covariance of modes*/
    ){
    dmat *U=0, *Vt=0, *S=0;
    dsvd(&U, &S, &Vt, cov);
    dmat *kl=0;
    dmm(&kl, 0, mod, U, "nn", 1);
    {
	//Drop modes with infinitesimal strength
	double ssmax=S->p[0];
	double thres=ssmax*1e-10;
	long count=0;
	for(long i=S->nx-1; i>0; i--){
	    if(S->p[i]<thres){
		count++;
	    }else{
		break;
	    }
	}
	if(count>0){
	    warning("Drop %ld last columns.\n", count);
	    dmat *kl2=dsub(kl, 0, 0, 0, kl->ny-count);
	    dfree(kl);
	    kl=kl2;
	}
    }
    dfree(U);
    dfree(Vt);
    dfree(S);
    return kl;
}

/**
   Create Karhunen-Loeve modes for which each mode is statistically independent
   in von Karman spectrum.
   
   The original method is to compute covariance of zernike modes in von Karman
   spectrum and diagnolize the covariance matrix. But this methods suffers for
   higher order systems because high order zernike modes are hard to generate on
   account of limited resolution of double precision floating point numbers. To
   make matters worse, to generate m KL modes, around 2*m Zernike modes are
   needed because each KL mode is a linear combination of all Zernike modes with
   same azimuthal (m) but higher radial (r) order. A work around of the
   precisiosn issue is to limit the absolute value of the zernike modes to
   within 100 to avoid excessive norm due to round off errors. 

   The new method is to use zonal modes (identity matrix). The advantage is that
   the zonal modes spans the whole range of possible modes and are
   straightforward to generate.

   In theory, the KL modes computes should be independent on the starting mode,
   as long as the modes span the whole vector space for the coordinate.
 */
dmat *KL_vonkarman_do(const loc_t *loc, double L0){
    dmat *modz=dnew(loc->nloc, loc->nloc);
    double val=sqrt(loc->nloc); //this ensures rms wfe is 1, or orthonormal.
    daddI(modz, val);
    dadds(modz, -val/loc->nloc);//this ensure every column sum to 0 (no piston)
    dmat *cov=cov_vonkarman(loc, modz, L0);
    dmat *kl=cov_diagnolize(modz, cov);
    dfree(modz);
    dfree(cov);
    return kl;
}
dmat *KL_vonkarman(const loc_t *loc, int nmod, double L0){
    if(!loc) return 0;
    uint32_t key=lochash(loc, 0);
    mymkdir("%s/.aos/cache", HOME);
    char fn[PATH_MAX];
    snprintf(fn, PATH_MAX, "%s/.aos/cache/KL_vonkarman_%u_%g_%ld.bin", HOME, key, L0, loc->nloc);
    char fnlock[PATH_MAX];
    snprintf(fnlock, PATH_MAX, "%s.lock", fn);
    dmat *kl=0;
  redo:
    if(loc->nloc<500){//fast. no need cache.
	kl=KL_vonkarman_do(loc, L0);
    }else if(!exist(fnlock) && zfexist(fn)){
	kl=dread("%s", fn);
    }else{
	info("trying to lock %s\n", fnlock);
	int fd=lock_file(fnlock, 0, 0);
	if(fd>=0){//start preparing
	    info2("locked\n");
	    kl=KL_vonkarman_do(loc, L0);
	    writebin(kl, "%s", fn);
	    close(fd); remove(fnlock);
	}else{
	    info("waiting to lock %s\n", fnlock);
	    fd=lock_file(fnlock, 1, 0);
	    info2("locked\n");
	    close (fd); remove(fnlock);
	    goto redo;
	}
    }
    if(nmod>0 && nmod<loc->nloc){
	//Take sub matrix
	dmat *klorig=kl;
	kl=dsub(klorig, 0, 0, 0, nmod);
	dfree(klorig);
    }
    return kl;
}
