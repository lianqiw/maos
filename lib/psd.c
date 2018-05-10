/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "../math/mathdef.h"
/**
   Compute the PSD from a sequence.
*/
//#define W_J(i,N2) (1-pow((double)(i-N2)/(double)N2, 2))
#define W_J(i,N2) 1
dmat *psd1d(const dmat *v, /**<[in] The data sequence*/
	    long nseg      /**<[in] Number of overlapping segments*/
    ){
    long nx;
    long ncol;
    if(v->nx==1){
	nx=v->ny;
	ncol=1;
    }else{
	nx=v->nx;
	ncol=v->ny;
    }
    if(nseg<=1) nseg=1;
    const int lseg2=nx/(nseg+1);
    const int lseg=lseg2*2;
    dmat *psd=dnew(lseg2+1, ncol);
    cmat *hat=cnew(lseg, 1);
    //cfft2plan(hat, -1);
    for(long icol=0; icol<ncol; icol++){
	double *ppsd=psd->p+icol*(lseg2+1);
	for(int iseg=0; iseg<nseg; iseg++){
	    double* p=v->p+icol*nx+iseg*lseg2;
	    for(int ix=0; ix<lseg; ix++){
		hat->p[ix]=p[ix]*W_J(ix, lseg2);
	    }
	    cfft2(hat, -1);
	    ppsd[0]+=cabs2(hat->p[0]);
	    for(int ix=1; ix<lseg2; ix++){
		ppsd[ix]+=cabs2(hat->p[ix])+cabs2(hat->p[lseg-ix]);
	    }
	    ppsd[lseg2]+=cabs2(hat->p[lseg2]);
	}
    }
    double sumwt=0;
    for(int ix=0; ix<lseg; ix++){
	sumwt+=pow(W_J(ix, lseg2), 2);
    }
    sumwt*=lseg*nseg;
    dscale(psd, 1./sumwt);
    cfree(hat);
    return psd;
}
/**
  Wrap of psd1d to put the frequency along the first column.
*/
dmat *psd1dt(const dmat *v, long nseg, double dt){
    dmat *psd=psd1d(v, nseg);
    dmat *psd2=dnew(psd->nx, psd->ny+1);
    int N=(psd->nx-1)*2;
    double df=1./(N*dt);
    for(int i=0; i<psd->nx; i++){
	psd2->p[i]=df*i;
    }
    dscale(psd, 1./df);//divide so the value is point, not integrated in a bin.
    memcpy(psd2->p+psd2->nx, psd->p, psd->nx*psd->ny*sizeof(double));
    dfree(psd);
    return psd2;
}

/*Interpolate psd onto new f. We interpolate in log space which is more linear.*/
dmat *psdinterp1(const dmat *psdin, const dmat *fnew, int uselog){
    dmat *f1=drefcols(psdin, 0, 1);
    dmat *psd1=dsub(psdin, 0, 0, 1, 1);//copy
    dmat *f2=dref(fnew);
    double t1=dtrapz(f1, psd1);
    double ydefault=1e-40;
    if(uselog){
	dcwlog(psd1);
	ydefault=log(ydefault);
    }
    dmat *psd2=dinterp1(f1, psd1, f2, ydefault);
    if(uselog){
	dcwexp(psd2,1);
    }
    double t2=dtrapz(f2, psd2);
    if(fabs(t1-t2)>fabs(0.5*(t1+t2)*2)){
	warning("psd interpolation failed. int_orig=%g, int_interp=%g\n", t1, t2);
    }
    //Don't scale psd2 as it may have overlapping frequency regions
    dfree(f1); dfree(f2); dfree(psd1);
    return psd2;
}
/**
   Find vibration peaks in the PSD by comparing the PSD against a LPF version plus noise.
 */
dmat *psd_vibid(const dmat *psdin){
    double *f=psdin->p;
    double *psd=psdin->p+psdin->nx;
    dmat *y=dsub(psdin, 0, 0, 1, 1);
    const double gain=0.1;
    const double gain2=0.05;
    int inpeak=0;
    double ylpf0=y->p[1];
    double dylpf0=fabs(y->p[1]-y->p[0]);
    double ylpf=0, dylpf=0;
    int nmaxp=100;
    dmat *res=dnew(4, nmaxp);
    double thres=25e-18;/*threshold: 5 nm*/
    double sumxy=0, sumy=0, sum=0;
    int count=0;
    for(long i=1; i<psdin->nx-1; i++){
	if(!inpeak){
	    //second order LPF
	    ylpf0=(1.-gain)*ylpf0+y->p[i]*gain; 
	    ylpf=(1.-gain)*ylpf+ylpf0*gain;
	    double diff=y->p[i]-y->p[i-1];
	    if(diff>0){
		dylpf0=(1.-gain2)*dylpf0+diff*gain2; 
		dylpf=(1.-gain2)*dylpf+dylpf0*gain2;
	    }
	    if(y->p[i+1]>ylpf+dylpf*5 && f[i]>1){//beginning of peak
		inpeak=1;
		if(count>0 && f[i] < f[(int)IND(res,3,count-1)] + 0.1){
		    //combine with last peak if within 1 Hz.
		    count--;
		}else{
		    IND(res,2,count)=i;
		    sumxy=f[i]*psd[i];//for CoG
		    sumy=psd[i];//for CoG
		    sum=0;//integration
		}
	    }
	}else{
	    //continuation of peak
	    sumxy+=f[i]*psd[i];
	    sumy+=psd[i];
	    sum+=(f[i]-f[i-1])*(psd[i]+psd[i-1]);
	    if(y->p[i]<ylpf+dylpf && y->p[i+1]<ylpf+dylpf){//end of peak
		inpeak=0;
		if(sum*0.5>thres){
		    IND(res,0,count)=sumxy/sumy;
		    IND(res,1,count)=sum*0.5;
		    IND(res,3,count)=i;
		    count++;
		    if(count==nmaxp){
			nmaxp*=2;
			dresize(res, 4, nmaxp);
		    }
		}
	    }
	}
    }
    dfree(y);
    dresize(res, 4, count);
    return res;
}

/**Convert temporal PSD to spatial*/
dmat *psdt2s(const dmat *psdt, double vmean){
    if(psdt->nx!=1 || psdt->ny>4){
	error("psdt needs to be 1 row and less than 4 cols\n");
    }
    double alpha=psdt->p[0];//power
    double beta=psdt->p[1];//strength
    double alpha2=alpha-1;
    //n is -alpha2
    double bfun=tgamma(0.5)*tgamma((-alpha2-1)*0.5)/tgamma(-alpha2*0.5);
    double beta2=beta/bfun*pow(vmean, 2+alpha2);
    dmat *psds=dnew(psdt->nx, psdt->ny);
    psds->p[0]=alpha2;
    psds->p[1]=beta2;
    if(psds->ny>2){
	psds->p[2]=psdt->p[2]/vmean;
    }
    if(psds->ny>3){
	psds->p[3]=psdt->p[3]/vmean;
    }
    return psds;
}

/**Convert special PSD to temporal*/
dmat *psds2t(const dmat *psds, double vmean){
    if(psds->nx!=1 || psds->ny>4){
	error("psds needs to be 1 row and less than 4 cols\n");
    }
    double alpha=psds->p[0];//power
    double beta=psds->p[1];//strength
    double alpha2=alpha+1;
    //n is -alpha
    double bfun=tgamma(0.5)*tgamma((-alpha-1)*0.5)/tgamma(-alpha*0.5);
    double beta2=beta*bfun*pow(vmean, -2-alpha);
    dmat *psdt=dnew(psds->nx, psds->ny);
    psdt->p[0]=alpha2;
    psdt->p[1]=beta2;
    if(psds->ny>2){
	psdt->p[2]=psds->p[2]*vmean;
    }
    if(psds->ny>3){
	psdt->p[3]=psds->p[3]*vmean;
    }
    return psdt;
}


/**
   Integrated a PSF that defines on linear or logrithmically spaced grid nu.
*/
double psd_inte(const double *nu, const double *psd, long n){
    double dnu=(nu[n-1]-nu[0])/(n-1);
    double dlognu=(log(nu[n-1])-log(nu[0]))/(n-1);
    double res_sig=0;
    if(fabs(nu[1]-nu[0]-dnu)<dnu*1.e-4){
	for(long i=0; i<n; i++){
	    res_sig+=psd[i];
	}
	res_sig*=dnu;
    }else if((log(nu[1])-log(nu[0])-dlognu)<dlognu*1.e-4){
        for(long i=0; i<n; i++){
	    res_sig+=psd[i]*nu[i];
	}	
	res_sig*=dlognu;
    }
    return res_sig;
}
/**
   wraps psd_inte
*/
double psd_inte2(const dmat *psdin){
    if(psdin->ny!=2){
	error("psdin  should have two columns\n");
    }
    double *nu=psdin->p;
    long n=psdin->nx;
    double *psd=nu+n;
    return psd_inte(nu, psd, n);
}

/**
   Convert PSD into time series.*/
dmat* psd2time(const dmat *psdin, rand_t *rstat, double dt, int nstepin){
    if(!psdin){
	error("psdin cannot be null\n");
    }
    long nstep=nextpow2(nstepin);
    double df=1./(dt*nstep);
    dmat *fs=dlinspace(0, df, nstep);
    dmat *psd=NULL;
    if(psdin->ny==1){//[alpha, beta, fmin, fmax] discribes power law with cut on/off freq.
	psd=dnew(nstep, 1);
	double alpha=psdin->p[0];
	double beta=psdin->p[1];
	long i0=1, imax=nstep;
	if(psdin->nx>2){
	    i0=(long)round(psdin->p[2]/df);
	    if(i0<1) i0=1;
	}
	if(psdin->nx>3){
	    imax=(long)round(psdin->p[3]/df);
	}
	dbg("fmin=%g, fmax=%g, df=%g, i0=%ld, imax=%ld\n", 
	     psdin->p[2], psdin->p[3], df, i0, imax);
	for(long i=i0; i<imax; i++){
	    psd->p[i]=beta*pow(i*df, alpha);
	}
    }else if(psdin->ny==2){
	if(psdin->nx<2){ 
	    error("Invalid PSD\n");
	}
	psd=dinterp1(psdin, 0, fs, 1e-40);
	psd->p[0]=0;/*disable pistion. */
    }else{
	error("psdin is invalid format.\n");
    }
    cmat *wshat=cnew(nstep, 1);
    //cfft2plan(wshat, -1);
    for(long i=0; i<nstep; i++){
	wshat->p[i]=sqrt(psd->p[i]*df)*COMPLEX(randn(rstat), randn(rstat));
    }
    cfft2(wshat, -1);
    dmat *out=NULL;
    creal2d(&out, 0, wshat, 1);
    cfree(wshat);
    dfree(psd);
    dfree(fs);
    dresize(out, nstepin, 1);
    return out;
}

/**
   Add two PSDs that doesn't have the same frequency. the first column of each
   dmat is the frequency nu, and the second column is PSD. Bug discovered on
   2013-03-24:only psd2 was added to to psd.*/
static dmat *add_psd_nomatch(const dmat *psd1,const dmat *psd2){
    dmat *nu1=dsub(psd1,0,psd1->nx,0,1);
    dmat *p2ynew=dinterp1(psd2, 0, nu1, 1e-40);
    dmat *psd=dnew(nu1->nx,2);
    double *py=psd->p+psd->nx;
    const double *p1y=psd1->p+psd1->nx;
    for(long i=0; i<psd->nx; i++){
	psd->p[i]=nu1->p[i];
	py[i]=p1y[i]+p2ynew->p[i];
    }
    dfree(nu1);
    dfree(p2ynew);
    return psd;
}

/**
   Add two PSDs. the first column of each dmat is the frequency nu, and the
   second column is PSD*/
dmat *add_psd(const dmat *psd1, const dmat *psd2){
    if(psd1->nx!=psd2->nx){
	//warning("The two PSDs have different length\n");
	return add_psd_nomatch(psd1, psd2);
    }
    dmat *psd=dnew(psd1->nx,2);
    double *restrict pp=psd->p+psd->nx;
    const double *p1=psd1->p+psd1->nx;
    const double *p2=psd2->p+psd2->nx;
    for(long i=0; i<psd->nx; i++){
	if(fabs(psd1->p[i]-psd2->p[i])>1.e-2){
	    warning("The two PSDs doesn't have the same freq.");
	    dfree(psd);
	    return add_psd_nomatch(psd1,psd2);
	    /*todo: apply interp1 to interpolate the second PSD. */
	}
	psd->p[i]=psd1->p[i];
	pp[i]=p1[i]+p2[i];
    }
    return psd;
}
/*
  Add a PSD to another.
*/
void add_psd2(dmat **out, const dmat *in){
    if(!*out){
	*out=ddup(in);
    }else{
	dmat *tmp=*out;
	*out=add_psd(tmp, in);
	dfree(tmp);
    }
}
/**
   Sum all the columns of PSDs (excluding first column), scaled by scale
*/
void psd_sum(dmat *psd, double scale){
    for(int ix=0; ix<psd->nx; ix++){
	double tmp=0;
	for(int iy=1; iy<psd->ny; iy++){
	    tmp+=IND(psd, ix, iy);
	}
	IND(psd, ix, 1)=tmp*scale;
    }
    dresize(psd, psd->nx, 2);
}
