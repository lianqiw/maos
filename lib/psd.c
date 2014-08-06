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
#include "../math/mathdef.h"
/**
   Compute the PSD from a sequence.
*/
//#define W_J(i,N2) (1-pow((double)(i-N2)/(double)N2, 2))
#define W_J(i,N2) 1
dmat *psd1d(dmat *v, /**<[in] The data sequence*/
	   long lseg /**<[in] The length of overlapping segments*/
	   ){
    long lseg2=lseg>>1;
    if(v->nx==1){
	v->nx=v->ny;
	v->ny=1;
    }
    int nseg=(v->nx/lseg2-1)>>1; /*number of segments*/
    if(nseg<1){
	nseg=1;
	lseg=v->nx;
	lseg2=lseg>>1;
    }
    long ncol=v->ny;
    dmat *psd=dnew(lseg2+1, ncol);
    cmat *hat=cnew(lseg, 1);
    cfft2plan(hat, -1);
    for(long icol=0; icol<ncol; icol++){
	double *ppsd=psd->p+icol*(lseg2+1);
	for(int iseg=0; iseg<nseg; iseg++){
	    double* p=v->p+icol*v->nx+iseg*lseg2;
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
dmat *psd1dt(dmat *v, long lseg, double dt){
    dmat *psd=psd1d(v, lseg);
    dmat *psd2=dnew(psd->nx, psd->ny+1);
    int N=(psd->nx-1)*2;
    double df=1./(N*dt);
    for(int i=0; i<psd->nx; i++){
	psd2->p[i]=df*i;
    }
    dscale(psd, 1./df);//divide so the value is point, not integrated in a bin.
    memcpy(psd2->p+psd->nx, psd->p, psd->nx*psd->ny*sizeof(double));
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
    double ym0=y->p[1];
    double yn0=fabs(y->p[1]-y->p[0]);
    double ym=0, yn=0;
    int nmaxp=100;
    dmat *res=dnew(4, nmaxp);
    PDMAT(res, pres);
    double thres=25e-18;/*threshold: 5 nm*/
    double sumxy, sumy, sum;
    int count=0;
    for(long i=1; i<psdin->nx-1; i++){
	if(!inpeak){
	    //second order LPF
	    ym0=(1.-gain)*ym0+y->p[i]*gain; 
	    ym=(1.-gain)*ym+ym0*gain;
	    double diff=y->p[i]-y->p[i-1];
	    if(diff>0){
		yn0=(1.-gain2)*yn0+diff*gain2; 
		yn=(1.-gain2)*yn+yn0*gain2;
	    }
	    if(y->p[i+1]>ym+yn*5 && f[i]>1){//beginning of peak
		inpeak=1;
		if(count>0 && f[i] < f[(int)pres[count-1][3]] + 0.1){
		    //combine with last peak if within 1 Hz.
		    count--;
		}else{
		    pres[count][2]=i;
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
	    if(y->p[i]<ym+yn && y->p[i+1]<ym+yn){//end of peak
		inpeak=0;
		if(sum*0.5>thres){
		    pres[count][0]=sumxy/sumy;
		    pres[count][1]=sum*0.5;
		    pres[count][3]=i;
		    count++;
		    if(count==nmaxp){
			nmaxp*=2;
			dresize(res, 4, nmaxp);
			pres=(void*)res->p;
		    }
		}
	    }
	}
    }
    dfree(y);
    dresize(res, 4, count);
    return res;
}
