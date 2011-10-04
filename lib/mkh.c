/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#include <math.h>
#include <alloca.h>
#include "common.h"
#include "mkh.h"
#include "misc.h"
#include "loc.h"

#define ONLY_FULL 0 /*1: comply with other propagators. */
/*
  ONLY_FULL==0: calc weights if not all points exist. preferred
  ONLY_FULL==1: calc weights only if all points exist
 */
/**
   \file mkh.c
   Contains functions that create ray tracing operator
*/
static dsp *mkhb_cubic(loc_t *locin, loc_t *locout, const double *ampout,
		       double displacex, double displacey, double scale,double cubic_iac);
/**
   Create ray tracing operator from coordinate locin to locout.  Locin is
   required to be evenly spaced.

   If vector Pin is defined on locin, Pout is defined on locout, H=mkh(locin,
   locout, ...), Pout=H*Pin does the bilinear interpolation.

   If cubic is 1, will call mkh_cubic to produce a cubical interpolation.
   
   A cubic influence function that can reproduce piston/tip/tilt is coined by
   Ellerbroek to model the piezostack DM actuator. The influence has the
   form\f[
   h(x;x_{i};\delta)=h_{0}((x-x_{i})/\delta)
   \f]
   where \f$\delta\f$ is the same as the grid spacing, and \f$h_{0}\f$ is the influence
   function defined in the normalized coordinates\f{eqnarray*}{
   h_{0}(x)=\frac{1}{1+2c}\begin{cases}
   1+(4c-\frac{5}{2})|x|^{2}+(\frac{3}{2}-3c)|x|^{3} & |x|\leq1\\
   (2c-\frac{1}{2})(2-|x|)^{2}+(\frac{1}{2}-c)(2-|x|)^{3} & 1<|x|\leq2\\ 0 &
   |x|>2\end{cases}\f} where c is the nearest neighbor coupling
   frequency. The leading coefficient is to normalize the influence function so
   that it sums to 1.

 */
dsp* mkh(loc_t *locin, loc_t *locout, const double *ampout,
	 double displacex, double displacey, double scale,
	 int cubic, double cubic_iac){
    dsp *Hb=mkhb(locin, locout, ampout,displacex, displacey, scale, cubic, cubic_iac);
    dsp *H=sptrans(Hb);
    spfree(Hb);
    return H;
}
/**
   Create transpose of mkh() result.
*/
dsp* mkhb(loc_t *locin, loc_t *locout, const double *ampout,
	  double displacex, double displacey, double scale,
	  int cubic, double cubic_iac){
    if(cubic){
	return mkhb_cubic(locin, locout, ampout, displacex, displacey, scale, cubic_iac);
    }
  
    loc_create_map_npad(locin,1);/*will only do once and save in locin. */
    dsp *hback;
    double dplocx, dplocy;
    int nplocx, nplocy, nplocx1, nplocy1;
    long iloc;
    /*int missing=0; */
    const int wrapx = locin->map->nx-1;
    const int wrapy = locin->map->ny-1;
    const double dx_in1 = 1./locin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-locin->map->ox)*dx_in1;
    displacey = (displacey-locin->map->oy)*dx_in1;
    const double *px=locout->locx;
    const double *py=locout->locy;
#if ONLY_FULL==1
    long iphi1,iphi2,iphi3,iphi4;
#else
    long iphi;
#endif
    /*-1 because we count from 1 in the map. */
    long (*map)[locin->map->nx] =(long(*)[locin->map->nx])(locin->map->p);
    /*transpose of hfor */
    long nzmax=locout->nloc*4;
    hback = spnew(locin->nloc, locout->nloc, nzmax);
    spint *bp=hback->p;
    spint *bi=hback->i;
    double *bx=hback->x;
    long count=0;
    double weight;
    /*double *phiin0=phiin-1; */
    for(iloc=0; iloc<locout->nloc; iloc++){
	bp[iloc]=count;/*column index */
	if(ampout && fabs(ampout[iloc])<EPS)
	    continue;
	if(count+5>nzmax){
	    nzmax*=2;
	    spsetnzmax(hback, nzmax);
	    bp=hback->p;
	    bi=hback->i;
	    bx=hback->x;
	}
	dplocy=myfma(py[iloc],dx_in2,displacey);
	dplocx=myfma(px[iloc],dx_in2,displacex);

	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);

	if(nplocx<0||nplocx>=wrapx||nplocy<0||nplocy>=wrapy){
	    continue;
	}
	nplocx1=nplocx+1;
	nplocy1=nplocy+1;
#if ONLY_FULL == 1 /*only proceed if all four points exist. preferred */
	iphi1=map[nplocy][nplocx];
	iphi2=map[nplocy][nplocx1];
	iphi3=map[nplocy1][nplocx];
	iphi4=map[nplocy1][nplocx1];
	if(iphi1-- && iphi2-- && iphi3-- && iphi4--){
	    if((weight=(1.-dplocx)*(1.-dplocy))>EPS){
		bi[count]=iphi1;
		bx[count]=weight;
		count++;
	    }
	    if((weight=(dplocx)*(1.-dplocy))>EPS){
		bi[count]=iphi2;
		bx[count]=weight;
		count++;
	    }
	    if((weight=(1.-dplocx)*(dplocy))>EPS){
		bi[count]=iphi3;
		bx[count]=weight;
		count++;
	    }
	    if((weight=(dplocx)*(dplocy))>EPS){
		bi[count]=iphi4;
		bx[count]=weight;
		count++;
	    }
	}
#else	
	if((iphi=map[nplocy][nplocx]) && (weight=(1.-dplocx)*(1.-dplocy))>EPS){
	    bi[count]=iphi-1;
	    bx[count]=weight;
	    count++;
	}
	if((iphi=map[nplocy][nplocx1]) && (weight=(dplocx)*(1.-dplocy))>EPS ){
	    bi[count]=iphi-1;
	    bx[count]=weight;
	    count++;
	}
	if((iphi=map[nplocy1][nplocx]) && (weight=(1.-dplocx)*(dplocy))>EPS){
	    bi[count]=iphi-1;
	    bx[count]=weight;
	    count++;
	}
	if((iphi=map[nplocy1][nplocx1]) && (weight=(dplocx)*(dplocy))>EPS){
	    bi[count]=iphi-1;
	    bx[count]=weight;
	    count++;
	}
#endif
    }
    bp[locout->nloc]=count;
    spsetnzmax(hback, count);
    spdroptol(hback,1e-12);
    return hback;
}
/**
   Create transpose of ray tracing operator from locin to locout using cubic
   influence function that can reproduce piston/tip/tilt.  */
static dsp *mkhb_cubic(loc_t *locin, loc_t *locout, const double *ampout,
		double displacex, double displacey, double scale,double cubic_iac){
    loc_create_map_npad(locin,2);
    dsp *hback;
    double dplocx, dplocy;
    int ix,iy;
    int nplocx, nplocy;
    long iloc;
    int missing=0;
    const double dx_in1 = 1./locin->dx;
    const double dx_in2 = scale*dx_in1;
    displacex = (displacex-locin->map->ox)*dx_in1;
    displacey = (displacey-locin->map->oy)*dx_in1;
    const double *px=locout->locx;
    const double *py=locout->locy;
    long iphi;
    /*-1 because we count from 1 in the map. */
    long (*map)[locin->map->nx]
	=(long(*)[locin->map->ny])(locin->map->p);
    /*cubic */
    double fx[4],fy[4];
    const double cubicn=1./(1.+2.*cubic_iac);
    const double c0=1*cubicn;
    const double c1=(4*cubic_iac-2.5)*cubicn;
    const double c2=(1.5-3*cubic_iac)*cubicn;
    const double c3=(2*cubic_iac-0.5)*cubicn;
    const double c4=(0.5-cubic_iac)*cubicn;
    double dplocx0, dplocy0;
    long nzmax=locout->nloc*16;
    hback = spnew(locin->nloc, locout->nloc, nzmax);
  
    spint *bp=hback->p;
    spint *bi=hback->i;
    double *bx=hback->x;
    long count=0;
    const int nmapx3=locin->map->nx-3;
    const int nmapy3=locin->map->ny-3;
    for(iloc=0; iloc<locout->nloc; iloc++){
	bp[iloc]=count;
	if(ampout && fabs(ampout[iloc])<EPS)
	    continue;
	if(count+17>nzmax){
	    nzmax*=2;
	    spsetnzmax(hback, nzmax);
	    bp=hback->p;
	    bi=hback->i;
	    bx=hback->x;
	}

	dplocy=myfma(py[iloc],dx_in2,displacey);
	dplocx=myfma(px[iloc],dx_in2,displacex);

	SPLIT(dplocx,dplocx,nplocx);
	SPLIT(dplocy,dplocy,nplocy);
	dplocy0=1.-dplocy;
	dplocx0=1.-dplocx;

	fx[0]=dplocx0*dplocx0*(c3+c4*dplocx0);
	fx[1]=c0+dplocx*dplocx*(c1+c2*dplocx);
	fx[2]=c0+dplocx0*dplocx0*(c1+c2*dplocx0);
	fx[3]=dplocx*dplocx*(c3+c4*dplocx);
	    
	fy[0]=dplocy0*dplocy0*(c3+c4*dplocy0);
	fy[1]=c0+dplocy*dplocy*(c1+c2*dplocy);
	fy[2]=c0+dplocy0*dplocy0*(c1+c2*dplocy0);
	fy[3]=dplocy*dplocy*(c3+c4*dplocy); 

	if(nplocy<1 || nplocy>nmapy3 || nplocx<1 || nplocx>nmapx3){
	    continue;
	}
	for(iy=nplocy-1; iy<nplocy+3; iy++){
	    for(ix=nplocx-1; ix<nplocx+3; ix++){
		iphi=map[iy][ix];
		double weight=fx[ix-nplocx+1]*fy[iy-nplocy+1];
		if(iphi && weight>EPS){
		    bi[count]=iphi-1;
		    bx[count]=weight;
		    count++;
		}
	    }
	}
    }/*for */
    bp[locout->nloc]=count;
    spsetnzmax(hback, count);
    spdroptol(hback,1e-12);
    if(missing>0){
	warning("%d points not covered by input screen\n", missing);
    }
    return hback;
}
#undef ONLY_FULL
/**
   Create a matrix to bin from coordinate xin to xout using bilinear interpolation. xin and
   xout should be 1-d arrays of coordinates. We require the coordinates to order
   incrementally monotonically, but do not require them to be evenly spaced.
 */
dsp *mkhbin1d(dmat *xin, dmat *xout){
    if(xin->ny!=1 || xout->ny!=1){
	error("We require both xin and xout to be only one column\n");
    }
    int iout=0;
    dsp *hbin=spnew(xout->nx, xin->nx, xin->nx*2);
    int count=0;
    for(int iin=0; iin<xin->nx; iin++){
	hbin->p[iin]=count;
	double ixin=xin->p[iin];
	while(iout+1<xout->nx && xout->p[iout+1]<ixin){
	    iout++;
	}
	if(xout->p[iout]>ixin || (iout+1==xout->nx && xout->p[iout]<ixin)){/*outside of the area. */
	    hbin->i[count]=iout;
	    hbin->x[count]=1;
	    count++;
	}else{/*within the area */
	    if(iout+1==xout->nx){
		error("This shouldn't happen\n");
	    }
	    double wt=(ixin-xout->p[iout])/(xout->p[iout+1]-xout->p[iout]);
	    hbin->i[count]=iout;
	    hbin->x[count]=1.-wt;
	    count++;
	    hbin->i[count]=iout+1;
	    hbin->x[count]=wt;
	    count++;
	}
    }
    hbin->p[xin->nx]=count;
    return hbin;
}

