/*
  Copyright 2009, 2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include <search.h>
#include <stdlib.h>
#include "laplacian.h"
double laplacian_coef(double r0, double weight, double dx){
    /*compute the factor for the laplacian*/
    double cf;
    cf=(2.*M_PI/0.5e-6)*sqrt(pow(r0/dx,5./3.)/3.44/weight)
	*pow(2.-pow(2,-1./3.)-pow(2,-1./6.),-1./2.);
    /*Fixme: need to study why is this.*/
    return cf;
}
static double laplacian_coef4(double r0, double weight, double dx){
    /*compute the factor for the laplacian*/
    double cf;
    cf=(2.*M_PI/0.5e-6)*sqrt(pow(r0/dx,5./3.)/3.44/weight)
	*pow(2.-pow(2,-4./3.)-pow(2,-1./6.),-1./2.);
    /*Fixme: need to study why is this.*/
    return cf;
}
static double laplacian_coef3(double r0, double weight, double dx){
    /*compute the factor for the laplacian*/
    double cf;
    cf=(2.*M_PI/0.5e-6)*sqrt(pow(r0/dx,5./3.)/3.44/weight)
	*pow(2.-pow(2,-1./6.),-1./2.);
    /*Fixme: need to study why is this.*/
    return cf;
}
void apply_laplacian_map(int nx, int ny, double dx, double r0, double weight, 
			 double *opd, double *opdout){
    /*Apply L2 directly to map with periodic condition.*/
    int ix,iy;
    double (*OPD)[nx]=(double(*)[nx])opd;
    double (*OPDout)[nx]=(double(*)[nx])opdout;
    double cf=laplacian_coef(r0, weight, dx);
    if(!opdout)
	error("opdout is not allocated\n");
    for(iy=0; iy<ny; iy++){
	for(ix=0; ix<nx; ix++){
	    OPDout[iy][ix]=cf*(-OPD[iy][ix]
			       +0.25*(OPD[iy==0?ny-1:iy-1][ix]
				      +OPD[iy][ix==0?nx-1:ix-1]
				      +OPD[iy][ix==nx-1?0:ix+1]
				      +OPD[iy==ny-1?0:iy+1][ix]));
	    
	}
    }
}
dsp* mklaplacian_map(int nx, int ny, double dx, double r0, double weight){
    /*build laplacian on square map using periodic conditions*/
    dsp *L2=spnew(nx*ny,nx*ny,nx*ny*5);
    int iy,ix;
    spint *pp=L2->p;
    spint *pi=L2->i;
    double *px=L2->x;
    double *px0=px;
    double cf=laplacian_coef(r0, weight, dx);
    for(iy=0; iy<ny; iy++){
	for(ix=0; ix<nx; ix++){
	    *(pp++)=px-px0;
	    *(pi++)=(iy==0?ny-1:iy-1)*nx+ix;
	    *(px++)=0.25*cf;
	    *(pi++)=iy*nx+(ix==0?nx-1:ix-1);
	    *(px++)=0.25*cf;
	    *(pi++)=iy*nx+ix;
	    *(px++)=-1*cf;
	    *(pi++)=iy*nx+(ix==nx-1?0:ix+1);
	    *(px++)=0.25*cf;
	    *(pi++)=(iy==ny-1?0:iy+1)*nx+ix;
	    *(px++)=0.25*cf;
	}
    }
    *(pp++)=px-px0;
    dsp *L2r=sptrans(L2);
    spfree(L2);
    return L2r;
}
 
dsp* mklaplacian_loc(loc_t *loc, double r0, double weight){
    /*
      The laplacian here is l=delta^2(u)/4;
      There is a scaling factor of 4.
      if I take an opd of x^2+y^2, the laplacian is 1. but the value should be 4.
      also, if the opd is cut off at some radius, the computed curvature along the edge is large number. this is not good feature.

      The result is not symmetrical for edge pixels.
      set USE_PARTIAL to 0 will result in reasonable boundary conditions.
     */

#define USE_PARTIAL 1 
    /*
      Must set USE_PARTIAL to 1. Otherwise, the SCAO, NGS, Phy, NF, case is worse than LAOS.
      USE_PARTIAL == 1 gives identical results to LAOS (upto round off error)
     */
    loc_create_map_npad(loc,1);//must have padding of at least 1.
    dsp *L2;
    L2=spnew(loc->nloc,loc->nloc,loc->nloc*5);
    int ix,iy;
    int nx1=loc->map->nx-1;
    int ny1=loc->map->ny-1;
    spint *pp=L2->p;
    spint *pi=L2->i;
    double  *px=L2->x;
    double  *px0=L2->x;
    long (*map)[loc->map->nx]=(long (*)[loc->map->nx])loc->map->p;
    double cf=laplacian_coef(r0, weight, loc->dx);
#if USE_PARTIAL == 1
    double cfs[5];
    cfs[0]=cfs[1]=0;
    cfs[2]=laplacian_coef3(r0, weight, loc->dx);
    cfs[3]=laplacian_coef4(r0, weight, loc->dx);
    cfs[4]=cf;
#endif
    for(iy=0; iy<loc->map->ny; iy++){
	for(ix=0; ix<loc->map->nx; ix++){
	    if (map[iy][ix]){
		*(pp++)=px-px0;
#if USE_PARTIAL == 1
		double *px1=px;
		double *p1, *p2;
		if(iy>0 && map[iy-1][ix]){
		    *(pi++)=map[iy-1][ix]-1;
		    *(px++)=0.25;
		    p2=px-1;
		}else{
		    p2=NULL;
		}
		if(ix>0 && map[iy][ix-1]){
		    *(pi++)=map[iy][ix-1]-1;
		    *(px++)=0.25;
		    p1=px-1;
		}else{
		    p1=NULL;
		}
		*(pi++)=map[iy][ix]-1;
		*(px++)=-1;
		if(ix<nx1 && map[iy][ix+1]){
		    *(pi++)=map[iy][ix+1]-1;
		    if(p1)
			*(px++)=0.25;
		    else
			*(px++)=0.5;
		}else{
		    if(p1)
			*p1=0.5;
		    else
			warning("Point is isolated");
		}
		if(iy<ny1 && map[iy+1][ix]){
		    *(pi++)=map[iy+1][ix]-1;
		    if(p2)
			*(px++)=0.25;
		    else
			*(px++)=0.5;
		}else{
		    if(p2)
			*p2=0.5;
		    else
			warning("Point is isolated");
		}
		double cfi=cfs[px-px1-1];
		for(double *px2=px1; px2<px; px2++){
		    *px2*=cfi;
		}
#else
		int has_up,has_rt;
		if(iy>0 && map[iy-1][ix] && iy<ny1 && map[iy+1][ix]){
		    *(pi++)=map[iy-1][ix]-1;
		    *(px++)=0.25*cf;
		    has_up=1;
		}else{
		    has_up=0;
		}
		if(ix>0 && map[iy][ix-1] && ix<nx1 && map[iy][ix+1]){
		    *(pi++)=map[iy][ix-1]-1;
		    *(px++)=0.25*cf;
		    has_rt=1;
		}else{
		    has_rt=0;
		}
		if(has_rt || has_up){
		    *(pi++)=map[iy][ix]-1;
		    *(px++)=(-0.5*has_rt-0.5*has_up)*cf;
		}
		if(has_rt){
		    *(pi++)=map[iy][ix+1]-1;
		    *(px++)=0.25*cf;
		}
		if(has_up){
		    *(pi++)=map[iy+1][ix]-1;
		    *(px++)=0.25*cf;
		}

#endif
	    }
	}
    }
    *pp=px-px0;
    if(px-px0 > L2->nzmax){
	error("Over flow happened\n");
    }
    spsetnzmax(L2,px-px0);
    //spcheck(L2);
    //L2->px may be reallocated. so scale before setnzmax.
    dsp *L2r=sptrans(L2);
    spfree(L2);
    loc_free_map(loc);
    //spclean(L2r);
    return L2r;
}
