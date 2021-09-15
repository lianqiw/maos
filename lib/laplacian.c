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


#include <search.h>

#include "laplacian.h"

/**
   compute the factor for the laplacian*/
real laplacian_coef(real r0, real weight, real dx){
	real cf;
	cf=(2.*M_PI/0.5e-6)*sqrt(pow(r0/dx, 5./3.)/3.44/weight)
		*pow(2.-pow(2, -1./3.)-pow(2, -1./6.), -1./2.);
		/*Fixme: need to study why is this.*/
	return cf;
}
static real laplacian_coef4(real r0, real weight, real dx){
	/*compute the factor for the laplacian*/
	real cf;
	cf=(2.*M_PI/0.5e-6)*sqrt(pow(r0/dx, 5./3.)/3.44/weight)
		*pow(2.-pow(2, -4./3.)-pow(2, -1./6.), -1./2.);
		/*Fixme: need to study why is this.*/
	return cf;
}
static real laplacian_coef3(real r0, real weight, real dx){
	/*compute the factor for the laplacian*/
	real cf;
	cf=(2.*M_PI/0.5e-6)*sqrt(pow(r0/dx, 5./3.)/3.44/weight)
		*pow(2.-pow(2, -1./6.), -1./2.);
		/*Fixme: need to study why is this.*/
	return cf;
}
/**
   Apply L2 directly to map with periodic condition.*/
void apply_laplacian_map(dmat* opdout, const dmat* opd, real dx, real r0, real weight){
	const long nx=NX(opd);
	const long ny=NY(opd);
	real cf=laplacian_coef(r0, weight, dx);
	if(!opdout)
		error("opdout is not allocated\n");
	for(long iy=0; iy<ny; iy++){
		for(long ix=0; ix<nx; ix++){
			P(opdout, ix, iy)=cf*(-P(opd, ix, iy)
				+0.25*(P(opd, ix, iy==0?ny-1:iy-1)
					+P(opd, ix==0?nx-1:ix-1, iy)
					+P(opd, ix==nx-1?0:ix+1, iy)
					+P(opd, ix, iy==ny-1?0:iy+1)));

		}
	}
}
/**
   build laplacian on square map using periodic conditions*/
dsp* mklaplacian_map(int nx, int ny, real dx, real r0, real weight){
	dsp* L2=dspnew(nx*ny, nx*ny, nx*ny*5);
	int iy, ix;
	spint* pp=L2->pp;
	spint* pi=L2->pi;
	real* px=L2->px;
	real* px0=px;
	real cf=laplacian_coef(r0, weight, dx);
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
	dsp* L2r=dsptrans(L2);
	dspfree(L2);
	return L2r;
}
/**
   Generate laplacian on loc_t
 */
dsp* mklaplacian_loc(loc_t* loc, real r0, real weight){
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
	loc_create_map(loc);
	map_t* map=loc->map;
	dsp* L2;
	L2=dspnew(loc->nloc, loc->nloc, loc->nloc*5);
	int ix, iy;
	spint* pp=L2->pp;
	spint* pi=L2->pi;
	real* px=L2->px;
	real* px0=L2->px;
	real cf=laplacian_coef(r0, weight, loc->dx);
#if USE_PARTIAL == 1
	real cfs[5];
	cfs[0]=cfs[1]=0;
	cfs[2]=laplacian_coef3(r0, weight, loc->dx);
	cfs[3]=laplacian_coef4(r0, weight, loc->dx);
	cfs[4]=cf;
#endif
	for(iy=0; iy<NY(map); iy++){
		for(ix=0; ix<NX(map); ix++){
			long iphi0;
			if((iphi0=loc_map_get(map, ix, iy))>0){
				*(pp++)=px-px0;
				long iphiL=loc_map_get(map, ix-1, iy);
				long iphiR=loc_map_get(map, ix+1, iy);
				long iphiD=loc_map_get(map, ix, iy-1);
				long iphiU=loc_map_get(map, ix, iy+1);

#if USE_PARTIAL == 1
				real* px1=px;
				real* p1, * p2;
				if(iphiD>0){
					*(pi++)=iphiD-1;
					*(px++)=0.25;
					p2=px-1;
				} else{
					p2=NULL;
				}
				if(iphiL>0){
					*(pi++)=iphiL-1;
					*(px++)=0.25;
					p1=px-1;
				} else{
					p1=NULL;
				}
				*(pi++)=iphi0-1;
				*(px++)=-1;
				if(iphiR>0){
					*(pi++)=iphiR-1;
					if(p1)
						*(px++)=0.25;
					else
						*(px++)=0.5;
				} else{
					if(p1)
						*p1=0.5;
					else
						warning("Point is isolated");
				}
				if(iphiU>0){
					*(pi++)=iphiU-1;
					if(p2)
						*(px++)=0.25;
					else
						*(px++)=0.5;
				} else{
					if(p2)
						*p2=0.5;
					else
						warning("Point is isolated");
				}
				real cfi=cfs[px-px1-1];
				for(real* px2=px1; px2<px; px2++){
					*px2*=cfi;
				}
#else
				int has_up, has_rt;
				if(iphiD>0&&iphiU>0){
					*(pi++)=iphiD-1;
					*(px++)=0.25*cf;
					has_up=1;
				} else{
					has_up=0;
				}
				if(iphiL>0&&iphiR>0){
					*(pi++)=iphiL-1;
					*(px++)=0.25*cf;
					has_rt=1;
				} else{
					has_rt=0;
				}
				if(has_rt||has_up){
					*(pi++)=iphi0-1;
					*(px++)=(-0.5*has_rt-0.5*has_up)*cf;
				}
				if(has_rt){
					*(pi++)=iphiR-1;
					*(px++)=0.25*cf;
				}
				if(has_up){
					*(pi++)=iphiU-1;
					*(px++)=0.25*cf;
				}

#endif
			}
		}
	}
	*pp=px-px0;
	if(px-px0>L2->nzmax){
		error("Over flow happened\n");
	}
	dspsetnzmax(L2, px-px0);
	dsp* L2r=dsptrans(L2);
	dspfree(L2);
	return L2r;
}
