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

#include "cure.h"
/**
   Cumulative Reconstructor

   Implementation is based on M.Rosensteiner, J. Opt. Soc. Am. A28, 2132-2138 (2011)

*/

void cure1d(dmat** pphix,   /**<Output: opd*/
	const dmat* gx,/**<x gradients, dimension is nx*ny. Invalid subapertures set to NaN*/
	const dmat* gy,/**<y gradients, dimension is ny*ny.*/
	real dx      /**<Size of subaperture*/
){
	const long nx=gx->nx;
	const long ny=gy->ny;
	const long maxseg=(nx+1)/2;
	dmat* lx=dnew(nx+1, ny);
	lmat* flag=lnew(1+maxseg*2, ny);//nseg, iseg[0].start, iseg[0].end, ...
	for(long iy=0; iy<ny; iy++){
	//Sum along first dimension (x)
		int on=0, iseg=0;
		for(long ix=0; ix<nx; ix++){
			int invalid=isnan(P(gx, ix, iy));
			if(!on&&!invalid){ //beginning of a segment
				on=1;
				iseg=P(flag, 0, iy);
				P(flag, 0, iy)++;
				if(P(flag, 0, iy)>maxseg){
					error("nseg=%ld > maxseg=%ld\n", P(flag, 0, iy), maxseg);
				}
				P(flag, iseg*2+1, iy)=ix;//mark beginning
				P(lx, ix, iy)=0;//set value at beginning of segment to 0
			} else if(on&&invalid){ //end of a segment
				on=0;
				P(flag, iseg*2+2, iy)=ix;//mark end (exclusive)
			}//else: continuation
			//integrating. Invalid subapertures set lx to NAN
			P(lx, ix+1, iy)=P(lx, ix, iy)+dx*P(gx, ix, iy);
		}
		if(on){//end of last segment (exclusive)
			P(flag, iseg*2+2, iy)=nx;
		}

		//Remove piston
		long nseg=P(flag, 0, iy);
		for(iseg=0; iseg<nseg; iseg++){
			long ix0=P(flag, iseg*2+1, iy);
			long ix1=P(flag, iseg*2+2, iy)+1;
			real lxs=0;
			//remove mean of each x-chain
			for(long ix=ix0; ix<ix1; ix++){
				lxs+=P(lx, ix, iy);
			}
			lxs/=-(ix1-ix0);
			//connect x-chains along y direction
			for(long ix=ix0; ix<ix1; ix++){
				P(lx, ix, iy)+=lxs;
			}
		}
	}
	//writebin(flag, "flag");
	//writebin(lx, "lx_1");

	//Connecting chains
	dmat* ty=dnew(maxseg, ny);
	long nseglast=P(flag, 0, 0);
	real piston=0; long pcount=0;
	//Include #elemenets in first chain in pcount
	for(long iy=0; iy<1; iy++){
		for(long iseg=0; iseg<nseglast; iseg++){
			pcount+=P(flag, iseg*2+2, iy)-P(flag, iseg*2+1, iy)+1;
		}
	}
	//Start alignment from second chain to previous chain
	for(long iy=1; iy<ny; iy++){
		long nseg=P(flag, 0, iy);
		for(long iseg=0; iseg<nseg; iseg++){
			const long ig0=P(flag, iseg*2+1, iy);
			const long ig1=P(flag, iseg*2+2, iy);
			real gym=0;
			long gyc=0;
			//compute mean gy along common boundary
			for(long ig=ig0; ig<ig1; ig++){
				real gyi=P(gy, ig, iy)+P(gy, ig, iy-1);
				if(!isnan(gyi)){//both valid
					gym+=gyi;
					gyc++;
				}
			}
			real lxm=0;
			long lxc=0;
			//compute x-chain mean difference along common boundary
			for(long ix=ig0; ix<=ig1; ix++){
				real lxd=P(lx, ix, iy)-P(lx, ix, iy-1);
				if(!isnan(lxd)){
					lxm+=lxd;
					lxc++;
				}
			}
			//skip if no common boundary
			if(!lxc||!gyc||isnan(gym)||isnan(lxm)){
				continue;
			}
			real tylast=0;
			//Handle mismatch in #segment
			if(nseglast!=nseg){//combine
				if(nseglast==1){//break up
					tylast=P(ty, 0, iy-1);
				} else{//combine
					tylast=0;
					for(long jseg=0; jseg<nseglast; jseg++){
						tylast+=P(ty, jseg, iy-1);
					}
					tylast/=nseglast;
				}
			} else{
				tylast=P(ty, iseg, iy-1);
			}
			P(ty, iseg, iy)=tylast+0.5*dx*(gym/gyc)-lxm/lxc;

			piston+=(ig1-ig0+1)*P(ty, iseg, iy);
			pcount+=(ig1-ig0+1);
		}
		nseglast=nseg;
	}
	piston/=pcount;
	//Adjust x chains. Separate from previous loop as lxd needs previous lx.
	for(long iy=0; iy<ny; iy++){
		long nseg=P(flag, 0, iy);
		for(long iseg=0; iseg<nseg; iseg++){
			long ix0=P(flag, iseg*2+1, iy);
			long ix1=P(flag, iseg*2+2, iy)+1;

			//connect x-chains along y direction and remove piston
			for(long ix=ix0; ix<ix1; ix++){
				P(lx, ix, iy)+=P(ty, iseg, iy)-piston;
			}
		}
	}


	//writebin(ty, "ty_1");
	// writebin(lx, "lx_2");
	//Extrapolate from x chain defined on edge to corner points
	if(*pphix){
		if((*pphix)->nx<nx+1||(*pphix)->ny<ny+1){
			warning("Output array is too small, resize\n");
			dresize(*pphix, nx+1, ny+1);
		}
	} else{
		*pphix=dnew(nx+1, ny+1);//x chains: integration of gx along x
	}
	dmat* phix=*pphix;

	//Interpolate inner points first
	for(long iy=0; iy<=ny; iy++){
		for(long ix=0; ix<=nx; ix++){
			real val=0;
			long count=0;
			if(ix>0&&iy>0&&!isnan(P(gy, ix-1, iy-1))){//BL
				val+=P(lx, ix, iy-1)+P(gy, ix-1, iy-1)*dx*0.5;
				count++;
			}
			if(iy>0&&ix<nx&&!isnan(P(gy, ix, iy-1))){//BR
				val+=P(lx, ix, iy-1)+P(gy, ix, iy-1)*dx*0.5;
				count++;
			}
			if(ix>0&&iy<ny&&!isnan(P(gy, ix-1, iy))){//TL
				val+=P(lx, ix, iy)-P(gy, ix-1, iy)*dx*0.5;
				count++;
			}
			if(ix<nx&&iy<ny&&!isnan(P(gy, ix, iy))){//TR
				val+=P(lx, ix, iy)-P(gy, ix, iy)*dx*0.5;
				count++;
			}

			if(count){
				P(phix, ix, iy)+=(val/count);
			}
		}
	}

	dfree(lx);
	dfree(ty);
	lfree(flag);
}
/**
   Implements 2d cure by calling cure1d twice.
*/
void cure(dmat** phi, const dmat* gx, const dmat* gy, real dx){
	dmat* gxt=dtrans(gx);
	dmat* gyt=dtrans(gy);
	dmat* phi2t=0;
	cure1d(phi, gx, gy, dx);
	cure1d(&phi2t, gyt, gxt, dx);
	dmat* phi2=dtrans(phi2t);
	dadd(phi, 0.5, phi2, 0.5);
	dfree(phi2);
	dfree(phi2t);
	dfree(gxt);
	dfree(gyt);
}
/**
   An alternative interface. grad is nsa*2 defined on saloc, instead of on square grid.
*/
void cure_loc(dmat** phix, const dmat* grad, const loc_t* saloc){
	if(grad->nx*grad->ny!=saloc->nloc*2){
		warning("Incorrect dimensions: grad %ldx%ld, saloc %ldx2",
			grad->nx, grad->ny, saloc->nloc);
		return;
	}
	loc_create_map(saloc);
	int nxm=saloc->map->nx;
	int npad=saloc->npad;
	int nx=saloc->map->nx-npad*2;
	int ny=saloc->map->ny-npad*2;
	dmat* gx=dnew(nx, ny);
	dmat* gy=dnew(nx, ny);
	long nsa=grad->nx/2;
	for(int iy=0; iy<ny; iy++){
		for(int ix=0; ix<nx; ix++){
			long ii=saloc->map->p[(ix+npad)+(iy+npad)*nxm];
			if(ii>0){
				gx->p[ix+iy*nx]=P(grad, ii-1, 0);
				gy->p[ix+iy*nx]=P(grad, ii-1+nsa, 0);
			} else{
				gx->p[ix+iy*nx]=NAN;
				gy->p[ix+iy*nx]=NAN;
			}
		}
	}
	cure(phix, gx, gy, saloc->dx);
	dfree(gx);
	dfree(gy);
}
/**
   An alternative interface. grad is nsa*2 defined on saloc, instead of on
   square grid. Return phix also on saloc grid.
*/
void cure_loc2(dmat** phix, const dmat* grad, const loc_t* saloc){
	dmat* phi2d=0;
	cure_loc(&phi2d, grad, saloc);
	if(!*phix){
		*phix=dnew(saloc->nloc, 1);
	} else if((*phix)->nx!=saloc->nloc){
		error("Mismatch: %ld vs %ld\n", (*phix)->nx, saloc->nloc);
	}
	dmat* phi=*phix;
	int npad=saloc->npad;
	int nxm=saloc->map->nx;
	for(int iy=0; iy<phi2d->ny-1; iy++){
		for(int ix=0; ix<phi2d->nx-1; ix++){
			long ii=saloc->map->p[(ix+npad)+(iy+npad)*nxm];
			if(ii>0){
				P(phi, ii-1)=P(phi2d, ix, iy);
			}
		}
	}
	dfree(phi2d);
}
