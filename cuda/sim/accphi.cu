/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "accphi.h"
#include "cudata.h"
/**
   Ray tracing from map to loc with boundary check. Real input
   This is memory bound. So increasing # of points processed does not help.
*/
__global__ void map2loc_linear(Real* restrict out, const Real* restrict in,
	const int nx, const int ny, KARG_COMMON){
	int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<nloc; i+=step){
		Real x=loc[i][0]*dxi+dispx;
		Real y=loc[i][1]*dyi+dispy;
		int ix=Z(floor)(x);
		int iy=Z(floor)(y);
		x=x-ix; y=y-iy;
		if(ix>=0&&ix<nx-1&&iy>=0&&iy<ny-1){
			Real tmp=((+in[iy*nx+ix]*(1.f-x)+in[iy*nx+ix+1]*x)*(1.f-y)
				+(+in[(iy+1)*nx+ix]*(1.f-x)+in[(iy+1)*nx+ix+1]*x)*y);
			add_valid(out[i], alpha, tmp);
		}
	}
}

/**
   Ray tracing from map to loc with boundary check. Complex input
*/
__global__ void map2loc_linear(Real* restrict out, const Comp* restrict in,
	const int nx, const int ny, KARG_COMMON){
	int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<nloc; i+=step){
		Real x=loc[i][0]*dxi+dispx;
		Real y=loc[i][1]*dyi+dispy;
		int ix=Z(floor)(x);
		int iy=Z(floor)(y);
		x=x-ix; y=y-iy;
		if(ix>=0&&ix<nx-1&&iy>=0&&iy<ny-1){
			Real tmp=((+in[iy*nx+ix].x*(1.f-x)+in[iy*nx+ix+1].x*x)*(1.f-y)
				+(+in[(iy+1)*nx+ix].x*(1.f-x)+in[(iy+1)*nx+ix+1].x*x)*y);
			add_valid(out[i], alpha, tmp);
		}
	}
}
/*
  Ray tracing from map to loc without boundary check. Real input.
*/
__global__ void map2loc_linear_nocheck(Real* restrict out, const Real* restrict in,
	const int nx, const int ny, KARG_COMMON){
	int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<nloc; i+=step){
		Real x=loc[i][0]*dxi+dispx;
		Real y=loc[i][1]*dyi+dispy;
		int ix=Z(floor)(x);
		int iy=Z(floor)(y);
		x=x-ix; y=y-iy;
		out[i]+=alpha*((in[iy*nx+ix]*(1-x)+in[iy*nx+ix+1]*x)*(1-y)
			+(in[(iy+1)*nx+ix]*(1-x)+in[(iy+1)*nx+ix+1]*x)*y);
	}
}
/*
  Ray tracing from map to loc with wrapping. Real input.
*/
__global__ void map2loc_linear_wrap(Real* restrict out, const Real* restrict in,
	const int nx, const int ny, KARG_COMMON){
	int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<nloc; i+=step){
		Real x=loc[i][0]*dxi+dispx;
		Real y=loc[i][1]*dyi+dispy;
		int ix=Z(floor)(x);
		int iy=Z(floor)(y);
		x=x-ix; y=y-iy;
		while(ix<0) ix+=nx;
		while(iy<0) iy+=ny;
		while(ix>nx-1) ix-=nx;
		while(iy>ny-1) iy-=ny;
		int ix1=(ix==nx-1)?0:(ix+1);
		int iy1=(iy==ny-1)?0:(iy+1);
		out[i]+=alpha*((in[iy*nx+ix]*(1.f-x)+in[iy*nx+ix1]*x)*(1.f-y)
			+(in[(iy1)*nx+ix]*(1.f-x)+in[(iy1)*nx+ix1]*x)*y);
	}
}

/*This is memory bound. So increasing # of points processed does not help. */
__global__ void map2loc_cubic(Real* restrict out, const Real* restrict in,
	const int nx, const int ny, KARG_COMMON, const Real* cc){
	int step=blockDim.x*gridDim.x;
	for(int i=blockIdx.x*blockDim.x+threadIdx.x; i<nloc; i+=step){
		Real x=loc[i][0]*dxi+dispx;
		Real y=loc[i][1]*dyi+dispy;
		int ix=Z(floor)(x); x=x-ix;
		int iy=Z(floor)(y); y=y-iy;
		Real fx[4], fy;
		Real sum=0;
		if(ix<1||ix>nx-3||iy<1||iy>ny-3){
			continue;/*out of range. */
		}
		/*cc need to be in device memory for sm_13 to work.*/
		fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));
		fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);
		fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));
		fx[3]=x*x*(cc[3]+cc[4]*x);

		fy=(1.f-y)*(1.f-y)*(cc[3]+cc[4]*(1.f-y));
#pragma unroll
		for(int kx=-1; kx<3; kx++){
			sum+=fx[kx+1]*fy*in[(iy-1)*nx+kx+ix];
		}

		fy=cc[0]+y*y*(cc[1]+cc[2]*y);
#pragma unroll
		for(int kx=-1; kx<3; kx++){
			sum+=fx[kx+1]*fy*in[iy*nx+kx+ix];
		}

		fy=cc[0]+(1.f-y)*(1.f-y)*(cc[1]+cc[2]*(1.f-y));
#pragma unroll
		for(int kx=-1; kx<3; kx++){
			sum+=fx[kx+1]*fy*in[(iy+1)*nx+kx+ix];
		}

		fy=y*y*(cc[3]+cc[4]*y);
#pragma unroll
		for(int kx=-1; kx<3; kx++){
			sum+=fx[kx+1]*fy*in[(iy+2)*nx+kx+ix];
		}
		add_valid(out[i], alpha, sum);
	}
}

/**
   Ray tracing of atm.
*/
void atm2loc(Real* phiout, const culoc_t& loc, Real hs, Real hc, Real thetax, Real thetay,
	Real misregx, Real misregy, Real dt, int isim, Real atmalpha0, cudaStream_t stream){
	cumapcell& cuatm=cudata->atm;
	Real atmalpha=atmalpha0;
	if(Z(fabs)(atmalpha)<EPS) return;
	const int atm_dtrat=cuglobal->atm_dtrat;
	if(cuglobal->atmscale&&!atm_dtrat){
		atmalpha*=cuglobal->atmscale->p[isim];
	}
	int nps=atm_dtrat?cuglobal->atm_nps:cudata->atm.N();
	real wt=1;
	for(int jps=0; jps<nps; jps++){
		int ips;
		if(atm_dtrat){
			ips=atm_interp(&wt, jps, isim, atm_dtrat, cudata->atm.N(), cuglobal->atm_interp);
			/*ips=wrap_seq(isim/atm_dtrat+jps, cudata->atm.N());
			Real wt2=0;
			if(nps>1&&atm_interp){
				wt2=(real)(isim%atm_dtrat)/atm_dtrat;
				if(atm_interp==2){
					wt2=pow(sin(wt2*M_PI/2), 2);//smoother interp with sin^2 function
				}
			}*/
			//dbg("isim=%d, atmalpha=%g\n", isim, atmalpha);
		}else{
			ips=jps;
		}
		const Real dx=cuatm[ips].dx;
		const Real dy=cuatm[ips].dy;
		const Real ht=cuatm[ips].ht;
		const Real vx=cuatm[ips].vx;
		const Real vy=cuatm[ips].vy;
		const Real scale=1.f-(ht-hc)/(hs-hc);
		if(scale<0) continue;
		const Real dispx=(ht*thetax+scale*misregx-vx*dt*isim-cuatm[ips].ox)/dx;
		const Real dispy=(ht*thetay+scale*misregy-vy*dt*isim-cuatm[ips].oy)/dy;
		const int nloc=loc.Nloc();

#define COMM loc(),loc.Nloc(),scale/dx,scale/dy, dispx, dispy, atmalpha*wt
		if(cuglobal->atm_full){
			map2loc_linear_wrap<<<DIM(nloc, 256), 0, stream>>>
				(phiout, cuatm[ips](), cuatm[ips].nx, cuatm[ips].ny, COMM);
		} else{/*we are gauranteed. */
			//check boundary
			if(loc.xmin/dx*scale+dispx>=0&&loc.ymin/dy*scale+dispy>=0
				&&loc.xmax/dx*scale+dispx+1<cuatm[ips].nx&&loc.ymax/dy*scale+dispy+1<cuatm[ips].ny){
				map2loc_linear_nocheck<<<DIM(nloc, 256), 0, stream>>>
					(phiout, cuatm[ips](), cuatm[ips].nx, cuatm[ips].ny, COMM);
			} else{
				warning("Unexpected: need to check boundary. min=(%g, %g), max=(%g, %g), map: (%ld, %ld)\n",
					loc.xmin/dx+dispx, loc.ymin/dx+dispy,
					loc.xmax/dx+dispx+1, loc.ymax/dx+dispy+1, cuatm[ips].nx, cuatm[ips].ny);
				print_backtrace();
				map2loc_linear<<<DIM(nloc, 256), 0, stream>>>
					(phiout, cuatm[ips](), cuatm[ips].nx, cuatm[ips].ny, COMM);
			}
		}
#undef COMM
	}
}
void map2loc(const cumap_t& map, const culoc_t& loc, Real* phiout,
	Real alpha, Real dispx, Real dispy, Real scale, int wrap, cudaStream_t stream){
	if(scale<0) return;
	dispx=(dispx-map.ox)/map.dx;
	dispy=(dispy-map.oy)/map.dy;
	const int nloc=loc.Nloc();
	if(map.cubic_cc){//128 is a good number for cubic. 
		if(wrap){
			error("Not supported\n");
		} else{
			map2loc_cubic<<<DIM(nloc, 128), 0, stream>>>
				(phiout, map(), map.nx, map.ny, loc(), loc.Nloc(), scale/map.dx, scale/map.dy, dispx, dispy, alpha,
					map.cubic_cc());
		}
	} else{
		if(wrap){
			map2loc_linear_wrap<<<DIM(nloc, 256), 0, stream>>>
				(phiout, map(), map.nx, map.ny, loc(), loc.Nloc(), scale/map.dx, scale/map.dy, dispx, dispy, alpha);
		} else{
			map2loc_linear<<<DIM(nloc, 256), 0, stream>>>
				(phiout, map(), map.nx, map.ny, loc(), loc.Nloc(), scale/map.dx, scale/map.dy, dispx, dispy, alpha);
		}
	}
}

/**
   Ray tracing of dm. use a different loc for each dm. so that distortion can be
   properly accounted for. Use the other version if no distortion.
*/
void mapcell2loc(Real* phiout, const Array<culoc_t>& locondm, const cumapcell& cudm, 
	Real hs, Real hc, Real thetax, Real thetay, Real misregx, Real misregy, Real alpha, cudaStream_t stream){
	const Real theta=RSS(thetax, thetay);
	for(int idm=0; idm<cudm.N(); idm++){
		assert(cudm[idm].ny>1);//prevent accidentally pass in a vector
		const Real ht=cudm[idm].ht;
		const Real scale=1.-(ht-hc)/(hs-hc);
		map2loc(cudm[idm], locondm[idm], phiout, alpha*cos(theta*cudm[idm].dratio), ht*thetax+scale*misregx, ht*thetay+scale*misregy, scale, 0, stream);
	}
}
/**
   Ray tracing of dm.
*/
void mapcell2loc(Real* phiout, const culoc_t& locout, const cumapcell& cudm, 
	Real hs, Real hc, Real thetax, Real thetay, Real misregx, Real misregy, Real alpha, cudaStream_t stream){
	const Real theta=RSS(thetax, thetay);
	for(int idm=0; idm<cudm.N(); idm++){
		assert(cudm[idm].ny>1);//prevent accidentally pass in a vector
		const Real ht=cudm[idm].ht;
		const Real scale=1.-(ht-hc)/(hs-hc);
		map2loc(cudm[idm], locout, phiout, alpha*cos(theta*cudm[idm].dratio), ht*thetax+scale*misregx, ht*thetay+scale*misregy, scale, 0, stream);
	}/*idm */
}
/**
   Convert NGS mode vector to aperture grid for science directions.  */
void ngsmod2loc(curmat& opd, Real(*restrict loc)[2],
	const ngsmod_t* ngsmod, const real* mod,
	real thetax, real thetay,
	real alpha, cudaStream_t stream){
	if(ngsmod->nmod==2){
		curaddptt(opd, loc, 0, mod[0]*alpha, mod[1]*alpha, stream);
	} else{
		const Real ht=ngsmod->hdm;
		const Real scale=ngsmod->scale;

		Real focus=0, ps1=0, ps2=0, ps3=0, astigx=0, astigy=0;
		if(ngsmod->indfocus){
			focus+=mod[ngsmod->indfocus];
		}
		if(ngsmod->indps){
			if(!ngsmod->ahst_focus){
				focus+=mod[ngsmod->indps]*(1.f-scale);
			}
			ps1=mod[ngsmod->indps];
			ps2=mod[ngsmod->indps+1];
			ps3=mod[ngsmod->indps+2];
		}
		if(ngsmod->indastig){
			astigx=mod[ngsmod->indastig];
			astigy=mod[ngsmod->indastig+1];
		}

		add_ngsmod_do<<<DIM(opd.N(), 256), 0, stream>>>
			(opd(), loc, opd.N(),
				mod[0], mod[1], ps1, ps2, ps3, astigx, astigy, focus,
				thetax, thetay, scale, ht, alpha);
	}
}
