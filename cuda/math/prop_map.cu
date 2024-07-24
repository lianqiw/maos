/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "prop_map.h"

/*
  One kernel that handles multiple layers/directions, linear or cubic to enable fully parallelization across the full GPU.

  Forward: Propagate from XLOC to WFS (Parallel across each WFS).
  Backward (Transpose): Propagate from WFS to XLOC (Parallel across each Layer)

  let XLOC be ps (phase screen)
  let WFS be dir (direction)

  alpha2: additional scaling defined on dimension dir.

  Do not separate the function branches because each layer/wfs combination may use different branches.

  The output grid aligns with the input grid along x/y. There cannot be rotation or higher order distortion effects.

*/
__global__ void
map2map_do(map2map_t* data, Real* const* pdirs, Real* const* ppss, int ndir, int nps, Real alpha1, const Real* alpha2, char trans){
	/*Using shared memory to reduce register spill */
	__shared__ map2map_shared_t shared;
	int& stepx=shared.stepx;
	int& stepy=shared.stepy;
	Real& dispx=shared.dispx;
	Real& dispy=shared.dispy;
	Real& xratio=shared.xratio;
	Real& yratio=shared.yratio;
	int& ndirx=shared.ndirx;
	int& ndiry=shared.ndiry;
	int& nx=shared.nx;
	int& ny=shared.ny;
	int& npsx=shared.npsx;
	int& nn=shared.nn;
	Real*& pps=shared.pps;
	Real*& pdir=shared.pdir;

	int& ix0=shared.ix0[threadIdx.x];
	int& iy0=shared.iy0[threadIdx.y];
	//only update shared variables by a single thread
	if(threadIdx.x==0&&threadIdx.y==0){
		if(ndir==0){//layer to layer. caching mechanism
			if(gridDim.z!=nps) return;
			nn=1;
		} else if(trans=='t'){
			if(gridDim.z!=nps) return;
			nn=ndir;
		} else{
			if(gridDim.z!=ndir) return;
			nn=nps;
		}
	}
	if(threadIdx.y==0){
		ix0=blockIdx.x*blockDim.x+threadIdx.x;
	}
	if(threadIdx.x==0){
		iy0=blockIdx.y*blockDim.y+threadIdx.y;
	}
	__syncthreads();//necessary here because otherwise different warps may modify the shared data.
	for(int ii=0; ii<nn; ii++){//number of screens to trace
		map2map_t* datai;
		int ips, idir;
		if(ndir==0){//plane to plane. no direction
			ips=idir=blockIdx.z;
			datai=data+blockIdx.z;
		} else{
			if(trans=='t'){
				ips=blockIdx.z;
				idir=ii;
			} else{
				idir=blockIdx.z;
				ips=ii;
			}
			datai=data+idir+ndir*ips;
		}
		const int match=Z(fabs)(datai->xratio-1.f)<EPS&&Z(fabs)(datai->yratio-1.f)<EPS;
		if(!datai->cc&&trans=='t'&&match){
			datai=datai->reverse;
		}
		if(datai->nx==0) continue;//skip empty wfs
		const Real alpha=alpha2?(alpha2[idir]*alpha1):alpha1;
		const Real* cc=datai->cc;
		__syncthreads();//necessary here because different warps may be doing different ii.
		if(threadIdx.x==0&&threadIdx.y==0){//only update shared variables by a single thread
			stepx=blockDim.x*gridDim.x;
			stepy=blockDim.y*gridDim.y;
			dispx=datai->dispx;
			dispy=datai->dispy;
			xratio=datai->xratio;
			yratio=datai->yratio;
			ndirx=datai->nxdir;
			ndiry=datai->nydir;
			nx=datai->nx;
			ny=datai->ny;
			npsx=datai->nxps;
			pdir=pdirs[idir]+datai->offdir;
			pps=ppss[ips]+datai->offps;
		}
		__syncthreads();
		if(!cc&&match){
			/* Matched bilinear propagation.  Reverse propagation is converted
			   into forward ray tracing. But atomic operation is necessary
			   because different thead blocks may be handling different
			   directions that have different offset and therefore have
			   concurrent memory access*/

			const Real dispx1=1.-dispx;
			const Real dispy1=1.-dispy;

			if(trans=='t'){
				for(int iy=iy0; iy<ny; iy+=stepy){
					for(int ix=ix0; ix<nx; ix+=stepx){
						atomicAdd(&pps[ix+iy*npsx],
							alpha*(+(pdir[ix+iy*ndirx]*dispx1
								+pdir[ix+1+iy*ndirx]*dispx)*dispy1
								+(pdir[ix+(iy+1)*ndirx]*dispx1
									+pdir[ix+1+(iy+1)*ndirx]*dispx)*dispy));
					}
				}
			} else{
				for(int iy=iy0; iy<ny; iy+=stepy){
					for(int ix=ix0; ix<nx; ix+=stepx){
						pdir[ix+iy*ndirx]+=
							alpha*(+(pps[ix+iy*npsx]*dispx1
								+pps[ix+1+iy*npsx]*dispx)*dispy1
								+(pps[ix+(iy+1)*npsx]*dispx1
									+pps[ix+1+(iy+1)*npsx]*dispx)*dispy);
					}
				}
			}
		} else if(cc){//cubic
			if(trans=='t'){
				if(Z(fabs)(xratio-0.5f)<EPS&&Z(fabs)(yratio-.5f)<EPS){
					//do without atomic operations.
					const int nxin=npsx;//ceil(nx*xratio)+1;
					const int nyin=datai->nyps;//ceil(ny*xratio)+1;
					const int xmaxdir=ndirx-datai->offdirx;
					const int ymaxdir=ndiry-datai->offdiry;
					const int xmindir=-datai->offdirx-1;
					const int ymindir=-datai->offdiry-1;
					int offx=0, offy=0;
					Real tcx=dispx;
					Real tcy=dispy;
					if(tcx>=0.5){
						offx=-1;
						tcx-=0.5f;
					}
					if(tcy>=0.5){
						offy=-1;
						tcy-=0.5f;
					}

					const Real tcx2=tcx+0.5f;
					const Real tcy2=tcy+0.5f;
					//Each thread has the same coefficients, so we use
					// shared memory to store them to avoid spill.
					Real* const& fx=(Real*)shared.fx;
					Real* const& fy=(Real*)shared.fy;
					if(threadIdx.x==0&&threadIdx.y==0){
						fy[0]=(cc[3]+cc[4]*tcy)*tcy*tcy;
						fy[1]=(cc[3]+cc[4]*tcy2)*tcy2*tcy2;
						fy[2]=(cc[1]+cc[2]*(1-tcy))*(1-tcy)*(1-tcy)+cc[0];
						fy[3]=(cc[1]+cc[2]*(1-tcy2))*(1-tcy2)*(1-tcy2)+cc[0];
						fy[4]=(cc[1]+cc[2]*tcy)*tcy*tcy+cc[0];
						fy[5]=(cc[1]+cc[2]*tcy2)*tcy2*tcy2+cc[0];
						fy[6]=(cc[3]+cc[4]*(1-tcy))*(1-tcy)*(1-tcy);
						fy[7]=(cc[3]+cc[4]*(1-tcy2))*(1-tcy2)*(1-tcy2);

						fx[0]=(cc[3]+cc[4]*tcx)*tcx*tcx;
						fx[1]=(cc[3]+cc[4]*tcx2)*tcx2*tcx2;
						fx[2]=(cc[1]+cc[2]*(1-tcx))*(1-tcx)*(1-tcx)+cc[0];
						fx[3]=(cc[1]+cc[2]*(1-tcx2))*(1-tcx2)*(1-tcx2)+cc[0];
						fx[4]=(cc[1]+cc[2]*tcx)*tcx*tcx+cc[0];
						fx[5]=(cc[1]+cc[2]*tcx2)*tcx2*tcx2+cc[0];
						fx[6]=(cc[3]+cc[4]*(1-tcx))*(1-tcx)*(1-tcx);
						fx[7]=(cc[3]+cc[4]*(1-tcx2))*(1-tcx2)*(1-tcx2);
					}
					__syncthreads();

					for(int my=iy0; my<nyin; my+=stepy){
						const int ycent=2*my+offy;
						for(int mx=ix0; mx<nxin; mx+=stepx){
							Real sum=0;
							const int xcent=2*mx+offx;
#pragma unroll
							for(int ky=-4; ky<4; ky++){
								const int ky2=ky+ycent;
								if(ky2>ymindir&&ky2<ymaxdir){
#pragma unroll
									for(int kx=-4; kx<4; kx++){
										const int kx2=kx+xcent;
										Real wt;
										if(kx2>xmindir&&kx2<xmaxdir&&(wt=fy[ky+4]*fx[kx+4])>EPS){
											sum+=wt*pdir[(kx2)+(ky2)*ndirx];
										}
									}
								}
							}
							//Need atomic because different layers have different offset.
							if(nn>1){
								atomicAdd(&pps[mx+my*npsx], sum*alpha);
							} else{
								pps[mx+my*npsx]+=sum*alpha;
							}
						}
					}
				} else{
					const int xmaxps=npsx-datai->offpsx;
					const int ymaxps=datai->nyps-datai->offpsy;
					const int xminps=-datai->offpsx;
					const int yminps=-datai->offpsy;
					Real fy[4]; Real fx[4];
					for(int my=iy0; my<ny; my+=stepy){
						Real jy;
						const Real y=Z(modf)(dispy+my*yratio, &jy);
						const int iy=(int)jy;
						fy[0]=(1.f-y)*(1.f-y)*(cc[3]+cc[4]*(1.f-y));
						fy[1]=cc[0]+y*y*(cc[1]+cc[2]*y);
						fy[2]=cc[0]+(1.f-y)*(1.f-y)*(cc[1]+cc[2]*(1.f-y));
						fy[3]=y*y*(cc[3]+cc[4]*y);
						for(int mx=ix0; mx<nx; mx+=stepx){
							Real jx;
							const Real x=Z(modf)(dispx+mx*xratio, &jx);
							const int ix=(int)jx;
							const Real value=pdir[mx+my*ndirx]*alpha;
							//cc need to be in device memory for sm_13 to work.
							//	if(threadIdx.x==0 && threadIdx.y==0){
							fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));
							fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);
							fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));
							fx[3]=x*x*(cc[3]+cc[4]*x);
							//	}
							//__syncthreads();
							const int ky0=(yminps-iy)>-1?(yminps-iy):-1;
							const int ky1=(ymaxps-iy)<3?(ymaxps-iy):3;
							for(int ky=ky0; ky<ky1; ky++){
								int kx0=(xminps-ix)>-1?(xminps-ix):-1;
								int kx1=(xmaxps-ix)<3?(xmaxps-ix):3;
								for(int kx=kx0; kx<kx1; kx++){
									Real wt;
									if((wt=fx[kx+1]*fy[ky+1])>EPS){
										atomicAdd(&pps[(iy+ky)*npsx+(kx+ix)], wt*value);
									}
								}
							}
						}
					}
				}
			} else{//CC, non trans
				const int xmaxps=npsx-datai->offpsx;
				const int ymaxps=datai->nyps-datai->offpsy;
				const int xminps=-datai->offpsx;
				const int yminps=-datai->offpsy;
				Real fy[4]; Real fx[4];
				for(int my=iy0; my<ny; my+=stepy){
					Real jy;
					const Real y=Z(modf)(dispy+my*yratio, &jy);
					const int iy=(int)jy;
					fy[0]=(1.f-y)*(1.f-y)*(cc[3]+cc[4]*(1.f-y));
					fy[1]=cc[0]+y*y*(cc[1]+cc[2]*y);
					fy[2]=cc[0]+(1.f-y)*(1.f-y)*(cc[1]+cc[2]*(1.f-y));
					fy[3]=y*y*(cc[3]+cc[4]*y);
					for(int mx=ix0; mx<nx; mx+=stepx){
						Real jx;
						const Real x=Z(modf)(dispx+mx*xratio, &jx);
						const int ix=(int)jx;
						Real sum=0;
						fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));
						fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);
						fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));
						fx[3]=x*x*(cc[3]+cc[4]*x);
						const int ky0=(yminps-iy)>-1?(yminps-iy):-1;
						const int ky1=(ymaxps-iy)<3?(ymaxps-iy):3;
						for(int ky=ky0; ky<ky1; ky++){
							int kx0=(xminps-ix)>-1?(xminps-ix):-1;
							int kx1=(xmaxps-ix)<3?(xmaxps-ix):3;
							for(int kx=kx0; kx<kx1; kx++){
								Real wt;
								if((wt=fx[kx+1]*fy[ky+1])>EPS)
									sum+=wt*pps[(iy+ky)*npsx+(kx+ix)];
							}
						}
						pdir[mx+my*ndirx]+=sum*alpha;
					}
				}
			}/*else trans*/
		} else if(trans=='t'){//linear, transpose
			for(int iy=iy0; iy<ny; iy+=stepy){
				Real jy;
				const Real fracy=Z(modf)(dispy+iy*yratio, &jy);
				const int ky=(int)jy;
				for(int ix=ix0; ix<nx; ix+=stepx){
					Real jx;
					const Real fracx=Z(modf)(dispx+ix*xratio, &jx);
					const int kx=(int)jx;
					const Real temp=pdir[ix+iy*ndirx]*alpha;
					Real wt;
					if((wt=(1.f-fracx)*(1.f-fracy))>EPS)
						atomicAdd(&pps[kx+ky*npsx], temp*wt);
					if((wt=fracx*(1.f-fracy))>EPS)
						atomicAdd(&pps[kx+1+ky*npsx], temp*wt);
					if((wt=(1.f-fracx)*fracy)>EPS)
						atomicAdd(&pps[kx+(ky+1)*npsx], temp*wt);
					if((wt=fracx*fracy)>EPS)
						atomicAdd(&pps[kx+1+(ky+1)*npsx], temp*wt);
				}
			}
		} else{//linear, regular
			for(int iy=iy0; iy<ny; iy+=stepy){
				Real jy;
				const Real fracy=Z(modf)(dispy+iy*yratio, &jy);
				const int ky=(int)jy;

				for(int ix=ix0; ix<nx; ix+=stepx){
					Real jx;
					const Real fracx=Z(modf)(dispx+ix*xratio, &jx);
					const int kx=(int)jx;
					Real tmp=0;
					Real wt;
					if((wt=1.f-fracy)>EPS){
						tmp+=(pps[kx+ky*npsx]*(1.f-fracx)+
							pps[kx+1+ky*npsx]*fracx)*wt;
					}
					if((wt=fracy)>EPS){
						tmp+=(pps[kx+(ky+1)*npsx]*(1.f-fracx)+
							pps[kx+1+(ky+1)*npsx]*fracx)*wt;
					}
					pdir[ix+iy*ndirx]+=alpha*tmp;
				}
			}
		}
	}/*for*/
}
/*
  Prepare data and copy to GPU so that one kernel (map2map_do) can handle multiple independent ray tracing without range checking.
  Forward: ps->dir (xloc -> wfs or dm->floc)
  Backward: dir->ps (wfs -> xloc or floc->dm)
*/

void map2map_prep(map2map_t* res, const cugrid_t& g_dir, const cugrid_t& g_ps,
	Real dispx, Real dispy, const curmat& cc){
	assert(g_ps.ny!=1);
	const Real uxi=1.f/g_ps.dx;
	const Real uyi=1.f/g_ps.dy;
	Real xratio=g_dir.dx*uxi;
	Real yratio=g_dir.dy*uyi;
	/*remove round off errors that causes trouble in nx, ny.*/
	if(Z(fabs)(xratio-1.f)<EPS&&Z(fabs)(yratio-1.f)<EPS){
		xratio=yratio=1.f;
		if(!cc&&!res->isreverse){
			res->reverse=new map2map_t;
			res->reverse->isreverse=1;
			map2map_prep(res->reverse, g_ps, g_dir, -dispx, -dispy, cc);
		}
	} else if(Z(fabs)(xratio-0.5f)<EPS&&Z(fabs)(yratio-0.5f)<EPS){
		xratio=yratio=0.5f;
	}
	const int nxdir=g_dir.nx;
	const int nydir=g_dir.ny;
	const int nxps=g_ps.nx;
	const int nyps=g_ps.ny;
	const Real xratio1=1.f/xratio;
	const Real yratio1=1.f/yratio;
	/*offset of origin in input grid spacing. */
	dispx=(dispx-g_ps.ox+g_dir.ox)*uxi;
	dispy=(dispy-g_ps.oy+g_dir.oy)*uyi;
	int offx1=0, offy1=0;/*for output. fine sampling. */
	/*if output is bigger than input, limit output grid. */
	if(dispx<0){
		offx1=(int)ceilf(-dispx*xratio1);
		dispx+=offx1*xratio;
	}
	if(dispy<0){
		offy1=(int)ceilf(-dispy*yratio1);
		dispy+=offy1*yratio;
	}
	int nd=cc?-1:1;//cubic allows extra points to be interpolated
	/*convert offset into input grid coordinate. -EPS to avoid laying on the last point. */
	int nx=(int)Z(floor)((nxps-nd-dispx-EPS)*xratio1+1);
	int ny=(int)Z(floor)((nyps-nd-dispy-EPS)*yratio1+1);
	//Sanity check.
	while((nx-1)*xratio+dispx+nd>=nxps){
		nx--;
	}
	while((nx)*xratio+dispx+nd<nxps){
		nx++;
	}
	while((ny-1)*yratio+dispy+nd>=nyps){
		ny--;
	}
	while((ny)*yratio+dispy+nd<nyps){
		ny++;
	}
	if(nx>nxdir-offx1) nx=nxdir-offx1;
	if(ny>nydir-offy1) ny=nydir-offy1;
	int offx2=(int)Z(floor)(dispx); 
	int offy2=(int)Z(floor)(dispy);
	if(cc){
		if(offx2>0) offx2--;
		if(offy2>0) offy2--;
	}
	dispx-=offx2;/*for input. coarse sampling. */
	dispy-=offy2;
	if(res->isreverse){//We un-reverse the naming
		res->offpsx=offx1;
		res->offpsy=offy1;
		res->offdirx=offx2;
		res->offdiry=offy2;
		res->offps=offy1*nxdir+offx1;
		res->offdir=offy2*nxps+offx2;
		res->nxps=nxdir;
		res->nyps=nydir;
		res->nxdir=nxps;
		res->nydir=nyps;
	} else{
		res->offdirx=offx1;
		res->offdiry=offy1;
		res->offpsx=offx2;
		res->offpsy=offy2;
		res->offdir=offy1*nxdir+offx1;
		res->offps=offy2*nxps+offx2;
		res->nxdir=nxdir;
		res->nydir=nydir;
		res->nxps=nxps;
		res->nyps=nyps;
	}
	res->dispx=dispx;
	res->dispy=dispy;
	res->xratio=xratio;
	res->yratio=yratio;
	res->nx=nx;
	res->ny=ny;
	res->cc=cc?cc():NULL;
	{
		Real nxps2, nyps2;
		Real fxps2=modff(offx2+dispx+(nx-1)*xratio, &nxps2);
		Real fyps2=modff(offy2+dispy+(ny-1)*yratio, &nyps2);
		if(offx1+nx>nxdir||offy1+ny>nydir||nxps2>=nxps-1||nyps2>=nyps-1){
			info("input: offset=%d, %d, to %g+%g, %g+%g, size %d, %d\n", offx2, offy2, nxps2, fxps2, nyps2, fyps2, nxps, nyps);
			info("output: offset=%d, %d, to %d, %d, size %d, %d\n", offx1, offy1, offx1+nx, offy1+ny, nxdir, nydir);
			error("overflow\n");
		}
	}
}

/**
   Initialize for ps to dir ray tracing, e.g., turbulence to WFS
*/
void map2map::init_l2d(const cugrid_t& out, const dir_t* dir, int _ndir, //output.
	const cugridcell& in,//input.
	Real delay){//directions and star height.
	if(hdata){
		deinit();
	}
	nlayer=in.N();
	ndir=_ndir;
	map2map_t* hdata_cpu=new map2map_t[nlayer*ndir];
	DO(cudaMalloc(&hdata, sizeof(map2map_t)*nlayer*ndir));
	for(int ilayer=0; ilayer<nlayer; ilayer++){
		const Real ht=in[ilayer].ht;
		for(int idir=0; idir<ndir; idir++){
			if(!dir[idir].skip){
				const Real scale=1.f-ht/dir[idir].hs;
				const Real dispx=dir[idir].thetax*ht+in[ilayer].vx*dir[idir].delay+dir[idir].misregx*scale;
				const Real dispy=dir[idir].thetay*ht+in[ilayer].vy*dir[idir].delay+dir[idir].misregy*scale;
				cugrid_t outscale=out.Scale(scale);
				map2map_prep(hdata_cpu+idir+ilayer*ndir, outscale, in[ilayer],
					dispx, dispy, in[ilayer].cubic_cc);
			}
			hdata_cpu[idir+ilayer*ndir].togpu(hdata+idir+ilayer*ndir);
		}
	}
	delete[] hdata_cpu;
}
/**
   Initialize for ps to ps ray tracing, e.g., for caching
*/
void map2map::init_l2l(const cugridcell& out, const cugridcell& in){//input. layers.
	if(nlayer) error("Already initialized\n");
	nlayer=in.N();
	ndir=0;//this is laye to layer.
	map2map_t* hdata_cpu=new map2map_t[nlayer];
	DO(cudaMalloc(&hdata, sizeof(map2map_t)*nlayer));
	for(int ilayer=0; ilayer<nlayer; ilayer++){
		if(Z(fabs)(out[ilayer].ht-in[ilayer].ht)>EPS){
			error("Layer height mis-match.\n");
		}
		map2map_prep(hdata_cpu+ilayer, out[ilayer], in[ilayer],
			0, 0, in[ilayer].cubic_cc);
		hdata_cpu[ilayer].togpu(hdata+ilayer);
	}
	delete[] hdata_cpu;
}
/**
   Rotate each OPD array by theta CCW (dir=-1) or CW (dir=1) around point (cx, cy)
 */
__global__ static void
map_rot_do(Real* const* outs, const Real* const* ins, Real cx, Real cy, long nx, long ny, const Real* wfsrot, int dir){
	const Real* in=ins[blockIdx.z];
	Real* out=outs[blockIdx.z];
	if(!in||!out) return;
	long nx1=nx-1;
	long ny1=ny-1;
	Real cs=wfsrot[blockIdx.z*2];
	Real ss=wfsrot[blockIdx.z*2+1]*dir;

	for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
		for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
			//Rotated index in the opposite direction
			Real ixnew=(ix-cx)*cs+(iy-cy)*ss+cx;
			Real iynew=-(ix-cx)*ss+(iy-cy)*cs+cy;
			int ix2=(int)floor(ixnew);
			int iy2=(int)floor(iynew);
			if(ix2>=0&&iy2>=0&&ix2<nx1&&iy2<ny1){
				ixnew-=ix2;
				iynew-=iy2;
				out[ix+iy*nx]+=
					+(in[ix2+iy2*nx]*(1-ixnew)+in[(1+ix2)+iy2*nx]*ixnew)*(1-iynew)
					+(in[ix2+(1+iy2)*nx]*(1-ixnew)+in[(1+ix2)+(1+iy2)*nx]*ixnew)*(iynew);
			}
		}
	}
}
/**
   Rotate each OPD array by theta CCW around point (cx, cy)
 */
void map_rot(curcell& out, const curcell& in, const curmat& wfsrot, int dir, stream_t& stream){
	long nx=0;
	long ny=0;
	long nb=in.N();
	for(long ib=0; ib<nb; ib++){
		if(in[ib].N()){
			if(!nx){
				nx=in[ib].Nx();
				ny=in[ib].Ny();
			} else if(nx!=in[ib].Nx()||ny!=in[ib].Ny()){
				error("Different cell has different dimensions (%ldx%ld) vs (%ldx%ld)\n",
					nx, ny, in[ib].Nx(), in[ib].Ny());
			}
		}
	}
	long cx=nx/2;//Fixed on 2024-07-18 
	long cy=ny/2;
	if(out.M()){
		out.M().Zero(stream);
	} else{
		out.Zero(stream);
	}
	map_rot_do<<<DIM3(nx, ny, 16, nb), 0, stream>>>
		(out.pm, in.pm, cx, cy, nx, ny, wfsrot, dir);

		/*static int count=-1; count++;
		if(count<10){
		cuwrite(in, stream, "in_%d_%d", count, dir);
		cuwrite(out, stream, "out_%d_%d", count, dir);
		}*/
}


