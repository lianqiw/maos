/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "utils.h"
#include "accphi.h"
#include "prop_wrap.h"

/*
  One kernel that handles multiple layers/directions, linear or cubic to enable fully parallelization across the full GPU.

  Forward: Propagate from XLOC to WFS (Parallel across each WFS).  
  Backward (Transpose): Propagate from WFS to XLOC (Parallel across each Layer)
  
  let XLOC be ps (phase screen)
  let WFS be dir (direction)

  alpha2: additional scaling defined on dimension dir.

  Do not separate the function branches because each layer/wfs combination may use different branches.
*/
__global__ void 
gpu_map2map_do(PROP_WRAP_T *data, Real *const*pdirs, Real *const*ppss, int ndir, int nps, Real alpha1, const Real *alpha2, char trans){
    /*Using shared memory to reduce register spill */
    __shared__ gpu_map2map_shared_t shared;
    int &stepx=shared.stepx;
    int &stepy=shared.stepy;
    Real &dispx=shared.dispx;
    Real &dispy=shared.dispy;
    Real &xratio=shared.xratio;
    Real &yratio=shared.yratio;
    int &ndirx=shared.ndirx;
    int &ndiry=shared.ndiry;
    int &nx=shared.nx;
    int &ny=shared.ny;
    int &npsx=shared.npsx;
    int &nn=shared.nn;
    Real *&pps=shared.pps;
    Real *&pdir=shared.pdir;
    
    int &ix0=shared.ix0[threadIdx.x];
    int &iy0=shared.iy0[threadIdx.y];

    if(threadIdx.x==0 && threadIdx.y==0){
	if(ndir==0){//layer to layer. caching mechanism
	    if(gridDim.z!=nps) return;
	    nn=1;
	}else if(trans=='t'){
	    if(gridDim.z!=nps) return;
	    nn=ndir;
	}else{
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
    __syncthreads();//necessary here because otherwise different wraps may modify the shared data.
    for(int ii=0; ii<nn; ii++){
	PROP_WRAP_T *datai;
	int ips, idir;
	if(ndir==0){//plane to plane. no direction
	    ips=idir=blockIdx.z;
	    datai=data+blockIdx.z;
	}else{
	    if(trans=='t'){
		ips=blockIdx.z;
		idir=ii;
	    }else{
		idir=blockIdx.z;
		ips=ii;
	    }
	    datai=data+idir+ndir*ips;
	}
	const int match=Z(fabs)(datai->xratio-1.f)<EPS && Z(fabs)(datai->yratio-1.f)<EPS;
	if(!datai->cc && trans=='t' && match){
	    datai=datai->reverse;
	}
	if(datai->nx==0) continue;//skip empty wfs
	const Real alpha=alpha2?(alpha2[idir]*alpha1):alpha1;
	const Real *cc=datai->cc;
	__syncthreads();//necessary here because different warps may be doing different ii.
	if(threadIdx.x==0 && threadIdx.y==0){
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
	if(!cc && match){
	    //Matched bilinear propagation. always forward prop. 
	    /*During reverse operation, for different idir, the offo is
	      different causing same thread to handle different memory in pdir
	      for different idir. This causes problem with pdir synchronization/atomic operation*/
	    const Real dispx1=1.-dispx;
	    const Real dispy1=1.-dispy;

	    if(trans=='t'){
		for(int iy=iy0; iy<ny; iy+=stepy){
		    for(int ix=ix0; ix<nx; ix+=stepx){
			atomicAdd(&pps[ix+iy*npsx],
				  alpha*(+(pdir[ix+    iy*ndirx]*dispx1
					   +pdir[ix+1+    iy*ndirx]*dispx)*dispy1
					 +(pdir[ix+(iy+1)*ndirx]*dispx1
					   +pdir[ix+1+(iy+1)*ndirx]*dispx)*dispy));
		    }
		}
	    }else{
		for(int iy=iy0; iy<ny; iy+=stepy){
		    for(int ix=ix0; ix<nx; ix+=stepx){
			pdir[ix+iy*ndirx]+=
			    alpha*(+(pps[ix+    iy*npsx]*dispx1
				     +pps[ix+1+    iy*npsx]*dispx)*dispy1
				   +(pps[ix+(iy+1)*npsx]*dispx1
				     +pps[ix+1+(iy+1)*npsx]*dispx)*dispy);
		    }
		}
	    }
	}else{//Generic
	    if(cc){
		/* Question: For cubic spline, don't we have to test whether pps
		   is within boundary?*/
		if(trans=='t'){
		    if(Z(fabs)(xratio-0.5f)<EPS && Z(fabs)(yratio-.5f)<EPS){
			//do without atomic operations.
			const int nxin=ceil(nx*xratio)+1;
			const int nyin=ceil(ny*xratio)+1;
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
			// shared memory to store them to avoid register spill.
			Real *const &fx=(Real*)shared.fx;
			Real *const &fy=(Real*)shared.fy;
			if(threadIdx.x==0 && threadIdx.y==0){
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
				    if(ky2>ymindir && ky2<ymaxdir){
#pragma unroll
					for(int kx=-4; kx<4; kx++){
					    const int kx2=kx+xcent;
					    Real wt;
					    if(kx2>xmindir && kx2<xmaxdir && (wt=fy[ky+4]*fx[kx+4])>EPS){
						sum+=wt*pdir[(kx2)+(ky2)*ndirx];
					    }
					}
				    }
				}
				//Need atomic because different layers have different offset.
				if(nn>1){
				    atomicAdd(&pps[mx+my*npsx],sum*alpha);
				}else{
				    pps[mx+my*npsx]+=sum*alpha;
				}
			    }
			}
		    }else{
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
				const int ky1=(ymaxps-iy)< 3?(ymaxps-iy): 3;
				for(int ky=ky0; ky<ky1; ky++){
				    int kx0=(xminps-ix)>-1?(xminps-ix):-1;
				    int kx1=(xmaxps-ix)< 3?(xmaxps-ix): 3;
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
		}else{//CC, non trans
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
			    const int ky1=(ymaxps-iy)< 3?(ymaxps-iy): 3;
			    for(int ky=ky0; ky<ky1; ky++){
				int kx0=(xminps-ix)>-1?(xminps-ix):-1;
				int kx1=(xmaxps-ix)< 3?(xmaxps-ix): 3;
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
	    }else if(trans=='t'){
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
			    atomicAdd(&pps[kx+      ky*npsx], temp*wt);
			if((wt=fracx*(1.f-fracy))>EPS)
			    atomicAdd(&pps[kx+1    +ky*npsx], temp*wt);
			if((wt=(1.f-fracx)*fracy)>EPS)
			    atomicAdd(&pps[kx+  (ky+1)*npsx], temp*wt);
			if((wt=fracx*fracy)>EPS)
			    atomicAdd(&pps[kx+1+(ky+1)*npsx], temp*wt);
		    }
		}
	    }else{
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
			    tmp+=(pps[kx+      ky*npsx]*(1.f-fracx)+
				  pps[kx+1+    ky*npsx]*fracx)*wt;
			}
			if((wt=fracy)>EPS){
			    tmp+=(pps[kx  +(ky+1)*npsx]*(1.f-fracx)+
				  pps[kx+1+(ky+1)*npsx]*fracx)*wt;
			}
			pdir[ix+iy*ndirx]+=alpha*tmp;
		    }
		}
	    }
	}
    }/*for*/
}
/*
  Prepare data and copy to GPU so that one kernel (gpu_map2map_do) can handle multiple independent ray tracing.
  Forward: ps->dir (xloc -> wfs or dm->floc)
  Backward: dir->ps (wfs -> xloc or floc->dm)
*/

void gpu_map2map_prep(PROP_WRAP_T*res, const cugrid_t &g_dir, const cugrid_t &g_ps,
		      Real dispx, Real dispy, const curmat &cc){
    assert(g_ps.ny!=1);
    const Real uxi=1.f/g_ps.dx;
    const Real uyi=1.f/g_ps.dy;
    Real xratio=g_dir.dx*uxi;
    Real yratio=g_dir.dy*uyi;
    /*remove round off errors that causes trouble in nx, ny.*/
    if(Z(fabs)(xratio-1.f)<EPS && Z(fabs)(yratio-1.f)<EPS){
	xratio=yratio=1.f;
	if(!cc && !res->isreverse){
	    res->reverse=new PROP_WRAP_T;
	    res->reverse->isreverse=1;
	    gpu_map2map_prep(res->reverse, g_ps, g_dir, -dispx, -dispy, cc);
	}
    }else if(Z(fabs)(xratio-0.5f)<EPS && Z(fabs)(yratio-0.5f)<EPS){
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
    /*if output is bigger than input. */
    if(dispx<0){
	offx1=(int)ceilf(-dispx*xratio1);
	dispx+=offx1*xratio;
    }
    if(dispy<0){
	offy1=(int)ceilf(-dispy*yratio1);
	dispy+=offy1*yratio;
    }
    /*convert offset into input grid coordinate. -EPS to avoid laying on the last point. */
    int nx=(int)Z(floor)((nxps-1-dispx-EPS)*xratio1+1);
    int ny=(int)Z(floor)((nyps-1-dispy-EPS)*yratio1+1);
    //Sanity check.
    while((nx-1)*xratio+dispx+1>=nxps){
	//warning("out point is at %g with ratio %g, and input %d\n",
	//	(nx-1)*xratio+dispx, xratio, nxps-1);
	nx--;
    }
    while((nx)*xratio+dispx+1<nxps){
	//warning("out point is at %g with ratio %g, and input %d\n",
	//	(nx-1)*xratio+dispx, xratio, nxps-1);
	nx++;
    }
    while((ny-1)*yratio+dispy+1>=nyps){
	//warning("out point is at %g with ratio %g, and input %d\n",
	//	(ny-1)*yratio+dispy, yratio, nyps-1);
	ny--;
    }
    while((ny)*yratio+dispy+1<nyps){
	//warning("out point is at %g with ratio %g, and input %d\n",
	//	(ny-1)*yratio+dispy, yratio, nyps-1);
	ny++;
    }
    if(nx>nxdir-offx1) nx=nxdir-offx1;
    if(ny>nydir-offy1) ny=nydir-offy1;
    int offx2=(int)Z(floor)(dispx); dispx-=offx2;/*for input. coarse sampling. */
    int offy2=(int)Z(floor)(dispy); dispy-=offy2;
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
    }else{
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
}

void gpu_map2map(cumap_t &out, const cumap_t &in, Real dispx, Real dispy, Real alpha, const curmat &cc, char trans){
    PROP_WRAP_T wrap;
    PROP_WRAP_T *wrap_gpu;
    cudaMalloc(&wrap_gpu, sizeof(PROP_WRAP_T));
    gpu_map2map_prep(&wrap, out, in, dispx, dispy, cc);
    wrap.togpu(wrap_gpu);
    Real **p;
    cudaMalloc(&p, sizeof(Real*)*2);
    const Real *tmp[2]={out(), in()};
    cudaMemcpy(p, tmp, sizeof(Real*)*2, cudaMemcpyHostToDevice);
    gpu_map2map_do<<<dim3(4,4,1),dim3(PROP_WRAP_TX,4),0,0>>>
	(wrap_gpu, p, p+1, 1, 1, alpha, 0, trans);
    cudaMemcpy(&wrap, wrap_gpu, sizeof(PROP_WRAP_T), cudaMemcpyDeviceToHost);
    if(wrap.reverse){
	cudaFree(wrap.reverse);
    }
    cudaFree(wrap_gpu);
}
