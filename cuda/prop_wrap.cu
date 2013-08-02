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
#include <cuda.h>
extern "C"
{
#include "gpu.h"
}
#include "utils.h"
#include "accphi.h"
#include "prop_wrap.h"
/*
  One kernel that handles multiple layers/directions.
  Forward: Propagate from XLOC to WFS (Parallel across each WFS).  
  Backward (Transpose): Propagate from WFS to XLOC (Parallel across each Layer)
  
  let XLOC be ps (phase screen)
  let WFS be dir (direction)

  alpha2: additional scaling defined on dimension dir.
*/

__global__ void gpu_prop_grid_do(PROP_WRAP_T *data, float **pdirs, float **ppss, int ndir, int nps, float alpha1, float *alpha2, char trans){
    int nn;
    if(ndir==1){
	assert(gridDim.z==nps);
	nn=1;
    }else if(trans=='t'){
	assert(gridDim.z==nps);
	nn=ndir;
    }else{
	assert(gridDim.z==ndir);
	nn=nps;
    }
    const int ix0=blockIdx.x*blockDim.x+threadIdx.x;
    const int iy0=blockIdx.y*blockDim.y+threadIdx.y;
    const int stepx=blockDim.x*gridDim.x;
    const int stepy=blockDim.y*gridDim.y;
  
    for(int ii=0; ii<nn; ii++){
	PROP_WRAP_T *datai;
	int ips, idir;
	if(ndir==1){//plane to plane. no direction
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
	const bool match=fabsf(datai->xratio-1.f)<EPS && fabsf(datai->yratio-1.f)<EPS;
	if(!datai->cc && trans=='t' && match){
	    datai=datai->reverse;
	}
	const int nx=datai->nx;
	const int ny=datai->ny;
	if(nx==0) continue;//skip empty wfs
	float *restrict pps=ppss[ips]+datai->offps;
	float *restrict pdir=pdirs[idir]+datai->offdir;
	const float alpha=alpha2?(alpha2[idir]*alpha1):alpha1;
	const float *cc=datai->cc;
	const int nxdir=datai->nxdir;
	const int nxps=datai->nxps;
	if(!cc && match){
	    //Matched bilinear propagation. always forward prop. 
	    const float fracx=datai->dispx;
	    const float fracy=datai->dispy;
	    const float fracx1=1.f-fracx;//this reduces register usage.
	    const float fracy1=1.f-fracy;
	    /*During reverse operation, for different idir, the offo is
	      different causing same thread to handle different memory in pdir
	      for different idir. This causes problem with pdir synchronization/atomic operation*/
	    if(trans=='t'){
		for(int iy=iy0; iy<ny; iy+=stepy){
		    for(int ix=ix0; ix<nx; ix+=stepx){
			atomicAdd(&pps[ix+iy*nxps],
				  alpha*(+(pdir[ix+    iy*nxdir]*fracx1+pdir[ix+1+    iy*nxdir]*fracx)*fracy1
					 +(pdir[ix+(iy+1)*nxdir]*fracx1+pdir[ix+1+(iy+1)*nxdir]*fracx)*fracy));
		    }
		}
	    }else{
		for(int iy=iy0; iy<ny; iy+=stepy){
		    for(int ix=ix0; ix<nx; ix+=stepx){
			pdir[ix+iy*nxdir]+=
			    alpha*(+(pps[ix+    iy*nxps]*fracx1+pps[ix+1+    iy*nxps]*fracx)*fracy1
				   +(pps[ix+(iy+1)*nxps]*fracx1+pps[ix+1+(iy+1)*nxps]*fracx)*fracy);
		    }
		}
	    }
	}else{//Generic
	    const float xratio=datai->xratio;
	    const float yratio=datai->yratio;
	    const float dispx=datai->dispx;
	    const float dispy=datai->dispy;
	    if(cc){
		/* Question: For cubic spline, don't we have to test whether pps
		   is within boundary?*/
		const bool match2=fabsf(datai->xratio-0.5f)<EPS && fabsf(datai->yratio-.5f)<EPS;
		if(trans=='t'){
		    if(match2){//do without atomic operations.
			const int nxin=ceil(nx*xratio)+1;
			const int nyin=ceil(ny*xratio)+1;
			//const int nxin=datai->nxps-datai->offpsx;
			//const int nyin=datai->nyps-datai->offpsy;
			const int xmaxdir=datai->nxdir-datai->offdirx;
			const int ymaxdir=datai->nydir-datai->offdiry;
			const int xmindir=-datai->offdirx-1;
			const int ymindir=-datai->offdiry-1;
			int offx=0, offy=0;
			float xc=dispx, yc=dispy;
			if(dispx>=0.5){
			    offx=-1;
			    xc=dispx-0.5f;
			}
			if(dispy>=0.5){
			    offy=-1;
			    yc=dispy-0.5f;
			}
		
			const float xc2=xc+0.5f;
			const float yc2=yc+0.5f;
			float fy[8], fx[8];

			fy[0]=(cc[3]+cc[4]*yc)*yc*yc;
			fy[1]=(cc[3]+cc[4]*yc2)*yc2*yc2;
			fy[2]=(cc[1]+cc[2]*(1-yc))*(1-yc)*(1-yc)+cc[0];
			fy[3]=(cc[1]+cc[2]*(1-yc2))*(1-yc2)*(1-yc2)+cc[0];
			fy[4]=(cc[1]+cc[2]*yc)*yc*yc+cc[0];
			fy[5]=(cc[1]+cc[2]*yc2)*yc2*yc2+cc[0];
			fy[6]=(cc[3]+cc[4]*(1-yc))*(1-yc)*(1-yc);
			fy[7]=(cc[3]+cc[4]*(1-yc2))*(1-yc2)*(1-yc2);

			fx[0]=(cc[3]+cc[4]*xc)*xc*xc;
			fx[1]=(cc[3]+cc[4]*xc2)*xc2*xc2;
			fx[2]=(cc[1]+cc[2]*(1-xc))*(1-xc)*(1-xc)+cc[0];
			fx[3]=(cc[1]+cc[2]*(1-xc2))*(1-xc2)*(1-xc2)+cc[0];
			fx[4]=(cc[1]+cc[2]*xc)*xc*xc+cc[0];
			fx[5]=(cc[1]+cc[2]*xc2)*xc2*xc2+cc[0];
			fx[6]=(cc[3]+cc[4]*(1-xc))*(1-xc)*(1-xc);
			fx[7]=(cc[3]+cc[4]*(1-xc2))*(1-xc2)*(1-xc2);

			for(int my=iy0; my<nyin; my+=stepy){
			    int ycent=2*my+offy;
			    for(int mx=ix0; mx<nxin; mx+=stepx){
				float sum=0;
				int xcent=2*mx+offx;
#pragma unroll
				for(int ky=-4; ky<4; ky++){
				    int ky2=ky+ycent;
				    if(ky2>ymindir && ky2<ymaxdir){
#pragma unroll
					for(int kx=-4; kx<4; kx++){
					    int kx2=kx+xcent;
					    if(kx2>xmindir && kx2<xmaxdir){
						sum+=fy[ky+4]*fx[kx+4]*pdir[(kx2)+(ky2)*nxdir];
					    }
					}
				    }
				}
				//Need atomic because different layers have different offset.
				if(nn>1){
				    atomicAdd(&pps[mx+my*nxps],sum*alpha);
				}else{
				    pps[mx+my*nxps]+=sum*alpha;
				}
			    }
			}
		    }else{
			for(int my=iy0; my<ny; my+=stepy){
			    float jy;
			    float y=modff(dispy+my*yratio, &jy);
			    int iy=(int)jy;	
			    float fy[4];
			    fy[0]=(1.f-y)*(1.f-y)*(cc[3]+cc[4]*(1.f-y));			
			    fy[1]=cc[0]+y*y*(cc[1]+cc[2]*y);			
			    fy[2]=cc[0]+(1.f-y)*(1.f-y)*(cc[1]+cc[2]*(1.f-y));			
			    fy[3]=y*y*(cc[3]+cc[4]*y);		
			    for(int mx=ix0; mx<nx; mx+=stepx){
				float jx;
				float x=modff(dispx+mx*xratio, &jx);
				int ix=(int)jx;
				float fx[4];
				float value=pdir[mx+my*nxdir]*alpha;
				/*cc need to be in device memory for sm_13 to work.*/
				fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));			
				fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);			
				fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));			
				fx[3]=x*x*(cc[3]+cc[4]*x);	
	
#pragma unroll
				for(int ky=-1; ky<3; ky++){
				    for(int kx=-1; kx<3; kx++){
					atomicAdd(&pps[(iy+ky)*nxps+(kx+ix)], fx[kx+1]*fy[ky+1]*value);
				    }
				}
			    }
			}
		    }
		}else{
		    for(int my=iy0; my<ny; my+=stepy){
			float jy;
			float y=modff(dispy+my*yratio, &jy);
			int iy=(int)jy;	
			float fy[4];
			fy[0]=(1.f-y)*(1.f-y)*(cc[3]+cc[4]*(1.f-y));			
			fy[1]=cc[0]+y*y*(cc[1]+cc[2]*y);			
			fy[2]=cc[0]+(1.f-y)*(1.f-y)*(cc[1]+cc[2]*(1.f-y));			
			fy[3]=y*y*(cc[3]+cc[4]*y);	

			for(int mx=ix0; mx<nx; mx+=stepx){
			    float jx;
			    float x=modff(dispx+mx*xratio, &jx);
			    int ix=(int)jx;
			    float fx[4];
			    float sum=0;
			    /*cc need to be in device memory for sm_13 to work.*/
			    fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));			
			    fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);			
			    fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));			
			    fx[3]=x*x*(cc[3]+cc[4]*x);		
			    
#pragma unroll
			    for(int ky=-1; ky<3; ky++){
				for(int kx=-1; kx<3; kx++){
				    sum+=fx[kx+1]*fy[ky+1]*pps[(iy+ky)*nxps+(kx+ix)];
				}
			    }
			    pdir[mx+my*nxdir]+=sum*alpha;
			}
		    }
		}/*else trans*/
	    }else if(trans=='t'){
		for(int iy=iy0; iy<ny; iy+=stepy){
		    float temp;
		    float fracy=modff(dispy+iy*yratio, &temp);
		    int ky=(int)temp;
		    for(int ix=ix0; ix<nx; ix+=stepx){
			float jx;
			float fracx=modff(dispx+ix*xratio, &jx);
			int kx=(int)jx;
			float temp=pdir[ix+iy*nxdir]*alpha;
			atomicAdd(&pps[kx+      ky*nxps], temp*(1.f-fracx)*(1.f-fracy));
			atomicAdd(&pps[kx+1    +ky*nxps], temp*fracx*(1.f-fracy));
			atomicAdd(&pps[kx+  (ky+1)*nxps], temp*(1.f-fracx)*fracy);
			atomicAdd(&pps[kx+1+(ky+1)*nxps], temp*fracx*fracy);
		    }
		}
	    }else{
		for(int iy=iy0; iy<ny; iy+=stepy){
		    float jy;
		    float fracy=modff(dispy+iy*yratio, &jy);
		    int ky=(int)jy;
		    
		    for(int ix=ix0; ix<nx; ix+=stepx){
			float fracx; int kx;
			
			float jx;
			fracx=modff(dispx+ix*xratio, &jx);
			kx=(int)jx;
			
			pdir[ix+iy*nxdir]+=
			    alpha*(+(pps[kx+      ky*nxps]*(1.f-fracx)+
				     pps[kx+1+    ky*nxps]*fracx)*(1.f-fracy)
				   +(pps[kx  +(ky+1)*nxps]*(1.f-fracx)+
				     pps[kx+1+(ky+1)*nxps]*fracx)*fracy);
		    }
		}
	    }
	}
    }/*for*/
}
/*
  Prepare data and copy to GPU so that one kernel (gpu_prop_grid_do) can handle multiple independent ray tracing.
  Forward: ps->dir (xloc -> wfs or dm->floc)
  Backward: dir->ps (wfs -> xloc or floc->dm)
*/

void gpu_prop_grid_prep(PROP_WRAP_T*res, 
			const cugrid_t &g_dir, const cugrid_t &g_ps,
			float dispx, float dispy, curmat *cc){
    assert(g_ps.ny!=1);
    const float uxi=1.f/g_ps.dx;
    const float uyi=1.f/g_ps.dy;
    float xratio=g_dir.dx*uxi;
    float yratio=g_dir.dy*uyi;
    /*remove round off errors that causes trouble in nx, ny.*/
    if(fabs(xratio-1.f)<EPS && fabs(yratio-1.f)<EPS){
	xratio=yratio=1.f;
	if(!cc && !res->isreverse){
	    res->reverse=new PROP_WRAP_T;
	    res->reverse->isreverse=1;
	    gpu_prop_grid_prep(res->reverse, g_ps, g_dir, -dispx, -dispy, NULL);
	}
    }else if(fabs(xratio-0.5f)<EPS && fabs(yratio-0.5f)<EPS){
	xratio=yratio=0.5f;
    }
    const int nxdir=g_dir.nx;
    const int nydir=g_dir.ny;
    const int nxps=g_ps.nx;
    const int nyps=g_ps.ny;
    const float xratio1=1.f/xratio;
    const float yratio1=1.f/yratio;
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
    int nx=(int)floorf((nxps-dispx-EPS)*xratio1);
    int ny=(int)floorf((nyps-dispy-EPS)*yratio1);
    
    if(nx>nxdir-offx1) nx=nxdir-offx1;
    if(ny>nydir-offy1) ny=nydir-offy1;
    int offx2=(int)floorf(dispx); dispx-=offx2;/*for input. coarse sampling. */
    int offy2=(int)floorf(dispy); dispy-=offy2;
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
    res->cc=cc?cc->p:NULL;
}
