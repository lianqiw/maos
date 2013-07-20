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
/*
  Ray tracing with matched spacing. Reverse, from out to in. out is xloc, in is ploc.
*/
__global__ static void prop_grid_match_do(float *restrict out, int nxo,
					  const float *restrict in, int nxi, 
					  float fracx, float fracy,
					  float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    float fracx1=1.f-fracx;
    float fracy1=1.f-fracy;
    for(int iy=blockIdx.y*blockDim.y+threadIdx.y; iy<ny; iy+=stepy){
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=stepx){
	    out[ix+iy*nxo]+=
		alpha*(+(in[ix+    iy*nxi]*fracx1+in[ix+1+    iy*nxi]*fracx)*fracy1
		       +(in[ix+(iy+1)*nxi]*fracx1+in[ix+1+(iy+1)*nxi]*fracx)*fracy);
	}
    }
}

__global__ static void prop_grid_nomatch_do(float *restrict out, int nxo, 
					    const float *restrict in, int nxi, 
					    float dispx, float dispy, float xratio, float yratio,
					    float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    for(int iy=blockIdx.y*blockDim.y+threadIdx.y; iy<ny; iy+=stepy){
	float jy;
	float fracy=modff(dispy+iy*yratio, &jy);
	int ky=(int)jy;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=stepx){
	    float jx;
	    float fracx=modff(dispx+ix*xratio, &jx);
	    int kx=(int)jx;
	    out[ix+iy*nxo]+=
		alpha*(+(in[kx+      ky*nxi]*(1.f-fracx)+
			 in[kx+1+    ky*nxi]*fracx)*(1.f-fracy)
		       +(in[kx  +(ky+1)*nxi]*(1.f-fracx)+
			 in[kx+1+(ky+1)*nxi]*fracx)*fracy);
	}
    }
}
__global__ static void prop_grid_nomatch_trans_do(const float *restrict out, int nxo, 
						  float *restrict in, int nxi,
						  float dispx, float dispy, float xratio, float yratio,
						  float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    for(int iy=blockIdx.y*blockDim.y+threadIdx.y; iy<ny; iy+=stepy){
	float jy;
	float fracy=modff(dispy+iy*yratio, &jy);
	int ky=(int)jy;
	for(int ix=blockIdx.x*blockDim.x+threadIdx.x; ix<nx; ix+=stepx){
	    float jx;
	    float fracx=modff(dispx+ix*xratio, &jx);
	    int kx=(int)jx;
	    float temp=out[ix+iy*nxo]*alpha;
	    atomicAdd(&in[kx+      ky*nxi], temp*(1.f-fracx)*(1.f-fracy));
	    atomicAdd(&in[kx+1    +ky*nxi], temp*fracx*(1.f-fracy));
	    atomicAdd(&in[kx+  (ky+1)*nxi], temp*(1.f-fracx)*fracy);
	    atomicAdd(&in[kx+1+(ky+1)*nxi], temp*fracx*fracy);
	}
    }
}

/*
  Ray tracing with over sampling. Reverse, from out to in. out is xloc, in is
  ploc. confirmed to agree with HXW'.  */
/*
__global__ static void prop_grid_os2_trans_old_do(const float *restrict out, int nxout,
						  float *restrict in, int nxin, 
						  float fracx, float fracy,
						  float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    int ax=fracx<0.5f?0:1;
    int ay=fracy<0.5f?0:1;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<(ny+1)/2; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<(nx+1)/2; ix+=stepx){

#pragma unroll 
	    for(int by=0; by<2; by++){
		int iy2=iy*2+by;
		int iy3=iy+ay*by;
		float fracy2=fracy+(0.5f-ay)*by;
		float fracy21=1.f-fracy2;
#pragma unroll 
		for(int bx=0; bx<2; bx++){
		    int ix2=ix*2+bx;
		    int ix3=ix+ax*bx;
		    float fracx2=fracx+(0.5f-ax)*bx;
		    float fracx21=1.f-fracx2;
		    if(ix2<nx && iy2<ny){
			float a=out[ix2+(iy2)*nxout]*alpha;
			atomicAdd(&in[ix3+    (iy3)*nxin], a*fracx21*fracy21);
			atomicAdd(&in[ix3+1+  (iy3)*nxin], a*fracx2*fracy21);
			atomicAdd(&in[ix3+  (iy3+1)*nxin], a*fracx21*fracy2);
			atomicAdd(&in[ix3+1+(iy3+1)*nxin], a*fracx2*fracy2);
		    }
		}
	    }
	}
    }
}*/
/**
   Do a single column (along x).
*/

__global__ static void prop_grid_os2_trans_col_do(float *restrict out, int nxout,
						  float *restrict in, int nxin, 
						  float fracx, float fracy2,
						  float alpha, int nx){
    int stepx=blockDim.x*gridDim.x;
    int ax=fracx<0.5f?0:1;
    const int iy3=0;
    float fracy21=1.f-fracy2;
    for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<(nx+1)/2; ix+=stepx){
	for(int bx=0; bx<2; bx++){
	    int ix2=ix*2+bx;
	    if(ix2<nx){
		int ix3=ix+ax*bx;
		float fracx2=fracx+(0.5f-ax)*bx;
		float fracx21=1.f-fracx2;
		float a = out[ix2]*alpha;
		atomicAdd(&in[ix3+    (iy3)*nxin], a*fracx21*fracy21);
		atomicAdd(&in[ix3+1+  (iy3)*nxin], a*fracx2*fracy21);
		atomicAdd(&in[ix3+  (iy3+1)*nxin], a*fracx21*fracy2);
		atomicAdd(&in[ix3+1+(iy3+1)*nxin], a*fracx2*fracy2);
	    }
	}
    }
}
/**
   Do a single row (along y).
*/
__global__ static void prop_grid_os2_trans_row_do(float *restrict out, int nxout,
						  float *restrict in, int nxin, 
						  float fracx2, float fracy,
						  float alpha, int ny){
    int stepy=blockDim.y*gridDim.y;
    int ay=fracy<0.5f?0:1;
    const int ix3=0;
    float fracx21=1.f-fracx2;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<(ny+1)/2; iy+=stepy){
	for(int by=0; by<2; by++){
	    int iy2=iy*2+by;
	    if(iy2<ny){
		int iy3=iy+ay*by;
		float fracy2=fracy+(0.5f-ay)*by;
		float fracy21=1.f-fracy2;
		float a = out[iy2*nxout]*alpha;
		atomicAdd(&in[ix3+    (iy3)*nxin], a*fracx21*fracy21);
		atomicAdd(&in[ix3+1+  (iy3)*nxin], a*fracx2*fracy21);
		atomicAdd(&in[ix3+  (iy3+1)*nxin], a*fracx21*fracy2);
		atomicAdd(&in[ix3+1+(iy3+1)*nxin], a*fracx2*fracy2);
	    }
	}
    }
}
/**
   Try to improve prop_grid_os2_trans_do using shared memory.

   2012-05-04: Problem found: When nx or ny is not multiple of block size,
   during the last set, the threads that are outside stop running because of
   test, iy<ny, ix<nx; Then the last col/row is cachein left not merged.

*/
__global__ static void prop_grid_os2_trans_share_do(float *restrict out, int nxout,
						    float *restrict in, int nxin, 
						    float fracx, float fracy,
						    float alpha, int nx, int ny){
    extern __shared__ float cachein[];/*caching. */
    const int ind=threadIdx.x+threadIdx.y*blockDim.x;
    cachein[ind]=0;
    __syncthreads();
    const int stepx=blockDim.x*gridDim.x;
    const int stepy=blockDim.y*gridDim.y;
    const int nyf=((ny+stepy-1)/stepy)*stepy;
    const int nxf=((nx+stepx-1)/stepx)*stepx;
    for(int iy=(blockIdx.y*blockDim.y+threadIdx.y); iy<nyf; iy+=stepy){
	for(int ix=(blockIdx.x*blockDim.x+threadIdx.x); ix<nxf; ix+=stepx){
	    int ix2=ix*2;
	    int iy2=iy*2;
	    if(iy<ny && ix<nx){
		float v_0_0=alpha*(+(+out[ix2+      iy2*nxout]*(1.f-fracx)
				     +out[ix2+1+    iy2*nxout]*(0.5f-fracx))*(1.f-fracy)
				   +(+out[ix2+  (iy2+1)*nxout]*(1.f-fracx)
				     +out[ix2+1+(iy2+1)*nxout]*(0.5f-fracx))*(0.5f-fracy));
		float v_1_0=alpha*(+(+out[ix2+      iy2*nxout]*(fracx)
				     +out[ix2+1+    iy2*nxout]*(0.5f+fracx))*(1.f-fracy)
				   +(+out[ix2+  (iy2+1)*nxout]*(fracx)
				     +out[ix2+1+(iy2+1)*nxout]*(0.5f+fracx))*(0.5f-fracy));
		float v_0_1=alpha*(+(+out[ix2+      iy2*nxout]*(1.f-fracx)
				     +out[ix2+1+    iy2*nxout]*(0.5f-fracx))*(fracy)
				   +(+out[ix2+  (iy2+1)*nxout]*(1.f-fracx)
				     +out[ix2+1+(iy2+1)*nxout]*(0.5f-fracx))*(0.5f+fracy));
		float v_1_1=alpha*(+(+out[ix2+      iy2*nxout]*(fracx)
				     +out[ix2+1+    iy2*nxout]*(0.5f+fracx))*(fracy)
				   +(+out[ix2+  (iy2+1)*nxout]*(fracx)
				     +out[ix2+1+(iy2+1)*nxout]*(0.5f+fracx))*(0.5f+fracy));
		atomicAdd(&cachein[ind], v_0_0);
		if(threadIdx.x+1==blockDim.x){/*edge*/
		    atomicAdd(&in[ix+1      +iy*nxin], v_1_0);
		    atomicAdd(&in[(ix+1)+(iy+1)*nxin], v_1_1);
		}else{
		    atomicAdd(&cachein[ind+1], v_1_0);
		    if(threadIdx.y+1==blockDim.y){
			atomicAdd(&in[(ix+1)+(iy+1)*nxin], v_1_1);
		    }else{
			atomicAdd(&cachein[ind+1+blockDim.x], v_1_1);
		    }
		}
		if(threadIdx.y+1==blockDim.y){
		    atomicAdd(&in[ix    +(iy+1)*nxin], v_0_1);
		}else{
		    atomicAdd(&cachein[ind+blockDim.x], v_0_1);
		}
	    }
	    /*atomicAdd(&in[ix        +iy*nxin], v_0_0); */
	    /*atomicAdd(&in[ix+1      +iy*nxin], v_1_0); */
	    /*atomicAdd(&in[ix    +(iy+1)*nxin], v_0_1); */
	    /*atomicAdd(&in[(ix+1)+(iy+1)*nxin], v_1_1); */
	    
	    __syncthreads();
	    atomicAdd(&in[ix+iy*nxin], cachein[ind]);
	    cachein[ind]=0;
	    __syncthreads();
	}
    }
}
/**
   Do the ray tracing
   from in to out if trans=='n'
   from out to in if trans=='t'
*/
void gpu_prop_grid(curmat *out, const cugrid_t &go,
		   curmat *in, const cugrid_t &gi,
		   float dispx, float dispy,
		   float alpha, char trans, cudaStream_t stream){
    assert(in->ny!=1);
    const float uxi=1.f/gi.dx;
    const float uyi=1.f/gi.dy;
    float xratio=go.dx*uxi;
    float yratio=go.dy*uyi;
    int match=0,match2=0;
    /*remove round off errors that causes trouble in nx, ny.*/
    if(fabs(xratio-1.f)<EPS && fabs(yratio-1.f)){
	xratio=yratio=1.f;
	match=1;
    }else if(fabs(xratio-0.5f)<EPS && fabs(yratio-0.5f)){
	xratio=yratio=0.5f;
	match2=1;
    }
    if(match && trans=='t'){
	gpu_prop_grid(in, gi, out, go, -dispx, -dispy, alpha,'n', stream);
	return;
    }
    const int nxo=out->nx;
    const int nyo=out->ny;
    const int nxi=in->nx;
    const int nyi=in->ny;
    const float xratio1=1.f/xratio;
    const float yratio1=1.f/yratio;
    /*offset of origin in input grid spacing. */
    dispx=(dispx-gi.ox+go.ox)*uxi;
    dispy=(dispy-gi.oy+go.oy)*uyi;
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
    int nx=(int)floorf((nxi-dispx-EPS)*xratio1);
    int ny=(int)floorf((nyi-dispy-EPS)*yratio1);
    
    if(nx>nxo-offx1) nx=nxo-offx1;
    if(ny>nyo-offy1) ny=nyo-offy1;
    int offx2=(int)floorf(dispx); dispx-=offx2;/*for input. coarse sampling. */
    int offy2=(int)floorf(dispy); dispy-=offy2;
    /*
      info("offx1=%d, offy1=%d, offx2=%d, offy2=%d\n", offx1, offy1, offx2, offy2);
      info("nxi=%d, nyi=%d, nxo=%d, nyo=%d\n", nxi, nyi, nxo, nyo);
      info("dispx=%g, dispy=%g\n", dispx, dispy);
      info("nx=%d, ny=%d\n", nx, ny);
    */
    if(trans!='t'){
	if(match){
	    prop_grid_match_do<<<DIM2(nx, ny, NTH2), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, alpha, nx, ny);
	}else{
	    prop_grid_nomatch_do<<<DIM2(nx, ny, NTH2), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, xratio, yratio, alpha, nx, ny);
	}
    }else{
	if(match){
	    error("Please revert the input/output and call with trans='n'\n");
	}else if(match2){
	    if(dispy>0.5f){/*do a single col first. */
		prop_grid_os2_trans_col_do<<<DIM2(nx,1,256),0,stream>>>
		    (out->p+offy1*nxo+offx1, nxo,
		     in->p+offy2*nxi+offx2, nxi, 
		     dispx, dispy, alpha, nx);
		dispy-=0.5f;
		offy1+=1; ny-=1;
		offy2+=1;
	    }
	 
	    if(dispx>0.5f){/*do a single row first */
		prop_grid_os2_trans_row_do<<<DIM2(1,ny,256),0,stream>>>
		    (out->p+offy1*nxo+offx1, nxo,
		     in->p+offy2*nxi+offx2, nxi, 
		     dispx, dispy, alpha, ny);
		dispx-=0.5f;
		offx1+=1; nx-=1;
		offx2+=1;
	    }
	    int nx2=nx>>1; 
	    int ny2=ny>>1;
	    
#define BS 16 // cannot be 32. /
	    prop_grid_os2_trans_share_do<<<DIM2(nx2, ny2, BS), (BS)*(BS)*sizeof(float), stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, alpha, nx2, ny2);
	    if(ny & 1 == 1){/*do the last col that is left over */
		prop_grid_os2_trans_col_do<<<DIM2(nx,1,256),0,stream>>>
		    (out->p+(offy1+ny2*2)*nxo+offx1, nxo,
		     in->p+(offy2+ny2)*nxi+offx2, nxi, 
		     dispx, dispy, alpha, nx);
	    }
	    if(nx & 1 == 1){/*do the last row that is left over */
		prop_grid_os2_trans_row_do<<<DIM2(1,ny,256),0,stream>>>
		    (out->p+offy1*nxo+(offx1+nx2*2), nxo,
		     in->p+offy2*nxi+(offx2+nx2), nxi, 
		     dispx, dispy, alpha, ny);
	    }
	}else{
	    prop_grid_nomatch_trans_do<<<DIM2(nx, ny, NTH2), 0, stream>>>
		(out->p+offy1*nxo+offx1, nxo,
		 in->p+offy2*nxi+offx2, nxi, 
		 dispx, dispy, xratio, yratio, alpha, nx, ny);
	}
    }
}
__global__ static void prop_grid_cubic_nomatch_do(float *restrict out, int nxo, 
						  const float *restrict in, int nxi, 
						  float dispx, float dispy, float xratio, float yratio, 
						  float *cc, float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    for(int my=blockIdx.y*blockDim.y+threadIdx.y; my<ny; my+=stepy){
	float jy;
	float y=modff(dispy+my*yratio, &jy);
	int iy=(int)jy;
	for(int mx=blockIdx.x*blockDim.x+threadIdx.x; mx<nx; mx+=stepx){
	    float jx;
	    float x=modff(dispx+mx*xratio, &jx);
	    int ix=(int)jx;
	    float fx[4],fy;
	    float sum=0;
	    /*cc need to be in device memory for sm_13 to work.*/
	    fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));			
	    fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);			
	    fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));			
	    fx[3]=x*x*(cc[3]+cc[4]*x);		
	    
	    fy=(1.f-y)*(1.f-y)*(cc[3]+cc[4]*(1.f-y)); 
#pragma unroll
	    for(int kx=-1; kx<3; kx++){
		sum+=fx[kx+1]*fy*in[(iy-1)*nxi+kx+ix];
	    }

	    fy=cc[0]+y*y*(cc[1]+cc[2]*y); 
#pragma unroll
	    for(int kx=-1; kx<3; kx++){
		sum+=fx[kx+1]*fy*in[iy*nxi+kx+ix];
	    }

	    fy=cc[0]+(1.f-y)*(1.f-y)*(cc[1]+cc[2]*(1.f-y)); 
#pragma unroll
	    for(int kx=-1; kx<3; kx++){
		sum+=fx[kx+1]*fy*in[(iy+1)*nxi+kx+ix];
	    }

	    fy=y*y*(cc[3]+cc[4]*y); 
#pragma unroll
	    for(int kx=-1; kx<3; kx++){
		sum+=fx[kx+1]*fy*in[(iy+2)*nxi+kx+ix];
	    }
	    out[mx+my*nxo]+=sum*alpha;
	}
    }
}
__global__ static void prop_grid_cubic_nomatch_trans_do(const float *restrict out, int nxo,
							float *restrict in, int nxi, 
							float dispx, float dispy, float xratio, float yratio,
							float *cc, float alpha, int nx, int ny){
    int stepx=blockDim.x*gridDim.x;
    int stepy=blockDim.y*gridDim.y;
    for(int my=blockIdx.y*blockDim.y+threadIdx.y; my<ny; my+=stepy){
	float jy;
	float y=modff(dispy+my*yratio, &jy);
	int iy=(int)jy;
	for(int mx=blockIdx.x*blockDim.x+threadIdx.x; mx<nx; mx+=stepx){
	    float jx;
	    float x=modff(dispx+mx*xratio, &jx);
	    int ix=(int)jx;
	    float fx[4],fy;
	    float value=out[mx+my*nxo]*alpha;
	    /*cc need to be in device memory for sm_13 to work.*/
	    fx[0]=(1.f-x)*(1.f-x)*(cc[3]+cc[4]*(1.f-x));			
	    fx[1]=cc[0]+x*x*(cc[1]+cc[2]*x);			
	    fx[2]=cc[0]+(1.f-x)*(1.f-x)*(cc[1]+cc[2]*(1.f-x));			
	    fx[3]=x*x*(cc[3]+cc[4]*x);		
	    
	    fy=(1.f-y)*(1.f-y)*(cc[3]+cc[4]*(1.f-y)); 
#pragma unroll
	    for(int kx=-1; kx<3; kx++){
		atomicAdd(&in[(iy-1)*nxi+kx+ix], fx[kx+1]*fy*value);
	    }

	    fy=cc[0]+y*y*(cc[1]+cc[2]*y); 
#pragma unroll
	    for(int kx=-1; kx<3; kx++){
		atomicAdd(&in[iy*nxi+kx+ix], fx[kx+1]*fy*value);
	    }

	    fy=cc[0]+(1.f-y)*(1.f-y)*(cc[1]+cc[2]*(1.f-y)); 
#pragma unroll
	    for(int kx=-1; kx<3; kx++){
		atomicAdd(&in[(iy+1)*nxi+kx+ix], fx[kx+1]*fy*value);
	    }

	    fy=y*y*(cc[3]+cc[4]*y); 
#pragma unroll
	    for(int kx=-1; kx<3; kx++){
		atomicAdd(&in[(iy+2)*nxi+kx+ix], fx[kx+1]*fy*value);
	    }
	}
    }
}
/**
   Do the ray tracing
   from in to out if trans=='n' (dm to ploc)
   from out to in if trans=='t' (ploc to dm)
*/
void gpu_prop_grid_cubic(curmat *out, const cugrid_t &go,
			 curmat *in, const cugrid_t &gi,
			 float dispx, float dispy, float *cc,
			 float alpha, char trans, cudaStream_t stream){
    assert(in->ny!=1);
    const float uxi=1.f/gi.dx;
    const float uyi=1.f/gi.dy;
    float xratio=go.dx*uxi;
    float yratio=go.dy*uyi;
    /*remove round off errors.*/
    if(fabs(xratio-1.f)<EPS && fabs(yratio-1.f)<EPS){
	xratio=yratio=1.f;
    }else if(fabs(xratio-0.5f)<EPS && fabs(yratio-0.5f)<EPS){
	xratio=yratio=0.5f;
    }
    const int nxo=out->nx;
    const int nyo=out->ny;
    const int nxi=in->nx;
    const int nyi=in->ny;
    const float xratio1=1.f/xratio;
    const float yratio1=1.f/yratio;
    /*offset of origin in input grid spacing. */
    dispx=(dispx-gi.ox+go.ox)*uxi;
    dispy=(dispy-gi.oy+go.oy)*uyi;
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
    int nx=(int)floorf((nxi-1-dispx-EPS)*xratio1)+1;
    int ny=(int)floorf((nyi-1-dispy-EPS)*yratio1)+1;

    if(nx>nxo-offx1) nx=nxo-offx1;
    if(ny>nyo-offy1) ny=nyo-offy1;
    int offx2=(int)floorf(dispx); dispx-=offx2;/*for input. coarse sampling. */
    int offy2=(int)floorf(dispy); dispy-=offy2;
    if(trans!='t'){
	prop_grid_cubic_nomatch_do<<<DIM2(nx, ny, NTH2), 0, stream>>>
	    (out->p+offy1*nxo+offx1, nxo,
	     in->p+offy2*nxi+offx2, nxi, 
	     dispx, dispy, xratio, yratio, cc, alpha, nx, ny);
    }else{
	prop_grid_cubic_nomatch_trans_do<<<DIM2(nx, ny, NTH2), 0, stream>>>
	    (out->p+offy1*nxo+offx1, nxo,
	     in->p+offy2*nxi+offx2, nxi, 
	     dispx, dispy, xratio, yratio, cc, alpha, nx, ny);
    }
}
