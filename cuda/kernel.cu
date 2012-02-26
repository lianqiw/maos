/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "kernel.h"
/**
   A few kernels.
*/


__global__ void set_do(float *a, float alpha, int n){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=alpha;
    }
}

__global__ void add_ptt_do(float *restrict opd, float (*restrict loc)[2], 
			   int n, float pis, float tx, float ty){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	opd[i]+=pis+loc[i][0]*tx+loc[i][1]*ty;
    }
}

__global__ void add_ngsmod_do(float *restrict opd, float (*restrict loc)[2], int n, 
			      float m0, float m1, float m2, float m3, float m4,
			      float thetax, float thetay, float scale, float ht, float MCC_fcp, float alpha
			      ){
    float scale1=1.f-scale;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	float x=loc[i][0];
	float y=loc[i][1];
	float xy=x*y;
	float x2=x*x;
	float y2=y*y;
	opd[i]+= alpha*(+x*m0
			+y*m1
			+m2*((x2+y2-MCC_fcp)*scale1-2*scale*ht*(thetax*x+thetay*y))
			+m3*((x2-y2)*scale1 - 2*scale*ht*(thetax*x-thetay*y))
			+m4*(xy*scale1-scale*ht*(thetay*x+thetax*y)));
    }
}
__global__ void inn_do(float *restrict res, const float *a, const float *b, const int n){
    extern __shared__ float sb[];
    sb[threadIdx.x]=0;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	sb[threadIdx.x]+=a[i]*b[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    sb[threadIdx.x]+=sb[threadIdx.x+step];
	}
    }
    if(threadIdx.x==0){
	atomicAdd(res, sb[0]);
    }
}
