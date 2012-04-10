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
/*somehow I must test both CUDA_ARCH existance and version.*/


__global__ void set_do(float *a, float alpha, int n){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=alpha;
    }
}
__global__ void scale_do(float *restrict in, int n, float alpha){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	in[i]*=alpha;
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

/**
   add a vector to another, scaled by alpha and beta. all in device memory.
   a=a*alpha+b*beta;
*/
__global__ void add_do(float *restrict a, float *alpha1, float alpha2, 
		       const float *restrict b, float *beta1, float beta2, int n){
    float alpha=alpha1?*alpha1*alpha2:alpha2;
    float beta=beta1?*beta1*beta2:beta2;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=a[i]*alpha+b[i]*beta;
    }
}

__global__ void add_do(float *restrict a, const float *restrict b, float *beta1, float beta2, int n){
    float beta=beta1?*beta1*beta2:beta2;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]+=b[i]*beta;
    }
}
__global__ void add_do(float *restrict a, float *alpha1, float alpha2, const float *restrict b,  int n){
    float alpha=alpha1?*alpha1*alpha2:alpha2;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=a[i]*alpha+b[i];
    }
}

/**
   add a beta to a vector. 
*/
__global__ void add_do(float *vec, float beta, int n){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	vec[i]+=beta;
    }
}

__global__ void addcabs2_do(float *restrict a, float alpha, 
			    const fcomplex *restrict b, float beta, int n){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=a[i]*alpha+CABS2(b[i])*beta;
    }
}

/*reduction routines*/

__global__ void sum_do(float *restrict res, const float *a, const int n){
    extern __shared__ float sb[];
    sb[threadIdx.x]=0;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	sb[threadIdx.x]+=a[i];
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
__global__ void max_do(float *restrict res, const float *a, const int n){
    extern __shared__ float sb[];
    sb[threadIdx.x]=0;
    int step=blockDim.x * gridDim.x ;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	if(sb[threadIdx.x]<a[i]) sb[threadIdx.x]=a[i];
    }
    for(step=(blockDim.x>>1);step>0;step>>=1){
	__syncthreads();
	if(threadIdx.x<step){
	    if(sb[threadIdx.x]<sb[threadIdx.x+step]){
		sb[threadIdx.x]=sb[threadIdx.x+step];
	    }
	}
    }
    if(threadIdx.x==0){
	atomicMax(res, sb[0]);
    }
}
/**
   tmp=sum(a.*b) res_add+=tmp

   2012-04-07: Bug found. The original implementation of res_rep does not work
   for multiple blocks where only the first block will be added to the final
   result.  */
__global__ void inn_do(float *res_add, const float *a, const float *b, const int n){
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
	atomicAdd(res_add, sb[0]);
    }
}
/* embed real to complex data.*/
__global__ void embed_do(fcomplex *out, float *in, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[ix]=make_cuComplex(in[ix], 0);
    }
}
/* extract real from complex data.*/
__global__ void extract_do(float *out, fcomplex *in, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[ix]=cuCrealf(in[ix]);
    }
}
__global__ void perm_f_do(fcomplex *restrict out, const fcomplex *restrict in, int *restrict perm, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[ix]=in[perm[ix]];
    }
}
__global__ void perm_i_do(fcomplex *restrict out, const fcomplex *restrict in, int *restrict perm, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[perm[ix]]=in[ix];
    }
}
__global__ void perm_f_do(float *restrict out, const float *restrict in, int *restrict perm, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[ix]=in[perm[ix]];
    }
}
__global__ void perm_i_do(float *restrict out, const float *restrict in, int *restrict perm, int nx){
    const int step=blockDim.x * gridDim.x;
    for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=step){
	out[perm[ix]]=in[ix];
    }
}

__global__ void embed_wvf_do(fcomplex *restrict wvf, 
					const float *restrict opd, const float *restrict amp, 
					const int *embed, const int nloc, const float wvl){
    const float pi2l=2.f*M_PI/wvl;
    for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nloc; ix+=blockDim.x*gridDim.x){
	float s,c;
	sincosf(pi2l*opd[ix], &s, &c);
	wvf[embed[ix]]=make_cuComplex(amp[ix]*c, amp[ix]*s);
    }
}

/**
   Embed or crop an array to another array. Preserve corner.
*/
__global__ void corner2center_do(fcomplex *restrict out, int noutx,  int nouty,
					    const fcomplex *restrict in, int ninx, int niny){
    int nx,ny;
    ny=MIN(niny, nouty)>>1;
    nx=MIN(ninx, noutx)>>1;
    int noutx2=noutx>>1;
    int nouty2=nouty>>1;
    for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
	for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
	    out[(iy+nouty2 )*noutx+(ix+noutx2)  ] = in[iy*ninx+ix];
	    out[(iy+nouty2 )*noutx+(noutx2-1-ix)] = in[iy*ninx+(ninx-1-ix)];
	    out[(nouty2-1-iy)*noutx+(noutx2-1-ix)] = in[(niny-1-iy)*ninx+(ninx-1-ix)];
	    out[(nouty2-1-iy)*noutx+(ix+noutx2)  ] = in[(niny-1-iy)*ninx+(ix)];
	}
    }
}

/**
   Embed or crop an array to another array. Preserve corner.
*/
__global__ void corner2center_abs2_do(float *restrict out, int noutx,  int nouty,
					    const fcomplex *restrict in, int ninx, int niny){
    int nx,ny;
    ny=MIN(niny, nouty)>>1;
    nx=MIN(ninx, noutx)>>1;
    int noutx2=noutx>>1;
    int nouty2=nouty>>1;
    for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
	for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
	    out[(iy+nouty2)*noutx+(ix+noutx2)]+=     CABS2(in[iy*ninx+ix]);
	    out[(iy+nouty2)*noutx+(noutx2-1-ix)]+=   CABS2(in[iy*ninx+(ninx-1-ix)]);
	    out[(nouty2-1-iy)*noutx+(noutx2-1-ix)]+= CABS2(in[(niny-1-iy)*ninx+(ninx-1-ix)]);
	    out[(nouty2-1-iy)*noutx+(ix+noutx2)]+=   CABS2(in[(niny-1-iy)*ninx+(ix)]);
	}
    }
}

/**
   Embed or crop an array to another array. Preserve corner.
*/
__global__ void corner2center_abs2_atomic_do(float *restrict out, int noutx,  int nouty,
					    const fcomplex *restrict in, int ninx, int niny){
    int nx,ny;
    ny=MIN(niny, nouty)>>1;
    nx=MIN(ninx, noutx)>>1;
    int noutx2=noutx>>1;
    int nouty2=nouty>>1;
    for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny; iy+=blockDim.y*gridDim.y){
	for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx; ix+=blockDim.x*gridDim.x){
	    atomicAdd(&out[(iy+nouty2)*noutx+(ix+noutx2)],     CABS2(in[iy*ninx+ix]));
	    atomicAdd(&out[(iy+nouty2)*noutx+(noutx2-1-ix)],   CABS2(in[iy*ninx+(ninx-1-ix)]));
	    atomicAdd(&out[(nouty2-1-iy)*noutx+(noutx2-1-ix)], CABS2(in[(niny-1-iy)*ninx+(ninx-1-ix)]));
	    atomicAdd(&out[(nouty2-1-iy)*noutx+(ix+noutx2)],   CABS2(in[(niny-1-iy)*ninx+(ix)]));
	}
    }
}
/**
   FFT Shift.
*/
__global__ void fftshift_do(fcomplex *wvf, const int nx, const int ny){
    int nx2=nx>>1;
    int ny2=ny>>1;
    for(int iy=threadIdx.y+blockDim.y*blockIdx.y; iy<ny2; iy+=blockDim.y*gridDim.y){
	for(int ix=threadIdx.x+blockDim.x*blockIdx.x; ix<nx2; ix+=blockDim.x*gridDim.x){
	    fcomplex tmp;
	    tmp=wvf[ix+iy*nx];
	    wvf[ix+iy*nx]=wvf[(ix+nx2)+(iy+ny2)*nx];
	    wvf[(ix+nx2)+(iy+ny2)*nx]=tmp;
	    tmp=wvf[ix+(iy+ny2)*nx];
	    wvf[ix+(iy+ny2)*nx]=wvf[(ix+nx2)+iy*nx];
	    wvf[(ix+nx2)+iy*nx]=tmp;
	}
    }
}
/**
   Add tip/tilt to the array. OPD=OPD+x*ttx+y*tty, where x=ix*dx+ox, y=iy*dy+oy; 
*/
__global__ void add_tilt_do(float *opd, int nx, int ny, float ox, float oy, float dx, float ttx, float tty){
    for(int iy=threadIdx.y; iy<ny; iy+=blockDim.y){
	float vty=(oy+iy*dx)*tty;
	for(int ix=threadIdx.x; ix<nx; ix+=blockDim.x){
	    opd[ix+iy*nx]+=vty+(ox+ix*dx)*ttx;
	}
    }
}

/**
   component wise multiply.
*/
__global__ void cwm_do(fcomplex *dest, float *from, int n){
    for(int i=threadIdx.x+threadIdx.y*blockDim.x; i<n; i+=blockDim.x*blockDim.y){
	dest[i].x*=from[i];
	dest[i].y*=from[i];
    }
}
/**
   unwrap the wvf back to opd. assume it within lambda/2 of the opd.
 */
__global__ void unwrap_phase_do(fcomplex *wvf, float *opd, int *embed, int n, float wvl){
    float kki=wvl/(2*M_PI);
    float wvlh=wvl*0.5;
    for(int i=threadIdx.x+threadIdx.y*blockDim.x; i<n; i+=blockDim.x*blockDim.y){
	float val=atan2(wvf[embed[i]].y, wvf[embed[i]].x)*kki;
	float diff=fmodf(val-opd[i]+wvlh, wvl);
	if(diff<0) diff+=wvl;
	opd[i]+=diff-wvlh;
    }
}
