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
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "wfs.h"
#include "recon.h"
#include "accphi.h"
#include "cucmat.h"

/*
  2012-08-01: Have tried the following with no improvement:
  Unroll with template.
  Ech 32 threads to handle each bs for bs<32.
*/

/* Each thread block is bs x 1. no need to lock.*/
/*
__global__ static void fdpcg_mul_block_sync(fcomplex *xin, fcomplex *M, int *restrict perm){
    int bs=blockDim.x;
    extern __shared__ fcomplex v[];
    fcomplex *vin=v+threadIdx.y*2*bs;
    fcomplex *vout=v+bs+threadIdx.y*2*bs;
    int ib=blockIdx.x*blockDim.y+threadIdx.y;
    perm+=ib*bs;
    M+=ib*bs*bs;
    int ix=threadIdx.x;
    vin[ix]=xin[perm[ix]];
    vout[ix]=make_cuComplex(0,0);
    __syncthreads();//wait for loading of vin
    for(int iy=0; iy<bs; iy++){
	vout[ix]=cuCfmaf(M[ix+iy*bs], vin[iy], vout[ix]);
    }
    xin[perm[ix]]=vout[ix];
    }*/
__global__ static void fdpcg_mul_block_sync_half(fcomplex *xin, fcomplex *M, int *restrict perm){
    int bs=blockDim.x;
    extern __shared__ fcomplex v[];
    fcomplex *vin=v+threadIdx.y*2*bs;
    fcomplex *vout=v+bs+threadIdx.y*2*bs;
    int ib=blockIdx.x*blockDim.y+threadIdx.y;
    M+=ib*bs*bs;
    int ix=threadIdx.x;
    int pm=perm[ib*bs+ix];
    vin[ix]=xin[abs(pm)];
    if(pm<0){
	vin[ix]=cuConjf(vin[ix]);
    }
    vout[ix]=make_cuComplex(0,0);
    __syncthreads();//wait for loading of vin
    for(int iy=0; iy<bs; iy++){
	vout[ix]=cuCfmaf(M[ix+iy*bs], vin[iy], vout[ix]);
    }
    if(pm<0){
	xin[-pm]=cuConjf(vout[ix]);
    }else{
	xin[pm]=vout[ix];
    }

}
__global__ static void fdpcg_scale(GPU_FDPCG_T *fddata, float **xall){
    int ips=blockIdx.z;
    int nx=fddata[ips].nx*fddata[ips].ny;
    int step=blockDim.x * gridDim.x; 
    float scale=fddata[ips].scale;
    float *restrict x=xall[ips];
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nx; i+=step){
	x[i]*=scale;
    }
}
__global__ static void fdpcg_scale(GPU_FDPCG_T *fddata, fcomplex **xall){
    int ips=blockIdx.z;
    int nx=fddata[ips].nx*fddata[ips].ny;
    int step=blockDim.x * gridDim.x; 
    float scale=fddata[ips].scale;
    fcomplex *restrict x=xall[ips];
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<nx; i+=step){
	cuCscalef(x[i], scale);
    }
}

/*
  We used FFT R2C for transforming from real to half part of complex array. Need
  to fill the other half of array. Scale the array in the mean time.*/
/*__global__ static void fdpcg_r2c(GPU_FDPCG_T *fddata, fcomplex **outall){
    const int ips=blockIdx.z;
    const int nx=fddata[ips].nx;
    const int ny=fddata[ips].ny;
    const int nx2=(nx>>1);
    const float scale=fddata[ips].scale;
    fcomplex *restrict out=outall[ips];
    const int ix0=threadIdx.x+blockDim.x*blockIdx.x;
    const int iy0=threadIdx.y+blockDim.y*blockIdx.y;
    for(int iy=iy0; iy<ny; iy+=blockDim.y*gridDim.y){
	if(ix0<2){
	    int ix=nx2*ix0;//0 and nx2 are real.
	    out[ix+iy*nx]=make_cuComplex(cuCrealf(out[ix+iy*nx])*scale, cuCimagf(out[ix+iy*nx])*scale);
	}
	for(int ix=1+ix0; ix<nx2; ix+=blockDim.x*gridDim.x){
	    float xr=cuCrealf(out[ix+iy*nx])*scale;
	    float xi=cuCimagf(out[ix+iy*nx])*scale;
	    out[ix+iy*nx]=make_cuComplex(xr, xi);
	    out[(nx-ix)+(iy==0?0:ny-iy)*nx]=make_cuComplex(xr, -xi);
	}
    }
    }*/

/**
   The positive (right half) base frequency (os=1) couples to both positive and
   negative frequencies to layers with os=2. Only apply the block matrix to
   positive frequencies maynot be good.

   If all lays have the same oversampling, either os=1 or os=2, then postive
   frequencies only couple with postive frequencies, so does negative
   frequencies. We can actually skip the multiplication of negative frequencies.

*/
#define TIMING 0
void gpu_Tomo_fdprecond(curcell **xout, const void *A, const curcell *xin, cudaStream_t stream){
#if TIMING==2
#define NEVENT 10
    static cudaEvent_t event[NEVENT]={0};
    static float times[NEVENT];
    if(!event[0]){
	for(int i=0; i<NEVENT; i++){
	    DO(cudaEventCreate(&event[i]));
	}
    }
#define RECORD(i) DO(cudaEventRecord(event[i], stream))
#else
#define RECORD(i)
#endif
    RECORD(0);
    curecon_t *curecon=cudata->recon;
    const RECON_T *recon=(const RECON_T *)A;
    if(!xin->m){
	error("xin is not continuous");
    }
    if(!*xout){
	*xout=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
    }else if(!(*xout)->m){
	error("xout is not continuous");
    }
    /*2012-07-11: Use real to complex fft*/
    cufdpcg_t *cufd=curecon->fdpcg;
    for(int ic=0; ic<cufd->fftnc; ic++){
	int ips=cufd->fftips[ic];
	CUFFTR2C(cufd->fft[ic], xin->p[ips]->p, cufd->xhat1->p[ips]->p);
    }
    RECORD(1);
    if(cufd->scale){
	fdpcg_scale<<<dim3(3,3,recon->npsr), dim3(16,16),0,stream>>>
	    (curecon->fddata, cufd->xhat1->pm);
    }
    RECORD(2);
    int bs=cufd->Mb->p[0]->nx;
    fdpcg_mul_block_sync_half<<<cufd->nbz, dim3(bs,cufd->nby), sizeof(fcomplex)*bs*2*cufd->nby, stream>>>
	(cufd->xhat1->m->p, cufd->Mb->m->p, cufd->perm);
    /*
      fdpcg_mul_block_sync<<<cufd->nbz, dim3(bs,cufd->nby), sizeof(fcomplex)*bs*2*cufd->nby, stream>>>
      (cufd->xhat1->m->p, cufd->Mb->m->p, cufd->perm);
    */
    RECORD(3);
    for(int ic=0; ic<cufd->fftnc; ic++){
	int ips=cufd->fftips[ic];
	CUFFTC2R(cufd->ffti[ic], cufd->xhat1->p[ips]->p, (*xout)->p[ips]->p);
    }
    RECORD(4);
    if(cufd->scale){
	fdpcg_scale<<<dim3(9,1,recon->npsr),dim3(256,1),0,stream>>>
	    (curecon->fddata, (*xout)->pm);
    }
    RECORD(5);

#if TIMING==2
    CUDA_SYNC_STREAM;
    for(int i=1; i<6; i++){
	DO(cudaEventElapsedTime(&times[i], event[i-1], event[i]));
	times[i]*=1e3;
    }
    info2("FDPCG: FFT %3.0f SC %3.0f MUL %3.0f FFTI %3.0f SC %3.0f \n", 
	  times[1], times[2], times[3], times[4], times[5]);
#endif
}
