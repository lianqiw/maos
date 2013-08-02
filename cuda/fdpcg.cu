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

#include "utils.h"
#include "tomo.h"
#include "recon.h"
#include "cucmat.h"
#define TIMING 0
namespace cuda_recon{
void cufdpcg_t::update(FDPCG_T *fdpcg){
    int nb=(fdpcg->nbx/2+1)*fdpcg->nby;//half frequency range
    //copy or update Mb. 
    int nxsave=fdpcg->Mbinv->nx;
    fdpcg->Mbinv->nx=nb;
    cp2gpu(&Mb, fdpcg->Mbinv);
    fdpcg->Mbinv->nx=nxsave;
}
cufdpcg_t::cufdpcg_t(FDPCG_T *fdpcg, curecon_geom *_grid)
    :grid(_grid),perm(0),Mb(0),fft(0),ffti(0),fftnc(0),fftips(0),
     xhat1(0),nby(0),nbz(0),scale(0),fddata(0){
    update(fdpcg);
    grid=_grid;
    scale=fdpcg->scale;
    int bs=fdpcg->bs;  
    int nb=(fdpcg->nbx/2+1)*fdpcg->nby;
    cp2gpu(&perm, fdpcg->permhf, nb*bs);
    int nps=grid->npsr;
    int count=0;
    int osi=-1;
    int start[nps];
    for(int ips=0; ips<nps; ips++){
	/*group layers with the same os together in a batch fft.*/
	if(osi != grid->xnx[ips]){
	    osi = grid->xnx[ips];
	    start[count]=ips;
	    count++;
	}
    }
    fft=(cufftHandle*)calloc(count, sizeof(cufftHandle));
    ffti=(cufftHandle*)calloc(count, sizeof(cufftHandle));
    fftnc=count;
    fftips=(int*)calloc(count+1, sizeof(int));
    for(int ic=0; ic<count; ic++){
	fftips[ic]=start[ic];
    }
    fftips[count]=nps;
    for(int ic=0; ic<count; ic++){
	int ncomp[2];
	/*Notice the reverse in specifying dimensions. THe first element is outmost rank.*/
	ncomp[0]=grid->xny[start[ic]];
	ncomp[1]=grid->xnx[start[ic]];

	int nembed[2];
	nembed[0]=grid->xnx[start[ic]]*grid->xny[start[ic]];
	nembed[1]=grid->xnx[start[ic]];
	DO(cufftPlanMany(&fft[ic], 2, ncomp, 
			 nembed, 1, ncomp[0]*ncomp[1], 
			 nembed, 1, ncomp[0]*ncomp[1], 
			 CUFFT_R2C, fftips[ic+1]-fftips[ic]));
	DO(cufftPlanMany(&ffti[ic], 2, ncomp, 
			 nembed, 1, ncomp[0]*ncomp[1], 
			 nembed, 1, ncomp[0]*ncomp[1],
			 CUFFT_C2R, fftips[ic+1]-fftips[ic]));
    }
    xhat1=cuccellnew(grid->npsr, 1, grid->xnx, grid->xny);
    {
	nby=256/bs;//number of blocks in each grid
	nbz=nb/nby;//number of grids to launch.
	while(nb!=nbz*nby){
	    nby--;
	    nbz=nb/nby;
	}
	nby=nby;
    }
    /* notice: performance may be improved by using
       R2C FFTs instead of C2C. Need to update perm
       and Mbinv to use R2C.*/
    GPU_FDPCG_T *FDDATA=new GPU_FDPCG_T[nps];
    for(int ips=0; ips<nps; ips++){
	FDDATA[ips].nx=grid->xnx[ips];
	FDDATA[ips].ny=grid->xny[ips];
	if(scale){
	    FDDATA[ips].scale=1.f/sqrtf((float)(grid->xnx[ips]*grid->xny[ips]));
	}else{
	    FDDATA[ips].scale=1.f;
	}
    }
    cudaMalloc(&fddata, sizeof(GPU_FDPCG_T)*nps);
    cudaMemcpy(fddata, FDDATA, sizeof(GPU_FDPCG_T)*nps, cudaMemcpyHostToDevice);
    delete [] FDDATA;
}
/*
  2012-08-01: Have tried the following with no improvement:
  Unroll with template.
  Ech 32 threads to handle each bs for bs<32.
*/

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

/**
   The positive (right half) base frequency (os=1) couples to both positive and
   negative frequencies to layers with os=2. Only apply the block matrix to
   positive frequencies maynot be good.

   If all lays have the same oversampling, either os=1 or os=2, then postive
   frequencies only couple with postive frequencies, so does negative
   frequencies. We can actually skip the multiplication of negative frequencies.

*/
void cufdpcg_t::P(curcell **xout, const curcell *xin, stream_t &stream){
#if TIMING
    EVENT_INIT(4)
#define RECORD(i) EVENT_TIC(i)
#else
#define RECORD(i)
#endif
    RECORD(0);
    if(!xin->m){
	error("xin is not continuous");
    }
    if(!*xout){
	*xout=curcellnew(grid->npsr, 1, grid->xnx, grid->xny);
    }else if(!(*xout)->m){
	error("xout is not continuous");
    }
    /*2012-07-11: Use real to complex fft*/
    for(int ic=0; ic<fftnc; ic++){
	int ips=fftips[ic];
	DO(cufftSetStream(fft[ic], stream));
	CUFFTR2C(fft[ic], xin->p[ips]->p, xhat1->p[ips]->p);
    }
    RECORD(1);
    if(scale){
	fdpcg_scale<<<dim3(3,3,grid->npsr), dim3(16,16),0,stream>>>
	    (fddata, xhat1->pm);
    }
    int bs=Mb->p[0]->nx;
    fdpcg_mul_block_sync_half<<<nbz, dim3(bs,nby), sizeof(fcomplex)*bs*2*nby, stream>>>
	(xhat1->m->p, Mb->m->p, perm);
    RECORD(2);
    for(int ic=0; ic<fftnc; ic++){
	int ips=fftips[ic];
	DO(cufftSetStream(ffti[ic], stream));
	CUFFTC2R(ffti[ic], xhat1->p[ips]->p, (*xout)->p[ips]->p);
    }
    if(scale){
	fdpcg_scale<<<dim3(9,1,grid->npsr),dim3(256,1),0,stream>>>
	    (fddata, (*xout)->pm);
    }
    RECORD(3);
#if TIMING
    EVENT_TOC;
    info2("FDPCG: FFT %3.0f MUL %3.0f FFTI %3.0f Total %3.0f\n", 
	  times[1], times[2], times[3], times[0]);
#endif
}
}//namespace
