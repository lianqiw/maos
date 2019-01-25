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
#include "curmat.h"
#include "cucmat.h"

/**
   \file mvmfull_pipe.cu

Copied from mvmfull_iwfs.cu to test another scheme for parallelizing. More
streams are used now with each stream doing tasks in series. This helps to increase the occupancy on Kepler GPUs with more lower speed cores.

*/

#define TIMING 0
typedef struct{
    curmat cumvm;//active mvm control matrix
    curmat cumvm_next;//inactive mvm control matrix.
    curmat cumvm1;
    curmat cumvm2;
    curmat mtch;
    curmat pix;//pixels. Each sa has 15x6=90 pixels.
    curmat grad;
    curmat act;
    Array<stream_t> stream;
    stream_t stream_mvm;//mvm
    int ism;//index of stream for mvm
    int count;
    int gpu;//Which GPU this data is for
    int istep;//Which time step we are in
    int copy_mvm;//1: need to copy mvm.
    int ic;//the column that we are copying.
    event_t *event_p;
    event_t *event_g;
    event_t *event_pall;

}GPU_DATA_T;
/*Does matched filter*/
static void __global__ mtch_do(const Real *mtch, const Real *pix, Real *grad, int pixpsa, int nsa){
    extern __shared__ Real cum[];//for cumulation and reduction
    Real *cumi=cum+threadIdx.y*blockDim.x;//2 padding for easy reduction
    int ig=threadIdx.y+blockDim.y*blockIdx.x;
    const Real *mtchi=mtch+ig*pixpsa;
    const Real *pixi=pix+ig/2*pixpsa;
    if(ig>nsa*2) return;//over range
    //sum 3 times for 90 pixels.
    cumi[threadIdx.x]=0;
    if(threadIdx.x<30){
	cumi[threadIdx.x]=mtchi[threadIdx.x]*pixi[threadIdx.x]
	    +mtchi[threadIdx.x+30]*pixi[threadIdx.x+30]
	    +mtchi[threadIdx.x+60]*pixi[threadIdx.x+60];
    }
    //reduction
    for(int step=16;step>0;step>>=1){
	if(threadIdx.x<step){
	    cumi[threadIdx.x]+=cumi[threadIdx.x+step];
	}
    }
    if(threadIdx.x==0){
	grad[ig]=cumi[0];
    }
}

/**
   A standalone routine that testes applying MVM for a single WFS and update mvm.*/
void mvmfull_pipe(char *fnmvm1, char *fnmvm2, char *fnpix1, char *fnpix2, char *fnmtch, 
	      int *gpus, int ngpu, int nstep){
    {
	DO(cudaFuncSetCacheConfig(mvm_do, cudaFuncCachePreferShared));
	struct cudaDeviceProp prop;
	DO(cudaGetDeviceProperties(&prop, 0));
    }
    const int nsm=6;

    dbg("Using %d gpus\n", ngpu);
    //warning("Notice that here we group x/y gradients together like xyxyxy instead of like"
    //	    "xxxyyy in matched filter and MVM here.\n");
    X(mat) *mvm1=X(read)("%s", fnmvm1);
    X(mat) *mvm2=X(read)("%s", fnmvm2);
    X(mat) *mvm=mvm1;
    
    /*X(mat) *grad1=X(read)("%s", fngrad1);
      X(mat) *grad2=X(read)("%s", fngrad2);*/
    X(mat) *pix1=X(read)("%s", fnpix1);
    X(mat) *pix2=X(read)("%s", fnpix2);
    /*Important: 
      1) Only page locked host memory can do async memcpy that overallps with computation
      2) Has to be Portable for multiple GPUs to do async memcpy concurrently.
     */
    X(mat) *pix=pix2;
    X(mat) *mtch=X(read)("%s", fnmtch);
    const int nsa=pix1->ny;
    const int ng=nsa*2;
    /*reduce the matrices to only a single wfs.*/
    X(resize)(mvm1, mvm1->nx, ng);
    X(resize)(mvm2, mvm2->nx, ng);
    X(cell) *dmres=X(cellnew)(ngpu, 1);
    X(pagelock)(pix1, pix2, mvm1, mvm2, mtch, NULL);
    const int pixpsa=90;//Change this need to change kernel mtch_do
    const int mtch_ngrid=40;//30;//can change to utilize GPU fully. 16 is good for cassiopeia
    const int mtch_dimx=32;//must launch 32 threads so that they belong to single wrap. use only 30 threads.
    const int mtch_dimy=12;//4 subapertures, 8 gradients
    const int sastep=mtch_dimy*mtch_ngrid/2;
    const int nact=mvm1->nx;

    GPU_DATA_T *data=(GPU_DATA_T*)calloc(ngpu, sizeof(GPU_DATA_T));
    const int sect_gpu=(nsa+sastep*ngpu-1)/(sastep*ngpu);
    for(int igpu=0; igpu<ngpu; igpu++){
	cudaSetDevice(gpus[igpu]);
	data[igpu].cumvm1=curmat(mvm1->nx, ng);
	data[igpu].cumvm2=curmat(mvm2->nx, ng);
	data[igpu].cumvm=data[igpu].cumvm1;
	data[igpu].cumvm_next=data[igpu].cumvm2;
	cp2gpu(data[igpu].cumvm1, mvm);
	data[igpu].pix=curmat(pixpsa, nsa);
	data[igpu].mtch=curmat(pixpsa, nsa*2);
	cp2gpu(data[igpu].mtch, mtch);
	data[igpu].grad=curmat(ng, 1);
	data[igpu].act=curmat(mvm1->nx, 1);
	data[igpu].stream=Array<stream_t>(nsm, 1);
	data[igpu].gpu=gpus[igpu];

	data[igpu].event_g=new event_t[sect_gpu];
	data[igpu].event_p=new event_t[sect_gpu];
	data[igpu].event_pall=new event_t;
	dmres->p[igpu]=X(new)(nact, 1);
	X(pagelock)(dmres->p[igpu], NULL);
    }
    X(mat) *timing=X(new)(nstep, 1);
    X(mat) *result=X(new)(nstep, 1);
    cudaProfilerStart();
    TIC;tic;
    for(int istep=0; istep<nstep; istep++){
	//if(istep%8000==7484)
	if(0)
	    {//need to update MVM
		if(mvm==mvm1){//switch mvm on host.
		    mvm=mvm2;
		}else{
		    mvm=mvm1;
		}
		for(int igpu=0; igpu<ngpu; igpu++){
		    data[igpu].copy_mvm=1;
		    if(data[igpu].ic!=0){
			warning("Sync error, skip update request at step %d\n", istep);
		    }
		}
	    }
	for(int igpu=0; igpu<ngpu; igpu++){
	    data[igpu].ism=-1;
	    data[igpu].count=0;
	    data[igpu].istep=istep;
	}
	if(pix==pix1){
	    pix=pix2;
	}else{
	    pix=pix1;
	}

	for(int isa=0, igpu=0; isa<nsa; isa+=sastep, igpu=((igpu+1)%ngpu)){
	    cudaSetDevice(gpus[igpu]); 
	    GPU_DATA_T *datai=&data[igpu];
	    int ism=datai->ism=(datai->ism+1)%nsm;
	    int nleft=(nsa-isa)<sastep?(nsa-isa):sastep;

	    DO(cudaMemcpyAsync(datai->pix()+isa*pixpsa, pix->p+isa*pixpsa, sizeof(Real)*nleft*pixpsa,
			       cudaMemcpyHostToDevice, datai->stream[ism]));
	    //Start matched filter in the same stream
	    mtch_do<<<mtch_ngrid, dim3(mtch_dimx, mtch_dimy), 
		mtch_dimx*mtch_dimy*sizeof(Real), datai->stream[ism]>>>
		(datai->mtch()+isa*2*pixpsa, datai->pix()+isa*pixpsa, 
		 datai->grad()+isa*2, pixpsa, nleft);
#if 0
	    Real one=1; Real zero=0; 
	    Real *pbeta;
	    if(!datai->count){
		pbeta=&zero;//initialize act.
	    }else{
		pbeta=&one;
	    }
	    DO(CUBL(gemv)(datai->stream[ism], CUBLAS_OP_N, nact, nleft*2, 
			  &one, datai->cumvm->p+nact*isa*2, nact, datai->grad->p+isa*2, 
			  1, pbeta, datai->act->p, 1));
#else
	    {
		const int naeach=128;
		const int nstream=10;
		const int nblock=(nact*nstream+naeach-1)/naeach;
		multimv_do<<<nblock, naeach, sizeof(Real)*naeach, datai->stream[ism]>>>
		    (datai->cumvm()+nact*isa*2, datai->act, datai->grad()+isa*2, 
		     nact, nleft*2);
	    }
#endif
	    datai->count++;
	}

	//Copy DM commands back to CPU
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=&data[igpu];
	    cudaSetDevice(gpus[igpu]); 
	    for(int ism=0; ism<nsm; ism++){
		datai->stream[ism].sync();
	    }
	    cudaMemcpyAsync(dmres->p[igpu]->p, datai->act, nact*sizeof(Real), 
			    cudaMemcpyDeviceToHost, datai->stream[0]);

	}
	//CPU sums them together
	for(int igpu=1; igpu<ngpu; igpu++){
	    cudaSetDevice(gpus[igpu]); 
	    data[igpu].stream[0].sync();
	    for(int iact=0; iact<nact; iact++){
		dmres->p[0]->p[iact]+=dmres->p[igpu]->p[iact];
	    }
	}
	result->p[istep]=dmres->p[0]->p[nact/2];
	timing->p[istep]=toc3;//do not tic.
	if(istep%1000==0 || timing->p[istep]>2.e-3){
	    info("Step %d takes %.0f us\n", istep, timing->p[istep]*1e6);
	}
	tic;
    }
    cudaProfilerStop();
    writebin(timing, "timing_%dgpu", ngpu);
    writebin(result, "result_%dgpu", ngpu);
    X(pageunlock)(pix1, pix2, mvm1, mvm2, NULL);
}
