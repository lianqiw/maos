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
#include "utils.h"
#include "curmat.h"
#include "cucmat.h"

/**
\file mvm_iwfs.cu

Test MVM for a single WFS. Using two GPUs, without networking.
This is not part of maos executable. Called by test/test_gpu executable.


two ways: 
1) broadcast gradients to both GPUs with each GPU handle part of actuators. (not tested yet)
2) partition gradients to GPUs with GPU handle all actuators. (selected)

Use schedtool -a 0x1 PID to let the exe only run one CPU 0. This prevents jitter.
*/

#define TIMING 0

#if TIMING 
unsigned int event_flag=cudaEventDefault;
#else
unsigned int event_flag=cudaEventDisableTiming;
#endif
typedef struct{
    curmat *cumvm;//active mvm control matrix
    curmat *cumvm_next;//inactive mvm control matrix.
    curmat *cumvm1;
    curmat *cumvm2;
    curmat *mtch;
    curmat *pix;//pixels. Each sa has 15x6=90 pixels.
    curmat *grad;
    curmat *act;
    stream_t *stream_p;//pixels
    stream_t *stream_g;//grads
    stream_t *stream_a;//act
    stream_t *stream_mvm;//mvm
    int count;
    int gpu;//Which GPU this data is for
    int istep;//Which time step we are in
    int copy_mvm;//1: need to copy mvm.
    int ic;//the column that we are copying.
    cudaEvent_t *event_p;
    cudaEvent_t *event_g;
    cudaEvent_t event_pall;
#if TIMING
    cudaEvent_t event0;
    cudaEvent_t *event0_p;
    cudaEvent_t *event0_g;
    cudaEvent_t event0_a;
    cudaEvent_t event_a;
    cudaEvent_t *event0_a2;
    cudaEvent_t *event_a2;
    cudaEvent_t event0_mvm;
    cudaEvent_t event_mvm;
#endif
}GPU_DATA_T;
/*Does matched filter*/
static void __global__ mtch_do(const float *mtch, const float *pix, float *grad, int pixpsa, int nsa){
    extern __shared__ float cum[];//for cumulation and reduction
    float *cumi=cum+threadIdx.y*blockDim.x;//2 padding for easy reduction
    int ig=threadIdx.y+blockDim.y*blockIdx.x;
    const float *mtchi=mtch+ig*pixpsa;
    const float *pixi=pix+ig/2*pixpsa;
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
void mvmfull_iwfs(char *fnmvm1, char *fnmvm2, char *fnpix1, char *fnpix2, char *fnmtch, 
	      int *gpus, int ngpu, int nstep){
    info("Using %d gpus\n", ngpu);
    warning2("Notice that here we group x/y gradients together like xyxyxy instead of like"
	    "xxxyyy in matched filter and MVM here.\n");
    smat *mvm1=sread("%s", fnmvm1);
    smat *mvm2=sread("%s", fnmvm2);
    smat *mvm=mvm1;
    
    /*smat *grad1=sread("%s", fngrad1);
      smat *grad2=sread("%s", fngrad2);*/
    smat *pix1=sread("%s", fnpix1);
    smat *pix2=sread("%s", fnpix2);
    /*Important: 
      1) Only page locked host memory can do async memcpy that overallps with computation
      2) Has to be Portable for multiple GPUs to do async memcpy concurrently.
     */
    smat *pix=pix2;
    smat *mtch=sread("%s", fnmtch);
    const int nsa=pix1->ny;
    const int ng=nsa*2;
    /*reduce the matrices to only a single wfs.*/
    sresize(mvm1, mvm1->nx, ng);
    sresize(mvm2, mvm2->nx, ng);
    scell *dmres=scellnew(ngpu, 1);
    spagelock(pix1, pix2, mvm1, mvm2, mtch, NULL);
    const int pixpsa=90;//Change this need to change kernel mtch_do
    const int mtch_ngrid=30;//30;//can change to utilize GPU fully. 16 is good for cassiopeia
    const int mtch_dimx=32;//must launch 32 threads so that they belong to single wrap. use only 30 threads.
    const int mtch_dimy=12;//4 subapertures, 8 gradients
    const int sastep=mtch_dimy*mtch_ngrid/2;
    const int nact=mvm1->nx;
    int nc=10;//each time copy nc column of mvm.
    GPU_DATA_T *data=(GPU_DATA_T*)calloc(ngpu, sizeof(GPU_DATA_T));
    const int sect_gpu=(nsa+sastep*ngpu-1)/(sastep*ngpu);
    for(int igpu=0; igpu<ngpu; igpu++){
	cudaSetDevice(gpus[igpu]);
	data[igpu].cumvm1=curnew(mvm1->nx, ng);
	data[igpu].cumvm2=curnew(mvm2->nx, ng);
	data[igpu].cumvm=data[igpu].cumvm1;
	data[igpu].cumvm_next=data[igpu].cumvm2;
	cp2gpu(&data[igpu].cumvm1, mvm);
	data[igpu].pix=curnew(pixpsa, nsa);
	data[igpu].mtch=curnew(pixpsa, nsa*2);
	cp2gpu(&data[igpu].mtch, mtch);
	data[igpu].grad=curnew(ng, 1);
	data[igpu].act=curnew(mvm1->nx, 1);
	data[igpu].stream_p=new stream_t;
	data[igpu].stream_g=new stream_t;
	data[igpu].stream_a=new stream_t;
	data[igpu].stream_mvm=new stream_t;
	data[igpu].gpu=gpus[igpu];
#if TIMING
	cudaEventCreateWithFlags(&data[igpu].event0, event_flag);
	data[igpu].event0_g=new cudaEvent_t[sect_gpu];
	data[igpu].event0_p=new cudaEvent_t[sect_gpu];
	data[igpu].event0_a2=new cudaEvent_t[sect_gpu];
	data[igpu].event_a2=new cudaEvent_t[sect_gpu];

	for(int i=0; i<sect_gpu; i++){
	    cudaEventCreateWithFlags(&data[igpu].event0_g[i],event_flag);
	    cudaEventCreateWithFlags(&data[igpu].event0_p[i],event_flag);
	    cudaEventCreateWithFlags(&data[igpu].event0_a2[i],event_flag);
	    cudaEventCreateWithFlags(&data[igpu].event_a2[i],event_flag);
	}
	cudaEventCreateWithFlags(&data[igpu].event0_mvm,event_flag);
	cudaEventCreateWithFlags(&data[igpu].event_mvm,event_flag);
	cudaEventCreateWithFlags(&data[igpu].event0_a,event_flag);
	cudaEventCreateWithFlags(&data[igpu].event_a,event_flag);
#endif
	data[igpu].event_g=new cudaEvent_t[sect_gpu];
	data[igpu].event_p=new cudaEvent_t[sect_gpu];
	for(int i=0; i<sect_gpu; i++){
	    cudaEventCreateWithFlags(&data[igpu].event_g[i],event_flag);
	    cudaEventCreateWithFlags(&data[igpu].event_p[i],event_flag);
	}
	cudaEventCreateWithFlags(&data[igpu].event_pall,event_flag);
	dmres->p[igpu]=snew(nact, 1);
	spagelock(dmres->p[igpu], NULL);
    }
    smat *timing=snew(nstep, 1);
    smat *timing2=snew(nstep, 1);
    smat *result=snew(nstep, 1);
    float one=1; float zero=0; float *pbeta;
    cudaProfilerStart();
    TIC;tic;
    for(int istep=0; istep<nstep; istep++){
#if TIMING
	if(istep%8000==0)
#else
	    //if(istep%8000==7484)
	    if(0)
#endif
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
	    data[igpu].count=0;
	    data[igpu].istep=istep;
#if TIMING
	    //beginning of each GPU operation.
	    DO(cudaEventRecord(data[igpu].event0, data[igpu].stream_a[0]));
#endif
	}
	if(pix==pix1){
	    pix=pix2;
	}else{
	    pix=pix1;
	}

	for(int isa=0, igpu=0; isa<nsa; isa+=sastep, igpu=((igpu+1)%ngpu)){
	    cudaSetDevice(gpus[igpu]); 
	    GPU_DATA_T *datai=&data[igpu];
	    int nleft=(nsa-isa)<sastep?(nsa-isa):sastep;
	    //One stream handling the memcpy
#if TIMING
	    DO(cudaEventRecord(datai->event0_p[datai->count], datai->stream_p[0]));
#endif
	    DO(cudaMemcpyAsync(datai->pix->p+isa*pixpsa, pix->p+isa*pixpsa, sizeof(float)*nleft*pixpsa,
			       cudaMemcpyHostToDevice, *datai->stream_p));
	    //Recored the event when the memcpy is finished
	    DO(cudaEventRecord(datai->event_p[datai->count], datai->stream_p[0]));
	    //Start matched filter when pixel transfer is done.
	    DO(cudaStreamWaitEvent(datai->stream_g[0], datai->event_p[datai->count], 0));
#if TIMING
	    DO(cudaEventRecord(datai->event0_g[datai->count], datai->stream_g[0]));    
#endif
	    mtch_do<<<mtch_ngrid, dim3(mtch_dimx, mtch_dimy), 
		mtch_dimx*mtch_dimy*sizeof(float), datai->stream_g[0]>>>
	       (datai->mtch->p+isa*2*pixpsa, datai->pix->p+isa*pixpsa, 
		datai->grad->p+isa*2, pixpsa, nleft);
	    //Record the event when matched filter is done
	    DO(cudaEventRecord(datai->event_g[datai->count], datai->stream_g[0]));
	    if(!datai->count){
		pbeta=&zero;//initialize act.
	    }else{
		pbeta=&one;
	    }
	    //Another stream does the matrix vector multiplication. Wait for the event before executing.
	    //The stream stream will wait only for the completion of the most recent host call to cudaEventRecord() on event
	    cudaStreamWaitEvent(datai->stream_a[0], datai->event_g[datai->count], 0);
#if TIMING
	    DO(cudaEventRecord(datai->event0_a2[datai->count], datai->stream_a[0]));    
#endif
	    DO(cublasSgemv(datai->stream_a[0], CUBLAS_OP_N, nact, nleft*2, &one, datai->cumvm->p+nact*isa*2, nact, datai->grad->p+isa*2, 1, pbeta, datai->act->p, 1));
#if TIMING
	    DO(cudaEventRecord(datai->event_a2[datai->count], datai->stream_a[0])); 
#endif
	    datai->count++;
	}
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=&data[igpu];
	    //Record an event when pixel tranporting is over. So we can start transporting mvm matrix.
	    DO(cudaEventRecord(datai->event_pall, datai->stream_p[0]));
	}
	//Queue copying MVM matrix to second slot.
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=&data[igpu];
	    if(datai->copy_mvm){
		int done=0, nleft;
		if(mvm->ny-datai->ic < nc){
		    done=1;
		    nleft=mvm->ny-datai->ic;
		}else{
		    nleft=nc;
		}
		//wait for mvm application to finish before copying.
		//DO(cudaStreamWaitEvent(datai->stream_mvm[0], datai->event_pall, 0));
#if TIMING
		DO(cudaEventRecord(datai->event0_mvm, datai->stream_mvm[0]));	
#endif
		DO(cudaMemcpyAsync(datai->cumvm_next->p+datai->ic*mvm->nx, 
				   mvm->p+datai->ic*mvm->nx, sizeof(float)*mvm->nx*nleft, 
				   cudaMemcpyHostToDevice, datai->stream_mvm[0]));
#if TIMING
		DO(cudaEventRecord(datai->event_mvm, datai->stream_mvm[0]));
#endif
		datai->ic+=nleft;
		if(done){
		    datai->ic=0;
		    datai->copy_mvm=0;
		    curmat *tmp=datai->cumvm;
		    datai->cumvm=datai->cumvm_next;
		    datai->cumvm_next=tmp;
		    info2("gpu %d switched over at step %d\n", datai->gpu, datai->istep);
		}
	    }
	}
	//Copy DM commands back to CPU
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=&data[igpu];
	    cudaSetDevice(gpus[igpu]); 
#if TIMING
	    DO(cudaEventRecord(datai->event0_a, datai->stream_a[0]));
#endif
	    cudaMemcpyAsync(dmres->p[igpu]->p, datai->act->p, nact*sizeof(float), cudaMemcpyDeviceToHost, datai->stream_a[0]);
#if TIMING
	    DO(cudaEventRecord(datai->event_a, datai->stream_a[0]));//record event when all act are copied so mvm can start.
#endif
	}
	data[0].stream_a->sync();
	//CPU sums them together
	for(int igpu=1; igpu<ngpu; igpu++){
	    cudaSetDevice(gpus[igpu]); 
	    data[igpu].stream_a->sync();
	    for(int iact=0; iact<nact; iact++){
		dmres->p[0]->p[iact]+=dmres->p[igpu]->p[iact];
	    }
	}
	result->p[istep]=dmres->p[0]->p[nact/2];
	timing->p[istep]=toc3;//do not tic.
	if(istep%1000==0 || timing->p[istep]>2.e-3){
	    info2("Step %d takes %.0f us\n", istep, timing->p[istep]*1e6);
	}

	//Wait for MVM matrix copy to finish and time.
	for(int igpu=0; igpu<ngpu; igpu++){
	    cudaSetDevice(data[igpu].gpu);
	    data[igpu].stream_mvm->sync();
	}
	timing2->p[istep]=toc3;tic;
#if TIMING 
	if(istep<100){
	    for(int igpu=0; igpu<ngpu; igpu++){
		cudaSetDevice(gpus[igpu]); 
		GPU_DATA_T *datai=data+igpu;
		const int count=datai->count;
		smat *tim=snew(count*6+4,2);
		PSMAT(tim,ptim);
		int ic;
		for(ic=0; ic<count; ic++){
		    cudaEventElapsedTime(&ptim[0][ic*6+0], datai->event0, datai->event0_p[ic]);//start of mtch
		    cudaEventElapsedTime(&ptim[0][ic*6+1], datai->event0, datai->event_p[ic]);//end of mtch
		    cudaEventElapsedTime(&ptim[0][ic*6+2], datai->event0, datai->event0_g[ic]);//start of g
		    cudaEventElapsedTime(&ptim[0][ic*6+3], datai->event0, datai->event_g[ic]);//end of g
		    cudaEventElapsedTime(&ptim[0][ic*6+4], datai->event0, datai->event0_a2[ic]);//start of mvm
		    cudaEventElapsedTime(&ptim[0][ic*6+5], datai->event0, datai->event_a2[ic]);//end of mvm
		    ptim[1][ic*6]=1;
		    ptim[1][ic*6+1]=1;
		    ptim[1][ic*6+2]=2;
		    ptim[1][ic*6+3]=2;
		    ptim[1][ic*6+4]=3;
		    ptim[1][ic*6+5]=3;
		}
		cudaEventElapsedTime(&ptim[0][ic*6+0], datai->event0, datai->event0_a);//start of a copy
		cudaEventElapsedTime(&ptim[0][ic*6+1], datai->event0, datai->event_a);//end of a copy
		cudaEventElapsedTime(&ptim[0][ic*6+2], datai->event0, datai->event0_mvm);//start of mvm copy
		cudaEventElapsedTime(&ptim[0][ic*6+3], datai->event0, datai->event_mvm);//end of mvm copy
		ptim[1][ic*6+0]=4;
		ptim[1][ic*6+1]=4;
		ptim[1][ic*6+2]=5;
		ptim[1][ic*6+3]=5;
		swrite(tim, "timing2_%dgpu%d_step%d", ngpu, igpu, istep);
		sfree(tim);
	    }
	}
	tic;
#endif
    }
    cudaProfilerStop();
    swrite(timing, "timing_%dgpu", ngpu);
    swrite(timing2, "timing2_%dgpu", ngpu);
    swrite(result, "result_%dgpu", ngpu);
    spageunlock(pix1, pix2, mvm1, mvm2, NULL);
}
