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
This is not part of maos executable.

two ways: 
1) broadcast gradients to both GPUs with each GPU handle part of actuators
2) partition gradients to GPUs with GPU handle all actuators.
*/
typedef struct{
    curmat *cumvm;//active mvm control matrix
    curmat *cumvm_next;//inactive mvm control matrix.
    curmat *cumvm1;
    curmat *cumvm2;
    curmat *grad;
    curmat *act;
    stream_t *stream_g;
    stream_t *stream_a;
    stream_t *stream_mvm;
    int count_g;
    int gpu;//Which GPU this data is for
    int istep;//Which time step we are in
    int copy_mvm;//1: need to copy mvm.
    int ic;//the column that we are copying.
    cudaEvent_t *event_g;
    cudaEvent_t event_a;
    cudaEvent_t event_mvm;
}GPU_DATA_T;
/**
   A standalone routine that testes applying MVM for a single WFS and update mvm.*/
void mvm_iwfs(char *fnmvm1, char *fnmvm2, char *fngrad1, char *fngrad2, int *gpus, int ngpu, int nstep){
    smat *mvm1=sread("%s", fnmvm1);
    smat *mvm2=sread("%s", fnmvm2);
    smat *mvm=mvm1;

    smat *grad1=sread("%s", fngrad1);
    smat *grad2=sread("%s", fngrad2);
    smat *grad=grad1;
    scell *dmres=scellnew(ngpu, 1);
    const int nact=mvm1->nx;
    const int ng=grad->nx;
    /*reduce the matrices to only a single wfs.*/
    sresize(mvm1, mvm1->nx, ng);
    sresize(mvm2, mvm2->nx, ng);
    int nc=10;//each time copy nc column of mvm.
    GPU_DATA_T *data=(GPU_DATA_T*)calloc(ngpu, sizeof(GPU_DATA_T));
    const int gstep=500;
    const int sect_gpu=(ng+gstep*ngpu-1)/(gstep*ngpu);
    for(int igpu=0; igpu<ngpu; igpu++){
	cudaSetDevice(gpus[igpu]);
	data[igpu].cumvm1=curnew(mvm1->nx, ng);
	data[igpu].cumvm2=curnew(mvm2->nx, ng);
	data[igpu].cumvm=data[igpu].cumvm1;
	data[igpu].cumvm_next=data[igpu].cumvm2;
	cp2gpu(&data[igpu].cumvm1, mvm);
	data[igpu].grad=curnew(ng, 1);
	data[igpu].act=curnew(mvm1->nx, 1);
	data[igpu].stream_g=new stream_t;
	data[igpu].stream_a=new stream_t;
	data[igpu].stream_mvm=new stream_t;
	data[igpu].event_g=new cudaEvent_t[sect_gpu];
	data[igpu].gpu=gpus[igpu];
	for(int i=0; i<sect_gpu; i++){
	    cudaEventCreateWithFlags(&data[igpu].event_g[i],cudaEventDisableTiming);
	}
	cudaEventCreateWithFlags(&data[igpu].event_mvm,cudaEventDisableTiming);
	cudaEventCreateWithFlags(&data[igpu].event_a,cudaEventDisableTiming);

	dmres->p[igpu]=snew(nact, 1);
	//pthread_create(&threads[igpu], NULL, (thread_fun)mvm_cp, &data[igpu]);
    }
    smat *timing=snew(nstep, 1);TIC;tic;
    smat *result=snew(nstep, 1);
    float one=1; float zero=0; float *pbeta;
    for(int istep=0; istep<nstep; istep++){
	if(istep%1000==999){//need to update MVM
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
	    data[igpu].count_g=0;
	    data[igpu].istep=istep;
	    cudaStreamWaitEvent(data[igpu].stream_a[0], data[igpu].event_mvm, 0);//wait for mvm to finish copying.
	}
	if(grad==grad1){
	    grad=grad2;
	}else{
	    grad=grad1;
	}
	for(int ig=0, igpu=0; ig<ng; ig+=gstep, igpu=((igpu+1)%ngpu)){
	    int nleft=ng-ig<gstep?ng-ig:gstep;
	    cudaSetDevice(gpus[igpu]); 
	    GPU_DATA_T *datai=&data[igpu];
	    //One stream handling the memcpy
	    DO(cudaMemcpyAsync(datai->grad->p+ig, grad->p+ig, sizeof(float)*nleft, 
			       cudaMemcpyHostToDevice, *datai->stream_g));
	    //Recored the event when the memcpy is finished
	    DO(cudaEventRecord(datai->event_g[datai->count_g], datai->stream_g[0]));
	    if(!datai->count_g){
		pbeta=&zero;
	    }else{
		pbeta=&one;
	    }
	    //Another stream does the matrix vector multiplication. Wait for the event before executing.
	    DO(cublasSgemv(datai->stream_a[0], CUBLAS_OP_N, nact, nleft, &one, datai->cumvm->p+nact*ig, nact, datai->grad->p+ig, 1, pbeta, datai->act->p, 1));
	    /*The stream stream will wait only for the completion of the most recent host call to cudaEventRecord() on event*/
	    cudaStreamWaitEvent(datai->stream_a[0], datai->event_g[datai->count_g], 0);
	    datai->count_g++;
	}
	//Copy DM commands back to CPU
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=&data[igpu];
	    cudaSetDevice(gpus[igpu]); 
	    cudaMemcpyAsync(dmres->p[igpu]->p, datai->act->p, nact*sizeof(float), cudaMemcpyDeviceToHost, datai->stream_a[0]);
	    DO(cudaEventRecord(datai->event_a, datai->stream_a[0]));//record event when all act are copied so mvm can start.
	    if(datai->copy_mvm){
		int done=0, nleft;
		if(mvm->ny-datai->ic < nc){
		    done=1;
		    nleft=mvm->ny-datai->ic;
		}else{
		    nleft=nc;
		}
		//wait for mvm application to finish before copying.
		DO(cudaStreamWaitEvent(datai->stream_mvm[0], datai->event_a, 0));
		DO(cudaMemcpyAsync(datai->cumvm_next->p+datai->ic*mvm->nx, mvm->p+datai->ic*mvm->nx, sizeof(float)*mvm->nx*nleft, 
				   cudaMemcpyHostToDevice, datai->stream_mvm[0]));
		DO(cudaEventRecord(datai->event_mvm, datai->stream_mvm[0]));
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
	//Wait for all to finish
	for(int igpu=0; igpu<ngpu; igpu++){
	    data[igpu].stream_a->sync();
	}
	//CPU sums them together
	for(int igpu=1; igpu<ngpu; igpu++){
	    for(int iact=0; iact<nact; iact++){
		dmres->p[0]->p[iact]+=dmres->p[igpu]->p[iact];
	    }
	}
	result->p[istep]=dmres->p[0]->p[nact/2];
	timing->p[istep]=toc3;tic;
	if(istep%1000==0 || timing->p[istep]>2e-3){
	    info2("Step %d takes %.0f us\n", istep, timing->p[istep]*1e6);
	}
    }
    swrite(timing, "timing");
    swrite(result, "result");
}
