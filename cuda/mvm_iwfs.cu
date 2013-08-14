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
    event_t *event_g;
    cudaEvent_t event_gall;
    cudaEvent_t event_mvm;
    event_t *event_a;
}GPU_DATA_T;
/*Only do MVM test. just test the bandwidth*/
void mvm_only(int *gpus, int ngpu, int nstep){
#if 1
    //const int nact=7673;//total
    const int nact=6981;//active
    const int nsa=2895;//total. all subaps transported to GPU.
#else
    const int nact=6981;//active
    const int nsa=2700;//active
#endif
    int ng=nsa*2/ngpu;
    cudaSetDevice(gpus[0]);//only do one GPU.
    smat *mvm1=snew(nact, ng);
    smat *grad1=snew(ng, 1);
    curmat *cumvm1=NULL, *cugrad1=NULL;
    cp2gpu(&cumvm1, mvm1);
    cp2gpu(&cugrad1, grad1);
    curmat *cuact=curnew(nact, 1);
    event_t event_mvm[nstep+1];
    stream_t stream_mvm;
    float one=1; 
    for(int istep=0; istep<nstep; istep++){
	DO(cudaEventRecord(event_mvm[istep], stream_mvm));
	DO(cublasSgemv(stream_mvm, CUBLAS_OP_N, nact, ng, &one,
		       cumvm1->p, nact, cugrad1->p, 1, &one, cuact->p, 1));
    }
    DO(cudaEventRecord(event_mvm[nstep], stream_mvm));
    stream_mvm.sync();
    smat *timing=snew(nstep,1);
    for(int istep=0; istep<nstep; istep++){
	cudaEventElapsedTime(timing->p+istep, event_mvm[istep], event_mvm[istep+1]);
    }
    sscale(timing, 1e-3);
    swrite(timing, "timing_mvmonly_%dgpu", ngpu);
}
//Timing only MVM.
/**
   A standalone routine that testes applying MVM for a single WFS and update mvm.*/
void mvm_iwfs(int *gpus, int ngpu, int nstep){
#if 1
    //const int nact=7673;//total
    const int nact=6981;//active
    const int nsa=2895;//total. all subaps transported to GPU.
#else
    const int nact=6981;//active
    const int nsa=2700;//active
#endif
    int ng=nsa*2;
    smat *mvm1=snew(nact, ng);
    smat *mvm2=snew(nact, ng);

    smat *mvm=mvm1;

    smat *grad1=snew(ng, 1);
    smat *grad2=snew(ng, 1);
    smat *grad=grad1;
    rand_t srand;
    seed_rand(&srand, 1);
    srandu(mvm1, 1, &srand);
    srandu(mvm2, 1,&srand);
    srandu(grad1,1, &srand);
    srandu(grad2,1, &srand);
    scell *dmres=scellnew(ngpu, 1);
    spagelock(mvm1, mvm2, grad1, grad2, dmres, NULL);
 
    int nc=10;//each time copy nc column of mvm.
    GPU_DATA_T *data=(GPU_DATA_T*)calloc(ngpu, sizeof(GPU_DATA_T));
    const int gstep=ng;

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
	data[igpu].event_g=new event_t[sect_gpu];
	data[igpu].event_a=new event_t[2];
	data[igpu].gpu=gpus[igpu];
	cudaEventCreateWithFlags(&data[igpu].event_mvm,cudaEventDisableTiming);
	cudaEventCreateWithFlags(&data[igpu].event_gall,cudaEventDisableTiming);

	dmres->p[igpu]=snew(nact, 1);
	//pthread_create(&threads[igpu], NULL, (thread_fun)mvm_cp, &data[igpu]);
    }
    smat *timing_mvm=snew(2, nstep);
    smat *timing=snew(nstep, 1);
    smat *result=snew(nstep, 1);
    float one=1; float zero=0; float *pbeta;
    TIC;
    for(int istep=0; istep<nstep; istep++){
	tic;
	if(istep%1000==999 && 0){//need to update MVM
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
	if(grad==grad1){//switch grad buffer.
	    grad=grad2;
	}else{
	    grad=grad1;
	}
	for(int ig=0, igpu=0; ig<ng; ig+=gstep, igpu=((igpu+1)%ngpu)){
	    int nleft=(ng-ig)<gstep?(ng-ig):gstep;
	    cudaSetDevice(gpus[igpu]); 
	    GPU_DATA_T *datai=&data[igpu];
	    //One stream handling the memcpy
	    DO(cudaMemcpyAsync(datai->grad->p+ig, grad->p+ig, sizeof(float)*nleft, 
			       cudaMemcpyHostToDevice, datai->stream_g[0]));
	    //Recored the event when the memcpy is finished
	    DO(cudaEventRecord(datai->event_g[datai->count_g], datai->stream_g[0]));
	    //Wait for gradient transport to complete before compute.
	    cudaStreamWaitEvent(datai->stream_a[0], datai->event_g[datai->count_g], 0);
	    DO(cudaEventRecord(datai->event_a[0], datai->stream_a[0]));
	    if(!datai->count_g){
		pbeta=&zero;
	    }else{
		pbeta=&one;
	    }
	    //Another stream does the matrix vector multiplication. Wait for the event before executing.
	    DO(cublasSgemv(datai->stream_a[0], CUBLAS_OP_N, nact, nleft, &one, datai->cumvm->p+nact*ig, nact, datai->grad->p+ig, 1, pbeta, datai->act->p, 1));
	    DO(cudaEventRecord(datai->event_a[1], datai->stream_a[0]));
	    datai->count_g++;
	}
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=&data[igpu];
	    datai->count_g=0;
	    cudaSetDevice(gpus[igpu]); 
	    //record event when all grads are copied so mvm copy can start.
	    DO(cudaEventRecord(datai->event_gall, datai->stream_g[0]));
	    //Copy DM commands back to CPU
	    cudaMemcpyAsync(dmres->p[igpu]->p, datai->act->p, nact*sizeof(float), cudaMemcpyDeviceToHost, datai->stream_a[0]);
	    if(datai->copy_mvm){
		int done=0, nleft;
		if(mvm->ny-datai->ic < nc){
		    done=1;
		    nleft=mvm->ny-datai->ic;
		}else{
		    nleft=nc;
		}
		//wait for gradient transport to finish before copying mvm.
		DO(cudaStreamWaitEvent(datai->stream_mvm[0], datai->event_gall, 0));
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
	    float tmp;
	    cudaEventElapsedTime(&tmp, data[igpu].event_a[0], data[igpu].event_a[1]);
	    timing_mvm->p[igpu+ngpu*istep]=tmp*1e-3;
	}
	//CPU sums them together
	for(int igpu=1; igpu<ngpu; igpu++){
	    for(int iact=0; iact<nact; iact++){
		dmres->p[0]->p[iact]+=dmres->p[igpu]->p[iact];
	    }
	}
	//wait for mvm transfer to finish
	for(int igpu=0; igpu<ngpu; igpu++){
	    data[igpu].stream_mvm->sync();
	}
	result->p[istep]=dmres->p[0]->p[nact/2];
	timing->p[istep]=toc3;
	if(istep%1000==0 || timing->p[istep]>2e-3){
	    info2("Step %d takes %.0f %.0fus\n", istep, timing->p[istep]*1e6, timing_mvm->p[istep]*1e6);
	}
    }
    swrite(timing, "timing_%dgpu", ngpu);
    swrite(timing_mvm, "timing_mvm_%dgpu", ngpu);
    swrite(result, "result_%dgpu", ngpu);
    spageunlock(mvm1, mvm2, grad1, grad2, dmres, NULL);
}
