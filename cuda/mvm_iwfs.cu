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
   \file mvm_iwfs.cu

   Test MVM for a single WFS. Using two GPUs, without networking.
   This is not part of maos executable.

   two ways: 
   1) broadcast gradients to both GPUs with each GPU handle part of actuators
   2) partition gradients to GPUs with GPU handle all actuators.
*/
typedef struct{
    curmat cumvm;//active mvm control matrix
    curmat cumvm_next;//inactive mvm control matrix.
    curmat cumvm1;
    curmat cumvm2;
    curmat grad;
    curmat act;
    stream_t stream_g;
    stream_t stream_a;
    stream_t stream_mvm;
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
    X(mat) *mvm1=X(new)(nact, ng);
    X(mat) *grad1=X(new)(ng, 1);
    curmat cumvm1, cugrad1;
    cp2gpu(cumvm1, mvm1);
    cp2gpu(cugrad1, grad1);
    curmat cuact(nact, 1);
    event_t *event_mvm=(event_t*)malloc(sizeof(event_t)*(nstep+1));
    stream_t stream_mvm;
    Real one=1; 
    for(int istep=0; istep<nstep; istep++){
	DO(cudaEventRecord(event_mvm[istep], stream_mvm));
	DO(CUBL(gemv)(stream_mvm, CUBLAS_OP_N, nact, ng, &one,
		      cumvm1, nact, cugrad1, 1, &one, cuact, 1));
    }
    DO(cudaEventRecord(event_mvm[nstep], stream_mvm));
    stream_mvm.sync();
    smat *timing=snew(nstep,1);
    for(int istep=0; istep<nstep; istep++){
	cudaEventElapsedTime(timing->p+istep, event_mvm[istep], event_mvm[istep+1]);
    }
    sscale(timing, 1e-3);
    writebin(timing, "timing_mvmonly_%dgpu", ngpu);
    free(event_mvm);
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
    X(mat) *mvm1=X(new)(nact, ng);
    X(mat) *mvm2=X(new)(nact, ng);

    X(mat) *mvm=mvm1;

    X(mat) *grad1=X(new)(ng, 1);
    X(mat) *grad2=X(new)(ng, 1);
    X(mat) *grad=grad1;
    rand_t srand;
    seed_rand(&srand, 1);
    X(randu)(mvm1, 1, &srand);
    X(randu)(mvm2, 1,&srand);
    X(randu)(grad1,1, &srand);
    X(randu)(grad2,1, &srand);
    X(cell) *dmres=X(cellnew)(ngpu, 1);
    X(pagelock)(mvm1, mvm2, grad1, grad2, dmres, NULL);
 
    int nc=10;//each time copy nc column of mvm.
    GPU_DATA_T **data=new GPU_DATA_T*[ngpu];
    const int gstep=ng;

    const int sect_gpu=(ng+gstep*ngpu-1)/(gstep*ngpu);
    for(int igpu=0; igpu<ngpu; igpu++){
	cudaSetDevice(gpus[igpu]);
	data[igpu]=new GPU_DATA_T;
	data[igpu]->cumvm1=curmat(mvm1->nx, ng);
	data[igpu]->cumvm2=curmat(mvm2->nx, ng);
	data[igpu]->cumvm=data[igpu]->cumvm1;
	data[igpu]->cumvm_next=data[igpu]->cumvm2;
	cp2gpu(data[igpu]->cumvm1, mvm);
	data[igpu]->grad=curmat(ng, 1);
	data[igpu]->act=curmat(mvm1->nx, 1);
	data[igpu]->event_g=new event_t[sect_gpu];
	data[igpu]->event_a=new event_t[2];
	data[igpu]->gpu=gpus[igpu];
	cudaEventCreateWithFlags(&data[igpu]->event_mvm,cudaEventDisableTiming);
	cudaEventCreateWithFlags(&data[igpu]->event_gall,cudaEventDisableTiming);

	dmres->p[igpu]=X(new)(nact, 1);
	//pthread_create(&threads[igpu], NULL, (thread_fun)mvm_cp, &data[igpu]);
    }
    X(mat) *timing_mvm=X(new)(2, nstep);
    X(mat) *timing=X(new)(nstep, 1);
    X(mat) *result=X(new)(nstep, 1);
    Real one=1; Real zero=0; Real *pbeta;
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
		data[igpu]->copy_mvm=1;
		if(data[igpu]->ic!=0){
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
	    GPU_DATA_T *datai=data[igpu];
	    //One stream handling the memcpy
	    DO(cudaMemcpyAsync(datai->grad.P()+ig, grad->p+ig, sizeof(Real)*nleft, 
			       cudaMemcpyHostToDevice, datai->stream_g));
	    //Recored the event when the memcpy is finished
	    DO(cudaEventRecord(datai->event_g[datai->count_g], datai->stream_g));
	    //Wait for gradient transport to complete before compute.
	    cudaStreamWaitEvent(datai->stream_a, datai->event_g[datai->count_g], 0);
	    DO(cudaEventRecord(datai->event_a[0], datai->stream_a));
	    if(!datai->count_g){
		pbeta=&zero;
	    }else{
		pbeta=&one;
	    }
	    //Another stream does the matrix vector multiplication. Wait for the event before executing.
	    DO(CUBL(gemv)(datai->stream_a, CUBLAS_OP_N, nact, nleft, &one, datai->cumvm.P()+nact*ig, nact, datai->grad.P()+ig, 1, pbeta, datai->act, 1));
	    DO(cudaEventRecord(datai->event_a[1], datai->stream_a));
	    datai->count_g++;
	}
	for(int igpu=0; igpu<ngpu; igpu++){
	    GPU_DATA_T *datai=data[igpu];
	    datai->count_g=0;
	    cudaSetDevice(gpus[igpu]); 
	    //record event when all grads are copied so mvm copy can start.
	    DO(cudaEventRecord(datai->event_gall, datai->stream_g));
	    //Copy DM commands back to CPU
	    cudaMemcpyAsync(dmres->p[igpu]->p, datai->act, nact*sizeof(Real), cudaMemcpyDeviceToHost, datai->stream_a);
	    if(datai->copy_mvm){
		int done=0, nleft;
		if(mvm->ny-datai->ic < nc){
		    done=1;
		    nleft=mvm->ny-datai->ic;
		}else{
		    nleft=nc;
		}
		//wait for gradient transport to finish before copying mvm.
		DO(cudaStreamWaitEvent(datai->stream_mvm, datai->event_gall, 0));
		DO(cudaMemcpyAsync(datai->cumvm_next.P()+datai->ic*mvm->nx, mvm->p+datai->ic*mvm->nx, sizeof(Real)*mvm->nx*nleft, 
				   cudaMemcpyHostToDevice, datai->stream_mvm));
		DO(cudaEventRecord(datai->event_mvm, datai->stream_mvm));
		datai->ic+=nleft;
		if(done){
		    datai->ic=0;
		    datai->copy_mvm=0;
		    curmat tmp=datai->cumvm;
		    datai->cumvm=datai->cumvm_next;
		    datai->cumvm_next=tmp;
		    info("gpu %d switched over at step %d\n", datai->gpu, datai->istep);
		}
	    }
	}
	//Wait for all to finish
	for(int igpu=0; igpu<ngpu; igpu++){
	    data[igpu]->stream_a.sync();
	    float tmp;
	    cudaEventElapsedTime(&tmp, data[igpu]->event_a[0], data[igpu]->event_a[1]);
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
	    data[igpu]->stream_mvm.sync();
	}
	result->p[istep]=dmres->p[0]->p[nact/2];
	timing->p[istep]=toc3;
	if(istep%1000==0 || timing->p[istep]>2e-3){
	    info("Step %d takes %.0f %.0fus\n", istep, timing->p[istep]*1e6, timing_mvm->p[istep]*1e6);
	}
    }
    writebin(timing, "timing_%dgpu", ngpu);
    writebin(timing_mvm, "timing_mvm_%dgpu", ngpu);
    writebin(result, "result_%dgpu", ngpu);
    X(pageunlock)(mvm1, mvm2, grad1, grad2, dmres, NULL);
}
