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
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "utils.h"
#include "wfs.h"
#include "recon.h"
#include "pcg.h"
#include "curmat.h"
#include "cucmat.h"
#include "accphi.h"

static int *gpu_avail=NULL;//implementes a single stack to manage available GPUs
static int gpu_pos=0;
PNEW(gpu_mutex);
typedef struct{
    const PARMS_T *parms;
    const RECON_T *recon;
    dmat *residual;
    dmat *residualfit;
    long ntotact;
    long ntotgrad;
    long ntotxloc;
    smat *mvmc;
    smat *mvmi;
}MVM_IGPU_T;
#define TIMING 0
static void mvm_direct_igpu(thread_t *info){
    int igpu=info->ithread;
    if(gpu_avail){
	LOCK(gpu_mutex);
	igpu=gpu_avail[--gpu_pos];
	UNLOCK(gpu_mutex);
    }
    gpu_set(igpu);
    info("Using GPU %d\n", igpu);
#if TIMING
#define RECORD(i) DO(cudaEventRecord(event[i], stream))
#define NEVENT 4
    cudaEvent_t event[NEVENT]={0};
    float times[NEVENT];
    for(int i=0; i<NEVENT; i++){
	DO(cudaEventCreate(&event[i]));
    }
#else
#define RECORD(i)
#endif

    MVM_IGPU_T *data=(MVM_IGPU_T*)info->data;
    const PARMS_T *parms=data->parms;
    const RECON_T *recon=data->recon;
    const long ntotact=data->ntotact;
    const long ntotgrad=data->ntotgrad;
    const long ntotxloc=data->ntotxloc;
    curcell *grad=curcellnew(parms->nwfsr, 1, recon->ngrad, (long*)NULL);//the I
    curcell *opdx=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);//right hand size
    curcell *opdr=NULL;//initialized later
    curcell *fitx=curcellnew(parms->ndm, 1, recon->anloc, (long*)NULL);
    curcell *fitr=curcellnew(parms->ndm, 1, recon->anloc, (long*)NULL, (float*)1);//skip data allocation.
    curmat *mvm=curnew(ntotact, info->end-info->start);
    float *eye2; cudaMalloc(&eye2, sizeof(float)*2);
    dmat *residual=data->residual;
    dmat *residualfit=data->residualfit;
    {
	float eye2c[2]={0,1.};
	cudaMemcpy(eye2, eye2c, sizeof(float)*2, cudaMemcpyHostToDevice);
    }
    curecon_t *curecon=cudata->recon;
    stream_t &stream=curecon->cgstream[0];
    if(parms->load.mvmf){
	cudaMemcpyAsync(mvm->p, data->mvmc->p+info->start*ntotact, 
			ntotact*(info->end-info->start)*sizeof(float), 
			cudaMemcpyHostToDevice, stream);
    }
    curmat *mvmi=NULL;
    if(parms->load.mvmi || parms->save.mvmi){
	info("Creating mvmi of size %ldx %ld\n", ntotxloc, info->end-info->start);
	mvmi=curnew(ntotxloc, info->end-info->start);
	if(parms->load.mvmi){
	    cudaMemcpyAsync(mvmi->p, data->mvmi->p+info->start*ntotxloc,
			    ntotxloc*(info->end-info->start)*sizeof(float),
			    cudaMemcpyHostToDevice, stream);
	}
	opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny, (float*)1);
    }else{
	opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
    }
    TIC;tic;
    for(int ig=info->start; ig<info->end; ig++){
	RECORD(0);
	if(info->ithread==0){
	    if(!detached){
		info2("%6d of %6ld\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", ig*NGPU, ntotgrad);
	    }else if(ig % 100==0){
		info2("%6d of %6ld\n", ig*NGPU, ntotgrad);
	    }
	}
	if(ig){
	    cudaMemcpyAsync(grad->m->p+ig-1, eye2, 2*sizeof(float), cudaMemcpyDeviceToDevice, stream);
	}else{
	    cudaMemcpyAsync(grad->m->p+ig, eye2+1, sizeof(float), cudaMemcpyDeviceToDevice, stream);
	}
	RECORD(1);
	if(mvmi){
	    opdr->replace(mvmi->p+(ig-info->start)*ntotxloc, 0, stream);
	}
	residual->p[ig]=gpu_tomo_do(parms, recon, opdr, opdx, grad, stream);
	RECORD(2);
	fitr->replace(mvm->p+(ig-info->start)*ntotact, 0, stream);
	residualfit->p[ig]=gpu_fit_do(parms, recon, fitr, fitx, opdr, stream);
	RECORD(3);
#if TIMING
	stream.sync();
	for(int i=1; i<NEVENT; i++){
	    DO(cudaEventElapsedTime(&times[i], event[i-1], event[i]));
	    times[i]*=1e3;//micro-second
	}
	info2("copy=%3.0f, Tomo=%3.0f, Fit=%3.0f\n", times[1], times[2], times[3]);
#endif	
    }
    DO(cudaMemcpyAsync(data->mvmc->p+info->start*ntotact, 
		       mvm->p, ntotact*(info->end-info->start)*sizeof(float), 
		       cudaMemcpyDeviceToHost, stream));
    if(parms->save.mvmi){
	DO(cudaMemcpyAsync(data->mvmi->p+info->start*ntotxloc,
			   mvmi->p, ntotxloc*(info->end-info->start)*sizeof(float),
			   cudaMemcpyDeviceToHost, stream));
    }
    stream.sync();
    toc2("Thread %ld mvm", info->ithread);
    curfree(mvm);
    curfree(mvmi);
    curcellfree(fitx);
    curcellfree(fitr);
    curcellfree(opdx);
    curcellfree(opdr);
    curcellfree(grad);
    if(gpu_avail){
	LOCK(gpu_mutex);
	gpu_avail[gpu_pos++]=igpu;
	UNLOCK(gpu_mutex);
    }
}
/**
   Assemble the MVM control matrix.
*/
void gpu_setup_recon_mvm_direct(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    TIC;tic;
    if(parms->recon.alg!=0){
	error("Please adapt to LSR\n");
    } 
    if(!parms->load.mvm){
	info2("Assembling MVR MVM (direct) in GPU\n");
	
	long ntotact=0;
	long ntotgrad=0;
	long ntotxloc=0;
	const int ndm=parms->ndm;
	for(int idm=0; idm<ndm; idm++){
	    ntotact+=recon->anloc[idm];
	} 
	for(int ips=0; ips<recon->npsr; ips++){
	    ntotxloc+=recon->xloc[ips]->nloc;
	}
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    ntotgrad+=recon->ngrad[iwfs];
	}
	smat *mvmc=NULL;//control matrix output to CPU
	if(parms->load.mvmf){
	    mvmc=sread("%s", parms->load.mvmf);
	}else{
	    mvmc=snew(ntotact, ntotgrad);
	}
	smat *mvmi=NULL;//intermediate result
	if(parms->load.mvmi){
	    TIC;tic; info2("Loading mvmi ...");
	    mvmi=sread("%s", parms->load.mvmi);
	    toc2("done");
	}else if(parms->save.mvmi){
	    mvmi=snew(ntotxloc, ntotgrad);
	}
	dmat *residual=dnew(ntotgrad,1);
	dmat *residualfit=dnew(ntotgrad, 1);
	MVM_IGPU_T data={parms, recon, residual, residualfit, ntotact, ntotgrad, ntotxloc, mvmc, mvmi};
	int nthread=NGPU;
	if(parms->load.mvmi || parms->save.mvmi){
	    /*Each GPU cannot handle all the mvmi if just divide to NGPU
	      runs. Do multiple pass to avoid memroy overflow. Assemes each GPU
	      has more than 2G free space.*/
	    int npass=iceil((double)ntotxloc*(double)ntotgrad*sizeof(float)/NGPU/2000000000);
	    info("mul=%ld\n", ntotxloc*ntotgrad*sizeof(float));
	    info("NGPU=%d\n", NGPU);
	    info("npass=%d\n", npass);
	    nthread=NGPU*npass;
	}
	thread_t info[nthread];
	thread_prep(info, 0, ntotgrad, nthread, mvm_direct_igpu, &data);
	if(nthread>NGPU){
	    THREAD_POOL_INIT(NGPU);//limit to only NGPU threads to avoid fighting
	    sleep(1);
	    gpu_avail=(int*)calloc(NGPU, sizeof(int));
	    for(int igpu=0; igpu<NGPU; igpu++){
		gpu_avail[gpu_pos++]=igpu;
	    }
	}
	CALL_THREAD(info, nthread, 1);
	if(nthread>NGPU){
	    THREAD_POOL_INIT(parms->sim.nthread);
	    free(gpu_avail); gpu_avail=NULL;
	}
	toc("Assembly");tic;
	dfree(residual);
	dfree(residualfit);
	if(parms->save.mvmf){
	    swrite(mvmc, "mvmf.bin");
	}
	if(parms->save.mvmi){
	    TIC;tic; info2("Saving mvmi ...");
	    swrite(mvmi, "mvmi.bin");
	    toc2("done");
	}
	sfree(mvmi);
	{
	    //first convert from smat to dmat
	    dmat *dmvmc=dnew(ntotact, ntotgrad);
	    for(long i=0; i<ntotact*ntotgrad; i++){
		dmvmc->p[i]=(double)mvmc->p[i];
	    }
	    int ndm=parms->ndm;
	    int nwfs=parms->nwfsr;
	    recon->MVM=dcellnew(ndm, nwfs);
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfsr[iwfs].powfs;
		if(!parms->powfs[ipowfs].skip){
		    for(int idm=0; idm<ndm; idm++){
			recon->MVM->p[idm+ndm*iwfs]=dnew(recon->anloc[idm], powfs[ipowfs].saloc->nloc*2);
		    }
		}
	    }
	    d2cell(&recon->MVM, dmvmc, NULL);
	    toc("Copy");
	    dfree(dmvmc);
	    sfree(mvmc);
	}
    }
}
