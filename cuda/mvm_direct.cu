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
#include "pcg.h"
#include "curmat.h"
#include "cucmat.h"
#include "accphi.h"
typedef struct{
    const PARMS_T *parms;
    const RECON_T *recon;
    dmat *residual;
    dmat *residualfit;
    int ntotact;
    int ntotgrad;
    smat *mvmc;
}MVM_IGPU_T;
#define TIMING 0
static void mvm_direct_igpu(thread_t *info){
    gpu_set(info->ithread);

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
    const int ntotact=data->ntotact;
    const int ntotgrad=data->ntotgrad;
    curcell *grad=curcellnew(parms->nwfs, 1, recon->ngrad, (long*)NULL);//the I
    curcell *opdx=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);//right hand size
    curcell *opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);//result
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
  
    for(int ig=info->start; ig<info->end; ig++){
	RECORD(0);
	if(info->ithread==0){
	    if(!detached){
		info2("%6d of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", ig*NGPU, ntotgrad);
	    }else if(ig % 100==0){
		info2("%6d of %6d\n", ig*NGPU, ntotact);
	    }
	}
	if(ig){
	    cudaMemcpyAsync(grad->m->p+ig-1, eye2, 2*sizeof(float), cudaMemcpyDeviceToDevice, stream);
	}else{
	    cudaMemcpyAsync(grad->m->p+ig, eye2+1, sizeof(float), cudaMemcpyDeviceToDevice, stream);
	}
	RECORD(1);
	residual->p[ig]=gpu_tomo_do(parms, recon, opdr, opdx, grad, stream);
	RECORD(2);
	fitr->replace(mvm->p+(ig-info->start)*ntotact, 0, stream);
	if(ig && parms->recon.warm_restart && 0){//no help
	    cudaMemcpyAsync(mvm->p+(ig-info->start)*ntotact, mvm->p+(ig-info->start-1)*ntotact, sizeof(float)*ntotact, cudaMemcpyDeviceToDevice,stream);
	}
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
    cudaMemcpyAsync(data->mvmc->p+info->start*ntotact, 
		    mvm->p, ntotact*(info->end-info->start)*sizeof(float), cudaMemcpyDeviceToHost, stream);
    stream.sync();
    curfree(mvm);
    curcellfree(fitx);
    curcellfree(fitr);
    curcellfree(opdx);
    curcellfree(opdr);
    curcellfree(grad);
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
	
	int ntotact=0;
	int ntotgrad=0;
	int ntotxloc=0;
	const int ndm=parms->ndm;
	for(int idm=0; idm<ndm; idm++){
	    ntotact+=recon->anloc[idm];
	} 
	for(int ips=0; ips<recon->npsr; ips++){
	    ntotxloc+=recon->xloc[ips]->nloc;
	}
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    ntotgrad+=recon->ngrad[iwfs];
	}
	smat *mvmc=snew(ntotact, ntotgrad);
	dmat *residual=dnew(ntotgrad,1);
	dmat *residualfit=dnew(ntotgrad, 1);
	MVM_IGPU_T data={parms, recon, residual, residualfit, ntotact, ntotgrad,mvmc};
	int nthread=NGPU;
	thread_t info[nthread];
	thread_prep(info, 0, ntotgrad, nthread, mvm_direct_igpu, &data);
	CALL_THREAD(info, nthread, 1);
	toc("Assembly");tic;
	dfree(residual);
	dfree(residualfit);

	{
	    //first convert from smat to dmat
	    dmat *dmvmc=dnew(ntotact, ntotgrad);
	    for(long i=0; i<ntotact*ntotgrad; i++){
		dmvmc->p[i]=(double)mvmc->p[i];
	    }
	    int ndm=parms->ndm;
	    int nwfs=parms->nwfs;
	    recon->MVM=dcellnew(ndm, nwfs);
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
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
