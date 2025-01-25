/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#define TIMING 0

//#include "wfs.h"
#include "recon.h"
#include "pcg.h"
#include "../sim/accphi.h"
#include "../sim/cudata.h"
/**
   \file mvm_direct.cu
   Compute the MVM control matrix in GPU colume by colume
*/
static int* gpu_avail=NULL;//implementes a single stack to manage available GPUs
static int gpu_pos=0;
PNEW(gpu_mutex);
typedef struct{
	const parms_t* parms;
	const recon_t* recon;
	dmat* residual;
	dmat* residualfit;
	long ntotact;
	long ntotgrad;
	long ntotxloc;
	X(mat)* mvmc;
	X(mat)* mvmi;
}mvm_igpu_t;

static void* mvm_direct_igpu(thread_t* info){
	int igpu=info->ithread;
	if(gpu_avail){
		LOCK(gpu_mutex);
		if(gpu_pos>0){
			igpu=gpu_avail[--gpu_pos];
		} else{
			igpu=-1;
			warning("error usage\n");
		}
		UNLOCK(gpu_mutex);
	}
	if(igpu==-1) return NULL;
	gpu_set(igpu);
	info("thread %ld is using GPU %d\n", info->ithread, igpu);

	mvm_igpu_t* data=(mvm_igpu_t*)info->data;
	const parms_t* parms=data->parms;
	const recon_t* recon=data->recon;
	const long ntotact=data->ntotact;
	const long ntotgrad=data->ntotgrad;
	const long ntotxloc=data->ntotxloc;
	curcell grad=curcell(parms->nwfsr, 1, P(recon->ngrad), (long*)NULL);//the I
	curcell opdx=curcell(recon->npsr, 1, P(recon->xnx), P(recon->xny));//right hand size
	curcell opdr;//initialized later
	//curcell *fitx=curcellnew(parms->ndm, 1, P(recon->anloc), (long*)NULL);
	curcell fitr=curcell(parms->ndm, 1, P(recon->anloc), (long*)NULL, (Real*)1);//skip data allocation.
	curmat mvm=curmat(ntotact, info->end-info->start);
	curmat eye2(2, 1);
	dmat* residual=data->residual;
	dmat* residualfit=data->residualfit;
	{
		Real eye2c[2]={0,1.};
		cudaMemcpy(eye2(), eye2c, sizeof(Real)*2, H2D);
	}
	curecon_t* curecon=cudata->recon;
	stream_t stream;
	if(parms->load.mvmf){
		cudaMemcpyAsync(mvm(), data->mvmc->p+info->start*ntotact,
			ntotact*(info->end-info->start)*sizeof(Real),
			H2D, stream);
	}
	curmat mvmi;
	if(parms->load.mvmi||parms->save.mvmi){
		dbg("Creating mvmi of size %ldx %ld\n", ntotxloc, info->end-info->start);
		mvmi=curmat(ntotxloc, info->end-info->start);
		if(parms->load.mvmi){
			cudaMemcpyAsync(mvmi(), data->mvmi->p+info->start*ntotxloc,
				ntotxloc*(info->end-info->start)*sizeof(Real),
				H2D, stream);
		}
		opdr=curcell(recon->npsr, 1, P(recon->xnx), P(recon->xny), (Real*)1);
	} else{
		opdr=curcell(recon->npsr, 1, P(recon->xnx), P(recon->xny));
	}
	TIC;tic;
	curcell tomo_rhs, fit_rhs;
	for(int ig=info->start; ig<info->end; ig++){
		ctoc_init(10);
		info_progress(ig, ntotgrad);
		if(ig){
			cudaMemcpyAsync(grad.M()()+ig-1, eye2(), 2*sizeof(Real), D2D, stream);
		} else{
			cudaMemcpyAsync(grad.M()()+ig, eye2()+1, sizeof(Real), D2D, stream);
		}
		ctoc("copy");
		if(mvmi){
			opdr.replace(mvmi()+(ig-info->start)*ntotxloc, stream);
		}
		curecon->RR->R(tomo_rhs, 0, grad, 1, stream);
		residual->p[ig]=curecon->RL->solve(opdr, tomo_rhs, stream);
		ctoc("Tomo");
		fitr.replace(mvm()+(ig-info->start)*ntotact, stream);
		curecon->FR->R(fit_rhs, 0, opdr, 1, stream);
		residualfit->p[ig]=curecon->FL->solve(fitr, fit_rhs, stream);
		ctoc("Fit");
		ctoc_final("mvm_direct");
	}
	DO(cudaMemcpyAsync(data->mvmc->p+info->start*ntotact,
		mvm(), ntotact*(info->end-info->start)*sizeof(Real),
		D2H, stream));
	if(parms->save.mvmi){
		DO(cudaMemcpyAsync(data->mvmi->p+info->start*ntotxloc,
			mvmi(), ntotxloc*(info->end-info->start)*sizeof(Real),
			D2H, stream));
	}
	stream.sync();
	toc2("Thread %ld mvm", info->ithread);
	if(gpu_avail){
		LOCK(gpu_mutex);
		gpu_avail[gpu_pos++]=igpu;
		UNLOCK(gpu_mutex);
	}
	dbg("thread %ld finish.\n", info->ithread);
	return NULL;
}
/**
   Assemble the MVM control matrix.
*/
void gpu_setup_recon_mvm_direct(const parms_t* parms, recon_t* recon){
	TIC;tic;
	if(parms->recon.alg!=0){
		error("Please adapt to LSR\n");
	}
	if(!parms->load.mvm){
		info("Assembling MVR MVM (direct) in GPU\n");

		long ntotact=0;
		long ntotgrad=0;
		long ntotxloc=0;
		const int ndm=parms->ndm;
		for(int idm=0; idm<ndm; idm++){
			ntotact+=recon->anloc->p[idm];
		}
		for(int ips=0; ips<recon->npsr; ips++){
			ntotxloc+=recon->xloc->p[ips]->nloc;
		}
		for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
			ntotgrad+=recon->ngrad->p[iwfs];
		}
		X(mat)* mvmc=NULL;//control matrix output to CPU
		if(parms->load.mvmf){
			mvmc=X(read)("%s", parms->load.mvmf);
		} else{
			mvmc=X(new)(ntotact, ntotgrad);
		}
		X(mat)* mvmi=NULL;//intermediate result
		if(parms->load.mvmi){
			tic; info("Loading mvmi ...");
			mvmi=X(read)("%s", parms->load.mvmi);
			toc2("done");
		} else if(parms->save.mvmi){
			mvmi=X(new)(ntotxloc, ntotgrad);
		}
		dmat* residual=dnew(ntotgrad, 1);
		dmat* residualfit=dnew(ntotgrad, 1);
		mvm_igpu_t data={parms, recon, residual, residualfit, ntotact, ntotgrad, ntotxloc, mvmc, mvmi};
		int nthread=NGPU;
		if(parms->load.mvmi||parms->save.mvmi){
			/*Each GPU cannot handle all the mvmi if just divide to NGPU
			  runs. Do multiple pass to avoid memroy overflow. Assemes each GPU
			  has more than 2G free space.*/
			int npass=iceil((real)ntotxloc*(real)ntotgrad*sizeof(Real)/NGPU/2000000000);
			dbg("mul=%ld\n", ntotxloc*ntotgrad*sizeof(Real));
			dbg("NGPU=%d\n", NGPU);
			dbg("npass=%d\n", npass);
			nthread=NGPU*npass;
		}
		thread_t *tdata=thread_prep(0, ntotgrad, nthread, mvm_direct_igpu, &data);
		if(nthread>NGPU){
			THREAD_POOL_INIT(NGPU);//limit to only NGPU threads to avoid fighting
			mysleep(1);
			gpu_avail=(int*)calloc(NGPU, sizeof(int));
			for(int igpu=0; igpu<NGPU; igpu++){
				gpu_avail[gpu_pos++]=igpu;
			}
		}
		CALL_THREAD(tdata, 1);
		if(nthread>NGPU){
			THREAD_POOL_INIT(NTHREAD);
			free(gpu_avail); gpu_avail=NULL;
		}
		toc2("Assembly");tic;
		free(tdata);
		dfree(residual);
		dfree(residualfit);
		if(parms->save.mvmf){
			writebin(mvmc, "mvmf.bin");
		}
		if(parms->save.mvmi){
			tic; info("Saving mvmi ...");
			writebin(mvmi, "mvmi.bin");
			toc2("done");
		}
		X(free)(mvmi);
		{
			dmat* dmvmc=dnew(ntotact, ntotgrad);
			for(long i=0; i<ntotact*ntotgrad; i++){
				dmvmc->p[i]=(real)mvmc->p[i];
			}
			recon->MVM=dmvmc;
			X(free)(mvmc);
		}
	}
}

