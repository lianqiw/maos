/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "wfs.h"
#include "recon.h"
#include "pcg.h"
#include "curmat.h"
#include "cucmat.h"
#include "accphi.h"
#include "cudata.h"

/**
   \file mvm_trans.cu
   Compute the transpose of MVM control matrix in GPU colume by colume.
*/
typedef struct MVM_IGPU_T{
    const PARMS_T *parms;
    RECON_T *recon;
    curcell &mvmig; /*intermediate TomoL result*/
    curcell &mvmfg; /*intermediate FitR result*/
    curmat &mvmt;     /*result: tranpose of MVM calculated by this GPU.*/
    Real *FLI;
    X(mat) *residual;
    X(mat) *residualfit;
    long (*curp)[2];
    int ntotact;
    int ntotgrad;
    int load_mvmf; /*intermediate FitR result is for 1) loading, 0) saving.*/
}MVM_IGPU_T;
#define MVM_DEBUG 0
static void mvm_trans_igpu(thread_t *info){
    TIC;tic;
    double tk_prep=0, tk_fitL=0, tk_fitR=0, tk_TomoL=0, tk_TomoR=0, tk_cp=0;
    MVM_IGPU_T *data=(MVM_IGPU_T*)info->data;
    const PARMS_T *parms=data->parms;
    RECON_T *recon=data->recon;
    X(mat) *residual=data->residual;
    X(mat) *residualfit=data->residualfit;
    long (*curp)[2]=data->curp;
    const int ntotact=data->ntotact;
    const int ntotgrad=data->ntotgrad;
    const int load_mvmf=data->load_mvmf;
    int igpu=info->ithread;
    gpu_set(igpu);
    cuda_recon::curecon_t *curecon=cudata->recon;
    cuda_recon::curecon_geom *grid=curecon->grid;
    curmat dummy;
    curmat &mvmi=data->mvmig?data->mvmig[igpu]:dummy;/*Tomography output, for warm restart*/
    curmat &mvmf=data->mvmfg?data->mvmfg[igpu]:dummy;/*loaded FitR output.*/

    curcell eyec;/* Only use eyec for CG.*/
    Real eye2c[2]={0,1.};
    curmat eye2(2,1);
    cudaMemcpy(eye2.P(), eye2c, sizeof(Real)*2, cudaMemcpyHostToDevice);
    //const int nwfs=parms->nwfsr;
    const int ndm=parms->ndm;
    /*fit*/
    const Real *FLI=data->FLI;
    if(!FLI && !load_mvmf){
	if(parms->fit.square){
	    eyec=curcell(ndm, 1, recon->anx->p, recon->any->p);
	}else{
	    eyec=curcell(ndm, 1, recon->anloc->p, (long*)0);
	}
    }
    curcell dmfit;
    if(!load_mvmf){
	if(parms->fit.square){
	    dmfit=curcell(grid->ndm, 1, grid->anx, grid->any);
	}else{
	    dmfit=curcell(grid->ndm, 1, recon->anloc->p, (long*)0);
	}
    }
    curcell opdx(recon->npsr, 1, recon->xnx->p, recon->xny->p, (Real*)(mvmf?1L:0L));
    curcell opdr(recon->npsr, 1, recon->xnx->p, recon->xny->p, (Real*)(mvmi?1L:0L));
    curcell grad(parms->nwfsr, 1, recon->ngrad->p, (long*)0, (Real*)1);
    if(ntotact==0){
	error("ntotact=0;\n");
    }
    curmat mvmt(ntotgrad, info->end-info->start);/*contains result*/
    tk_prep+=toc3;tic;
    stream_t stream;
    for(int iact=info->start; iact<info->end; iact++){
	int curdm=curp[iact][0];
	int curact=curp[iact][1];
	if(info->ithread==0){
	    if(!detached){
		info2("%6d of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", iact*NGPU, ntotact);
	    }else if(iact % 100==0){
		info2("%6d of %6d\n", iact*NGPU, ntotact);
	    }
	}
	if(eyec){
	    if(iact){
		cudaMemcpyAsync(eyec.M().P()+iact-1, eye2.P(), 2*sizeof(Real),
				cudaMemcpyDeviceToDevice, stream);
	    }else{
		cudaMemcpyAsync(eyec.M().P()+iact, eye2.P()+1, sizeof(Real), 
				cudaMemcpyDeviceToDevice, stream);
	    }
#if MVM_DEBUG == 1
		    cuwrite(eyec, "mvm_eyec_%d", iact);
#endif
	}
	if(!recon->actcpl || recon->actcpl->p[curdm]->p[curact]>EPS){
	    if(mvmf){
		opdx.replace(mvmf.P()+(iact-info->start)*mvmf.Nx(), stream);
	    }
	    if(!load_mvmf){
		if(eyec){ /*Fitting operator*/
		    cuzero(dmfit, stream);//temp
		    residualfit->p[iact]=curecon->FL->solve(dmfit, eyec, stream);
#if MVM_DEBUG == 1
		    cuwrite(dmfit, "mvm_dmfit_%d", iact);
#endif
		}else{
		    cudaMemcpyAsync(dmfit.M().P(), FLI+iact*ntotact, sizeof(Real)*ntotact, 
				    cudaMemcpyHostToDevice, stream);
		}
    		tk_fitL+=toc3; tic;
		/*Transpose of fitting operator*/
		curecon->FR->Rt(opdx, 0.f, dmfit, 1.f, stream);
#if MVM_DEBUG == 1
		cuwrite(opdx, "mvm_opdx_%d", iact);
#endif
	    }
	    tk_fitR+=toc3; tic;
	    if(mvmi){
		opdr.replace(mvmi.P()+(iact-info->start)*mvmi.Nx(),  stream);
	    }else{
		cuzero(opdr.M(), stream);
	    }
	    residual->p[iact]=curecon->RL->solve(opdr, opdx, stream);
#if MVM_DEBUG == 1
	    cuwrite(opdx, "mvm_opdr_%d", iact);
#endif
	    tk_TomoL+=toc3; tic;
	    /*Right hand side. output directly to mvmt*/
	    grad.replace(mvmt.P()+(iact-info->start)*ntotgrad, stream);
	    curecon->RR->Rt(grad, 0, opdr, 1, stream);
#if MVM_DEBUG == 1
	    cuwrite(grad, "mvm_grad_%d", iact);
	    if(iact>10){
		_Exit(0);
	    }
#endif
	    tk_TomoR+=toc3; tic;
	}
    }//for iact
    int nn=ntotgrad*(info->end-info->start)*sizeof(Real);
    Real *mvmtc=data->mvmt.P()+info->start*ntotgrad;
    cudaMemcpyAsync(mvmtc, mvmt.P(), nn, cudaMemcpyDeviceToDevice, stream);
    stream.sync();
    tk_cp+=toc3;tic;
    info2("GPU %d: Prep %.2f FitL %.2f FitR %.2f TomoL %.1f TomoR %.1f cp %.2f\n", 
	  igpu, tk_prep, tk_fitL, tk_fitR, tk_TomoL, tk_TomoR, tk_cp);
}

void gpu_setup_recon_mvm_trans(const PARMS_T *parms, RECON_T *recon){
    TIC;tic;
    if(parms->recon.alg!=0){
	error("Please adapt to LSR\n");
    } 
    if(!parms->load.mvm){
	info2("Assembling MVR MVM (transpose) in GPU\n");
	int ntotact=0;
	int ntotgrad=0;
	int ntotxloc=0;
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
	

	long (*curp)[2]=(long(*)[2])malloc(ntotact*2*sizeof(long));
	int nact=0;
	for(int idm=0; idm<ndm; idm++){
	    for(int iact=0; iact<recon->anloc->p[idm]; iact++){
		curp[nact+iact][0]=idm;
		curp[nact+iact][1]=iact;
	    }
	    nact+=recon->anloc->p[idm];
	}   

	X(mat) *residual=NULL;
	X(mat) *residualfit=NULL;
	if(parms->tomo.alg==1){
	    residual=X(new)(ntotact, 1);
	}
	if(parms->fit.alg==1){
	    residualfit=X(new)(ntotact, 1);
	}
	dmat *FLId=NULL; /* MI is inv(FL) for direct methods*/
	Real *FLI=NULL;

	/* Loading or saving intermediate TomoL result. */
	X(mat) *mvmi=NULL; 
	if(parms->load.mvmi){
	    mvmi=X(read)("%s", parms->load.mvmi);
	    if(mvmi->nx!=ntotxloc || mvmi->ny!=ntotact){
		error("loaded mvmi has dimension (%ld, %ld) but we expect (%d, %d)",
		      mvmi->nx, mvmi->ny, ntotxloc, ntotact);
	    }
       	}else if(parms->save.mvmi){
	    mvmi=X(new)(ntotxloc, ntotact);
	}
	curcell mvmig;
	if(mvmi){
	    mvmig=curcell(NGPU, 1);
	}

	/* Loading or saving intermediate FitR Result.*/
	X(mat) *mvmf=NULL;
	if(parms->load.mvmf){
	    /*Load FitR FitL results from file. Resembling warm restart case
	      where mvmf is kept in memory*/
	    mvmf=X(read)("%s", parms->load.mvmf);
	    if(mvmf->nx!=ntotxloc || mvmf->ny!=ntotact){
		error("loaded mvmf has dimension (%ld, %ld) but we expect (%d, %d)",
		      mvmf->nx, mvmf->ny, ntotxloc, ntotact);
	    }
	}else if(parms->save.mvmf){
	    /*save FitR FitL resutls to file, for later loading.*/
	    mvmf=X(new)(ntotxloc, ntotact);
	}
	curcell mvmfg;
	if(mvmf){
	    mvmfg=curcell(NGPU, 1);
	}
	if(!parms->load.mvmf){
	    /*Prepare FitR, FitL is don't load fitting results using mvmf*/
	    switch(parms->fit.alg){
	    case 0:{//Use CPU to handle CBS.
		dmat *eye=dnew(ntotact, ntotact);
		daddI(eye, 1);
		FLId=dnew(ntotact, ntotact);
		muv_direct_solve_mat(&FLId, &recon->FL, eye);
		dfree(eye);
		toc("Fit CBS");tic;
	    }
		break;
	    case 1://Use GPU.
		break;
	    case 2:
		FLId=dref(recon->FL.MI);
		break;
	    default:
		error("Invalid fit.alg=%d\n", parms->fit.alg);
	    }
	    if(FLId){
		FLI=(Real*)malloc4async(sizeof(Real)*ntotact*ntotact);
		for(long i=0; i<ntotact*ntotact; i++){
		    FLI[i]=(Real)FLId->p[i];
		}
		dfree(FLId);
	    }
	}
	gpu_set(cudata_t::recongpu);
    	curmat mvmt=curmat(ntotgrad, ntotact);
	MVM_IGPU_T data={parms, recon, mvmig, mvmfg, mvmt, FLI, residual, residualfit, curp, ntotact, ntotgrad, parms->load.mvmf?1:0};
	int nthread=NGPU;
	thread_t info[nthread];
	thread_prep(info, 0, ntotact, nthread, mvm_trans_igpu, &data);

	/*Initialyze intermediate TomoL result array in GPU. Send intermediate
	  TomoL results to GPU if load.mvmi is set.*/
	if(mvmi){
	    tic;
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		mvmig[i]=curmat(ntotxloc, info[i].end-info[i].start);
		if(parms->load.mvmi){
		    cudaMemcpy(mvmig[i].P(), mvmi->p+info[i].start*ntotxloc, 
			       sizeof(Real)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyHostToDevice);
		}
	    }
	    if(parms->load.mvmi){
		toc2("copy mvmi to gpu");
	    }
	}
	/*Initialyze intermediate FitL/FitR result array in GPU. Send
	  intermediate FitL/FitR results to GPU if load.mvmf is set.*/
	if(mvmf){
	    tic;
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		mvmfg[i]=curmat(ntotxloc, info[i].end-info[i].start);
		if(parms->load.mvmf){
		    cudaMemcpy(mvmfg[i].P(), mvmf->p+info[i].start*ntotxloc, 
			       sizeof(Real)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyHostToDevice);
		}
	    }
	    if(parms->load.mvmf){
		toc2("copy mvmf to gpu");
	    }
	}
	/*Do real MVM control matrix assemble in multiply CPU/GPU*/
	CALL_THREAD(info, 1);
	/*Copy MVM control matrix results back*/
	{
	    gpu_set(cudata_t::recongpu);
	    tic;
	    stream_t stream;
	    curmat mvmtt=mvmt.trans(stream);
	    stream.sync();
	    toc2("MVM Reshape in GPU");tic;
	    cp2cpu(&recon->MVM, mvmtt, stream);
	    stream.sync();
	    toc2("MVM copy to CPU");
	}
	if(parms->save.setup){
	    writebin(residual, "MVM_RL_residual");
	    writebin(residualfit, "MVM_FL_residual");
	}
	if(parms->save.mvmi){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		cudaMemcpy(mvmi->p+info[i].start*ntotxloc, mvmig[i].P(),  
			   sizeof(Real)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyDeviceToHost);
	    }
	    writebin(mvmi, "MVM_Tomo.bin");
	}
	if(parms->save.mvmf){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		cudaMemcpy(mvmf->p+info[i].start*ntotxloc, mvmfg[i].P(),  
			   sizeof(Real)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyDeviceToHost);
	    }
	    writebin(mvmf, "MVM_FitL.bin");
	}

	X(free)(mvmi);
	X(free)(mvmf);
	X(free)(residual);
	X(free)(residualfit);
	free(curp);
	if(FLI) free4async(FLI);
    }//if assemble in gpu
}

