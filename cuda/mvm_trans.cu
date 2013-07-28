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
#include "cudata.h"

typedef struct MVM_IGPU_T{
    const PARMS_T *parms;
    RECON_T *recon;
    POWFS_T *powfs;
    curcell *mvmig; /*intermediate TomoL result*/
    curcell *mvmfg; /*intermediate FitR result*/
    curmat *mvmt;     /*result: tranpose of MVM calculated by this GPU.*/
    float *FLI;
    smat *residual;
    smat *residualfit;
    long (*curp)[2];
    int ntotact;
    int ntotgrad;
    int load_mvmf; /*intermediate FitR result is for 1) loading, 0) saving.*/
}MVM_IGPU_T;
static void mvm_trans_igpu(thread_t *info){
    TIC;tic;
    double tk_prep=0, tk_fitL=0, tk_fitR=0, tk_TomoL=0, tk_TomoR=0, tk_cp=0;
    MVM_IGPU_T *data=(MVM_IGPU_T*)info->data;
    const PARMS_T *parms=data->parms;
    RECON_T *recon=data->recon;
    smat *residual=data->residual;
    smat *residualfit=data->residualfit;
    long (*curp)[2]=data->curp;
    const int ntotact=data->ntotact;
    const int ntotgrad=data->ntotgrad;
    const int load_mvmf=data->load_mvmf;
    int igpu=info->ithread;
    gpu_set(igpu);
    curecon_t *curecon=cudata->recon;
    curecon_geom *grid=curecon->grid;
    curmat *mvmi=data->mvmig?data->mvmig->p[igpu]:NULL;/*Tomography output, for warm restart*/
    curmat *mvmf=data->mvmfg?data->mvmfg->p[igpu]:NULL;/*loaded FitR output.*/

    curcell *eyec=NULL;/* Only use eyec for CG.*/
    float eye2c[2]={0,1.};
    float *eye2;
    cudaMalloc(&eye2, sizeof(float)*2);
    cudaMemcpy(eye2, eye2c, sizeof(float)*2, cudaMemcpyHostToDevice);
    //const int nwfs=parms->nwfsr;
    const int ndm=parms->ndm;
    /*fit*/
    const float *FLI=data->FLI;
    if(!FLI && !load_mvmf){
	if(parms->fit.square){
	    eyec=curcellnew(ndm, 1, recon->anx, recon->any);
	}else{
	    eyec=curcellnew(ndm, 1, recon->anloc, (long*)0);
	}
    }
 
    curcell *dmfit=load_mvmf?NULL:curcellnew(grid->ndm, 1, grid->anx, grid->any);
    curcell *opdx=curcellnew(recon->npsr, 1, recon->xnx, recon->xny, (float*)(mvmf?1L:0L));
    curcell *opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny, (float*)(mvmi?1L:0L));
    curcell *grad=curcellnew(parms->nwfsr, 1, recon->ngrad, (long*)0, (float*)1);
    if(ntotact==0){
	error("ntotact=0;\n");
    }
    curmat *mvmt=curnew(ntotgrad, info->end-info->start);/*contains result*/
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
		cudaMemcpyAsync(eyec->m->p+iact-1, eye2, 2*sizeof(float),
				cudaMemcpyDeviceToDevice, stream);
	    }else{
		cudaMemcpyAsync(eyec->m->p+iact, eye2+1, sizeof(float), 
				cudaMemcpyDeviceToDevice, stream);
	    }
	}
	if(!recon->actcpl || recon->actcpl->p[curdm]->p[curact]>EPS){
	    if(mvmf) opdx->replace(mvmf->p+(iact-info->start)*mvmf->nx, 0, stream);
	    if(!load_mvmf){
		if(eyec){ /*Fitting operator*/
		    curcellzero(dmfit, stream);//temp
		    residualfit->p[iact]=curecon->FL->solve(&dmfit, eyec, NULL, stream);
		}else{
		    cudaMemcpyAsync(dmfit->m->p, FLI+iact*ntotact, sizeof(float)*ntotact, 
				    cudaMemcpyHostToDevice, stream);
		}
    		tk_fitL+=toc3; tic;
		/*Transpose of fitting operator*/
		curecon->FR->Rt(&opdx, 0.f, dmfit, 1.f, stream);
	    }
	    tk_fitR+=toc3; tic;
	    if(mvmi){
		opdr->replace(mvmi->p+(iact-info->start)*mvmi->nx, 0, stream);
	    }
	    residual->p[iact]=curecon->RL->solve(&opdr, opdx, NULL, stream);
	    tk_TomoL+=toc3; tic;
	    /*Right hand side. output directly to mvmt*/
	    grad->replace(mvmt->p+(iact-info->start)*ntotgrad, 0, stream);
	    curecon->RR->Rt(&grad, 0, opdr, 1, stream);
	    tk_TomoR+=toc3; tic;
	}
    }//for iact
    int nn=ntotgrad*(info->end-info->start)*sizeof(float);
    float *mvmtc=data->mvmt->p+info->start*ntotgrad;
    cudaMemcpyAsync(mvmtc, mvmt->p, nn, cudaMemcpyDeviceToDevice, stream);
    stream.sync();
    curcellfree(dmfit);
    curcellfree(opdx);
    curcellfree(opdr);
    curcellfree(grad);
    curcellfree(eyec);
    curfree(mvmt);
    cudaFree(eye2);
    tk_cp+=toc3;tic;
    info2("GPU %d: Prep %.2f FitL %.2f FitR %.2f TomoL %.1f TomoR %.1f cp %.2f\n", 
	  igpu, tk_prep, tk_fitL, tk_fitR, tk_TomoL, tk_TomoR, tk_cp);
}

void gpu_setup_recon_mvm_trans(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
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
	    ntotact+=recon->anloc[idm];
	} 
	for(int ips=0; ips<recon->npsr; ips++){
	    ntotxloc+=recon->xloc[ips]->nloc;
	}
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    ntotgrad+=recon->ngrad[iwfs];
	}
	

	long (*curp)[2]=(long(*)[2])malloc(ntotact*2*sizeof(long));
	int nact=0;
	for(int idm=0; idm<ndm; idm++){
	    for(int iact=0; iact<recon->anloc[idm]; iact++){
		curp[nact+iact][0]=idm;
		curp[nact+iact][1]=iact;
	    }
	    nact+=recon->anloc[idm];
	}   

	smat *residual=NULL;
	smat *residualfit=NULL;
	if(parms->tomo.alg==1){
	    residual=snew(ntotact, 1);
	}
	if(parms->fit.alg==1){
	    residualfit=snew(ntotact, 1);
	}
	dmat *FLId=NULL; /* MI is inv(FL) for direct methods*/
	float *FLI=NULL;

	/* Loading or saving intermediate TomoL result. */
	smat *mvmi=NULL; 
	if(parms->load.mvmi){
	    mvmi=sread("%s", parms->load.mvmi);
	    if(mvmi->nx!=ntotxloc || mvmi->ny!=ntotact){
		error("loaded mvmi has dimension (%ld, %ld) but we expect (%d, %d)",
		      mvmi->nx, mvmi->ny, ntotxloc, ntotact);
	    }
       	}else if(parms->save.mvmi){
	    mvmi=snew(ntotxloc, ntotact);
	}
	curcell *mvmig=NULL;
	if(mvmi){
	    mvmig=curcellnew(NGPU, 1);
	}

	/* Loading or saving intermediate FitR Result.*/
	smat *mvmf=NULL;
	if(parms->load.mvmf){
	    /*Load FitR FitL results from file. Resembling warm restart case
	      where mvmf is kept in memory*/
	    mvmf=sread("%s", parms->load.mvmf);
	    if(mvmf->nx!=ntotxloc || mvmf->ny!=ntotact){
		error("loaded mvmf has dimension (%ld, %ld) but we expect (%d, %d)",
		      mvmf->nx, mvmf->ny, ntotxloc, ntotact);
	    }
	}else if(parms->save.mvmf){
	    /*save FitR FitL resutls to file, for later loading.*/
	    mvmf=snew(ntotxloc, ntotact);
	}
	curcell *mvmfg=NULL;
	if(mvmf){
	    mvmfg=curcellnew(NGPU, 1);
	}
	if(!parms->load.mvmf){
	    /*Prepare FitR, FitL is don't load fitting results using mvmf*/
	    switch(parms->fit.alg){
	    case 0:{//Use CPU to handle CBS.
		dmat *eye=dnew(ntotact, ntotact);
		daddI(eye, 1);
		FLId=dnew(ntotact, ntotact);
		muv_direct_solve(&FLId, &recon->FL, eye);
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
		FLI=(float*)malloc4async(sizeof(float)*ntotact*ntotact);
		for(long i=0; i<ntotact*ntotact; i++){
		    FLI[i]=(float)FLId->p[i];
		}
		dfree(FLId);
	    }
	}
	gpu_set(gpu_recon);
    	curmat *mvmt=curnew(ntotgrad, ntotact);
	MVM_IGPU_T data={parms, recon, powfs, mvmig, mvmfg, mvmt, FLI, residual, residualfit, curp, ntotact, ntotgrad, parms->load.mvmf?1:0};
	int nthread=NGPU;
	thread_t info[nthread];
	thread_prep(info, 0, ntotact, nthread, mvm_trans_igpu, &data);

	/*Initialyze intermediate TomoL result array in GPU. Send intermediate
	  TomoL results to GPU if load.mvmi is set.*/
	if(mvmi){
	    TIC;tic;
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		mvmig->p[i]=curnew(ntotxloc, info[i].end-info[i].start);
		if(parms->load.mvmi){
		    cudaMemcpy(mvmig->p[i]->p, mvmi->p+info[i].start*ntotxloc, 
			       sizeof(float)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyHostToDevice);
		}
	    }
	    if(parms->load.mvmi){
		toc2("copy mvmi to gpu");
	    }
	}
	/*Initialyze intermediate FitL/FitR result array in GPU. Send
	  intermediate FitL/FitR results to GPU if load.mvmf is set.*/
	if(mvmf){
	    TIC;tic;
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		mvmfg->p[i]=curnew(ntotxloc, info[i].end-info[i].start);
		if(parms->load.mvmf){
		    cudaMemcpy(mvmfg->p[i]->p, mvmf->p+info[i].start*ntotxloc, 
			       sizeof(float)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyHostToDevice);
		}
	    }
	    if(parms->load.mvmf){
		toc2("copy mvmf to gpu");
	    }
	}
	/*Do real MVM control matrix assemble in multiply CPU/GPU*/
	CALL_THREAD(info, nthread, 1);
	/*Copy MVM control matrix results back*/
	{
	    gpu_set(gpu_recon);
	    TIC;tic;
	    stream_t stream;
	    curmat *mvmtt=mvmt->trans(stream);
	    stream.sync();
	    toc2("MVM Reshape in GPU");tic;
	    cp2cpu(&recon->MVM, 0, mvmtt, 1, stream, NULL);
	    stream.sync();
	    curfree(mvmtt);
	    toc2("MVM copy to CPU");
	}
	curfree(mvmt);
	swrite(residual, "MVM_RL_residual");
	swrite(residualfit, "MVM_FL_residual");
	
	if(parms->save.mvmi){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		cudaMemcpy(mvmi->p+info[i].start*ntotxloc, mvmig->p[i]->p,  
			   sizeof(float)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyDeviceToHost);
	    }
	    swrite(mvmi, "MVM_Tomo.bin");
	}
	if(parms->save.mvmf){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		cudaMemcpy(mvmf->p+info[i].start*ntotxloc, mvmfg->p[i]->p,  
			   sizeof(float)*ntotxloc*(info[i].end-info[i].start), cudaMemcpyDeviceToHost);
	    }
	    swrite(mvmf, "MVM_FitL.bin");
	}
	if(mvmig){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		curfree(mvmig->p[i]);
	    }
	    curcellfree(mvmig);
	}
	if(mvmfg){
	    for(int i=0; i<NGPU; i++){
		gpu_set(i);
		curfree(mvmfg->p[i]);
	    }
	    curcellfree(mvmfg);
	}
	sfree(mvmi);
	sfree(mvmf);
	sfree(residual);
	sfree(residualfit);
	free(curp);
	if(FLI) free4async(FLI);
    }//if assemble in gpu
}
