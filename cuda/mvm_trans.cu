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


typedef struct MVM_IGPU_T{
    const PARMS_T *parms;
    RECON_T *recon;
    POWFS_T *powfs;
    curcell *mvmig; /*intermediate TomoL result*/
    curcell *mvmfg; /*intermediate FitR result*/
    smat *mvmt;     /*result: tranpose of MVM calculated by this GPU.*/
    float *FLI;
    smat *residual;
    long (*curp)[2];
    int ntotact;
    int ntotgrad;
    int load_mvmf; /*intermediate FitR result is for 1) loading, 0) saving.*/
}MVM_IGPU_T;
static void gpu_setup_recon_mvmt_igpu(thread_t *info){
    TIC;tic;
    double tk_prep=0, tk_fitL=0, tk_fitR=0, tk_TomoL=0, tk_TomoR=0, tk_cp=0;
    MVM_IGPU_T *data=(MVM_IGPU_T*)info->data;
    const PARMS_T *parms=data->parms;
    RECON_T *recon=data->recon;
    smat *residual=data->residual;
    long (*curp)[2]=data->curp;
    const int ntotact=data->ntotact;
    const int ntotgrad=data->ntotgrad;
    const int load_mvmf=data->load_mvmf;
    int igpu=info->ithread;
    gpu_set(igpu);
    curecon_t *curecon=cudata->recon;
    curmat *mvmi=data->mvmig?data->mvmig->p[igpu]:NULL;/*Tomography output, for warm restart*/
    curmat *mvmf=data->mvmfg?data->mvmfg->p[igpu]:NULL;/*loaded FitR output.*/
    /*Tomography*/
    G_PREFUN prefun=NULL;
    void *predata=NULL;
    if(parms->tomo.precond==1){
	prefun=gpu_Tomo_fdprecond;
	predata=(void*)recon;
    }
    G_CGFUN cg_fun;
    void *cg_data;
    curcell *eyec=NULL;/* Only use eyec for CG.*/
    float eye2[2]={0,1.};
    float eye1[1]={1.};
    //const int nwfs=parms->nwfs;
    const int ndm=parms->ndm;
    if(parms->gpu.fit==1){//sparse matrix
	cg_fun=(G_CGFUN) cumuv;
	cg_data=&curecon->FL;
    }else{
	cg_fun=(G_CGFUN) gpu_FitL;
	cg_data=recon;
    }
    const float *FLI=data->FLI;
    if(!FLI && !load_mvmf){
	if(parms->fit.square){
	    eyec=curcellnew(ndm, 1, recon->anx, recon->any);
	}else{
	    eyec=curcellnew(ndm, 1, recon->anloc, (long*)0);
	}
    }
 
    curcell *dmfit=load_mvmf?NULL:curcellnew(curecon->dmfit);
    curcell *opdx=curcellnew(recon->npsr, 1, recon->xnx, recon->xny, (float*)(mvmf?1:0));
    curcell *opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny, (float*)(mvmi?1:0));
    curcell *grad=curcellnew(parms->nwfs, 1, recon->ngrad, (long*)0, (float*)1);
    if(ntotact==0){
	error("ntotact=0;\n");
    }
    curmat *mvmt=curnew(ntotgrad, info->end-info->start);/*contains result*/
    tk_prep+=toc3;tic;
    stream_t &stream=curecon->cgstream[0];
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
				cudaMemcpyHostToDevice, stream);
	    }else{
		cudaMemcpyAsync(eyec->m->p+iact, eye1, sizeof(float), 
				cudaMemcpyHostToDevice, stream);
	    }
	}
	if(!recon->actcpl || recon->actcpl->p[curdm]->p[curact]>EPS){
	    if(mvmf) opdx->replace(mvmf->p+(iact-info->start)*mvmf->nx, 0, stream);
	    if(!load_mvmf){
		if(eyec){
		    /*Fitting operator*/
		    curcellzero(dmfit, stream);//temp
		    if(gpu_pcg(&dmfit, (G_CGFUN)cg_fun, cg_data, NULL, NULL, eyec, &curecon->cgtmp_fit,
			       parms->recon.warm_restart, parms->fit.maxit, stream)>1.){
			warning("Fit CG not converge.\n");
		    }
		}else{
		    cudaMemcpyAsync(dmfit->m->p, FLI+iact*ntotact, sizeof(float)*ntotact, 
				    cudaMemcpyHostToDevice, stream);
		}
		cudaStreamSynchronize(stream);
    		tk_fitL+=toc3; tic;
		/*Transpose of fitting operator*/
		if(parms->gpu.fit==1){//sparse matrix
		    cumuv_trans(&opdx, 0, &curecon->FR, dmfit, 1);
		}else{
		    gpu_FitRt(&opdx, 0, recon, dmfit, 1);
		}
	    }
	    tk_fitR+=toc3; tic;
	    switch(parms->tomo.alg){
	    case 0:
		if(!opdr->m || !opdx->m){
		    error("opdr and opdx must be continuous\n");
		}
		cuchol_solve(opdr->m->p, curecon->RCl, curecon->RCp, opdx->m->p, stream);
		if(curecon->RUp){
		    curmat *tmp=curnew(curecon->RVp->ny, 1);
		    curmv(tmp->p, 0, curecon->RVp, opdx->m->p, 't', -1, stream);
		    curmv(opdr->m->p, 1, curecon->RUp, tmp->p, 'n', 1, stream);
		    curfree(tmp);
		}
		break;
	    case 1:{
		if(mvmi){
		    opdr->replace(mvmi->p+(iact-info->start)*mvmi->nx, 0, stream);
		}
		if(parms->recon.mvm==2){//disable warm restart in CG using neighboring act.
		    curcellzero(opdr, stream); 
		}
		/*disable the t/t removal lrt in split tomo that creats problem in fdpcg mode*/
		if((residual->p[iact]=gpu_pcg(&opdr, gpu_TomoL, recon, prefun, predata, opdx, &curecon->cgtmp_tomo,
					      parms->recon.warm_restart, parms->tomo.maxit,
					      stream, parms->tomo.cgthres))>0.5){
		    warning2("Tomo CG residual is %.2f for %d\n", residual->p[iact], iact);
		}
	    }
		break;
	    case 2:
		curmv(opdr->m->p, 0, curecon->RMI, opdx->m->p, 'n', 1, stream);
		break;
	    default:
		error("Invalid");
	    }
	    tk_TomoL+=toc3; tic;
	    /*Right hand side. output directly to mvmt*/
	    grad->replace(mvmt->p+(iact-info->start)*ntotgrad, 0, stream);
	    gpu_TomoRt(&grad, 0, recon, opdr, 1, stream);
	    tk_TomoR+=toc3; tic;
	}
    }//for iact
    int nn=ntotgrad*(info->end-info->start)*sizeof(float);
    float *mvmtc=data->mvmt->p+info->start*ntotgrad;
    cudaMemcpyAsync(mvmtc, mvmt->p, nn, cudaMemcpyDeviceToHost, curecon->cgstream[0]);
    cudaStreamSynchronize(curecon->cgstream[0]);
    curcellfree(dmfit);
    curcellfree(opdx);
    curcellfree(opdr);
    curcellfree(grad);
    curcellfree(eyec);
    curfree(mvmt);
    tk_cp+=toc3;tic;
    info2("GPU %d: Prep %.2f FitL %.2f FitR %.2f TomoL %.1f TomoR %.1f cp %.2f\n", 
	  igpu, tk_prep, tk_fitL, tk_fitR, tk_TomoL, tk_TomoR, tk_cp);
}

void gpu_setup_recon_mvm_trans(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    TIC;tic;
    if(parms->recon.alg!=0){
	error("Please adept to LSR\n");
    } 
    if(!parms->load.mvm){
	info2("Assembling MVR MVM in GPU\n");
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
	long (*curp)[2]=(long(*)[2])malloc(ntotact*2*sizeof(long));
	int nact=0;
	for(int idm=0; idm<ndm; idm++){
	    for(int iact=0; iact<recon->anloc[idm]; iact++){
		curp[nact+iact][0]=idm;
		curp[nact+iact][1]=iact;
	    }
	    nact+=recon->anloc[idm];
	}   
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    ntotgrad+=recon->ngrad[iwfs];
	}
	
	smat *residual=NULL;
	if(parms->tomo.alg==1){
	    residual=snew(ntotact, 1);
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
	    case 0:{
		dmat *eye=dnew(ntotact, ntotact);
		daddI(eye, 1);
		FLId=dnew(ntotact, ntotact);
		muv_direct_solve(&FLId, &recon->FL, eye);
		dfree(eye);
		toc("Fit CBS");tic;
	    }
		break;
	    case 1:
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
		dwrite(FLId, "FLId");
		dfree(FLId);
	    }
	}
    	smat *mvmt=snew(ntotgrad, ntotact);
	MVM_IGPU_T data={parms, recon, powfs, mvmig, mvmfg, mvmt, FLI, residual, curp, ntotact, ntotgrad, parms->load.mvmf?1:0};
	int nthread=NGPU;
	thread_t info[NGPU];
	thread_prep(info, 0, ntotact, nthread, gpu_setup_recon_mvmt_igpu, &data);

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
	    TIC;tic;
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
	    dmat *mvmtt=dnew(mvmt->ny, mvmt->nx);
	    for(int iy=0; iy<mvmtt->ny; iy++){
		for(int ix=0; ix<mvmtt->nx; ix++){
		    mvmtt->p[ix+iy*mvmtt->nx]=(double)mvmt->p[iy+ix*mvmt->nx];
		}
	    }
	    toc2("MVM Reshape in CPU 1");
	    sfree(mvmt);
	    d2cell(&recon->MVM, mvmtt, NULL);
	    dfree(mvmtt);
	    toc2("MVM Reshape in CPU 2");
	}
	if(parms->save.setup || parms->save.mvm){
	    dcellwrite(recon->MVM, "MVM.bin");
	}
	swrite(residual, "MVM_RL_residual");
	
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

	free(curp);
	if(FLI) free4async(FLI);
    }//if assemble in gpu
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	gpu_recon_free_do();
	CUDA_SYNC_DEVICE;
    }///for GPU
    if(!parms->sim.mvmport){
	gpu_set(gpu_recon);
	curecon_t *curecon=cudata->recon;
	cp2gpu(&curecon->MVM, recon->MVM);
	dcellfree(recon->MVM);
    }
    toc("MVM Final");
    gpu_print_mem("MVM");
}
