/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

//#include "wfs.h"
#include "recon.h"
#include "pcg.h"
#include "../sim/accphi.h"
#include "../sim/cudata.h"

/**
   \file mvm_trans.cu
   Compute the transpose of MVM control matrix in GPU colume by colume.
*/
typedef struct mvm_igpu_t{
	const parms_t* parms;
	recon_t* recon;
	curcell& mvmig; /*intermediate TomoL result*/
	curcell& mvmfg; /*intermediate FitR result*/
	curmat& mvmt;   /*control matrix transposed.*/
	Real* FLI;
	X(mat)* residual;
	X(mat)* residualfit;
	long(*curp)[2];
	int ntotact;
	int ntotgrad;
	int load_mvmf; /*intermediate FitR result is for 1) loading, 0) saving.*/
}mvm_igpu_t;
#define MVM_DEBUG 0
static void* mvm_trans_igpu(thread_t* info){
	TIC;tic;
	real tk_prep=0, tk_fitL=0, tk_fitR=0, tk_TomoL=0, tk_TomoR=0, tk_cp=0;
	mvm_igpu_t* data=(mvm_igpu_t*)info->data;
	const parms_t* parms=data->parms;
	recon_t* recon=data->recon;
	X(mat)* residual=data->residual;
	X(mat)* residualfit=data->residualfit;
	long(*curp)[2]=data->curp;
	const int ntotact=data->ntotact;
	const int ntotgrad=data->ntotgrad;
	const int load_mvmf=data->load_mvmf;
	int igpu=info->ithread;
	gpu_set(igpu);
	curecon_t* curecon=cudata->recon;
	curecon_geom* grid=curecon->grid;
	curmat dummy;
	curmat& mvmi=data->mvmig?data->mvmig[igpu]:dummy;/*Tomography output, for warm restart*/
	curmat& mvmf=data->mvmfg?data->mvmfg[igpu]:dummy;/*loaded FitR output.*/

	curcell eyec;/* Only use eyec for CG.*/
	Real eye2c[2]={0,1.};
	curmat eye2(2, 1);
	cudaMemcpy(eye2(), eye2c, sizeof(Real)*2, H2D);
	//const int nwfs=parms->nwfsr;
	const int ndm=parms->ndm;
	const lmat *anmod=parms->recon.modal?recon->anmod:recon->anloc;
	/*fit*/
	const Real* FLI=data->FLI;
	if(!FLI&&!load_mvmf){
		if(parms->fit.square && !parms->recon.modal){
			eyec=curcell(ndm, 1, P(recon->anx), P(recon->any));
		}else{
			eyec=curcell(ndm, 1, P(anmod), (long*)0);
		}
	}
	curcell dmrecon;
	if(!load_mvmf){
		if(parms->fit.square&&!parms->recon.modal){
			dmrecon=curcell(grid->ndm, 1, grid->anx, grid->any);
		}else{
			dmrecon=curcell(grid->ndm, 1, P(anmod), (long*)0);
		}
	}
	curcell opdx(recon->npsr, 1, P(recon->xnx), P(recon->xny), (Real*)(mvmf?1L:0L));
	curcell opdr(recon->npsr, 1, P(recon->xnx), P(recon->xny), (Real*)(mvmi?1L:0L));
	curcell grad(parms->nwfsr, 1, P(recon->ngrad), (long*)0, (Real*)1);
	if(ntotact==0){
		error("ntotact=0;\n");
	}
	curmat mvmt(ntotgrad, info->end-info->start);/*contains result*/
	tk_prep+=toc3;tic;
	stream_t stream;
	for(int iact=info->start; iact<info->end; iact++){
		int curdm=curp[iact][0];
		int curact=curp[iact][1];
		info_progress(iact, ntotact);
		if(eyec){
			if(iact){
				cudaMemcpyAsync(eyec.M()()+iact-1, eye2(), 2*sizeof(Real),
					D2D, stream);
			} else{
				cudaMemcpyAsync(eyec.M()()+iact, eye2()+1, sizeof(Real),
					D2D, stream);
			}
#if MVM_DEBUG == 1
			cuwrite(eyec, stream, "mvm_eyec_%d", iact);
#endif
		}
		if(!recon->actcpl||recon->actcpl->p[curdm]->p[curact]>EPS){
			if(mvmf){
				opdx.replace(mvmf()+(iact-info->start)*mvmf.Nx(), stream);
			}
			if(!load_mvmf){
				if(eyec){ /*Fitting operator*/
					cuzero(dmrecon, stream);//temp
					residualfit->p[iact]=curecon->FL->solve(dmrecon, eyec, stream);
#if MVM_DEBUG == 1
					cuwrite(dmrecon, stream, "mvm_dmrecon_%d", iact);
#endif
				} else{
					cudaMemcpyAsync(dmrecon.M()(), FLI+iact*ntotact, sizeof(Real)*ntotact,
						H2D, stream);
				}
				tk_fitL+=toc3; tic;
			/*Transpose of fitting operator*/
				curecon->FR->Rt(opdx, 0.f, dmrecon, 1.f, stream);
#if MVM_DEBUG == 1
				cuwrite(opdx, stream, "mvm_opdx_%d", iact);
#endif
			}
			tk_fitR+=toc3; tic;
			if(mvmi){
				opdr.replace(mvmi()+(iact-info->start)*mvmi.Nx(), stream);
			} else{
				cuzero(opdr.M(), stream);
			}
			residual->p[iact]=curecon->RL->solve(opdr, opdx, stream);
#if MVM_DEBUG == 1
			cuwrite(opdx, stream, "mvm_opdr_%d", iact);
#endif
			tk_TomoL+=toc3; tic;
			/*Right hand side. output directly to mvmt*/
			grad.replace(mvmt()+(iact-info->start)*ntotgrad, stream);
			curecon->RR->Rt(grad, 0, opdr, 1, stream);
#if MVM_DEBUG == 1
			cuwrite(grad, stream, "mvm_grad_%d", iact);
			if(iact>10){
				_Exit(0);
			}
#endif
			tk_TomoR+=toc3; tic;
		}
	}//for iact
	int nn=ntotgrad*(info->end-info->start)*sizeof(Real);
	Real* mvmtc=data->mvmt()+info->start*ntotgrad;
	cudaMemcpyAsync(mvmtc, mvmt(), nn, D2D, stream);
	stream.sync();
	tk_cp+=toc3;tic;
	info("GPU %d: Prep %.2f FitL %.2f FitR %.2f TomoL %.1f TomoR %.1f cp %.2f\n",
		igpu, tk_prep, tk_fitL, tk_fitR, tk_TomoL, tk_TomoR, tk_cp);
	return NULL;
}

void gpu_setup_recon_mvm_trans(const parms_t* parms, recon_t* recon){
	TIC;tic;
	if(parms->recon.alg!=0){
		error("Please adapt to LSR\n");
	}
	if(parms->load.mvm){
		return;
	}
	gpu_set(cuglobal->recongpu);

	info("Assembling MVR MVM (transpose) in GPU\n");
	int ntotact=0;
	int ntotgrad=0;
	int ntotxloc=0;
	const int ndm=parms->ndm;
	const lmat *anmod=parms->recon.modal?recon->anmod:recon->anloc;
	for(int idm=0; idm<ndm; idm++){
		ntotact+=anmod->p[idm];
	}
	for(int ips=0; ips<recon->npsr; ips++){
		ntotxloc+=recon->xloc->p[ips]->nloc;
	}
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		ntotgrad+=recon->ngrad->p[iwfs];
	}
	//Store control matrix in cuglobal.
	curmat mvmt=curmat(ntotgrad, ntotact);
	{
		long(*curp)[2]=(long(*)[2])malloc(ntotact*2*sizeof(long));
		int nact=0;
		for(int idm=0; idm<ndm; idm++){
			for(int iact=0; iact<anmod->p[idm]; iact++){
				curp[nact+iact][0]=idm;
				curp[nact+iact][1]=iact;
			}
			nact+=anmod->p[idm];
		}

		X(mat)* residual=NULL;
		X(mat)* residualfit=NULL;
		if(parms->tomo.alg==1){
			residual=X(new)(ntotact, 1);
		}
		if(parms->fit.alg==1){
			residualfit=X(new)(ntotact, 1);
		}
		//Real *FLI=NULL;
		Array<Real, Pinned>FLI;
		/* Loading or saving intermediate TomoL result. */
		X(mat)* mvmi=NULL;
		if(parms->load.mvmi){
			mvmi=X(read)("%s", parms->load.mvmi);
			if(mvmi->nx!=ntotxloc||mvmi->ny!=ntotact){
				error("loaded mvmi has dimension (%ld, %ld) but we expect (%d, %d)",
					mvmi->nx, mvmi->ny, ntotxloc, ntotact);
			}
		} else if(parms->save.mvmi){
			mvmi=X(new)(ntotxloc, ntotact);
		}
		curcell mvmig;
		if(mvmi){
			mvmig=curcell(NGPU, 1);
		}

		/* Loading or saving intermediate FitR Result.*/
		X(mat)* mvmf=NULL;
		if(parms->load.mvmf){
			/*Load FitR FitL results from file. Resembling warm restart case
			  where mvmf is kept in memory*/
			mvmf=X(read)("%s", parms->load.mvmf);
			if(mvmf->nx!=ntotxloc||mvmf->ny!=ntotact){
				error("loaded mvmf has dimension (%ld, %ld) but we expect (%d, %d)",
					mvmf->nx, mvmf->ny, ntotxloc, ntotact);
			}
		} else if(parms->save.mvmf){
			/*save FitR FitL resutls to file, for later loading.*/
			mvmf=X(new)(ntotxloc, ntotact);
		}
		curcell mvmfg;
		if(mvmf){
			mvmfg=curcell(NGPU, 1);
		}
		if(!parms->load.mvmf){
			dmat* FLId=NULL; /* MI is inv(FL) for direct methods*/
			/*Prepare FitR, FitL is don't load fitting results using mvmf*/
			switch(parms->fit.alg){
			case 0:{//Use CPU to handle CBS.
				dmat* eye=dnew(ntotact, ntotact);
				daddI(eye, 1);
				FLId=dnew(ntotact, ntotact);
				muv_direct_solve_mat(&FLId, &recon->fit->FL, eye);
				dfree(eye);
				toc2("Fit CBS");tic;
			}
				  break;
			case 1://Use GPU.
				break;
			case 2:
				FLId=dref(recon->fit->FL.MI);
				break;
			default:
				error("Invalid fit.alg=%d\n", parms->fit.alg);
			}
			if(FLId){
				FLI.init(ntotact*ntotact, 1);
				for(int i=0; i<ntotact*ntotact; i++){
					FLI[i]=(Real)FLId->p[i];
				}
				dfree(FLId);
			}
		}
		mvm_igpu_t data={parms, recon, mvmig, mvmfg, mvmt, FLI(), residual, residualfit, curp, ntotact, ntotgrad, parms->load.mvmf?1:0};
		int nthread=NGPU;
		thread_t *thdata=thread_prep(0, ntotact, nthread, mvm_trans_igpu, &data);

		/*Initialyze intermediate TomoL result array in GPU. Send intermediate
		  TomoL results to GPU if load.mvmi is set.*/
		if(mvmi){
			tic;
			for(int i=0; i<NGPU; i++){
				gpu_set(i);
				mvmig[i]=curmat(ntotxloc, thdata[i].end-thdata[i].start);
				if(parms->load.mvmi){
					cudaMemcpy(mvmig[i](), mvmi->p+thdata[i].start*ntotxloc,
						sizeof(Real)*ntotxloc*(thdata[i].end-thdata[i].start), H2D);
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
				mvmfg[i]=curmat(ntotxloc, thdata[i].end-thdata[i].start);
				if(parms->load.mvmf){
					cudaMemcpy(mvmfg[i](), mvmf->p+thdata[i].start*ntotxloc,
						sizeof(Real)*ntotxloc*(thdata[i].end-thdata[i].start), H2D);
				}
			}
			if(parms->load.mvmf){
				toc2("copy mvmf to gpu");
			}
		}
		/*Do real MVM control matrix assemble in multiply CPU/GPU*/
		tic;
		CALL_THREAD(thdata, 1);
		toc2("MVM Assembly in GPU");
		

		if(parms->save.setup){
			writebin(residual, "MVM_RL_residual");
			writebin(residualfit, "MVM_FL_residual");
		}
		if(parms->save.mvmi){
			for(int i=0; i<NGPU; i++){
				gpu_set(i);
				cudaMemcpy(mvmi->p+thdata[i].start*ntotxloc, mvmig[i](),
					sizeof(Real)*ntotxloc*(thdata[i].end-thdata[i].start), D2H);
			}
			writebin(mvmi, "MVM_Tomo.bin");
		}
		if(parms->save.mvmf){
			for(int i=0; i<NGPU; i++){
				gpu_set(i);
				cudaMemcpy(mvmf->p+thdata[i].start*ntotxloc, mvmfg[i](),
					sizeof(Real)*ntotxloc*(thdata[i].end-thdata[i].start), D2H);
			}
			writebin(mvmf, "MVM_FitL.bin");
		}
		free(thdata);
		X(free)(mvmi);
		X(free)(mvmf);
		X(free)(residual);
		X(free)(residualfit);
		free(curp);
	}
	{
		for(int i=0; i<NGPU; i++){
			gpu_set(i);
			if(cudata->recon){
				delete cudata->recon;
				cudata->recon=NULL;
			}
		}
	}
	//Clear memory used to do the assembly to avoid out of memory error below.
	{
		gpu_set(cuglobal->recongpu);
		tic;
		gpu_print_mem("before trans");
		long mem=gpu_get_free_mem();
		if((long)(mvmt.N()*sizeof(Real)+100000L)>mem){
			warning("Not enough memory to do the transpose in GPU.\n");
			dmat *MVMt=0;
			cp2cpu(&MVMt, mvmt);
			toc2("Copy to CPU"); tic;
			recon->MVM=dtrans(MVMt);
			toc2("Transpose in CPU"); tic;
			dfree(MVMt);
		}else{
			stream_t stream;
			cuglobal->mvm=mvmt.trans(stream);
			stream.sync();
			toc2("MVM Transpose in GPU");tic;
			cp2cpu(&recon->MVM, cuglobal->mvm, stream);
			stream.sync();
			toc2("MVM copy to CPU");
		}
	}
}

