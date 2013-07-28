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
#include "curmat.h"
#include "cucmat.h"
#include "utils.h"
#include "tomo.h"
#include "fit.h"
#include "recon.h"
#include "cudata.h"
#undef TIMING
#define TIMING 0
#if !TIMING
#define TIC_test
#define tic_test
#define toc_test(A)
#else
#define TIC_test TIC
#define tic_test tic
#define toc_test(A) toc2(A);tic
#endif
#define CATCH_ERR 0
/*__global__ static void saloc2ptr_do(int (*restrict saptr)[2], float (*restrict saloc)[2], 
  int nsa, float ox, float oy, float dx, float dy){
  const int step=blockDim.x * gridDim.x;
  const float dx1=1./dx;
  const float dy1=1./dy;
  for(int isa=blockIdx.x * blockDim.x + threadIdx.x; isa<nsa; isa+=step){
  saptr[isa][0]=(int)roundf((saloc[isa][0]-ox)*dx1);
  saptr[isa][1]=(int)roundf((saloc[isa][1]-oy)*dy1);
  }
  }*/

/*
  The caller must specify current GPU.
*/
void curecon_t::init(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    if(parms->gpu.fit || parms->recon.mvm){
	if(parms->fit.square){
	    dmfit=curcellnew(parms->ndm, 1, recon->anx, recon->any);
	}else{
	    dmfit=curcellnew(parms->ndm, 1, recon->anloc, (long*)NULL);
	}
	dmfit_vec=curcellnew(parms->ndm, 1);
	for(int idm=0; idm<parms->ndm; idm++){
	    dmfit_vec->p[idm]=dmfit->p[idm]->ref(1);
	}
    }
    if(parms->recon.mvm){
	if(!parms->gpu.tomo || !parms->gpu.fit){
	    return; /*Use CPU to assemble MVM*/
	}
	if(recon->MVM){//MVM already exists
	    MVM=new cusolve_mvm(recon->MVM);
	}
    }
    grid=new curecon_geom(parms, recon);
    if(parms->recon.alg!=0){
	error("Only MVR is implemented in GPU\n");
    }
    if(parms->gpu.fit){
	switch(parms->gpu.fit){
	case 1:
	    FR=new cusolve_sparse(parms->fit.maxit, parms->recon.warm_restart,
				  &recon->FR, &recon->FL);
	    break;
	case 2:
	    FR=new cufit_grid(parms, recon, grid);
	    break;
	default:error("Invalid");
	}
	switch(parms->fit.alg){
	case 0:
	    FL=new cusolve_cbs(recon->FL.C, recon->FL.Up, recon->FL.Vp);
	    break;
	case 1:
	    FL=dynamic_cast<cusolve_l*>(FR);
	    break;
	case 2:
	    FL=new cusolve_svd(recon->FL.MI);
	    break;
	default:
	    error("Invalid");
	}
    }
    if(parms->gpu.tomo){
	RR=new cutomo_grid(parms, recon, powfs, grid);
	switch(parms->tomo.alg){
	case 0:
	    RL=new cusolve_cbs(recon->RL.C, recon->RL.Up, recon->RL.Vp);
	    break;
	case 1:
	    RL=dynamic_cast<cusolve_l*>(RR);
	    break;
	case 2:
	    RL=new cusolve_svd(recon->RL.MI);
	    break;
	default:
	    error("Invalid");
	}
    }
    if((parms->gpu.tomo || parms->gpu.fit) && !parms->sim.idealfit){
	opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
	opdr_vec=curcellnew(recon->npsr, 1);
	for(int ips=0; ips<recon->npsr; ips++){
	    opdr_vec->p[ips]=opdr->p[ips]->ref(1);
	}
    }
 
    if(parms->dbg.deltafocus){
	cp2gpu(&RFdfx, recon->RFdfx);
    }
    if(parms->recon.split==2){
	cp2gpu(&GXL, recon->GXL);
    }
    gpu_print_mem("recon init");
}
void gpu_setup_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    for(int igpu=0; igpu<NGPU; igpu++){
	if((parms->recon.mvm && parms->gpu.tomo && parms->gpu.fit && !parms->load.mvm)
	   || igpu==gpu_recon){
	    gpu_set(igpu);
	    if(!cudata->recon){
		cudata->recon=new curecon_t(parms, powfs, recon);
	    }
	}
    }
}
void gpu_recon_free(){
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	if(cudata->recon){
	    delete cudata->recon;
	    cudata->recon=NULL;
	}
    }
}
void gpu_setup_recon_mvm(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    /*The following routine assemble MVM and put in recon->MVM*/
    if(!parms->load.mvm){
	if(parms->recon.mvm==1){
	    gpu_setup_recon_mvm_trans(parms, recon, powfs);
	}else{
	    gpu_setup_recon_mvm_direct(parms, recon, powfs);
	}
    }
    //free existing data
    gpu_recon_free();
    if(!parms->sim.mvmport){
	gpu_set(gpu_recon);
	//recreate curecon_t that uses MVM.
	cudata->recon=new curecon_t(parms, powfs, recon);
    }
    gpu_print_mem("MVM");
}
void gpu_update_recon(const PARMS_T *parms, RECON_T *recon){
    if(!parms->gpu.tomo) return;
    for(int igpu=0; igpu<NGPU; igpu++){
	if((parms->recon.mvm && parms->gpu.tomo && parms->gpu.fit && !parms->load.mvm)
	   || igpu==gpu_recon){
	    gpu_set(igpu);

	    curecon_t *curecon=cudata->recon;
	    curecon_geom *grid=curecon->grid;
	    if(parms->tomo.predict){
		if(!grid->vx){
		    grid->vx=new float[recon->npsr];
		    grid->vy=new float[recon->npsr];
		}
		for(int ips=0; ips<recon->npsr; ips++){ 
		    int ips0=parms->atmr.indps[ips]; 
		    grid->vx[ips]=cudata->atm[ips0].vx;
		    grid->vy[ips]=cudata->atm[ips0].vy;
		}
	    }
	    cutomo_grid *RR=dynamic_cast<cutomo_grid*>(curecon->RR);
	    cutomo_grid *RL=dynamic_cast<cutomo_grid*>(curecon->RL);
	    RL->init_hxdata(parms, recon);
	    if(RR && RR !=RL){
		RR->init_hxdata(parms, recon);
	    }
	    RL->init_fdpcg(recon->fdpcg, grid);
	}
    }
}

void curecon_t::reset(const PARMS_T *parms){
    curcellzero(opdr, 0);
    curcellzero(dmfit, 0);
    if(dm_wfs){
	for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
	    curcellzero(dm_wfs[iwfs], 0);
	}
    }
    if(dm_evl){
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    curcellzero(dm_evl[ievl], 0);
	}
    }
   
}
void gpu_recon_reset(const PARMS_T *parms){/*reset warm restart.*/
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	curecon_t *curecon=cudata->recon;
	if(curecon){
	    curecon->reset(parms);
	}
	if(cudata->dm_wfs){
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		cudata->dm_wfs[iwfs].p.zero();
	    }
	}
	if(cudata->dm_evl){
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		cudata->dm_evl[ievl].p.zero();
	    }
	}
	CUDA_SYNC_DEVICE;
    }
}
#define DBG_RECON 1
float curecon_t::tomo(dcell **_opdr, dcell **_gngsmvst, dcell **_deltafocus,
		      const dcell *_gradin){
    cp2gpu(&gradin, _gradin);
#ifdef DBG_RECON
    curcellcp(&opdr_save, opdr, cgstream);
#endif
    float cgres=RL->solve(&opdr, gradin, RR, cgstream);
#ifdef DBG_RECON
    int isimr=grid->isimr;
    if(isimr>5){
	float omax=curmax(opdr->m, cgstream);
	if(omax>4e-6){
	    dcellwrite(_gradin, "tomo_gradin_%d", isimr);
	    curcellwrite(opdr_save, "tomo_opdrlast_%d", isimr);
	    curcellwrite(opdr, "tomo_opdr_%d", isimr);
	    curcellcp(&opdr, opdr_save, cgstream);
	    float newres=RL->solve(&opdr, gradin, RR, cgstream);
	    info2("tomo[%d]: oldres=%g. newres=%g\n", isimr, cgres, newres);
	    curcellwrite(opdr, "tomo_opdrredo_%d", isimr);
	    cgres=newres;
	}
    }
#endif
    if(_opdr){
	cp2cpu(_opdr, 0, opdr_vec, 1, cgstream);
    }
    if(GXL){
	curcellmm(&gngsmvst, 0, GXL, opdr_vec, "nn", 1, cgstream);
	add2cpu(_gngsmvst, gngsmvst, cgstream);
    }
    if(RFdfx){
	curcellmm(&deltafocus, 0, RFdfx, opdr_vec, "nn", 1, cgstream);
	cp2cpu(_deltafocus, deltafocus, cgstream);
    }
    cgstream.sync();
    return cgres;
}
float curecon_t::fit(dcell **_dmfit, dcell *_opdr){
    if(_opdr){
	cp2gpu(&opdr_vec, _opdr);
    }
#ifdef DBG_RECON
    curcellcp(&dmfit_save, dmfit, cgstream);
#endif
    float cgres=FL->solve(&dmfit, opdr, FR, cgstream);
#ifdef DBG_RECON
    int isimr=grid->isimr;
    if(isimr > 5 && cgres>1e-2){
        curcellwrite(opdr, "fit_opdr_%d", isimr);
	curcellwrite(dmfit_save, "fit_dmfitlast_%d", isimr);
	curcellwrite(dmfit, "fit_dmfit_%d", isimr);
	curcellcp(&dmfit, dmfit_save, cgstream);
	float newres=FL->solve(&dmfit, opdr, FR, cgstream);
	info2("fit[%d]: oldres=%g. newres=%g\n", isimr, cgres, newres);
	curcellwrite(dmfit, "fit_dmfitredo_%d", isimr);
	cgres=newres;
    }
#endif
    cp2cpu(_dmfit, 0, dmfit_vec, 1, cgstream);
    cgstream.sync();
    return cgres;
}
void curecon_t::mvm(dcell **_dmerr, dcell *_gradin){
    cp2gpu(&gradin, _gradin);
    MVM->solve(&dmfit, gradin, NULL, cgstream);
    cp2cpu(_dmerr, 0, dmfit_vec, 1, cgstream);
    cgstream.sync();
}
void gpu_tomo(SIM_T *simu){
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    curecon->grid->isim=simu->isim;
    curecon->grid->isimr=simu->reconisim;
    TIC_test;tic_test;
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
#if 0
    gpu_tomo_test(simu);
#endif
    toc_test("Before gradin");
    int copy2cpu=(!parms->gpu.fit || parms->save.opdr || (recon->moao && !parms->gpu.moao));
    simu->cgres->p[0]->p[simu->reconisim]=
	curecon->tomo(copy2cpu?&simu->opdr:NULL, &simu->gngsmvst, &simu->deltafocus,
		      parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
    toc_test("Tomo");
}

void gpu_fit(SIM_T *simu){
    gpu_set(gpu_recon);
    TIC_test;tic_test;
    curecon_t *curecon=cudata->recon;
    curecon->grid->isim=simu->isim;
    curecon->grid->isimr=simu->reconisim;
    const PARMS_T *parms=simu->parms;
#if 0
    curecon->fit_test(simu);
#endif
    toc_test("Before FitR");

    simu->cgres->p[1]->p[simu->reconisim]=
	curecon->fit(&simu->dmfit, parms->gpu.tomo?NULL:simu->opdr);
    /*Don't free opdr. Needed for warm restart in tomo.*/
    toc_test("Fit");
}
void gpu_recon_mvm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    curecon->mvm(&simu->dmerr, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
}

void curecon_t::tomo_test(SIM_T *simu){
    gpu_set(gpu_recon);
    const PARMS_T *parms=simu->parms;
    stream_t stream;
    RECON_T *recon=simu->recon;
    /*Debugging. */
    dcell *rhsc=NULL;
    dcell *lc=NULL;
    dcell *rtc=NULL;
    curcell *rhsg=NULL;
    curcell *lg=NULL;
    curcell *rtg=NULL;
    muv(&rhsc, &recon->RR, simu->gradlastol, 1);
    dcellwrite(rhsc, "CPU_TomoR");
    muv_trans(&rtc, &recon->RR, rhsc, 1);
    dcellwrite(rtc, "CPU_TomoRt");
    muv(&lc, &recon->RL, rhsc, 1);
    dcellwrite(lc, "CPU_TomoL");
    muv(&lc, &recon->RL, rhsc, -1);
    dcellwrite(lc, "CPU_TomoL2");
    if(parms->tomo.alg==1 && parms->tomo.precond==1){
	dcell *lp=NULL;
	fdpcg_precond(&lp, recon, lc);
	dcellwrite(lp, "CPU_TomoP");
	fdpcg_precond(&lp, recon, lc);
	dcellwrite(lp, "CPU_TomoP2");
    }
    dcellzero(lc);
    for(int i=0; i<2; i++){
	muv_solve(&lc, &recon->RL, NULL, rhsc);
	dcellwrite(lc, "CPU_Tomo_%d", i);
    }
	
    cp2gpu(&gradin, simu->gradlastol);
    RR->R(&rhsg, 0, gradin, 1, stream);
    curcellwrite(rhsg, "GPU_TomoR");
    RR->Rt(&rtg, 0, rhsg, 1, stream);
    curcellwrite(rtg, "GPU_TomoRt");
    if(parms->tomo.alg==1){
	cucg_t *RL2=dynamic_cast<cucg_t*>(RL);
	RL2->L(&lg, 0, rhsg, 1,stream);
	curcellwrite(lg, "GPU_TomoL1"); 
	RL2->L(&lg, 1, rhsg, -1,stream);
	curcellwrite(lg, "GPU_TomoL2");
    }
    /*if(parms->tomo.alg==1 && parms->tomo.precond==1){
      curcell *lp=NULL;
      gpu_Tomo_fdprecond(&lp, recon, lg, stream);
      curcellwrite(lp, "GPU_TomoP");
      gpu_Tomo_fdprecond(&lp, recon, lg, stream);
      curcellwrite(lp, "GPU_TomoP2");
      }*/
    for(int i=0; i<2; i++){
	curcellzero(lg, stream);
	RL->solve(&lg, rhsg, NULL, stream);
	curcellwrite(lg, "GPU_Tomo_%d", i);
    }
    CUDA_SYNC_DEVICE;
    exit(0);
}

void curecon_t::fit_test(SIM_T *simu){	/*Debugging. */
    gpu_set(gpu_recon);
    stream_t stream;
    const RECON_T *recon=simu->recon;
    dcell *rhsc=NULL;
    dcell *lc=NULL;	
    if(!simu->opdr){
	cp2cpu(&simu->opdr, 0, opdr_vec, 1, 0);
    }
    dcellwrite(simu->opdr, "opdr");
    muv(&rhsc, &recon->FR, simu->opdr, 1);
    dcellwrite(rhsc, "CPU_FitR");
    muv(&lc, &recon->FL, rhsc, 1);
    dcellwrite(lc, "CPU_FitL");
    muv(&lc, &recon->FL, rhsc, -1);
    dcellwrite(lc, "CPU_FitL2");
    dcellzero(lc);
    muv_solve(&lc, &recon->FL, NULL, rhsc);
    dcellwrite(lc, "CPU_FitCG");
    dcell *lhs=NULL;
    muv_trans(&lhs, &recon->FR, rhsc, 1);
    dcellwrite(lhs, "CPU_FitRt");

    curcell *rhsg=NULL;
    curcell *lg=NULL;
    FR->R(&rhsg, 0.f, opdr, 1.f, stream);
    curcellwrite(rhsg, "GPU_FitR");
    cucg_t *FL2=dynamic_cast<cucg_t*>(FL);
    if(FL2){
	FL2->L(&lg, 0, rhsg, 1, stream);
	curcellwrite(lg, "GPU_FitL");
	FL2->L(&lg, 1, rhsg, -1, stream);
	curcellwrite(lg, "GPU_FitL2");
    }
    curcell *lhsg=NULL;
    FR->Rt(&lhsg, 0, rhsg, 1, stream);
    curcellwrite(lhsg, "GPU_FitRt");
    curcellzero(lg, stream);
    FL->solve(&lg, rhsg, NULL, stream);
    curcellwrite(lg, "GPU_FitCG");
    CUDA_SYNC_DEVICE;
    exit(0);
}
