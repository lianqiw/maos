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

#include <cuda.h>

#include "curmat.h"
#include "cucmat.h"
#include "utils.h"
#include "tomo.h"
#include "fit.h"
#include "recon.h"
#include "cudata.h"
#include "perf.h"
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

/*
  The caller must specify current GPU.
*/
namespace cuda_recon{
curecon_t::curecon_t(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon)
    :grid(0), FR(0), FL(0), RR(0), RL(0), MVM(0),
     gradin(0),opdr(0),opdr_vec(0),opdr_save(0),tomo_rhs(0),
     dmfit(0), dmfit_vec(0),dmfit_save(0),fit_rhs(0),
     RFdfx(0),GXL(0),gngsmvst(0), deltafocus(0),
     moao(0), dm_wfs(0),dm_evl(0){
    if(!parms) return;
    cgstream = new stream_t;
    if(parms->gpu.fit || parms->recon.mvm || parms->gpu.moao){
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
	    return;//done.
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
	    FL=new cusolve_mvm(recon->FL.MI);
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
	    RL=new cusolve_mvm(recon->RL.MI);
	    break;
	default:
	    error("Invalid");
	}
    }
    if((parms->gpu.tomo || parms->gpu.fit) && !parms->sim.idealfit){
	if(parms->tomo.square){
	    opdr=curcellnew(recon->npsr, 1, recon->xnx, recon->xny);
	}else{
	    opdr=curcellnew(recon->npsr, 1, recon->xnloc, (long*)NULL);
	}
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
    if(parms->nmoao){
	nmoao=parms->nmoao;
	const int nwfs=parms->nwfs;
	const int nevl=parms->evl.nevl;
	moao=(cumoao_t**)calloc(nmoao, sizeof(cumoao_t*));
	dm_moao=(curcell***)calloc(nmoao, sizeof(curcell**));
	moao_gwfs=X(new)(nwfs, 1);
	moao_gevl=X(new)(nevl, 1);
	/*Pre-allocate moao output and assign to wfs or evl*/
	for(int imoao=0; imoao<parms->nmoao; imoao++){
	    if(!parms->moao[imoao].used) continue;
	    int ntot=nwfs+nevl;
	    int count=0;
	    dir_t dir[ntot];
	    dm_moao[imoao]=(curcell**)calloc(ntot, sizeof(curcell*));
	    int anx=recon->moao[imoao].amap->nx;
	    int any=recon->moao[imoao].amap->ny;
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].moao==imoao){
		    if(!dm_wfs){
			dm_wfs=curcellnew(nwfs, 1);
		    }
		    dm_moao[imoao][count]=curcellnew(1,1,anx,any);
		    dm_wfs->p[iwfs]=dm_moao[imoao][count]->p[0]->ref();
		    moao_gwfs->p[iwfs]=parms->moao[imoao].gdm;
		    dir[count].thetax=parms->wfs[iwfs].thetax;
		    dir[count].thetay=parms->wfs[iwfs].thetay;
		    dir[count].hs=parms->powfs[ipowfs].hs;
		    count++;
		}
	    }
	    if(parms->evl.moao==imoao){
		if(!dm_evl){
		    dm_evl=curcellnew(nevl, 1);
		}
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		    dm_moao[imoao][count]=curcellnew(1,1, anx, any);
		    dm_evl->p[ievl]=dm_moao[imoao][count]->p[0]->ref();
		    moao_gevl->p[ievl]=parms->moao[imoao].gdm;
		    dir[count].thetax=parms->evl.thetax[ievl];
		    dir[count].thetay=parms->evl.thetay[ievl];
		    dir[count].hs=parms->evl.hs[ievl];
		    count++;
		}
	    }
	    dm_moao[imoao]=(curcell**)realloc(dm_moao[imoao], count*sizeof(curcell*));
	    moao[imoao]=new cumoao_t(parms, recon->moao+imoao, dir, count, grid);
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].moao==imoao){
		    gpu_set(cudata_t::wfsgpu[iwfs]);
		    if(!cudata->dm_wfs){
			cudata->dm_wfs=(cumap_t**)calloc(nwfs, sizeof(cumap_t*));
		    }
		    cudata->dm_wfs[iwfs]=new cumap_t(recon->moao[imoao].amap); 
		}
	    }
	    if(parms->evl.moao==imoao){
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		    gpu_set(cudata_t::evlgpu[ievl]);
		    if(!cudata->dm_evl){
			cudata->dm_evl=(cumap_t**)calloc(parms->evl.nevl, sizeof(cumap_t*));
		    }
		    cudata->dm_evl[ievl]=new cumap_t(recon->moao[imoao].amap); 
		}
	    }
	    gpu_set(gpu_recon);
	}
    }
    gpu_print_mem("recon init");
}

void curecon_t::update(const PARMS_T *parms, RECON_T *recon){
    if(parms->tomo.predict){
	for(int ips=0; ips<recon->npsr; ips++){ 
	    int ips0=parms->atmr.indps[ips]; 
	    grid->xmap[ips].vx=cudata->atm[ips0].vx;
	    grid->xmap[ips].vy=cudata->atm[ips0].vy;
	}
    }
   cutomo_grid *_RR=dynamic_cast<cutomo_grid*>(RR);
   cutomo_grid *_RL=dynamic_cast<cutomo_grid*>(RL);
   if(_RL){
       _RL->init_hx(parms, recon);
   }
   if(_RR && _RR != _RL){
       _RR->init_hx(parms, recon);
   }
   if(parms->tomo.precond==1){
       _RL->update_fdpcg(recon->fdpcg);
   }
}
void curecon_t::reset(const PARMS_T *parms){
    curcellzero(opdr, 0);
    curcellzero(dmfit, 0);
    if(dm_wfs){
	dm_wfs->zero(0);
    }
    if(dm_evl){
	dm_evl->zero(0);
    }
   
}

#define DBG_RECON 1
Real curecon_t::tomo(dcell **_opdr, dcell **_gngsmvst, dcell **_deltafocus,
		      const dcell *_gradin){
    cp2gpu(&gradin, _gradin);
#ifdef DBG_RECON
    curcellcp(&opdr_save, opdr, *cgstream);
#endif
    RR->R(&tomo_rhs, 0, gradin, 1, *cgstream);
    Real cgres=RL->solve(&opdr, tomo_rhs, *cgstream);
#ifdef DBG_RECON
    static Real cgres_last=INFINITY;
    int isimr=grid->isimr;
    if(isimr>5 || cgres>cgres_last*5){
	Real omax=curmax(opdr->m, *cgstream);
	static Real omax_last=INFINITY;
	if(omax>omax_last*5 || cgres>cgres_last*5){
	    info("tomo cgres=%g cgres_last=%g. omax=%g, omax_last=%g\n", cgres, cgres_last, omax, omax_last);
	    if(!disable_save){
		dcellwrite(_gradin, "tomo_gradin_%d", isimr);
		curcellwrite(opdr_save, "tomo_opdrlast_%d", isimr);
		curcellwrite(opdr, "tomo_opdr_%d", isimr);
		curcellwrite(tomo_rhs, "tomo_rhs_%d", isimr);
	    }
	    curcellcp(&opdr, opdr_save, *cgstream);
	    Real newres=RL->solve(&opdr, tomo_rhs, *cgstream);
	    info2("tomo[%d]: omax=%g, oldres=%g. newres=%g\n", 
		  isimr, omax, cgres, newres);
	    if(!disable_save){
		curcellwrite(opdr, "tomo_opdrredo_%d", isimr);
	    }
	    cgres=newres;
	}
	omax_last=omax;
    }
    cgres_last=cgres;
#endif
    if(_opdr){
	cp2cpu(_opdr, opdr_vec, *cgstream);
    }
    if(GXL){
	info("computing ngsmvst\n");
	curcellmm(&gngsmvst, 0, GXL, opdr_vec, "nn", 1, *cgstream);
	add2cpu(_gngsmvst, 1, gngsmvst, 1, *cgstream);
    }
    if(RFdfx){
	info2("computing deltafocus\n");
	curcellmm(&deltafocus, 0, RFdfx, opdr_vec, "nn", 1, *cgstream);
	cp2cpu(_deltafocus, deltafocus, *cgstream);
    }
    cgstream->sync();
return cgres;
}
Real curecon_t::fit(dcell **_dmfit, dcell *_opdr){
    if(_opdr){
	cp2gpu(&opdr_vec, _opdr);
    }
#ifdef DBG_RECON
    curcellcp(&dmfit_save, dmfit, *cgstream);
#endif
    FR->R(&fit_rhs, 0, opdr, 1, *cgstream);
    Real cgres=FL->solve(&dmfit, fit_rhs, *cgstream);
#ifdef DBG_RECON
    static Real cgres_last=INFINITY;
    int isimr=grid->isimr;
    if(cgres>cgres_last*5){
	info("fit cgres=%g\n", cgres);
	if(!disable_save){
	    curcellwrite(opdr, "fit_opdr_%d", isimr);
	    curcellwrite(dmfit_save, "fit_dmfitlast_%d", isimr);
	    curcellwrite(dmfit, "fit_dmfit_%d", isimr);
	    curcellwrite(fit_rhs, "fit_rhs_%d", isimr);
	}
	curcellcp(&dmfit, dmfit_save, *cgstream);
	Real newres=FL->solve(&dmfit, fit_rhs, *cgstream);
	info2("fit[%d]: oldres=%g. newres=%g\n", isimr, cgres, newres);
	if(!disable_save){
	    curcellwrite(dmfit, "fit_dmfitredo_%d", isimr);
	}
	cgres=newres;
    }
    cgres_last=cgres;
#endif
    add2cpu(_dmfit, 0, dmfit_vec, 1, *cgstream);
    cgstream->sync();
    return cgres;
}
Real curecon_t::moao_recon(dcell *_dmfit, dcell *_opdr){
    if(_dmfit){
	cp2gpu(&dmfit_vec, _dmfit);
    }
    if(_opdr){
	cp2gpu(&opdr_vec, _opdr);
    }
    for(int imoao=0; imoao<nmoao; imoao++){
	if(!moao[imoao]) continue;
	moao[imoao]->moao_solve(dm_moao[imoao], opdr, dmfit, *cgstream);
    }
    return 0;
}

void curecon_t::moao_filter(dcell *_dm_wfs, dcell *_dm_evl){
    warning_once("MOAO temporal filter implemented with LPF\n");
    const int *wfsgpu=cudata_t::wfsgpu;
    if(dm_wfs){
	int nwfs=dm_wfs->nx;
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    if(!dm_wfs->p[iwfs]) continue;
	    curmat *temp=NULL;
	    if(wfsgpu) gpu_set(wfsgpu[iwfs]);
	    stream_t stream;
	    if(wfsgpu && wfsgpu[iwfs] != gpu_recon){
		curcp(&temp, dm_wfs->p[iwfs], stream);
	    }else{
		temp=dm_wfs->p[iwfs]->ref();
	    }
	    Real g=moao_gwfs->p[iwfs];
	    curadd(&cudata->dm_wfs[iwfs]->p, 1-g, temp, g, stream);
	    if(!wfsgpu){//gpu.moao implies fit.square=1
		cp2cpu(&_dm_wfs->p[iwfs], cudata->dm_wfs[iwfs]->p, stream);
	    }
	    delete temp;
	}
    }

    if(dm_evl){
	const int nevl=dm_evl->nx;
	for(int ievl=0; ievl<nevl; ievl++){
	    if(cudata_t::evlgpu) gpu_set(cudata_t::evlgpu[ievl]);
	    stream_t stream;
	    curmat *temp=0;
	    if(cudata_t::evlgpu && cudata_t::evlgpu[ievl]!=gpu_recon){
		curcp(&temp, dm_evl->p[ievl], stream);
	    }else{
		temp=dm_evl->p[ievl]->ref();
	    }
	    Real g=moao_gevl->p[ievl];
	    curadd(&cudata->dm_evl[ievl]->p, 1-g, temp, g, stream);
	    if(!cudata_t::evlgpu){
		cp2cpu(&_dm_evl->p[ievl], cudata->dm_evl[ievl]->p, stream);
	    }
	    delete temp;
	}
    }
}
void curecon_t::mvm(dcell **_dmerr, dcell *_gradin){
    cp2gpu(&gradin, _gradin);
    MVM->solve(&dmfit, gradin, *cgstream);
    cp2cpu(_dmerr, dmfit_vec, *cgstream);
    cgstream->sync();
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
    if(parms->tomo.alg==1){
	muv(&lc, &recon->RL, rhsc, 1);
	dcellwrite(lc, "CPU_TomoL");
	muv(&lc, &recon->RL, rhsc, -1);
	dcellwrite(lc, "CPU_TomoL2");
	if(parms->tomo.precond==1){
	    dcell *lp=NULL;
	    fdpcg_precond(&lp, recon, rhsc);
	    dcellwrite(lp, "CPU_TomoP");
	    fdpcg_precond(&lp, recon, rhsc);
	    dcellwrite(lp, "CPU_TomoP2");
	    dcellfree(lp);
	}
    }
    dcellzero(lc);
    for(int i=0; i<5; i++){
	muv_solve(&lc, &recon->RL, NULL, rhsc);
	dcellwrite(lc, "CPU_TomoCG%d", i);
    }
	
    cp2gpu(&gradin, simu->gradlastol);
    RR->R(&rhsg, 0, gradin, 1, stream);
    curcellwrite(rhsg, "GPU_TomoR");
    RR->Rt(&rtg, 0, rhsg, 1, stream);
    curcellwrite(rtg, "GPU_TomoRt");
    if(parms->tomo.alg==1){
	cucg_t *RL2=dynamic_cast<cucg_t*>(RL);
	RL2->L(&lg, 0, rhsg, 1,stream);
	curcellwrite(lg, "GPU_TomoL"); 
	RL2->L(&lg, 1, rhsg, -1,stream);
	curcellwrite(lg, "GPU_TomoL2");
	if(parms->tomo.precond==1){
	    cucg_t *RL2=dynamic_cast<cucg_t*>(RL);
	    curcell *lp=NULL;
	    RL2->P(&lp, rhsg, stream);
	    curcellwrite(lp, "GPU_TomoP");
	    RL2->P(&lp, rhsg, stream);
	    curcellwrite(lp, "GPU_TomoP2");
	    delete lp;
	}
    }
    curcellzero(lg, stream);
    for(int i=0; i<5; i++){
	RL->solve(&lg, rhsg, stream);
	curcellwrite(lg, "GPU_TomoCG%d", i);
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
	cp2cpu(&simu->opdr, opdr_vec, 0);
    }
    if(!simu->parms->gpu.tomo){
	cp2gpu(&opdr_vec, simu->opdr);
    }
    dcellwrite(simu->opdr, "opdr");
    muv(&rhsc, &recon->FR, simu->opdr, 1);
    dcellwrite(rhsc, "CPU_FitR");
    muv(&lc, &recon->FL, rhsc, 1);
    dcellwrite(lc, "CPU_FitL");
    muv(&lc, &recon->FL, rhsc, -1);
    dcellwrite(lc, "CPU_FitL2");
    dcellzero(lc);
    for(int i=0; i<5; i++){
	muv_solve(&lc, &recon->FL, NULL, rhsc);
	dcellwrite(lc, "CPU_FitSolve%d", i);
    }
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
    curcellzero(lg, stream);
    for(int i=0; i<5; i++){
	FL->solve(&lg, rhsg, stream);
	curcellwrite(lg, "GPU_FitSolve%d", i);
    }
    curcell *lhsg=NULL;
    FR->Rt(&lhsg, 0, rhsg, 1, stream);
    curcellwrite(lhsg, "GPU_FitRt");
    CUDA_SYNC_DEVICE;
    exit(0);
}
}//namespace
using cuda_recon::curecon_t;
void gpu_setup_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    for(int igpu=0; igpu<NGPU; igpu++){
	if((parms->recon.mvm && parms->gpu.tomo && parms->gpu.fit && !parms->load.mvm)
	   || igpu==gpu_recon){
	    gpu_set(igpu);
	    if(cudata->recon){
		error("Already initialized\n");
	    }
	    cudata->recon=new curecon_t(parms, powfs, recon);
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
	    curecon->update(parms, recon);
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
		gpu_set(cudata_t::wfsgpu[iwfs]);
		if(cudata->dm_wfs[iwfs]){
		    cudata->dm_wfs[iwfs]->p->zero();
		}
	    }
	}
	if(cudata->dm_evl){
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		gpu_set(cudata_t::evlgpu[ievl]);
		if(cudata->dm_evl[ievl]){
		    cudata->dm_evl[ievl]->p->zero();
		}
	    }
	}
	CUDA_SYNC_DEVICE;
    }
}
void gpu_tomo(SIM_T *simu){
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    curecon->grid->isim=simu->isim;
    curecon->grid->isimr=simu->reconisim;
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    if(parms->dbg.tomo){
	curecon->tomo_test(simu);
    }else{
	int copy2cpu=(parms->plot.run>1 || !parms->gpu.fit || parms->save.opdr || (recon->moao && !parms->gpu.moao));
	simu->cgres->p[0]->p[simu->reconisim]=
	    curecon->tomo(copy2cpu?&simu->opdr:NULL, &simu->gngsmvst, &simu->deltafocus,
			  parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
    }
}

void gpu_fit(SIM_T *simu){
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    curecon->grid->isim=simu->isim;
    curecon->grid->isimr=simu->reconisim;
    const PARMS_T *parms=simu->parms;
    if(parms->dbg.fit){
	curecon->fit_test(simu);
    }else{
	simu->cgres->p[1]->p[simu->reconisim]=
	    curecon->fit(&simu->dmfit, parms->gpu.tomo?NULL:simu->opdr);
    }
    /*Don't free opdr. Needed for warm restart in tomo.*/
}
void gpu_recon_mvm(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    curecon->mvm(&simu->dmerr, parms->tomo.psol?simu->gradlastol:simu->gradlastcl);
}

void gpu_moao_recon(SIM_T *simu){
    gpu_set(gpu_recon);
    const PARMS_T *parms=simu->parms;
    curecon_t *curecon=cudata->recon;
    curecon->moao_recon(parms->gpu.fit?NULL:simu->dmfit, (parms->gpu.tomo || parms->gpu.fit)?NULL:simu->opdr);
}

void gpu_moao_filter(SIM_T *simu){
    gpu_set(gpu_recon);
    curecon_t *curecon=cudata->recon;
    curecon->moao_filter(simu->dm_wfs, simu->dm_evl);
}
