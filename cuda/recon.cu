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

#include "curmat.h"
#include "cucmat.h"
#include "utils.h"
#include "tomo.h"
#include "fit.h"
#include "recon.h"
#include "cudata.h"
#include "perf.h"
#if !USE_CPP
extern "C"{
#endif
#include "../maos/fdpcg.h"
#if !USE_CPP
}
#endif
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
	:grid(0),FR(0), FL(0), RR(0), RL(0), MVM(0), nmoao(0),moao(0){
    if(!parms) return;
    if((parms->recon.alg==0 && parms->gpu.fit || parms->recon.mvm ||parms->gpu.moao)
       || (parms->recon.alg==1 && parms->gpu.lsr)
	){
	if(parms->recon.alg==0 && parms->fit.square){
	    dmfit=curcell(parms->ndm, 1, recon->anx->p, recon->any->p);
	}else if(parms->recon.modal){
	    dmfit=curcell(parms->ndm, 1, recon->anmod->p, (long*)NULL);
	}else{
	    dmfit=curcell(parms->ndm, 1, recon->anloc->p, (long*)NULL);
	}
	dmfit_vec=curcell(parms->ndm, 1);
	for(int idm=0; idm<parms->ndm; idm++){
	    dmfit_vec[idm]=dmfit[idm].Vector();
	}
    }

    if(parms->recon.alg==0 && (parms->gpu.tomo || parms->gpu.fit)
       && !parms->sim.idealfit && !parms->load.mvm){
	if(parms->tomo.square){
	    opdr=curcell(recon->npsr, 1, recon->xnx->p, recon->xny->p);
	}else{
	    opdr=curcell(recon->npsr, 1, recon->xnloc->p, (long*)NULL);
	}
	opdr_vec=curcell(recon->npsr, 1);
	for(int ips=0; ips<recon->npsr; ips++){
	    opdr_vec[ips]=opdr[ips].Vector();
	}
    }
    if(parms->recon.alg==0 && !(parms->recon.mvm && parms->load.mvm)){
	grid=new curecon_geom(parms, recon);//does not change
    }
    update(parms, powfs, recon);
    gpu_print_mem("recon init");
}
void curecon_t::reset_config(){
    if(FL && FL!=dynamic_cast<cusolve_l*>(FR)) delete FL; FL=0;
    delete FR; FR=0;
    if(RL && RL!=dynamic_cast<cusolve_l*>(RR)) delete RL; RL=0;
    delete RR; RR=0;//problematic.
    delete MVM; MVM=0;
    if(moao){
	for(int im=0; im<nmoao; im++){
	    delete moao[im];
	}
	delete [] moao;
	moao=0; nmoao=0;
    }
}
void curecon_t::update(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    reset_config();
    if(parms->recon.mvm){
	if(recon->MVM){//MVM already exists
	    MVM=new cusolve_mvm(recon->MVM);
	    return;//done.
	}else if(!parms->gpu.tomo || !parms->gpu.fit){
	    return; //Use CPU to assemble MVM
	}
    }
    if(parms->gpu.fit){
	switch(parms->gpu.fit){
	case 1:
	    FR=new cusolve_sparse(parms->fit.maxit, parms->recon.warm_restart,
				  &recon->fit->FR, &recon->fit->FL);
	    break;
	case 2:
	    FR=new cufit_grid(parms, recon, grid);
	    break;
	default:
	    error("Invalid");
	}
	switch(parms->fit.alg){
	case 0:
	    FL=new cusolve_cbs(recon->fit->FL.C, recon->fit->FL.Up, recon->fit->FL.Vp);
	    break;
	case 1:
	    FL=dynamic_cast<cusolve_l*>(FR);
	    break;
	case 2:
	    FL=new cusolve_mvm(recon->fit->FL.MI);
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
 
    if(parms->recon.split==2){
	cp2gpu(GXL, recon->GXL);
    }
    if(parms->nmoao){
	nmoao=parms->nmoao;
	const int nwfs=parms->nwfs;
	const int nevl=parms->evl.nevl;
	moao=new cumoao_t*[nmoao];
	dm_moao=curcccell(nmoao, 1);
	moao_gwfs=X(new)(nwfs, 1);
	moao_gevl=X(new)(nevl, 1);
	//Pre-allocate moao output and assign to wfs or evl
	for(int imoao=0; imoao<parms->nmoao; imoao++){
	    if(!parms->moao[imoao].used) continue;
	    int ntot=nwfs+nevl;
	    int count=0;
	    dir_t* dir=new dir_t[ntot];
	    dm_moao[imoao]=curccell(ntot, 1);
	    int anx=recon->moao[imoao].amap->p[0]->nx;
	    int any=recon->moao[imoao].amap->p[0]->ny;
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].moao==imoao){
		    if(!dm_wfs){
			dm_wfs=curcell(nwfs, 1);
		    }
		    dm_wfs[iwfs]=dm_moao[imoao][count][0];
		    dm_moao[imoao][count]=curcell(1,1,anx,any);
		    moao_gwfs->p[iwfs]=parms->moao[imoao].gdm;
		    dir[count].thetax=parms->wfs[iwfs].thetax;
		    dir[count].thetay=parms->wfs[iwfs].thetay;
		    dir[count].hs=parms->wfs[iwfs].hs;
		    dir[count].skip=0;
		    count++;
		}
	    }
	    if(parms->evl.moao==imoao){
		if(!dm_evl){
		    dm_evl=curcell(nevl, 1);
		}
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		    dm_moao[imoao][count]=curcell(1,1, anx, any);
		    dm_evl[ievl]=dm_moao[imoao][count][0];
		    moao_gevl->p[ievl]=parms->moao[imoao].gdm;
		    dir[count].thetax=parms->evl.thetax->p[ievl];
		    dir[count].thetay=parms->evl.thetay->p[ievl];
		    dir[count].hs=parms->evl.hs->p[ievl];
		    dir[count].skip=0;
		    count++;
		}
	    }
	    moao[imoao]=new cumoao_t(parms, recon->moao+imoao, dir, count, grid);
	    delete[] dir;
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs=parms->wfs[iwfs].powfs;
		if(parms->powfs[ipowfs].moao==imoao){
		    gpu_set(cuglobal->wfsgpu[iwfs]);
		    if(!cudata->dm_wfs){
			cudata->dm_wfs=Array<cumapcell>(nwfs, 1);
		    }
		    cudata->dm_wfs[iwfs]=cumapcell(1,1);
		    cudata->dm_wfs[iwfs][0]=(recon->moao[imoao].amap->p[0]); 
		}
	    }
	    if(parms->evl.moao==imoao){
		for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		    gpu_set(cuglobal->evlgpu[ievl]);
		    if(!cudata->dm_evl){
			cudata->dm_evl=Array<cumapcell>(parms->evl.nevl, 1);
		    }
		    cudata->dm_evl[ievl]=cumapcell(1,1);
		    cudata->dm_evl[ievl][0]=(recon->moao[imoao].amap->p[0]); 
		}
	    }
	    gpu_set(cuglobal->recongpu);
	}
    }
}
void curecon_t::update_cn2(const PARMS_T *parms, RECON_T *recon){
    if(parms->tomo.predict){
	for(int ips=0; ips<recon->npsr; ips++){ 
	    int ips0=parms->atmr.indps->p[ips]; 
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
void curecon_t::reset_runtime(){
    opdr.zero(0);
    dmfit.zero(0);
    dm_wfs.zero(0);
    dm_evl.zero(0);
}
    
#define DBG_RECON 0
Real curecon_t::tomo(dcell **_opdr, dcell **_gngsmvst, 
		     const dcell *_gradin){
    cp2gpu(gradin, _gradin);

    RR->R(tomo_rhs, 0, gradin, 1, cgstream);
    Real cgres=RL->solve(opdr, tomo_rhs, cgstream);

    if(_opdr){
	cp2cpu(_opdr, opdr_vec, cgstream);
    }
    if(GXL){
	dbg("computing ngsmvst\n");
	curcellmm(gngsmvst, 0, GXL, opdr_vec, "nn", 1, cgstream);
	add2cpu(_gngsmvst, 1, gngsmvst, 1, cgstream);
    }
    cgstream.sync();
    return cgres;
}
Real curecon_t::fit(dcell **_dmfit, dcell *_opdr){
    if(_opdr){
	cp2gpu(opdr_vec, _opdr);
    }

    FR->R(fit_rhs, 0, opdr, 1, cgstream);
    Real cgres=FL->solve(dmfit, fit_rhs, cgstream);

    add2cpu(_dmfit, 0, dmfit_vec, 1, cgstream);
    cgstream.sync();
    return cgres;
}
Real curecon_t::moao_recon(dcell *_dmfit, dcell *_opdr){
    if(_dmfit){
	cp2gpu(dmfit_vec, _dmfit);
    }
    if(_opdr){
	cp2gpu(opdr_vec, _opdr);
    }
    for(int imoao=0; imoao<nmoao; imoao++){
	if(!moao[imoao]) continue;
	moao[imoao]->moao_solve(dm_moao[imoao], opdr, dmfit, cgstream);
    }
    return 0;
}

void curecon_t::moao_filter(dcell *_dm_wfs, dcell *_dm_evl){
    warning_once("MOAO temporal filter implemented with LPF\n");
    const int *wfsgpu=cuglobal->wfsgpu();
    if(dm_wfs){
	int nwfs=dm_wfs.Nx();
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    if(!dm_wfs[iwfs]) continue;
	    curmat temp;
	    if(wfsgpu) gpu_set(wfsgpu[iwfs]);
	    stream_t stream;
	    if(wfsgpu && wfsgpu[iwfs] != cuglobal->recongpu){
		curcp(temp, dm_wfs[iwfs], stream);
	    }else{
		temp=dm_wfs[iwfs];
	    }
	    Real g=moao_gwfs->p[iwfs];
	    curadd(cudata->dm_wfs[iwfs][0], 1-g, temp, g, stream);
	    if(!wfsgpu || (_dm_wfs && _dm_wfs->p[iwfs])){//gpu.moao implies fit.square=1
		cp2cpu(&_dm_wfs->p[iwfs], cudata->dm_wfs[iwfs][0], stream);
	    }
	}
    }

    if(dm_evl){
	const int nevl=dm_evl.Nx();
	for(int ievl=0; ievl<nevl; ievl++){
	    if(cuglobal->evlgpu) gpu_set(cuglobal->evlgpu[ievl]);
	    stream_t stream;
	    curmat temp;
	    if(cuglobal->evlgpu && cuglobal->evlgpu[ievl]!=cuglobal->recongpu){
		curcp(temp, dm_evl[ievl], stream);
	    }else{
		temp=dm_evl[ievl];
	    }
	    Real g=moao_gevl->p[ievl];
	    curadd(cudata->dm_evl[ievl][0], 1-g, temp, g, stream);
	    if(!cuglobal->evlgpu || (_dm_evl && _dm_evl->p[ievl])){
		cp2cpu(&_dm_evl->p[ievl], cudata->dm_evl[ievl][0], stream);
	    }
	}
    }
}
void curecon_t::mvm(dcell **_dmerr, dcell *_gradin){
    cp2gpu(gradin, _gradin);
    MVM->solve(dmfit, gradin, cgstream);
    cp2cpu(_dmerr, dmfit_vec, cgstream);
    cgstream.sync();
}

void curecon_t::tomo_test(SIM_T *simu){
    gpu_set(cuglobal->recongpu);
    const PARMS_T *parms=simu->parms;
    stream_t stream;
    RECON_T *recon=simu->recon;
    //Debugging. 
    dcell *rhsc=NULL;
    dcell *lc=NULL;
    dcell *rtc=NULL;
    curcell rhsg;
    curcell lg;
    curcell rtg;
    muv(&rhsc, &recon->RR, simu->gradlastol, 1);
    writebin(rhsc, "CPU_TomoR");
    muv_trans(&rtc, &recon->RR, rhsc, 1);
    writebin(rtc, "CPU_TomoRt");
    if(parms->tomo.alg==1){
	muv(&lc, &recon->RL, rhsc, 1);
	writebin(lc, "CPU_TomoL");
	muv(&lc, &recon->RL, rhsc, -1);
	writebin(lc, "CPU_TomoL2");
	if(parms->tomo.precond==1){
	    dcell *lp=NULL;
	    fdpcg_precond(&lp, recon, rhsc);
	    writebin(lp, "CPU_TomoP");
	    fdpcg_precond(&lp, recon, rhsc);
	    writebin(lp, "CPU_TomoP2");
	    dcellfree(lp);
	}
    }
    dcellzero(lc);
    for(int i=0; i<5; i++){
	muv_solve(&lc, &recon->RL, NULL, rhsc);
	writebin(lc, "CPU_TomoCG%d", i);
    }
	
    cp2gpu(gradin, simu->gradlastol);
    RR->R(rhsg, 0, gradin, 1, stream);
    cuwrite(rhsg, "GPU_TomoR");
    RR->Rt(rtg, 0, rhsg, 1, stream);
    cuwrite(rtg, "GPU_TomoRt");
    if(parms->tomo.alg==1){
	cucg_t *RL2=dynamic_cast<cucg_t*>(RL);
	RL2->L(lg, 0, rhsg, 1,stream);
	cuwrite(lg, "GPU_TomoL"); 
	RL2->L(lg, 1, rhsg, -1,stream);
	cuwrite(lg, "GPU_TomoL2");
	if(parms->tomo.precond==1){
	    RL2=dynamic_cast<cucg_t*>(RL);
	    curcell lp;
	    RL2->Pre(lp, rhsg, stream);
	    cuwrite(lp, "GPU_TomoP");
	    RL2->Pre(lp, rhsg, stream);
	    cuwrite(lp, "GPU_TomoP2");
	}
    }
    cuzero(lg, stream);
    for(int i=0; i<5; i++){
	RL->solve(lg, rhsg, stream);
	cuwrite(lg, "GPU_TomoCG%d", i);
    }
    CUDA_SYNC_DEVICE;
    exit(0);
}
void curecon_t::fit_test(SIM_T *simu){	//Debugging. 
    gpu_set(cuglobal->recongpu);
    stream_t stream;
    const RECON_T *recon=simu->recon;
    dcell *rhsc=NULL;
    dcell *lc=NULL;	
    if(!simu->opdr && opdr_vec){
	cp2cpu(&simu->opdr, opdr_vec, 0);
    }
    if(!simu->parms->gpu.tomo && simu->opdr){
	cp2gpu(opdr_vec, simu->opdr);
    }
    writebin(simu->opdr, "opdr");
    muv(&rhsc, &recon->fit->FR, simu->opdr, 1);
    writebin(rhsc, "CPU_FitR");
    muv(&lc, &recon->fit->FL, rhsc, 1);
    writebin(lc, "CPU_FitL");
    muv(&lc, &recon->fit->FL, rhsc, -1);
    writebin(lc, "CPU_FitL2");
    dcellzero(lc);
    for(int i=0; i<5; i++){
	muv_solve(&lc, &recon->fit->FL, NULL, rhsc);
	writebin(lc, "CPU_FitSolve%d", i);
    }
    dcell *lhs=NULL;
    if(recon->fit->FR.M){
	muv_trans(&lhs, &recon->fit->FR, rhsc, 1);
	writebin(lhs, "CPU_FitRt");
    }
    curcell rhsg;
    curcell lg;
    FR->R(rhsg, 0.f, opdr, 1.f, stream);
    cuwrite(rhsg, "GPU_FitR");
    cucg_t *FL2=dynamic_cast<cucg_t*>(FL);
    if(FL2){
	FL2->L(lg, 0, rhsg, 1, stream);
	cuwrite(lg, "GPU_FitL");
	FL2->L(lg, 1, rhsg, -1, stream);
	cuwrite(lg, "GPU_FitL2");
    }
    cuzero(lg, stream);
    for(int i=0; i<5; i++){
	FL->solve(lg, rhsg, stream);
	cuwrite(lg, "GPU_FitSolve%d", i);
    }
    //Start from the same RHS.
    cp2gpu(rhsg, rhsc);
    cuzero(lg, stream);
    for(int i=0; i<5; i++){
	FL->solve(lg, rhsg, stream);
	cuwrite(lg, "GPU_FitSolveCPU%d", i);
    }
    curcell lhsg;
    FR->Rt(lhsg, 0, rhsg, 1, stream);
    cuwrite(lhsg, "GPU_FitRt");
    CUDA_SYNC_DEVICE;
    exit(0);
}
}//namespace

typedef cuda_recon::curecon_t curecon_t;
void gpu_setup_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    for(int igpu=0; igpu<NGPU; igpu++){
	if((parms->recon.mvm && parms->recon.alg==0 && parms->gpu.tomo && parms->gpu.fit && !parms->load.mvm)
	   || igpu==cuglobal->recongpu){
	    gpu_set(igpu);
	    if(cudata->recon){
		delete cudata->recon;
	    }
	    cudata->recon=new curecon_t(parms, powfs, recon);
	}
    }
}
void gpu_setup_recon_mvm(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs){
    //The following routine assemble MVM and put in recon->MVM
    if(!recon->MVM){
	if(parms->recon.mvm==1){
	    gpu_setup_recon_mvm_trans(parms, recon);
	}else{
	    gpu_setup_recon_mvm_direct(parms, recon);
	}
    }
    //free existing data
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	if(cudata->recon){
	    delete cudata->recon;
	    cudata->recon=NULL;
	}
    }

    if(!parms->sim.mvmport){
	gpu_set(cuglobal->recongpu);
	//recreate curecon_t that uses MVM.
	cudata->recon=new curecon_t(parms, powfs, recon);
    }
    gpu_print_mem("MVM");
}
void gpu_update_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon){
    if(!parms->gpu.tomo) return;
    for(int igpu=0; igpu<NGPU; igpu++){
	if((parms->recon.mvm && parms->gpu.tomo && parms->gpu.fit && !parms->load.mvm)
	   || igpu==cuglobal->recongpu){
	    gpu_set(igpu);
	    
	    curecon_t *curecon=cudata->recon;
	    curecon->update(parms, powfs, recon);
	}
    }
}
void gpu_update_recon_cn2(const PARMS_T *parms, RECON_T *recon){
    if(!parms->gpu.tomo) return;
    for(int igpu=0; igpu<NGPU; igpu++){
	if((parms->recon.mvm && parms->gpu.tomo && parms->gpu.fit && !parms->load.mvm)
	   || igpu==cuglobal->recongpu){
	    gpu_set(igpu);
	    
	    curecon_t *curecon=cudata->recon;
	    curecon->update_cn2(parms, recon);
	}
    }
}
void gpu_recon_reset(const PARMS_T *parms){//reset warm restart.
    for(int igpu=0; igpu<NGPU; igpu++){
	gpu_set(igpu);
	curecon_t *curecon=cudata->recon;
	if(curecon){
	    curecon->reset_runtime();
	}
	if(cudata->dm_wfs){
	    for(int iwfs=0; iwfs<parms->nwfs; iwfs++){
		gpu_set(cuglobal->wfsgpu[iwfs]);
		if(cudata->dm_wfs[iwfs]){
		    cuzero(cudata->dm_wfs[iwfs][0].p);
		}
	    }
	}
	if(cudata->dm_evl){
	    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
		gpu_set(cuglobal->evlgpu[ievl]);
		if(cudata->dm_evl[ievl]){
		    cuzero(cudata->dm_evl[ievl][0].p);
		}
	    }
	}
	CUDA_SYNC_DEVICE;
    }
}
void gpu_tomo(SIM_T *simu, dcell *gradin){
    gpu_set(cuglobal->recongpu);
    curecon_t *curecon=cudata->recon;
    curecon->grid->reconisim=simu->reconisim;
    const PARMS_T *parms=simu->parms;
    RECON_T *recon=simu->recon;
    if(parms->dbg.tomo){
	curecon->tomo_test(simu);
    }else{
	int copy2cpu=(parms->plot.run>1 
		      || !parms->gpu.fit 
		      || parms->save.opdr 
		      || (recon->moao && !parms->gpu.moao)
		      || parms->evl.tomo);
	simu->cgres->p[0]->p[simu->reconisim]=
	    curecon->tomo(copy2cpu?&simu->opdr:NULL, &simu->gngsmvst, gradin);
    }
}

void gpu_fit(dcell **dmout, SIM_T *simu){
    gpu_set(cuglobal->recongpu);
    curecon_t *curecon=cudata->recon;
    curecon->grid->reconisim=simu->reconisim;
    const PARMS_T *parms=simu->parms;
    if(parms->dbg.fit){
	curecon->fit_test(simu);
    }else{
	simu->cgres->p[1]->p[simu->reconisim]=
	    curecon->fit(dmout, parms->gpu.tomo?NULL:simu->opdr);
    }
    //Don't free opdr. Needed for warm restart in tomo.
}
void gpu_recon_mvm(dcell **dmout, dcell *gradin){
    gpu_set(cuglobal->recongpu);
    curecon_t *curecon=cudata->recon;
    curecon->mvm(dmout, gradin);
}

void gpu_moao_recon(SIM_T *simu){
    gpu_set(cuglobal->recongpu);
    const PARMS_T *parms=simu->parms;
    curecon_t *curecon=cudata->recon;
    curecon->moao_recon(parms->gpu.fit?NULL:simu->dmfit, (parms->gpu.tomo || parms->gpu.fit)?NULL:simu->opdr);
}

void gpu_moao_filter(SIM_T *simu){
    gpu_set(cuglobal->recongpu);
    curecon_t *curecon=cudata->recon;
    curecon->moao_filter(simu->dm_wfs, simu->dm_evl);
}

