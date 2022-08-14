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

#include <signal.h>
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
curecon_t::curecon_t(const parms_t* parms, recon_t* recon)
	:grid(0), FR(0), FL(0), RR(0), RL(0), MVM(0), nmoao(0), moao(0){
	if(!parms) return;
	dbg("initialize reconstructor in GPU\n");
	if(parms->recon.alg==RECON_MVR&&!(parms->recon.mvm&&parms->load.mvm)){
		grid=new curecon_geom(parms, recon);//does not change
	}
	if(parms->recon.split==2){//mvst
		cp2gpu(GXL, recon->GXL);
	}
	init_mvm(parms, recon);
	init_fit(parms, recon);
	init_tomo(parms, recon);
	init_moao(parms, recon);

	if((parms->recon.alg==RECON_MVR&&parms->gpu.fit||parms->recon.mvm||parms->gpu.moao)
		||(parms->recon.alg==RECON_LSR&&parms->gpu.lsr)
		){
		if(parms->recon.alg==RECON_MVR&&parms->fit.square){
			dmfit=curcell(parms->ndm, 1, P(recon->anx), P(recon->any));
		} else if(parms->recon.modal){
			dmfit=curcell(parms->ndm, 1, P(recon->anmod), (long*)NULL);
		} else{
			dmfit=curcell(parms->ndm, 1, P(recon->anloc), (long*)NULL);
		}
		dmfit_vec.init(parms->ndm, 1);
		for(int idm=0; idm<parms->ndm; idm++){
			dmfit_vec[idm]=dmfit[idm].Vector();
		}
	}

	if(parms->recon.alg==RECON_MVR&&(parms->gpu.tomo||parms->gpu.fit)
		&&!parms->sim.idealfit&&!parms->load.mvm&&!recon->MVM&&!cuglobal->mvm){
		if(parms->tomo.square){
			opdr=curcell(recon->npsr, 1, P(recon->xnx), P(recon->xny));
		} else{
			opdr=curcell(recon->npsr, 1, P(recon->xnloc), (long*)NULL);
		}
		opdr_vec.init(recon->npsr, 1);
		for(int ips=0; ips<recon->npsr; ips++){
			opdr_vec[ips]=opdr[ips].Vector();
		}
	}
	gpu_print_mem("recon init");
}

void curecon_t::reset_fit(){
	if(FL || FR){
		dbg("reset DM fitting in GPU %d.\n", current_gpu());
		if(FL!=dynamic_cast<cusolve_l*>(FR)) delete FL; 
		FL=0;
		delete FR; FR=0;
	}
	dmfit.Zero();
}
void curecon_t::reset_tomo(){
	if(RL || RR){
		dbg("reset tomography in GPU %d.\n", current_gpu());
		if(RL!=dynamic_cast<cusolve_l*>(RR)) delete RL; 
		RL=0;
		delete RR; RR=0;
	}
	//opdr.Zero();//no need here. first_run is automatically set in tomo pcg.
}
void curecon_t::reset_mvm(){
	if(MVM){
		dbg("reset MVM in GPU %d.\n", current_gpu());
		delete MVM; MVM=0;
	}
}
void curecon_t::reset_moao(){
	if(moao){
		dbg("reset MOAO in GPU %d.\n", current_gpu());
		for(int im=0; im<nmoao; im++){
			delete moao[im];
		}
		delete[] moao;
		moao=0; nmoao=0;
	}
}
//this is only for initialization. updating mvm requires calling gpu_setup_recon_mvm and perhaps more changes
void curecon_t::init_mvm(const parms_t* parms, recon_t* recon){
	if(parms->recon.mvm){
		reset_mvm();
		dbg("initialize MVM in GPU %d.\n", current_gpu());
		if(cuglobal->mvm){
			MVM=new cusolve_mvm(cuglobal->mvm);
		} else if(recon->MVM){//MVM already exists
			MVM=new cusolve_mvm(recon->MVM);
		}//else: waiting for future call
	}
}
//There is usually no need to update DM fitting
void curecon_t::init_fit(const parms_t* parms, recon_t* recon){
	int skip_tomofit=parms->recon.mvm&&(cuglobal->mvm||recon->MVM);
	if(parms->gpu.fit&&!skip_tomofit){
		reset_fit();
		dbg("Initialize DM fit in GPU %d with gpu.fit=%d and fit.alg=%d.\n", current_gpu(), parms->gpu.fit, parms->fit.alg);
		switch(parms->gpu.fit){
		case 1:
			FR=new cusolve_sparse(parms->fit.maxit, parms->fit.cgwarm,
				&recon->fit->FR, &recon->fit->FL);
			break;
		case 2:
			FR=new cufit_grid(parms, recon, grid);
			break;
		default:
			error("Invalid");
		}
		switch(parms->fit.alg){
		case ALG_CBS:
			FL=new cusolve_cbs(recon->fit->FL.C, recon->fit->FL.Up, recon->fit->FL.Vp);
			break;
		case ALG_CG:
			FL=dynamic_cast<cusolve_l*>(FR);
			break;
		case ALG_SVD:
			FL=new cusolve_mvm(recon->fit->FL.MI);
			break;
		default:
			error("Invalid");
		}
	}
}

void curecon_t::init_tomo(const parms_t*parms, recon_t*recon){
	int skip_tomofit=parms->recon.mvm&&(cuglobal->mvm||recon->MVM);
	if(parms->gpu.tomo&&!skip_tomofit){
		reset_tomo();
		dbg("initialize tomography in GPU %d.\n", current_gpu());
		RR=new cutomo_grid(parms, recon, grid);
		switch(parms->tomo.alg){
		case ALG_CBS:
			RL=new cusolve_cbs(recon->RL.C, recon->RL.Up, recon->RL.Vp);
			break;
		case ALG_CG:
			RL=dynamic_cast<cusolve_l*>(RR);
			break;
		case ALG_SVD:
			RL=new cusolve_mvm(recon->RL.MI);
			break;
		default:
			error("Invalid");
		}
	}
}
void curecon_t::init_moao(const parms_t*parms, recon_t*recon){
	if(!parms->nmoao) return;
	dbg("initialize MOAO in GPU %d.\n", current_gpu());
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
				dm_moao[imoao][count]=curcell(1, 1, anx, any);
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
				dm_moao[imoao][count]=curcell(1, 1, anx, any);
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
					cudata->dm_wfs.init(nwfs, 1);
				}
				cudata->dm_wfs[iwfs].init(1, 1);
				cudata->dm_wfs[iwfs][0]=(recon->moao[imoao].amap->p[0]);
			}
		}
		if(parms->evl.moao==imoao){
			for(int ievl=0; ievl<parms->evl.nevl; ievl++){
				gpu_set(cuglobal->evlgpu[ievl]);
				if(!cudata->dm_evl){
					cudata->dm_evl.init(parms->evl.nevl, 1);
				}
				cudata->dm_evl[ievl].init(1, 1);
				cudata->dm_evl[ievl][0]=(recon->moao[imoao].amap->p[0]);
			}
		}
		gpu_set(cuglobal->recongpu);
	}
}
void curecon_t::update_cn2(const parms_t* parms, recon_t* recon){
	dbg("update turbulence profile.\n");
	if(parms->tomo.predict){
		for(int ips=0; ips<recon->npsr; ips++){
			int ips0=parms->atmr.indps->p[ips];
			grid->xmap[ips].vx=cudata->atm[ips0].vx;
			grid->xmap[ips].vy=cudata->atm[ips0].vy;
		}
	}
	cutomo_grid* _RR=dynamic_cast<cutomo_grid*>(RR);
	cutomo_grid* _RL=dynamic_cast<cutomo_grid*>(RL);
	if(_RL){
		_RL->init_hx(parms, recon);
	}
	if(_RR&&_RR!=_RL){
		_RR->init_hx(parms, recon);
	}
	if(parms->tomo.precond==1){
		_RL->update_fdpcg(recon->fdpcg);
	}
}
void curecon_t::reset_runtime(){
	dbg("zero out runtime data.\n");
	opdr.Zero();
	dmfit.Zero();
	dm_wfs.Zero();
	dm_evl.Zero();
}

Real curecon_t::tomo(dcell** _opdr, dcell** _gngsmvst,
	const dcell* _gradin){
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
Real curecon_t::fit(dcell** _dmfit, dcell* _opdr){
	if(_opdr){
		cp2gpu(opdr_vec, _opdr);
	}

	FR->R(fit_rhs, 0, opdr, 1, cgstream);
	Real cgres=FL->solve(dmfit, fit_rhs, cgstream);

	add2cpu(_dmfit, 0, dmfit_vec, 1, cgstream);
	cgstream.sync();
	return cgres;
}
Real curecon_t::moao_recon(dcell* _dmfit, dcell* _opdr){
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

void curecon_t::moao_filter(dcell* _dm_wfs, dcell* _dm_evl){
	info_once("MOAO temporal filter implemented with LPF\n");
	const int* wfsgpu=cuglobal->wfsgpu();
	if(dm_wfs){
		int nwfs=dm_wfs.Nx();
		for(int iwfs=0; iwfs<nwfs; iwfs++){
			if(!dm_wfs[iwfs]) continue;
			curmat temp;
			if(wfsgpu) gpu_set(wfsgpu[iwfs]);
			stream_t stream;
			if(wfsgpu&&wfsgpu[iwfs]!=cuglobal->recongpu){
				cucp(temp, dm_wfs[iwfs], stream);
			} else{
				temp=dm_wfs[iwfs];
			}
			Real g=moao_gwfs->p[iwfs];
			curadd(cudata->dm_wfs[iwfs][0], 1-g, temp, g, stream);
			if(!wfsgpu||(_dm_wfs&&_dm_wfs->p[iwfs])){//gpu.moao implies fit.square=1
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
			if(cuglobal->evlgpu&&cuglobal->evlgpu[ievl]!=cuglobal->recongpu){
				cucp(temp, dm_evl[ievl], stream);
			} else{
				temp=dm_evl[ievl];
			}
			Real g=moao_gevl->p[ievl];
			curadd(cudata->dm_evl[ievl][0], 1-g, temp, g, stream);
			if(!cuglobal->evlgpu||(_dm_evl&&_dm_evl->p[ievl])){
				cp2cpu(&_dm_evl->p[ievl], cudata->dm_evl[ievl][0], stream);
			}
		}
	}
}
void curecon_t::mvm(dcell** _dmerr, dcell* _gradin){
	cp2gpu(gradin, _gradin);
	MVM->solve(dmfit, gradin, cgstream);
	cp2cpu(_dmerr, dmfit_vec, cgstream);
	cgstream.sync();
}

void curecon_t::tomo_test(sim_t* simu){
	gpu_set(cuglobal->recongpu);
	const parms_t* parms=simu->parms;
	stream_t stream;
	recon_t* recon=simu->recon;
	//Debugging. 
	dcell* rhsc=NULL;
	dcell* lc=NULL;
	dcell* rtc=NULL;
	curcell rhsg;
	curcell lg;
	curcell rtg;
	muv(&rhsc, &recon->RR, simu->gradlastol, 1);
	writebin(rhsc, "CPU_TomoR");
	muv_trans(&rtc, &recon->RR, rhsc, 1);
	writebin(rtc, "CPU_TomoRt");
	if(parms->tomo.alg==ALG_CG){
		muv(&lc, &recon->RL, rhsc, 1);
		writebin(lc, "CPU_TomoL");
		muv(&lc, &recon->RL, rhsc, -1);
		writebin(lc, "CPU_TomoL2");
		if(parms->tomo.precond==1){
			dcell* lp=NULL;
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
	cuwrite(rhsg, stream, "GPU_TomoR");
	CUDA_SYNC_STREAM;//check for errors
	RR->Rt(rtg, 0, rhsg, 1, stream);
	cuwrite(rtg, stream, "GPU_TomoRt");
	CUDA_SYNC_STREAM;//check for errors
	if(parms->tomo.alg==ALG_CG){
		cusolve_cg* RL2=dynamic_cast<cusolve_cg*>(RL);
		RL2->L(lg, 0, rhsg, 1, stream);
		cuwrite(lg, stream, "GPU_TomoL");
		CUDA_SYNC_STREAM;//check for errors
		RL2->L(lg, 1, rhsg, -1, stream);
		cuwrite(lg, stream, "GPU_TomoL2");
		if(parms->tomo.precond==1){
			RL2=dynamic_cast<cusolve_cg*>(RL);
			curcell lp;
			RL2->Pre(lp, rhsg, stream);
			cuwrite(lp, stream, "GPU_TomoP");
			CUDA_SYNC_STREAM;//check for errors
			RL2->Pre(lp, rhsg, stream);
			cuwrite(lp, stream, "GPU_TomoP2");
		}
	}
	cuzero(lg, stream);
	for(int i=0; i<5; i++){
		RL->solve(lg, rhsg, stream);
		cuwrite(lg, stream, "GPU_TomoCG%d", i);
		CUDA_SYNC_STREAM;//check for errors
	}
	CUDA_SYNC_DEVICE;
	exit(0);
}
void curecon_t::fit_test(sim_t* simu){	//Debugging. 
	gpu_set(cuglobal->recongpu);
	stream_t stream;
	const recon_t* recon=simu->recon;
	dcell* rhsc=NULL;
	dcell* lc=NULL;
	if(!simu->opdr&&opdr_vec){
		cp2cpu(&simu->opdr, opdr_vec, 0);
	}
	if(!simu->parms->gpu.tomo&&simu->opdr){
		cp2gpu(opdr_vec, simu->opdr);
	}
	writebin(simu->opdr, "opdr");
	muv(&rhsc, &recon->fit->FR, simu->opdr, 1);
	writebin(rhsc, "CPU_FitR");
	muv(&lc, &recon->fit->FL, rhsc, 1);
	writebin(lc, "CPU_FitL");
	/*muv(&lc, &recon->fit->FL, rhsc, -1);
	writebin(lc, "CPU_FitL2");*/
	dcellzero(lc);
	for(int i=0; i<5; i++){
		muv_solve(&lc, &recon->fit->FL, NULL, rhsc);
		writebin(lc, "CPU_FitSolve%d", i);
	}
	/*dcell* lhs=NULL;
	if(recon->fit->FR.M){
		muv_trans(&lhs, &recon->fit->FR, rhsc, 1);
		writebin(lhs, "CPU_FitRt");
	}*/
	curcell rhsg;
	curcell lg;
	FR->R(rhsg, 0.f, opdr, 1.f, stream);
	cuwrite(rhsg, stream, "GPU_FitR");
	cusolve_cg* FL2=dynamic_cast<cusolve_cg*>(FL);
	if(FL2){
		FL2->L(lg, 0, rhsg, 1, stream);
		cuwrite(lg, stream, "GPU_FitL");
		/*FL2->L(lg, 1, rhsg, -1, stream);
		cuwrite(lg, stream, "GPU_FitL2");*/
	}
	cuzero(lg, stream);
	for(int i=0; i<5; i++){
		FL->solve(lg, rhsg, stream);
		cuwrite(lg, stream, "GPU_FitSolve%d", i);
	}
	//Start from the same RHS.
	cp2gpu(rhsg, rhsc);
	cuzero(lg, stream);
	for(int i=0; i<5; i++){
		FL->solve(lg, rhsg, stream);
		cuwrite(lg, stream, "GPU_FitSolveCPU%d", i);
	}
	/*curcell lhsg;
	FR->Rt(lhsg, 0, rhsg, 1, stream);
	cuwrite(lhsg, stream, "GPU_FitRt");*/ //fix crash
	CUDA_SYNC_DEVICE;
	exit(1);
}
}//namespace

typedef cuda_recon::curecon_t curecon_t;
void gpu_setup_recon(const parms_t* parms, recon_t* recon){
	gpu_set(cuglobal->recongpu);
	if(cudata->recon){
		delete cudata->recon;
	}
	cudata->recon=new curecon_t(parms, recon);
}
/**Assembles MVM in gpu and put in recon->MVM. */
void gpu_setup_recon_mvm(const parms_t* parms, recon_t* recon){
	if(!recon->MVM){
		const int tomo_cgwarm=parms->tomo.cgwarm;
		const int fit_cgwarm=parms->fit.cgwarm;
		//disable warm_restart to disable solve.cu checking
		((parms_t *)parms)->tomo.cgwarm=0;
		((parms_t *)parms)->fit.cgwarm=0;
		for(int igpu=0; igpu<NGPU; igpu++){
			gpu_set(igpu);
			if(cudata->recon){
				delete cudata->recon;
			}
			cudata->recon=new curecon_t(parms, recon);
		}

		if(parms->recon.mvm==1){
			gpu_setup_recon_mvm_trans(parms, recon);
		} else{
			gpu_setup_recon_mvm_direct(parms, recon);
		}

		//free existing data
		for(int igpu=0; igpu<NGPU; igpu++){
			gpu_set(igpu);
			if(cudata->recon){
				delete cudata->recon;
				cudata->recon=NULL;
			}
		}
		((parms_t *)parms)->tomo.cgwarm=tomo_cgwarm;
		((parms_t *)parms)->fit.cgwarm=fit_cgwarm;
		gpu_print_mem("MVM");
	}
}
void gpu_update_recon_control(const parms_t* parms, recon_t* recon){
	gpu_set(cuglobal->recongpu);
	for(int igpu=0; igpu<NGPU; igpu++){
		if((parms->recon.mvm&&parms->gpu.tomo&&parms->gpu.fit)||igpu==cuglobal->recongpu){
			gpu_set(igpu);
			curecon_t* curecon=cudata->recon;
			if(parms->recon.mvm){
				error("Please implement\n");
			}else{
				//curecon->init_fit(parms, recon);//temp
				curecon->init_tomo(parms, recon);
			}
		}
	}
}
void gpu_update_recon_cn2(const parms_t* parms, recon_t* recon){
	if(!parms->gpu.tomo) return;
	for(int igpu=0; igpu<NGPU; igpu++){
		if((parms->recon.mvm&&parms->gpu.tomo&&parms->gpu.fit)||igpu==cuglobal->recongpu){
			gpu_set(igpu);
			curecon_t* curecon=cudata->recon;
			if(parms->recon.mvm){
				error("Please implement\n");
			} else{
				curecon->update_cn2(parms, recon);
			}
		}
	}
}
void gpu_recon_reset(const parms_t* parms){//reset warm restart.
	for(int igpu=0; igpu<NGPU; igpu++){
		gpu_set(igpu);
		curecon_t* curecon=cudata->recon;
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
void gpu_tomo(sim_t* simu, dcell* gradin){
	gpu_set(cuglobal->recongpu);
	curecon_t* curecon=cudata->recon;
	curecon->grid->reconisim=simu->reconisim;
	const parms_t* parms=simu->parms;
	recon_t* recon=simu->recon;
	if(parms->dbg.tomo){
		curecon->tomo_test(simu);
	} else{
		int copy2cpu=((parms->plot.run&&draw_current("opdr", NULL))
			||!parms->gpu.fit
			||parms->save.opdr
			||(recon->moao&&!parms->gpu.moao)
			||parms->evl.tomo);
		simu->cgres->p[0]->p[simu->reconisim]=
			curecon->tomo(copy2cpu?&simu->opdr:NULL, &simu->gngsmvst, gradin);
	}
}

void gpu_fit(dcell** dmout, sim_t* simu){
	gpu_set(cuglobal->recongpu);
	curecon_t* curecon=cudata->recon;
	curecon->grid->reconisim=simu->reconisim;
	const parms_t* parms=simu->parms;
	if(parms->dbg.fit){
		curecon->fit_test(simu);
	} else{
		simu->cgres->p[1]->p[simu->reconisim]=
			curecon->fit(dmout, parms->gpu.tomo?NULL:simu->opdr);
	}
	//Don't free opdr. Needed for warm restart in tomo.
}
void gpu_recon_mvm(dcell** dmout, dcell* gradin){
	gpu_set(cuglobal->recongpu);
	curecon_t* curecon=cudata->recon;
	curecon->mvm(dmout, gradin);
}

void gpu_moao_recon(sim_t* simu){
	gpu_set(cuglobal->recongpu);
	const parms_t* parms=simu->parms;
	curecon_t* curecon=cudata->recon;
	curecon->moao_recon(parms->gpu.fit?NULL:simu->dmfit, (parms->gpu.tomo||parms->gpu.fit)?NULL:simu->opdr);
}

void gpu_moao_filter(sim_t* simu){
	gpu_set(cuglobal->recongpu);
	curecon_t* curecon=cudata->recon;
	curecon->moao_filter(simu->dm_wfs, simu->dm_evl);
}

