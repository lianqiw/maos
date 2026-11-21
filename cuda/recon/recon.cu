/*
  Copyright 2009-2026 Lianqi Wang
  
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
#include "../sim/cudata.h"
#include "tomo.h"
#include "fit.h"
#include "recon.h"
//#include "perf.h"
#if !USE_CPP
extern "C"{
#endif
#include "../../maos/fdpcg.h"
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
/**
   Create a curcell for tomography output on the xloc layers, as 2-d grids if
   tomo.square is set, otherwise as vectors. If p is not NULL, the cell wraps
   the existing memory p.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data.
   @param[in] p Existing device memory to wrap, or NULL to allocate.
   @return The new curcell.
*/
curcell new_opdr(const parms_t *parms, const recon_t* recon, Real *p){
	if(parms->tomo.square){
		return curcell(recon->npsr, 1, P(recon->xnx), P(recon->xny), p);
	}else{
		return curcell(recon->npsr, 1, P(recon->xnloc), (long*)NULL, p);
	}
}
/**
   Create a curcell for DM commands: modal vectors if recon.modal is set and
   zonal is 0, 2-d grids for MVR with fit.square, otherwise actuator vectors.
   If p is not NULL, the cell wraps the existing memory p.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data.
   @param[in] zonal Force zonal DM commands even if recon.modal is set.
   @param[in] p Existing device memory to wrap, or NULL to allocate.
   @return The new curcell.
*/
curcell new_dmr(const parms_t *parms, const recon_t* recon, int zonal, Real *p){
	if(parms->recon.modal && !zonal){
		return curcell(parms->ndm, 1, P(recon->anmod), (long*)NULL, p);
	}else if(parms->recon.alg==RECON_MVR && parms->fit.square){
		return curcell(parms->ndm, 1, P(recon->anx), P(recon->any), p);
	}else{
		return curcell(parms->ndm, 1, P(recon->anloc), (long*)NULL, p);
	}
}
/**
   Set up the GPU reconstructor: geometry, MVM, DM fitting, tomography and
   MOAO, depending on parms->gpu settings. The caller must specify current GPU.

   @param[in] parms Simulation parameters. If NULL, an empty object is created.
   @param[in] recon Reconstruction data.
*/
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

	if(((parms->recon.alg==RECON_MVR&&parms->gpu.fit)||parms->recon.mvm||parms->gpu.moao)
		||(parms->recon.alg==RECON_LSR&&parms->gpu.lsr)
		){
		dmrecon=new_dmr(parms, recon);
		dmrecon_vec=dmrecon.Vector();
	}
	if(parms->recon.modal){
		cp2gpu(amod, recon->amod);
	}
	if(parms->recon.alg==RECON_MVR&&(parms->gpu.tomo||parms->gpu.fit)
		&&!parms->sim.idealtomo&&!parms->load.mvm&&!recon->MVM&&!cuglobal->mvm){
		opdr=new_opdr(parms, recon);
		opdr_vec=opdr.Vector();
		if(parms->tomo.square){
			opdr_map=cumapcell(recon->npsr, 1);
			for(int i=0; i<recon->npsr; i++){
				opdr_map[i]=cumap_t(grid->xmap[i], opdr[i]);
			}
			cudata->opdr=opdr_map;
		}
	}
	gpu_print_mem("recon init");
}

/**
   Set up the MVM solver in GPU from cuglobal->mvm or recon->MVM if
   recon.mvm is set. This is only for initialization. Updating mvm requires
   calling gpu_setup_recon_mvm and perhaps more changes.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data.
*/
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
/**
   Set up the DM fitting right hand side (FR) and left hand side (FL) solvers
   in GPU according to gpu.fit and fit.alg. Skipped if MVM is used and already
   available. There is usually no need to update DM fitting.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data.
*/
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

/**
   Set up the tomography right hand side (RR) and left hand side (RL) solvers
   in GPU according to tomo.alg (CBS, CG with optional FDPCG, or SVD). Skipped
   if MVM is used and already available. Also called to update the control
   parameters.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data.
*/
void curecon_t::init_tomo(const parms_t*parms, recon_t*recon){
	int skip_tomofit=parms->recon.mvm&&(cuglobal->mvm||recon->MVM);
	if(parms->gpu.tomo&&!skip_tomofit){
		reset_tomo();
		dbg("initialize tomography in GPU %d.\n", current_gpu());
		if(parms->tomo.square && !parms->tomo.assemble){
			RR=new cutomo_grid(parms, recon, grid);
		}else{
			RR=new cusolve_sparse(parms->tomo.maxit, parms->tomo.cgwarm, 
				&recon->RR, &recon->RL);
		}
		switch(parms->tomo.alg){
		case ALG_CBS:
			RL=new cusolve_cbs(recon->RL.C, recon->RL.Up, recon->RL.Vp);
			break;
		case ALG_CG:
			RL=dynamic_cast<cusolve_cg*>(RR);
			if(parms->tomo.precond==1){
				if(parms->tomo.square){
					(dynamic_cast<cusolve_cg*>(RR))->precond=new cufdpcg_t(recon->fdpcg, grid);
				}else{
					error("To implement: FDPCG with tomo.square=0");
				}
			}
			break;
		case ALG_SVD:
			RL=new cusolve_mvm(recon->RL.MI);
			break;
		default:
			error("Invalid");
		}
	}
}
/**
   Set up the MOAO solvers in GPU and allocate the MOAO DM commands for the WFS
   and evaluation directions that use each MOAO DM. The MOAO DM grids are also
   registered in the WFS and evaluation GPUs.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data.
*/
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
				dm_moao[imoao][count]=curcell(1, 1, anx, any);
				dm_wfs[iwfs]=dm_moao[imoao][count][0];
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
	dmrecon_zonal=new_dmr(parms, recon, 1);
	dmrecon_zonal_vec=dmrecon_zonal.Vector();
}
/**
   Update the tomography operators after the turbulence profile (cn2) is
   updated: ray tracing, laplacian regularization, and FDPCG preconditioner.
   For predictive tomography, layer wind velocities are also updated.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data with updated operators.
*/
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
	if(parms->tomo.precond==1&&_RL){
		_RL->update_fdpcg(recon->fdpcg);
	}
}
/**
   Zero runtime data (tomography and DM fitting results, MOAO DM commands) used
   for warm restart.
*/
void curecon_t::reset_runtime(){
	dbg("zero out runtime data.\n");
	opdr.Zero();
	dmrecon.Zero();
	dm_wfs.Zero();
	dm_evl.Zero();
}

/**
   Tomography from gradients _gradin in GPU. The result is copied to *_opdr if
   _opdr is not NULL. For MVST split tomography, GXL*opdr is accumulated to
   *_gngsmvst. Synchronizes cgstream.

   @param[in,out] _opdr If not NULL, receives a copy of the tomography result.
   Allocated if NULL.
   @param[in,out] _gngsmvst For MVST, GXL*opdr is accumulated into it.
   @param[in] _gradin Gradients.
   @return The CG residual.
*/
Real curecon_t::tomo(dcell** _opdr, dcell** _gngsmvst, const dcell* _gradin){
	cp2gpu(gradin, _gradin);

	RR->R(tomo_rhs, 0, gradin, 1, cgstream);
	Real cgres=RL->solve(opdr, tomo_rhs, cgstream);

	if(_opdr){
		cp2cpu(_opdr, opdr_vec, cgstream);
	}
	if(GXL){
		curcellmm(gngsmvst, 0, GXL, opdr_vec, "nn", 1, cgstream);
		add2cpu(_gngsmvst, 1, gngsmvst, 1, cgstream);
	}
	cgstream.sync();
	return cgres;
}
/**
   Copy tomography result from the reconstruction GPU to all other GPUs (for
   evl.tomo).
*/
void curecon_t::copy_opdr(){
	for(int igpu=0; igpu<NGPU; igpu++){
		if(igpu!=cuglobal->recongpu){
			gpu_set(igpu);
			Copy(cudata->opdr, opdr_map);
		}
	}
	gpu_set(cuglobal->recongpu);
}
/**
   DM fitting in GPU. Uses _opdr as input if not NULL, otherwise the GPU
   tomography result. The DM commands are copied to *_dmrecon. Synchronizes
   cgstream.

   @param[in,out] _dmrecon Receives the DM commands. Allocated if NULL.
   @param[in] _opdr Tomography result from CPU, or NULL to use the GPU result.
   @return The CG residual.
*/
Real curecon_t::fit(dcell** _dmrecon, dcell* _opdr){
	if(_opdr){
		cp2gpu(opdr_vec, _opdr);
	}

	FR->R(fit_rhs, 0, opdr, 1, cgstream);
	Real cgres=FL->solve(dmrecon, fit_rhs, cgstream);

	add2cpu(_dmrecon, 0, dmrecon_vec, 1, cgstream);
	cgstream.sync();
	return cgres;
}
/**
   MOAO reconstruction for all MOAO DMs. _dmrecon and _opdr, if not NULL,
   override the GPU DM fitting and tomography results. Modal DM commands are
   converted to zonal first.

   @param[in] _dmrecon DM fitting result from CPU, or NULL.
   @param[in] _opdr Tomography result from CPU, or NULL.
   @return Always 0.
*/
Real curecon_t::moao_recon(dcell* _dmrecon, dcell* _opdr){
	if(_dmrecon){
		cp2gpu(dmrecon_vec, _dmrecon);
	}
	if(_opdr){
		cp2gpu(opdr_vec, _opdr);
	}
	curcell dmin=dmrecon;
	if(amod){
		curcellmm(dmrecon_zonal_vec, 0, amod, dmrecon, "nn", 1, cgstream);
		/*for(int idm=0; idm<dmrecon.Nx(); idm++){
			curmv(dmrecon_zonal_vec(idm), 0., amod(idm,idm), dmrecon(idm), 'n', 1, cgstream);
		}*/
		dmin=dmrecon_zonal;
	}
	for(int imoao=0; imoao<nmoao; imoao++){
		if(!moao[imoao]) continue;
		moao[imoao]->moao_solve(dm_moao[imoao], opdr, dmin, cgstream);
	}
	return 0;
}

/**
   Apply a low pass filter with gain gdm to the MOAO DM commands for WFS and
   evaluation directions: dm=(1-g)*dm+g*dm_new, in the GPU that owns each
   direction. The results are copied to _dm_wfs and _dm_evl if needed by CPU.

   @param[in,out] _dm_wfs MOAO DM commands of WFS on CPU. Elements are updated
   if not NULL or if WFS is not simulated in GPU.
   @param[in,out] _dm_evl MOAO DM commands of evaluation directions on CPU.
   Elements are updated if not NULL or if evaluation is not done in GPU.
*/
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
				Copy(temp, dm_wfs[iwfs], stream);
			} else{
				temp=dm_wfs[iwfs];
			}
			Real g=moao_gwfs->p[iwfs];
			Add(cudata->dm_wfs[iwfs][0].p, (Real)(1.-g), temp, g, stream);
			if(_dm_wfs&&(!wfsgpu||_dm_wfs->p[iwfs])){//gpu.moao implies fit.square=1
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
				Copy(temp, dm_evl[ievl], stream);
			} else{
				temp=dm_evl[ievl];
			}
			Real g=moao_gevl->p[ievl];
			Add(cudata->dm_evl[ievl][0].p, 1-g, temp, g, stream);
			if(_dm_evl&&(!cuglobal->evlgpu||_dm_evl->p[ievl])){
				cp2cpu(&_dm_evl->p[ievl], cudata->dm_evl[ievl][0], stream);
			}
		}
	}
}
/**
   Matrix vector multiply reconstruction: *_dmerr=MVM*_gradin. Synchronizes
   cgstream.

   @param[in,out] _dmerr DM error signal. Allocated if NULL, otherwise overwritten.
   @param[in] _gradin Gradients.
*/
void curecon_t::mvm(dcell** _dmerr, dcell* _gradin){
	cp2gpu(gradin, _gradin);
	MVM->solve(dmrecon, gradin, cgstream);
	cp2cpu(_dmerr, dmrecon_vec, cgstream);
	cgstream.sync();
}

/**
   Debugging: compare tomography operators (R, Rt, L, preconditioner, solve)
   between CPU and GPU by saving their results to files. Exits the program.

   @param[in] simu Simulation data. simu->gradlastol is used as input.
*/
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
		writebin(lc, "CPU_TomoSolve%d", i);
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
		cuwrite(lg, stream, "GPU_TomoSolve%d", i);
		CUDA_SYNC_STREAM;//check for errors
	}
	CUDA_SYNC_DEVICE;
	exit(0);
}
/**
   Debugging: compare DM fitting operators (R, L, solve) between CPU and GPU by
   saving their results to files. Exits the program.

   @param[in,out] simu Simulation data. simu->opdr is used as input and
   filled from GPU if empty.
*/
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

/**
   Create the GPU reconstructor in the reconstruction GPU, replacing any
   existing one.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data.
*/
void gpu_setup_recon(const parms_t* parms, recon_t* recon){
	gpu_set(cuglobal->recongpu);
	if(cudata->recon){
		delete cudata->recon;
	}
	cudata->recon=new curecon_t(parms, recon);
}
/**
   Assembles MVM in gpu and put in recon->MVM if it does not exist yet. A
   temporary GPU reconstructor is created in each GPU and deleted afterwards.

   @param[in] parms Simulation parameters. tomo.cgwarm and fit.cgwarm are
   temporarily set to 0 and restored.
   @param[in,out] recon Reconstruction data. recon->MVM is set.
*/
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
/**
   Update the GPU tomography solver after control parameters change. Not
   implemented for MVM.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data.
*/
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
/**
   Update the GPU tomography operators after the turbulence profile is updated.
   Not implemented for MVM.

   @param[in] parms Simulation parameters.
   @param[in] recon Reconstruction data with updated turbulence profile.
*/
void gpu_update_recon_cn2(const parms_t* parms, recon_t* recon){
	if(!parms->gpu.tomo) return;
	for(int igpu=0; igpu<NGPU; igpu++){
		if((parms->recon.mvm&&parms->gpu.tomo&&parms->gpu.fit)||igpu==cuglobal->recongpu){
			gpu_set(igpu);
			curecon_t* curecon=cudata->recon;
			if(parms->recon.mvm){
				warning("Not yet implemented, skip.\n");
			} else{
				curecon->update_cn2(parms, recon);
			}
		}
	}
}
/**
   Reset GPU reconstruction runtime data and MOAO DM commands in all GPUs, for
   example between seeds. Synchronizes each device.

   @param[in] parms Simulation parameters.
*/
void gpu_recon_reset(const parms_t* parms){//reset warm restart.
	for(int igpu=0; igpu<NGPU; igpu++){
		gpu_set(igpu);
		curecon_t* curecon=cudata->recon;
		if(curecon){
			curecon->reset_runtime();
		}
		//MOAO DM commands of each direction are registered only in the GPU that owns it.
		if(cudata->dm_wfs){
			for(int iwfs=0; iwfs<parms->nwfs&&iwfs<cudata->dm_wfs.N(); iwfs++){
				if(cudata->dm_wfs[iwfs]){
					cuzero(cudata->dm_wfs[iwfs][0].p);
				}
			}
		}
		if(cudata->dm_evl){
			for(int ievl=0; ievl<parms->evl.nevl&&ievl<cudata->dm_evl.N(); ievl++){
				if(cudata->dm_evl[ievl]){
					cuzero(cudata->dm_evl[ievl][0].p);
				}
			}
		}
		CUDA_SYNC_DEVICE;
	}
}
/**
   Tomography in GPU from gradin. Stores the CG residual in simu->cgres and
   copies the result to simu->opdr when needed by CPU (plotting, CPU fitting or
   MOAO, saving, or evl.tomo). Runs tomo_test() instead if dbg.tomo is set.

   @param[in,out] simu Simulation data.
   @param[in] gradin Gradients.
*/
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
		if(parms->evl.tomo){
			curecon->copy_opdr();//copy to all gpus
		}
	}
}

/**
   DM fitting in GPU. Uses simu->opdr as input if tomography is done in CPU.
   The result is copied to *dmout and the CG residual stored in simu->cgres.
   Runs fit_test() instead if dbg.fit is set.

   @param[in,out] dmout DM commands. Allocated if NULL, otherwise overwritten.
   @param[in,out] simu Simulation data.
*/
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
/**
   MVM reconstruction in GPU: *dmout=MVM*gradin.

   @param[in,out] dmout DM error signal. Allocated if NULL, otherwise overwritten.
   @param[in] gradin Gradients.
*/
void gpu_recon_mvm(dcell** dmout, dcell* gradin){
	gpu_set(cuglobal->recongpu);
	curecon_t* curecon=cudata->recon;
	curecon->mvm(dmout, gradin);
}

/**
   MOAO reconstruction in GPU. Uses CPU tomography or fitting results from
   simu when those steps are not done in GPU.

   @param[in] simu Simulation data.
*/
void gpu_moao_recon(sim_t* simu){
	gpu_set(cuglobal->recongpu);
	const parms_t* parms=simu->parms;
	curecon_t* curecon=cudata->recon;
	curecon->moao_recon(parms->gpu.fit?NULL:simu->dmrecon, (parms->gpu.tomo||parms->gpu.fit)?NULL:simu->opdr);
}

/**
   Temporal filtering of MOAO DM commands in GPU, results copied to
   simu->dm_wfs and simu->dm_evl when needed.

   @param[in,out] simu Simulation data. simu->dm_wfs and simu->dm_evl are updated.
*/
void gpu_moao_filter(sim_t* simu){
	gpu_set(cuglobal->recongpu);
	curecon_t* curecon=cudata->recon;
	curecon->moao_filter(simu->dm_wfs, simu->dm_evl);
}

