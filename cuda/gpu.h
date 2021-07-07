/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_CUDA_GPU_H
#define AOS_CUDA_GPU_H
#if defined(__cplusplus) && !USE_CPP
extern "C"{
#endif
#include "../lib/aos.h"
#include "../maos/parms.h"
#include "../maos/types.h"

	void gpu_dbg(void);
	int  gpu_init(const parms_t* parms, int* gpus, int ngpu);
	void gpu_cleanup(void);
	void gpu_atm2gpu(const mapcell* atm, const dmat* atmscale, const parms_t* parms, int iseed, int isim);
	void gpu_dmreal2gpu(mapcell* dmreal);
	void gpu_dmproj2gpu(mapcell* dmproj);

	void gpu_wfsgrad_init(const parms_t* parms, const powfs_t* powfs);
	void gpu_wfs_init_sim(const parms_t* parms, powfs_t* powfs);
	void gpu_wfsgrad_update_etf(const parms_t* parms, const powfs_t* powfs, int ipowfs);
	void gpu_wfsgrad_update_mtche(const parms_t* parms, const powfs_t* powfs, int ipowfs);
	void gpu_wfsgrad_seeding(const parms_t* parms, const powfs_t* powfs, rand_t* rstat);
	void gpu_wfsgrad_queue(thread_t* info);
	void gpu_wfsgrad_sync(sim_t* simu, int iwfs);
	dmat* gpu_pywfs_mkg(const struct pywfs_t* pywfs, const loc_t* locin, const loc_t* locfft, const dmat* mod, real displacex, real displacey);
	void gpu_save_pistat(sim_t* simu);
	void gpu_wfssurf2gpu(const parms_t* parms, powfs_t* powfs);
	void gpu_perfevl_init(const parms_t* parms, aper_t* aper);
	void gpu_perfevl_ngsr(sim_t* simu, real* cleNGSm);
	void gpu_perfevl_queue(thread_t* info);
	void gpu_perfevl_sync(thread_t* info);
	void gpu_perfevl_save(sim_t* simu);
	void gpu_perfevl_init_sim(const parms_t* parms, aper_t* aper);
	void gpu_setup_recon(const parms_t* parms, recon_t* recon);
	void gpu_setup_recon_mvm(const parms_t* parms, recon_t* recon);
	void gpu_update_recon(const parms_t* parms, recon_t* recon);
	void gpu_update_recon_cn2(const parms_t* parms, recon_t* recon);
	void gpu_recon_reset(const parms_t* parms);
	void gpu_recon_free(void);
	void gpu_tomo(sim_t* simu, dcell* gradin);
	void gpu_fit(dcell** dmout, sim_t* simu);
	void gpu_recon_mvm(dcell** dmout, dcell* gradin);
	void gpu_moao_recon(sim_t* simu);
	void gpu_moao_filter(sim_t* simu);
	void gpu_moao_2gpu(sim_t* simu);
	void gpu_mvm_daemon(int port);
	void mvm_iwfs(int* gpus, int ngpu, int nstep);
	void mvm_only(int* gpus, int ngpu, int nstep);
	void mvmfull_iwfs(int* gpus, int ngpu, int nstep);
	void mvmfull_real(int* gpus, int ngpu, int nstep);
	void mvmfull_pipe(const char* mvm1, const char* mvm2, const char* pix1, const char* pix2, const char* mtch, int* gpus, int ngpu, int nstep);
	void mvm_test(int igpu);
#if defined(__cplusplus) && !USE_CPP
}
#endif
#endif
