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
#ifndef AOS_CUDA_GPU_H
#define AOS_CUDA_GPU_H
#ifdef __cplusplus
extern "C"{
#endif
#include "../lib/aos.h"
#include "../maos/parms.h"
#include "../maos/types.h"
void gpu_info(void);
int  gpu_init(const PARMS_T *parms, int *gpus, int ngpu);
void gpu_cleanup(void);
void gpu_atm2gpu(const mapcell *atm, const dmat *atmscale, const PARMS_T *parms, int iseed, int isim);
void gpu_dmreal2gpu(mapcell *dmreal);
void gpu_dmproj2gpu(mapcell *dmproj);

void gpu_wfsgrad_init(const PARMS_T *parms, const POWFS_T *powfs);
void gpu_wfs_init_sim(const PARMS_T *parms, POWFS_T *powfs);
void gpu_wfsgrad_update_etf(const PARMS_T *parms, const POWFS_T *powfs);
void gpu_wfsgrad_update_mtche(const PARMS_T *parms, const POWFS_T *powfs);
void gpu_wfsgrad_seeding(const PARMS_T *parms, const POWFS_T *powfs, rand_t *rstat);
void gpu_wfsgrad_queue(thread_t *info);
void gpu_wfsgrad_sync(SIM_T *simu, int iwfs);
dmat *gpu_pywfs_mkg(const struct PYWFS_T *pywfs, const loc_t* locin, const dmat *mod, double displacex, double displacey, double scale);
void gpu_save_gradstat(SIM_T *simu);
void gpu_wfssurf2gpu(const PARMS_T *parms, POWFS_T *powfs);
void gpu_perfevl_init(const PARMS_T *parms, APER_T *aper);
void gpu_perfevl_ngsr(SIM_T *simu, double *cleNGSm);
void gpu_perfevl_queue(thread_t *info);
void gpu_perfevl_sync(thread_t *info);
void gpu_perfevl_save(SIM_T *simu);
void gpu_perfevl_init_sim(const PARMS_T *parms, APER_T *aper);
void gpu_setup_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon);
void gpu_setup_recon_mvm(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs);
void gpu_update_recon(const PARMS_T *parms, POWFS_T *powfs, RECON_T *recon);
void gpu_update_recon_cn2(const PARMS_T *parms, RECON_T *recon);
void gpu_recon_reset(const PARMS_T *parms);
void gpu_recon_free(void);
void gpu_tomo(SIM_T *simu);
void gpu_fit(SIM_T *simu);
void gpu_recon_mvm(dcell **dmout, dcell *gradin);
void gpu_moao_recon(SIM_T *simu);
void gpu_moao_filter(SIM_T *simu);
void gpu_moao_2gpu(SIM_T *simu);
void gpu_mvm_daemon(int port);
void mvm_iwfs(int *gpus, int ngpu, int nstep);
void mvm_only(int *gpus, int ngpu, int nstep);
void mvmfull_iwfs(int *gpus, int ngpu, int nstep);
void mvmfull_real(int *gpus, int ngpu, int nstep);
void mvmfull_pipe(char *mvm1, char *mvm2, char *pix1, char *pix2, char *mtch, int *gpus, int ngpu, int nstep);
void mvm_test(int igpu);
#ifdef __cplusplus
}
#endif
#endif

