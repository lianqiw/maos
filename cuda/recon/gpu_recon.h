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
#ifndef AOS_CUDA_GPU_RECON_H
#define AOS_CUDA_GPU_RECON_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
/**
 * \file gpu_math.h
 * Routines that can be used by CPU for simulation.
 * */
#if defined(__cplusplus) && !USE_CPP
extern "C"{
#endif
#include "../../lib/aos.h"
#include "../../maos/parms.h"
#include "../../maos/types.h"

	void gpu_setup_recon(const parms_t* parms, recon_t* recon);
	void gpu_setup_recon_mvm(const parms_t* parms, recon_t* recon);
	void gpu_update_recon_control(const parms_t* parms, recon_t* recon);
	void gpu_update_recon_cn2(const parms_t* parms, recon_t* recon);
	void gpu_recon_reset(const parms_t* parms);
	void gpu_recon_free(void);
	void gpu_tomo(sim_t* simu, dcell* gradin);
	void gpu_fit(dcell** dmout, sim_t* simu);
	void gpu_recon_mvm(dcell** dmout, dcell* gradin);
	void gpu_moao_recon(sim_t* simu);
	void gpu_moao_filter(sim_t* simu);


#if defined(__cplusplus) && !USE_CPP
}
#endif
#endif
