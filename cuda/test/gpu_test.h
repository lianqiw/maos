/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#ifndef AOS_CUDA_GPU_TEST_H
#define AOS_CUDA_GPU_TEST_H
/**
 * \file gpu_test.h
 * Routines that can be used by CPU.
 * */
#if defined(__cplusplus) && !USE_CPP
extern "C"{
#endif
#include "../../lib/aos.h"
	void gpu_mvm_daemon(int port);
	void mvm_iwfs(int *gpus, int ngpu, int nstep);
	void mvm_only(int *gpus, int ngpu, int nstep);
	void mvmfull_iwfs(int *gpus, int ngpu, int nstep);
	void mvmfull_real(int *gpus, int ngpu, int nstep);
	void mvmfull_pipe(const char *mvm1, const char *mvm2, const char *pix1, const char *pix2, const char *mtch, int *gpus, int ngpu, int nstep);
	void mvm_test(int igpu);


#if defined(__cplusplus) && !USE_CPP
}
#endif
#endif
