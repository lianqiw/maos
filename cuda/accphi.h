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
#ifndef AOS_CUDA_ACCPHI_H
#define AOS_CUDA_ACCPHI_H
#include "curmat.h"
void gpu_map2loc(const cumap_t &map, const culoc_t &loc, Real *phiout,
		 Real alpha, Real dispx, Real dispy, Real scale, int wrap, cudaStream_t stream);
void gpu_atm2loc(Real *phiout, const culoc_t &loc, Real hs, Real thetax,Real thetay,
		 Real mispx, Real mispy, Real dt, int isim, Real atmalpha, cudaStream_t stream);
void gpu_dm2loc(Real *phiout, const cuarray<culoc_t> &locarr, const cumapcell &cudm, int ndm,
		Real hs, Real thetax, Real thetay,
		Real mispx, Real mispy, Real dmalpha, cudaStream_t stream);
void gpu_dm2loc(Real *phiout, const culoc_t &locout, const cumapcell &cudm, int ndm,
		Real hs, Real thetax, Real thetay,
		Real mispx, Real mispy, Real dmalpha, cudaStream_t stream);

void gpu_ngsmod2science(curmat &opd, Real (*restrict loc)[2],
			const NGSMOD_T *ngsmod, const double *mod, 
			double thetax, double thetay, 
			double alpha, cudaStream_t stream);
#define KARG_COMMON const Real (*restrict loc)[2], const int nloc, const Real dxi, const Real dyi, const Real dispx, const Real dispy, const Real alpha
__global__ void prop_linear(Real *restrict out, const Real *restrict in, const int nx, const int ny, KARG_COMMON);
__global__ void prop_linear(Real *restrict out, const Comp *restrict in, const int nx, const int ny, KARG_COMMON);

#endif
