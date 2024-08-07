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

#ifndef AOS_CUDA_ACCPHI_H
#define AOS_CUDA_ACCPHI_H
#include "../math/math.h"
#include "gpu_sim.h"
void map2loc(const cumap_t& map, const culoc_t& loc, Real* phiout,
	Real alpha, Real dispx, Real dispy, Real scale, int wrap, cudaStream_t stream);
void atm2loc(Real* phiout, const culoc_t& loc, Real hs, Real hc, Real thetax, Real thetay,
	Real misregx, Real misregy, Real dt, int isim, Real atmalpha, cudaStream_t stream);
void mapcell2loc(Real* phiout, const Array<culoc_t>& locarr, const cumapcell& cudm, 
	Real hs, Real hc, Real thetax, Real thetay,
	Real misregx, Real misregy, Real dmalpha, cudaStream_t stream);
void mapcell2loc(Real* phiout, const culoc_t& locout, const cumapcell& cudm, 
	Real hs, Real hc, Real thetax, Real thetay,
	Real misregx, Real misregy, Real dmalpha, cudaStream_t stream);

void ngsmod2loc(curmat& opd, Real(*restrict loc)[2],
	const ngsmod_t* ngsmod, const real* mod,
	real thetax, real thetay,
	real alpha, cudaStream_t stream);
#define KARG_COMMON const Real (*restrict loc)[2], const int nloc, const Real dxi, const Real dyi, const Real dispx, const Real dispy, const Real alpha
__global__ void map2loc_linear(Real* restrict out, const Real* restrict in, const int nx, const int ny, KARG_COMMON);
__global__ void map2loc_linear(Real* restrict out, const Comp* restrict in, const int nx, const int ny, KARG_COMMON);


#endif
