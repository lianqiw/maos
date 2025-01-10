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

#ifndef AOS_CUDA_PERF_H
#define AOS_CUDA_PERF_H
#include "../math/cumath.h"
#include "gpu_sim.h"

/**
   Data per GPU
 */
class cuperf_t{
public:
	culoc_t locs;
	Array<Array<culoc_t> > locs_dm;
	curmat imcc;
	curmat amp;
	Array<cuimat> embed;

	curcell psfol;
	curmat  opdcovol;
	curmat  opdmeanol;
	cuperf_t(){}
	~cuperf_t(){};
};
/**
   Global data independent of GPU.
*/
class cuperf_g{
public:
	int nevl;
	Array<int> nembed;
	Array<int> psfsize;
	Array<Real> wvls;
	Array<cufftHandle>    plan;
	curcell surf;
	curcell opd;
	curcell psfcl;
	curcell psfcl_ngsr;
	curcell opdcov;
	curcell opdcov_ngsr;
	curcell opdmean;
	curcell opdmean_ngsr;
	curcell cc_ol, cc_cl, coeff;
	Array<Array<Real, Pinned> >ccb_ol;
	Array<Array<Real, Pinned> >ccb_cl;
};
#endif
