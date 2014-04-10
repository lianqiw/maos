/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_CUDA_PROP_WRAP_H
#define AOS_CUDA_PROP_WRAP_H
#include "common.h"
/*data to be used by kernel */
struct PROP_WRAP_T{
    int offdir, offdirx, offdiry;
    int offps, offpsx, offpsy;
    Real *cc;
    int nxdir,nydir;
    int nxps,nyps;
    Real dispx, dispy;
    Real xratio, yratio;
    int nx, ny;
    int isreverse;/*Indicate this is an reverse*/
    PROP_WRAP_T *reverse;/*store pointer for the reversely prepared data.*/
    PROP_WRAP_T(){ /*This is not good. zeros out already initialized childs.*/
	memset(this, 0, sizeof(*this));
    }
    void togpu(PROP_WRAP_T *pgpu){
	if(reverse){
	    PROP_WRAP_T *gpureverse;
	    DO(cudaMalloc(&gpureverse, sizeof(PROP_WRAP_T)));
	    reverse->togpu(gpureverse);
	    delete reverse;
	    reverse=gpureverse;
	}
	DO(cudaMemcpy(pgpu, this, sizeof(PROP_WRAP_T),cudaMemcpyHostToDevice));  
    }
};
__global__ void 
gpu_prop_grid_do(PROP_WRAP_T *data, Real **pdirs, Real **ppss, 
		 int ndir, int nps, Real alpha1, Real *alpha2, char trans);

void gpu_prop_grid_prep(PROP_WRAP_T*res, 
			const cugrid_t &g_dir, const cugrid_t &gi,
			Real dispx, Real dispy, curmat *cc);
#endif
