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
#ifndef AOS_CUDA_ACCPHI_H
#define AOS_CUDA_ACCPHI_H
#include "curmat.h"

void gpu_atm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, const float hs, const float thetax,const float thetay,
		 const float mispx, const float mispy, const float dtisim, const float atmalpha, cudaStream_t stream);
void gpu_dm2loc(float *phiout, const float (*restrict loc)[2], const int nloc, cumap_t **cudm, int ndm,
		const float hs, const float thetax, const float thetay,
		const float mispx, const float mispy, const float dmalpha, cudaStream_t stream);

void prop_grid_match(curmat *out, float oxo, float oyo,
		     const curmat *in, float oxi, float oyi, float dxi,
		     float dispx, float dispy,
		     float alpha, cudaStream_t stream);
void gpu_prop_grid(curmat *out, float oxo, float oyo, float dxo,
		   curmat *in, float oxi, float oyi, float dxi,
		   float dispx, float dispy,
		   float alpha, char trans, cudaStream_t stream);
float* gpu_dmcubic_cc(float iac);

void gpu_prop_grid_cubic(curmat *out, float oxo, float oyo, float dxo,
			 curmat *in, float oxi, float oyi, float dxi,
			 float dispx, float dispy, float *cc,
			 float alpha, char trans, cudaStream_t stream);

void gpu_prop_grid_prep(GPU_PROP_GRID_T*res, curmat *out, float oxo, float oyo, float dxo,
			curmat *in, float oxi, float oyi, float dxi,
			float dispx, float dispy, char trans);

void gpu_ngsmod2science(curmat *opd, float (*restrict loc)[2],
			const NGSMOD_T *ngsmod, const double *mod, 
			double thetax, double thetay, 
			double alpha, cudaStream_t stream);
#endif

