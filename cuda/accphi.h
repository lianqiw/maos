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

void gpu_atm2loc(Real *phiout, culoc_t *loc, Real hs, Real thetax,Real thetay,
		 Real mispx, Real mispy, Real dt, int isim, Real atmalpha, cudaStream_t stream);
void gpu_dm2loc(Real *phiout, culoc_t **locarr, cumap_t *cudm, int ndm,
		Real hs, Real thetax, Real thetay,
		Real mispx, Real mispy, Real dmalpha, cudaStream_t stream);

void gpu_ngsmod2science(curmat *opd, Real (*restrict loc)[2],
			const NGSMOD_T *ngsmod, const double *mod, 
			double thetax, double thetay, 
			double alpha, cudaStream_t stream);
#endif
