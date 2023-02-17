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

/**
 * \file gpu.cu
 * 
 * Wraps cuda routines for CPU data type.
 * */
#include "cublas.h"
#include "utils.h"
#include "gpu_math.h"
#if CUDA_VERSION>10000
void gpu_dsvd(dmat **U_, dmat **S_, dmat **Vt_, dmat *A_){
  NumArray<real, Gpu> U, S, Vt, A;
	cp2gpu(A, A_);
	cusvd(U, S, Vt, A); 
	cp2cpu(U_, U);
	cp2cpu(S_, S);
  //dmat *V_=NULL, *Vt2=NULL;
	//cp2cpu(&V_, V);
  //Vt2=dtrans(V_); dfree(V_);
  cp2cpu(Vt_, Vt);
  //dfree(Vt2);
}
#endif
