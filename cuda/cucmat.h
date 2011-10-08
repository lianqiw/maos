/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_CUDA_CUCMAT_H
#define AOS_CUDA_CUCMAT_H

#include "utils.h"
#include "types.h"
cucmat *cucnew(int nx, int ny);
cucmat *cucnew(int nx, int ny, float *p, int own);
cucmat *cucnew(int nx, int ny, cudaStream_t stream);
cucmat *cucref(cucmat *A);
void cucfree(cucmat *A);
cuccell *cuccellnew(int nx, int ny);
cuccell *cuccellnew(int nx, int ny, int mx, int my);
void curcellfree(curcell *A);
void cucwrite(const cucmat *A, const char *format, ...);

inline void cuczero(cucmat *A, cudaStream_t stream){
    if(A && A->p){
	DO(cudaMemsetAsync(A->p, 0, A->nx*A->ny*sizeof(fcomplex), stream));
    }
}
#endif
