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
#ifndef AOS_CUDA_PCG_H
#define AOS_CUDA_PCG_H
typedef void (*G_CGFUN)(curcell**, float, const void*, const curcell*, float);
typedef void (*G_PREFUN)(curcell**, const void*, const curcell*, cudaStream_t stream);

int gpu_pcg(curcell **px, 
	    G_CGFUN Amul, const void *A, 
	    G_PREFUN Mmul, const void *M, 
	    const curcell *b, int warm, int maxiter,
	    cudaStream_t stream);
#endif
