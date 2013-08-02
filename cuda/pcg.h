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
#ifndef AOS_CUDA_PCG_H
#define AOS_CUDA_PCG_H
namespace cuda_recon{
/**
   hold data struct for temporary data used for CG to avoid alloc/free at every call to CG.
*/
typedef struct CGTMP_T{
    curcell *r0;
    curcell *z0;
    curcell *p0;
    curcell *Ap;
    float *store;
    float *diff;
    int count_fail;
    CGTMP_T(){
	memset(this, 0, sizeof(CGTMP_T));
    }
    ~CGTMP_T(){
	if(!this) return;
	delete r0;r0=NULL;
	delete z0;z0=NULL;
	delete p0;p0=NULL;
	delete Ap;Ap=NULL;
	if(store) cudaFree(store); store=NULL;
	if(diff) cudaFreeHost(diff); diff=NULL;
    }
}CGTMP_T;
//typedef void (*G_CGFUN)(curcell**, float, const curcell*, float, stream_t &stream);
//typedef void (*G_PREFUN)(curcell**, const curcell*, stream_t &stream);
class cucg_t;
class cucgpre_t;
float gpu_pcg(curcell **x0, cucg_t *Amul, cucgpre_t *Mmul,
	      const curcell *b, CGTMP_T *cg_data, int warm, int maxiter, 
	      stream_t &stream, double cgthres=-1);
}//namespace
#endif
