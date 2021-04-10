/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "curmat.h"
namespace cuda_recon{
/**
   hold data struct for temporary data used for CG to avoid alloc/free at every call to CG.
*/
class cgtmp_t{
public:
	curcell r0;
	curcell z0;
	curcell p0;
	curcell Ap;
	curmat store;
	Array<Real, Pinned> diff;
	int count_fail, count;

	cgtmp_t():count_fail(0), count(0){}
	~cgtmp_t(){}
};
//typedef void (*G_cgfun_t)(curcell**, Real, const curcell*, Real, stream_t &stream);
//typedef void (*G_prefun_t)(curcell**, const curcell*, stream_t &stream);
class cusolve_cg;
class cusolve_cgpre;
Real pcg(curcell& x0, cusolve_cg* Amul, cusolve_cgpre* Mmul,
	const curcell& b, cgtmp_t& cg_data, int warm, int maxiter,
	stream_t& stream, Real cgthres=-1);
}//namespace
#endif
