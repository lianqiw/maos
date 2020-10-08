/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>

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
#ifndef AOS_CUDA_MOAO_H
#define AOS_CUDA_MOAO_H
#include "types.h"
#include "solve.h"
#include "recon_geom.h"
namespace cuda_recon{
class cumoao_t:private cusolve_cg{
private:
	curecon_geom* grid;
	curmat NW, dotNW;
	cugridcell amap;
	cusp actslave;
	curcell opdfit;
	curcell opdfit2;
	map2map ha;
	Array<map2map>hxp;
	Array<map2map>hap;
	int ndir;
	curcell rhs;
public:
	friend class curecon_t;
	virtual void L(curcell& xout, Real beta, const curcell& xin, Real alpha, stream_t& stream);
	Real moao_solve(curccell& xout, const curcell& xin, const curcell& ain, stream_t& stream);
	cumoao_t(const PARMS_T* parms, MOAO_T* moao, dir_t* dir, int _ndir, curecon_geom* _grid);
	~cumoao_t(){}
	operator bool(){
		return amap?1:0;
	}
};
}//namespace
#endif
