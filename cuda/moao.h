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
#ifndef AOS_CUDA_MOAO_H
#define AOS_CUDA_MOAO_H
#include "types.h"
#include "solve.h"
#include "recon_base.h"
namespace cuda_recon{
class cumoao_l:public cucg_t{
protected:
    curecon_geom *grid;
    curmat *NW, *dotNW;
    cugrid_t *amap;
    cusp *actslave;
    curcell *opdfit;
    curcell *opdfit2;
    map_ray *ha;
public:
    cumoao_l(const PARMS_T *parms, MOAO_T *moao, curecon_geom *_grid);
    virtual ~cumoao_l(){
	delete NW;
	delete dotNW;
	delete amap;
	delete actslave;
	delete opdfit;
	delete opdfit2;
	delete ha;
    }
    virtual void L(curcell **xout, Real beta, const curcell *xin, Real alpha, stream_t &stream);
};

class cumoao_t:public cumoao_l{
    map_ray **hxp;
    map_ray **hap;
    int ndir;
    curcell *rhs;
public:
    cumoao_t(const PARMS_T *parms, MOAO_T *moao, dir_t *dir, int _ndir, curecon_geom *_grid);
    Real moao_solve(curcell **xout, const curcell *xin, const curcell *ain, stream_t &stream);
    ~cumoao_t(){
	for(int idir=0; idir<ndir; idir++){
	    delete hxp[idir];
	    delete hap[idir];
	}
	delete hxp;
	delete hap;
	delete rhs;
    }
};
}//namespace
#endif
