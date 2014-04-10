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
#ifndef AOS_CUDA_FIT_H
#define AOS_CUDA_FIT_H
#include "gpu.h"
#include "solve.h"
#include "recon_base.h"
#include "prop_wrap.h"
/*cufit_grid implements both RHS and LHS. LHS is through cucg*/
namespace cuda_recon{
class cufit_grid:public cusolve_r,public cucg_t{/*grid based*/
protected:
    /*the following data are input or intermediate data for the operation*/
    curecon_geom *grid;
    int nfit;
    cugrid_t *acmap;
    curcell *dmcache, *xcache;
    curcell *opdfit; /**<OPDs defined on ploc for fitting.*/
    curcell *opdfit2;/**<OPDs defined on ploc for fitting.*/
    curmat *opdfitv;/**<Concatenated version of opdfit. 1 column for each direction*/
    curmat *opdfit2v;/**<Concatenated version of opdfit2. 1 column for each direction*/
    curmat *fitwt, *fitNW, *dotNW;
    culoc_t *floc;
    dir_t *dir;
    cusp *actslave;
    map_ray *hxp, *hxp0, *hxp1, *ha, *ha0, *ha1;
    /*PROP_WRAP_T *hxpdata,*hxp0data,*hxp1data;
    PROP_WRAP_T *hapdata;//for moao
    PROP_WRAP_T *hadata,*ha0data,*ha1data;*/
public:
    cufit_grid(const PARMS_T *parms=0, const RECON_T *recon=0, curecon_geom *_grid=0);
    void do_hxp(const curcell *xin, stream_t &stream);
    void do_hxpt(const curcell *xout, Real alpha, stream_t &stream);
    void do_ha(const curcell *xin, stream_t &stream);
    void do_hat(curcell *xout,  Real alpha, stream_t &stream);
    virtual ~cufit_grid(){
	info2("cufit_grid::destructor\n");
	if(!this) return;
	delete[] acmap;
	delete dmcache;
	delete xcache;
	delete opdfit;
	delete opdfit2;
	delete opdfitv;
	delete opdfit2v;
	delete fitwt;
	delete fitNW;
	delete dotNW;
	delete dir;
	delete floc;
	delete actslave;
	delete [] ha;
	delete [] ha0;
	delete [] ha1;

	delete [] hxp;
	delete [] hxp0;
	delete [] hxp1;
    }
    virtual void R(curcell **out, Real beta, 
		   const curcell *xin, Real alpha, stream_t &stream);
    virtual void L(curcell **out, Real beta, 
		   const curcell *xin, Real alpha, stream_t &stream);
    virtual void Rt(curcell **out, Real beta, 
		    const curcell *xin, Real alpha, stream_t &stream);
};

}//namespace
#endif
