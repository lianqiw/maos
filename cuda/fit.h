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

#ifndef AOS_CUDA_FIT_H
#define AOS_CUDA_FIT_H
#include "solve.h"
#include "recon_geom.h"
/*cufit_grid implements both RHS and LHS. LHS is through cucg*/
namespace cuda_recon{
class cufit_grid:public cusolve_r,public cusolve_cg{/*grid based*/
protected:
    /*the following data are input or intermediate data for the operation*/
    const curecon_geom *grid;
    int nfit;
    int idealfit;   /**<Use input from real turbulence instead of tomography output*/
    cugridcell acmap;
    curcell dmcache;
    curcell xcache;
    curcell opdfit; /**<OPDs defined on ploc for fitting.*/
    curcell opdfit2;/**<OPDs defined on ploc for fitting.*/
    curmat opdfitv;/**<Concatenated version of opdfit. 1 column for each direction*/
    curmat opdfit2v;/**<Concatenated version of opdfit2. 1 column for each direction*/
    curmat fitwt, fitNW, dotNW;
    culoc_t floc;
    dir_t *dir;
    cusp actslave;
    map2map hxp, hxp0, hxp1, ha, ha0, ha1;
    /*PROP_WRAP_T *hxpdata,*hxp0data,*hxp1data;
    PROP_WRAP_T *hapdata;//for moao
    PROP_WRAP_T *hadata,*ha0data,*ha1data;*/
    void do_hxp(const curcell &xin, stream_t &stream);
    void do_hxpt(const curcell &xout, Real alpha, stream_t &stream);
    void do_ha(const curcell &xin, stream_t &stream);
    void do_hat(curcell &xout,  Real alpha, stream_t &stream);
public:
    cufit_grid(const parms_t *parms=0, const recon_t *recon=0, const curecon_geom *_grid=0);
    virtual ~cufit_grid(){
	delete [] dir;
    }
    virtual void R(curcell &out, Real beta, 
		   curcell &xin, Real alpha, stream_t &stream);
    virtual void L(curcell &out, Real beta, 
		   const curcell &xin, Real alpha, stream_t &stream);
    virtual void Rt(curcell &out, Real beta, 
		    const curcell &xin, Real alpha, stream_t &stream);
};

}//namespace
#endif
