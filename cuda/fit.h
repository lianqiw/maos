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
#include "recon_geom.h"
#include "prop_wrap.h"
/*cufit_grid implements both RHS and LHS. LHS is through cucg*/
class curecon_geom;
class cufit_grid:public cusolve_r,public cucg_t{/*grid based*/
    /*the following data are input or intermediate data for the operation*/
    curecon_geom *grid;
    int nfit;
    curcell *dmcache, *xcache;
    curcell *opdfit; /**<OPDs defined on ploc for fitting.*/
    curcell *opdfit2;/**<OPDs defined on ploc for fitting.*/
    curmat *opdfitv;/**<Concatenated version of opdfit. 1 column for each direction*/
    curmat *opdfit2v;/**<Concatenated version of opdfit2. 1 column for each direction*/
    curmat *pis;     /**<contains result of W1'*opdfit*/
    curmat *fitwt, *fitNW, *dotNW;
    culoc_t *floc;
    float *fit_thetax, *fit_thetay;
    cusp *actslave;
    PROP_WRAP_T *hxpdata,*hadata,*ha0data,*ha1data,*hxp0data,*hxp1data;
public:
    void init(const PARMS_T *parms, const RECON_T *recon);
    cufit_grid(const PARMS_T *parms=0, const RECON_T *recon=0, curecon_geom *_grid=0):
	cucg_t(parms?parms->fit.maxit:0, parms?parms->recon.warm_restart:0),grid(_grid),
	nfit(0),dmcache(0),xcache(0),opdfit(0),opdfit2(0), opdfitv(0),opdfit2v(0),
	pis(0),fitwt(0),fitNW(0),dotNW(0),floc(0),fit_thetax(0),fit_thetay(0),
	actslave(0), 
	hxpdata(0),hadata(0),ha0data(0),ha1data(0),hxp0data(0),hxp1data(0){
	if(!parms) return;
	init(parms, recon);
    }
    void do_hxp(const curcell *xin, stream_t &stream);
    void do_hxpt(const curcell *xout, float alpha, stream_t &stream);
    void do_w(stream_t &stream);
    void do_ha(const curcell *xin, stream_t &stream);
    void do_hat(curcell *xout,  float alpha, stream_t &stream);
    virtual ~cufit_grid(){
	info2("cufit_grid::destructor\n");
	if(!this) return;
	delete dmcache;
	delete xcache;
	delete opdfit;
	delete opdfit2;
	delete opdfitv;
	delete opdfit2v;
	delete pis;
	delete fitwt;
	delete fitNW;
	delete dotNW;
	delete [] fit_thetax;
	delete [] fit_thetay;
	delete floc;
	delete actslave;
	delete [] hadata;
	delete [] ha0data;
	delete [] ha1data;

	delete [] hxpdata;
	delete [] hxp0data;
	delete [] hxp1data;
    }
    virtual void R(curcell **out, float beta, 
		   const curcell *xin, float alpha, stream_t &stream);
    virtual void L(curcell **out, float beta, 
		   const curcell *xin, float alpha, stream_t &stream);
    virtual void Rt(curcell **out, float beta, 
		    const curcell *xin, float alpha, stream_t &stream);
};

class cufit_sparse:public cusolve_sparse{
public:
    cufit_sparse(const PARMS_T *parms, const RECON_T *recon);
};

#endif
