/*
  Copyright 2009-2015 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#ifndef AOS_CUDA_SOLVE_H
#define AOS_CUDA_SOLVE_H
#include "utils.h"
#include "curmat.h"
#include "pcg.h"
namespace cuda_recon{
class cusolve_r{/*Interface for RHS*/
public:
    virtual void R(curcell **xout, Real beta, const curcell *xin, Real alpha, stream_t &stream)=0;
    virtual void Rt(curcell **xout, Real beta, const curcell *xin, Real alpha, stream_t &stream)=0;
    virtual ~cusolve_r(){}
};

class cusolve_l{/*Interface for LHS*/
public:
    virtual Real solve(curcell **xout, const curcell *xin, stream_t &stream)=0;
    virtual ~cusolve_l(){}
};

class cucgpre_t{/*Interface for preconditioner*/
public:
    virtual void P(curcell **xout, const curcell *xin, stream_t &stream)=0;
    void operator()(curcell **xout, const curcell *xin, stream_t &stream){
	P(xout, xin, stream);
    }
    virtual ~cucgpre_t(){}
};

class cucg_t:public cusolve_l{/*Implementes LHS with cg algorithm*/
    int maxit, warm_restart;
    CGTMP_T cgtmp;
protected:
    cucgpre_t *precond;
public:
    cucg_t(int _maxit=0, int _warm_restart=0)
	:maxit(_maxit),warm_restart(_warm_restart),precond(NULL){}
    virtual ~cucg_t(){
	delete precond;
    }
    /*Left hand side forward operaton*/
    virtual void L(curcell **xout, Real beta, const curcell *xin, Real alpha, stream_t &stream)=0;
    void operator()(curcell **xout, Real beta, const curcell *xin, Real alpha, stream_t &stream){
	L(xout, beta, xin, alpha, stream);
    }
    virtual Real solve(curcell **xout, const curcell *xin, stream_t &stream);
    void P(curcell **xout, const curcell *xin, stream_t &stream){
	if(precond){
	    (*precond)(xout, xin, stream);
	}
    }
};

class cumuv_t{
    cusp *M;
    curmat *U;
    curmat *V;
    curmat *Vx;
    int nx, ny, *nxs, *nys;
 public:
    cumuv_t(const MUV_T *in=0);
    ~cumuv_t();
    void Forward(curcell **out, Real beta, const curcell *in, Real alpha, stream_t &stream);
    void operator()(curcell **out, Real beta, const curcell *in, Real alpha, stream_t &stream){
	Forward(out, beta, in, alpha, stream);
    }
    void Trans(curcell **out, Real beta, const curcell *in, Real alpha, stream_t &stream);
};

class cusolve_sparse:public cusolve_r,public cucg_t{
protected:
    cumuv_t *CR, *CL;
public:
    cusolve_sparse(int _maxit=0, int _warm_restart=0, MUV_T *_R=0, MUV_T *_L=0);
    virtual ~cusolve_sparse(){
	info2("cusolve_sparse::destructor\n");
	delete CR;
	delete CL;
    }
    virtual void R(curcell **out, Real beta, 
		   const curcell *xin, Real alpha, stream_t &stream){
	CR->Forward(out, beta, xin, alpha, stream);
    }
    virtual void L(curcell **out, Real beta, 
		   const curcell *xin, Real alpha, stream_t &stream){
	CL->Forward(out, beta, xin, alpha, stream);
    }
    virtual void Rt(curcell **out, Real beta, 
		    const curcell *xin, Real alpha, stream_t &stream){
	CR->Trans(out, beta, xin, alpha, stream);
    }
};

class cusolve_cbs:public cusolve_l{
protected:
    cusp *Cl;
    int  *Cp;
    curmat *Up;
    curmat *Vp;
    curmat *y;
    curmat *Vr;
public:
    cusolve_cbs(spchol *_C=0, dmat *_Up=0, dmat *_Vp=0);
    virtual ~cusolve_cbs(){
	delete Vr;
	delete Cl;
	delete Up;
	delete Vp;
	delete y;
	info2("cusolve_cbs::destructor\n");
    }
    void chol_solve(Real *out, const Real *in,  stream_t &stream);
    virtual Real solve(curcell **xout, const curcell *xin, stream_t &stream);
};

class cusolve_mvm:public cusolve_l{
    curmat *M;
public:
    cusolve_mvm(dmat *_M=0):M(0){
	cp2gpu(&M, _M);
    }
    ~cusolve_mvm(){
	delete M;
    }
    virtual Real solve(curcell **xout, const curcell *xin, stream_t &stream){
	if(!*xout) *xout=curcellnew(xin);
	curmv((*xout)->m->p, 0., M, xin->m->p, 'n', 1., stream);
	return 0;
    }
};
}//namespace
#endif
