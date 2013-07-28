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
#ifndef AOS_CUDA_SOLVE_H
#define AOS_CUDA_SOLVE_H
#include "utils.h"
#include "curmat.h"
#include "pcg.h"

class cusolve_r{/*Interface for RHS*/
public:
    virtual void R(curcell **xout, float beta, const curcell *xin, float alpha, stream_t &stream)=0;
    virtual void Rt(curcell **xout, float beta, const curcell *xin, float alpha, stream_t &stream)=0;
};

class cusolve_l{/*Interface for LHS*/
public:
    virtual float solve(curcell **xout, const curcell *xin, cusolve_r *Rfun, stream_t &stream)=0;
};

class cucgpre_t{/*Interface for preconditioner*/
public:
    virtual void P(curcell **xout, const curcell *xin, stream_t &stream)=0;
    void operator()(curcell **xout, const curcell *xin, stream_t &stream){
	P(xout, xin, stream);
    }
};

class cucg_t:public cusolve_l{/*Implementes LHS with cg algorithm*/
    int maxit, warm_restart;
    CGTMP_T cgtmp;
    curcell *rhs;
protected:
    cucgpre_t *precond;
public:
    cucg_t(int _maxit=0, int _warm_restart=0)
	:maxit(_maxit),warm_restart(_warm_restart),rhs(NULL),precond(NULL){}
    /*Left hand side forward operaton*/
    virtual void L(curcell **xout, float beta, const curcell *xin, float alpha, stream_t &stream)=0;
    void operator()(curcell **xout, float beta, const curcell *xin, float alpha, stream_t &stream){
	L(xout, beta, xin, alpha, stream);
    }
    virtual float solve(curcell **xout, const curcell *xin, cusolve_r *Rfun, stream_t &stream);
};

class cumuv_t{
    cusp *M;
    curmat *U;
    curmat *V;
    curmat *Vx;
    int nx, ny, *nxs, *nys;
    
public:
    void init(const MUV_T *in);
    cumuv_t(MUV_T *in=0):M(0),U(0),V(0),Vx(0),nx(0),ny(0),nxs(0),nys(0){
	if(in) init(in);
    }
    void Forward(curcell **out, float beta, const curcell *in, float alpha, stream_t &stream);
    void operator()(curcell **out, float beta, const curcell *in, float alpha, stream_t &stream){
	Forward(out, beta, in, alpha, stream);
    }
    void Trans(curcell **out, float beta, const curcell *in, float alpha, stream_t &stream);
};

class cusolve_sparse:public cusolve_r,public cucg_t{
protected:
    cumuv_t CR, CL;
public:
    cusolve_sparse(int _maxit=0, int _warm_restart=0, MUV_T *_R=0, MUV_T *_L=0);
    virtual ~cusolve_sparse(){
	info2("cusolve_sparse::destructor\n");
    }
    virtual void R(curcell **out, float beta, 
		   const curcell *xin, float alpha, stream_t &stream){
	CR(out, beta, xin, alpha, stream);
    }
    virtual void L(curcell **out, float beta, 
		   const curcell *xin, float alpha, stream_t &stream){
	CL(out, beta, xin, alpha, stream);
    }
    virtual void Rt(curcell **out, float beta, 
		    const curcell *xin, float alpha, stream_t &stream){
	CR.Trans(out, beta, xin, alpha, stream);
    }
};

class cusolve_cbs:public cusolve_l{
protected:
    cusp *Cl;
    int  *Cp;
    curmat *Up;
    curmat *Vp;
    curcell *rhs; 
    curmat *y;
    curmat *Vr;
public:
    cusolve_cbs(spchol *_C=0, dmat *_Up=0, dmat *_Vp=0);
    virtual ~cusolve_cbs(){
	if(!this) return;
	delete rhs;
	delete Vr;
	delete Cl;
	delete Up;
	delete Vp;
	delete y;
	info2("cusolve_cbs::destructor\n");
    }
    void chol_solve(float *out, const float *in,  stream_t &stream);
    virtual float solve(curcell **xout, const curcell *xin, cusolve_r *Rfun, stream_t &stream);
};

class cusolve_svd:public cusolve_l{
    curcell *rhs;
    curmat *LI;//Inversion of left hand operator
public:
    cusolve_svd(dmat *_LI=0):rhs(NULL),LI(NULL){
	if(_LI) cp2gpu(&LI, _LI);
    }
    virtual ~cusolve_svd(){
	info2("cusolve_svd::destructor\n");
    }
    virtual float solve(curcell **xout, const curcell *xin, cusolve_r *Rfun, stream_t &stream){
	Rfun->R(&rhs, 0, xin, 1, stream);
	if(!*xout) *xout=curcellnew(xin);
	curmv((*xout)->m->p, 0, LI, rhs->m->p, 'n', 1, stream);
	return 0;
    }
}; 

class cusolve_mvm:public cusolve_l{
    curmat *MVM;
public:
    cusolve_mvm(dmat *mvm=0):MVM(0){
	if(!mvm) cp2gpu(&MVM, mvm);
    }
    ~cusolve_mvm(){
	if(!this) return;
	delete MVM;
    }
    virtual float solve(curcell **xout, const curcell *xin, cusolve_r *Rfun, stream_t &stream){
	assert(*xout && (*xout)->m);
	(void)Rfun;
	curmv((*xout)->m->p, 0., MVM, xin->m->p, 'n', 1., stream);
	return 0;
    }
};
#endif
