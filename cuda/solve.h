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

#ifndef AOS_CUDA_SOLVE_H
#define AOS_CUDA_SOLVE_H
#include "curmat.h"
#include "pcg.h"
namespace cuda_recon{
class cusolve_r{/*Interface for RHS*/
public:
    virtual void R(curcell &xout, Real beta, curcell &xin, Real alpha, stream_t &stream)=0;
    virtual void Rt(curcell &xout, Real beta, const curcell &xin, Real alpha, stream_t &stream)=0;
    virtual ~cusolve_r(){}
};

class cusolve_l{/*Interface for LHS*/
public:
    virtual Real solve(curcell &xout, const curcell &xin, stream_t &stream)=0;
    virtual ~cusolve_l(){}
};

class cusolve_cgpre{/*Interface for preconditioner*/
public:
    virtual void Pre(curcell &xout, const curcell &xin, stream_t &stream)=0;
    void operator()(curcell &xout, const curcell &xin, stream_t &stream){
	Pre(xout, xin, stream);
    }
    virtual ~cusolve_cgpre(){}
};

class cusolve_cg:public cusolve_l,nonCopyable{/*Implementes LHS with cg algorithm*/
    int maxit;       //number of iterations
    int warm_restart;//warm restart when first_run is 0.
    int first_run;   //indicate first run to temporarily use cold start rather than warm-restart
    cgtmp_t cgtmp;
    int id;
    static int counter;
protected:
    cusolve_cgpre *precond;
private:
    
public:
    cusolve_cg(int _maxit=0, int _warm_restart=0):maxit(_maxit), warm_restart(_warm_restart), precond(0), first_run(_warm_restart?1:0){
        id=counter; counter++;
    }
    ~cusolve_cg(){
    	delete precond;
    }
    /*Left hand side forward operaton*/
    virtual void L(curcell &xout, Real beta, const curcell &xin, Real alpha, stream_t &stream)=0;
    void operator()(curcell &xout, Real beta, const curcell &xin, Real alpha, stream_t &stream){
	    L(xout, beta, xin, alpha, stream);
    }
    virtual Real solve(curcell &xout, const curcell &xin, stream_t &stream);
    void Pre(curcell &xout, const curcell &xin, stream_t &stream){
	    if(precond){
	        (*precond)(xout, xin, stream);
	    }
    }
};

class cusolve_muv:public nonCopyable{
    cusp M;
    curmat U;
    curmat V;
    curmat Vx;
    int nx, ny, *nxs, *nys;
 public:
    cusolve_muv():nx(0),ny(0),nxs(0),nys(0){};
    void init(const muv_t *in);
    ~cusolve_muv(){
	delete[] nxs;
	delete[] nys;
    }
    void Forward(curcell &out, Real beta, const curcell &in, Real alpha, stream_t &stream);
    void operator()(curcell &out, Real beta, const curcell &in, Real alpha, stream_t &stream){
	Forward(out, beta, in, alpha, stream);
    }
    void Trans(curcell &out, Real beta, const curcell &in, Real alpha, stream_t &stream);
};

class cusolve_sparse:public cusolve_r,public cusolve_cg{
protected:
    cusolve_muv CR, CL;
public:
    cusolve_sparse(int _maxit, int _warm_restart, muv_t *_R, muv_t *_L);
    virtual void R(curcell &out, Real beta, 
		   curcell &xin, Real alpha, stream_t &stream){
	CR.Forward(out, beta, xin, alpha, stream);
    }
    virtual void L(curcell &out, Real beta, 
		   const curcell &xin, Real alpha, stream_t &stream){
	CL.Forward(out, beta, xin, alpha, stream);
    }
    virtual void Rt(curcell &out, Real beta, 
		    const curcell &xin, Real alpha, stream_t &stream){
	CR.Trans(out, beta, xin, alpha, stream);
    }
};

class cusolve_cbs:public cusolve_l,nonCopyable{
protected:
    cusp Cl;      //Lower diagonal sparse matrix L. M=LL';
    cuimat Cp;//Permutation matrix.
    curmat Up;   //Low rank
    curmat Vp;   //Low rank
    curmat y;    //Temporary data.
    curmat Vr;   //Temporary data.
public:
    cusolve_cbs(spchol *_C, dmat *_Up, dmat *_Vp);
    void chol_solve(Real *out, const Real *in,  stream_t &stream);
    virtual Real solve(curcell &xout, const curcell &xin, stream_t &stream);
};

class cusolve_mvm:public cusolve_l,nonCopyable{
    curmat M;
public:
    cusolve_mvm(dmat *_M);
    cusolve_mvm(curmat &_M){
	M=_M;
    }
    virtual Real solve(curcell &xout, const curcell &xin, stream_t &stream){
	if(!xout) xout=xin.New();
	curmv(xout.M()(), 0., M, xin.M()(), 'n', 1., stream);
	return 0;
    }
};
};//namespace
#endif
