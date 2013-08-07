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
#include "solve.h"
namespace cuda_recon{
float cucg_t::solve(curcell **xout, const curcell *xin, stream_t &stream){
    float ans;
    if((ans=gpu_pcg(xout, this, precond, xin, &cgtmp,
		    warm_restart, maxit, stream))>1){
	cgtmp.count_fail++;
	warning2("CG %5d does not converge. ans=%g. maxit=%d\n", 
		 cgtmp.count_fail, ans, maxit);
    }
    return ans;
}

void cumuv_t::Forward(curcell **out, float beta, const curcell *in, float alpha, stream_t &stream){
    if(!M) error("M Can not be empty\n");
    if(!*out){
	*out=curcellnew(nx, 1, nxs, (int*)NULL);
    }else{
	curscale((*out)->m, beta, stream);
    }
    cuspmul((*out)->m->p, M, in->m->p, 1, 'n', alpha, stream);
    if(U && V){
	curmv(Vx->p, 0, V, in->m->p, 't', 1, stream);
	curmv((*out)->m->p, 1, U, Vx->p, 'n', -alpha, stream);
    }
}
void cumuv_t::Trans(curcell **out, float beta, const curcell *in, float alpha, stream_t &stream){
    if(!M) error("M Can not be empty\n");
    if(!*out){
	*out=curcellnew(ny, 1, nys, (int*)NULL);
    }else{
	curscale((*out)->m, beta, stream);
    }
    
    curscale((*out)->m, beta, stream);
    cuspmul((*out)->m->p, M, in->m->p, 1, 't', alpha, stream);
    if(U && V){
	curmv(Vx->p, 0, U, in->m->p, 't', 1, stream);
	curmv((*out)->m->p, 1, V, Vx->p, 'n', -alpha, stream);
    }
}
cumuv_t::cumuv_t(const MUV_T *in)
    :M(0),U(0),V(0),Vx(0),nx(0),ny(0),nxs(0),nys(0){
    if(!in) return;
    if(M || !in->M) error("in->M should not be NULL and M should be NULL\n");
    dsp *Mc=spcell2sp(in->M);
    dmat *Uc=dcell2m(in->U);
    dmat *Vc=dcell2m(in->V);
    nx=in->M->nx;
    ny=in->M->ny;
    nxs=(int*)malloc(sizeof(int)*nx);
    nys=(int*)malloc(sizeof(int)*ny);
    for(int i=0; i<nx; i++){
	nxs[i]=in->M->p[i]->m;
    }
    for(int i=0; i<ny; i++){
	nys[i]=in->M->p[i*in->M->nx]->n;
    }
    M=new cusp(Mc, 1);
    cp2gpu(&U, Uc);
    cp2gpu(&V, Vc);
    spfree(Mc); dfree(Uc); dfree(Vc);
    Vx=curnew(V->ny, 1);
}
cusolve_sparse::cusolve_sparse(int _maxit, int _warm_restart, MUV_T *_R, MUV_T *_L)
    :cucg_t(_maxit, _warm_restart),CR(NULL),CL(NULL){
    if(!_R) return;
    CR=new cumuv_t(_R);
    CL=new cumuv_t(_L);
}
cusolve_cbs::cusolve_cbs(spchol *_C, dmat *_Up, dmat *_Vp)
    :Vr(0),Cl(0),Cp(0),Up(0),Vp(0),y(0){
    if(!_C) return;
    chol_convert(_C, 0);
    Cl=new cusp(_C->Cl, 0);
    cp2gpu(&Cp, _C->Cp, _C->Cl->m);
    if(_Up){
	cp2gpu(&Up, _Up);
	cp2gpu(&Vp, _Vp);
    }
}
float cusolve_cbs::solve(curcell **xout, const curcell *xin, stream_t &stream){
    if(!*xout) *xout=curcellnew(xin);
    if(Cl->type==SP_CSC){
	chol_solve((*xout)->m->p, xin->m->p, stream);
    }else{
	error("To implemente\n");
    }
    if(Up){
	if(!Vr){
	    Vr=curnew(Vp->ny, 1);
	}
	curmv(Vr->p, 0, Vp, xin->m->p, 't', -1, stream);
	curmv((*xout)->m->p, 1, Up, Vr->p, 'n', 1, stream);
    }
    return 0;
}
/*solve in place*/
static __global__ void cuchol_solve_lower_do(float *restrict y, float *Cx, int *Cp, int *Ci, int n){
    int id=threadIdx.x;
    int nd=blockDim.x;
    extern __shared__ float sb[];
    __shared__ float val;
    /*first solve L\y*/
    
    for(int icol=0; icol<n; icol++){
	if(id==0){
	    y[icol]/=Cx[Cp[icol]];//divide by diagonal.
	    val=-y[icol];
	}
	__syncthreads();//this is necessary!
	for(int irow=Cp[icol]+1+id; irow<Cp[icol+1]; irow+=nd){
	    y[Ci[irow]]+=val*Cx[irow];
	}
	__syncthreads();
    }
    /*Next solve L'\y. Use reduction algorithm instead of atomic add.*/
    for(int icol=n-1; icol>-1; icol--){
	sb[id]=0;
	for(int irow=Cp[icol]+1+id; irow<Cp[icol+1]; irow+=nd){
	    sb[id]+=Cx[irow]*y[Ci[irow]];
	}
	for(int step=(blockDim.x>>1);step>0;step>>=1){
	    __syncthreads();
	    if(id<step){
		sb[id]+=sb[id+step];
	    }
	}
	if(id==0){
	    y[icol]=(y[icol]-sb[0])/Cx[Cp[icol]];
	}
	__syncthreads();//this is necessary!
    }
}
void cusolve_cbs::chol_solve(float *out, const float *in, stream_t &stream){
    if(!Cl || !Cp) error("Invalid\n");
    int n=Cl->nx;
    if(!y){
	y=curnew(Cl->nx, 1);
    }
    perm_f_do<<<DIM(n, 256),0,stream>>>(y->p, in, Cp, n);
    //only 1 block for synchronization.
    const int NTH=256;
    cuchol_solve_lower_do<<<1,NTH, NTH*sizeof(float),stream>>>(y->p, Cl->x, Cl->p, Cl->i, n); 
    perm_i_do<<<DIM(n, 256),0,stream>>>(out, y->p, Cp, n);
    cudaStreamSynchronize(stream);
}
}
