/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


int cusolve_cg::counter=0;
Real cusolve_cg::solve(curcell& xout, const curcell& xin, stream_t& stream){
	Real ans=0;
	cgtmp.count++;
	Real thres=MAX(cgtmp.residual*5,0.05);
	int restarted=0;
repeat:
	if(first_run){
		dbg("CG %5d(%d) %s with maxit scaled by %d.\n", cgtmp.count, id, first_run==2?"restart":"cold start", first_run+1);
	}
	ans=pcg(xout, this, precond, xin, cgtmp, first_run?0:cgwarm, 
			(first_run&&cgwarm)?(maxit*(1+first_run)):maxit, stream);
	first_run=0;
	if(cgwarm){
		if(restarted){
			info("CG %5d(%d) after restart %d times, residual=%.5f threshold is %.5f.\n", cgtmp.count, id, restarted, ans, thres);
		}
	
		if(ans>thres){//only check during warm_restart
			cgtmp.count_fail++;
			//failure is usually caused by a rare temporary memory corruption.
			info("CG %5d(%d) does not converge: residual=%.5f, threshold is %.5f.\n", 
				cgtmp.count, id, ans, thres);
			if(!restarted){
				info("CG %5d(%d) is restarted\n", cgtmp.count, id);
				restarted++;
				first_run=2;
				goto repeat;
			}
		}
		
		if(!disable_save&&cgtmp.count_fail<10&&ans>thres*5){
			info("CG result saved.\n");
			cuwrite(xin, stream, "cucg_solve_xin_%d", cgtmp.count_fail);
			cuwrite(xout, stream, "cucg_solve_xout_%d", cgtmp.count_fail);
			cuwrite(cgtmp.r0, stream, "cucg_solve_r0_%d", cgtmp.count_fail);
			cuwrite(cgtmp.z0, stream, "cucg_solve_z0_%d", cgtmp.count_fail);
			cuwrite(cgtmp.p0, stream, "cucg_solve_p0_%d", cgtmp.count_fail);
			cuwrite(cgtmp.Ap, stream, "cucg_solve_Ap_%d", cgtmp.count_fail);
			cuwrite(cgtmp.store, stream, "cucg_solve_store_%d", cgtmp.count_fail);
			cuwrite(cgtmp.diff, stream, "cucg_solve_diff_%d", cgtmp.count_fail);
		}
		
		//if(ans<thres){
		cgtmp.residual=cgtmp.residual*0.9+ans*0.1;
		//}
	}
	return ans;
}

void cusolve_muv::Forward(curcell& out, Real beta, const curcell& in, Real alpha, stream_t& stream){
	if(!out){
		out=curcell(nx, 1, nxs, (int*)NULL);
	} else{
		Scale(out.M(), beta, stream);
	}
	if(M){
		cuspmul(out.M(), M, in.M(), 1, 'n', alpha, stream);
	}else{
		cugemm(out.M(), (Real)1, Md, in.M(), "nn", alpha, stream);
	}
	if(U&&V){
		curmv(Vx(), 0, V, in.M()(), 't', 1, stream);
		curmv(out.M()(), 1, U, Vx(), 'n', -alpha, stream);
	}
}
void cusolve_muv::Trans(curcell& out, Real beta, const curcell& in, Real alpha, stream_t& stream){
	if(!out){
		out=curcell(ny, 1, nys, (int*)NULL);
	} else{
		Scale(out.M(), beta, stream);
	}

	Scale(out.M(), beta, stream);
	if(M){
		cuspmul(out.M(), M, in.M(), 1, 't', alpha, stream);
	}else{
		cugemm(out.M(), (Real)1, Md, in.M(), "tn", alpha, stream);
	}
	if(U&&V){
		curmv(Vx(), 0, U, in.M()(), 't', 1, stream);
		curmv(out.M()(), 1, V, Vx(), 'n', -alpha, stream);
	}
}
void cusolve_muv::init(const muv_t* in){
	if(!in||!in->M) {
		dbg("in or in.M is empty, do nothing\n");
		return;
	}
	dspcell *inM=NULL;
	if((inM=dspcell_cast(in->M))){
		dsp* Mc=dspcell2sp(inM);
		M=cusp(Mc);
		dspfree(Mc); 
	}else{
		dmat *MM=dcell2m(in->M);
		cp2gpu(Md, MM); dfree(MM);
	}
	dmat *Uc=dcell2m(in->U);
	cp2gpu(U, Uc);dfree(Uc); 
	dmat* Vc=dcell2m(in->V);
	cp2gpu(V, Vc);dfree(Vc);
	Vx=curmat(V.Ny(), 1);
	
	nx=in->M->nx;
	ny=in->M->ny;
	nxs=new int[nx];
	nys=new int[ny];
	for(int i=0; i<nx; i++){
		nxs[i]=NX(in->M,i,0);
	}
	for(int i=0; i<ny; i++){
		nys[i]=NY(in->M,0,i);
	}
}

cusolve_sparse::cusolve_sparse(int _maxit, int _cgwarm, muv_t* _R, muv_t* _L)
	:cusolve_cg(_maxit, _cgwarm){
	CR.init(_R);
	CL.init(_L);
}
cusolve_cbs::cusolve_cbs(spchol* _C, const dmat* _Up, const dmat* _Vp){
	if(!_C){
		error("C cannot be empty\n");
	}
	chol_convert(_C, 1);
	if(_C->Cl){
		Cl=cusp(_C->Cl, 0, 0);
	} else{
		Cl=cusp(_C->Cu, 0, 1);
	}
	cp2gpu(Cp, _C->Cp, Cl.Nx(), 1);
	if(_Up){
		cp2gpu(Up, _Up);
		cp2gpu(Vp, _Vp);
	}
}
/**
 * @brief Apply the following: xout=solve(xin)-Up*(Vp'*xin). Notice that, in
 * single precision mode, the error is large due to accumulation.
 *
 * @param xout 
 * @param xin 
 * @param stream 
 * @return Real 
 */
Real cusolve_cbs::solve(curcell& xout, const curcell& xin, stream_t& stream){
	if(!xout) xout=xin.New();
	if(Up){//Apply Vp first in case xin and xout is the same
		if(!Vr)	Vr=curmat(Vp.Ny(), 1);
		curmv(Vr(), 0, Vp, xin.M()(), 't', -1, stream);
	}
	if(Cl.Type()==SP_CSC){
		chol_solve(xout.M()(), xin.M()(), stream);
	} else{
		error("To implement.\n");
	}
	if(Up){
		curmv(xout.M()(), 1, Up, Vr(), 'n', 1, stream);
	}
	return 0;
}
/*solve in place*/
static __global__ void cuchol_solve_lower_do(Real* restrict y, Real* Cx, Spint* Cp, Spint* Ci, int n){
	int id=threadIdx.x;
	int nd=blockDim.x;
	extern __shared__ Real sb[];
	__shared__ Real val;
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
void cusolve_cbs::chol_solve(Real* out, const Real* in, stream_t& stream){
	if(!Cl||!Cp) error("Invalid\n");
	int n=Cl.Nx();
	if(!y){
		y=curmat(Cl.Nx(), 1);
	}
	perm_f_do<<<DIM(n, 256), 0, stream>>>(y(), in, Cp(), n);
	//only 1 block for synchronization. //todo: improve the implementation.
	const int NTH=256;
	cuchol_solve_lower_do<<<1, NTH, NTH*sizeof(Real), stream>>>(y(), Cl.Px(), Cl.Pp(), Cl.Pi(), n);
	perm_i_do<<<DIM(n, 256), 0, stream>>>(out, y(), Cp(), n);
	//CUDA_SYNC_STREAM;
}
cusolve_mvm::cusolve_mvm(dmat* _M){
	cp2gpu(M, _M);
}
