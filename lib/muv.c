/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "muv.h"
/**
   Decomposes a matrix into sparse and low rank terms and do math
   operations.

   The forward operations can be represented by a matrix M and two low rank terms as
   y = ( M - u*v' ) *x , or by a function Mfun and corresponding argument Mdata.

   This can also be used to solve linear equations using muv_solve() or muv_bgs_solve().

*/
/*
  First, forward operations.
*/

/**
   Apply the sparse plus low rank compuation to xin with scaling of alpha:
   \f$xout=(A.M-A.U*A.V')*xin*alpha\f$; U,V are low rank.  */
void muv(dcell** xout, const void* A_, const dcell* xin, const real alpha){
	const muv_t* A=(const muv_t*)A_;//A_ is declared void for cg to use without casting.
	if(A->M&&xin){
		dcellmm_cell(xout, A->M, xin, "nn", alpha);
		if(A->U&&A->V){
			dcell* tmp=NULL;
			dcellmm(&tmp, A->V, xin, "tn", -1.);
			dcellmm(xout, A->U, tmp, "nn", alpha);
			dcellfree(tmp);
		}
		if(A->extra){
			A->exfun(xout, A->extra, xin, alpha, -1, -1);
		}
	} else{
		A->Mfun(xout, A->Mdata, xin, alpha);
	}
}
/**
   Apply the transpose of operation muv(); Apply the sparse plug low rand
   compuation to xin with scaling of alpha: \f$xout=(A.M-A.V*A.U')*xin*alpha\f$;
   U,V are low rank.  */
void muv_trans(dcell** xout, const void* A_, const dcell* xin, const real alpha){
	const muv_t* A=(const muv_t*)A_;
	if(A->M){
		if(!xin) return;
		dcellmm_cell(xout, A->M, xin, "tn", alpha);
		if(A->U&&A->V){
			dcell* tmp=NULL;
			dcellmm(&tmp, A->U, xin, "tn", -1.);
			dcellmm(xout, A->V, tmp, "nn", alpha);
			dcellfree(tmp);
		}
		if(A->extra){
			error("Need to implement this\n");
			A->exfun(xout, A->extra, xin, alpha, -1, -1);
		}
	} else if(A->Mtfun){
		A->Mtfun(xout, A->Mdata, xin, alpha);
	} else{
		error("Please assign Mtfun for M' operation\n");
	}
}
/*
  Wrap of data for call with muv_ib.
 */
typedef struct{
	const muv_t* A;
	int xb;
	int yb;
}MUV_IB_T;
/**
   Apply the sparse plus low rand compuation to xin with scaling of alpha for block (xb, yb):
   \f$xout_x=(A.M_xy-A.U_xi*A.V_yi')*xin_y*alpha\f$; U,V are low rank. */
void muv_ib(dcell** xout, const void* B, const dcell* xin, const real alpha){
	const MUV_IB_T* C=(const MUV_IB_T*)B;
	int xb=C->xb;
	int yb=C->yb;
	const muv_t* A=C->A;
	assert(NY(xin)==1);/*if this is not true, make a loop here. */
	if(xb==-1||yb==-1){
		muv(xout, A, xin, alpha);
		return;
	}
	assert(A->M);/*Don't support A->Mfun yet.  */
	if(!*xout){
		*xout=dcellnew(NX(A->M), NY(xin));
	}
	if(P(A->M, xb, yb)->id==M_REAL){
		dmm(&P(*xout,xb), 1, dmat_cast(P(A->M, xb, yb)), P(xin,yb), "nn", alpha);
	} else{
		dspmm(&P(*xout,xb), dsp_cast(P(A->M, xb, yb)), P(xin,yb), "nn", alpha);
	}
	if(A->V&&A->U){
		dcell* AV=A->V;
		dcell* AU=A->U;
		for(int jb=0; jb<NY(A->U); jb++){
			dmat* tmp=NULL;
			dmm(&tmp, 0, P(AV, yb, jb), P(xin,yb), "tn", -1);
			dmm(&P(*xout,xb), 1, P(AU, xb, jb), tmp, "nn", alpha);
			dfree(tmp);
		}
	}
	if(A->extra){
		A->exfun(xout, A->extra, xin, alpha, xb, yb);
	}
}

/*
  A few routines to do Cholesky deocomposition or SVD.
*/
/**
   Prepare the low rank terms:
   \f[
   Up=M^{-1}*U,
   Vp=M^{-1}*[V*(I-Up'*V)^{-1}]
   \f]
   Where M^{-1} is either using MI (svd) or C (CBS)
*/
static void muv_direct_prep_lowrank(dmat** Up, dmat** Vp, spchol* C, dmat* MI, dmat* U, dmat* V){
	if(MI){
		dmm(Up, 1, MI, U, "nn", 1);
	} else{
		chol_solve(Up, C, U);
	}
	dmat* UpV=NULL;
	/*UpV=(I-Up'*V)^{-1} */
	dmm(&UpV, 0, *Up, V, "tn", -1);
	for(long ii=0; ii<NY(UpV); ii++){
		P(UpV, ii, ii)+=1;
	}
	dinv_inplace(UpV);
	dmat* VI=NULL;
	/*VI=V*UpV */
	dmm(&VI, 0, V, UpV, "nn", -1);
	dfree(UpV);
	if(MI){
		dmm(Vp, 1, MI, VI, "nn", 1);
	} else{
		chol_solve(Vp, C, VI);
	}
	dfree(VI);
}

/**
   Cholesky factorization (svd=0) or svd (svd>0) on the sparse matrix and its low
   rank terms.  if svd is less than 1, it is also used as the threshold in SVD
   inversion.

   When there are low rank terms,
   \f[
	   Up=M^{-1}*U,
	   Vp=M^{-1}*[V*(I-Up'*V)^{-1}]
	\f]

*/
void muv_direct_prep(muv_t* A, real svd){
	int use_svd=fabs(svd)>0;
	if(!A->M) error("M has to be none NULL\n");
	if(A->extra) error("Direct solutions does not support extra/exfun\n");
	TIC;tic;
	muv_direct_free(A);
	if(use_svd){/*Do SVD */
		dfree(A->MI);
		A->MI=dcell2m_any(A->M);
		info("muv_direct_prep: (%s) on %ldx%ld array ", use_svd?"svd":"chol", NX(A->MI), NY(A->MI));
		dsvd_pow(A->MI, -1, svd<1?svd:2e-4, 0);
	} else{/*Do Cholesky decomposition. */
		dsp* muvM=dspcell2sp((const dspcell*)A->M);
		info("muv_direct_prep: (%s) on %ldx%ld array ", use_svd?"svd":"chol", NX(muvM), NY(muvM));
		A->C=chol_factorize(muvM);
		dspfree(muvM);
	}
	if(A->U&&A->V){
		dmat* U=dcell2m(A->U);
		dmat* V=dcell2m(A->V);
		if(U&&V&&NX(U)&&NY(U)){
			muv_direct_prep_lowrank(&A->Up, &A->Vp, A->C, A->MI, U, V);
			if(use_svd){
				dmm(&A->MI, 1, A->Up, A->Vp, "nt", -1);
				dfree(A->Up);
				dfree(A->Vp);
			}
		}
		dfree(V);
		dfree(U);
	}
	toc2("done.");
}
/**
  Cholesky factorization (svd=0) or svd (svd>0) on each diagonal element of the
  block sparse matrix and its low rank terms.  if svd is less than 1, it is also
  used as the threshold in SVD inversion.

  When there are low rank terms,
	  \f[
	  Up=M^{-1}U,
	  Vp=M^{-1}[V*(I-Up'*V)^{-1}]
	  \f]
  For each diagonal block.
*/
void muv_direct_diag_prep(muv_t* A, real svd){
	int use_svd=fabs(svd)>0;
	if(!A->M) error("M has to be none NULL\n");
	if(A->extra) error("Direct solutions does not support extra/exfun\n");
	assert(NX(A->M)==NY(A->M));
	info("muv_direct_prep_diag: (%s) ", use_svd?"svd":"chol");
	TIC;tic;
	muv_direct_diag_free(A);

	int nb=NX(A->M);
	A->nb=nb;
	if(use_svd){
		A->MIB=dcellnew(nb, 1);
	} else{
		A->CB=mycalloc(nb, spchol*);
	}
	for(int ib=0; ib<nb; ib++){/*Invert each diagonal block. */
		if(use_svd){
			if(P(A->M,ib,ib)->id==M_REAL){
				dcp(&P(A->MIB,ib), dmat_cast(P(A->M,ib,ib)));
			} else{
				dspfull(&P(A->MIB,ib), dsp_cast(P(A->M,ib,ib)), 'n', 1);
			}
			dsvd_pow(P(A->MIB,ib), -1, svd<1?svd:2e-4, 0);
		} else{
			A->CB[ib]=chol_factorize(dsp_cast(P(A->M,ib,ib)));
		}
	}
	if(A->U&&A->V){/*deal with low rank terms */
	/*First reduce U, V to column block vectors.  */
		dcell* U2=dcellreduce(A->U, 2);
		dcellfree(A->U); A->U=U2;
		dcell* V2=dcellreduce(A->V, 2);
		dcellfree(A->V); A->V=V2;

		A->UpB=dcellnew(NX(A->U), NY(A->U));
		A->VpB=dcellnew(NX(A->V), NY(A->V));
		for(int ib=0; ib<nb; ib++){
			muv_direct_prep_lowrank(&P(A->UpB,ib),
				&P(A->VpB,ib),
				A->CB?A->CB[ib]:NULL,
				A->MIB?P(A->MIB,ib):NULL,
				P(A->U,ib),
				P(A->V,ib));
		}
	}
	toc2("done.");
}
/*
  A few routines to apply CBS or SVD solver.
*/
/**
   Apply CBS or SVD multiply to xin to get xout. Create xout is not exist
   already.  xout = A^-1 * xin; xin and xout can be the same for in place operation.*/

void muv_direct_solve_mat(dmat** xout, const muv_t* A, dmat* xin){
	dmat* dotpt=NULL;
	if(A->Up&&A->Vp){
		dmm(&dotpt, 0, A->Vp, xin, "tn", -1);
	}
	if(A->MI){
		if(*xout==xin){
			dmat* tmp=NULL;
			dmm(&tmp, 0, A->MI, xin, "nn", 1);
			dcp(xout, tmp);
			dfree(tmp);
		} else{
			dmm(xout, 0, A->MI, xin, "nn", 1);
		}
	} else{
		chol_solve(xout, A->C, xin);
	}
	if(dotpt){
		dmm(xout, 1, A->Up, dotpt, "nn", 1);
		dfree(dotpt);
	}
}
/**
   Apply CBS or SVD multiply to square sparse matrix.
	   xout = A^-1 * xin *
	   (A^-1)^T.
   The output is dmat if using SVD, and dsp if using CBS.
 */
void* muv_direct_spsolve(const muv_t* A, const dsp* xin){
	void* xout=NULL;
	if(A->MI){
		dmat* x1=NULL;
		dmulsp(&x1, A->MI, xin, "nn", 1);
		dmm((dmat**)&xout, 0, x1, A->MI, "nt", 1);
		dfree(x1);
	} else{
		dsp* x1=chol_spsolve(A->C, xin);
		dsp* x1t=dsptrans(x1);
		dspfree(x1);
		dsp* x2=chol_spsolve(A->C, x1t);
		dspfree(x1t);
		xout=dsptrans(x2);
		dspfree(x2);
	}
	return xout;
}
/**
   Apply cholesky backsubstitution solve to xin to get xout. Create xout is not
   exist already.  xin and xout can be the same for in place operation. */

void muv_direct_diag_solve(dmat** xout, const muv_t* A, dmat* xin, int ib){
	if(!xin) return;
	dmat* dotpt=NULL;
	if(A->UpB&&A->VpB){
		dmm(&dotpt, 0, P(A->VpB,ib), xin, "tn", -1);
	}
	if(A->MIB){
		if(*xout==xin){
			dmat* tmp=NULL;
			dmm(&tmp, 0, P(A->MIB,ib), xin, "nn", 1);
			dcp(xout, tmp);
			dfree(tmp);
		} else{
			dmm(xout, 0, P(A->MIB,ib), xin, "nn", 1);
		}
	} else{
		chol_solve(xout, A->CB[ib], xin);
	}
	if(dotpt){
		dmm(xout, 1, P(A->UpB,ib), dotpt, "nn", 1);
		dfree(dotpt);
	}
}
/**
   convert the data from dcell to dmat and apply muv_direct_solve() . May be done in place.*/
void muv_direct_solve(dcell** xout, const muv_t* A, dcell* xin){
	if(!xin) return;
	if(NX(xin)*NY(xin)==1){/*there is only one cell. */
		if(!*xout) *xout=dcellnew(1, 1);
		muv_direct_solve_mat(&P(*xout,0), A, P(xin,0));
	} else{
		dmat* xin2=dcell2m(xin);
		muv_direct_solve_mat(&xin2, A, xin2);/*in place solve. */
		d2cell(xout, xin2, xin);/*xin is the reference for dimensions. copy data into xout. */
		dfree(xin2);
	}
}
/**
   solve A*x=b using the Block Gauss Seidel algorithm. The matrix is always
   assembled. We use routines from muv. to do the computation.  */
void muv_bgs_solve(dcell** px,    /**<[in,out] The output vector. input for warm restart.*/
	const muv_t* A,/**<[in] Contain info about the The left hand side matrix A*/
	const dcell* b /**<[in] The right hand side vector to solve*/
){
	if(!b) return;
	int nb=NX(b);
	assert(NY(b)==1);
	if(!*px){
		*px=dcellnew2(b); /*initialize the output array. */
	} else if(!A->warm){
		dcellzero(*px);
	}
	dcell* x0=*px; /*Just a convenient pointer. */
	dcell* c0=dcellnew2(b);/*Auxillary array */
	MUV_IB_T B={A,-1,-1};
	dcell* x0b=dcellnew(nb, 1);
	dcell* c0b=dcellnew(nb, 1);
	if(A->pfun){
		warning("We don't handle preconditioner yet\n");
	}
	for(int iter=0; iter<A->bgs; iter++){
		for(int ib=0; ib<nb; ib++){
			dcp(&P(c0,ib), P(b,ib));
			for(int jb=0; jb<nb; jb++){
				if(jb==ib) continue;
				B.xb=ib;
				B.yb=jb;
				muv_ib(&c0, &B, x0, -1);
			}
			switch(A->alg){
			case 0:
			case 2:
				muv_direct_diag_solve(&P(x0,ib), A, P(c0,ib), ib);
				break;
			case 1:/*Have to call CG, for this block. Embed this block into an empty block matrix. */
				P(x0b,ib)=P(x0,ib);
				P(c0b,ib)=P(c0,ib);
				B.xb=ib;
				B.yb=ib;
				pcg(&x0b, muv_ib, &B, NULL, NULL, c0b, A->warm, A->maxit);
				P(x0b,ib)=0;
				P(c0b,ib)=0;
				break;
			default:
				error("Invalid alg=%d\n", A->alg);
			}
		}
	}
	dcellfree(c0);
	dcellfree(c0b);
	dcellfree(x0b);
}
/**
   solve x=A^-1*B*b using algorithms depend on algorithm, wrapper.
*/
real muv_solve(dcell** px,    /**<[in,out] The output vector. input for warm restart.*/
	const muv_t* L,/**<[in] Contain info about the left hand side matrix A*/
	const muv_t* R,/**<[in] Contain info about the right hand side matrix B*/
	dcell* b       /**<[in] The right hand side vector to solve. null is special case*/
){
	dcell* rhs=NULL;
	real res=0;
	if(R){
		muv(&rhs, R, b, 1);
	} else{
		rhs=b;
	}
	if(L->bgs){
		muv_bgs_solve(px, L, rhs);
	} else{
		switch(L->alg){
		case 0:/*CBS */
		case 2:/*SVD */
			muv_direct_solve(px, L, rhs);
			break;
		case 1:/*CG */
			res=pcg(px, L->M?muv:L->Mfun, L->M?L:L->Mdata, L->pfun, L->pdata, rhs, L->warm, L->maxit);
			break;
		default:
			error("Invalid alg=%d\n", L->alg);
		}
	}
	if(rhs!=b) dcellfree(rhs);
	return res;
}

/**
   Free cholesky decompositions or SVD.
 */
void muv_direct_free(muv_t* A){
	if(A->C){
		chol_free(A->C);
		A->C=NULL;
	}
	if(A->MI){
		dfree(A->MI); A->MI=NULL;
	}
	dfree(A->Up);
	dfree(A->Vp);
}
/**
   Free cholesky decompositions or SVD.
 */
void muv_direct_diag_free(muv_t* A){
	if(A->CB){
		for(int i=0; i<A->nb; i++){
			chol_free(A->CB[i]);
		}
		free(A->CB); A->CB=NULL;
	}
	if(A->MIB){
		dcellfree(A->MIB);
	}
	dcellfree(A->UpB);
	dcellfree(A->VpB);
}
/**
   Free muv_t struct
*/
void muv_free(muv_t* A){
	cellfree(A->M);
	dcellfree(A->U);
	dcellfree(A->V);
	if(A->extra){
		free(A->extra);/*a struct contains pointers. */
	}
	muv_direct_free(A);
	muv_direct_diag_free(A);
	memset(A, 0, sizeof(muv_t));
}

