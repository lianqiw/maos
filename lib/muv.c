/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "sys/thread.h"
/**
   \file lib/muv.c 

   Decomposes a matrix into sparse and low rank terms and do math
   operations. 
   
   The forward operations can be represented by a matrix M and two low rank terms as
   y = ( M - u*v' ) *x , or by a function Mfun and corresponding argument Mdata.
   
   This can also be used to solve linear equations using muv_solve() or muv_bgs_solve().
   
*/


/**
   Apply the sparse plug low rand compuation to xin with scaling of alpha:
   \f$xout=(A.M-A.U*A.V')*xin*alpha\f$; U,V are low rank.  */
void muv(dcell **xout, const void *B, const dcell *xin, const double alpha){
    if(!xin) return;
    const MUV_T *A = B;
    if(A->M){
	spcellmulmat(xout, A->M, xin, alpha);
	dcell *tmp=NULL;
	dcellmm(&tmp,A->V, xin, "tn", -1.);
	dcellmm(xout,A->U, tmp, "nn", alpha);
	dcellfree(tmp);
	if(A->extra){
	    A->exfun(xout, A->extra, xin, alpha, -1, -1);
	}
    }else{
	A->Mfun(xout, A->Mdata, xin, alpha);
    }
}
/*
  Wrap of data for call with muv_ib.
 */
typedef struct{
    const MUV_T *A;
    int xb;
    int yb;
}MUV_IB_T;
/**
   Apply the sparse plug low rand compuation to xin with scaling of alpha for block (xb, yb):
   \f$xout_x=(A.M_xy-A.U_xi*A.V_yi')*xin_y*alpha\f$; U,V are low rank. */
void muv_ib(dcell **xout, const void *B, const dcell *xin, const double alpha){
    const MUV_IB_T *C=B;
    int xb=C->xb;
    int yb=C->yb;
    const MUV_T *A=C->A;
    assert(xin->ny==1);//if this is not true, make a loop here.
    if(xb==-1 || yb==-1){
	muv(xout, A, xin, alpha);
	return;
    }
    assert(A->M);//Don't support A->Mfun yet. 
    PDSPCELL(A->M, AM);
    if(!*xout){
	*xout=dcellnew(A->M->nx, xin->ny);
    }
    spmulmat(&(*xout)->p[xb], AM[yb][xb], xin->p[yb], alpha);

    PDCELL(A->V, AV);
    PDCELL(A->U, AU);
    for(int jb=0; jb<A->U->ny; jb++){
	dmat *tmp=NULL;
	dmm(&tmp, AV[jb][yb], xin->p[yb], "tn", -1);
	dmm(&(*xout)->p[xb], AU[jb][xb], tmp, "nn", alpha);
	dfree(tmp);
    }
    if(A->extra){
	A->exfun(xout, A->extra, xin, alpha, xb, yb);
    }
}

/**
   Prepare the low rank terms:
   Up=M^{-1}*U, 
   Vp=M^{-1}*[V*(I-Up'*V)^{-1}]
   Where M^{-1} is either using MI (svd) or C (CBS)
*/
static void muv_direct_prep_lowrank(dmat **Up, dmat **Vp, spchol *C, dmat *MI, dmat *U, dmat *V){
    if(MI){
	dmm(Up, MI, U, "nn", 1);
    }else{
	chol_solve(Up,C,U);
    }
    dmat *UpV=NULL;
    //UpV=(I-Up'*V)^{-1}
    dmm(&UpV,*Up,V,"tn",-1);
    for(unsigned long ii=0; ii<UpV->ny; ii++){
	UpV->p[ii+ii*UpV->nx]+=1;
    }
    dinv_inplace(UpV);
    dmat *VI=NULL;
    //VI=V*UpV
    dmm(&VI,V,UpV,"nn",-1);
    dfree(UpV);
    if(MI){
	dmm(Vp, MI, VI, "nn", 1);
    }else{
	chol_solve(Vp,C,VI);
    }
    dfree(VI);
}

/**
  Cholesky factorization (svd=0) or svd (svd>0) on the sparse matrix and its low
  rank terms.  if svd is less than 1, it is also used as the threshold in SVD
  inversion.

  When there are low rank terms, 
  Up=M^{-1}*U, 
  Vp=M^{-1}*[V*(I-Up'*V)^{-1}]

*/
void muv_direct_prep(MUV_T *A, double svd){ 
    if(!A->M) error("M has to be none NULL\n");
    if(A->extra) error("Direct solutions does not support extra/exfun\n");
    info2("muv_direct_prep: (%s) ", svd?"svd":"chol");
    TIC;tic;
    muv_direct_free(A);
    dsp *muvM=spcell2sp(A->M);
    if(svd>0){//Do SVD
	spfull(&A->MI, muvM, 1);
	if(svd<1){//use svd as threashold
	    dsvd_pow(A->MI, -1, 1, svd);
	}else{//use a threshold good for lsr.
	    dsvd_pow(A->MI, -1, 1, 2e-4);
	}
    }else{//Do Cholesky decomposition.
	A->C=chol_factorize(muvM);
    }
    spfree(muvM);
    if(A->U && A->V){
	dmat *U=dcell2m(A->U);
	dmat *V=dcell2m(A->V);
	muv_direct_prep_lowrank(&A->Up, &A->Vp, A->C, A->MI, U, V);
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
  Up=M\U, 
  Vp=M\[V*(I-Up'*V)^{-1}]
  For each diagonal block.
*/
void muv_direct_diag_prep(MUV_T *A, double svd){ 
    if(!A->M) error("M has to be none NULL\n");
    if(A->extra) error("Direct solutions does not support extra/exfun\n");
    assert(A->M->nx == A->M->ny);
    info2("muv_direct_prep_diag: (%s) ", svd?"svd":"chol");
    TIC;tic;
    muv_direct_diag_free(A);
    
    int nb=A->M->nx;
    A->nb=nb;
    if(svd>0){
	A->MIB=dcellnew(nb,1);
    }else{
	A->CB=calloc(nb, sizeof(spchol*));
    }
    for(int ib=0; ib<nb; ib++){//Invert each diagonal block.
	dsp *muvM=A->M->p[ib+ib*nb];
	if(svd>0){
	    spfull(&A->MIB->p[ib], muvM, 1);
	    if(svd<1){
		dsvd_pow(A->MIB->p[ib], -1, 1, svd);
	    }else{
		dsvd_pow(A->MIB->p[ib], -1, 1, 2e-4);
	    }
	}else{
	    A->CB[ib]=chol_factorize(muvM);
	}
	//Don't free muvM.
    }
    if(A->U && A->V){//deal with low rank terms
	//First reduce U, V to column block vectors. 
	dcell *U2=dcellreduce(A->U, 2);
	dcellfree(A->U); A->U=U2;
	dcell *V2=dcellreduce(A->V, 2);
	dcellfree(A->V); A->V=V2;

	A->UpB=dcellnew(A->U->nx, A->U->ny);	
	A->VpB=dcellnew(A->V->nx, A->V->ny);	
	for(int ib=0; ib<nb; ib++){
	    muv_direct_prep_lowrank(&A->UpB->p[ib], 
				    &A->VpB->p[ib], 
				    A->CB ? A->CB[ib] : NULL, 
				    A->MIB ? A->MIB->p[ib] : NULL, 
				    A->U->p[ib],
				    A->V->p[ib]);
	}
    }
    toc2("done.");
}

/**
   Apply cholesky backsubstitution solve to xin to get xout. Create xout is not
   exist already.  */

void muv_direct_solve(dmat **xout, const MUV_T *A, const dmat *xin){
    dzero(*xout);
    if(A->MI){
	dmm(xout, A->MI, xin, "nn", 1);
    }else{
	chol_solve(xout,A->C,xin);
    }
    if(A->Up && A->Vp){
	dmat *tmp=NULL;
	dmm(&tmp,A->Vp,xin,"tn",-1);
	dmm(xout,A->Up,tmp,"nn",1);
	dfree(tmp);
    }
}
/**
   Apply cholesky backsubstitution solve to xin to get xout. Create xout is not
   exist already.  */

void muv_direct_diag_solve(dmat **xout, const MUV_T *A, const dmat *xin, int ib){
    dzero(*xout);
    if(A->MIB){
	dmm(xout, A->MIB->p[ib], xin, "nn", 1);
    }else{
	chol_solve(xout,A->CB[ib],xin);
    }
    if(A->UpB && A->VpB){
	dmat *tmp=NULL;
	dmm(&tmp,A->VpB->p[ib],xin,"tn",-1);
	dmm(xout,A->UpB->p[ib],tmp,"nn",1);
	dfree(tmp);
    }
}
/**
   convert the data from dcell to dmat and apply muv_direct_solve() */
void muv_direct_solve_cell(dcell **xout, const MUV_T *A, const dcell *xin){
    if(xin->nx*xin->ny==1){//there is only one cell.
	if(!*xout) *xout=dcellnew(1,1);
	muv_direct_solve(&((*xout)->p[0]), A, xin->p[0]);
    }else{
	dmat *xin2=dcell2m(xin);
	dmat *xout2=NULL;
	muv_direct_solve(&xout2, A, xin2);
	dfree(xin2);
	d2cell(xout,xout2,xin);//xin is the reference for dimensions. copy data into xout.
	dfree(xout2);
    }
}
/**
   solve A*x=b using the Block Gauss Seidel algorithm. The matrix is always
   assembled. We use routines from muv. to do the computation.  */
void muv_bgs_solve(dcell **px,    /**<[in,out] The output vector. input for warm restart.*/
		   const MUV_T *A,/**<[in] Contain info about the The left hand side matrix A*/
		   const dcell *b/**<[in] The right hand side vector to solve*/ 
		   ){
    if(!b) return;
    int nb=b->nx;
    assert(b->ny==1);
    if(!*px){
	*px=dcellnew2(b); //initialize the output array.
    }else if(!A->warm){
	dcellzero(*px);
    }
    dcell *x0=*px; //Just a convenient pointer.
    dcell *c0=dcellnew2(b);//Auxillary array
    MUV_IB_T B={A,-1,-1};
    dcell *x0b=dcellnew(nb,1);
    dcell *c0b=dcellnew(nb,1);
    if(A->pfun){
	warning("We don't handle preconditioner yet\n");
    }
    for(int iter=0; iter<A->bgs; iter++){
	for(int ib=0; ib<nb; ib++){
	    dcp(&c0->p[ib], b->p[ib]);
	    for(int jb=0; jb<nb; jb++){
		if(jb==ib) continue;
		B.xb=ib;
		B.yb=jb;
		muv_ib(&c0, &B, x0, -1);
	    }
	    switch(A->alg){
	    case 0:
	    case 2:
		muv_direct_diag_solve(&x0->p[ib], A, c0->p[ib],ib);
		break;
	    case 1://Have to call CG, for this block. Embed this block into an empty block matrix.
		x0b->p[ib]=x0->p[ib];
		c0b->p[ib]=c0->p[ib];
		B.xb=ib;
		B.yb=ib;
		pcg(&x0b, muv_ib, &B, NULL, NULL, c0b, A->warm, A->maxit);
		x0b->p[ib]=0;
		c0b->p[ib]=0;
		break;
	    default:
		error("Invalid alg=%d\n", A->alg);
	    }
	}
    }
    dcellfree(c0);
}
/**
   solve A*x=b using algorithms depend on algorithm
*/
void muv_solve(dcell **px,    /**<[in,out] The output vector. input for warm restart.*/
	       const MUV_T *A,/**<[in] Contain info about the The left hand side matrix A*/
	       const dcell *b/**<[in] The right hand side vector to solve*/ 
	       ){
    if(A->bgs){
	muv_bgs_solve(px, A, b);
    }else{
	switch(A->alg){
	case 0://CBS
	case 2://SVD
	    muv_direct_solve_cell(px, A, b);
	    break;
	case 1://CG
	    pcg(px, A->M?muv:A->Mfun, A->M?A:A->Mdata, A->pfun, A->pdata, b, A->warm, A->maxit);
	    break;
	default:
	    error("Invalid alg=%d\n", A->alg);
	}
    }
}

/**
   Free cholesky decompositions or SVD.
 */
void muv_direct_free(MUV_T *A){
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
void muv_direct_diag_free(MUV_T *A){
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
   Free MUV_T struct
*/
void muv_free(MUV_T *A){
    spcellfree(A->M);
    dcellfree(A->U);
    dcellfree(A->V);
    if(A->extra){
	free(A->extra);//a struct contains pointers.
    }
    muv_direct_free(A);
    muv_direct_diag_free(A);
    memset(A, 0, sizeof(MUV_T));
}

