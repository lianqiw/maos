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
/**
   \file lib/muv.c 

   Decomposes a matrix into sparse and low rank terms and do math
   operations.  Apply M - u*v'*/


/**
   Apply the sparse plug low rand compuation to xin without scaling of alpha:
   \f$xout=(muv.M-muv.U*muv.V')*xin*alpha\f$; U,V are low rank.  */
void muv(dcell **xout, const MUV_T *A, const dcell *xin, const double alpha){
    spcellmulmat(xout, A->M, xin, alpha);
    dcell *tmp=NULL;
    dcellmm(&tmp,A->V, xin, "tn", -1.);
    dcellmm(xout,A->U, tmp, "nn", alpha);
    dcellfree(tmp);
    if(A->extra){
	A->exfun(xout, A->extra, xin, alpha);
    }
}

/**
   Apply cholesky backsubstitution solve to xin to get xout. Create xout is not
   exist already. The cholesky factors are prepared in setup_recon.c and stored
   in ((MUV_T*)A)->C */

void muv_direct_solve(dmat **xout, const MUV_T *A, const dmat *xin){
    dzero(*xout);
    if(A->MI){
	dmm(xout, A->MI, xin, "nn", 1);
    }else{
	chol_solve(xout,A->C,xin);
    }
    dmat *tmp=NULL;
    dmm(&tmp,A->Vp,xin,"tn",-1);
    dmm(xout,A->Up,tmp,"nn",1);
    dfree(tmp);
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
  Cholesky factorization (svd=0) or svd (svd>0) on the sparse matrix and its low
  rank terms.  if svd is less than 1, it is also used as the threshold in SVD
  inversion.*/
void muv_direct_prep(MUV_T *A, double svd){ 
    if(!A->M) error("M has to be none NULL\n");
    info2("muv_direct_prep: (svd=%g)", svd);    
    TIC;tic;
    muv_direct_free(A);
    dsp *muvM=spcell2sp(A->M);
    if(svd>0){
	spfull(&A->MI, muvM, 1);
	if(svd<1){
	    dsvd_pow(A->MI, -1, 1, svd);
	}else{
	    dsvd_pow(A->MI, -1, 1, 2e-4);
	}
    }else{
	A->C=chol_factorize(muvM);
    }
    spfree(muvM);
    if(A->U){
	dmat *U=dcell2m(A->U);
	A->Up=NULL;
	if(A->MI){
	    dmm(&A->Up, A->MI, U, "nn", 1);
	}else{
	    chol_solve(&(A->Up),A->C,U);
	}
	dfree(U);
	dmat *V=dcell2m(A->V);
	dmat *UpV=NULL;//UpV=I-Up'*V
	dmm(&UpV,A->Up,V,"tn",-1);
	for(unsigned long ii=0; ii<UpV->ny; ii++){
	    UpV->p[ii+ii*UpV->nx]+=1;
	}
	dinv_inplace(UpV);//Upv=(Upv)^-1
	dmat *VI=NULL;//VI=V*UpV
	dmm(&VI,V,UpV,"nn",-1);
	dfree(UpV);
	dfree(V);
	if(A->MI){
	    dmm(&A->Vp, A->MI, VI, "nn", 1);
	}else{
	    chol_solve(&(A->Vp),A->C,VI);
	}
	dfree(VI);
    }
    toc2("done.");
}

/**
   Free cholesky decompositions.
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
}
