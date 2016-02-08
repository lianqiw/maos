/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_LIB_MUV_H
#define AOS_LIB_MUV_H
#include "../math/mathdef.h"
#include "pcg.h"
/**
   \file muv.h

   Decomposes a matrix into sparse and low rank terms and do math
   operations. 
*/
/**
   Avoid function casting. It will hide data safety check and hide bugs.
*/
typedef void (*EXFUN) (dcell **xout, const void *A, const dcell *xin, const double alpha, int xb, int yb);
/**
   Decompose an operator into a sparse operator and optional low rank
   terms. M-U*V'; M is usually symmetrical.
 */
typedef struct MUV_T{
    /*The forward computation can be done by the following 3 matrices, together with exfun, extra (data) */
    cell *M;       /**<block sparse matrix or dense matrix*/
    dcell  *U;     /**<low rank terms U*/
    dcell  *V;     /**<low rank terms V*/
    EXFUN exfun;   /**<Optionally attach an extra function that applies extra data*/
    void *extra;   /**<Data used by fun to apply: (*exfun)(y,extra,x,alpha) to compute*/
    /*And the following are the inversion of above */
    spchol *C;     /**<Cholesky factor.*/
    dmat  *Up;     /**<M\U*/
    dmat  *Vp;     /**<M\[V*(I-Up'*V)^{-1}]*/
    dmat  *MI;     /**<Inverse of M via svd*/
    /*The following are for Block Gaussian Seidel Method. */
    dcell *MIB;    /**<Inverse of each diagonal element of M via svd*/
    spchol**CB;     /**<Cholesky factor of each diagonal element of M*/
    dcell *UpB;    /**<Up for each diagonal component*/
    dcell *VpB;    /**<Vp for each diagonal component*/
    int nb;        /**<Number of blocks in CB;*/
    /*The operation of M can be done with the folloing function and data */
    CGFUN Mfun;    /**<Do M*x with a function*/
    CGFUN Mtfun;   /**<Do M'*x with a function*/
    void *Mdata;   /**<Parameter to Mfun*/
    /*For CG purpose, the preconditioner function and data */
    PREFUN pfun;   /**<The preconditioner function*/
    void *pdata;   /**<The precondtioner data*/
    /*For solving the linear problem, a few arguments. */
    int alg;       /**<The algorithm: 0: CBS, 1: CG, 2: SVD*/
    int bgs;       /**<Whether use BGS (Block Gauss Seidel) method, and how many iterations*/
    int warm;      /**<Whether use warm restart*/
    int maxit;     /**<How many iterations*/
}MUV_T;

void muv(dcell **xout, const void *A, const dcell *xin, const double alpha);
void muv_trans(dcell **xout, const void *A, const dcell *xin, const double alpha);
void muv_ib(dcell **xout, const void *A, const dcell *xin, const double alpha);
void muv_direct_solve(dcell **xout, const MUV_T *A, dcell *xin);
void muv_direct_solve_mat(dmat **xout, const MUV_T *A, dmat *xin);
void muv_direct_prep(MUV_T *muv, double svd);
void muv_direct_free(MUV_T *muv);
void muv_direct_diag_solve(dmat **xout, const MUV_T *A, dmat *xin, int ib);
void muv_bgs_solve(dcell **px, const MUV_T *A, const dcell *b);
double muv_solve(dcell **px, const MUV_T *L, const MUV_T *R, dcell *b);
void* muv_direct_spsolve(const MUV_T *A, const dsp *xin);
void muv_direct_diag_prep(MUV_T *muv, double svd);
void muv_direct_diag_free(MUV_T *muv);
void muv_free(MUV_T *A);
#endif
