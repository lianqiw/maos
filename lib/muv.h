/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_LIB_muv_H
#define AOS_LIB_muv_H
#include "dmat.h"
#include "cmat.h"
#include "cell.h"
#include "dsp.h"
#include "chol.h"
#include "pcg.h"
/**
   Decompose an operator into a sparse operator and optional low rank
   terms. M-U*V';
 */
typedef struct MUV_T{
    spcell *M;     /**<block sparse matrix*/
    dcell  *U;     /**<low rank terms U*/
    dcell  *V;     /**<low rank terms V*/
    //dcell *UP;     /**<*/
    //dcell *VP;
    spchol *C;     //Cholesky factor.
    dsp   *Cs;       //Converted cholesky factor.
    long  *Cp;    //permutation vector for Cs.
    dmat  *Up;
    dmat  *Vp;

    CGFUN exfun;     /**<Optionally attach an extra function that applies extra data*/
    void *extra;     /**<Data used by fun to apply: (*exfun)(y,extra,x,alpha) to compute*/
}MUV_T;

void muv(dcell **xout, const MUV_T *A, const dcell *xin, const double alpha);
void muv_chol_solve_cell(dcell **xout, const MUV_T *A, const dcell *xin);
void muv_chol_solve(dmat **xout, const MUV_T *A, const dmat *xin);
void muv_chol_prep(MUV_T *muv);
void muv_chol_free(MUV_T *muv);
void muv_free(MUV_T *A);
#endif
