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

#ifndef AOS_CELL_H
#define AOS_CELL_H
#include "type.h"

#define AOS_CELL_DEF(X,Y,T)\
X(cell) *X(cellnew)(const long nx, const long ny);\
X(cell) *X(cellnew2)(const X(cell) *A);\
void X(cellzero)(X(cell) *dc);\
void X(cellset)(X(cell)*dc, T alpha);\
void X(cellfree_do)(X(cell) *dc);\
X(cell) *X(celltrans)(const X(cell) *A);\
X(cell) *X(cellref)(const X(cell) *in);\
X(cell) *X(celldup)(const X(cell) *in);\
void X(cellcp)(X(cell)** out0, const X(cell) *in);\
double X(cellnorm2)(const X(cell) *in);\
void X(cellscale)(X(cell) *A, double w);\
X(cell) *X(cellreduce)(const X(cell) *A, int dim);\
X(cell) *X(cellcat)(const X(cell) *A, const X(cell) *B, int dim);\
X(cell) *X(cellcat_each)(const X(cell) *A, const X(cell) *B, int dim);\
void X(celldropempty)(X(cell) **A0, int dim);\
void X(celladd)(X(cell) **B0, double bc, const X(cell) *A,const double ac);\
T X(cellinn)(const X(cell)*A, const X(cell)*B);\
void X(cellcwm)(X(cell) *B, const X(cell) *A);\
void X(cellmm)(X(cell) **C0, const X(cell) *A, const X(cell) *B, const char trans[2], const double alpha);\
X(cell)* X(cellinvspd)(X(cell) *A);\
X(cell)* X(cellinv)(X(cell) *A);\
X(cell)* X(cellinvspd_each)(X(cell) *A);\
X(mat) *X(cell2m)(const X(cell) *A);\
X(cell)* X(2cellref)(const X(mat) *A, int*dims, int ndim);\
void X(2cell)(X(cell) **B, const X(mat) *A, const X(cell) *ref);\
void X(celldropzero)(X(cell) *B, double thres);\
double X(celldiff)(const X(cell) *A, const X(cell) *B);\
int X(cellclip)(X(cell) *Ac, double min, double max);\
void X(celltikcr)(X(cell) *A, double thres);\
X(cell) *X(cellpinv)(const X(cell) *A, const X(cell) *wt, const Y(spcell) *Wsp);\
void X(cellmulsp)(X(cell) **C0, const X(cell) *A, const Y(spcell) *B, double alpha);\
void X(celladdI)(X(cell) *A, double a);\
void X(cellsvd_pow)(X(cell) *A, double power);\
void X(cellcwpow)(X(cell)*A, double power);

#endif
