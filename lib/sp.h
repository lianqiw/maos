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

#ifndef AOS_SP_H
#define AOS_SP_H
#include "random.h"
#include "type.h"

#define spfree(A)      {spfree_do(A); A=NULL;}
#define spcellfree(A)  {spcellfree_do(A); A=NULL;}
#define PSPCELL(M,P)   dsp* (*restrict P)[(M)->nx]=(dsp*(*)[(M)->nx])(M)->p
#define PDSPCELL(M,P)  dsp* (*restrict P)[(M)->nx]=(dsp*(*)[(M)->nx])(M)->p

#define cspfree(A)     {cspfree_do(A); A=NULL;}
#define cspcellfree(A) {cspcellfree_do(A); A=NULL;}
#define PCSPCELL(M,P)  csp* (*restrict P)[(M)->nx]=(csp*(*)[(M)->nx])(M)->p

#define AOS_SP_DEF(X,Y,T) \
X(sp)* Y(spnew)(long nx, long ny, long nzmax) CHECK_UNUSED_RESULT;\
X(sp) *Y(spref)(X(sp) *A) CHECK_UNUSED_RESULT;\
X(sp) *Y(spdup)(const X(sp) *A) CHECK_UNUSED_RESULT;\
void Y(spmove)(X(sp) *A, X(sp) *res);\
X(sp) *Y(spnew2)(const X(sp) *A) CHECK_UNUSED_RESULT;\
X(sp)* Y(spnewrandu)(int nx, int ny, const T mean, double fill,rand_t *rstat) CHECK_UNUSED_RESULT;\
void Y(spsetnzmax)(X(sp) *sp, long nzmax);\
void Y(spfree_do)(X(sp) *sp);\
void Y(sparrfree)(X(sp) **sparr, int n);\
void Y(spdisp)(const X(sp) *sp);\
void Y(spcheck)(const X(sp) *sp);\
void Y(spscale)(X(sp) *A, const T beta);\
void Y(spcellscale)(Y(spcell) *A, const T beta);\
X(sp)* Y(spnewdiag)(long N, T *vec, T alpha) CHECK_UNUSED_RESULT;\
X(mat) *Y(spdiag)(const X(sp) *A) CHECK_UNUSED_RESULT;\
void Y(spmuldiag)(X(sp) *restrict A, const T* w, T alpha);\
void Y(spmulvec_thread)(T *restrict y, const X(sp) *A, const T * restrict x, T alpha, int nthread);	\
void Y(spmulvec)(T *restrict y, const X(sp) *A, const T * restrict x, T alpha);\
void Y(spmulcreal)(T *restrict y, const X(sp) *A, const dcomplex * restrict x, T alpha);\
void Y(sptmulvec)(T *restrict y, const X(sp) *A, const T * restrict x,const T alpha);\
void Y(sptmulvec_thread)(T *restrict y, const X(sp) *A, const T * restrict x,const T alpha, int nthread); \
void Y(spmulmat)(X(mat) **yout, const X(sp) *A, const X(mat) *x, const T alpha);\
void Y(sptmulmat)(X(mat) **yout, const X(sp) *A, const X(mat) *x, const T alpha);\
T Y(spwdinn)(const X(mat) *y, const X(sp) *A, const X(mat) *x) CHECK_UNUSED_RESULT;\
T Y(spcellwdinn)(const X(cell) *y, const Y(spcell) *A, const X(cell) *x) CHECK_UNUSED_RESULT;\
void Y(spcellmulmat)(X(cell) **C, const Y(spcell)*A, const X(cell)*B, const T alpha);\
void Y(sptcellmulmat)(X(cell) **C, const Y(spcell)*A, const X(cell)*B, const T alpha);\
void Y(spcellmulmat_thread)(X(cell) **C, const Y(spcell)*A, const X(cell)*B, const T alpha, const int nthread);\
void Y(sptcellmulmat_thread)(X(cell) **C, const Y(spcell)*A, const X(cell)*B, const T alpha, const int nthread);\
void Y(spcellmulmat_each)(X(cell) **xout, Y(spcell) *A, X(cell) *xin, T alpha, int trans, int nthread); \
void Y(spfull)(X(mat) **out0, const X(sp) *A, const T f);\
void Y(sptfull)(X(mat) **out0, const X(sp) *A, const T f);\
void Y(spcellfull)(X(cell) **out0, const Y(spcell) *A, const T f);\
void Y(sptcellfull)(X(cell) **out0, const Y(spcell) *A, const T f);\
X(sp) *Y(spadd2)(X(sp) *A,X(sp)*B,T a,T b) CHECK_UNUSED_RESULT;\
void Y(spadd)(X(sp) **A0, const X(sp) *B);\
void Y(spcelladd)(Y(spcell) **A0, const Y(spcell) *B);\
void Y(spaddI)(X(sp) **A0, double alpha);\
void Y(spcelladdI)(Y(spcell) *A0, double alpha);\
X(sp) *Y(sptrans)(const X(sp) *A) CHECK_UNUSED_RESULT;\
X(sp) *Y(spmulsp)(const X(sp) *A, const X(sp) *B) CHECK_UNUSED_RESULT;\
X(sp) *Y(sptmulsp)(const X(sp) *A, const X(sp) *B) CHECK_UNUSED_RESULT;\
void Y(spmulsp2)(X(sp) **C0, const X(sp) *A, const X(sp) *B, const T scale);\
Y(spcell) *Y(spcellmulspcell)(const Y(spcell) *A, const Y(spcell) *B, const T scale) CHECK_UNUSED_RESULT;\
Y(spcell) *Y(spcellnew)(const long nx, const long ny) CHECK_UNUSED_RESULT;\
Y(spcell) *Y(spcelltrans)(const Y(spcell) *spc) CHECK_UNUSED_RESULT;\
void Y(spcellfree_do)(Y(spcell) *spc);\
X(sp) *Y(spcat)(const X(sp) *A, const X(sp) *B, int type) CHECK_UNUSED_RESULT;\
X(sp) *Y(spcell2sp)(const Y(spcell) *A) CHECK_UNUSED_RESULT;\
X(mat) *Y(spsum)(const X(sp) *A, int col) CHECK_UNUSED_RESULT;\
X(mat) *Y(spsumabs)(const X(sp) *A, int col) CHECK_UNUSED_RESULT;\
void Y(spclean)(X(sp) *A);\
void Y(spcellmulvec)(T *restrict yc, const Y(spcell) *Ac, const T * restrict xc, T alpha);\
void Y(spdropeps)(X(sp) *A);\
void Y(spcelldropeps)(Y(spcell) *A);\
void Y(spsort)(X(sp) *A);\
void Y(spcellsort)(Y(spcell) *A);\
void Y(spsym)(X(sp) *A);\
void Y(spcellsym)(Y(spcell) *A);\
X(sp) *Y(spconvolvop)(X(mat) *A) CHECK_UNUSED_RESULT;\
X(sp) *Y(spperm)(X(sp) *A, int reverse, long *pcol, long *prow) CHECK_UNUSED_RESULT;\
X(sp) *Y(spinvbdiag)(const X(sp) *A, long bs) CHECK_UNUSED_RESULT;\
X(cell) *Y(spblockextract)(const X(sp) *A, long bs) CHECK_UNUSED_RESULT;

#endif
