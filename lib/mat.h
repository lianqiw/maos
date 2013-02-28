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

#ifndef AOS_LIB_MAT_H
#define AOS_LIB_MAT_H
#include "random.h"
#include "type.h"
void dembed(dmat *restrict A, dmat *restrict B, const double theta);
#define PDMAT(M,P)   double (*restrict P)[(M)->nx]=(double(*)[(M)->nx])(M)->p
#define PDCELL(M,P)  dmat* (*restrict P)[(M)->nx]=(dmat*(*)[(M)->nx])(M)->p
#define dfree(A)     ({dfree_do((A),0);(A)=NULL;})
#define dcp2(A,B)    memcpy(A->p,B->p,sizeof(double)*A->nx*A->ny)
#define dcellfree(A) ({dcellfree_do(A);A=NULL;})
#define dcellfreearr(A,n) ({for(int in=0; A&&in<n; in++){dcellfree(A[in]);};free(A);A=NULL;})
#define dzero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(double))
#define dhash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(double), key)

#define PSMAT(M,P)   float (*restrict P)[(M)->nx]=(float(*)[(M)->nx])(M)->p
#define PSCELL(M,P)  smat* (*restrict P)[(M)->nx]=(smat*(*)[(M)->nx])(M)->p
#define sfree(A)     ({sfree_do((A),0);(A)=NULL;})
#define scp2(A,B)    memcpy(A->p,B->p,sizeof(float)*A->nx*A->ny)
#define scellfree(A) ({scellfree_do(A);A=NULL;})
#define scellfreearr(A,n) ({for(int in=0; A&&in<n; in++){scellfree(A[in]);};free(A);A=NULL;})
#define szero(A) if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(float))
#define shash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(float), key)

#define PCMAT(M,P)   dcomplex (*restrict P)[(M)->nx]=(dcomplex(*)[(M)->nx])(M)->p
#define PCCELL(M,P)  cmat* (*restrict P)[(M)->nx]=(cmat*(*)[(M)->nx])(M)->p
#define cfree(A)     ({cfree_do(A,0);A=NULL;})
#define ccellfree(A) ({ccellfree_do(A);A=NULL;})
#define ccellfreearr(A,n) ({for(int in=0; A&&in<n; in++){ccellfree(A[in]);};free(A);A=NULL;})
#define cabs2(A)     (pow(creal(A),2)+pow(cimag(A),2))
#define czero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(dcomplex))
#define chash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(dcomplex), key)

#define PZMAT(M,P)   fcomplex (*restrict P)[(M)->nx]=(fcomplex(*)[(M)->nx])(M)->p
#define PZCELL(M,P)  zmat* (*restrict P)[(M)->nx]=(zmat*(*)[(M)->nx])(M)->p
#define zfree(A)     ({zfree_do(A,0);A=NULL;})
#define zcellfree(A) ({zcellfree_do(A);A=NULL;})
#define zabs2(A)     (pow(crealf(A),2)+pow(cimagf(A),2))
#define zzero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(fcomplex))
#define zhash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(fcomplex), key)

#define AOS_MAT_DEF(X,XR,Y,T,R)					\
X(mat) *X(new_ref)(long nx, long ny, T *p) CHECK_UNUSED_RESULT; \
X(mat) *X(new_data)(long nx, long ny, T *p) CHECK_UNUSED_RESULT; \
X(mat) *X(new)(long nx, long ny) CHECK_UNUSED_RESULT;\
void X(init)(X(mat)**A, long nx, long ny) ; \
void X(free_keepdata)(X(mat) *A);\
void X(free_do)(X(mat) *A, int keepdata);\
void X(resize)(X(mat) *A, long nx, long ny);\
X(mat) *X(ref)(X(mat) *in) CHECK_UNUSED_RESULT;\
X(mat) *X(ref_reshape)(X(mat) *in, long nx, long ny) CHECK_UNUSED_RESULT;\
X(mat) *X(refcols)(X(mat) *in, long icol, long ncol) CHECK_UNUSED_RESULT;\
X(mat) *X(sub)(const X(mat) *in, long sx, long nx, long sy, long ny) CHECK_UNUSED_RESULT;\
int X(isnan)(const X(mat)*A);\
X(mat) *X(cat)(const X(mat) *in1, const X(mat) *in2, int dim) CHECK_UNUSED_RESULT;\
void X(arrfree)(X(mat) **As, int n);\
X(mat) *X(dup)(const X(mat) *in) CHECK_UNUSED_RESULT;\
void X(cp)(X(mat) **out0, const X(mat) *in);\
X(mat) *X(trans)(const X(mat) *A) CHECK_UNUSED_RESULT;\
void X(set)(X(mat) *A, const T val);\
R X(max)(const X(mat) *A) CHECK_UNUSED_RESULT;\
R X(min)(const X(mat) *A) CHECK_UNUSED_RESULT;\
R X(norm2)(const X(mat) *in) CHECK_UNUSED_RESULT;\
void X(randu)(X(mat) *A, const T mean, rand_t *rstat);\
void X(randn)(X(mat) *A, const T sigma, rand_t *rstat);\
void X(show)(const X(mat) *A, const char *format,...) CHECK_ARG(2);\
void X(scale)(X(mat) *A, T w);\
T X(sum)(const X(mat) *A) CHECK_UNUSED_RESULT;\
void X(add)(X(mat) **B0, T bc,const X(mat) *A, const T ac);\
void X(adds)(X(mat*)A, const T ac);\
T X(inn)(const X(mat)*A, const X(mat) *B);			\
T X(wdot)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT;\
T X(wdot2)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT;\
T X(wdot3)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT;\
void X(cwm)(X(mat) *B, const X(mat) *A);\
void X(cwdiv)(X(mat) *B, const X(mat) *A, T value);			\
void X(mulvec)(T *restrict y, const X(mat) * restrict A, const T *restrict x, const T alpha);\
void X(mm)(X(mat)**C0, const X(mat) *A, const X(mat) *B, const char trans[2], const T alpha);\
void X(invspd_inplace)(X(mat) *A);\
X(mat)* X(invspd)(const X(mat) *A) CHECK_UNUSED_RESULT;\
void X(inv_inplace)(X(mat)*A);\
X(mat)* X(inv)(const X(mat) *A) CHECK_UNUSED_RESULT;\
X(mat) *X(chol)(const X(mat) *A) CHECK_UNUSED_RESULT; \
X(mat) *X(mcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT;\
X(mat) *X(imcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT;\
X(mat) *X(tmcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT;\
X(mat) *X(pinv)(const X(mat) *A, const X(mat) *wt, const X(sp) *Wsp) CHECK_UNUSED_RESULT;\
T X(diff)(const X(mat) *A, const X(mat) *B) CHECK_UNUSED_RESULT;\
void X(circle)(X(mat) *A, double cx, double cy, double r, T val);\
void X(circle_mul)(X(mat) *A, double cx, double cy, double r, T val);\
void X(circle_symbolic)(X(mat) *A, double cx, double cy, double r);\
void X(fftshift)(X(mat) *A);\
void X(cpcorner2center)(X(mat) *A, const X(mat)*B);\
void X(shift)(X(mat) **B0, const X(mat) *A, int sx, int sy);\
void X(rotvec)(X(mat) *A, const double theta);\
void X(rotvect)(X(mat) *A, const double theta);\
void X(rotvecnn)(X(mat) **B0, const X(mat) *A, double theta);\
void X(mulvec3)(T *y, const X(mat) *A, const T *x);\
void X(cog)(double *grad,const X(mat) *i0,double offsetx, double offsety, double thres, double bkgrnd);\
void X(shift2center)(X(mat) *A, double offsetx, double offsety);\
int X(clip)(X(mat) *A, double min, double max);\
void X(gramschmidt)(X(mat) *Mod, R *amp);	\
void X(muldiag)(X(mat) *A, X(mat) *s);\
void X(cwpow)(X(mat) *A, double power);\
void X(svd)(X(mat) **U, XR(mat) **Sdiag, X(mat) **VT, const X(mat) *A); \
void X(evd)(X(mat) **U, XR(mat) **Sdiag, const X(mat) *A); \
void X(svd_pow)(X(mat) *A, double power, int issym, double thres);\
void X(addI)(X(mat) *A, T val);\
void X(tikcr)(X(mat) *A, T thres);\
void X(mulsp)(X(mat) **yout, const X(mat) *x, const X(sp) *A, const T alpha);\
X(mat)* X(logspace)(double emin, double emax, long n) CHECK_UNUSED_RESULT;\
X(mat)* X(linspace)(double min, double dx, long n) CHECK_UNUSED_RESULT;\
X(mat)* X(interp1)(X(mat) *xin, X(mat) *yin, X(mat) *xnew) CHECK_UNUSED_RESULT;\
X(mat)* X(interp1linear)(X(mat) *xin, X(mat) *yin, X(mat) *xnew) CHECK_UNUSED_RESULT;\
X(mat)* X(interp1log)(X(mat) *xin, X(mat) *yin, X(mat) *xnew) CHECK_UNUSED_RESULT;\
void X(blend)(X(mat) *restrict A, X(mat) *restrict B, int overlap);\
void X(histfill)(X(mat) **out, const X(mat)* A, double center, double spacing, int n);\
X(mat) *X(spline_prep)(X(mat) *x, X(mat) *y);\
X(mat)* X(spline_eval)(X(mat) *coeff, X(mat)* x, X(mat)*xnew);\
X(mat)* X(spline)(X(mat) *x,X(mat) *y,X(mat) *xnew);\
X(cell)* X(bspline_prep)(X(mat)*x, X(mat)*y, X(mat) *z);\
X(mat) *X(bspline_eval)(X(cell)*coeff, X(mat) *x, X(mat) *y, X(mat) *xnew, X(mat) *ynew);\
void X(cwlog10)(X(mat) *A);\
void X(embed_locstat)(X(mat) **out, double alpha, loc_t *loc, R *oin, double beta, int reverse);\
long X(fwhm)(X(mat) *A);\
void X(sort)(X(mat) *A, int ascend);
#endif
