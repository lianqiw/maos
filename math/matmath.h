/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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



#ifndef AOS_LIB_MATMATH_H
#define AOS_LIB_MATMATH_H
#ifndef AOS_LIB_TYPE
#error "Don't include this file directly"
#endif
#define AOS_MATMATH_DEF(X,XR,T,R)					\
    int X(isnan)(const X(mat)*A);					\
    R X(norm)(const X(mat) *in) CHECK_UNUSED_RESULT;			\
    R X(std)(const X(mat) *in) CHECK_UNUSED_RESULT;			\
    void X(randu)(X(mat) *A, const T mean, rand_t *rstat);		\
    void X(randn)(X(mat) *A, const T sigma, rand_t *rstat);		\
    void X(scale)(X(mat) *A, R w);					\
    T X(inn)(const X(mat)*A, const X(mat) *B);				\
    T X(wdot)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT; \
    void X(cwm)(X(mat) *A, const X(mat) *B);				\
    void X(cwm2)(X(mat) *A, const X(mat) *B1, R wt1, const X(mat)*B2, R wt2);\
    void X(cwm3)(X(mat) *restrict A, const X(mat) *restrict W , const X(mat) *restrict B1, const R wt1, const X(mat) *restrict B2, const R wt2); \
    void X(cwmcol)(X(mat) *restrict A, const X(mat) *restrict B);	\
    void X(cwm3col)(X(mat) *restrict A,const X(mat) *restrict W, const X(mat) *restrict B1, const R wt1, const X(mat) *restrict B2, const R wt2); \
    void X(cwmrow)(X(mat) *restrict A, const X(mat) *restrict B);	\
    void X(cwmcol2)(X(mat) *restrict A,	const X(mat) *restrict B1, const R wt1, const X(mat) *restrict B2, const R wt2);		\
    void X(cwmrow2)(X(mat) *restrict A,	const X(mat) *restrict B1, const R wt1, const X(mat) *restrict B2, const R wt2);		\
    void X(cwdiv)(X(mat) *B, const X(mat) *A, T value);			\
    void X(mulvec)(T *restrict y, const X(mat) * restrict A, const T *restrict x, const T alpha); \
    void X(mm)(X(mat)**C0, const T beta, const X(mat) *A, const X(mat) *B, const char trans[2], const T alpha); \
    void X(invspd_inplace)(X(mat) *A);					\
    X(mat)* X(invspd)(const X(mat) *A) CHECK_UNUSED_RESULT;		\
    void X(inv_inplace)(X(mat)*A);					\
    X(mat)* X(inv)(const X(mat) *A) CHECK_UNUSED_RESULT;		\
    X(mat) *X(chol)(const X(mat) *A) CHECK_UNUSED_RESULT;		\
    X(mat) *X(mcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT; \
    X(mat) *X(imcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT; \
    X(mat) *X(tmcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT; \
    X(mat) *X(pinv2)(const X(mat) *A, const_anyarray W, R thres) CHECK_UNUSED_RESULT; \
    X(mat) *X(pinv)(const X(mat) *A, const_anyarray W) CHECK_UNUSED_RESULT; \
    T X(diff)(const X(mat) *A, const X(mat) *B) CHECK_UNUSED_RESULT;	\
    int X(circle)(X(mat) *A, R cx, R cy, R dx, R dy, R r, T val); \
    int X(circle_symbolic)(X(mat) *A, R cx, R cy, R dx, R dy, R r);	\
    void X(rotvec)(X(mat) *A, const R theta);				\
    void X(rotvect)(X(mat) *A, const R theta);				\
    void X(rotvecnn)(X(mat) **B0, const X(mat) *A, R theta);		\
    void X(mulvec3)(T *y, const X(mat) *A, const T *x);			\
    void X(corr)(X(mat) **corr, const X(mat) *A, const X(mat) *B);	\
    void X(para3)(R *grad, const X(mat) *corr);				\
    void X(cog)(R *grad,const X(mat) *i0,R offsetx, R offsety, R thres, R bkgrnd, R flux); \
    void X(shift2center)(X(mat) *A, R offsetx, R offsety);		\
    int X(clip)(X(mat) *A, R min, R max);				\
    void X(gramschmidt)(X(mat) *Mod, R *amp);				\
    void X(muldiag)(X(mat) *A, const X(mat) *s);			\
    void X(cwpow)(X(mat) *A, R power);					\
    void X(cwexp)(X(mat) *A, R alpha);					\
    void X(cwpow_thres)(X(mat) *A, R power, R thres);			\
    void X(svd)(X(mat) **U, XR(mat) **Sdiag, X(mat) **VT, const X(mat) *A); \
    void X(svd_cache)(X(mat) **U, XR(mat) **Sdiag, X(mat) **VT, const X(mat) *A); \
    void X(svd_pow)(X(mat) *A, R power, R thres);			\
    void X(expm)(X(mat) **out, R alpha, const X(mat) *A, R beta);	\
    void X(polyval)(X(mat) *A, XR(mat)*p);				\
    void X(addI)(X(mat) *A, T val);					\
    void X(add)(X(mat) **B0, T bc,const X(mat) *A, const T ac);		\
    void X(add_relax)(X(mat) **B0, T bc,const X(mat) *A, const T ac);		\
    void X(addcol)(X(mat) *B, long icol, T bc, const X(mat) *A, const T ac);\
    void X(adds)(X(mat*)A, const T ac);	\
    X(mat)* X(logspace)(R emin, R emax, long n) CHECK_UNUSED_RESULT;	\
    X(mat)* X(linspace)(R min, R dx, long n) CHECK_UNUSED_RESULT;	\
    X(mat)* X(interp1)(const X(mat) *xin, const X(mat) *yin, const X(mat) *xnew, T y0) CHECK_UNUSED_RESULT; \
    X(mat)* X(interp1_2)(const X(mat) *xyin, const X(mat) *xnew) CHECK_UNUSED_RESULT; \
    X(mat)* X(interp1linear)(const X(mat) *xin, const X(mat) *yin, const X(mat) *xnew, T y0) CHECK_UNUSED_RESULT; \
    X(mat)* X(interp1log)(const X(mat) *xin, const X(mat) *yin, const X(mat) *xnew, T y0) CHECK_UNUSED_RESULT; \
    void X(blend)(X(mat) *restrict A, X(mat) *restrict B, int overlap);	\
    void X(histfill)(X(mat) **out, const X(mat)* A, R center, R spacing, int n); \
    X(mat) *X(spline_prep)(X(mat) *x, X(mat) *y);			\
    X(mat)* X(spline_eval)(X(mat) *coeff, X(mat)* x, X(mat)*xnew);	\
    X(mat)* X(spline)(X(mat) *x,X(mat) *y,X(mat) *xnew);		\
    void X(cwlog10)(X(mat) *A);						\
    void X(cwlog)(X(mat) *A);						\
    void X(embed)(X(mat) *restrict A, const X(mat) *restrict B, const R theta); \
    void X(embedd)(X(mat) *restrict A, const XR(mat) *restrict B, const R theta); \
    R X(fwhm)(X(mat) *A);						\
    void X(gauss_fit)(R*mr, R*ma, R*mb, R*theta, X(mat)*A, R thres);\
    R X(gauss_width)(X(mat)*A, R thres);\
    R X(fwhm_gauss)(X(mat) *A);\
    X(mat) *X(enc)(const X(mat) *A, const X(mat) *dvec, int type, int nthread);	typedef T (*X(minsearch_fun))(const T *x, void *info);		\
    int X(minsearch)(T *x, int nmod, T ftol, int nmax, X(minsearch_fun) fun, void *info); \
    void X(bessik)(T x, T xnu, T *ri, T *rk, T *rip, T *rkp);		\
    T X(trapz)(const X(mat)*x, const X(mat)*y);				\
    R X(cellnorm)(const X(cell) *in);					\
    void X(cellscale)(anyarray A, R w);					\
    void X(celldropempty)(X(cell) **A0, int dim);			\
    T X(cellinn)(const X(cell)*A, const X(cell)*B);			\
    void X(cellcwm)(X(cell) *B, const X(cell) *A);			\
    X(cell)* X(cellinvspd)(X(cell) *A);					\
    X(cell)* X(cellinv)(X(cell) *A);					\
    X(cell)* X(cellinvspd_each)(X(cell) *A);				\
    X(cell)* X(cellpinv2)(const X(cell) *A, const_anyarray W, R thres);		\
    X(cell)* X(cellpinv)(const X(cell) *A, const X(spcell) *W);		\
    void X(cellsvd_pow)(X(cell) *A, R power, R thres);		\
    void X(cellcwpow)(X(cell)*A, R power);				\
    void X(celldropzero)(X(cell) *B, R thres);				\
    R X(celldiff)(const X(cell) *A, const X(cell) *B);			\
    int X(cellclip)(X(cell) *Ac, R min, R max);				\
    X(cell) *X(cellsub)(const X(cell) *in, long sx, long nx, long sy, long ny);	\
    X(cell) *X(bspline_prep)(X(mat)*x, X(mat)*y, X(mat) *z);		\
    X(mat) *X(bspline_eval)(X(cell)*coeff, X(mat) *x, X(mat) *y, X(mat) *xnew, X(mat) *ynew);\
	void X(maprot)(anyarray A_, real theta);

/*The following are only useful for cmat */
#define AOS_CMATMATH_DEF(X,XR,T,R)					\
    void X(cwmc)(X(mat) *restrict A, const X(mat) *restrict B, const R alpha); \
    void X(cwmd)(X(mat) *restrict A, const XR(mat) *restrict B, const R alpha); \
    void X(embed_wvf)(X(mat) *restrict A, const R *opd, const R *amp,	const int nopdx, const int nopdy, const R wvl, const R theta);			\
    void X(embedc_flag)(X(mat) *restrict A, const X(mat) *restrict B, const R theta,CEMBED flag); \
    void X(cpcorner)(X(mat) *A, const X(mat) *restrict B, CEMBED flag);	\
    void X(abstoreal)(X(mat) *A);					\
    void X(abs2toreal)(X(mat) *A, R scale);					\
    void X(cpd)(X(mat)**restrict A, const XR(mat) *restrict B);		\
    void X(real2d)(XR(mat)**restrict A0, R alpha,const X(mat) *restrict B, R beta); \
    void X(imag2d)(XR(mat)**restrict A0, R alpha,const X(mat) *restrict B, R beta); \
    void X(abs22d)(XR(mat)**restrict A0, R alpha,const X(mat) *restrict B, R beta); \
    void X(tilt2)(X(mat) *otf, X(mat) *otfin, R sx, R sy, int peak_corner); \
    void X(tilt)(X(mat) *otf, R sx, R sy, int peak_corner);		\
    void X(cogreal)(R *grad,const X(mat) *i0,R offsetx,	R offsety,R thres, R bkgrnd);	\
    void X(cogabs)(R *grad,const X(mat) *i0,R offsetx,  R offsety,R thres, R bkgrnd);			

#endif
