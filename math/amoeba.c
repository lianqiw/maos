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
#include <stdlib.h>
#include "../sys/sys.h"
#include "defs.h"
#include "mathdef.h"
#define NRANSI
#define NR_END 1
#define FREE_ARG char*
#ifdef USE_COMPLEX
#error "Do not use for complex numbers"
#endif
/**
   2-D mininum search using downhill simplex method.
 */

/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

/**
   Extrapolates by a factor fac through the face of the simplex across from the
   high point, tries it, and replaces the high point if the new point is better.
*/
static T amotry(T **p, T y[], T psum[], int ndim,
		     T (*funk)(T [], void *data), void *data, int ihi, T fac){
    int j;
    T fac1,fac2,ytry;
    T ptry[ndim];
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=0;j<ndim;j++) {
	ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    }
    ytry=(*funk)(ptry,data);
    if (ytry < y[ihi]) {
	y[ihi]=ytry;
	for (j=0;j<ndim;j++) {
	    psum[j] += ptry[j]-p[ihi][j];
	    p[ihi][j]=ptry[j];
	}
    }
    return ytry;
}

#define NMAX 5000
#define GET_PSUM					\
    for (j=0;j<ndim;j++) {				\
	for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];	\
	psum[j]=sum;					\
    }
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

/**
   Multidimensional minimization of the function funk(x) where x[1..ndim] is a
   vector in ndim dimensions, by the downhill simplex method of Nelder and
   Mead. 

   The matrix p[1..ndim+1] [1..ndim] is input. Its ndim+1 rows are
   ndim-dimensional vectors which are the vertices of the starting simplex. Also
   input is the vector y[1..ndim+1], whose components must be preinitialized to
   the values of funk evaluated at the ndim+1 vertices (rows) of p; and ftol the
   fractional convergence tolerance to be achieved in the function value
   (n.b.!). 

   On output, p and y will have been reset to ndim+1 new points all within ftol
   of a minimum function value, and nfunk gives the number of function
   evaluations taken.
 */
static void amoeba(T **p, T y[], int ndim, T ftol,
		   T (*funk)(T [], void *data), void *data, int *nfunk){
    int i,ihi,ilo,inhi,j,mpts=ndim+1;
    T rtol,sum,swap,ysave,ytry;
    T psum[ndim];
    *nfunk=0;
    GET_PSUM;
    for (;;) {
	ilo=0;
	ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
	for (i=0;i<mpts;i++) {
	    if (y[i] <= y[ilo]) {
		ilo=i;
	    }
	    if (y[i] > y[ihi]) {
		inhi=ihi;
		ihi=i;
	    } else if (y[i] > y[inhi] && i != ihi) {
		inhi=i;
	    }
	}
	/* We use ftol as absolute instead of relative tolerance. relative tol does not work near zero.*/
	/*rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));*/
	rtol=2.0*fabs(y[ihi]-y[ilo]);
	if (rtol < ftol) {//solution found at ilo.
	    SWAP(y[0],y[ilo]);
	    for (i=0;i<ndim;i++) {
		SWAP(p[0][i],p[ilo][i]);
	    }
	    break;
	}
	if (*nfunk >= NMAX) {
	    warning("NMAX exceeded\n");
	    SWAP(y[0],y[ilo]);
	    for (i=0;i<ndim;i++) {
		SWAP(p[0][i],p[ilo][i]);
	    }
	    break;
	}
	*nfunk += 2;
	ytry=amotry(p,y,psum,ndim,funk,data,ihi,-1.0);
	if (ytry <= y[ilo])
	    ytry=amotry(p,y,psum,ndim,funk,data,ihi,2.0);
	else if (ytry >= y[inhi]) {
	    ysave=y[ihi];
	    ytry=amotry(p,y,psum,ndim,funk,data,ihi,0.5);
	    if (ytry >= ysave) {
		for (i=0;i<mpts;i++) {
		    if (i != ilo) {
			for (j=0;j<ndim;j++)
			    p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
			y[i]=(*funk)(psum,data);
		    }
		}
		*nfunk += ndim;
		GET_PSUM;
	    }
	} else --(*nfunk);
    }
}
#undef SWAP
#undef GET_PSUM
#undef NMAX
#undef NRANSI
/**
  Search minimum along multiple dimenstions. scale is the contains the scale of each dimension. x
  contains initial warm restart values.
*/
int X(minsearch)(T *x, T *scale, int nmod, T ftol, X(minsearch_fun) fun, void *info){
    T pinit[nmod+1][nmod];
    T *pinit2[nmod+1];
    T yinit[nmod+1];
    for(int i=0; i<nmod+1; i++){
	pinit2[i]=pinit[i];
	for(int j=0; j<nmod; j++){
	    pinit[i][j]=x[j];
	}
	if(i>0){
	    pinit[i][i-1]+=scale[i-1];
	}
	yinit[i]=fun(pinit[i], info);
    }
    int ncall=0;
    amoeba(pinit2, yinit, nmod, ftol, fun, info, &ncall);
    for(int j=0; j<nmod; j++){
	x[j]=pinit[0][j];
    }
    return ncall;
}
