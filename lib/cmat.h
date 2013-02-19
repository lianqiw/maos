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

#ifndef AOS_LIB_CMAT_H
#define AOS_LIB_CMAT_H
/**
   \file cmat.h Contains the mathematically functions regarding to cmat and ccell object
*/
#include "mat.h"
#include "cell.h"
#include "matbin.h"
AOS_MAT_DEF(AOS_CMAT,AOS_DMAT,AOS_CSP,dcomplex,double)
AOS_CELL_DEF(AOS_CMAT,AOS_CSP,dcomplex)
AOS_MATBIN_DEF(AOS_CMAT,AOS_CSP,dcomplex)
/*The following are only useful for cmat */
double cmaxabs(const cmat *A);
double cminabs(const cmat *A);
double csumabs(const cmat *A);

void ccwm3(cmat *restrict A, const cmat *restrict B, const cmat *restrict C);
void ccwm2(cmat *restrict A, const cmat *restrict B, const double alpha);
void ccwmcol(cmat *restrict A, const cmat *restrict B);
void ccwm3col(cmat *restrict A,const cmat *restrict W,const cmat *restrict B);
void ccwmrow(cmat *restrict A, const cmat *restrict B);
void ccwmcol2(cmat *restrict A, 
	     const dcomplex *restrict B1, const double wt1,
	     const dcomplex *restrict B2, const double wt2);
void ccwmrow2(cmat *restrict A, 
	     const dcomplex *restrict B1, const double wt1,
	     const dcomplex *restrict B2, const double wt2);
void ccwmc(cmat *restrict A, const cmat *restrict B, const double alpha);
void ccwmd(cmat *restrict A, const dmat *restrict B, const double alpha);
void cembed_wvf(cmat *restrict A, const double *opd, const double *amp,
	       const int nopdx, const int nopdy, 
		const double wvl, const double theta);
void cembed(cmat *restrict A, const cmat *restrict B, const double theta,CEMBED flag);
void cembedd(cmat *restrict A, dmat *restrict B, const double theta);
void cembedscaleout(cmat *restrict A, const cmat * in, 
		    double xoutscale,double youtscale,
		    const double theta, CEMBED flag);
void ccpcorner(cmat *A, const cmat *restrict B, CEMBED flag);
void cabstoreal(cmat *A);
void cabs2toreal(cmat *A);
void ccpd(cmat**restrict A, const dmat *restrict B);
void creal2d(dmat**restrict A0, double alpha,const cmat *restrict B, double beta);
void cabs22d(dmat**restrict A0, double alpha,const cmat *restrict B, double beta);
void ccp(cmat**restrict A0, const cmat *restrict B);
void cscale(cmat *A, dcomplex alpha);
void ctilt2(cmat *otf, cmat *otfin, double sx, double sy, int pinct);
void ctilt(cmat *otf, double sx, double sy, int pinct);
void ccogreal(double *grad,const cmat *i0,double offsetx,
	      double offsety,double thres, double bkgrnd);
void ccogabs(double *grad,const cmat *i0,double offsetx,
	      double offsety,double thres, double bkgrnd);
void cinv_inplace(cmat*A);
void cinvspd_inplace(cmat *A);
void cmulvec(dcomplex *restrict y, const cmat * restrict A,
	     const dcomplex *restrict x, const dcomplex alpha);

#endif
