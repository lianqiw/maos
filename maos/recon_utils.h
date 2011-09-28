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

#ifndef AOS_RECON_UTILS_H
#define AOS_RECON_UTILS_H
#include "maos.h"

void apply_L2(dcell **xout, const spcell *L2, const dcell *xin, double alpha, int nthread);
void apply_invpsd(dcell **xout, const void *extra, const dcell *xin, double alpha, int xb, int yb);
void apply_fractal(dcell **xout, const void *extra, const dcell *xin, double alpha, int xb, int yb);
void TTFR(dcell* x, const dcell *TTF, const dcell *PTTF);
void applyW(dcell *xin, const dsp *W0, const dmat *W1, const double *wt);
dcell* calcWmcc(const dcell *A, const dcell *B, const dsp *W0, 
		const dmat *W1, const dmat *wt);
void focus_tracking(SIM_T*simu);

void TomoR(dcell **xout, const void *A, 
	   const dcell *xin, const double alpha);
void TomoL(dcell **xout, const void *A, 
	   const dcell *xin, const double alpha);

void FitL(dcell **xout, const void *A, 
	  const dcell *xin, const double alpha);
void FitR(dcell **xout, const void *A, 
	  const dcell *xin, const double alpha);

dsp *nea2sp(dmat **nea, long nsa);
void psfr_calc(SIM_T *simu, dcell *opdr, dcell *dmpsol, dcell *dmerr_hi, dcell *dmerr_lo);
void shift_grad(SIM_T *simu);
imat* act_coord2ind(loc_t *aloc, const char *fndead);
#endif
