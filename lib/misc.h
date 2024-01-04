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


/**
   \file lib/misc.h
*/
#ifndef AOS_LIB_MISC_H
#define AOS_LIB_MISC_H
#include "../math/mathdef.h"
void addnoise(dmat *A, rand_t* rstat, 
	      const real bkgrnd, const real bkgrndc, 
	      const dmat *bkgrnd2, const dmat *bkgrnd2c,
	      const dmat* qe, real rne, real excess);
void addnoise_grad(dmat *grad, const dmat *neal, rand_t *srand);
dmat *poly2fit(const dmat *in, const dmat *out, int maxorder);
dmat *loc_calib(const dsp *GA, const loc_t *aloc, const loc_t *saloc,
		real dispx, real dispy, real scale, int maxorder);
dmat *polyfit(const dmat *x, const dmat *y, int maxorder);
dmat *polyval(const dmat *x, const dmat *coeff, int separate);
real calc_dither_amp(dmat **res, const dmat *signal, long dtrat, long npoint, int detrend, int combine);
#endif
