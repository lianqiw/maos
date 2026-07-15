/*
  Copyright 2009-2026 Lianqi Wang
  
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
   \file lib/libmisc.h
*/
#ifndef AOS_LIB_MISC_H
#define AOS_LIB_MISC_H
#include "../math/mathdef.h"
#include "accphi.h"
extern const real RAD2AS;
extern const real RAD2MAS;
extern const real MAS2RAD;
extern const real AS2RAD;

void addnoise(dmat *A, rand_t* rstat, 
	      const real bkgrnd, const real bkgrndc, 
	      const dmat *bkgrnd2, const dmat *bkgrnd2c,
	      const dmat* qe, real rne, real excess);
void addnoise_grad(dmat *grad, const dmat *neal, rand_t *srand);
int cog_multi(dmat** cg, const dmat* im,	const dmat* loc,int wsize);
void cog_nea(real* nea, const dmat* ints, const dmat* cogmask, real cogthres, real cogoff, int ntry,
  rand_t* rstat, real bkgrnd, real bkgrndc, const dmat* bkgrnd2i, const dmat* bkgrnd2ic, real rne
);
dmat *poly2fit(const dmat *in, const dmat *out, int maxorder);
dmat *loc_calib(const dsp *GA, propdata_t *propdata, int maxorder);
dmat *polyfit(const dmat *x, const dmat *y, int maxorder);
dmat *polyval(const dmat *x, const dmat *coeff, int separate);
int wrap_seq(long index, long n);
real wrap2range(real val, real low, real high);
real pchip_wt(int i, float u);
int atm_interp(real *wt, int ips, int isim, int dtrat, int natm, int interp);
dccell *wfe_fov_prep(const dmat *thetax, const dmat *thetay, real unit, int npath, int nwvl);
void wfe_fov_fill(dccell *res, const dcell *clep, int ipath, int isim);
#endif
