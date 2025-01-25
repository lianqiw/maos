/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "common.h"

void apply_L2(dcell **xout, const dspcell *L2, const dcell *xin, real alpha);
void apply_invpsd(dcell **xout, const void *extra, const dcell *xin, real alpha, int xb, int yb);
void apply_fractal(dcell **xout, const void *extra, const dcell *xin, real alpha, int xb, int yb);
void remove_mode(dcell *x, const dcell *M, const dcell *PM);
void applyW(dcell *xin, const dsp *W0, const dmat *W1, const real *wt);
dcell* calcWmcc(const dcell* A, const dcell* B, const dsp* W0, const dmat* W1, const dmat* wt);

dsp *nea2sp(dmat *nea, int ll, int ur, int ng);
void nea_chol(dmat **pout, const dmat *in, int ng);
void nea_inv(dmat **pout, const dmat *in, int ng, real scale);
void nea_mm(dmat **pout, const dmat *in, int ng);
void nea_check(dmat *nea, int nsa, int ng);

void psfr_calc(sim_t *simu, dcell *opdr, dcell *dmpsol, dcell *dmerr, dcell *dmerr_lo);
cn2est_t* cn2est_prepare(const parms_t *parms, const powfs_t *powfs);
void cn2est_isim(dcell *cn2res, recon_t *recon, const parms_t *parms, const dcell *grad, int *tomo_update);

#endif
