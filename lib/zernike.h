/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file zernike.h
*/
#ifndef AOS_LIB_ZERNIKE_H
#define AOS_LIB_ZERNIKE_H
#include "../math/mathdef.h"
dmat *zernike_Rnm(const dmat *locr, int ir, int im);
dmat* zernike(const loc_t *loc, real D, int rstart, int rend, int flag);
dmat *zernike_cov_kolmogorov(int nr);
dmat *cov_vonkarman(const loc_t *loc, const dmat *modz, real L0);
dmat *cov_diagnolize(const dmat *mz, const dmat *cov);
dmat *KL_vonkarman(const loc_t *loc, int nmod, real L0);
dmat *fft_mode(const loc_t *loc, real D, real px, real py);
#endif
