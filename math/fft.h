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


#ifndef AOS_LIB_FFT_H
#define AOS_LIB_FFT_H
#ifndef AOS_LIB_TYPE
#error "Don't include this file directly"
#endif
#include "type.h"
/**
   Routines to do FFT on X(mat) or dcell.
*/
typedef struct fft_t fft_t;
#define AOS_FFT_DEF(X, XR)						\
    void X(fft_free_plan)(struct fft_t *fft);			\
    void X(fft2)(X(mat) *A, int dir);				\
    void X(fft2i)(X(mat) *A, int dir);				\
    void X(fft2s)(X(mat) *A, int dir);				\
    void X(fft2partial)(X(mat) *A, int ncomp, int dir);		\
    X(mat) *X(ffttreat)(X(mat) *A);				\
    void XR(fft_free_plan)(struct fft_t *fft);			\
    void XR(cell_fft2)(XR(cell) *dc, int dir);			\
    void XR(fft1plan_r2hc)(XR(mat) *out, int dir);		\
    void XR(fft1)(XR(mat) *A, int dir);
#endif