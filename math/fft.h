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
#ifndef AOS_LIB_FFT_H
#define AOS_LIB_FFT_H
#ifndef AOS_LIB_MATH_H
#error "Don't include this file directly"
#endif
#include "type.h"
/**
   \file fft.h
   Routines to do FFT on X(mat) or dcell.
*/
typedef struct fft_t fft_t;
#define AOS_FFT_DEF(X)						\
    void X(fft_free_plan)(struct fft_t *fft);			\
    void X(fft2plan)(X(mat) *A, int dir);			\
    void X(fft2)(X(mat) *A, int dir);				\
    void X(fft2i)(X(mat) *A, int dir);				\
    void X(fft2s)(X(mat) *A, int dir);				\
    void X(fft2partialplan)(X(mat) *A, int ncomp, int dir);	\
    void X(fft2partial)(X(mat) *A, int ncomp, int dir);		\
    X(mat) *X(ffttreat)(X(mat) *A);				\
    void X(cell_fft2plan)(X(cell) *dc, int dir, int nthreads);	\
    void X(cell_fft2)(X(cell) *dc, int dir);			\
    void X(fft1plan_r2hc)(X(mat) *out, int dir);		\
    void X(fft1)(X(mat) *A, int dir);
#endif
