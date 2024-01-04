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
   \file cmath.h
*/
#ifndef AOS_MATH_CMATH_H
#define AOS_MATH_CMATH_H
#ifndef AOS_LIB_TYPE
#define AOS_LIB_TYPE
#include "type.h"
#include "mat.h"
#include "matmath.h"
#include "sp.h"
#include "fft.h"
#include "matbin.h"
#include "spbin.h"
#endif
#define AOS_CMAT(A) c##A
#define AOS_DMAT(A) d##A
//Real, which can be double or float
AOS_MAT_DEF(AOS_CMAT, comp, real);//; to bypass doxygen problem
AOS_MATMATH_DEF(AOS_CMAT, AOS_DMAT, comp, real);
AOS_CMATMATH_DEF(AOS_CMAT, AOS_DMAT, comp, real);
AOS_MATBIN_DEF(AOS_CMAT, comp);
AOS_SP_DEF(AOS_CMAT, comp, real, comp);
AOS_SPBIN_DEF(AOS_CMAT, comp);
AOS_FFT_DEF(AOS_CMAT, AOS_DMAT);
#endif

