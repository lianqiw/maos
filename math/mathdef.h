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

#ifndef AOS_LIB_MATH_H
#define AOS_LIB_MATH_H
#include "../sys/sys.h"
#include "type.h"
#include "mat.h"
#include "matbin.h"
#include "cell.h"
#include "sp.h"
#include "spbin.h"
#include "fft.h"
#include "chol.h"

#include "imat.h"
#include "cellarr.h"
#include "mathmisc.h"
//#ifdef USE_SINGLE
//Single
AOS_MAT_DEF(AOS_SMAT,AOS_SMAT,AOS_SSP,float,float)
AOS_MAT_DEF(AOS_ZMAT,AOS_SMAT,AOS_ZSP,fcomplex,float)
AOS_CMAT_DEF(AOS_ZMAT,AOS_SMAT,AOS_ZSP,fcomplex,float)

AOS_MATBIN_DEF(AOS_SMAT,AOS_SSP,float)
AOS_MATBIN_DEF(AOS_ZMAT,AOS_ZSP,fcomplex)

AOS_CELL_DEF(AOS_SMAT,AOS_SSP,float,float)
AOS_CELL_DEF(AOS_ZMAT,AOS_ZSP,fcomplex,float)

AOS_SP_DEF(AOS_SMAT,AOS_SSP,float,float,fcomplex)
AOS_SP_DEF(AOS_ZMAT,AOS_ZSP,fcomplex,float,fcomplex)

AOS_SPBIN_DEF(AOS_SMAT,AOS_SSP,float)
AOS_SPBIN_DEF(AOS_ZMAT,AOS_ZSP,fcomplex)

AOS_FFT_DEF(AOS_SMAT)
AOS_FFT_DEF(AOS_ZMAT)
//#else
//Double
AOS_MAT_DEF(AOS_DMAT,AOS_DMAT,AOS_DSP,double,double)
AOS_MAT_DEF(AOS_CMAT,AOS_DMAT,AOS_CSP,dcomplex,double)
AOS_MATBIN_DEF(AOS_DMAT,AOS_DSP,double)
AOS_MATBIN_DEF(AOS_CMAT,AOS_CSP,dcomplex)
AOS_CMAT_DEF(AOS_CMAT,AOS_DMAT,AOS_CSP,dcomplex,double)

AOS_CELL_DEF(AOS_DMAT,AOS_DSP,double,double)
AOS_CELL_DEF(AOS_CMAT,AOS_CSP,dcomplex,double)

AOS_SP_DEF(AOS_DMAT,AOS_DSP,double,double,dcomplex)
AOS_SP_DEF(AOS_CMAT,AOS_CSP,dcomplex,double,dcomplex)

AOS_SPBIN_DEF(AOS_DMAT,AOS_DSP,double)
AOS_SPBIN_DEF(AOS_CMAT,AOS_CSP,dcomplex)

AOS_FFT_DEF(AOS_DMAT)
AOS_FFT_DEF(AOS_CMAT)

//#endif
#endif
