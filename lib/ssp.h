/*
  Copyright 2009-2012 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_LIB_SSP_H
#define AOS_LIB_SSP_H
/**
   \file dsp.h
   Contains functions for double sparse dsp
*/
#include "sp.h"
#include "spbin.h"
AOS_SP_DEF   (AOS_SMAT,AOS_SSP,float)
AOS_SPBIN_DEF(AOS_SMAT,AOS_SSP,float)
#endif
