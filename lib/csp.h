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

#ifndef AOS_LIB_CSP_H
#define AOS_LIB_CSP_H
/**
   \file csp.h
   Contains functions for complex sparse csp
*/
#include "sp.h"
#include "spbin.h"
AOS_SP_DEF(AOS_CMAT,AOS_CSP,dcomplex)
AOS_SPBIN_DEF(AOS_CMAT,AOS_CSP,dcomplex)
#endif
