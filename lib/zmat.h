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
#ifndef AOS_LIB_ZMAT_H
#define AOS_LIB_ZMAT_H
/**
   \file dmat.h Contains the mathematically functions regarding to dmat and dcell object
*/
#include "mat.h"
#include "cell.h"
#include "matbin.h"
AOS_MAT_DEF   (AOS_ZMAT,AOS_SMAT,AOS_SSP,fcomplex,float)
AOS_CELL_DEF  (AOS_ZMAT,AOS_ZSP,fcomplex)
AOS_MATBIN_DEF(AOS_ZMAT,AOS_ZSP,fcomplex)

/*The following are only for dmat. */
/*dmat *denc(dmat *A, dmat *dvec, int type, int nthread); */
#endif
