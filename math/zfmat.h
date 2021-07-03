/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_LIB_ZFMAT_H
#define AOS_LIB_ZFMAT_H
#include "type.h"
/**
   \file zfmat.h
   zfmatr is an object used to write dmat into file incrementally.
*/

typedef struct zfmat zfmat;
zfmat* zfmat_init(long nx, long ny, const char*format,...) CHECK_ARG(3);
void zfmat_push(zfmat *ca, long count, real *p);
void zfmat_close(zfmat *ca);
void zfmat_close_n(zfmat **ca, int nc);
#endif
