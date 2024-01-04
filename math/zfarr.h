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



#ifndef AOS_LIB_ZFARR_H
#define AOS_LIB_ZFARR_H
#include "type.h"
/**
   \file zfarr.h
   zfarr is an object used to write arrays of dcell or ccell into file.
   Mainly used to output PSF into files.
*/
/*
   used to save array of dmat, cmat, ccell or dcell. mainly used to save
psfout. No user modifiable entries.  */
typedef struct zfarr zfarr;
zfarr* zfarr_init2(long nx, long ny, const char *keywords, const char* format,...) CHECK_ARG(4);
zfarr *zfarr_init(long nx, long ny, const char *format, ...) CHECK_ARG(3);
void zfarr_push(zfarr *ca, int i, const_anyarray A);
void zfarr_close(zfarr *ca);
void zfarr_close_n(zfarr **ca, int nc);
#endif
