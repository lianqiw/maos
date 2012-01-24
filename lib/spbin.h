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

#ifndef AOS_SPBIN_H
#define AOS_SPBIN_H
#include "bin.h"
#include "type.h"

#define AOS_SPBIN_DEF(X,Y,T)\
void Y(spwritedata)(file_t *fp, const X(sp) *sp);\
X(sp) *Y(spreaddata)(file_t *fp, header_t *header);		     \
void Y(spwrite)(const X(sp) *sp, const char *format,...) CHECK_ARG(2); \
void Y(spcellwrite)(const Y(spcell) *spc, const char *format,...) CHECK_ARG(2);\
X(sp)* Y(spread)(const char *format,...) CHECK_ARG(1);\
Y(spcell) *Y(spcellread)(const char *format,...) CHECK_ARG(1);
#endif

