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

#ifndef AOS_MATBIN_H
#define AOS_MATBIN_H
#include "bin.h"
#include "type.h"
#define AOS_MATBIN_DEF(X,Y,T)\
void X(writedata)(file_t *fp, const X(mat) *A);	\
void X(cellwritedata)(file_t *fp, const X(cell) *dc);\
X(mat) *X(readdata)(file_t *fp, uint32_t magic);	\
X(cell)* X(cellreaddata)(file_t *fp, uint32_t magic);	\
void Y(spwritedata)(file_t *fp, const X(sp) *sp);\
X(sp) *Y(spreaddata)(file_t *fp, uint32_t magic);		     \
void X(write)(const X(mat) *A, const char *format,...) CHECK_ARG(2); \
void X(cellwrite)(const X(cell) *dc, const char *format,...) CHECK_ARG(2); \
X(mat)* X(read)(const char *format,...) CHECK_ARG(1);\
X(cell)* X(cellread)(const char *format,...) CHECK_ARG(1);\
void Y(spwrite)(const X(sp) *sp, const char *format,...) CHECK_ARG(2); \
void Y(spcellwrite)(const Y(spcell) *spc, const char *format,...) CHECK_ARG(2);\
X(sp)* Y(spread)(const char *format,...) CHECK_ARG(1);\
Y(spcell) *Y(spcellread)(const char *format,...) CHECK_ARG(1);\
X(mat) *X(new_mmap)(long nx, long ny, const char *format,...) CHECK_ARG(3);\
X(cell)* X(cellnew_mmap)(long nx, long ny, long *nnx,long *nny,const char *format,...) CHECK_ARG(5);\
X(cell)* X(cellnewsame_mmap)(long nx,long ny,long mx,long my,const char *format, ...) CHECK_ARG(5);\
void X(cellswrite)(X(cell) *A, double scale, const char *format, ...) CHECK_ARG(3); \
void X(swrite)(X(mat) *A, double scale, const char *format, ...) CHECK_ARG(3);
#endif
