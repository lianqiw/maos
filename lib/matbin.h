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

#ifndef AOS_MATBIN_H
#define AOS_MATBIN_H
#include "type.h"

#define AOS_MATBIN_DEF(X,Y,T)\
void X(writedata)(file_t *fp, const X(mat) *A);	\
void X(cellwritedata)(file_t *fp, const X(cell) *dc);\
X(mat) *X(readdata)(file_t *fp, header_t *header); \
X(cell)* X(cellreaddata)(file_t *fp, header_t *header); \
void X(write)(const X(mat) *A, const char *format,...) CHECK_ARG(2); \
void X(cellwrite)(const X(cell) *dc, const char *format,...) CHECK_ARG(2); \
X(mat)* X(read)(const char *format,...) CHECK_ARG(1);\
X(cell)* X(cellread)(const char *format,...) CHECK_ARG(1);\
X(cell)**X(cellreadarr)(long *nxout, long *nyout, const char *format,...) CHECK_ARG(3); \
void X(cellwritearr)(X(cell)**A,long nxin, long nyin, const char *format,...) CHECK_ARG(4); \
void X(cellswrite)(X(cell) *A, double scale, const char *format, ...) CHECK_ARG(3); \
void X(swrite)(X(mat) *A, double scale, const char *format, ...) CHECK_ARG(3);\
X(mat) *X(new_mmap)(long nx, long ny, char *header, const char *format,...) CHECK_ARG(4); \
X(cell)* X(cellnew_mmap)(long nx,long ny,long *nnx,long *nny, char *header1,char**header2,const char *format,...) CHECK_ARG(7); \
X(cell)* X(cellnewsame_mmap)(long nx,long ny,long mx,long my, char *header, const char *format,...) CHECK_ARG(6); \
X(mat*) X(read_mmap)(const char *format, ...);\
X(cell*) X(cellread_mmap)(const char *format, ...);

#endif
