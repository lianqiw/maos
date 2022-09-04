/*
  Copyright 2009-2022 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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


#ifndef AOS_LIB_MATBIN_H
#define AOS_LIB_MATBIN_H
#ifndef AOS_LIB_TYPE
#error "Don't include this file directly"
#endif

#define AOS_MATBIN_DEF(X,T)						\
    void X(writedata)(file_t *fp, const X(mat) *A, long ncol);			\
    X(mat) *X(readdata)(file_t *fp, header_t *header);			\
    X(mat) *X(new_mmap)(long nx, long ny, const char *header, const char *format,...) CHECK_ARG(4); \
    X(cell)* X(cellnew_mmap)(long nx,long ny,long *nnx,long *nny, const char *header, const char *format,...) CHECK_ARG(6); \
    X(cell)* X(cellnewsame_mmap)(long nx,long ny,long mx,long my, const char *header, const char *format,...) CHECK_ARG(6); \
    X(mat)* X(read_mmap)(const char *format, ...) CHECK_ARG(1);			\
    X(cell)* X(cellread_mmap)(const char *format, ...) CHECK_ARG(1);

#endif
