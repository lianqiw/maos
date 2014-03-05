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
#ifndef AOS_LIB_IMAT_H
#define AOS_LIB_IMAT_H
/**
   \file imat.h
   Routines for imat.
*/
typedef struct imat{
    long nx;
    long ny;
    long *p;
}imat;

typedef struct icell{
    long nx;
    long ny;
    imat **p;
}icell;

imat* inew(long nx, long ny);
icell* icellnew(long nx, long ny);
void ifree(imat *A);
void icellfree(icell *A);
void iwrite(const imat *A, const char *format, ...);
void icellwrite(const icell *A, const char *format, ...);
#endif
