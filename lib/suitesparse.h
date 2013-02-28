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

/*Taken from suitesparse package from UFL. */
#ifndef AOS_LIB_SUITE_SPARSE_H
#define AOS_LIB_SUITE_SPARSE_H
#include "type.h"
dsp *cs_multiply (const dsp *A, const dsp *B);
dsp *cs_add (const dsp *A, const dsp *B, double alpha, double beta);
spint cs_dropzeros (dsp *A);
spint cs_droptol (dsp *A, double tol);
dsp *cs_transpose (const dsp *A, spint values);


csp *ccs_multiply (const csp *A, const csp *B);
csp *ccs_add (const csp *A, const csp *B, dcomplex alpha, dcomplex beta);
spint ccs_dropzeros (csp *A);
spint ccs_droptol (csp *A, double tol);
csp *ccs_transpose (const csp *A, spint values);
#endif
