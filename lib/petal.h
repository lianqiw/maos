/*
  Copyright 2009-2026 Lianqi Wang
  
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
#ifndef AOS_LIB_PETAL_H
#define AOS_LIB_PETAL_H
#include "../math/mathdef.h"
dsp *petal_mkh(long nx, long ny, real cx, real cy, long npetal, real theta);
dsp *petal_mkh_loc(loc_t *loc, long npetal, real theta);
void petal_opd(anydmat opd, real cx, real cy, long npetal, real theta, const dmat *mode);
dmat *petal_mkopd(long nx, long ny, real cx, real cy, long npetal, real theta, const dmat *mode);
typedef struct petal_t petal_t;
void petal_free(petal_t *p);
void petal_free_arr(petal_t **p, int np);
petal_t *petal_setup(const loc_t *saloc, real dx, const dmat *amp, real pdtheta, real pixblur, real theta, int npsf, int widthtt);
void petal_save(petal_t *petal, const char *format, ...) CHECK_ARG(2);
void petal_solve(dcell **phi1, dmat **mphi1, const petal_t *petal, const dcell *ints, const dmat *phi1b, int nrep);
void petal_solve_wfs(dcell **phi1, dmat **mphi1, const dcell *ints, const dmat *phi1b, const loc_t *saloc, const dmat *amp, real pdtheta, real pixblur, real theta, int nrep, int withtt);
#endif
