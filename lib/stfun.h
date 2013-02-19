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
#ifndef AOS_LIB_STFUN
#define AOS_LIB_STFUN
typedef struct stfun_t stfun_t;
stfun_t* stfun_init(long nx, long ny, double *amp);
void stfun_push(stfun_t *A, dmat *opd);
dmat *stfun_finalize(stfun_t *A);
dmat* stfun_kolmogorov(loc_t *loc, double r0);
dmat *vkcov(long nx, double dx, double r0, double L0);
#endif
