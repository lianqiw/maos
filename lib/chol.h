/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_CHOL_H
#define AOS_CHOL_H
#ifndef DLONG
#define DLONG
#endif
#include "dsp.h"
#include "dmat.h"
#include "cmat.h"

/* assume large file support.  If problems occur, compile with -DNLARGEFILE */
#include "../external/cholmod/Include/cholmod_io64.h"
/* define UF_long */
#include "../external/cholmod/Include/UFconfig.h"
#include "../external/cholmod/Include/cholmod_config.h"
#include "../external/cholmod/Include/cholmod_core.h"

typedef struct spchol{
    cholmod_factor *L;
    cholmod_common c;
}spchol;
spchol* chol_factorize(dsp *A_in);
void chol_solve(dmat **x, spchol *A, const dmat *y);
void chol_free_do(spchol *A);
#define chol_free(A) ({chol_free_do(A);A=NULL;})
//dsp* chol_sp(spchol *A, int keep);
void chol_save(spchol *A, const char *format,...) CHECK_ARG(2);
void chol_convert(dsp **Cs, long **Cp, spchol *A, int keep);
void chol_solve_lower(dmat **x, dsp *A, long *p, const dmat *y);
void chol_solve_upper(dmat **x, dsp *A, long *perm, const dmat *y);
#define M_CHOL 0x6501
#endif
