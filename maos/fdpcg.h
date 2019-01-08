/*
  Copyright 2009-2019 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef AOS_RECON_FDPCG
#define AOS_RECON_FDPCG
#include "common.h"
FDPCG_T *fdpcg_prepare(const PARMS_T *parms, const RECON_T *recon, 
		       const POWFS_T *powfs, mapcell *atm);
void fdpcg_precond(dcell **xout, const void *A, const dcell *xin);
void fdpcg_free(FDPCG_T *fdpcg);
#endif
