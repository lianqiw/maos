/*
  Copyright 2009-2025 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
/**
   \file maos/fdpcg.h

   Fourier Domain Preconditioner for Tomography step.

*/
#ifndef AOS_RECON_FDPCG
#define AOS_RECON_FDPCG
#include "common.h"
fdpcg_t *fdpcg_prepare(const parms_t *parms, const recon_t *recon, 
		       const powfs_t *powfs, mapcell *atm);
void fdpcg_precond(dcell **xout, const void *A, const dcell *xin);
void fdpcg_free(fdpcg_t *fdpcg);
#endif
