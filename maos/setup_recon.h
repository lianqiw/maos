/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#ifndef SETUP_RECON_H
#define SETUP_RECON_H
#include "common.h"
//Called by maos.c
RECON_T *setup_recon_prep(const PARMS_T *parms, const APER_T *aper, const POWFS_T *powfs);
void setup_recon_prep2(RECON_T *recon, const PARMS_T *parms, const APER_T *aper, const POWFS_T *powfs);
void setup_recon_fit(RECON_T *recon, const PARMS_T *parms);
void setup_powfs_fit(POWFS_T *powfs, const RECON_T *recon, const PARMS_T *parms);
void free_fit(FIT_T *fit);
//Called by setup_recon.c
void setup_recon_tomo_prep(RECON_T *recon, const PARMS_T *parms);
void setup_recon_tomo_matrix(RECON_T *recon, const PARMS_T *parms);
void setup_recon_tomo_update(RECON_T *recon, const PARMS_T *parms);
void setup_recon_tomo(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs);
void setup_recon_lsr(RECON_T *recon, const PARMS_T *parms);
void setup_recon_mvm(const PARMS_T *parms, RECON_T *recon, POWFS_T *powfs);
void setup_recon(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs);
void setup_recon_update(RECON_T *recon, const PARMS_T *parms, POWFS_T *powfs);
void setup_recon_post(RECON_T *recon, const PARMS_T *parms, const APER_T *aper);
void free_recon(const PARMS_T *parms, RECON_T *recon);
void free_recon_unused(const PARMS_T *parms, RECON_T *recon);

void test_recon_GX(RECON_T *recon, const PARMS_T *parms, 
		   const POWFS_T *powfs);
void test_recon_GA(RECON_T *recon, const PARMS_T *parms, 
		   const POWFS_T *powfs);
void setup_recon_mvst(RECON_T *recon, const PARMS_T *parms);
void setup_recon_dmttr(RECON_T *recon, const PARMS_T *parms);
void setup_recon_dither_dm(RECON_T *recon, const POWFS_T *powfs, const PARMS_T *parms);
#endif
