/*
  Copyright 2009-2021 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
recon_t *setup_recon_prep(const parms_t *parms, const aper_t *aper, const powfs_t *powfs);
void setup_recon_prep2(recon_t *recon, const parms_t *parms, const aper_t *aper, const powfs_t *powfs);
void setup_recon_fit(recon_t *recon, const parms_t *parms);
void setup_powfs_fit(powfs_t *powfs, const recon_t *recon, const parms_t *parms);
void free_powfs_fit(powfs_t *powfs, const parms_t *parms);
void free_fit(fit_t *fit, int nfit);
//Called by setup_recon.c
void setup_recon_tomo_prep(recon_t *recon, const parms_t *parms);
void setup_recon_tomo_matrix(recon_t *recon, const parms_t *parms);
void setup_recon_update_cn2(recon_t *recon, const parms_t *parms);
void setup_recon_tomo(recon_t *recon, const parms_t *parms, powfs_t *powfs);
void setup_recon_lsr(recon_t *recon, const parms_t *parms);
void setup_recon_mvm(const parms_t *parms, recon_t *recon, powfs_t *powfs);
void setup_recon_control(recon_t *recon, const parms_t *parms, powfs_t *powfs);
void setup_recon_post(recon_t *recon, const parms_t *parms, const aper_t *aper);
void free_recon(const parms_t *parms, recon_t *recon);
void free_recon_unused(const parms_t *parms, recon_t *recon);

void setup_recon_mvst(recon_t *recon, const parms_t *parms);
void setup_recon_dmttr(recon_t *recon, const parms_t *parms);
void setup_recon_dither_dm(recon_t *recon, const powfs_t *powfs, const parms_t *parms);
#endif
