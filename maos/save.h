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
#ifndef AOS_MAOS_SAVE_H
#define AOS_MAOS_SAVE_H
#include "common.h"
void save_wfsgrad(sim_t *simu);
void save_pistat(sim_t *simu);
void save_gradol(sim_t *simu);
void save_recon(sim_t *simu);
void save_dmproj(sim_t *simu);
void save_dmreal(sim_t *simu);
void draw_dm(const parms_t *parms, const recon_t *recon, const dcell *ac, int modal, const char *title, const char *type);
void draw_dm_lo(sim_t *simu, dcell *merr, const char *title, const char *type);
void plot_gradoff(sim_t *simu, int iwfs);
void plot_psf(ccell* psf2s, const char* psfname, int type, int ievl, dmat* wvl, int zlog, real psfmin);
#endif
