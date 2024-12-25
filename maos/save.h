/*
  Copyright 2009-2024 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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



#include "common.h"
void save_wfsgrad(sim_t *simu);
void save_pistat(sim_t *simu);
void save_gradol(sim_t *simu);
void save_recon(sim_t *simu);
void save_dmproj(sim_t *simu);
void save_dmreal(sim_t *simu);
void draw_dm(const parms_t *parms, const recon_t *recon, const dcell *ac, int modal, const char *title, const char *type);
void draw_dm_lo(sim_t *simu, dcell *merr, const char *title, const char *type);
