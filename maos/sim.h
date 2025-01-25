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
#ifndef AOS_SIM_H
#define AOS_SIM_H
/**
   \file sim.h

   Contains main simulation blocks.
*/
void *perfevl_ievl(thread_t *info);
void *perfevl(sim_t *simu);
void plot_psf(ccell *psf2s, const char *psfname, int type, int ievl, dmat *wvl, int zlog, real psfmin);
void prep_cachedm(sim_t *simu);
void calc_cachedm(sim_t *simu);
void filter_dm(sim_t *simu);
void update_dm(sim_t *simu);
void* wfsgrad(sim_t *simu);
void* wfsints(thread_t *thread_data);
void wfs_ideal_atm(sim_t *simu, dmat *opd, int iwfs, real alpha);
void* wfsgrad_iwfs(thread_t *info);
void* wfsgrad_post(thread_t *info);
void addlow2dm(dcell **dmval, const sim_t *simu, const dcell *low_val, real gain);
void filter_fsm(sim_t *simu);
#endif
