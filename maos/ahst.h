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



#ifndef AOS_AHST_H
#define AOS_AHST_H

void ngsmod_prep(const parms_t *parms, recon_t *recon,  const aper_t *aper);
void ngsmod_setup(const parms_t *parms, recon_t *recon);
void ngsmod_dot(real *pttr_out, real *pttrcoeff_out, real *ngsmod_out, const parms_t *parms, const ngsmod_t *ngsmod, const aper_t *aper, const real *opd, int ievl);
int  ngsmod_dot_post(real *pttr_out, real *pttrcoeff_out, real *ngsmod_out, real tot, const real *coeff,  const ngsmod_t *ngsmod, const aper_t *aper,real thetax, real thetay);
void ngsmod_free(ngsmod_t *ngsmod);
void ngsmod_split(dcell **Merr, sim_t *simu, dcell *dmerr);
void ngsmod_remove(sim_t *simu, dcell *dmerr);
void ngsmod2science(dmat *iopd, const loc_t *loc, const ngsmod_t *ngsmod, real thetax, real thetay, const real *mod, real alpha);
#endif
