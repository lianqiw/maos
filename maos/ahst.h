/*
  Copyright 2009-2020 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
   \file ahst.h

   Contains functions to setup NGS modes and reconstructor
   using AHST for one or more DMs.  Use parms->wfsr instead of parms->wfs for wfs
   information, which hands GLAO mode correctly.

   Notice that update of this file may require GPU code update accordingly
*/
#ifndef AOS_AHST_H
#define AOS_AHST_H

void ngsmod2science(dmat *iopd, const loc_t *loc, const NGSMOD_T *ngsmod, 
		    real thetax, real thetay,
		    const real *mod, real alpha);
void setup_ngsmod_prep(const PARMS_T *parms, RECON_T *recon, 
		       const APER_T *aper, const POWFS_T* powfs);

void setup_ngsmod_recon(const PARMS_T *parms, RECON_T *recon);

void calc_ngsmod_dot(real *pttr_out, real *pttrcoeff_out,
		     real *ngsmod_out,
		     const PARMS_T *parms, const NGSMOD_T *ngsmod, const APER_T *aper, 
		     const real *opd, int ievl);
void calc_ngsmod_post(real *pttr_out, real *pttrcoeff_out, real *ngsmod_out,
		      real tot, const real *coeff,  const NGSMOD_T *ngsmod, 
		      const APER_T *aper,real thetax, real thetay);
void ngsmod_free(NGSMOD_T *ngsmod);
void remove_dm_ngsmod(SIM_T *simu, dcell *dmerr);
#endif
