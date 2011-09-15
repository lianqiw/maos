/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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


dcell *ngsmod_hm_ana(const PARMS_T *parms, RECON_T *recon, APER_T *aper);

dcell *ngsmod_hm_accphi(const PARMS_T *parms, RECON_T *recon, APER_T *aper);

void ngsmod2dm(dcell **dmc, const RECON_T *recon, const dcell *M, double gain);

void ngsmod2science(dmat *iopdevl, const PARMS_T *parms,
		    const RECON_T *recon, const APER_T *aper,
		    const double *mod, int ievl, double alpha);

void setup_ngsmod(const PARMS_T *parms, RECON_T *recon, 
		  APER_T *aper, POWFS_T* powfs);

void calc_ngsmod_dot(double *pttr_out, double *pttrcoeff_out,
		     double *ngsmod_out,
		     const PARMS_T *parms,
		     const RECON_T *recon, const APER_T *aper, 
		     const double *opd, int ievl);
void remove_dm_ngsmod(SIM_T *simu, dcell *dmerr);
void remove_dm_tt(SIM_T *simu, dcell *dmerr);
#endif
